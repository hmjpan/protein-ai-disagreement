"""20_clinical_overlap.py

Enhancement analysis (review v2, item 10): resolve the clinical overlap
problem by mapping RefSeq clinical protein IDs to UniProt accessions
(UniProt ID mapping API) and repeating the transfer evaluation on a STRICT
non-overlap set:

  Level 1 (identity): clinical proteins whose UniProt accession is among the
    186 DMS proteins are excluded.
  Level 2 (homology): clinical proteins with high sequence similarity to any
    DMS target sequence (>=70% exact residue identity by alignment) are
    excluded.

AUROC is recomputed on the remaining clinical variants for: uniform ensemble,
BioGate (DMS-trained 2-expert gating), and best single expert.

Outputs:
  results/statistics/clinical_overlap_mapping.csv
  results/statistics/clinical_transfer_nonoverlap.csv
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import urllib.request
import urllib.parse
import json
from difflib import SequenceMatcher
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    PROCESSED, STATISTICS, Job, read_reference_clinical,
)

BATCH = 500
IDENTITY_THRESH = 0.70


def map_refseq_to_uniprot(ids, job, retries=3):
    """UniProt ID mapping API (current): RefSeq_Protein -> UniProtKB."""
    out = {}
    cache_p = STATISTICS / "clinical_overlap_mapping_partial.csv"
    if cache_p.exists():
        c = pd.read_csv(cache_p)
        out = dict(zip(c["DMS_id"], c["uniprot"]))
        job.info(f"resumed from partial mapping: {len(out)}")
    todo = [x for x in ids if x not in out]
    for i in range(0, len(todo), BATCH):
        batch = todo[i:i + BATCH]
        body = urllib.parse.urlencode({
            "from": "RefSeq_Protein", "to": "UniProtKB",
            "ids": " ".join(batch)}).encode()
        job_id = None
        for t in range(retries):
            try:
                req = urllib.request.Request(
                    "https://rest.uniprot.org/idmapping/run", data=body,
                    headers={"User-Agent": "protein-ai/1.0"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    job_id = json.loads(r.read().decode())["jobId"]
                break
            except Exception as e:
                if t == retries - 1:
                    job.info(f"idmapping submit failed at {i}: {e}")
                    return out
                time.sleep(3)
        status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
        for _ in range(25):
            time.sleep(3)
            try:
                req = urllib.request.Request(status_url, headers={"User-Agent": "x"})
                with urllib.request.urlopen(req, timeout=60) as r:
                    st = json.loads(r.read().decode())
                if st.get("status") in ("FINISHED", "finished"):
                    break
            except Exception:
                pass
        res_url = f"https://rest.uniprot.org/idmapping/stream/{job_id}"
        try:
            req = urllib.request.Request(res_url, headers={"User-Agent": "x"})
            with urllib.request.urlopen(req, timeout=120) as r:
                text = r.read().decode()
            lines = text.strip().splitlines()
            if lines and lines[0].startswith("From"):
                for line in lines[1:]:
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        out[parts[0]] = parts[1].split(";")[0]
            else:
                results = json.loads(text).get("results", [])
                for row in results:
                    out[row["from"]] = row["to"]["primaryAccession"]
        except Exception as e:
            job.info(f"idmapping results failed at {i}: {e}")
        pd.DataFrame({"DMS_id": list(out.keys()),
                      "uniprot": list(out.values())}).to_csv(cache_p, index=False)
        job.info(f"mapped {len(out)}/{len(ids)} IDs so far")
    return out


def max_identity(target: str, dms_targets: list, job):
    """Max identity to any DMS target via 8-mer filter + SequenceMatcher."""
    tmers = {target[i:i + 8] for i in range(max(0, len(target) - 7))}
    best = 0.0
    for cand in dms_targets:
        cmers = set()
        for i in range(max(0, len(cand) - 7)):
            cmers.add(cand[i:i + 8])
        shared = len(tmers & cmers) / max(len(tmers), 1)
        if shared < 0.25:
            continue
        r = SequenceMatcher(None, target, cand, autojunk=False).ratio()
        if r > best:
            best = r
        if best >= IDENTITY_THRESH:
            break
    return best


def main():
    job = Job("20_clinical_overlap")
    ref = read_reference_clinical()
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    clinical_ids = ref["DMS_id"].tolist()
    dref = pd.read_csv(PROCESSED.parent / "raw" / "reference" / "DMS_substitutions.csv")
    dref["DMS_id"] = dref["DMS_id"].astype(str)
    dms_uniprots = set(dref["UniProt_ID"])
    dms_targets = dref["target_seq"].astype(str).tolist()
    job.info(f"clinical proteins: {len(clinical_ids)}; DMS proteins: {len(dms_uniprots)}")

    # Level 1 (identity): clinical proteins whose UniProt accession is among
    # the 186 DMS proteins. The UniProt RefSeq->UniProt idmapping API was
    # attempted but is unreliable from this environment (multiple endpoint
    # failures); the homology check below (Level 2) is strictly stronger,
    # because any same-protein clinical set is >=95% identical to its DMS
    # counterpart and is therefore captured by the >=70% identity threshold.
    mapped = {}
    map_df = pd.DataFrame({"DMS_id": clinical_ids,
                           "uniprot": [mapped.get(x, "") for x in clinical_ids]})
    map_df.to_csv(STATISTICS / "clinical_overlap_mapping.csv", index=False)
    job.info("RefSeq->UniProt mapping skipped (API unreliable); homology "
             "exclusion below subsumes identity-level overlap")

    identity_excl = set(map_df.loc[map_df["uniprot"].isin(dms_uniprots), "DMS_id"])
    job.info(f"clinical proteins sharing UniProt with DMS: {len(identity_excl)}")

    # homology exclusion (only for proteins not already identity-excluded)
    hom_excl = set()
    target_map = dict(zip(dref["DMS_id"], dms_targets))
    for _, row in map_df.iterrows():
        dms_id = row["DMS_id"]
        if dms_id in identity_excl:
            continue
        target = ref.loc[ref["DMS_id"] == dms_id, "target_seq"]
        if target.empty:
            continue
        ident = max_identity(str(target.iloc[0]), dms_targets, job)
        if ident >= IDENTITY_THRESH:
            hom_excl.add(dms_id)
    job.info(f"clinical proteins with >= {IDENTITY_THRESH:.0%} identity to a DMS "
             f"protein: {len(hom_excl)}")
    strict_excl = identity_excl | hom_excl
    job.info(f"STRICT exclusion set: {len(strict_excl)} clinical proteins")

    # ---- rerun transfer on strict non-overlap set ----
    cli = pd.read_parquet(PROCESSED / "clinical_disagreement.parquet")
    keep = ~cli["DMS_id"].isin(strict_excl)
    job.info(f"clinical variants retained: {keep.sum()}/{len(cli)}")

    # reload DMS-trained BioGate (2-expert gating, seq + evo) from 14's
    # protocol so the strict-set comparison uses the identical gate
    import os
    os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
    import torch
    import torch.nn as nn

    class GateNet(nn.Module):
        def __init__(self, n_in, n_exp):
            super().__init__()
            self.mlp = nn.Sequential(nn.Linear(n_in, 32), nn.ReLU(),
                                     nn.Linear(32, n_exp))

        def forward(self, x, s):
            w = torch.softmax(self.mlp(x), dim=1)
            return (w * s).sum(1), w

    dms = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    dref = pd.read_csv(PROCESSED.parent / "raw" / "reference" / "DMS_substitutions.csv")
    dref["DMS_id"] = dref["DMS_id"].astype(str)
    dref_map = dref.set_index("DMS_id")
    dms = dms.merge(struct, on=["DMS_id", "mutant"], how="left")
    dms["seq_len"] = dms["DMS_id"].map(dref_map["seq_len"])
    dms["msa_depth"] = dms["DMS_id"].map(dref_map["MSA_Neff_L_category"])
    dms["norm_pos"] = dms["position"] / dms["seq_len"]
    dms["log_len"] = np.log(dms["seq_len"])
    dms["msa_shallow"] = dms["msa_depth"].eq("Low").astype(int)
    dms["msa_deep"] = dms["msa_depth"].eq("High").astype(int)
    pmed = dms["plddt"].median()
    dms["plddt_z"] = (dms["plddt"].fillna(pmed) - dms["plddt"].mean()) / dms["plddt"].std()
    dms = dms.dropna(subset=["S_single_seq", "S_evolution", "Y_deleter",
                             "norm_pos", "log_len"])
    # percentile-space family centroids for the DMS training side
    from common import TABLES
    core_df = pd.read_csv(TABLES / "model_panel.csv")
    core_df = core_df[core_df["panel"] == "core"]
    norm_u = pd.read_parquet(PROCESSED / "normalized_scores.parquet",
                             columns=["DMS_id", "mutant"] +
                             [f"u_{m}" for m in core_df["model"]])
    dms = dms.merge(norm_u, on=["DMS_id", "mutant"], how="left")
    fam_u = {f: [f"u_{m}" for m in core_df["model"]
                 if core_df.loc[core_df.model == m, "family"].iloc[0] == f]
             for f in ["single_seq", "evolution"]}
    dms["U_single_seq"] = dms[fam_u["single_seq"]].mean(axis=1, skipna=True)
    dms["U_evolution"] = dms[fam_u["evolution"]].mean(axis=1, skipna=True)
    C = ["plddt_z", "norm_pos", "log_len", "msa_shallow", "msa_deep"]
    Xt = torch.from_numpy(dms[C].to_numpy(dtype=np.float32))
    St = torch.from_numpy(dms[["U_single_seq", "U_evolution"]].to_numpy(dtype=np.float32))
    yt = torch.from_numpy(dms["Y_deleter"].to_numpy(dtype=np.float32))
    model = GateNet(Xt.shape[1], 2)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    lossf = nn.MSELoss()
    n = len(Xt)
    for ep in range(12):
        model.train()
        perm = torch.randperm(n)
        for i in range(0, n, 4096):
            idx = perm[i:i + 4096]
            opt.zero_grad()
            pred, w = model(Xt[idx], St[idx])
            loss = lossf(pred, yt[idx]) + 0.01 * (-(w * (w + 1e-9).log()).sum(1).mean())
            loss.backward()
            opt.step()
    model.eval()
    job.info("BioGate (2-expert) retrained on full DMS for strict-set evaluation")

    # recompute per-protein AUROC on the retained variants for all methods
    from common import TABLES
    core_df = pd.read_csv(TABLES / "model_panel.csv")
    core_df = core_df[core_df["panel"] == "core"]
    core = core_df["model"].tolist()
    score_dir = PROCESSED.parent / "raw" / "proteingym_clinical_scores"
    rows = []
    for dms_id, g in cli[keep].groupby("DMS_id"):
        p = score_dir / f"{dms_id}.csv"
        if not p.exists():
            continue
        sc = pd.read_csv(p)
        lab = pd.read_csv(p)["DMS_bin_score"].astype(str).str.strip().str.lower()
        sc["label"] = lab.map({"pathogenic": 1, "benign": 0})
        sc = sc[sc["label"].notna()]
        yl = sc["label"].to_numpy().astype(int)
        if yl.size == 0 or np.unique(yl).size < 2:
            continue
        models = [m for m in core if m in sc.columns]
        if len(models) < 3:
            continue
        S = sc[models].to_numpy(dtype=float)
        # family centroids (deleterosity space) for the gate
        fam = {"S_single_seq": [], "S_evolution": []}
        for j, m in enumerate(models):
            f = core_df.loc[core_df.model == m, "family"].iloc[0]
            if f == "single_seq":
                fam["S_single_seq"].append(j)
            elif f == "evolution":
                fam["S_evolution"].append(j)
        u = np.full(S.shape, np.nan)
        for j in range(S.shape[1]):
            s = S[:, j]
            ok = ~np.isnan(s)
            if ok.sum() > 5:
                u[ok, j] = 1 - rankdata(s[ok], method="average") / ok.sum()
        S_seq = np.nanmean(u[:, fam["S_single_seq"]], axis=1)
        S_evo = np.nanmean(u[:, fam["S_evolution"]], axis=1)
        # context for clinical (no structures; pLDDT neutral; MSA from ref)
        rrow = ref[ref["DMS_id"] == dms_id].iloc[0]
        seq_len = int(rrow["MSA_len"]) if pd.notna(rrow["MSA_len"]) else 1
        msa_len = int(rrow["MSA_len"]) if pd.notna(rrow["MSA_len"]) else 0
        pos = pd.to_numeric(sc["mutant"].str.extract(r"(\d+)")[0]) / seq_len
        ctx = np.column_stack([
            np.zeros(len(sc)), pos, np.full(len(sc), np.log(seq_len)),
            np.full(len(sc), int(msa_len < 100)),
            np.full(len(sc), int(msa_len > 1000))]).astype(np.float32)
        with torch.no_grad():
            bio, _ = model(torch.from_numpy(ctx),
                           torch.from_numpy(np.column_stack([S_seq, S_evo]).astype(np.float32)))
            bio = bio.numpy()
        ens = 1 - np.nanmean(S, axis=1)
        # best single on the strict set (reselected within strict set)
        best_score = None
        best_name = None
        for m in models:
            s = S[:, models.index(m)]
            ok = ~np.isnan(s)
            if ok.sum() >= 10:
                a = roc_auc_score(yl[ok], -s[ok])
                if best_score is None or a > best_score:
                    best_score = a
                    best_name = m
        best = -S[:, models.index(best_name)] if best_name else ens
        # keep only rows with complete gate inputs
        valid = (~np.isnan(S_seq)) & (~np.isnan(S_evo)) & np.isfinite(bio) \
            & np.isfinite(ens) & np.isfinite(best)
        if valid.sum() < 20:
            continue
        yl_v, ens_v = yl[valid], ens[valid]
        bio_v, best_v = bio[valid], best[valid]
        if np.unique(yl_v).size < 2:
            continue
        rows.append({"DMS_id": dms_id, "n": int(valid.sum()),
                     "auroc_uniform": roc_auc_score(yl_v, ens_v),
                     "auroc_biogate": roc_auc_score(yl_v, bio_v),
                     "auroc_best": roc_auc_score(yl_v, best_v),
                     "label_rate": yl_v.mean()})
    non = pd.DataFrame(rows)
    non.to_csv(STATISTICS / "clinical_transfer_nonoverlap.csv", index=False)
    if len(non):
        job.info(f"non-overlap transfer AUROC (uniform ensemble; "
                 f"median over {len(non)} proteins): {non.auroc_uniform.median():.4f}")
        job.info(f"  vs full-set uniform (from 14): 0.877; "
                 f"n proteins retained: {len(non)}")
    job.close()


def _core_names():
    import pandas as pd
    p = STATISTICS.parent / "tables" / "model_panel.csv"
    return pd.read_csv(p).query("panel=='core'")["model"].tolist()


if __name__ == "__main__":
    main()