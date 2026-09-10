"""14_clinical_validation.py

Phase 7 -- external clinical validation.

Uses the ProteinGym clinical substitution benchmark (2,525 proteins, ~63K
variants, Pathogenic/Benign labels) and its precomputed zero-shot model scores.

Clinical proteins are identified by RefSeq DMS_id; the clinical panel
available in the released score files is {TranceptEVE_L, GEMME, EVE, ESM1b,
PoET} -- no structure-based models. Analyses are therefore restricted to:
  1. Clinical disagreement (D_std over the available panel) and its
     association with classification difficulty and with pathogenicity.
  2. Transfer of a DMS-trained BioGate (2-expert gating: seq + evo) to
     clinical variants without retraining; AUROC vs uniform ensemble and best
     single expert.
  3. Honest limitations: no clinical structures (pLDDT feature neutralised),
     no UniProt functional annotations (flags set to 0), no RefSeq->UniProt
     overlap exclusion (reported as limitation).

Outputs:
  data/processed/clinical_disagreement.parquet
  results/statistics/clinical_disagreement_summary.csv
  results/statistics/clinical_transfer_performance.csv
  figures/main/Fig7A_clinical_disagreement.png
  figures/main/Fig7B_clinical_transfer.png
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from scipy.stats import norm, rankdata
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    PROCESSED, STATISTICS, FIGURES_MAIN, Job, load_core_panel,
    read_reference_clinical,
)

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
SEED = 2026
CLINICAL_SCORES_DIR = PROCESSED.parent / "raw" / "proteingym_clinical_scores"

HIDDEN = 32
EPOCHS = 12
BATCH = 4096
LR = 1e-3
ENT_PEN = 0.01
CONTEXT = ["plddt_z", "norm_pos", "log_len", "msa_shallow", "msa_deep"]
EXPERTS = ["S_single_seq", "S_evolution", "S_structure"]


class GateNet(nn.Module):
    def __init__(self, n_in, n_exp):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(n_in, HIDDEN), nn.ReLU(), nn.Linear(HIDDEN, n_exp))
        self.n_exp = n_exp

    def forward(self, x, s):
        w = torch.softmax(self.mlp(x), dim=1)
        return (w * s).sum(1), w


def main():
    job = Job("14_clinical_validation")
    ref = read_reference_clinical()
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    ref_map = ref.set_index("DMS_id")
    job.info(f"clinical reference proteins: {len(ref)}")

    score_files = sorted(CLINICAL_SCORES_DIR.glob("*.csv")) \
        if CLINICAL_SCORES_DIR.exists() else []
    job.info(f"clinical score files: {len(score_files)}")
    core = load_core_panel(job)
    core_names = set(core["model"])

    frames = []
    for p in score_files:
        dms_id = p.stem
        if dms_id not in ref_map.index:
            continue
        sc = pd.read_csv(p)
        keep = [c for c in sc.columns
                if c in {"mutant", "DMS_bin_score"} or c in core_names]
        sc = sc[keep].copy()
        sc["DMS_id"] = dms_id
        frames.append(sc)
    cli = pd.concat(frames, ignore_index=True)
    model_list = [c for c in cli.columns
                  if c in core_names and c != "mutant"]
    job.info(f"clinical variants: {len(cli)}; proteins: {cli.DMS_id.nunique()}; "
             f"clinical core panel: {model_list}")

    lbl = cli["DMS_bin_score"].astype(str).str.strip().str.lower()
    cli["label"] = lbl.map({"pathogenic": 1, "benign": 0})
    cli = cli[cli["label"].notna()].copy()
    cli["label"] = cli["label"].astype(int)
    job.info(f"labeled: {len(cli)} (pathogenic={int((cli.label == 1).sum())}, "
             f"benign={int((cli.label == 0).sum())})")

    # ---------------- clinical disagreement ----------------
    rows = []
    for dms_id, g in cli.groupby("DMS_id"):
        sub = g[model_list].to_numpy(dtype=float)
        if sub.shape[1] < 3:
            continue
        u = np.full(sub.shape, np.nan)
        for j in range(sub.shape[1]):
            s = sub[:, j]
            ok = ~np.isnan(s)
            if ok.sum() > 5:
                r = rankdata(s[ok], method="average") / ok.sum()
                u[ok, j] = 1 - r  # released scores: higher = fitness
        z = norm.ppf(np.clip(u, 0.001, 0.999))
        D = np.nanstd(z, axis=1)
        ens = np.nanmean(u, axis=1)  # higher = more deleterious
        rows.append(pd.DataFrame({"DMS_id": dms_id,
                                  "mutant": g["mutant"].to_numpy(),
                                  "D_clin": D, "ens_clin": ens,
                                  "label": g["label"].to_numpy()}))
    cdisc = pd.concat(rows, ignore_index=True)
    cdisc.to_parquet(PROCESSED / "clinical_disagreement.parquet", index=False)
    job.info(f"clinical disagreement computed on {len(cdisc)} variants; "
             f"median D = {cdisc.D_clin.median():.4f}")

    diff_rows = []
    for dms_id, g in cdisc.groupby("DMS_id"):
        if g.label.nunique() < 2 or len(g) < 20:
            continue
        thr = _best_threshold(g["ens_clin"].to_numpy(), g["label"].to_numpy())
        pred = (g["ens_clin"].to_numpy() >= thr).astype(int)
        correct = pred == g["label"].to_numpy()
        diff_rows.append({"DMS_id": dms_id,
                          "D_correct": g.loc[correct, "D_clin"].median(),
                          "D_incorrect": g.loc[~correct, "D_clin"].median(),
                          "n": len(g), "n_incorrect": int((~correct).sum())})
    diff = pd.DataFrame(diff_rows)
    diff.to_csv(STATISTICS / "clinical_disagreement_summary.csv", index=False)
    if len(diff):
        job.info(f"median D incorrect={diff.D_incorrect.median():.4f} vs "
                 f"correct={diff.D_correct.median():.4f}; n_proteins={len(diff)}; "
                 f"frac D_inc>D_cor: {(diff.D_incorrect > diff.D_correct).mean():.3f}")
    med_p = cdisc.loc[cdisc.label == 1, "D_clin"].median()
    med_b = cdisc.loc[cdisc.label == 0, "D_clin"].median()
    job.info(f"median D pathogenic={med_p:.4f} benign={med_b:.4f}")

    # ---------------- BioGate transfer (train on DMS, test clinical) ----------------
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
    dms = dms.dropna(subset=EXPERTS + ["Y_deleter", "norm_pos", "log_len"])

    # clinical experts: seq + evo only (percentile u scale)
    fam_avail = [c for c in EXPERTS
                 if any(core.loc[core.model == m, "family"].iloc[0] ==
                        c.split("_", 1)[1] for m in model_list)]
    job.info(f"clinical expert families: {fam_avail}")
    if len(fam_avail) < 2:
        job.info("need >=2 expert families for transfer -- abort")
        job.close()
        return

    # DMS training side: percentile-space family centroids (U in [0,1],
    # same scale as Y), computed from released u scores
    core_names = [c for c in EXPERTS]
    norm_u = pd.read_parquet(PROCESSED / "normalized_scores.parquet",
                             columns=["DMS_id", "mutant"] +
                             [f"u_{m}" for m in core["model"]])
    dms = dms.merge(norm_u, on=["DMS_id", "mutant"], how="left")
    fam_u = {f: [f"u_{m}" for m in core["model"]
                 if core.loc[core.model == m, "family"].iloc[0] == f]
             for f in ["single_seq", "evolution", "structure"]}
    for f, cols in fam_u.items():
        dms[f"U_{f}"] = dms[cols].mean(axis=1, skipna=True)
    fam_avail = [f"U_{f}" for f in ["single_seq", "evolution"]
                 if f"U_{f}" in dms.columns]

    X = dms[CONTEXT].to_numpy(dtype=np.float32)
    S = dms[fam_avail].to_numpy(dtype=np.float32)
    y = dms["Y_deleter"].to_numpy(dtype=np.float32)
    model = GateNet(X.shape[1], S.shape[1])
    opt = torch.optim.Adam(model.parameters(), lr=LR)
    lossf = nn.MSELoss()
    Xt = torch.from_numpy(X); St = torch.from_numpy(S); yt = torch.from_numpy(y)
    n = len(X)
    for ep in range(EPOCHS):
        model.train()
        perm = np.random.default_rng(SEED + ep).permutation(n)
        for i in range(0, n, BATCH):
            idx = perm[i:i + BATCH]
            opt.zero_grad()
            pred, w = model(Xt[idx], St[idx])
            loss = lossf(pred, yt[idx]) + ENT_PEN * (-(w * (w + 1e-9).log()).sum(1).mean())
            loss.backward()
            opt.step()
    model.eval()
    job.info("BioGate trained on full DMS set (2-expert gating)")

    # clinical context (build from cli so model scores are present)
    cli_ctx = []
    for dms_id, g in cli.groupby("DMS_id"):
        rrow = ref_map.loc[dms_id]
        g = g.copy()
        g["seq_len"] = int(rrow["MSA_len"]) if pd.notna(rrow["MSA_len"]) else 1
        msa_len = int(rrow["MSA_len"]) if pd.notna(rrow["MSA_len"]) else 0
        g["norm_pos"] = pd.to_numeric(g["mutant"].str.extract(r"(\d+)")[0]) / g["seq_len"]
        g["log_len"] = np.log(g["seq_len"])
        g["msa_shallow"] = int(msa_len < 100)
        g["msa_deep"] = int(msa_len > 1000)
        g["plddt_z"] = 0.0  # no clinical structures available
        cli_ctx.append(g)
    cctx = pd.concat(cli_ctx, ignore_index=True)
    cctx = cctx.dropna(subset=["norm_pos", "log_len"])
    job.info(f"clinical transfer set: {len(cctx)} variants, "
             f"{cctx.DMS_id.nunique()} proteins")

    m_map = {m: i for i, m in enumerate(model_list)}
    cS_raw = cctx[model_list].to_numpy(dtype=float)
    # family centroids on clinical (deleter) space
    for f in fam_avail:
        fam = f.split("_", 1)[1]
        idx = [m_map[m] for m in model_list
               if core.loc[core.model == m, "family"].iloc[0] == fam]
        cS_u = np.full(cS_raw.shape, np.nan)
        for j in range(cS_raw.shape[1]):
            s = cS_raw[:, j]
            ok = ~np.isnan(s)
            if ok.sum() > 5:
                cS_u[ok, j] = 1 - rankdata(s[ok], method="average") / ok.sum()
        cctx[f] = np.nanmean(cS_u[:, idx], axis=1)
    cctx = cctx.dropna(subset=fam_avail)

    cX = cctx[CONTEXT].to_numpy(dtype=np.float32)
    cS2 = cctx[fam_avail].to_numpy(dtype=np.float32)
    with torch.no_grad():
        bio_pred, _ = model(torch.from_numpy(cX), torch.from_numpy(cS2))
        bio_pred = bio_pred.numpy()

    ens_pred = 1 - cctx[model_list].mean(axis=1)  # flip to deleter space
    best = None
    best_name = None
    for m in model_list:
        r = []
        for dms_id, g in cctx.groupby("DMS_id"):
            if g.label.nunique() == 2 and len(g) >= 10:
                s = g[m].to_numpy()
                ok = ~np.isnan(s)
                if ok.sum() >= 10:
                    r.append(roc_auc_score(g["label"].to_numpy()[ok], -s[ok]))
        if r and (best is None or np.median(r) > best):
            best = np.median(r)
            best_name = m
    best_pred = 1 - cctx[best_name].to_numpy() if best_name else ens_pred
    job.info(f"clinical best single model: {best_name}")

    trans_rows = []
    for dms_id, g in cctx.groupby("DMS_id"):
        if g.label.nunique() < 2 or len(g) < 10:
            continue
        yl = g["label"].to_numpy()
        trans_rows.append({"DMS_id": dms_id, "n": len(g),
                           "auroc_uniform": _auroc(ens_pred[g.index], yl),
                           "auroc_biogate": _auroc(bio_pred[g.index], yl),
                           "auroc_best": _auroc(best_pred[g.index], yl)})
    trans = pd.DataFrame(trans_rows)
    trans.to_csv(STATISTICS / "clinical_transfer_performance.csv", index=False)
    if len(trans):
        job.info(f"clinical transfer AUROC (protein-level median, n={len(trans)}):")
        for c in ["auroc_uniform", "auroc_biogate", "auroc_best"]:
            job.info(f"  {c}: {trans[c].median():.4f}")

    # ---------------- figures ----------------
    if len(diff):
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(diff["D_correct"], diff["D_incorrect"], s=12, alpha=0.6)
        lim = [0, max(diff[["D_correct", "D_incorrect"]].max().max(), 1)]
        ax.plot(lim, lim, "r--", lw=1)
        ax.set_xlabel("median disagreement (correctly classified)")
        ax.set_ylabel("median disagreement (misclassified)")
        fig.tight_layout()
        fig.savefig(FIGURES_MAIN / "Fig7A_clinical_disagreement.png", dpi=200)
        plt.close(fig)
    if len(trans):
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.boxplot([trans["auroc_uniform"], trans["auroc_biogate"],
                    trans["auroc_best"]], labels=["uniform", "biogate", "best single"])
        ax.set_ylabel("protein-level AUROC")
        fig.tight_layout()
        fig.savefig(FIGURES_MAIN / "Fig7B_clinical_transfer.png", dpi=200)
        plt.close(fig)
    job.info("figures written: Fig7A, Fig7B")
    job.close()


def _best_threshold(pred, y):
    from sklearn.metrics import roc_curve
    fpr, tpr, thr = roc_curve(y, pred)
    j = tpr - fpr
    return thr[np.argmax(j)]


def _auroc(pred, y):
    try:
        return roc_auc_score(y, pred)
    except Exception:
        return np.nan


if __name__ == "__main__":
    main()