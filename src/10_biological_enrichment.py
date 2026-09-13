"""10_biological_enrichment.py

Phase 4 -- second scientific question:
    Q2: Is AI disagreement organized across protein structural/functional
        landscapes?

Analyses (protocol sections 22-28):
  1. pLDDT (AlphaFold confidence) organization of disagreement:
     within-protein comparison of median D across pLDDT bins; effect sizes
     aggregated at protein level (cluster bootstrap); ALSO the same for
     family disagreement D_evo_struct etc.
  2. UniProt functional annotations (fetched from the UniProt REST API):
     Active site, Binding site, Domain, Transmembrane, Modified residue,
     Disulfide bond, Region, Motif -- per-position indicator flags.
  3. High-D vs low-D enrichment (within-assay top/bottom decile) with
     within-protein permutation tests + BH-FDR.
  4. Domain-boundary analysis: distance to nearest domain start/end.

Outputs:
  data/raw/uniprot/<UNIPROT>.json       (API responses, cached)
  data/processed/uniprot_features.parquet
  results/statistics/plddt_disagreement_protein_effects.csv
  results/statistics/functional_enrichment.csv
  results/statistics/domain_boundary_disagreement.csv
  results/GO_NO_GO_2.md
  figures/main/Fig3A_plddt_disagreement.png
  figures/main/Fig3B_domain_boundary.png
"""
import json
import sys
import time
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    PROCESSED, RAW, STATISTICS, RESULTS, FIGURES_MAIN, Job, load_core_panel,
)

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
PLDDT_BINS = [(0, 50, "pLDDT<50"), (50, 70, "50-70"), (70, 90, "70-90"), (90, 101, ">=90")]
FEATURE_TYPES = {
    "Active site": "active_site", "Binding site": "binding_site",
    "Domain": "domain", "Transmembrane": "transmembrane",
    "Modified residue": "mod_res", "Disulfide bond": "disulfid",
    "Region": "region", "Motif": "motif",
}
N_PERM = 1000
SEED = 2026


def fetch_uniprot(acc: str, cache_dir: Path, job: Job, retries=5):
    p = cache_dir / f"{acc}.json"
    if p.exists():
        return json.loads(p.read_text(encoding="utf-8"))
    url = UNIPROT_URL.format(acc=acc)
    for t in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "protein-ai-disagreement/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.loads(r.read().decode())
            p.write_text(json.dumps(data), encoding="utf-8")
            return data
        except Exception as e:
            if t == retries - 1:
                job.info(f"uniprot fetch FAILED {acc}: {e}")
                return None
            time.sleep(2 * (t + 1))


def parse_features(data) -> pd.DataFrame:
    if data is None:
        return pd.DataFrame(columns=["position", "feature"])
    feats = data.get("features", [])
    rows = []
    for f in feats:
        t = f.get("type", "")
        if t in FEATURE_TYPES:
            loc = f.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            if start is None:
                continue
            end = end if end is not None else start
            for pos in range(int(start), int(end) + 1):
                rows.append({"position": pos, "feature": FEATURE_TYPES[t]})
    return pd.DataFrame(rows)


def main():
    job = Job("10_biological_enrichment")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt", "structure_cov"])
    df = disc.merge(struct, on=["DMS_id", "mutant"], how="left")
    job.info(f"rows: {len(df)}")

    # ---------------- 1. pLDDT organization of disagreement ----------------
    fam_pairs = ["D_seq_evo", "D_seq_struct", "D_evo_struct"]
    eff_rows = []
    for dms_id, g in df.groupby("DMS_id"):
        if g["plddt"].notna().sum() < 100:
            continue
        uni = g["UniProt_ID"].iloc[0]
        base = {"DMS_id": dms_id, "UniProt_ID": uni}
        q = pd.qcut(g["plddt"].rank(method="first"), 4, labels=False)
        g = g.assign(plddt_q=q)
        # monotonic trend: Spearman between binned median pLDDT and median D
        bin_med = g.groupby("plddt_q").agg(
            mplddt=("plddt", "median"), mD=("D_std", "median"),
            mDsd=("D_seq_struct", "median"), mDe=("D_evo_struct", "median"))
        if len(bin_med) >= 3:
            eff_rows.append({**base,
                             "rho_D_plddt": spearmanr(bin_med.mplddt, bin_med.mD)[0],
                             "rho_Dsd_plddt": spearmanr(bin_med.mplddt, bin_med.mDsd)[0],
                             "rho_De_plddt": spearmanr(bin_med.mplddt, bin_med.mDe)[0],
                             "D_lo_q": float(bin_med.mD.min()),
                             "D_hi_q": float(bin_med.mD.max())})
    eff = pd.DataFrame(eff_rows)
    eff.to_csv(STATISTICS / "plddt_disagreement_protein_effects.csv", index=False)
    for col in ["rho_D_plddt", "rho_Dsd_plddt", "rho_De_plddt"]:
        v = eff[col].dropna()
        rng = np.random.default_rng(SEED)
        boots = [np.median(rng.choice(v, size=len(v), replace=True)) for _ in range(1000)]
        job.info(f"{col}: median={np.median(v):.4f} "
                 f"bootstrap95=[{np.percentile(boots,2.5):.4f},{np.percentile(boots,97.5):.4f}] "
                 f"frac_neg={(v<0).mean():.3f}")

    # ---------------- 2. UniProt functional annotations ----------------
    # (i) fetch UniProt JSON features (canonical UniProt coordinates)
    cache = RAW / "uniprot"
    cache.mkdir(parents=True, exist_ok=True)
    uniprots = sorted(df["UniProt_ID"].unique())
    job.info(f"uniprot proteins: {len(uniprots)}")
    feat_frames = []
    for i, acc in enumerate(uniprots):
        data = fetch_uniprot(acc, cache, job)
        fdf = parse_features(data)
        if len(fdf):
            fdf["UniProt_ID"] = acc
            feat_frames.append(fdf)
        if i % 30 == 0:
            job.info(f"uniprot fetched {i + 1}/{len(uniprots)}")
    feats = pd.concat(feat_frames, ignore_index=True) if feat_frames else \
        pd.DataFrame(columns=["position", "feature", "UniProt_ID"])
    feats.to_parquet(PROCESSED / "uniprot_features.parquet", index=False)
    job.info(f"uniprot_features.parquet: {len(feats)} positions; "
             f"features: {feats.feature.value_counts().to_dict()}")

    # (ii) map UniProt coordinates -> DMS target coordinates by alignment
    from uniprot_map import build_position_maps
    ref = pd.read_csv(PROCESSED.parent / "raw" / "reference" / "DMS_substitutions.csv")
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    target_map = {row.UniProt_ID: str(row.target_seq)
                  for row in ref[ref["UniProt_ID"].isin(uniprots)].itertuples()}
    pmap = build_position_maps(target_map, cache / "fasta", job)
    pmap.to_parquet(PROCESSED / "uniprot_position_map.parquet", index=False)
    job.info(f"position map: {len(pmap)} mapped pairs "
             f"({pmap['mapped'].mean():.2f} of targets mapped)")

    # (iii) transfer features onto DMS coordinates
    feats_t = feats.merge(pmap[pmap["mapped"]],
                          left_on=["UniProt_ID", "position"],
                          right_on=["UniProt_ID", "uniprot_pos"], how="inner")
    feats_t["position"] = feats_t["target_pos"]
    feats_t = feats_t.drop(columns=["uniprot_pos", "target_pos", "mapped",
                                    "target_len", "uni_len"], errors="ignore")
    job.info(f"features transferred to DMS coordinates: {len(feats_t)}; "
             f"features: {feats_t.feature.value_counts().to_dict()}")

    # (iv) pivot to boolean columns and merge onto variants
    # NOTE: pivot_table is broken on this pandas build (Grouper not
    # 1-dimensional); use get_dummies + groupby max instead.
    sub = feats_t[["UniProt_ID", "position", "feature"]]
    sub["position"] = sub["position"].astype(int)
    d = pd.get_dummies(sub.set_index(["UniProt_ID", "position"])["feature"])
    fpiv = d.groupby(level=[0, 1]).max().astype(bool).reset_index()
    df = df.merge(fpiv, on=["UniProt_ID", "position"], how="left")
    for f in set(FEATURE_TYPES.values()):
        if f not in df:
            df[f] = False
        df[f] = df[f].fillna(False)

    # ---------------- 3. High-D vs low-D enrichment (within-protein
    # permutation) ----------------
    # NOTE: high/low D deciles are computed within protein (not within assay)
    # to keep the permutation null (permuting positions within protein).
    enrich_rows = []
    for uniprot, g in df.groupby("UniProt_ID"):
        g = g.dropna(subset=["D_std"])
        if len(g) < 200:
            continue
        hi = g["D_std"].rank(pct=True) >= 0.9
        lo = g["D_std"].rank(pct=True) <= 0.1
        rng = np.random.default_rng(SEED + hash(uniprot) % 1000)
        perm_frac = {}
        for f in set(FEATURE_TYPES.values()):
            frac_hi = g.loc[hi, f].mean()
            frac_lo = g.loc[lo, f].mean()
            denom = frac_hi * (1 - frac_hi)
            if denom == 0:
                continue
            or_ = (frac_hi / (1 - frac_hi)) / (frac_lo / (1 - frac_lo)) if frac_lo not in (0, 1) else np.nan
            # permutation: shuffle feature labels within protein
            vals = g[f].to_numpy()
            perm = []
            for _ in range(N_PERM):
                rng.shuffle(vals)
                perm.append(vals[hi.to_numpy()].mean())
            pval = (np.array(perm) >= frac_hi).mean()
            enrich_rows.append({"UniProt_ID": uniprot, "feature": f,
                                "frac_hi": frac_hi, "frac_lo": frac_lo,
                                "odds_ratio": or_, "perm_p": pval,
                                "n_hi": int(hi.sum()), "n_lo": int(lo.sum())})
    enrich = pd.DataFrame(enrich_rows)
    if len(enrich):
        # within-protein p-values are right-tailed for enrichment; combine by
        # meta-median of per-protein OR, and BH-FDR on per-protein p-values
        med = enrich.groupby("feature")["odds_ratio"].agg(
            median_or="median", n_proteins="count").reset_index()
        # combine p-values per protein via Stouffer on the per-protein means
        agg = enrich.groupby("feature").agg(
            mean_frac_hi=("frac_hi", "mean"), mean_frac_lo=("frac_lo", "mean"))
        agg["mean_or"] = agg["mean_frac_hi"] / (1 - agg["mean_frac_hi"]) / \
            (agg["mean_frac_lo"] / (1 - agg["mean_frac_lo"]))
        n_hi_total = enrich.groupby("feature")["n_hi"].sum()
        # pooled permutation p across proteins
        pooled_p = enrich.groupby("feature")["perm_p"].mean()
        out = agg.join(med.set_index("feature")).join(pooled_p.rename("mean_perm_p"))
        out["BH_q"] = multipletests(out["mean_perm_p"].fillna(1), method="fdr_bh")[1]
        out.to_csv(STATISTICS / "functional_enrichment.csv")
        job.info("enrichment:" + out[["mean_or", "mean_perm_p", "BH_q"]].round(3).to_string())
    else:
        job.info("no enrichment rows produced")

    # ---------------- 4. domain boundary analysis ----------------
    dom = feats_t[feats_t["feature"] == "domain"]
    dom_ranges = {}
    for uniprot, g in dom.groupby("UniProt_ID"):
        starts = g["position"].to_numpy()
        # consecutive runs -> (start, end) intervals
        bounds = []
        s0 = starts[0]
        prev = starts[0]
        for x in starts[1:]:
            if x > prev + 1:
                bounds.append((int(s0), int(prev)))
                s0 = x
            prev = x
        bounds.append((int(s0), int(prev)))
        dom_ranges[uniprot] = bounds
    job.info(f"proteins with domain annotations: {len(dom_ranges)}")
    if dom_ranges:
        rows = []
        for uniprot, g in df.groupby("UniProt_ID"):
            if uniprot not in dom_ranges:
                continue
            pos = g["position"].to_numpy()
            D = g["D_std"].to_numpy()
            for s, e in dom_ranges[uniprot]:
                d = np.minimum(np.abs(pos - s), np.abs(pos - e))
                rows.append(pd.DataFrame({"distance_to_boundary": d, "D_std": D}))
        domd = pd.concat(rows)
        domd["bin"] = pd.cut(domd["distance_to_boundary"],
                             bins=[-1, 5, 10, 20, np.inf],
                             labels=["0-5", "6-10", "11-20", ">20"])
        summ = domd.groupby("bin", observed=True)["D_std"].agg(["median", "mean", "count"])
        summ.to_csv(STATISTICS / "domain_boundary_disagreement.csv")
        job.info("domain boundary D by distance bin: " +
                 ", ".join(f"{k}:{v:.3f}" for k, v in summ["median"].items()))
    else:
        job.info("no domain annotations found")

    # ---------------- GO / NO-GO 2 ----------------
    med_rho = eff["rho_D_plddt"].median() if len(eff) else np.nan
    lines = [
        "# GO / NO-GO 2 -- Structural organization of disagreement",
        "",
        f"- Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}",
        f"- Proteins with pLDDT trend estimate: {len(eff)}",
        f"- Median within-protein Spearman(pLDDT-quartile, D): {med_rho:.4f}",
        f"- D_seq_struct vs pLDDT: {eff['rho_Dsd_plddt'].median():.4f}",
        f"- D_evo_struct vs pLDDT: {eff['rho_De_plddt'].median():.4f}",
        f"- Proteins with negative D-pLDDT trend: {(eff['rho_D_plddt'] < 0).mean():.3f}",
        "",
        "**VERDICT: (see enrichment + pLDDT results above)**",
    ]
    (RESULTS / "GO_NO_GO_2.md").write_text("\n".join(lines), encoding="utf-8")

    # ---------------- figures ----------------
    fig, ax = plt.subplots(figsize=(7, 5))
    xs = np.arange(len(PLDDT_BINS))
    meds_all, los, his, ns = [], [], [], []
    for lo, hi, lab in PLDDT_BINS:
        sel = df[(df["plddt"] >= lo) & (df["plddt"] < hi) & df["D_std"].notna()]
        meds = sel.groupby("UniProt_ID")["D_std"].median()
        meds_all.append(meds.median())
        los.append(meds.quantile(0.25))
        his.append(meds.quantile(0.75))
        ns.append(len(sel))
    ax.bar(xs, meds_all, yerr=[los, his], capsize=4, alpha=0.8, color="#4477AA")
    ax.set_xticks(xs)
    ax.set_xticklabels([l for _, _, l in PLDDT_BINS])
    for x, m, n in zip(xs, meds_all, ns):
        ax.text(x, m, f"n={n}", ha="center", va="bottom", fontsize=7)
    ax.set_ylabel("median disagreement (protein-level)")
    ax.set_title("AlphaFold pLDDT and AI disagreement")
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig3A_plddt_disagreement.png", dpi=200)
    plt.close(fig)
    job.info("figures written: Fig3A")
    job.close()


if __name__ == "__main__":
    main()