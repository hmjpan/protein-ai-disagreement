"""41_twofamily_control.py

Reviewer-response: is the clinical failure attributable to the DOMAIN or to
the two-family (sequence + evolution) model composition available in the
clinical benchmark?

Matched control: recompute the disagreement-error association in DMS using
EXACTLY the clinical configuration -- disagreement D2 = |U_seq - U_evo|
(percentile family centroids, the two families present in the clinical
score files), ensemble error E2 = |mean(U_seq, U_evo) - Y| -- per assay,
same aggregation and bootstrap as Section 2.6. Then compute the identical
statistic on the clinical benchmark per protein.

Outputs: results/statistics/twofamily_dms_vs_clinical.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, TABLES, STATISTICS, Job, read_reference_clinical  # noqa: E402

SEED = 2026
RNG = np.random.default_rng(SEED)


def boot_ci(vals, reps=10000):
    v = np.asarray(vals, dtype=float)
    meds = np.array([np.median(RNG.choice(v, size=len(v), replace=True))
                     for _ in range(reps)])
    return float(np.nanpercentile(meds, 2.5)), float(np.nanpercentile(meds, 97.5))


def main():
    job = Job("41_twofamily_control")
    panel = pd.read_csv(TABLES / "model_panel.csv")
    core = panel[panel["panel"] == "core"]
    evo = core[core["family"] == "evolution"]["model"].tolist()
    sq = core[core["family"] == "single_seq"]["model"].tolist()
    job.info(f"families present in clinical: {len(sq)} single-seq, {len(evo)} evolution")

    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet")
    need = [f"u_{m}" for m in evo + sq]
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet",
                           columns=["DMS_id", "mutant", "Y_deleter", "D_std",
                                    "S_evolution", "S_single_seq", "S_structure"])
    d = disc.merge(norm[["DMS_id", "mutant"] + need], on=["DMS_id", "mutant"],
                   how="left")
    d["U_evo2"] = d[[f"u_{m}" for m in evo]].mean(axis=1)
    d["U_seq2"] = d[[f"u_{m}" for m in sq]].mean(axis=1)
    d["D2"] = (d["U_seq2"] - d["U_evo2"]).abs()
    d["E2"] = ((d["U_seq2"] + d["U_evo2"]) / 2 - d["Y_deleter"]).abs()
    # anchor: three-family configuration as in Section 2.6
    d["U3"] = (d["S_structure"].rank(pct=True) * 0)  # placeholder not used
    rows = []
    sub = d.dropna(subset=["D2", "E2"])
    rhos2 = [spearmanr(g["D2"], g["E2"])[0] for _, g in sub.groupby("DMS_id")]
    lo, hi = boot_ci(rhos2)
    rows.append({"dataset": "DMS (two-family, matched to clinical)",
                 "n_units": len(rhos2), "median_rho": float(np.median(rhos2)),
                 "ci_low": lo, "ci_high": hi,
                 "frac_positive": float(np.mean(np.array(rhos2) > 0))})

    # anchor: three-family total disagreement vs ensemble error (Section 2.6)
    norm2 = pd.read_parquet(PROCESSED / "normalized_scores.parquet")
    allcore = core["model"].tolist()
    s3 = disc.copy()
    s3["S_ens"] = np.mean(np.nan_to_num(norm2[[f"u_{m}" for m in allcore]]
                                        .to_numpy(), nan=np.nan), axis=1)
    s3["E3"] = (s3["S_ens"] - s3["Y_deleter"]).abs()
    s3 = s3.dropna(subset=["D_std", "E3"])
    rhos3 = [spearmanr(g["D_std"], g["E3"])[0] for _, g in s3.groupby("DMS_id")]
    lo3, hi3 = boot_ci(rhos3)
    rows.append({"dataset": "DMS (three-family total, Section 2.6 anchor)",
                 "n_units": len(rhos3), "median_rho": float(np.median(rhos3)),
                 "ci_low": lo3, "ci_high": hi3,
                 "frac_positive": float(np.mean(np.array(rhos3) > 0))})

    # clinical: identical two-family statistic per protein
    cli = pd.read_parquet(PROCESSED / "clinical_disagreement.parquet")
    ref = read_reference_clinical()
    cli["UniProt"] = cli["DMS_id"].map(dict(zip(ref["DMS_id"], ref["UniProt"]))) \
        if "UniProt" in ref.columns else cli["DMS_id"]
    rhosc = []
    for _, g in cli.groupby("DMS_id"):
        if g["label"].nunique() < 2 or len(g) < 25:
            continue
        e = (g["ens_clin"] - g["label"]).abs()
        r = spearmanr(g["D_clin"], e)[0]
        if np.isfinite(r):
            rhosc.append(r)
    lo_c, hi_c = boot_ci(rhosc)
    rows.append({"dataset": "Clinical (two-family, same definition)",
                 "n_units": len(rhosc), "median_rho": float(np.median(rhosc)),
                 "ci_low": lo_c, "ci_high": hi_c,
                 "frac_positive": float(np.mean(np.array(rhosc) > 0))})

    # scale-matched discriminative comparison: AUROC of disagreement for
    # "high error" (DMS, within-assay top half) vs misclassification (clinical)
    from sklearn.metrics import roc_auc_score
    aucs_d, aucs_c = [], []
    for _, g in sub.groupby("DMS_id"):
        e_med = g["E2"].median()
        y_bin = (g["E2"] > e_med).astype(int)
        if 0 < y_bin.mean() < 1:
            aucs_d.append(roc_auc_score(y_bin, g["D2"]))
    for _, g in cli.groupby("DMS_id"):
        if g["label"].nunique() < 2 or len(g) < 25:
            continue
        mis = ((g["ens_clin"] - g["label"]).abs() > 0.5).astype(int)
        if 0 < mis.mean() < 1:
            aucs_c.append(roc_auc_score(mis, -g["D_clin"]))
    lo_d, hi_d = boot_ci(aucs_d)
    lo_cc, hi_cc = boot_ci(aucs_c)
    rows.append({"dataset": "DMS: AUROC of two-family disagreement for "
                            "above-median error",
                 "n_units": len(aucs_d), "median_rho": float(np.median(aucs_d)),
                 "ci_low": lo_d, "ci_high": hi_d,
                 "frac_positive": float(np.mean(np.array(aucs_d) > 0.5))})
    rows.append({"dataset": "Clinical: AUROC of low two-family disagreement "
                            "for misclassification",
                 "n_units": len(aucs_c),
                 "median_rho": float(np.median(aucs_c)),
                 "ci_low": lo_cc, "ci_high": hi_cc,
                 "frac_positive": float(np.mean(np.array(aucs_c) > 0.5))})

    out = pd.DataFrame(rows)
    out.to_csv(STATISTICS / "twofamily_dms_vs_clinical.csv", index=False)
    print(out.to_string(index=False))
    job.close()


if __name__ == "__main__":
    main()
