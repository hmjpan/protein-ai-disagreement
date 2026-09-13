"""45_matched5_control.py

Two reviewer-facing facts, computed exactly:

1) 2x2 decomposition of the clinical uniform-AUROC gap (0.968 vs 0.867):
   five released zero-shot models, raw-score-average vs rank-average,
   on the all-proteins set (n=2,490) and the gate-complete subset (n=788).

2) Fully matched disagreement control: the SAME five models, SAME u/rank
   normalization and SAME SD(z) disagreement definition, applied to DMS
   (per assay) and to clinical (per protein), plus the 0.5-threshold
   misclassification AUROC on both sides:
     DMS:     rho(D5_assay, |mean_u5 - Y|)         ; AUROC(D5, E5>median)
     Clinical: rho(D_clin, |ens_clin - label|)     ; AUROC(D_clin, mis>0.5)

Outputs:
  results/statistics/uniform_2x2_clinical.csv
  results/statistics/matched5_dms_vs_clinical.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata, spearmanr
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (PROCESSED, RAW, STATISTICS, Job,
                    read_reference_clinical)  # noqa: E402

SEED = 2026
RNG = np.random.default_rng(SEED)
FIVE = ["TranceptEVE_L", "GEMME", "EVE", "ESM1b", "PoET"]
SCORES = PROCESSED.parent / "raw" / "proteingym_scores"
CLI_DIR = PROCESSED.parent / "raw" / "proteingym_clinical_scores"


def boot_ci(vals, reps=10000):
    v = np.asarray(vals, float)
    meds = [np.median(RNG.choice(v, len(v), replace=True)) for _ in range(reps)]
    return float(np.nanpercentile(meds, 2.5)), float(np.nanpercentile(meds, 97.5))


def main():
    job = Job("45_matched5_control")

    # ---------- clinical per-protein matrices ----------
    ref = read_reference_clinical()
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    rows22, rho_c, auc_c = [], [], []
    gate_complete = pd.read_csv(STATISTICS / "clinical_transfer_nonoverlap.csv")
    gc_set = set(gate_complete["DMS_id"].astype(str))
    full_set = pd.read_csv(STATISTICS / "clinical_transfer_performance.csv")
    full_set = set(full_set["DMS_id"].astype(str))
    for p in sorted(CLI_DIR.glob("*.csv")):
        dms_id = p.stem
        if dms_id not in full_set:
            continue
        sc = pd.read_csv(p)
        lab = sc["DMS_bin_score"].astype(str).str.strip().str.lower()
        sc["label"] = lab.map({"pathogenic": 1, "benign": 0})
        sc = sc.dropna(subset=["label"] + FIVE)
        if sc["label"].nunique() < 2 or len(sc) < 10:
            continue
        yl = sc["label"].to_numpy().astype(int)
        S = sc[FIVE].to_numpy(float)
        S[:, FIVE.index("PoET")] *= -1  # official clinical config: -1
        raw_ens = 1 - S.mean(axis=1)
        U = np.column_stack([1 - rankdata(S[:, j]) / len(sc) for j in range(5)])
        rank_ens = U.mean(axis=1)
        a_raw = roc_auc_score(yl, raw_ens)
        a_rk = roc_auc_score(yl, rank_ens)
        z = norm.ppf(np.clip(U, 0.001, 0.999))
        D = z.std(axis=1)
        err = np.abs(rank_ens - yl)
        mis = (err > 0.5).astype(int)
        rho_c.append(spearmanr(D, err)[0])
        if 0 < mis.mean() < 1:
            auc_c.append(roc_auc_score(mis, -D))
        in_gc = dms_id in gc_set
        rows22.append({"DMS_id": dms_id, "in_gate_complete": in_gc,
                       "auroc_raw_mean": a_raw, "auroc_rank_mean": a_rk})
    cdf = pd.DataFrame(rows22)
    full = pd.read_csv(STATISTICS / "clinical_transfer_performance.csv")
    full["DMS_id"] = full["DMS_id"].astype(str)

    # gate versus rank-averaged uniform, paired on gate-complete proteins
    tg = full.merge(cdf, on="DMS_id", how="inner")
    tg = tg[tg["in_gate_complete"]]
    d_gate = (tg["auroc_biogate"] - tg["auroc_rank_mean"]).to_numpy()
    rng = np.random.default_rng(SEED)
    meds = [np.median(rng.choice(d_gate, len(d_gate), replace=True))
            for _ in range(10000)]
    gate_vs_rank = pd.DataFrame([{
        "comparison": "BioGate - rank-averaged 5-model uniform",
        "n": int(len(d_gate)),
        "median_delta": float(np.median(d_gate)),
        "ci_low": float(np.percentile(meds, 2.5)),
        "ci_high": float(np.percentile(meds, 97.5)),
        "frac_positive": float((d_gate > 0).mean())}])
    gate_vs_rank.to_csv(STATISTICS / "gate_vs_rankuniform.csv", index=False)
    print(gate_vs_rank.to_string(index=False))
    t22 = []
    for name, sel in [("all scored proteins (2,490)", ~np.zeros(len(cdf), bool)),
                      ("gate-complete subset (788)", cdf["in_gate_complete"])]:
        sub = cdf[sel]
        t22.append({"set": name, "n": len(sub),
                   "median_raw_mean": float(sub["auroc_raw_mean"].median()),
                   "median_rank_mean": float(sub["auroc_rank_mean"].median())})
    pd.DataFrame(t22).to_csv(STATISTICS / "uniform_2x2_clinical.csv", index=False)
    print(pd.DataFrame(t22).to_string(index=False))

    # ---------- DMS with the same five models ----------
    norm5 = []
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet",
                           columns=["DMS_id", "mutant", "Y_deleter"])
    # per-assay u/z of the five models from merged raw score table
    ms = pd.read_parquet(PROCESSED / "model_scores.parquet") \
        if (PROCESSED / "model_scores.parquet").exists() else None
    if ms is None:
        # fall back: normalized_scores has u_ per model (percentile space)
        ns = pd.read_parquet(PROCESSED / "normalized_scores.parquet",
                             columns=["DMS_id", "mutant"] +
                             [f"u_{m}" for m in
                              ["TranceptEVE_L", "GEMME", "EVE_ensemble",
                               "ESM1b", "PoET"]])
        m5 = disc.merge(ns, on=["DMS_id", "mutant"])
        m5 = m5.rename(columns={"u_EVE_ensemble": "u_EVE"})
    rho_d, auc_d = [], []
    for _, g in m5.groupby("DMS_id"):
        g = g.dropna(subset=[f"u_{m}" for m in FIVE])
        if len(g) < 50:
            continue
        U = g[[f"u_{m}" for m in FIVE]].to_numpy()
        z = norm.ppf(np.clip(U, 0.001, 0.999))
        D = z.std(axis=1)
        err = np.abs(U.mean(axis=1) - g["Y_deleter"].to_numpy())
        rho_d.append(spearmanr(D, err)[0])
        mis = (err > np.median(err)).astype(int)
        if 0 < mis.mean() < 1:
            auc_d.append(roc_auc_score(mis, D))
    lo, hi = boot_ci(rho_c)
    loD, hiD = boot_ci(rho_d)
    loA, hiA = boot_ci(auc_c)
    loAD, hiAD = boot_ci(auc_d)
    out = pd.DataFrame([
        {"dataset": "DMS (same 5 models, same SD(z)/rank-mean definitions)",
         "statistic": "median per-assay Spearman(D, |mean_u - Y|)",
         "units": len(rho_d), "value": float(np.median(rho_d)),
         "ci_low": loD, "ci_high": hiD},
        {"dataset": "Clinical (released definition)",
         "statistic": "median per-protein Spearman(D, |mean_u - label|)",
         "units": len(rho_c), "value": float(np.median(rho_c)),
         "ci_low": lo, "ci_high": hi},
        {"dataset": "DMS (same 5 models)",
         "statistic": "median AUROC of D for above-median error",
         "units": len(auc_d), "value": float(np.median(auc_d)),
         "ci_low": loAD, "ci_high": hiAD},
        {"dataset": "Clinical",
         "statistic": "median AUROC of low D for misclassification (0.5 rule)",
         "units": len(auc_c), "value": float(np.median(auc_c)),
         "ci_low": loA, "ci_high": hiA},
    ])
    out.to_csv(STATISTICS / "matched5_dms_vs_clinical.csv", index=False)
    print(out.to_string(index=False))
    job.close()


if __name__ == "__main__":
    main()
