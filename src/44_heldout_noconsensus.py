"""44_heldout_noconsensus.py

Reviewer-response: is the held-out regime-profile Spearman of 0.94 driven
mainly by the two easy consensus regimes? Recompute held-out concordance
after dropping the consensus-damaging and consensus-tolerate regimes, and
report per-regime training vs held-out experimental Y medians with
protein-bootstrap CIs.

Same fold protocol as 33_heldout_regime.py (GroupKFold by protein, frozen
per-fold GMM K=6).

Outputs: results/statistics/heldout_noconsensus.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job  # noqa: E402

SEED = 2026
K = 6
S_COLS = ["S_single_seq", "S_evolution", "S_structure"]


def main():
    job = Job("44_heldout_noconsensus")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    keep = ~disc[S_COLS].isna().any(axis=1)
    d = disc[keep].reset_index(drop=True)
    X = d[S_COLS].to_numpy(dtype=float)
    prot = d["UniProt_ID"].to_numpy()
    y = d["Y_deleter"].to_numpy()

    recs = []
    gkf = GroupKFold(n_splits=5)
    for fold, (tr, te) in enumerate(gkf.split(X, y, prot)):
        gmm = GaussianMixture(n_components=K, random_state=SEED + fold,
                              covariance_type="full", max_iter=300, n_init=3)
        gmm.fit(X[tr])
        lab_tr = gmm.predict(X[tr])
        prof_tr = pd.Series(y[tr]).groupby(lab_tr).median()
        # consensus regimes by training profile
        cons = set(prof_tr.index[(prof_tr > 0.70) | (prof_tr < 0.30)])
        lab_te = gmm.predict(X[te])
        g = pd.DataFrame({"prot": prot[te], "reg": lab_te, "y": y[te]})
        for uni, gg in g.groupby("prot"):
            if len(gg) < 100:
                continue
            prof_te = gg.groupby("reg")["y"].median()
            common = prof_tr.index.intersection(prof_te.index)
            if len(common) < 4:
                continue
            rho_all = spearmanr(prof_tr[common], prof_te[common])[0]
            nc = [r for r in common if r not in cons]
            rho_nc = (spearmanr(prof_tr[nc], prof_te[nc])[0]
                      if len(nc) >= 3 else np.nan)
            recs.append({"UniProt_ID": uni, "fold": fold,
                         "rho_all": rho_all, "rho_noconsensus": rho_nc,
                         "n_noncons": len(nc),
                         "n": len(gg)})
            for r_ in common:
                tag = "consensus" if r_ in cons else "nonconsensus"
                recs[-1][f"trY_r{r_}"] = prof_tr[r_]
        # per-regime pooled held-out Y + train profile
        for r_ in range(K):
            sub = g[g.reg == r_]
            if len(sub) >= 50:
                recs.append({"UniProt_ID": f"__REGIME__{r_}", "fold": fold,
                             "reg": r_, "heldout_n": len(sub),
                             "heldout_Y": float(sub["y"].median()),
                             "train_Y": float(prof_tr.get(r_, np.nan)),
                             "rho_all": np.nan, "rho_noconsensus": np.nan,
                             "n": len(sub)})
    rep = pd.DataFrame(recs)

    per_prot = rep[~rep["UniProt_ID"].astype(str).str.startswith("__REGIME__")]
    med_all = per_prot["rho_all"].median()
    nc_vals = per_prot["rho_noconsensus"].dropna()
    med_nc = nc_vals.median()
    rng = np.random.default_rng(SEED)
    boot_nc = np.median(rng.choice(nc_vals, (2000, len(nc_vals)),
                                   replace=True), axis=1)

    reg = rep[rep["UniProt_ID"].astype(str).str.startswith("__REGIME__")]
    reg_summary = reg.groupby("reg").agg(
        heldout_n=("heldout_n", "sum"),
        heldout_Y=("heldout_Y", "median"),
        train_Y=("train_Y", "median")).reset_index()

    summary = pd.DataFrame([{
        "metric": "median held-out regime-profile Spearman (all regimes)",
        "value": float(med_all), "ci_low": float(np.nanpercentile(
            per_prot["rho_all"], 25)),
        "ci_high": float(np.nanpercentile(per_prot["rho_all"], 75)),
        "n": len(per_prot)},
        {"metric": "median held-out Spearman excluding the two consensus regimes",
         "value": float(med_nc),
         "ci_low": float(np.percentile(nc_vals, 25)),
         "ci_high": float(np.percentile(nc_vals, 75)),
         "n": int(len(nc_vals))},
        {"metric": "fraction of proteins positive (no-consensus)",
         "value": float((nc_vals > 0).mean()), "ci_low": np.nan,
         "ci_high": np.nan, "n": int(len(nc_vals))}])
    summary.to_csv(STATISTICS / "heldout_noconsensus_summary.csv", index=False)
    reg_summary.to_csv(STATISTICS / "heldout_noconsensus_regimes.csv",
                       index=False)
    print(summary.to_string(index=False))
    print(reg_summary.to_string(index=False))
    job.close()


if __name__ == "__main__":
    main()
