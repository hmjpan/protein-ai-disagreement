"""29_regime_phenotype_replication.py

Review-response analysis (item 2): replace the 'regime assignments agree
across assays' claim with true EXPERIMENTAL replication -- does the
relationship between regime and the measured DMS phenotype replicate across
independent assays of the same protein?

For each multi-assay protein:
  * assign regime labels from model predictions (K = 6, as in the main
    analysis) on the merged variant set
  * per assay, compute the median experimental Y within each regime
  * correlate regime-level Y across assay pairs (Pearson/Spearman over
    regimes with >= 20 variants in both assays)
  * additionally: consensus-regime direction consistency
      - consensus-damaging Y > 0.5 in both assays?
      - consensus-tolerant Y < 0.5 in both assays?
  * structure-dissenting regime: does DMS Y exceed the structure-family
    prediction-implied expectation in both assays (i.e., phenotype aligns
    with seq/evo)?

Outputs:
  results/statistics/regime_phenotype_replication.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job  # noqa: E402

SEED = 2026
K = 6
MIN_SHARED = 20
S_COLS = ["S_single_seq", "S_evolution", "S_structure"]


def main():
    job = Job("29_regime_phenotype_replication")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    X = disc[S_COLS].to_numpy(dtype=float)
    keep = ~np.isnan(X).any(axis=1)
    Xk = X[keep]
    gmm = GaussianMixture(n_components=K, random_state=SEED,
                          covariance_type="full", max_iter=300, n_init=3)
    gmm.fit(Xk)
    lab = gmm.predict(Xk)
    ids = disc.loc[keep, ["DMS_id", "UniProt_ID", "mutant", "Y_deleter"]].copy()
    ids["regime"] = lab
    job.info(f"variants assigned: {len(ids)}")

    rows = []
    for uni, g in ids.groupby("UniProt_ID"):
        assays = sorted(g["DMS_id"].unique())
        if len(assays) < 2:
            continue
        for a in range(len(assays) - 1):
            for b in range(a + 1, len(assays)):
                ga = g[g.DMS_id == assays[a]]
                gb = g[g.DMS_id == assays[b]]
                ya = ga.groupby("regime")["Y_deleter"].median()
                yb = gb.groupby("regime")["Y_deleter"].median()
                na = ga.groupby("regime").size()
                nb = gb.groupby("regime").size()
                common = [r for r in ya.index
                          if r in yb.index and na[r] >= MIN_SHARED and nb[r] >= MIN_SHARED]
                if len(common) < 3:
                    continue
                rho = spearmanr(ya[common].to_numpy(), yb[common].to_numpy())[0]
                # consensus-regime direction consistency
                cent = pd.DataFrame(gmm.means_, columns=S_COLS)
                dam_regs = cent[cent[S_COLS].min(axis=1) > 0.5].index.tolist()
                tol_regs = cent[cent[S_COLS].max(axis=1) < -0.5].index.tolist()
                diss_regs = cent[(cent.S_single_seq > 0.3) & (cent.S_evolution > 0.3)
                                 & (cent.S_structure < 0.3)].index.tolist()
                dam_ok = all(ya.get(r, np.nan) > 0.5 and yb.get(r, np.nan) > 0.5
                             for r in dam_regs if r in ya and r in yb)
                tol_ok = all(ya.get(r, np.nan) < 0.5 and yb.get(r, np.nan) < 0.5
                             for r in tol_regs if r in ya and r in yb)
                diss_y = [(ya.get(r, np.nan), yb.get(r, np.nan)) for r in diss_regs
                          if r in ya and r in yb]
                diss_ok = all(not np.isnan(x) and not np.isnan(y_) and x > 0.5 and y_ > 0.5
                              for x, y_ in diss_y) if diss_y else None
                rows.append({"UniProt_ID": uni, "assay_a": assays[a], "assay_b": assays[b],
                             "n_regimes_shared": len(common),
                             "rho_regime_Y": rho,
                             "consensus_damaging_both": dam_ok,
                             "consensus_tolerant_both": tol_ok,
                             "structure_dissenting_damaging_both": diss_ok})
    rep = pd.DataFrame(rows)
    rep.to_csv(STATISTICS / "regime_phenotype_replication.csv", index=False)
    if len(rep):
        job.info(f"assay pairs with >=3 shared regimes: {len(rep)}")
        job.info(f"median rho(regime-level Y across assays) = "
                 f"{rep.rho_regime_Y.median():+.3f}")
        job.info(f"consensus-damaging consistent (both assays Y>0.5): "
                 f"{rep.consensus_damaging_both.mean():.3f}")
        job.info(f"consensus-tolerant consistent (both assays Y<0.5): "
                 f"{rep.consensus_tolerant_both.mean():.3f}")
        d = rep["structure_dissenting_damaging_both"].dropna()
        if len(d):
            job.info(f"structure-dissenting damaging in both assays (where present): "
                     f"{d.mean():.3f} ({len(d)} pairs)")
    else:
        job.info("no analysable pairs")
    job.close()


if __name__ == "__main__":
    main()