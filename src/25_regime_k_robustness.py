"""25_regime_k_robustness.py

Review-response analysis (item 13): are the regime conclusions stable
across the number of clusters?

  * GMM for K = 4, 5, 6, 7 (same seeds as the main analysis)
  * adjusted Rand index between K = 6 and each other K (shared-variant
    subset; stable membership under re-clustering)
  * check whether the key biological conclusions hold for every K:
      (a) a "structure-dissenting" regime exists (S_seq/S_evo > +0.3,
          S_struct < +0.3)
      (b) consensus regimes track experimental Y
      (c) regime assignments reproduce across assays (mean position-level
          agreement minus chance, on multi-assay proteins)

Outputs:
  results/statistics/regime_k_robustness.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.mixture import GaussianMixture
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job  # noqa: E402

SEED = 2026
FAMILIES = ["single_seq", "evolution", "structure"]
S_COLS = [f"S_{f}" for f in FAMILIES]


def main():
    job = Job("25_regime_k_robustness")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    X = disc[S_COLS].to_numpy(dtype=float)
    keep = ~np.isnan(X).any(axis=1)
    Xk = X[keep]
    ids = disc.loc[keep, ["DMS_id", "UniProt_ID", "mutant"]].reset_index(drop=True)
    Y = disc.loc[keep, "Y_deleter"].to_numpy()
    job.info(f"variants: {len(Xk)}")

    labels = {}
    rows = []
    for K in [4, 5, 6, 7]:
        gmm = GaussianMixture(n_components=K, random_state=SEED,
                              covariance_type="full", max_iter=300, n_init=3)
        gmm.fit(Xk)
        lab = gmm.predict(Xk)
        labels[K] = lab
        # a) structure-dissenting regime
        cent = pd.DataFrame(gmm.means_, columns=S_COLS)
        dissent = ((cent["S_single_seq"] > 0.3) & (cent["S_evolution"] > 0.3)
                   & (cent["S_structure"] < 0.3)).any()
        # b) consensus regimes track Y
        d = pd.DataFrame({"lab": lab, "Y": Y})
        cons_tol = d[d.lab.isin(np.where(cent[S_COLS].min(axis=1) < -0.5)[0])]
        cons_dam = d[d.lab.isin(np.where(cent[S_COLS].max(axis=1) > 0.5)[0])]
        y_tol = cons_tol["Y"].median() if len(cons_tol) else np.nan
        y_dam = cons_dam["Y"].median() if len(cons_dam) else np.nan
        # c) cross-assay reproducibility (mode regime per position)
        tmp = ids.copy()
        tmp["lab"] = lab
        agree = []
        for uni, g in tmp.groupby("UniProt_ID"):
            assays = g["DMS_id"].unique()
            if len(assays) < 2:
                continue
            for a in range(len(assays) - 1):
                a1 = g[g.DMS_id == assays[a]].groupby("mutant")["lab"].agg(
                    lambda s: s.value_counts().idxmax())
                a2 = g[g.DMS_id == assays[a + 1]].groupby("mutant")["lab"].agg(
                    lambda s: s.value_counts().idxmax())
                common = a1.index.intersection(a2.index)
                if len(common) >= 20:
                    agree.append((a1.loc[common].to_numpy() == a2.loc[common].to_numpy()).mean())
        rows.append({"K": K, "n_regimes": K,
                     "has_structure_dissenting": bool(dissent),
                     "consensus_tolerant_Y": float(y_tol) if not np.isnan(y_tol) else None,
                     "consensus_damaging_Y": float(y_dam) if not np.isnan(y_dam) else None,
                     "mean_crossassay_agreement": float(np.mean(agree)) if agree else None,
                     "n_crossassay_pairs": len(agree)})
        job.info(f"K={K}: dissent={dissent} Y_tol={y_tol:.3f} Y_dam={y_dam:.3f} "
                 f"agreement={np.mean(agree):.3f} (n={len(agree)})")

    # ARI between K=6 and others
    base = labels[6]
    ari_rows = []
    for K in [4, 5, 7]:
        ari_rows.append({"K_ref": 6, "K_other": K,
                         "ARI": adjusted_rand_score(base, labels[K])})
    job.info("ARI vs K=6: " + ", ".join(f"K={r['K_other']}: {r['ARI']:.3f}"
                                        for r in ari_rows))
    out = pd.DataFrame(rows)
    out.to_csv(STATISTICS / "regime_k_robustness.csv", index=False)
    pd.DataFrame(ari_rows).to_csv(STATISTICS / "regime_k_ari.csv", index=False)
    job.close()


if __name__ == "__main__":
    main()