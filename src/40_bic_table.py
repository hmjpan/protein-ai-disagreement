"""40_bic_table.py

Review-response audit: expose the GMM model-selection BIC values for
K = 2..8 (identical settings to 12_disagreement_regimes.py: covariance
full, seed 2026, n_init 3) so the "BIC-optimal K = 6" claim is auditable.

Outputs:
  results/statistics/gmm_bic_table.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job  # noqa: E402

SEED = 2026
FAMILIES = ["single_seq", "evolution", "structure"]


def main():
    job = Job("40_bic_table")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    S_cols = [f"S_{f}" for f in FAMILIES]
    X = disc[S_cols].to_numpy(dtype=float)
    X = X[~np.isnan(X).any(axis=1)]
    job.info(f"fit rows: {len(X)}")
    out = pd.DataFrame({"K": [], "BIC": []})
    rows = []
    for k in range(2, 9):
        m = GaussianMixture(n_components=k, random_state=SEED,
                            covariance_type="full", max_iter=300, n_init=3)
        m.fit(X)
        rows.append({"K": k, "BIC": float(m.bic(X))})
    out = pd.DataFrame(rows)
    out["selected"] = out["K"] == int(out.loc[out["BIC"].idxmin(), "K"])
    out.to_csv(STATISTICS / "gmm_bic_table.csv", index=False)
    job.info("selected K = " + str(out.loc[out['selected'], 'K'].iloc[0]))



if __name__ == "__main__":
    main()
