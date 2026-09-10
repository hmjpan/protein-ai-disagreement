"""24_disagreement_definition_robustness.py

Review-response analysis (item 11): show the pLDDT organization finding does
not depend on the disagreement definition.

Disagreement definitions compared (all at residue level, per protein):
  1. z-std      : SD of z-scores across core models (primary, D_std)
  2. u-std      : SD of percentile scores
  3. D_MAD      : median absolute deviation of z-scores
  4. mean-pair  : mean pairwise absolute z difference
  5. Kendall    : Kendall's W (concordance) transformed to disagreement 1-W

Outputs:
  results/statistics/disagreement_definition_plddt.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job, load_core_panel  # noqa: E402

MIN_RESIDUES = 50
SEED = 2026


def main():
    job = Job("24_disagreement_definition_robustness")
    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    core = load_core_panel(job)
    z_cols = [f"z_{m}" for m in core["model"]]
    u_cols = [f"u_{m}" for m in core["model"]]

    df = disc.merge(struct, on=["DMS_id", "mutant"], how="left")
    df = df.merge(norm[["DMS_id", "mutant"] + u_cols], on=["DMS_id", "mutant"],
                  how="left")
    Z = norm[z_cols].to_numpy(dtype=float)
    U = df[u_cols].to_numpy(dtype=float)
    # per-variant definitions (z-basis except u-std)
    D_z = np.nanstd(Z, axis=1)
    D_u = np.nanstd(U, axis=1)
    D_mad = df["D_mad"].to_numpy()
    with np.errstate(invalid="ignore"):
        n = Z.shape[1]
        # mean pairwise absolute difference over the M(M-1)/2 unordered pairs,
        # chunked over variants (constant 2 not needed for Spearman, but the
        # true mean includes it)
        D_pair = np.full(len(Z), np.nan)
        chunk = 50000
        npairs = n * (n - 1) / 2
        for i in range(0, len(Z), chunk):
            Zc = Z[i:i + chunk]
            d3 = np.abs(Zc[:, :, None] - Zc[:, None, :])  # (c,40,40)
            idx = np.arange(d3.shape[1])
            d3[:, idx, idx] = 0
            D_pair[i:i + chunk] = d3.sum(axis=(1, 2)) / (2 * npairs)
        D_pair[~np.isfinite(D_pair)] = np.nan

    # Four independent disagreement definitions (review: D_rankvar was
    # mathematically redundant with D_u (var(u) = SD(u)^2) and removed)
    defs = {"D_std": D_z, "D_u": D_u, "D_MAD": D_mad, "D_pair": D_pair}
    df = df.assign(**{k: v for k, v in defs.items()})

    rows = []
    for uni, g in df.groupby("UniProt_ID"):
        rl = g.groupby("position").agg(
            **{k: (k, "median") for k in defs},
            plddt=("plddt", "median")).dropna()
        if len(rl) < MIN_RESIDUES:
            continue
        for k in defs:
            rows.append({"UniProt_ID": uni, "definition": k, "n": len(rl),
                         "rho": spearmanr(rl["plddt"], rl[k])[0]})
    out = pd.DataFrame(rows)
    out.to_csv(STATISTICS / "disagreement_definition_plddt.csv", index=False)

    rng = np.random.default_rng(SEED)
    summ = []
    for k, g in out.groupby("definition"):
        v = g["rho"].dropna().to_numpy()
        boots = np.array([np.median(rng.choice(v, size=len(v), replace=True))
                          for _ in range(2000)])
        summ.append({"definition": k, "n_proteins": len(v),
                     "median_rho": np.median(v),
                     "ci_low": np.percentile(boots, 2.5),
                     "ci_high": np.percentile(boots, 97.5),
                     "frac_negative": (v < 0).mean()})
    s = pd.DataFrame(summ)
    s.to_csv(STATISTICS / "disagreement_definition_plddt_summary.csv", index=False)
    job.info("\n" + s.round(4).to_string(index=False))
    job.close()


if __name__ == "__main__":
    main()