"""43_struct_subset_sensitivity.py

Reviewer-response: the structure family has only three released models
(ESM-IF1, MIF, ProteinMPNN), so all structure-family analyses fix that
panel. Sensitivity: recompute the residue-level pLDDT association for
every structure-subset configuration.

Configurations (structure centroid = mean percentile-space u over the
subset; evolution centroid fixed at the full 12-model family centroid):
  each single model; each pair; leave-one-out (= pairs); full trio.

Statistic (primary, Section 2.2): within-protein Spearman(pLDDT,
|U_evo - U_struct_subset|) aggregated by median over proteins with
>= 50 mapped residues; fraction negative; binomial sign p.

Output: results/statistics/struct_subset_sensitivity.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job  # noqa: E402

STRUCT = ["ESM-IF1", "MIF", "ProteinMPNN"]
EVO = None  # filled from family centroid column
MIN_RES = 50


def main():
    job = Job("43_struct_subset_sensitivity")
    norm = pd.read_parquet(
        PROCESSED / "normalized_scores.parquet",
        columns=["DMS_id", "mutant"] + [f"u_{m}" for m in STRUCT])
    disc = pd.read_parquet(
        PROCESSED / "disagreement_scores.parquet",
        columns=["DMS_id", "UniProt_ID", "mutant", "position", "S_evolution"])
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    d = disc.merge(struct, on=["DMS_id", "mutant"]).merge(
        norm, on=["DMS_id", "mutant"])
    d = d.dropna(subset=["plddt"])
    job.info(f"rows: {len(d)}")

    # S_evolution is z-space family centroid; convert to percentile within
    # assay to match the U (percentile) scale of the subset centroids
    d["U_evo"] = d.groupby("DMS_id")["S_evolution"].rank(pct=True)

    def run(cols, name):
        d["U_sub"] = d[cols].mean(axis=1)
        d["D_sub"] = (d["U_evo"] - d["U_sub"]).abs()
        rhos = []
        for _, g in d.groupby(["UniProt_ID"]):
            rl = g.groupby("position").agg(plddt=("plddt", "median"),
                                           D=("D_sub", "median")).dropna()
            if len(rl) >= MIN_RES:
                r = spearmanr(rl["plddt"], rl["D"])[0]
                if np.isfinite(r):
                    rhos.append(r)
        rhos = np.array(rhos)
        med = float(np.median(rhos))
        boot = np.median(rhos[np.random.default_rng(2026).integers(
            0, len(rhos), (2000, len(rhos)))], axis=1)
        return {"config": name, "n_struct_models": len(cols) if len(cols) > 1
                else 1, "n_proteins": int(len(rhos)),
                "median_rho": med,
                "ci_low": float(np.percentile(boot, 2.5)),
                "ci_high": float(np.percentile(boot, 97.5)),
                "frac_negative": float(np.mean(rhos < 0)),
                "sign_p": float(binomtest(int((rhos < 0).sum()), len(rhos),
                                          0.5).pvalue)}

    rows = []
    u = [f"u_{m}" for m in STRUCT]
    rows.append(run(u, "full trio (ESM-IF1 + MIF + ProteinMPNN)"))
    for m in STRUCT:
        rows.append(run([f"u_{m}"], f"single: {m}"))
    for i in range(3):
        sub = [u[j] for j in range(3) if j != i]
        rows.append(run(sub, f"pair: {' + '.join(STRUCT[j] for j in range(3) if j != i)}"))
    out = pd.DataFrame(rows)
    out.to_csv(STATISTICS / "struct_subset_sensitivity.csv", index=False)
    print(out.to_string(index=False))
    job.close()


if __name__ == "__main__":
    main()
