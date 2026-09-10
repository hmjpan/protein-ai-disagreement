"""23_residue_level_plddt.py

Review-response analysis (item 2): replace the quartile-mean (4-point)
Spearman as the PRIMARY statistical evidence with protein-wise residue-level
correlation.

For each protein:
  * aggregate D_std / D_evo_struct to the residue level (median over
    mutations at the same position; pLDDT is per-position)
  * compute Spearman(pLDDT, disagreement) directly over residues
  * also compute Kendall tau as a second association measure

Aggregation: median across proteins; protein-level bootstrap 95% CI;
Wilcoxon signed-rank test. The quartile visualization is retained but
demoted to a figure.

Outputs:
  results/statistics/residue_level_plddt_effects.csv
  results/statistics/residue_level_plddt_summary.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job  # noqa: E402

SEED = 2026
MIN_RESIDUES = 50


def main():
    job = Job("23_residue_level_plddt")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    df = disc.merge(struct, on=["DMS_id", "mutant"], how="left")
    job.info(f"rows: {len(df)}")

    rows = []
    for uni, g in df.groupby("UniProt_ID"):
        # aggregate to residue level
        rl = g.groupby("position").agg(
            D_std=("D_std", "median"),
            D_evo_struct=("D_evo_struct", "median"),
            D_seq_struct=("D_seq_struct", "median"),
            plddt=("plddt", "median"),
            n_mut=("D_std", "size")).dropna()
        if len(rl) < MIN_RESIDUES:
            continue
        rho_D = spearmanr(rl["plddt"], rl["D_std"])[0]
        rho_De = spearmanr(rl["plddt"], rl["D_evo_struct"])[0]
        rho_Ds = spearmanr(rl["plddt"], rl["D_seq_struct"])[0]
        tau_De = kendalltau(rl["plddt"], rl["D_evo_struct"])[0]
        rows.append({"UniProt_ID": uni, "n_residues": len(rl),
                     "n_mutations": int(rl.n_mut.sum()),
                     "rho_D_plddt": rho_D, "rho_De_plddt": rho_De,
                     "rho_Ds_plddt": rho_Ds, "tau_De_plddt": tau_De})
    eff = pd.DataFrame(rows)
    eff.to_csv(STATISTICS / "residue_level_plddt_effects.csv", index=False)
    job.info(f"proteins analysed: {len(eff)}")

    rng = np.random.default_rng(SEED)
    summ = []
    for col in ["rho_D_plddt", "rho_De_plddt", "rho_Ds_plddt", "tau_De_plddt"]:
        v = eff[col].dropna().to_numpy()
        boots = np.array([np.median(rng.choice(v, size=len(v), replace=True))
                          for _ in range(2000)])
        w, p = wilcoxon(v)  # one-sample: H0 median = 0
        summ.append({"metric": col, "median": np.median(v),
                     "ci_low": np.percentile(boots, 2.5),
                     "ci_high": np.percentile(boots, 97.5),
                     "frac_negative": (v < 0).mean(),
                     "wilcoxon_p": p})
    sdf = pd.DataFrame(summ)
    sdf.to_csv(STATISTICS / "residue_level_plddt_summary.csv", index=False)
    job.info("\n" + sdf.round(4).to_string(index=False))

    # effect size in disagreement units: low vs high pLDDT tertile
    eff2 = []
    for uni, g in df.groupby("UniProt_ID"):
        rl = g.groupby("position").agg(
            D_evo_struct=("D_evo_struct", "median"),
            D_std=("D_std", "median"),
            plddt=("plddt", "median")).dropna()
        if len(rl) < MIN_RESIDUES:
            continue
        t = rl["plddt"].rank(pct=True)
        lo = rl[t <= 0.34]
        hi = rl[t >= 0.66]
        if len(lo) < 10 or len(hi) < 10:
            continue
        eff2.append({"UniProt_ID": uni,
                     "dD_evo_struct": lo.D_evo_struct.median() - hi.D_evo_struct.median(),
                     "dD_std": lo.D_std.median() - hi.D_std.median()})
    e2 = pd.DataFrame(eff2)
    e2.to_csv(STATISTICS / "residue_level_plddt_tertile_effects.csv", index=False)
    for col in ["dD_evo_struct", "dD_std"]:
        v = e2[col].to_numpy()
        boots = np.array([np.median(rng.choice(v, size=len(v), replace=True))
                          for _ in range(2000)])
        job.info(f"{col}: median low-tertile-minus-high = {np.median(v):+.4f} "
                 f"[{np.percentile(boots,2.5):+.4f}, {np.percentile(boots,97.5):+.4f}], "
                 f"frac positive {(v > 0).mean():.3f}")
    job.close()


if __name__ == "__main__":
    main()