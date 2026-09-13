"""26_effect_size_plddt.py

Review-response analysis (item 1): provide intuitive effect sizes for the
residue-level pLDDT-disagreement association, beyond Spearman:

  * within-protein comparison of the lowest vs highest pLDDT decile:
    median D_evo_struct difference, and Cohen's d (standardized)
  * enrichment: in low-pLDDT decile residues, the odds of falling in the
    protein-internal top disagreement decile vs all other residues
  * pooled enrichment OR with 95% CI (protein-level bootstrap)

Outputs:
  results/statistics/plddt_effect_sizes.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job  # noqa: E402

SEED = 2026
MIN_RESIDUES = 50


def main():
    job = Job("26_effect_size_plddt")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    df = disc.merge(struct, on=["DMS_id", "mutant"], how="left")
    job.info(f"rows: {len(df)}")

    rows = []
    for uni, g in df.groupby("UniProt_ID"):
        rl = g.groupby("position").agg(
            D_evo_struct=("D_evo_struct", "median"),
            D_std=("D_std", "median"),
            plddt=("plddt", "median")).dropna()
        if len(rl) < MIN_RESIDUES:
            continue
        pct = rl["plddt"].rank(pct=True)
        lo = rl[pct <= 0.10]
        hi = rl[pct >= 0.90]
        if len(lo) < 5 or len(hi) < 5:
            continue
        dD_evo = lo.D_evo_struct.median() - hi.D_evo_struct.median()
        dD_tot = lo.D_std.median() - hi.D_std.median()
        pooled_sd = np.sqrt((lo.D_evo_struct.std() ** 2 + hi.D_evo_struct.std() ** 2) / 2)
        cohens_d = dD_evo / pooled_sd if pooled_sd > 0 else np.nan
        # enrichment: low-pLDDT decile vs top-disagreement decile
        topD = rl["D_evo_struct"].rank(pct=True) >= 0.90
        a = int(((pct <= 0.10) & topD).sum())     # low pLDDT & high D
        b = int(((pct <= 0.10) & ~topD).sum())
        c = int(((pct > 0.10) & topD).sum())
        dcell = int(((pct > 0.10) & ~topD).sum())
        if a > 0 and b > 0 and c > 0 and dcell > 0:
            or_ = (a * dcell) / (b * c)
        else:
            or_ = np.nan
        rows.append({"UniProt_ID": uni, "n_residues": len(rl),
                     "dD_evo_struct_decile": dD_evo, "dD_std_decile": dD_tot,
                     "cohens_d": cohens_d, "enrich_OR": or_,
                     "n_low_hi": a, "n_low_other": b,
                     "n_other_hi": c, "n_other_other": dcell})
    eff = pd.DataFrame(rows)
    eff.to_csv(STATISTICS / "plddt_effect_sizes.csv", index=False)
    job.info(f"proteins analysed: {len(eff)}")

    rng = np.random.default_rng(SEED)
    for col in ["dD_evo_struct_decile", "dD_std_decile", "cohens_d", "enrich_OR"]:
        v = eff[col].dropna().to_numpy()
        boots = np.array([np.median(rng.choice(v, size=len(v), replace=True))
                          for _ in range(2000)])
        job.info(f"{col}: median = {np.median(v):+.4f} "
                 f"[{np.percentile(boots,2.5):+.4f}, {np.percentile(boots,97.5):+.4f}], "
                 f"n_proteins = {len(v)}")
    # pooled enrichment OR (marginal)
    a, b, c, dcell = (eff["n_low_hi"].sum(), eff["n_low_other"].sum(),
                      eff["n_other_hi"].sum(), eff["n_other_other"].sum())
    or_pooled = (a * dcell) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / dcell)
    job.info(f"pooled enrichment OR = {or_pooled:.2f} "
             f"(95% CI [{np.exp(np.log(or_pooled)-1.96*se):.2f}, "
             f"{np.exp(np.log(or_pooled)+1.96*se):.2f}])")
    job.close()


if __name__ == "__main__":
    main()