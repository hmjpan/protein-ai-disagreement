"""12b_regime_reproducibility.py

Reviewer-defense analysis: are disagreement regimes reproducible rather than
assay-specific noise?

For every protein with >= 2 DMS assays in ProteinGym, compare the regime
assignment of mutations at SHARED positions across assay pairs (regime labels
come from 12_disagreement_regimes.py). Report:
  * mean position-level agreement within protein (identical regime label)
  * chance expectation (proportional agreement from the marginal regime
    distribution per assay pair)
  * fraction of proteins with agreement > chance

Outputs: results/statistics/regime_reproducibility.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job  # noqa: E402


def main():
    job = Job("12b_regime_reproducibility")
    assign = pd.read_parquet(STATISTICS / "regime_assignment.parquet")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet",
                           columns=["DMS_id", "UniProt_ID", "mutant", "position"])
    df = assign.merge(disc[["DMS_id", "UniProt_ID", "mutant", "position"]],
                      on=["DMS_id", "mutant"], how="left")
    job.info(f"assigned variants: {len(df)}; proteins: {df.UniProt_ID.nunique()}")

    rows = []
    n_multi = 0
    for uni, g in df.groupby("UniProt_ID"):
        assays = g["DMS_id"].unique()
        if len(assays) < 2:
            continue
        n_multi += 1
        # pair-wise agreement at shared positions
        pairs = list(zip(assays, assays[1:])) + \
            ([(assays[0], a) for a in assays[2:]] if len(assays) > 2 else [])
        for a1, a2 in pairs:
            g1 = g[g.DMS_id == a1].groupby("position")["regime_name"].agg(
                lambda s: s.value_counts().idxmax())
            g2 = g[g.DMS_id == a2].groupby("position")["regime_name"].agg(
                lambda s: s.value_counts().idxmax())
            common = g1.index.intersection(g2.index)
            if len(common) < 20:
                continue
            lab1 = g1.loc[common].to_numpy()
            lab2 = g2.loc[common].to_numpy()
            obs = (lab1 == lab2).mean()
            # chance: product of marginals summed
            p1 = pd.Series(lab1).value_counts(normalize=True)
            p2 = pd.Series(lab2).value_counts(normalize=True)
            chance = sum(p1.get(k, 0) * p2.get(k, 0) for k in set(p1.index) | set(p2.index))
            rows.append({"UniProt_ID": uni, "assay_a": a1, "assay_b": a2,
                         "n_shared_positions": len(common),
                         "agreement": obs, "chance_agreement": chance,
                         "excess": obs - chance})
    rep = pd.DataFrame(rows)
    rep.to_csv(STATISTICS / "regime_reproducibility.csv", index=False)
    if len(rep):
        job.info(f"protein-assay pairs with >=20 shared positions: {len(rep)} "
                 f"(from {n_multi} multi-assay proteins)")
        job.info(f"mean position-level agreement: {rep.agreement.mean():.3f} "
                 f"(chance {rep.chance_agreement.mean():.3f}); "
                 f"excess {rep.excess.mean():.3f}")
        job.info(f"fraction of pairs with agreement > chance: "
                 f"{(rep.agreement > rep.chance_agreement).mean():.3f}")
    else:
        job.info("no analysable pairs")
    job.close()


if __name__ == "__main__":
    main()