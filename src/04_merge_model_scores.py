"""04_merge_model_scores.py

Merge the zero-shot model prediction scores onto the primary single-mutation
set (is_single & ~sequence_mismatch). Produces one wide parquet:
  assay_id, UniProt_ID, mutant, wt_aa, position, mut_aa, DMS_score,
  DMS_score_bin, <model_1> ... <model_N>

Also writes results/tables/TableS3b_variant_level_model_coverage.csv with the
fraction of primary variants scored by each model per assay.

Inputs : data/processed/processed_single_mutations.parquet
         data/raw/proteingym_scores/zero_shot_substitutions_scores/*.csv
Outputs: data/intermediate/merged_scores_wide.parquet
         results/tables/TableS3b_variant_level_model_coverage.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import INTERMEDIATE, PROCESSED, SCORES_DIR, TABLES, Job  # noqa: E402


def main():
    job = Job("04_merge_model_scores")
    prim = pd.read_parquet(PROCESSED / "processed_single_mutations.parquet")
    job.info(f"primary set: {len(prim)} rows, {prim.DMS_id.nunique()} assays")

    score_files = {p.stem: p for p in SCORES_DIR.glob("*.csv")}
    present = prim["DMS_id"].isin(score_files)
    job.info(f"assays with score files: {present.sum()} / {prim.DMS_id.nunique()}")
    prim = prim[present].copy()

    frames = []
    coverage_rows = []
    model_union = set()
    for dms_id, group in prim.groupby("DMS_id"):
        p = score_files[dms_id]
        sc = pd.read_csv(p)
        sc["mutant"] = sc["mutant"].astype(str)
        merged = group.merge(sc, on="mutant", how="left", suffixes=("", "_dup"))
        frames.append(merged)
        n = len(group)
        for c in sc.columns:
            if c == "mutant":
                continue
            model_union.add(c)
            coverage_rows.append({"DMS_id": dms_id, "model": c,
                                  "n_primary": n,
                                  "n_scored": int(merged[c].notna().sum())})
        job.info(f"{dms_id}: merged {n} variants, {len(sc.columns)-1} model columns")

    merged_all = pd.concat(frames, ignore_index=True)
    merged_all.to_parquet(INTERMEDIATE / "merged_scores_wide.parquet", index=False)
    job.info(f"merged_scores_wide.parquet: {len(merged_all)} rows x {merged_all.shape[1]} cols")

    cov = pd.DataFrame(coverage_rows)
    cov["frac_scored"] = cov["n_scored"] / cov["n_primary"]
    cov.to_csv(TABLES / "TableS3b_variant_level_model_coverage.csv", index=False)
    job.info(f"coverage table: {len(cov)} rows; {len(model_union)} distinct models")
    job.close()


if __name__ == "__main__":
    main()