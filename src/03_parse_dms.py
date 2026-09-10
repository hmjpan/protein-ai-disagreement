"""03_parse_dms.py

Parse every DMS assay file, keep ALL mutations for downstream sensitivity
analyses, and flag:
  * is_single        -- single amino-acid substitution (regex ^[A-Z]\\d+[A-Z]$)
  * sequence_mismatch-- target_seq[position-1] != wt_aa  (never silently fixed;
                        also true when position is out of range)
  * is_duplicate     -- duplicate mutant within the same assay

Primary analysis set (Phase 1) = is_single & ~sequence_mismatch.

Inputs : data/raw/proteingym_dms/DMS_ProteinGym_substitutions/*.csv
         data/raw/reference/DMS_substitutions.csv
Outputs: data/processed/processed_single_mutations.parquet  (primary set)
         data/intermediate/all_mutations.parquet           (full, flagged)
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DMS_DIR, INTERMEDIATE, PROCESSED, Job, read_reference_dms  # noqa: E402

SINGLE_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def main():
    job = Job("03_parse_dms")
    ref = read_reference_dms()
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    seq_map = dict(zip(ref["DMS_id"], ref["target_seq"].astype(str)))

    frames = []
    n_seq_mismatch = 0
    n_dup = 0
    for p in sorted(DMS_DIR.glob("*.csv")):
        dms_id = p.stem
        df = pd.read_csv(p)
        df["DMS_id"] = dms_id
        df["mutant"] = df["mutant"].astype(str)
        m = df["mutant"].str.extract(SINGLE_RE)
        df["wt_aa"] = m[0]
        df["position"] = pd.to_numeric(m[1], errors="coerce")
        df["mut_aa"] = m[2]
        df["is_single"] = df["mutant"].str.match(SINGLE_RE)

        target = seq_map.get(dms_id, "")
        tlist = list(target)
        n = len(df)
        mismatch = np.zeros(n, dtype=bool)
        single_idx = (df["is_single"] & df["position"].notna()).to_numpy()
        pos = df.loc[single_idx, "position"].astype(int).to_numpy()
        wt = df.loc[single_idx, "wt_aa"].to_numpy()
        out_of_range = (pos < 1) | (pos > len(tlist))
        mismatch[single_idx] = out_of_range
        in_range = single_idx.copy()
        in_range[in_range] = ~out_of_range
        if in_range.any():
            ok_pos = pos[~out_of_range] - 1
            mismatch[in_range] = np.array([tlist[i] for i in ok_pos]) != wt[~out_of_range]
        df["sequence_mismatch"] = mismatch

        df["is_duplicate"] = df.duplicated(subset=["mutant"], keep=False)
        n_seq_mismatch += int(df["sequence_mismatch"].sum())
        n_dup += int(df["is_duplicate"].sum())
        frames.append(df)
        job.info(f"{dms_id}: variants={len(df)} single={int(df.is_single.sum())} "
                 f"mismatch={int(df.sequence_mismatch.sum())} dup={int(df.is_duplicate.sum())}")

    all_df = pd.concat(frames, ignore_index=True)
    all_df.to_parquet(INTERMEDIATE / "all_mutations.parquet", index=False)
    job.info(f"all_mutations.parquet: {len(all_df)} rows, {all_df.DMS_id.nunique()} assays")

    primary = all_df[all_df["is_single"] & ~all_df["sequence_mismatch"]].copy()
    primary["UniProt_ID"] = primary["DMS_id"].map(
        dict(zip(ref["DMS_id"], ref["UniProt_ID"])))
    primary.to_parquet(PROCESSED / "processed_single_mutations.parquet", index=False)
    job.info(f"processed_single_mutations.parquet: {len(primary)} rows, "
             f"{primary.DMS_id.nunique()} assays, {primary.UniProt_ID.nunique()} proteins")
    job.info(f"totals: seq_mismatch={n_seq_mismatch} duplicate={n_dup}")
    job.close()


if __name__ == "__main__":
    main()