"""36_aligner_mapping.py

Reviewer-response (item 11-2): replace difflib SequenceMatcher UniProt->DMS
coordinate mapping with a standard Needleman-Wunsch global alignment
(Biopython PairwiseAligner, C implementation), and report mapping quality.

Outputs:
  data/processed/uniprot_position_map_pa.parquet   (new mapping)
  results/statistics/mapping_quality.csv           (identity/coverage/gap rate)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, RAW, STATISTICS, Job, read_reference_dms  # noqa: E402

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2.0
aligner.mismatch_score = -1.0
aligner.open_gap_score = -5.0
aligner.extend_gap_score = -0.5


def main():
    job = Job("36_aligner_mapping")
    ref = read_reference_dms()
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    # cache fasta from 10's fetch
    fasta_cache = RAW / "uniprot" / "fasta"
    rows, qual = [], []
    for uni, g in ref.groupby("UniProt_ID"):
        acc = uni.split("_")[0]
        f = fasta_cache / f"{uni}.fasta"
        if not f.exists():
            f = fasta_cache / f"{acc}.fasta"
        if not f.exists():
            # resolve from cached JSON
            jp = RAW / "uniprot" / f"{uni}.json"
            if jp.exists():
                import json
                try:
                    pa = json.loads(jp.read_text(encoding="utf-8")).get("primaryAccession", acc)
                except Exception:
                    pa = acc
                f = fasta_cache / f"{pa}.fasta"
                if not f.exists():
                    job.info(f"{uni}: fasta missing")
                    continue
            else:
                continue
        uni_seq = "".join(f.read_text().splitlines()[1:])
        target = str(g["target_seq"].iloc[0])
        if not uni_seq or not target:
            job.info(f"{uni}: empty sequence skipped")
            continue
        aln = aligner.align(target, uni_seq)[0]
        tstr, ustr = str(aln[0]), str(aln[1])
        ti = ui = 0
        mapped = 0
        identity_pos = 0
        for tc, uc in zip(tstr, ustr):
            if tc != "-":
                ti += 1
            if uc != "-":
                ui += 1
            if tc != "-" and uc != "-":
                rows.append({"UniProt_ID": uni, "target_pos": ti,
                             "uniprot_pos": ui, "mapped": True,
                             "target_len": len(target), "uni_len": len(uni_seq)})
                mapped += 1
                if tc == uc:
                    identity_pos += 1
        qual.append({"UniProt_ID": uni, "target_len": len(target),
                     "uni_len": len(uni_seq), "n_mapped": mapped,
                     "coverage": mapped / max(len(target), 1),
                     "identity": identity_pos / max(mapped, 1)})
    pm = pd.DataFrame(rows)
    pm.to_parquet(PROCESSED / "uniprot_position_map_pa.parquet", index=False)
    q = pd.DataFrame(qual)
    q.to_csv(STATISTICS / "mapping_quality.csv", index=False)
    job.info(f"proteins aligned: {len(q)}; mean coverage {q.coverage.mean():.3f}; "
             f"mean identity {q.identity.mean():.3f}")

    # agreement with difflib mapping
    old = pd.read_parquet(PROCESSED / "uniprot_position_map.parquet")
    o = old[old["mapped"]][["UniProt_ID", "target_pos", "uniprot_pos"]].rename(
        columns={"uniprot_pos": "old_u"})
    n = pm.rename(columns={"uniprot_pos": "new_u"})
    cmp_ = o.merge(n, on=["UniProt_ID", "target_pos"], how="inner")
    agree = (cmp_["old_u"] == cmp_["new_u"]).mean()
    job.info(f"position-mapping agreement with difflib: {agree:.4f} "
             f"(n={len(cmp_)} shared positions)")
    job.close()


if __name__ == "__main__":
    main()
