"""09_structure_annotations.py

Extract per-residue AlphaFold pLDDT for every assay protein from the
ProteinGym-provided AlphaFold2 structures (PDB format; pLDDT stored in the
B-factor column) and attach it to every mutation.

Alignment protocol (pre-registered):
  1. Use the chain with the most CA atoms.
  2. Build a residue-number -> pLDDT map.
  3. Map DMS position (1-based in the assay target sequence) -> pLDDT by the
     SAME residue number.
  4. QC: report per-protein sequence alignment pass rate (fraction of target
     residues with a matching residue number in the PDB). If a protein has a
     PDB but < 50% of residues map, flag it; positions without pLDDT are left
     as NaN (structure-coverage sensitivity analysis handles the rest).
  5. pLDDT is a model confidence, NOT experimental disorder (wording rule).

Inputs : data/raw/structures/ProteinGym_AF2_structures/*.pdb
         data/processed/processed_single_mutations.parquet
Outputs: data/processed/structure_features.parquet
         results/statistics/structure_alignment_QC.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    PROCESSED, STRUCTURES_DIR, STATISTICS, Job, read_reference_dms,
)


def plddt_from_pdb(path: Path) -> dict:
    """Return {residue_number: pLDDT} for the best chain."""
    best = {}
    current = {}
    counts = {}
    chain = None
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM"):
                ch = line[21]
                resnum = int(line[22:26])
                if chain is None or ch != chain:
                    chain = ch
                    current = {}
                    counts[ch] = 0
                counts[ch] = counts.get(ch, 0) + 1
                if line[12:16].strip() == "CA":
                    current[resnum] = float(line[60:66])
        # pick the chain with the most CA atoms seen
        best_chain = max(counts, key=counts.get) if counts else None
        return dict(current) if best_chain is not None else {}


def main():
    job = Job("09_structure_annotations")
    mut = pd.read_parquet(PROCESSED / "processed_single_mutations.parquet")
    mut = mut.copy()
    mut["plddt"] = np.nan
    mut["structure_cov"] = np.nan  # per-assay fraction of residues with pLDDT

    pdb_files = {p.stem: p for p in STRUCTURES_DIR.glob("*.pdb")}
    job.info(f"PDB files: {len(pdb_files)}")

    ref = read_reference_dms()
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    target_len_map = {row.DMS_id: len(str(row.target_seq))
                      for row in ref.itertuples()}
    qc_rows = []
    n_with_structure = 0
    for dms_id, g in mut.groupby("DMS_id"):
        uniprot = g["UniProt_ID"].iloc[0]
        target = target_len_map.get(dms_id, 0)
        pdb = pdb_files.get(uniprot)
        if pdb is None:
            qc_rows.append({"UniProt_ID": uniprot, "DMS_id": dms_id,
                            "has_pdb": False, "n_residues_mapped": 0,
                            "target_len": target, "align_pass_rate": 0.0,
                            "median_plddt": np.nan})
            continue
        plddt_map = plddt_from_pdb(pdb)
        n_with_structure += 1
        pos = g["position"].astype(int).to_numpy()
        vals = np.array([plddt_map.get(p, np.nan) for p in pos])
        mut.loc[g.index, "plddt"] = vals
        mapped_residues = sum(1 for p in range(1, target + 1) if p in plddt_map)
        pass_rate = mapped_residues / target
        mut.loc[g.index, "structure_cov"] = pass_rate
        qc_rows.append({"UniProt_ID": uniprot, "DMS_id": dms_id,
                        "has_pdb": True, "n_residues_mapped": mapped_residues,
                        "target_len": target, "align_pass_rate": pass_rate,
                        "median_plddt": float(np.nanmedian(vals)) if len(vals) else np.nan})
        job.info(f"{dms_id} ({uniprot}): mapped {mapped_residues}/{target} "
                 f"residues, pass_rate={pass_rate:.2f}, "
                 f"median_plddt={np.nanmedian(vals):.1f}")

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(STATISTICS / "structure_alignment_QC.csv", index=False)

    n_assays = mut["DMS_id"].nunique()
    frac_struct = mut["plddt"].notna().mean()
    job.info(f"assays with PDB: {n_with_structure}/{n_assays}")
    job.info(f"variants with pLDDT: {frac_struct:.3f}")
    low = qc[(qc.has_pdb) & (qc.align_pass_rate < 0.5)]
    job.info(f"proteins with pass_rate < 0.5: {len(low)} -> {low.UniProt_ID.tolist()}")

    mut.to_parquet(PROCESSED / "structure_features.parquet", index=False)
    job.info(f"structure_features.parquet: {len(mut)} rows")
    job.close()


if __name__ == "__main__":
    main()