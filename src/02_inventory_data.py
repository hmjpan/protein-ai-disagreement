"""02_inventory_data.py

Produce the three Phase-1 inventory tables.

Table S1 : DMS assay inventory (DMS_id, UniProt_ID, protein, taxon, selection
           type, number of variants, number of single mutants, sequence length,
           MSA depth category)
Table S2 : AI model inventory (model, family, modality, assays covered, variant
           coverage, missing rate)
Table S3 : model x assay coverage matrix (heatmap-ready)

Inputs : data/raw/reference/DMS_substitutions.csv
         data/raw/proteingym_dms/DMS_ProteinGym_substitutions/*.csv
         data/raw/proteingym_scores/zero_shot_substitutions_scores/*.csv
Outputs: results/tables/TableS1_assay_inventory.csv
         results/tables/TableS2_model_inventory.csv
         results/tables/TableS3_model_assay_coverage.csv
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    DMS_DIR, PROCESSED, SCORES_DIR, TABLES, Job, read_reference_dms,
)

SINGLE_RE = re.compile(r"^[A-Z]\d+[A-Z]$")

# Fallback family map (only used if official metadata is unavailable)
FAMILY_FALLBACK = {
    "EVE": "evolution", "GEMME": "evolution", "EVmutation": "evolution",
    "DeepSequence": "evolution", "MSA_Transformer": "evolution", "MSA-VAE": "evolution",
    "ARD": "evolution", "Bayes_ML": "evolution", "Site_Independent": "evolution",
    "Wavenet": "evolution", "Protriever": "evolution", "PoET": "evolution",
    "TranceptEVE": "evolution",
    "ESM-1v": "single_seq", "ESM1b": "single_seq", "ESM2": "single_seq",
    "ESM3": "single_seq", "ESM C": "single_seq", "ProGen2": "single_seq",
    "Progen": "single_seq", "Progen3": "single_seq", "CARP": "single_seq",
    "RITA": "single_seq", "UniRep": "single_seq", "LSTM": "single_seq",
    "ProtGPT2": "single_seq", "VESPA": "single_seq", "VespaG": "single_seq",
    "Wavenet": "evolution",
    "ESM-IF1": "structure", "ProteinMPNN": "structure", "MIF": "structure",
    "MIF-ST": "structure", "MIFST": "structure",
    "SaProt": "hybrid", "ProtSSN": "hybrid", "S2F": "hybrid", "S3F": "hybrid",
    "VenusREM": "hybrid", "Escott": "hybrid", "ESCOTT": "hybrid",
    "Tranception": "hybrid", "RSALOR": "hybrid", "MULAN": "hybrid",
    "ProSST": "hybrid", "AIDO": "hybrid",
}


def load_family_map():
    meta = PROCESSED / "model_metadata.csv"
    if meta.exists():
        df = pd.read_csv(meta)
        return dict(zip(df["model_name"], df["family"]))
    return {}


FAMILY_OFFICIAL = load_family_map()


def family_of(name: str) -> str:
    if name in FAMILY_OFFICIAL:
        return FAMILY_OFFICIAL[name]
    for key, fam in FAMILY_FALLBACK.items():
        if key.lower() in name.lower():
            return fam
    return "other"


def main():
    job = Job("02_inventory_data")
    ref = read_reference_dms()
    job.info(f"reference rows: {len(ref)}; columns: {list(ref.columns)}")
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    ref_map = ref.set_index("DMS_id")

    dms_files = sorted(DMS_DIR.glob("*.csv"))
    score_files = sorted(SCORES_DIR.glob("*.csv"))
    job.info(f"DMS files: {len(dms_files)}; score files: {len(score_files)}")
    dms_ids = {p.stem for p in dms_files}
    score_ids = {p.stem for p in score_files}
    common = sorted(dms_ids & score_ids)
    job.info(f"assay overlap (DMS n scores): {len(common)}")

    # ---- Table S1: assay inventory (variant counts from DMS files) ----
    rows1 = []
    for p in dms_files:
        dms_id = p.stem
        meta = ref_map.loc[dms_id] if dms_id in ref_map.index else None
        df = pd.read_csv(p)
        n_var = len(df)
        n_single = df["mutant"].astype(str).str.match(SINGLE_RE).sum()
        seq_len = None
        msa = None
        taxon = None
        selection = None
        if meta is not None:
            seq_len = len(str(meta.get("target_seq", ""))) if pd.notna(meta.get("target_seq", None)) else None
            msa = meta.get("MSA_depth", None)
            taxon = meta.get("taxon", None)
            for c in ("selection_type", "assay_type", "phenotype"):
                if c in meta.index and pd.notna(meta.get(c, None)):
                    selection = meta[c]
                    break
        rows1.append({
            "DMS_id": dms_id,
            "UniProt_ID": meta["UniProt_ID"] if meta is not None else None,
            "taxon": taxon,
            "selection_type": selection,
            "n_variants": n_var,
            "n_single_mutants": n_single,
            "sequence_length": seq_len,
            "MSA_depth": msa,
            "has_scores": dms_id in score_ids,
        })
    t1 = pd.DataFrame(rows1)
    t1.to_csv(TABLES / "TableS1_assay_inventory.csv", index=False)
    job.info(f"TableS1: {len(t1)} assays; single-mutant fraction overall "
             f"{t1.n_single_mutants.sum() / t1.n_variants.sum():.3f}")

    # ---- Table S2: model inventory + Table S3: coverage matrix ----
    model_cols = []
    coverage = {}
    for p in score_files:
        header = pd.read_csv(p, nrows=1).columns.tolist()
        model_cols.append(header[1:])
    all_models = []
    for cols in model_cols:
        for c in cols:
            if c not in all_models:
                all_models.append(c)
    job.info(f"distinct model columns found: {len(all_models)}")

    cov = pd.DataFrame(0.0, index=sorted(common), columns=all_models)
    for p in score_files:
        dms_id = p.stem
        if dms_id not in cov.index:
            continue
        df = pd.read_csv(p, usecols=lambda c: c == "mutant" or c in all_models)
        n = len(df)
        for c in all_models:
            if c in df.columns:
                cov.loc[dms_id, c] = df[c].notna().sum() / n
    t3 = cov.copy()
    t3.insert(0, "DMS_id", cov.index)
    t3.to_csv(TABLES / "TableS3_model_assay_coverage.csv", index=False)
    job.info(f"TableS3 written: {t3.shape}")

    n_assays_with_any = (cov > 0).sum(axis=0)
    median_cov = cov.where(cov > 0).median(axis=0)
    rows2 = []
    for c in all_models:
        rows2.append({
            "model_name": c,
            "model_family": family_of(c),
            "n_assays_with_scores": int(n_assays_with_any[c]),
            "median_variant_coverage_among_scored": float(median_cov[c]) if pd.notna(median_cov[c]) else 0.0,
            "missing_rate": float(1 - (cov[c] > 0).mean()),
        })
    t2 = pd.DataFrame(rows2).sort_values("model_family")
    t2.to_csv(TABLES / "TableS2_model_inventory.csv", index=False)
    job.info(f"TableS2 written: {len(t2)} models; families: "
             + str(t2.model_family.value_counts().to_dict()))
    job.close()


if __name__ == "__main__":
    main()