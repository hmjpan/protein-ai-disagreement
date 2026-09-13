"""07_compute_disagreement.py

Compute the core scientific variable -- AI model disagreement -- for every
variant in the primary analysis set.

  D_std   = sqrt( 1/M * sum_m (z_im - mean_z_i)^2 )   over core-panel models
  D_MAD   = median_m |z_im - median(z_i)|             robust version
  family centroids: S_seq, S_evo, S_struct, S_hybrid (mean z within family)
  D_seq_evo / D_seq_struct / D_evo_struct = |centroid_a - centroid_b|

Model panel comes from results/tables/model_panel.csv (panel == 'core').
Family labels come from model_panel.csv.

Inputs : data/processed/normalized_scores.parquet
         results/tables/model_panel.csv
Outputs: data/processed/disagreement_scores.parquet
         results/tables/TableS4_disagreement_summary.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, TABLES, Job  # noqa: E402

BASE_COLS = ["DMS_id", "UniProt_ID", "mutant", "wt_aa", "position", "mut_aa",
             "DMS_score", "DMS_score_bin", "Y_deleter"]
FAMILIES = ["evolution", "single_seq", "structure", "hybrid"]


def main():
    job = Job("07_compute_disagreement")
    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet")
    panel = pd.read_csv(TABLES / "model_panel.csv")
    core = panel[panel["panel"] == "core"].copy()
    job.info(f"core panel: {len(core)} models -> {core['family'].value_counts().to_dict()}")

    z_cols = {m: f"z_{m}" for m in core["model"]}
    missing = [m for m, c in z_cols.items() if c not in norm.columns]
    if missing:
        job.info(f"dropping core models missing from normalized data: {missing}")
        z_cols = {m: c for m, c in z_cols.items() if m not in missing}

    Z = norm[list(z_cols.values())].to_numpy(dtype=float)
    m_models = Z.shape[1]

    with np.errstate(invalid="ignore"):
        D_std = np.nanstd(Z, axis=1)
        med = np.nanmedian(Z, axis=1)
        D_mad = np.nanmedian(np.abs(Z - med[:, None]), axis=1)

    out = norm[BASE_COLS].copy()
    out["D_std"] = D_std
    out["D_mad"] = D_mad

    fam_centroids = {}
    for fam in FAMILIES:
        cols = [f"z_{m}" for m in core["model"] if core.loc[core["model"] == m, "family"].iloc[0] == fam]
        cols = [c for c in cols if c in norm.columns]
        if cols:
            fam_centroids[f"S_{fam}"] = np.nanmean(norm[cols].to_numpy(dtype=float), axis=1)
            out[f"S_{fam}"] = fam_centroids[f"S_{fam}"]

    pairs = [("D_seq_evo", "S_single_seq", "S_evolution"),
             ("D_seq_struct", "S_single_seq", "S_structure"),
             ("D_evo_struct", "S_evolution", "S_structure"),
             ("D_seq_hybrid", "S_single_seq", "S_hybrid"),
             ("D_evo_hybrid", "S_evolution", "S_hybrid"),
             ("D_struct_hybrid", "S_structure", "S_hybrid")]
    for out_name, a, b in pairs:
        if a in out.columns and b in out.columns:
            out[out_name] = np.abs(out[a] - out[b])

    out.to_parquet(PROCESSED / "disagreement_scores.parquet", index=False)
    job.info(f"disagreement_scores.parquet: {len(out)} rows; "
             f"mean D_std={out['D_std'].mean():.4f}; non-null D_std={out['D_std'].notna().mean():.3f}")

    summ = out.groupby("DMS_id").agg(
        n_variants=("D_std", "size"),
        median_D=("D_std", "median"),
        p95_D=("D_std", lambda s: s.quantile(0.95)),
        frac_highD=("D_std", lambda s: (s >= s.quantile(0.9)).mean()),
    ).reset_index()
    summ.to_csv(TABLES / "TableS4_disagreement_summary.csv", index=False)
    job.info(f"TableS4: {len(summ)} assays")
    job.close()


if __name__ == "__main__":
    main()