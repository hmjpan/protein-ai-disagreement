"""06_normalize_scores.py

Normalize every model prediction to a per-assay *deleteriousness percentile*.

Direction is taken ONLY from model_direction.csv (official ProteinGym config
where available; curated fallback otherwise). It is NEVER inferred from DMS
labels. For each assay and model:
    r_im = percentile rank of raw score within the assay
    u_im = 1 - r_im        if higher raw score = higher fitness
         = r_im            if higher raw score = higher deleteriousness
    z_im = Phi^-1(clip(u_im, 0.001, 0.999))

Experimental deleteriousness for validation:
    Y_i  = 1 - rank(DMS_score_i within assay)      (ProteinGym: higher DMS_score
                                                    = higher fitness, official)

A direction-diagnostic table (median within-assay Spearman between each model's
rank and DMS rank) is written for QC purposes only -- never used to flip scores.

Inputs : data/processed/panel_scores.parquet
         data/processed/model_direction.csv
Outputs: data/processed/normalized_scores.parquet
         results/statistics/direction_diagnostic.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job  # noqa: E402

SCORE_COLS_EXCLUDE = {"DMS_id", "UniProt_ID", "mutant", "wt_aa", "position",
                      "mut_aa", "DMS_score", "DMS_score_bin",
                      "n_core_models_scored", "frac_core_models_scored"}


def match_direction(model: str, dir_df: pd.DataFrame):
    """Exact match first, then substring match against curated names."""
    hit = dir_df[dir_df["model_name"] == model]
    if len(hit):
        return hit.iloc[0]
    for _, row in dir_df.iterrows():
        if row["model_name"].lower() in model.lower():
            return row
    return None


def main():
    job = Job("06_normalize_scores")
    df = pd.read_parquet(PROCESSED / "panel_scores.parquet")
    job.info(f"panel scores: {len(df)} rows, {df.DMS_id.nunique()} assays")

    dir_df = pd.read_csv(PROCESSED / "model_direction.csv")
    meta = pd.read_csv(PROCESSED / "model_metadata.csv")
    official_models = set(meta["model_name"])
    model_cols = [c for c in df.columns
                  if c in official_models and c not in SCORE_COLS_EXCLUDE]
    job.info(f"model columns (official list, {len(model_cols)}): "
             + ", ".join(model_cols))

    # apply directions; collect diagnostics
    diag_rows = []
    norm_frames = []
    used_directions = {}
    for assay_id, g in df.groupby("DMS_id", sort=False):
        g = g.sort_values("mutant").reset_index(drop=True)
        dms_rank = rankdata(g["DMS_score"], method="average") / len(g)
        out = g[["DMS_id", "UniProt_ID", "mutant", "wt_aa", "position",
                 "mut_aa", "DMS_score", "DMS_score_bin"]].copy()
        out["Y_deleter"] = 1 - dms_rank
        for m in model_cols:
            s = g[m].to_numpy(dtype=float)
            ok = ~np.isnan(s)
            r = np.full(len(g), np.nan)
            if ok.any():
                r[ok] = rankdata(s[ok], method="average") / ok.sum()
            row = match_direction(m, dir_df)
            if row is None:
                job.info(f"NO DIRECTION for {m} -- excluding from u/z (NaN)")
                used_directions[m] = "unknown"
                out[f"u_{m}"] = np.nan
                out[f"z_{m}"] = np.nan
                continue
            higher_means = row.get("released_higher_means", "fitness")
            # value semantics: mutational_direction = higher_is_better -> higher=fitness
            flip = (row["mutational_direction"] == "higher_is_better")
            u = (1 - r) if flip else r
            z = np.full(len(g), np.nan)
            valid = ~np.isnan(u)
            uu = np.clip(u[valid], 0.001, 0.999)
            from scipy.stats import norm
            z[valid] = norm.ppf(uu)
            out[f"u_{m}"] = u
            out[f"z_{m}"] = z
            used_directions[m] = row["mutational_direction"]
            if ok.sum() > 50:
                rho, _ = spearmanr(s[ok], g.loc[ok, "DMS_score"].to_numpy())
                diag_rows.append({"model": m, "assay": assay_id,
                                  "n": int(ok.sum()), "spearman_vs_dms": rho})
        norm_frames.append(out)

    norm = pd.concat(norm_frames, ignore_index=True)
    norm.to_parquet(PROCESSED / "normalized_scores.parquet", index=False)
    job.info(f"normalized_scores.parquet: {len(norm)} rows")

    diag = pd.DataFrame(diag_rows)
    diag_summary = diag.groupby("model").agg(
        n_assays=("assay", "nunique"),
        median_rho=("spearman_vs_dms", "median"),
        mean_rho=("spearman_vs_dms", "mean")).reset_index()
    diag_summary.to_csv(STATISTICS / "direction_diagnostic.csv", index=False)
    job.info("direction_diagnostic.csv written (QC only)")
    neg = diag_summary[diag_summary["median_rho"] < -0.05]
    if len(neg):
        job.info(f"WARNING: {len(neg)} models show negative median DMS correlation "
                 f"-- check direction table: {neg['model'].tolist()}")
    job.info(f"direction usage: {sorted(set(used_directions.values()))}")
    job.close()


if __name__ == "__main__":
    main()