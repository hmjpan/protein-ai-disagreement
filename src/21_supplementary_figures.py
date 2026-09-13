"""21_supplementary_figures.py

Create the remaining supplementary figures from existing statistics tables:
  FigS11  leave-one-model-out sensitivity (LOMO median rhos)
  FigS12  u-basis vs z-basis disagreement-error association
  FigS13  per-assay disagreement-error rho by assay selection type and MSA
          depth
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    PROCESSED, STATISTICS, FIGURES_SUPP, Job, load_core_panel,
)


def main():
    job = Job("21_supplementary_figures")

    # ---- FigS11: LOMO sensitivity ----
    lomo = pd.read_csv(STATISTICS / "disagreement_lomo_sensitivity.csv")
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(lomo["median_rho"], np.arange(len(lomo)), s=25)
    ax.axvline(lomo["median_rho"].median(), color="red", ls="--",
               label=f"median {lomo['median_rho'].median():.4f}")
    ax.set_yticks(range(len(lomo)))
    ax.set_yticklabels(lomo["left_out"], fontsize=6)
    ax.set_xlabel("median assay-level Spearman(D_std, error) with model left out")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES_SUPP / "FigS6_lomo_sensitivity.png", dpi=300)
    plt.close(fig)
    job.info("FigS11 written")

    # ---- FigS12: u-basis vs z-basis ----
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet")
    core = load_core_panel(job)
    u_cols = [f"u_{m}" for m in core["model"]]
    z_cols = [f"z_{m}" for m in core["model"]]
    m_df = disc.merge(norm[["DMS_id", "mutant"] + u_cols + z_cols],
                      on=["DMS_id", "mutant"], how="left")
    S = m_df[u_cols].mean(axis=1, skipna=True)
    m_df["error"] = (S - m_df["Y_deleter"]).abs()
    Z = m_df[z_cols].to_numpy(dtype=float)
    U = m_df[u_cols].to_numpy(dtype=float)
    rows = []
    for dms_id, g in m_df.groupby("DMS_id"):
        if g["D_std"].notna().sum() < 50:
            continue
        rho_z = spearmanr(g["D_std"], g["error"])[0]
        rho_u = spearmanr(np.nanstd(U[g.index], axis=1), g["error"])[0]
        rows.append({"DMS_id": dms_id, "rho_z": rho_z, "rho_u": rho_u})
    rdf = pd.DataFrame(rows)
    rdf.to_csv(STATISTICS / "disagreement_error_scale_comparison.csv", index=False)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.boxplot([rdf["rho_z"], rdf["rho_u"]], labels=["z-basis (primary)", "u-basis"])
    ax.set_ylabel("assay-level Spearman(disagreement, error)")
    fig.tight_layout()
    fig.savefig(FIGURES_SUPP / "FigS7_scale_comparison.png", dpi=300)
    plt.close(fig)
    job.info(f"FigS12 written (median z={np.median(rdf.rho_z):.4f} "
             f"u={np.median(rdf.rho_u):.4f})")

    # ---- FigS13: rho by assay type / MSA depth ----
    ref = pd.read_csv(PROCESSED.parent / "raw" / "reference" / "DMS_substitutions.csv")
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    rdf = rdf.merge(ref[["DMS_id", "coarse_selection_type", "MSA_Neff_L_category"]],
                    on="DMS_id", how="left")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, col in zip(axes, ["coarse_selection_type", "MSA_Neff_L_category"]):
        groups = []
        labs = []
        for k, g in rdf.groupby(col):
            groups.append(g["rho_z"].to_numpy())
            labs.append(f"{k}\n(n={len(g)})")
        ax.boxplot(groups, labels=labs)
        ax.set_ylabel("Spearman(D, error)")
        ax.set_title(col)
    fig.tight_layout()
    fig.savefig(FIGURES_SUPP / "FigS8_rho_by_assay_context.png", dpi=300)
    plt.close(fig)
    job.info("FigS13 written")
    job.close()


if __name__ == "__main__":
    main()