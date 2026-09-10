"""08_global_analysis.py

Phase 2 -- first scientific question:
    Q1: Do protein AI models learn the same thing?

Analyses (all with the core 40-model panel, normalized z-scores):
  1. Pairwise model correlation: per-assay Spearman -> median across assays
     (pooled variant-level Spearman reported as a sensitivity check)
  2. Hierarchical clustering (distance = 1 - median Spearman)
  3. PCA of the variant x model prediction matrix (visualisation only)

Outputs:
  results/statistics/Fig1D_model_correlation_matrix.csv   (median Spearman)
  results/statistics/Fig1E_model_clustering.csv           (linkage + order)
  results/statistics/Fig1C_model_pca.csv                  (PC1/PC2 loadings)
  results/statistics/model_correlations_assay_rho.csv     (per-assay rho)
  figures/main/Fig1D_model_correlation_heatmap.png
  figures/main/Fig1E_model_clustering_dendrogram.png
  figures/main/Fig1C_model_pca.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata, spearmanr
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    PROCESSED, STATISTICS, TABLES, FIGURES_MAIN, Job, load_core_panel,
)

CLR = {"evolution": "#E69F00", "single_seq": "#56B4E9",
       "structure": "#009E73", "hybrid": "#CC79A7"}


def main():
    job = Job("08_global_analysis")
    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet")
    core = load_core_panel(job)
    z_cols = [f"z_{m}" for m in core["model"]]
    job.info(f"core panel: {len(core)} models; variants: {len(norm)}")

    # ---- per-assay Spearman correlation matrices ----
    rho_sum = None
    count = 0
    assay_rows = []
    model_list = core["model"].tolist()
    for assay_id, g in norm.groupby("DMS_id", sort=False):
        Z = g[z_cols].to_numpy(dtype=float)
        keep = ~np.isnan(Z).any(axis=1)
        if keep.sum() < 50:
            continue
        Z = Z[keep]
        ranks = np.apply_along_axis(rankdata, 0, Z)
        R = np.corrcoef(ranks, rowvar=False)
        R = np.nan_to_num(R)
        rho_sum = R if rho_sum is None else rho_sum + R
        count += 1
        for i, mi in enumerate(model_list):
            for j, mj in enumerate(model_list):
                if j > i:
                    assay_rows.append({"assay": assay_id, "model_a": mi,
                                       "model_b": mj, "spearman": R[i, j]})
    rho_median = rho_sum / count
    job.info(f"assays contributing to correlation matrix: {count}")

    corr_df = pd.DataFrame(rho_median, index=model_list, columns=model_list)
    corr_df.to_csv(STATISTICS / "Fig1D_model_correlation_matrix.csv")
    pd.DataFrame(assay_rows).to_csv(
        STATISTICS / "model_correlations_assay_rho.csv", index=False)

    # summary: within vs between family
    fam_map = dict(zip(core["model"], core["family"]))
    vals = pd.DataFrame(assay_rows)
    vals["fam_a"] = vals["model_a"].map(fam_map)
    vals["fam_b"] = vals["model_b"].map(fam_map)
    vals["same_family"] = vals["fam_a"] == vals["fam_b"]
    med_same = vals.loc[vals["same_family"], "spearman"].median()
    med_diff = vals.loc[~vals["same_family"], "spearman"].median()
    job.info(f"median Spearman within family: {med_same:.3f}; "
             f"between families: {med_diff:.3f}")

    # ---- hierarchical clustering ----
    dist = 1 - np.clip(rho_median, 0, 1)
    dist = (dist + dist.T) / 2
    np.fill_diagonal(dist, 0)
    cond = squareform(dist)
    Zlink = linkage(cond, method="average")
    from scipy.cluster.hierarchy import leaves_list
    order = leaves_list(Zlink)
    clust_df = pd.DataFrame({"model": model_list,
                             "cluster_order": order,
                             "family": [fam_map[m] for m in model_list]})
    clust_df.to_csv(STATISTICS / "Fig1E_model_clustering.csv", index=False)
    job.info(f"clustering leaves: {order[:10]} ...")

    # ---- PCA (visualisation only) ----
    Zall = norm[z_cols].to_numpy(dtype=float)
    keep = ~np.isnan(Zall).any(axis=1)
    X = Zall[keep]
    pca = PCA(n_components=5, random_state=2026)
    scores = pca.fit(X)
    pca_df = pd.DataFrame({
        "model": model_list,
        "family": [fam_map[m] for m in model_list],
        "PC1": scores.components_[0],
        "PC2": scores.components_[1],
        "PC3": scores.components_[2],
    })
    pca_df.to_csv(STATISTICS / "Fig1C_model_pca.csv", index=False)
    job.info(f"PCA explained variance: {scores.explained_variance_ratio_[:3].round(4)}")
    job.info(f"PC1 |PC1| loadings top: "
             + str(pca_df.reindex(pca_df.PC1.abs().sort_values(ascending=False).index).head(5)[["model", "PC1"]].values.tolist()))

    # ---- figures ----
    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(rho_median[np.ix_(order, order)], cmap="RdBu_r", vmin=0, vmax=1)
    ax.set_xticks(range(len(model_list)), [model_list[i] for i in order], rotation=90, fontsize=6)
    ax.set_yticks(range(len(model_list)), [model_list[i] for i in order], fontsize=6)
    fig.colorbar(im, label="median Spearman")
    ax.set_title("Model prediction correlation (assay-median Spearman)")
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig1D_model_correlation_heatmap.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = [CLR[fam_map[m]] for m in [model_list[i] for i in order]]
    dendrogram(Zlink, labels=model_list, leaf_rotation=90, leaf_font_size=7,
               link_color_func=lambda k: colors[k] if k < len(colors) else "#333333")
    ax.set_ylabel("1 - Spearman")
    ax.set_title("Hierarchical clustering of protein AI models")
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig1E_model_clustering_dendrogram.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 7))
    for fam, sub in pca_df.groupby("family"):
        ax.scatter(sub["PC1"], sub["PC2"], label=fam, color=CLR[fam], s=60, alpha=0.85)
        for _, r in sub.iterrows():
            ax.annotate(r["model"], (r["PC1"], r["PC2"]), fontsize=6, alpha=0.7)
    ax.axhline(0, color="grey", lw=0.5)
    ax.axvline(0, color="grey", lw=0.5)
    ax.set_xlabel(f"PC1 ({scores.explained_variance_ratio_[0]:.1%})")
    ax.set_ylabel(f"PC2 ({scores.explained_variance_ratio_[1]:.1%})")
    ax.legend()
    ax.set_title("PCA of model predictions (variant-level)")
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig1C_model_pca.png", dpi=200)
    plt.close(fig)

    job.info("figures written: Fig1C, Fig1D, Fig1E")
    job.close()


if __name__ == "__main__":
    main()