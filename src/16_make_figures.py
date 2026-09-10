"""16_make_figures.py

Recreate ALL main figures from machine-readable results tables (v1.3 numbering:
figure 1-6 main; supplementary S1-S10). One-command reproducible; numeric
content is never hand-edited.

Usage:  python src/16_make_figures.py
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
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import STATISTICS, PROCESSED, FIGURES_MAIN, FIGURES_SUPP, Job, load_core_panel  # noqa: E402

CLR = {"evolution": "#E69F00", "single_seq": "#56B4E9",
       "structure": "#009E73", "hybrid": "#CC79A7", "other": "#999999"}


def save(fig, name, supp=False):
    fig.savefig((FIGURES_SUPP if supp else FIGURES_MAIN) / name, dpi=300)
    plt.close(fig)


def main():
    job = Job("16_make_figures")

    # ---------- figure 1: model landscape ----------
    corr = pd.read_csv(STATISTICS / "Fig1D_model_correlation_matrix.csv", index_col=0)
    clust = pd.read_csv(STATISTICS / "Fig1E_model_clustering.csv")
    pca = pd.read_csv(STATISTICS / "Fig1C_model_pca.csv")
    order = clust.sort_values("cluster_order")["model"].tolist()
    R = corr.loc[order, order].to_numpy()

    fig, ax = plt.subplots(figsize=(10, 9))
    im = ax.imshow(R, cmap="RdBu_r", vmin=0, vmax=1)
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=6)
    ax.set_yticklabels(order, fontsize=6)
    fig.colorbar(im, ax=ax, label="median assay-level Spearman")
    fig.tight_layout()
    save(fig, "Fig1D_model_correlation_heatmap.png")

    dist = 1 - np.clip(R, 0, 1)
    dist = (dist + dist.T) / 2
    np.fill_diagonal(dist, 0)
    Z = linkage(squareform(dist), method="average")
    fig, ax = plt.subplots(figsize=(10, 6))
    fam = dict(zip(clust["model"], clust["family"]))
    cols = [CLR[fam[m]] for m in order]
    dendrogram(Z, labels=order, leaf_rotation=90, leaf_font_size=7,
               link_color_func=lambda k: cols[k % len(cols)])
    ax.set_ylabel("1 - median Spearman")
    fig.tight_layout()
    save(fig, "Fig1E_model_clustering_dendrogram.png")

    fig, ax = plt.subplots(figsize=(8, 7))
    for f, sub in pca.groupby("family"):
        ax.scatter(sub["PC1"], sub["PC2"], label=f, color=CLR[f], s=70,
                   alpha=0.85, edgecolor="k", linewidth=0.3)
        for _, r in sub.iterrows():
            ax.annotate(r["model"], (r["PC1"], r["PC2"]), fontsize=7.5, alpha=0.7)
    ax.axhline(0, color="grey", lw=0.5)
    ax.axvline(0, color="grey", lw=0.5)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "Fig1C_model_pca.png")

    # ---------- figure 2: structural organization ----------
    eff = pd.read_csv(STATISTICS / "plddt_disagreement_protein_effects.csv")
    rde = pd.read_csv(STATISTICS / "residue_level_plddt_effects.csv")
    fig, ax = plt.subplots(figsize=(6, 5))
    v = rde["rho_De_plddt"].dropna()
    ax.hist(v, bins=35, color="#4477AA", alpha=0.85)
    ax.axvline(v.median(), color="red", ls="--", lw=2,
               label=f"median = {v.median():.3f}")
    ax.set_xlabel("within-protein residue-level Spearman\n(pLDDT, evolution-vs-structure disagreement)")
    ax.set_ylabel("proteins")
    ax.annotate(f"{(v < 0).mean() * 100:.1f}% of {len(v)}\nproteins negative",
                xy=(0.03, 0.92), xycoords="axes fraction", fontsize=9)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    save(fig, "Fig2A_plddt_disagreement.png")

    ez = pd.read_csv(STATISTICS / "plddt_effect_sizes.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    v = ez["dD_evo_struct_decile"].dropna()
    axes[0].hist(v, bins=40, color="#4477AA", alpha=0.85)
    axes[0].axvline(v.median(), color="red", ls="--", lw=1.5)
    axes[0].set_xlabel("low-minus-high pLDDT decile D_evo_struct")
    v2 = ez["enrich_OR"].dropna().clip(upper=10)
    axes[1].hist(v2, bins=40, color="#88CCEE", alpha=0.85)
    axes[1].axvline(ez["enrich_OR"].median(), color="red", ls="--", lw=1.5)
    axes[1].set_xlabel("enrichment OR (low-pLDDT decile, top-D decile)")
    fig.tight_layout()
    save(fig, "Fig2B_effect_sizes.png")

    ds = pd.read_csv(STATISTICS / "disagreement_definition_plddt_summary.csv")
    meds = np.load(STATISTICS / "panel_residue_level_B_medians.npy")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    y = np.arange(len(ds))
    axes[0].errorbar(ds["median_rho"], y,
                     xerr=[ds["median_rho"] - ds["ci_low"], ds["ci_high"] - ds["median_rho"]],
                     fmt="o", ms=7, capsize=4)
    axes[0].axvline(0, color="grey", lw=1)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(ds["definition"], fontsize=8)
    axes[0].set_xlabel("median residue-level rho(D, pLDDT)")
    vals, edges = np.histogram(meds, bins=np.arange(-0.65, 0.05, 0.05))
    axes[1].bar(edges[:-1], vals, width=0.045, align="edge", color="#4477AA")
    axes[1].axvline(-0.05, color="red", ls="--", lw=1.5)
    axes[1].set_xlabel("median per balanced resampling")
    axes[1].set_ylabel("iterations")
    fig.tight_layout()
    save(fig, "Fig2C_robustness.png")

    idr = pd.read_csv(STATISTICS / "disorder_disagreement_effects.csv")
    fig, ax = plt.subplots(figsize=(6, 5))
    for i, col in enumerate(["D_evo_struct", "D_seq_struct", "D_std"]):
        ax.boxplot([idr[f"{col}_idr"] - idr[f"{col}_ord"]], positions=[i], widths=0.5)
    ax.axhline(0, color="red", ls="--", lw=1)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["D_evo_struct", "D_seq_struct", "D_std"])
    ax.set_ylabel("median disagreement (IDR - ordered), per protein")
    fig.tight_layout()
    save(fig, "Fig2D_disorder_disagreement.png")

    # ---------- figure 3: regimes ----------
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    cols = ["DMS_id", "S_single_seq", "S_evolution", "S_structure"]
    full = disc[[c for c in cols if disc[c].nunique() < 1000 or True]].loc[:, [c for c in dict.fromkeys(cols)]]
    rng = np.random.default_rng(2026)
    idx = rng.choice(len(full), size=6000, replace=False)
    samp = full.iloc[idx]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, (a, b) in zip(axes, [("S_evolution", "S_single_seq"),
                                 ("S_evolution", "S_structure")]):
        ax.scatter(samp[a].iloc[::40], samp[b].iloc[::40], s=6, alpha=0.4, c="#555555")
        ax.axhline(0, color="grey", lw=0.4)
        ax.axvline(0, color="grey", lw=0.4)
        ax.set_xlabel(a)
        ax.set_ylabel(b)
    fig.tight_layout()
    save(fig, "Fig3A_regime_scatter.png")

    rc = pd.read_csv(STATISTICS / "regime_centroids.csv")
    fig, ax = plt.subplots(figsize=(8, 5))
    order_r = rc.sort_values("D_std", key=lambda s: -s)["name"] if "D_std" in rc else rc["name"]
    summ = pd.read_csv(STATISTICS / "regime_summary.csv")
    sm = summ.set_index("regime")
    x = np.arange(len(sm))
    ax.bar(x - 0.15, sm["median_plddt"], 0.3, label="median pLDDT", color="#4477AA")
    ax.bar(x + 0.15, sm["median_Y"] * 100, 0.3, label="experimental Y x100", color="#CC6677")
    ax.set_xticks(x)
    ax.set_xticklabels(sm.index, rotation=25, fontsize=8)
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "Fig3B_regime_plddt.png")

    fig, ax = plt.subplots(figsize=(6, 4.5))
    for _, r in rc.iterrows():
        ax.scatter(r["S_evolution"], r["S_structure"], s=140,
                   c=r["S_single_seq"], cmap="coolwarm")
        ax.annotate(r["name"], (r["S_evolution"], r["S_structure"]), fontsize=8)
    ax.set_xlabel("S_evolution")
    ax.set_ylabel("S_structure")
    ax.axhline(0, color="grey", lw=0.5)
    ax.axvline(0, color="grey", lw=0.5)
    fig.tight_layout()
    save(fig, "Fig3C_regime_biology.png")

    # ---------- figure 4: reproducibility ----------
    rep = pd.read_csv(STATISTICS / "regime_reproducibility.csv")
    rep2 = pd.read_csv(STATISTICS / "regime_phenotype_replication.csv")
    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.scatter(rep["chance_agreement"], rep["agreement"], s=20, alpha=0.7, c="#4477AA")
    ax.plot([0, 1], [0, 1], "r--", lw=1)
    ax.set_xlabel("chance agreement (marginal product)")
    ax.set_ylabel("observed position-level regime agreement")
    if len(rep2):
        axins = ax.inset_axes([0.12, 0.10, 0.42, 0.38])
        v = rep2["rho_regime_Y"].dropna()
        axins.hist(v, bins=25, color="#CC6677", alpha=0.85)
        axins.axvline(v.median(), color="black", ls="--", lw=1)
        axins.set_xlabel("rho(regime-level Y)", fontsize=7)
        axins.tick_params(labelsize=6)
        axins.set_title(f"phenotype replication (median {v.median():.2f})", fontsize=7)
    fig.tight_layout()
    save(fig, "Fig4A_regime_reproducibility.png")

    # ---------- figure 5: expertise ----------
    adv = pd.read_csv(STATISTICS / "family_advantage_plddt.csv")
    labs = ["seq", "evo", "struct"]
    bins = ["pLDDT<50", "50-70", "70-90", ">=90"]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for i, lab in enumerate(labs):
        sub = adv[adv.family == lab]
        vals = [sub.loc[sub.plddt_bin == b, "median_adv"].iloc[0] for b in bins]
        ax.bar(np.arange(len(bins)) + i * 0.25, vals, 0.25, label=lab,
               color=CLR[{"seq": "single_seq", "evo": "evolution", "struct": "structure"}[lab]])
    ax.set_xticks(np.arange(len(bins)) + 0.25)
    ax.set_xticklabels(bins)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("median family advantage")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "Fig5A_family_advantage_plddt.png")

    madv = pd.read_csv(STATISTICS / "mechanics_family_advantage.csv")
    ctxs = ["core", "surface", "helix", "sheet", "loop"]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for i, lab in enumerate(labs):
        sub = madv[madv.family == lab]
        vals = [sub.loc[sub.context == c, "median_adv"].iloc[0] for c in ctxs]
        ax.bar(np.arange(len(ctxs)) + i * 0.25, vals, 0.25, label=lab)
    ax.set_xticks(np.arange(len(ctxs)) + 0.25)
    ax.set_xticklabels(ctxs)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("median family advantage")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "Fig5B_family_advantage_mechanics.png")

    dfa = pd.read_csv(STATISTICS / "disorder_family_advantage.csv")
    ctxs = ["IDR", "ordered", "lowconf", "highconf"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for i, lab in enumerate(labs):
        sub = dfa[dfa.family == lab]
        vals = [sub.loc[sub.context == c, "median_adv"].iloc[0] for c in ctxs]
        ax.bar(np.arange(len(ctxs)) + i * 0.25, vals, 0.25, label=lab)
    ax.set_xticks(np.arange(len(ctxs)) + 0.25)
    ax.set_xticklabels(ctxs)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("median family advantage")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "Fig5C_family_advantage_disorder.png")

    # ---------- figure 6: clinical ----------
    diff = pd.read_csv(STATISTICS / "clinical_disagreement_summary.csv")
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(diff["D_correct"], diff["D_incorrect"], s=12, alpha=0.6, c="#882255")
    lim = [0, max(diff[["D_correct", "D_incorrect"]].max().max(), 1)]
    ax.plot(lim, lim, "r--", lw=1)
    ax.set_xlabel("median disagreement (correctly classified)")
    ax.set_ylabel("median disagreement (misclassified)")
    fig.tight_layout()
    save(fig, "Fig6B_clinical_disagreement.png")

    tr = pd.read_csv(STATISTICS / "clinical_transfer_performance.csv")
    non = pd.read_csv(STATISTICS / "clinical_transfer_nonoverlap.csv")
    fig, ax = plt.subplots(figsize=(9, 5.5))
    data = [tr["auroc_uniform"], tr["auroc_biogate"], tr["auroc_best"],
            non["auroc_uniform"], non["auroc_biogate"], non["auroc_best"]]
    labs = [f"Uniform\nfull (n={len(tr)})", f"BioGate\nfull", f"Best single\nfull",
            f"Uniform\nstrict (n={len(non)})", f"BioGate\nstrict", f"Best single\nstrict"]
    ax.boxplot(data, labels=labs)
    ax.set_ylabel("protein-level AUROC")
    ax.set_title("Clinical transfer: full vs strict non-overlap set")
    fig.tight_layout()
    save(fig, "Fig6A_clinical_transfer.png")

    # ---------- control figures: exp-structure gain + SSEmb replication ----------
    eg = pd.read_csv(STATISTICS / "expstruct_input_gain.csv")
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.scatter(eg["mean_plddt"], eg["gain"], s=22, alpha=0.8, c="#E69F00")
    ax.axhline(0, color="grey", lw=1)
    ax.set_xlabel("assay mean pLDDT (AlphaFold)")
    ax.set_ylabel("ESM-IF1 Spearman gain\n(experimental - AlphaFold2 structures)")
    ax.set_title("Experimental-structure control (assay level)")
    fig.tight_layout()
    save(fig, "FigS10_expstruct_gain.png", supp=True)

    sr = pd.read_csv(STATISTICS / "ssemb_residue_replication.csv")
    fig, ax = plt.subplots(figsize=(6, 4.5))
    v2 = sr["rho"].dropna()
    ax.hist(v2, bins=35, color="#009E73", alpha=0.85)
    ax.axvline(v2.median(), color="red", ls="--", lw=2,
               label=f"median = {v2.median():.3f}")
    ax.set_xlabel("within-protein residue-level Spearman\n(pLDDT, evolution-vs-SSEmb disagreement)")
    ax.set_ylabel("proteins")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "FigS11_ssemb_replication.png", supp=True)

    # ---------- supplementary ----------
    rho = pd.read_csv(STATISTICS / "disagreement_error_assay_rho.csv")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(rho["rho"], bins=50, color="#555555", alpha=0.85)
    ax.axvline(rho["rho"].median(), color="red", ls="--", label=f"median {rho['rho'].median():.3f}")
    ax.set_xlabel("assay-level Spearman(D_std, error)")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "FigS1_disagreement_distribution.png", supp=True)

    dec = pd.read_csv(STATISTICS / "disagreement_deciles.csv", index_col=0)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(dec.index + 1, dec["mean"], color="#CC6677")
    ax.set_xlabel("within-assay disagreement decile")
    ax.set_ylabel("mean |ensemble - DMS| error")
    fig.tight_layout()
    save(fig, "FigS2_disagreement_deciles.png", supp=True)

    dom = pd.read_csv(STATISTICS / "domain_boundary_disagreement.csv", index_col=0)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(dom.index.astype(str), dom["median"], color="#88CCEE")
    ax.set_xlabel("distance to domain boundary")
    ax.set_ylabel("median disagreement")
    fig.tight_layout()
    save(fig, "FigS3_domain_boundary.png", supp=True)

    # BioGate supplementary
    perf = pd.read_csv(STATISTICS / "biogate_performance.csv")
    fig, ax = plt.subplots(figsize=(7, 5))
    order_m = list(perf.sort_values("median_rho", ascending=False)["method"].unique())[:5]
    box = {}
    cv = pd.read_parquet(STATISTICS / "biogate_cv_predictions.parquet")
    methods = ["uniform_ensemble", "linear_stacking", "xgboost", "biogate"]
    for uni, g in cv.groupby("UniProt_ID"):
        if len(g) < 20:
            continue
        for m in methods:
            box.setdefault(m, []).append(spearmanr(g[m], g["y"])[0])
    data = [np.array(box[m]) for m in order_m if m in box]
    order_m = [m for m in order_m if m in box]
    bp = ax.boxplot(data, labels=order_m, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("#DDDDDD")
    ax.tick_params(axis="x", rotation=20)
    ax.set_ylabel("protein-level Spearman rho")
    fig.tight_layout()
    save(fig, "FigS4_biogate_performance.png", supp=True)

    ab = pd.read_csv(STATISTICS / "biogate_ablation.csv")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(ab["ablation"], ab["median_rho"], color="#88CCEE")
    ax.tick_params(axis="x", rotation=20)
    ax.set_ylabel("median protein-level rho")
    fig.tight_layout()
    save(fig, "FigS5_biogate_ablation.png", supp=True)

    meff = pd.read_csv(STATISTICS / "mechanics_disagreement_effects.csv")
    fig, ax = plt.subplots(figsize=(8, 5))
    cols = ["D_evo_struct_core", "D_evo_struct_surface", "D_evo_struct_helix",
            "D_evo_struct_sheet", "D_evo_struct_loop"]
    for i, c in enumerate(cols):
        ax.boxplot([meff[c].dropna()], positions=[i], widths=0.6)
    ax.set_xticks(range(5))
    ax.set_xticklabels(["core", "surface", "helix", "sheet", "loop"])
    ax.set_ylabel("median evolution-vs-structure disagreement")
    fig.tight_layout()
    save(fig, "FigS9_mechanics_null.png", supp=True)

    job.info("all figures regenerated (v1.3 numbering: figures 1-6 + supp S1-S10)")
    job.close()


if __name__ == "__main__":
    main()