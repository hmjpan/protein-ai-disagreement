"""12_disagreement_regimes.py

Phase 5 -- fourth & fifth scientific questions:
    Q4: Do distinct disagreement patterns correspond to distinct biological
        regimes?
    Q5: Do protein AI models exhibit context-dependent expertise?

Analysis 1 (regimes):
  * Input: per-variant family centroids (z-basis): S_seq, S_evo, S_struct
  * Gaussian Mixture Model (K = 2..8, BIC selection), pre-registered
  * Regime naming follows the actual centroids (protocol s30) -- the template
    names (consensus tolerant / consensus damaging / evolution-dominant /
    structure-dominant / high-disagreement) are used ONLY if the data support
    them.
  * Regime characterization: pLDDT, assay selection type, functional
    annotations, experimental error, experimental effect.

Analysis 2 (context-dependent expertise):
  * Per-variant family prediction error E_fam = |S_fam - Y|
  * Advantage_fam = E_others - E_fam  (positive = family more accurate)
  * Association of family advantage with pLDDT bin, functional annotation,
    assay selection type (protein-level effect aggregation).

Outputs:
  results/statistics/regime_assignment.parquet
  results/statistics/regime_summary.csv
  results/statistics/family_advantage_plddt.csv
  results/statistics/family_advantage_annotation.csv
  figures/main/Fig4A_regime_ternary.png
  figures/main/Fig4B_regime_plddt.png
  figures/main/Fig5A_family_advantage.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    PROCESSED, STATISTICS, FIGURES_MAIN, Job, load_core_panel,
)

SEED = 2026
FAMILIES = ["single_seq", "evolution", "structure"]
PLDDT_BINS = [(0, 50, "pLDDT<50"), (50, 70, "50-70"), (70, 90, "70-90"), (90, 101, ">=90")]


def main():
    job = Job("12_disagreement_regimes")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    feats = pd.read_parquet(PROCESSED / "uniprot_features.parquet")
    ref = pd.read_csv(PROCESSED.parent / "raw" / "reference" / "DMS_substitutions.csv")
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    sel_type = dict(zip(ref["DMS_id"], ref["selection_type"]))
    core = load_core_panel(job)

    df = disc.merge(struct, on=["DMS_id", "mutant"], how="left")
    # percentile-space family centroids (U scale, same [0,1] space as Y)
    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet",
                           columns=["DMS_id", "mutant"] +
                           [f"u_{m}" for m in core["model"]])
    fam_u = {f: [f"u_{m}" for m in core["model"]
                 if core.loc[core.model == m, "family"].iloc[0] == f]
             for f in ["single_seq", "evolution", "structure"]}
    for f, cols in fam_u.items():
        df[f"U_{f}"] = norm[cols].mean(axis=1, skipna=True)
    # also attach z-scale family centroids (used for regime clustering)
    df = df.merge(norm[["DMS_id", "mutant"]], on=["DMS_id", "mutant"], how="left")
    # transfer UniProt features to DMS target coordinates (map from 10)
    pmap = pd.read_parquet(PROCESSED / "uniprot_position_map.parquet")
    pm = pmap[pmap["mapped"]].copy()
    feats_t = feats.merge(pm, left_on=["UniProt_ID", "position"],
                          right_on=["UniProt_ID", "uniprot_pos"], how="inner")
    feats_t["position"] = feats_t["target_pos"].astype(int)
    feats_t = feats_t[["UniProt_ID", "position", "feature"]]
    # pivot per-position features into boolean columns BEFORE merging
    # (pivot_table is unreliable on this pandas build; get_dummies + max)
    d = pd.get_dummies(feats_t.set_index(["UniProt_ID", "position"])["feature"])
    feats_pivot = d.groupby(level=[0, 1]).max().astype(bool).reset_index()
    df = df.merge(feats_pivot, on=["UniProt_ID", "position"], how="left")
    for f in feats["feature"].unique():
        if f not in df:
            df[f] = False
        df[f] = df[f].fillna(False)
    df["selection_type"] = df["DMS_id"].map(sel_type)
    job.info(f"rows: {len(df)}")

    S_cols = [f"S_{f}" for f in FAMILIES]
    X = df[S_cols].to_numpy(dtype=float)
    keep = ~np.isnan(X).any(axis=1)
    job.info(f"variants with all three family centroids: {keep.sum()}")

    # ---------------- Analysis 1: GMM regimes ----------------
    Xk = X[keep]
    bics = []
    models = {}
    for k in range(2, 9):
        m = GaussianMixture(n_components=k, random_state=SEED, covariance_type="full",
                            max_iter=300, n_init=3)
        m.fit(Xk)
        bics.append(m.bic(Xk))
        models[k] = m
    best_k = range(2, 9)[int(np.argmin(bics))]
    job.info(f"BIC per K: " + ", ".join(f"{k}:{b:.0f}" for k, b in zip(range(2, 9), bics)))
    job.info(f"best K = {best_k}")

    gmm = models[best_k]
    labels = gmm.predict(Xk)
    assign = pd.DataFrame({"DMS_id": df.loc[keep, "DMS_id"].values,
                           "mutant": df.loc[keep, "mutant"].values,
                           "regime": labels})
    # map regime -> template name based on actual centroids
    cent = pd.DataFrame(gmm.means_, columns=S_cols)
    cent["size"] = np.bincount(labels, minlength=best_k)
    cent["name"] = ""
    for i in range(best_k):
        row = cent.iloc[i]
        s = {f: row[f"S_{f}"] for f in FAMILIES}
        if max(s.values()) < -0.3:
            cent.loc[i, "name"] = "consensus_tolerant"
        elif min(s.values()) > 0.3:
            cent.loc[i, "name"] = "consensus_damaging"
        elif s["evolution"] > 0.3 and s["structure"] < -0.3:
            cent.loc[i, "name"] = "evolution_dominant"
        elif s["structure"] > 0.3 and s["evolution"] < -0.3:
            cent.loc[i, "name"] = "structure_dominant"
        else:
            spread = max(s.values()) - min(s.values())
            if spread > 1.5:
                cent.loc[i, "name"] = "high_disagreement"
            else:
                cent.loc[i, "name"] = f"regime_{i}"
    cent.to_csv(STATISTICS / "regime_centroids.csv")
    job.info("regime centroids:\n" + cent.round(2).to_string())

    name_map = dict(zip(cent.index, cent["name"]))
    assign["regime_name"] = assign["regime"].map(name_map)
    assign.to_parquet(STATISTICS / "regime_assignment.parquet", index=False)
    job.info("regime sizes: " + str(assign["regime_name"].value_counts().to_dict()))

    # regime characterization table
    full = df.loc[keep].copy()
    full["regime_name"] = labels
    full["regime_name"] = full["regime_name"].map(name_map)
    summ_rows = []
    for name, g in full.groupby("regime_name"):
        summ_rows.append({
            "regime": name, "n": len(g),
            "median_plddt": g["plddt"].median(),
            "frac_low_plddt": (g["plddt"] < 50).mean(),
            "median_error": g["error"].median() if "error" in g else np.nan,
            "median_Y": g["Y_deleter"].median(),
            "median_D": g["D_std"].median(),
            "frac_active": g.get("active_site", pd.Series(False)).mean(),
            "frac_binding": g.get("binding_site", pd.Series(False)).mean(),
            "frac_domain": g.get("domain", pd.Series(False)).mean(),
            "frac_transmem": g.get("transmembrane", pd.Series(False)).mean(),
        })
    summ = pd.DataFrame(summ_rows)
    summ.to_csv(STATISTICS / "regime_summary.csv", index=False)
    job.info("regime summary:\n" + summ.round(3).to_string(index=False))

    # ---------------- Analysis 2: family advantage ----------------
    # NOTE (review): all error/advantage quantities are computed in the
    # percentile (u) space so that predictions and the experimental outcome
    # Y in [0,1] share the same scale.
    adv_rows = []
    fam_labels = {"U_single_seq": "seq", "U_evolution": "evo", "U_structure": "struct"}
    full["error"] = (full[["U_single_seq", "U_evolution", "U_structure"]]
                     .mean(axis=1) - full["Y_deleter"]).abs()
    for fam, lab in fam_labels.items():
        other = [c for c in fam_labels if c != fam]
        E_fam = (full[fam] - full["Y_deleter"]).abs()
        E_oth = full[other].sub(full["Y_deleter"], axis=0).abs().mean(axis=1)
        full[f"adv_{lab}"] = E_oth - E_fam
    adv_plddt = []
    for fam, lab in fam_labels.items():
        for lo, hi, bname in PLDDT_BINS:
            sel = full[(full["plddt"] >= lo) & (full["plddt"] < hi)]
            if len(sel) < 500:
                continue
            meds = sel.groupby("UniProt_ID")[f"adv_{lab}"].median()
            adv_plddt.append({"family": lab, "plddt_bin": bname, "n": len(sel),
                              "median_adv": meds.median(),
                              "frac_pos": (meds > 0).mean()})
    advp = pd.DataFrame(adv_plddt)
    advp.to_csv(STATISTICS / "family_advantage_plddt.csv", index=False)
    job.info("family advantage by pLDDT:\n" + advp.round(3).to_string(index=False))

    adv_ann = []
    for fam, lab in fam_labels.items():
        for feat in ["active_site", "binding_site", "domain", "transmembrane"]:
            if feat not in full:
                continue
            for ann in [True, False]:
                sel = full[full[feat] == ann]
                if len(sel) < 200:
                    continue
                meds = sel.groupby("UniProt_ID")[f"adv_{lab}"].median()
                adv_ann.append({"family": lab, "annotation": feat,
                                "in_annotation": ann, "n": len(sel),
                                "median_adv": meds.median(),
                                "frac_pos": (meds > 0).mean()})
    adva = pd.DataFrame(adv_ann)
    adva.to_csv(STATISTICS / "family_advantage_annotation.csv", index=False)
    job.info("family advantage by annotation (median adv, frac pos):")
    for _, r in adva.iterrows():
        job.info(f"  {r.family} {r.annotation}={r.in_annotation}: "
                 f"{r.median_adv:.4f} ({r.frac_pos:.2f}, n={r.n})")

    # ---------------- figures ----------------
    # Fig4A: ternary-ish scatter (seq vs evo vs struct via 2D projection)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    cmap = plt.get_cmap("tab10")
    regime_colors = {rn: cmap(i) for i, rn in enumerate(full["regime_name"].unique())}
    for ax, (a, b) in zip(axes, [("S_evolution", "S_single_seq"),
                                 ("S_evolution", "S_structure")]):
        for rn in full["regime_name"].unique():
            sel = full[full["regime_name"] == rn]
            idx = np.random.default_rng(SEED).choice(len(sel), size=min(len(sel), 3000),
                                                     replace=False)
            ax.scatter(sel[a].to_numpy()[idx], sel[b].to_numpy()[idx],
                       s=2, alpha=0.3, label=rn, color=regime_colors[rn])
        ax.axhline(0, color="grey", lw=0.4)
        ax.axvline(0, color="grey", lw=0.4)
        ax.set_xlabel(a)
        ax.set_ylabel(b)
        ax.legend(fontsize=6, markerscale=4)
    fig.suptitle("AI disagreement regimes (GMM, K=%d)" % best_k)
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig4A_regime_ternary.png", dpi=200)
    plt.close(fig)

    # Fig4B: regime x pLDDT
    fig, ax = plt.subplots(figsize=(7, 4.5))
    full["plddt_bin"] = pd.cut(full["plddt"], bins=[0, 50, 70, 90, 101],
                               labels=[l for _, _, l in PLDDT_BINS])
    tab = pd.crosstab(full["regime_name"], full["plddt_bin"], normalize="index")
    tab.plot(kind="bar", stacked=True, ax=ax, colormap="viridis")
    ax.set_ylabel("fraction of variants")
    ax.legend(title="pLDDT", fontsize=7)
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig4B_regime_plddt.png", dpi=200)
    plt.close(fig)

    # Fig5A: family advantage by pLDDT
    fig, ax = plt.subplots(figsize=(8, 4.5))
    width = 0.25
    for i, (fam, lab) in enumerate(fam_labels.items()):
        sub = advp[advp["family"] == lab]
        ax.bar(np.arange(len(sub)) + i * width, sub["median_adv"], width=width,
               label=fam)
    ax.set_xticks(np.arange(len(PLDDT_BINS)) + width)
    ax.set_xticklabels([l for _, _, l in PLDDT_BINS])
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("median family advantage (protein-level)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig5A_family_advantage.png", dpi=200)
    plt.close(fig)
    job.info("figures written: Fig4A, Fig4B, Fig5A")
    job.close()


if __name__ == "__main__":
    main()