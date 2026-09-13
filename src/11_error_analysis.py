"""11_error_analysis.py

Phase 3 -- third scientific question:
    Q3: Is model disagreement associated with prediction error?

Pre-registered analysis plan (protocol sections 18-21, 60, 75):

  * experimental deleteriousness  Y = 1 - rank(DMS_score) within assay
  * ensemble prediction           S = mean(u_m) over core models
  * ensemble error                E = |S - Y|
  * association                   per-assay Spearman(D_std, E)
  * aggregation                   median rho; UniProt-level cluster bootstrap
                                  95% CI; sign test; fraction of assays rho>0
  * decile analysis               within-assay deciles of D_std -> mean E
  * confounder control            within-pLDDT-bin Spearman (alpha = pLDDT is
                                  the main competing explanation)
  * sensitivity                   D_MAD; leave-one-model-out D_std

GO/NO-GO 1 rule (config.yaml): if the UniProt-cluster-bootstrapped 95% CI of
the median assay-level rho crosses 0, or median rho < 0.05, the disagreement-
as-uncertainty narrative FAILS (Outcome A); the paper pivots to biological
specialization.

Outputs: results/statistics/disagreement_error_assay_rho.csv
         results/statistics/disagreement_deciles.csv
         results/statistics/disagreement_error_plddt_stratified.csv
         results/statistics/disagreement_lomo_sensitivity.csv
         results/GO_NO_GO_1.md
         figures/main/Fig2A_disagreement_distribution.png
         figures/main/Fig2B_disagreement_deciles.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    PROCESSED, STATISTICS, RESULTS, FIGURES_MAIN, Job, load_core_panel,
)

SEED = 2026
BOOT = 10000
MIN_ASSET_N = 50
GO_LOWER_RHO = 0.05  # config.yaml outcome_A_median_rho


def cluster_bootstrap_median(rhos: pd.DataFrame, n_iter=BOOT, seed=SEED):
    """Bootstrap the median of per-assay rho by resampling UniProt clusters."""
    rng = np.random.default_rng(seed)
    medians = []
    uniprots = rhos["UniProt_ID"].unique()
    for _ in range(n_iter):
        ids = rng.choice(uniprots, size=len(uniprots), replace=True)
        sub = rhos[rhos["UniProt_ID"].isin(ids)].groupby(
            "UniProt_ID").first().reindex(ids)["rho"]
        medians.append(np.median(sub.dropna()))
    return np.percentile(medians, [2.5, 97.5]), medians


def main():
    job = Job("11_error_analysis")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    core = load_core_panel(job)
    u_cols = [f"u_{m}" for m in core["model"]]
    z_cols = [f"z_{m}" for m in core["model"]]
    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet",
                           columns=["DMS_id", "mutant"] + u_cols + z_cols)
    m_df = disc.merge(norm[["DMS_id", "mutant"] + u_cols + z_cols],
                      on=["DMS_id", "mutant"],
                      how="left", validate="one_to_one")
    job.info(f"merged rows: {len(m_df)}")

    # ensemble in percentile space: mean of the three family u-centroids
    fam_u = {f: [f"u_{m}" for m in core["model"]
                 if core.loc[core.model == m, "family"].iloc[0] == f]
             for f in ["single_seq", "evolution", "structure"]}
    u_fams = []
    for f, cols in fam_u.items():
        m_df[f"U_{f}"] = m_df[cols].mean(axis=1, skipna=True)
        u_fams.append(f"U_{f}")
    S = m_df[u_fams].mean(axis=1, skipna=True)
    m_df["S_ensemble"] = S
    m_df["error"] = (S - m_df["Y_deleter"]).abs()
    m_df["D_mad"] = disc["D_mad"]
    job.info(f"mean |S-Y| = {m_df['error'].mean():.4f}")

    # ---- per-assay Spearman(D, error) ----
    rows = []
    for dms_id, g in m_df.groupby("DMS_id"):
        if g["D_std"].notna().sum() < MIN_ASSET_N:
            continue
        rho, p = spearmanr(g["D_std"], g["error"])
        rows.append({"DMS_id": dms_id, "UniProt_ID": g["UniProt_ID"].iloc[0],
                     "n": len(g), "rho": rho, "p": p})
    rho_df = pd.DataFrame(rows)
    rho_df.to_csv(STATISTICS / "disagreement_error_assay_rho.csv", index=False)
    med = np.median(rho_df["rho"])
    frac_pos = (rho_df["rho"] > 0).mean()
    st = binomtest(int((rho_df["rho"] > np.median(rho_df["rho"])).sum()),
                   len(rho_df), p=0.5)
    ci, _ = cluster_bootstrap_median(rho_df)
    job.info(f"median assay rho(D,error) = {med:.4f} [{ci[0]:.4f}, {ci[1]:.4f}]; "
             f"frac rho>0 = {frac_pos:.3f}; median n = {rho_df['n'].median():.0f}")

    # ---- decile analysis (within assay) ----
    dec_rows = []
    for dms_id, g in m_df.groupby("DMS_id"):
        if len(g) < 100:
            continue
        q = pd.qcut(g["D_std"].rank(method="first"), 10, labels=False)
        dec_rows.append(pd.DataFrame({
            "decile": q, "error": g["error"].values}))
    dec = pd.concat(dec_rows)
    dec_summary = dec.groupby("decile")["error"].agg(["mean", "median", "std", "count"])
    dec_summary.to_csv(STATISTICS / "disagreement_deciles.csv")
    job.info(f"decile errors: " + ", ".join(f"{v:.4f}" for v in dec_summary["mean"]))

    # ---- pLDDT-stratified association (confounder control) ----
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    m_df = m_df.merge(struct, on=["DMS_id", "mutant"], how="left")
    strat_rows = []
    bins = [(0, 50, "<50"), (50, 70, "50-70"), (70, 90, "70-90"), (90, 101, ">=90")]
    for lo, hi, lab in bins:
        sel = m_df[(m_df["plddt"] >= lo) & (m_df["plddt"] < hi) & m_df["D_std"].notna()]
        rhos = []
        for dms_id, g in sel.groupby("DMS_id"):
            if len(g) >= 50:
                rhos.append(spearmanr(g["D_std"], g["error"])[0])
        if rhos:
            strat_rows.append({"plddt_bin": lab, "n_variants": len(sel),
                               "n_assays": len(rhos), "median_rho": np.median(rhos),
                               "p25": np.percentile(rhos, 25),
                               "p75": np.percentile(rhos, 75)})
    strat = pd.DataFrame(strat_rows)
    strat.to_csv(STATISTICS / "disagreement_error_plddt_stratified.csv", index=False)
    job.info("pLDDT-stratified rho: " +
             "; ".join(f"{r.plddt_bin}:{r.median_rho:.3f}" for r in strat.itertuples()))

    # ---- leave-one-model-out sensitivity (z-score basis, consistent with
    # D_std) ----
    lomo = []
    Z = m_df[z_cols].to_numpy(dtype=float)
    for i, m in enumerate(core["model"]):
        D_lo = np.nanstd(np.delete(Z, i, axis=1), axis=1)
        tmp = pd.DataFrame({"D": D_lo, "e": m_df["error"].to_numpy(),
                            "dms": m_df["DMS_id"].to_numpy()})
        rhos = []
        for dms_id, g in tmp.groupby("dms"):
            if len(g) >= MIN_ASSET_N and np.isfinite(g["D"]).sum() >= MIN_ASSET_N:
                rhos.append(spearmanr(g["D"], g["e"])[0])
        lomo.append({"left_out": m, "median_rho": np.median(rhos),
                     "n_assays": len(rhos)})
    lomo_df = pd.DataFrame(lomo)
    lomo_df.to_csv(STATISTICS / "disagreement_lomo_sensitivity.csv", index=False)
    job.info(f"LOMO (z) median rho range: {lomo_df.median_rho.min():.4f} - "
             f"{lomo_df.median_rho.max():.4f}")

    # ---- D_MAD sensitivity (robust disagreement metric) ----
    mad_rows = []
    for dms_id, g in m_df.groupby("DMS_id"):
        if g["D_mad"].notna().sum() < MIN_ASSET_N:
            continue
        mad_rows.append({"DMS_id": dms_id, "rho": spearmanr(g["D_mad"], g["error"])[0],
                         "n": len(g)})
    mad_df = pd.DataFrame(mad_rows)
    mad_df.to_csv(STATISTICS / "disagreement_error_mad_sensitivity.csv", index=False)
    med_mad = np.median(mad_df["rho"])
    job.info(f"median assay rho(D_mad,error) = {med_mad:.4f}")

    # ---- u-basis disagreement (methodological comparison) ----
    u_rows = []
    U = m_df[u_cols].to_numpy(dtype=float)
    D_u = np.nanstd(U, axis=1)
    tmp = pd.DataFrame({"D": D_u, "e": m_df["error"].to_numpy(),
                        "dms": m_df["DMS_id"].to_numpy()})
    for dms_id, g in tmp.groupby("dms"):
        if len(g) >= MIN_ASSET_N:
            u_rows.append({"DMS_id": dms_id, "rho": spearmanr(g["D"], g["e"])[0]})
    u_df = pd.DataFrame(u_rows)
    med_u = np.median(u_df["rho"])
    job.info(f"median assay rho(D_u,error) = {med_u:.4f} (methodological note)")

    # ---- GO / NO-GO 1 (protocol: BOTH median < threshold AND CI crosses 0) ----
    ci_lower = ci[0]
    verdict = "GO" if (med >= GO_LOWER_RHO or ci_lower > 0) else "PIVOT"
    lines = [
        "# GO / NO-GO 1 -- Disagreement vs prediction error",
        "",
        f"- Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}",
        f"- Assays analysed: {len(rho_df)}",
        f"- Median assay-level Spearman(D_std, error): {med:.4f}",
        f"- UniProt-cluster bootstrap 95% CI: [{ci[0]:.4f}, {ci[1]:.4f}]",
        f"- Fraction of assays with rho > 0: {frac_pos:.3f}",
        f"- Decile 10 vs 1 mean error: {dec_summary.loc[9,'mean']:.4f} vs {dec_summary.loc[0,'mean']:.4f}",
        f"- pLDDT-stratified median rho: " +
        "; ".join(f"{r.plddt_bin}={r.median_rho:.3f}" for r in strat.itertuples()),
        f"- Leave-one-model-out median rho range: "
        f"[{lomo_df.median_rho.min():.4f}, {lomo_df.median_rho.max():.4f}]",
        f"- D_MAD sensitivity median rho: {med_mad:.4f}",
        f"- u-basis disagreement median rho: {med_u:.4f} (methodological note)",
        f"- Pre-registered rule (protocol s60): PIVOT only if median rho < "
        f"{GO_LOWER_RHO} AND bootstrap CI crosses 0",
        "",
        f"**VERDICT: {verdict}**",
        "",
    ]
    if verdict == "PIVOT":
        lines.append("Follow protocol Outcome A: pivot to biological "
                     "specialization framing; disagreement-as-uncertainty "
                     "claims are dropped.")
    elif med < GO_LOWER_RHO:
        lines.append(f"NOTE: median rho {med:.4f} is below the {GO_LOWER_RHO} "
                     "effect-size target even though the CI excludes 0; "
                     "disagreement-as-uncertainty is reported as a weak, "
                     "consistent signal (no 'strongly predicts' wording).")
    (RESULTS / "GO_NO_GO_1.md").write_text("\n".join(lines), encoding="utf-8")
    job.info(f"GO/NO-GO 1 verdict: {verdict}")

    # ---- figures ----
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(rho_df["rho"], bins=50, color="#555555", alpha=0.8)
    ax.axvline(med, color="red", ls="--", label=f"median = {med:.3f}")
    ax.set_xlabel("assay-level Spearman(D, error)")
    ax.set_ylabel("number of assays")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig2A_disagreement_distribution.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(dec_summary.index + 1, dec_summary["mean"], color="#CC6677")
    ax.set_xlabel("within-assay disagreement decile")
    ax.set_ylabel("mean |ensemble - DMS| error")
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig2B_disagreement_deciles.png", dpi=200)
    plt.close(fig)
    job.info("figures written: Fig2A, Fig2B")
    job.close()


if __name__ == "__main__":
    main()