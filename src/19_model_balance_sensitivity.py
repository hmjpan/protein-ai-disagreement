"""19_model_balance_sensitivity.py

Enhancement analysis (review v2, item 4): are the family-disagreement
results robust to panel composition?

Analysis A -- strict one-per-architecture panel:
  {ESM1v_ensemble, GEMME, EVE_ensemble, MSA_Transformer_ensemble, ESM-IF1,
   ProteinMPNN, SaProt_650M_AF2, Tranception_L}
  -> recompute per-protein Spearman(pLDDT-quartile, D_evo_struct) and
     D_std error-association.

Analysis B -- balanced family sampling:
  1000 iterations; in each, draw 3 models per family from
  {evolution, single_seq, structure}; recompute family centroids and the
  per-protein pLDDT trend of D_evo_struct. Report the distribution of
  median within-protein Spearman.

Outputs:
  results/statistics/panel_sensitivity_analysisA.csv
  results/statistics/panel_sensitivity_analysisB.csv
  results/statistics/panel_sensitivity_summary.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job, load_core_panel  # noqa: E402

SEED = 2026
ARCH_PANEL = ["ESM1v_ensemble", "GEMME", "EVE_ensemble", "MSA_Transformer_ensemble",
              "ESM-IF1", "ProteinMPNN", "SaProt_650M_AF2", "Tranception_L"]
N_BALANCED = 1000
N_PER_FAMILY = 3
FAMILIES = ["evolution", "single_seq", "structure"]


def trend_of_pair(zcols: list, fam_evo: list, fam_struct: list, Z: np.ndarray,
                  plddt: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Per-group Spearman between pLDDT-quartile median and |S_evo - S_struct|."""
    idx_evo = [i for i, c in enumerate(zcols) if c in fam_evo]
    idx_struct = [i for i, c in enumerate(zcols) if c in fam_struct]
    if not idx_evo or not idx_struct:
        return np.array([])
    out = []
    for g in np.unique(groups):
        m = groups == g
        if m.sum() < 100:
            continue
        d = np.abs(np.nanmean(Z[m][:, idx_evo], axis=1) -
                   np.nanmean(Z[m][:, idx_struct], axis=1))
        p = plddt[m]
        ok = ~np.isnan(d) & ~np.isnan(p)
        if ok.sum() < 50:
            continue
        q = pd.qcut(pd.Series(p[ok]).rank(method="first"), 4, labels=False)
        t = pd.DataFrame({"q": q, "p": p[ok], "d": d[ok]})
        med = t.groupby("q").agg(mp=("p", "median"), md=("d", "median"))
        if len(med) >= 3:
            out.append(spearmanr(med["mp"], med["md"])[0])
    return np.array(out)


def main():
    job = Job("19_model_balance_sensitivity")
    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    core = load_core_panel(job)
    fam_map = dict(zip(core["model"], core["family"]))
    df = norm[["DMS_id", "UniProt_ID", "mutant"]].merge(
        struct, on=["DMS_id", "mutant"], how="left")
    groups = df["UniProt_ID"].to_numpy()
    plddt = df["plddt"].to_numpy()
    job.info(f"rows: {len(df)}; proteins: {len(np.unique(groups))}")

    # ---- Analysis A: strict one-per-architecture panel ----
    zcols = [f"z_{m}" for m in ARCH_PANEL if f"z_{m}" in norm.columns]
    Z = norm[zcols].to_numpy(dtype=float)
    fam_evo_a = [f"z_{m}" for m in ARCH_PANEL if fam_map.get(m) == "evolution"]
    fam_struct_a = [f"z_{m}" for m in ARCH_PANEL if fam_map.get(m) == "structure"]
    job.info(f"Analysis A panel: {len(zcols)} models: {[c[2:] for c in zcols]}")
    rhoA = trend_of_pair(zcols, fam_evo_a, fam_struct_a, Z, plddt, groups)
    a_sum = pd.DataFrame({
        "analysis": "A_arch_panel",
        "median_rho_D_evo_struct": [float(np.median(rhoA))],
        "n_proteins": [int(len(rhoA))],
        "frac_negative": [float((rhoA < 0).mean())],
    })
    a_sum.to_csv(STATISTICS / "panel_sensitivity_analysisA.csv", index=False)
    job.info(f"Analysis A (D_evo_struct): median within-protein rho = "
             f"{np.median(rhoA):.4f} (n={len(rhoA)} proteins, "
             f"frac negative { (rhoA < 0).mean():.3f})")

    # ---- Analysis B: balanced family sampling ----
    fam_models = {f: [m for m in core["model"] if fam_map[m] == f]
                  for f in FAMILIES}
    job.info(f"family sizes: { {f: len(fam_models[f]) for f in FAMILIES} }")
    if any(len(v) < N_PER_FAMILY for v in fam_models.values()):
        job.info("cannot balance: a family has < 3 members (structure has 3) "
                 "-- structure sampled with replacement")
    rng = np.random.default_rng(SEED)
    medians = []
    frac_neg = []
    for it in range(N_BALANCED):
        chosen = []
        for f in FAMILIES:
            pool = fam_models[f]
            if len(pool) <= N_PER_FAMILY:
                chosen.extend(pool)
            else:
                chosen.extend(rng.choice(pool, size=N_PER_FAMILY, replace=False))
        cols = [f"z_{m}" for m in chosen if f"z_{m}" in norm.columns]
        if len(cols) < 9:
            continue
        # family centroids
        fam_cols = {}
        for f in FAMILIES:
            fc = [c for c in cols if fam_map[c[2:]] == f]
            if fc:
                fam_cols[f] = fc
        if len(fam_cols) < 3:
            continue
        Zb = norm[cols].to_numpy(dtype=float)
        S_evo = np.nanmean(Zb[:, [cols.index(c) for c in fam_cols["evolution"]]], axis=1)
        S_str = np.nanmean(Zb[:, [cols.index(c) for c in fam_cols["structure"]]], axis=1)
        De = np.abs(S_evo - S_str)
        # per-protein pLDDT-quartile trend of De
        rho = []
        for g in np.unique(groups):
            m = groups == g
            if m.sum() < 100:
                continue
            p = plddt[m]
            ok = ~np.isnan(De[m]) & ~np.isnan(p)
            if ok.sum() < 50:
                continue
            q = pd.qcut(pd.Series(p[ok]).rank(method="first"), 4, labels=False)
            t = pd.DataFrame({"q": q, "p": p[ok], "d": De[m][ok]})
            med = t.groupby("q").agg(mp=("p", "median"), md=("d", "median"))
            if len(med) >= 3:
                rho.append(spearmanr(med["mp"], med["md"])[0])
        if rho:
            medians.append(np.median(rho))
            frac_neg.append((np.array(rho) < 0).mean())
        if (it + 1) % 200 == 0:
            job.info(f"balanced sampling {it + 1}/{N_BALANCED}: "
                     f"current median of medians = {np.median(medians):.4f}")
    med = np.array(medians)
    b_sum = pd.DataFrame({
        "analysis": "B_balanced_sampling",
        "n_iterations": [len(med)],
        "median_of_medians": [float(np.median(med))],
        "p2.5": [float(np.percentile(med, 2.5))],
        "p97.5": [float(np.percentile(med, 97.5))],
        "frac_iterations_negative_median": [float((med < 0).mean())],
    })
    b_sum.to_csv(STATISTICS / "panel_sensitivity_analysisB.csv", index=False)
    job.info(f"Analysis B: median of per-iteration medians = {np.median(med):.4f} "
             f"[{np.percentile(med, 2.5):.4f}, {np.percentile(med, 97.5):.4f}]; "
             f"{(med < 0).mean():.3f} of iterations negative")

    np.save(STATISTICS / "panel_sensitivity_B_medians.npy", med)
    pd.concat([a_sum, b_sum], ignore_index=True).to_csv(
        STATISTICS / "panel_sensitivity_summary.csv", index=False)
    job.close()


if __name__ == "__main__":
    main()