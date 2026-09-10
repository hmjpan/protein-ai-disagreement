"""27_panel_residue_level.py

Verify that the panel-composition claims in Section 2.2 hold at the
residue level (the primary analysis scale):

  A. strict 8-model architecture panel: per-protein residue-level
     Spearman(pLDDT, D_evo_struct)
  B. 1,000 balanced resamplings (3 models per family; structure fixed by
     availability): per-protein residue-level Spearman, distribution of
     medians

Outputs:
  results/statistics/panel_residue_level_A.csv
  results/statistics/panel_residue_level_B_medians.npy
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job, load_core_panel  # noqa: E402

SEED = 2026
ARCH_PANEL = ["ESM1v_ensemble", "GEMME", "EVE_ensemble", "MSA_Transformer_ensemble",
              "ESM-IF1", "ProteinMPNN", "SaProt_650M_AF2", "Tranception_L"]
N_BALANCED = 1000
N_PER_FAMILY = 3
FAMILIES = ["evolution", "single_seq", "structure"]
MIN_RESIDUES = 50


def residue_trend_evo_struct(zmat, plddt, pos, groups):
    """Per-protein residue-level Spearman(pLDDT, |S_evo - S_struct|).

    Mutations at the same position are aggregated to the residue median;
    pLDDT is per-position.
    """
    n_mod = zmat.shape[1]
    evo = [i for i in range(n_mod) if fam_of_col[i] == "evolution"]
    struct = [i for i in range(n_mod) if fam_of_col[i] == "structure"]
    if not evo or not struct:
        return np.array([])
    out = []
    for g in np.unique(groups):
        m = groups == g
        if m.sum() < 100:
            continue
        d = np.abs(np.nanmean(zmat[m][:, evo], axis=1) -
                   np.nanmean(zmat[m][:, struct], axis=1))
        p = plddt[m]
        ps = pos[m]
        ok = ~np.isnan(d) & ~np.isnan(p) & (ps > 0)
        if ok.sum() < MIN_RESIDUES:
            continue
        rl = pd.DataFrame({"d": d[ok], "p": p[ok], "pos": ps[ok]}).groupby(
            "pos").agg(d=("d", "median"), p=("p", "median"))
        if len(rl) >= MIN_RESIDUES:
            out.append(spearmanr(rl["p"], rl["d"])[0])
    return np.array(out)


def main():
    global fam_of_col
    job = Job("27_panel_residue_level")
    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    core = load_core_panel(job)
    fam_map = dict(zip(core["model"], core["family"]))
    df = norm[["DMS_id", "UniProt_ID", "mutant", "position"]].merge(
        struct, on=["DMS_id", "mutant"], how="left")
    groups = df["UniProt_ID"].to_numpy()
    plddt = df["plddt"].to_numpy()
    pos = df["position"].to_numpy(dtype=float)
    job.info(f"rows: {len(df)}; proteins: {len(np.unique(groups))}")

    # ---- Analysis A: strict 8-model architecture panel ----
    zcols = [f"z_{m}" for m in ARCH_PANEL if f"z_{m}" in norm.columns]
    Z = norm[zcols].to_numpy(dtype=float)
    fam_of_col = [fam_map.get(c[2:], "other") for c in zcols]
    job.info(f"Analysis A panel: {zcols}")
    rhoA = residue_trend_evo_struct(Z, plddt, pos, groups)
    a_sum = pd.DataFrame({"median_rho": np.median(rhoA),
                          "n_proteins": len(rhoA),
                          "frac_negative": (rhoA < 0).mean()}, index=[0])
    a_sum.to_csv(STATISTICS / "panel_residue_level_A.csv", index=False)
    job.info(f"Analysis A (residue-level): median rho = {np.median(rhoA):.4f} "
             f"(n={len(rhoA)}, frac neg {(rhoA < 0).mean():.3f})")

    # ---- Analysis B: balanced resampling, residue level ----
    fam_models = {f: [m for m in core["model"] if fam_map[m] == f]
                  for f in FAMILIES}
    rng = np.random.default_rng(SEED)
    medians = []
    t0 = time.time()
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
        Zb = norm[cols].to_numpy(dtype=float)
        fam_of_col = [fam_map.get(c[2:], "other") for c in cols]
        rho = residue_trend_evo_struct(Zb, plddt, pos, groups)
        if len(rho):
            medians.append(np.median(rho))
        if (it + 1) % 200 == 0:
            el = time.time() - t0
            job.info(f"balanced {it + 1}/{N_BALANCED}: "
                     f"median of medians = {np.median(medians):.4f} "
                     f"({el:.0f}s elapsed)")
    med = np.array(medians)
    np.save(STATISTICS / "panel_residue_level_B_medians.npy", med)
    job.info(f"Analysis B (residue-level): median = {np.median(med):.4f} "
             f"[{np.percentile(med, 2.5):.4f}, {np.percentile(med, 97.5):.4f}], "
             f"frac negative = {(med < 0).mean():.3f}")
    job.close()


if __name__ == "__main__":
    main()