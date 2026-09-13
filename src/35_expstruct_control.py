"""35_expstruct_control.py

Reviewer-response (item 1): quantifying the shared-input confound with
independent data from gitter-lab (Sharma & Gitter, ICLR-W 2025,
Zenodo 10.5281/zenodo.13821399).

Part A (assay-level): for assays with BOTH ESM-IF1 with experimental
structures (spearman) and with AlphaFold2 structures
(predicted_struct_score), gain G = exp - af2 measures the input-quality
penalty. Relate G to assay mean pLDDT (our mapping).

Part B (cross-architecture, variant-level): SSEmb (equivariant GNN,
216 assays, different architecture & training pipeline than ESM-IF1/
ProteinMPNN/MIF). Recompute residue-level Spearman(pLDDT, |U_evo -
U_ssemb|) and compare to the AF2 inverse-folding result.

Outputs:
  results/statistics/expstruct_input_gain.csv
  results/statistics/ssemb_residue_replication.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, RAW, STATISTICS, Job, load_core_panel  # noqa: E402

EXP_CSV = RAW / "exp_struct" / "results" / "experimental_struct_scores.csv"
SSEMB_CSV = RAW / "exp_struct" / "results" / "df_total_proteingym_1.csv"
MIN_RESIDUES = 50


def main():
    job = Job("35_expstruct_control")

    # ---------- Part A: assay-level input-quality gain ----------
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    d = disc.merge(struct, on=["DMS_id", "mutant"], how="left")
    plddt_by_assay = d.groupby("DMS_id")["plddt"].median()

    es = pd.read_csv(EXP_CSV)
    es["mean_plddt"] = es["DMS_id"].map(plddt_by_assay)
    es = es.dropna(subset=["spearman", "predicted_struct_score", "mean_plddt"])
    es["gain"] = es["spearman"] - es["predicted_struct_score"]
    es.to_csv(STATISTICS / "expstruct_input_gain.csv", index=False)
    rho_g, p_g = spearmanr(es["mean_plddt"], es["gain"])
    low = es[es.mean_plddt < es.mean_plddt.median()]
    high = es[es.mean_plddt >= es.mean_plddt.median()]
    job.info(f"Part A: n={len(es)} assays with exp-structure scores")
    job.info(f"  mean gain (exp - AF2) = {es.gain.mean():+.4f}")
    job.info(f"  Spearman(mean_pLDDT, gain) = {rho_g:+.3f} (p={p_g:.3g})")
    job.info(f"  gain in low-confidence assays: {low.gain.mean():+.4f}; "
             f"high-confidence: {high.gain.mean():+.4f}")

    # ---------- Part B: SSEmb cross-architecture replication ----------
    core = load_core_panel(job)
    evo_models = core[core.family == "evolution"]["model"].tolist()
    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet",
                           columns=["DMS_id", "mutant"] + [f"u_{m}" for m in evo_models])
    u_evo = norm[[f"u_{m}" for m in evo_models]].mean(axis=1, skipna=True)
    norm = norm[["DMS_id", "mutant"]].assign(U_evolution=u_evo)

    pos = d[["DMS_id", "mutant", "position", "UniProt_ID", "plddt"]].merge(
        norm, on=["DMS_id", "mutant"], how="left")

    sm = pd.read_csv(SSEMB_CSV)
    sm = sm.rename(columns={"dms_id": "DMS_id", "variant_set": "mutant",
                            "score_ml": "ssemb"})
    sm = sm.dropna(subset=["ssemb"])
    # per-assay rank normalize SSEmb into u (deleterosity): released
    # convention check via DMS correlation sign
    merged = pos.merge(sm[["DMS_id", "mutant", "ssemb"]],
                       on=["DMS_id", "mutant"], how="inner")
    job.info(f"SSEmb matched variants: {len(merged)}; "
             f"proteins: {merged.UniProt_ID.nunique()}")
    # per-assay correlation of raw ssemb vs DMS_score determines the
    # deleterosity orientation (sanity check; released scores may keep
    # either sign)
    signs = []
    for dms, g in merged.groupby("DMS_id"):
        yy = d[(d.DMS_id == dms)].set_index("mutant")["DMS_score"]
        common = g["mutant"].isin(yy.index)
        if common.sum() > 50:
            signs.append(spearmanr(g.loc[common, "ssemb"],
                                   yy.loc[g.loc[common, "mutant"]])[0])
    sign = 1 if np.median(signs) > 0 else -1
    job.info(f"SSEmb median rho with DMS (raw) = {np.median(signs):+.3f} "
             f"-> deleterosity sign = {sign}")
    merged["ssemb_signed"] = sign * merged["ssemb"]
    u_s = []
    for dms, g in merged.groupby("DMS_id"):
        r = rankdata(g["ssemb_signed"].to_numpy()) / len(g)
        u_s.append(pd.DataFrame({"DMS_id": dms, "mutant": g["mutant"].to_numpy(),
                                 "U_ssemb": 1 - r}))
    usdf = pd.concat(u_s)
    merged = merged.merge(usdf, on=["DMS_id", "mutant"])
    merged["D_evo_ssemb"] = (merged["U_evolution"] - merged["U_ssemb"]).abs()

    rows = []
    for uni, g in merged.groupby("UniProt_ID"):
        rl = g.groupby("position").agg(
            D=("D_evo_ssemb", "median"), p=("plddt", "median")).dropna()
        if len(rl) < MIN_RESIDUES:
            continue
        rows.append({"UniProt_ID": uni, "n": len(rl),
                     "rho": spearmanr(rl["p"], rl["D"])[0]})
    rep = pd.DataFrame(rows)
    rep.to_csv(STATISTICS / "ssemb_residue_replication.csv", index=False)
    v = rep["rho"].dropna()
    job.info(f"Part B SSEmb: n={len(v)} proteins; median rho(pLDDT, D_evo-SSEmb) "
             f"= {v.median():+.4f}; frac negative {(v < 0).mean():.3f}")
    # compare with AF2 inverse-folding primary (-0.042, 67.5%)
    job.close()


if __name__ == "__main__":
    main()