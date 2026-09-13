"""28_confound_controls.py

Review-response analysis (item 1): partial controls for the shared-input
confound (structure-conditioned models and pLDDT both derive from the same
AlphaFold structures).

Control 1 (asymmetry): if disagreement were purely a technical artifact of
degraded structural inputs at low pLDDT, then ALL structure-vs-other family
disagreements should rise at low pLDDT. We test whether the
sequence-vs-structure disagreement (which involves the same structure
models) shows the trend. The primary result is evolution-vs-structure
specific; a null seq-vs-structure trend argues against a pure input-quality
artifact.

Control 2 (structure-model accuracy by pLDDT stratum): within each pLDDT
bin, the median per-assay Spearman of each structure-conditioned model with
DMS. If structure-model scores remain informative at low pLDDT, the
disagreement is not simply 'structure models are garbage there'.

Control 3 (coverage): whether structure-model score availability varies
with pLDDT.

Outputs:
  results/statistics/confound_controls.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job, load_core_panel  # noqa: E402

BINS = [(0, 50, "pLDDT<50"), (50, 70, "50-70"), (70, 90, "70-90"), (90, 101, ">=90")]


def main():
    job = Job("28_confound_controls")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet")
    core = load_core_panel(job)
    struct_models = core[core.family == "structure"]["model"].tolist()
    df = disc.merge(struct, on=["DMS_id", "mutant"], how="left")
    job.info(f"rows: {len(df)}; structure models: {struct_models}")

    # ---- Control 1: asymmetry of family-pair trends (residue level) ----
    rows = []
    for uni, g in df.groupby("UniProt_ID"):
        rl = g.groupby("position").agg(
            D_evo_struct=("D_evo_struct", "median"),
            D_seq_struct=("D_seq_struct", "median"),
            plddt=("plddt", "median")).dropna()
        if len(rl) < 50:
            continue
        rows.append({"UniProt_ID": uni,
                     "rho_evo_struct": spearmanr(rl.plddt, rl.D_evo_struct)[0],
                     "rho_seq_struct": spearmanr(rl.plddt, rl.D_seq_struct)[0]})
    c1 = pd.DataFrame(rows)
    job.info(f"Control 1: n_proteins = {len(c1)}")
    job.info(f"  rho_evo_struct median = {c1.rho_evo_struct.median():+.4f} "
             f"(frac neg {(c1.rho_evo_struct < 0).mean():.3f})")
    job.info(f"  rho_seq_struct median = {c1.rho_seq_struct.median():+.4f} "
             f"(frac neg {(c1.rho_seq_struct < 0).mean():.3f})")
    job.info(f"  asymmetry (evo-struct minus seq-struct) median = "
             f"{(c1.rho_evo_struct - c1.rho_seq_struct).median():+.4f}")

    # ---- Control 2: structure-model DMS accuracy by pLDDT bin ----
    u_cols = [f"u_{m}" for m in struct_models]
    df2 = df.merge(norm[["DMS_id", "mutant"] + u_cols], on=["DMS_id", "mutant"],
                   how="left")
    acc_rows = []
    for lo, hi, lab in BINS:
        sel = df2[(df2.plddt >= lo) & (df2.plddt < hi)]
        if len(sel) < 500:
            continue
        for m in struct_models:
            rhos = []
            for dms_id, g in sel.groupby("DMS_id"):
                if len(g) >= 50:
                    rhos.append(spearmanr(g[f"u_{m}"], g["DMS_score"])[0])
            if rhos:
                acc_rows.append({"plddt_bin": lab, "model": m,
                                 "n_variants": len(sel),
                                 "n_assays": len(rhos),
                                 "median_rho_dms": np.median(rhos)})
    c2 = pd.DataFrame(acc_rows)
    job.info("Control 2: structure-model DMS correlation by pLDDT bin")
    for _, r in c2.iterrows():
        job.info(f"  {r.model} @ {r.plddt_bin}: median rho(DMS) = {r.median_rho_dms:+.4f}")

    # ---- Control 3: coverage by pLDDT ----
    cov_rows = []
    for lo, hi, lab in BINS:
        sel = df2[(df2.plddt >= lo) & (df2.plddt < hi)]
        for m in struct_models:
            cov_rows.append({"plddt_bin": lab, "model": m,
                             "frac_scored": sel[f"u_{m}"].notna().mean()})
    c3 = pd.DataFrame(cov_rows)
    job.info("Control 3: structure-model coverage by pLDDT bin")
    for _, r in c3.iterrows():
        job.info(f"  {r.model} @ {r.plddt_bin}: coverage = {r.frac_scored:.4f}")

    # save
    with pd.ExcelWriter(STATISTICS / "confound_controls.xlsx") as w:
        c1.to_excel(w, sheet_name="asymmetry", index=False)
        c2.to_excel(w, sheet_name="accuracy_by_plddt", index=False)
        c3.to_excel(w, sheet_name="coverage_by_plddt", index=False)
    c1.to_csv(STATISTICS / "confound_control_asymmetry.csv", index=False)
    c2.to_csv(STATISTICS / "confound_control_accuracy.csv", index=False)
    c3.to_csv(STATISTICS / "confound_control_coverage.csv", index=False)
    job.close()


if __name__ == "__main__":
    main()