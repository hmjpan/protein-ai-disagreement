"""30_multivariable_regression.py

Review-response analysis (item 8): is the pLDDT-disagreement association
robust to adjustment for correlated residue-level covariates?

Protein-fixed-effects regression (within-protein demeaning, equivalent to
protein dummies) of residue-level D_evo_struct on:
  pLDDT, IDR (UniProt disorder), contact density, secondary structure
  (helix/sheet), normalized position, MSA depth, structural coverage

Report: standardized beta for pLDDT with and without covariates; cluster-
robust standard errors at the protein level (statsmodels).

Outputs:
  results/statistics/multivariable_regression.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job, read_reference_dms  # noqa: E402


def main():
    job = Job("30_multivariable_regression")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt", "structure_cov"])
    mech = pd.read_parquet(PROCESSED / "structure_mechanics.parquet")
    disf = pd.read_parquet(PROCESSED / "disorder_features.parquet",
                           columns=["DMS_id", "mutant", "disorder"])
    ref = read_reference_dms()
    ref["DMS_id"] = ref["DMS_id"].astype(str)

    df = disc.merge(struct, on=["DMS_id", "mutant"], how="left")
    df = df.merge(mech[["UniProt_ID", "position", "ss3", "contact8", "burial"]],
                  on=["UniProt_ID", "position"], how="left")
    df = df.merge(disf[["DMS_id", "mutant", "disorder"]], on=["DMS_id", "mutant"],
                  how="left")
    df["seq_len"] = df["DMS_id"].map(dict(zip(ref["DMS_id"], ref["seq_len"])))
    df["msa_depth"] = df["DMS_id"].map(dict(zip(ref["DMS_id"], ref["MSA_Neff_L_category"])))
    df["norm_pos"] = df["position"] / df["seq_len"]
    df["disorder"] = df["disorder"].fillna(False).astype(int)
    df["is_helix"] = df["ss3"].eq("helix").astype(int)
    df["is_sheet"] = df["ss3"].eq("sheet").astype(int)
    df["msa_low"] = df["msa_depth"].eq("Low").astype(int)
    df["msa_high"] = df["msa_depth"].eq("High").astype(int)

    # residue-level aggregation
    rl = df.groupby(["UniProt_ID", "position"]).agg(
        D_evo_struct=("D_evo_struct", "median"),
        plddt=("plddt", "median"),
        disorder=("disorder", "max"),
        contact8=("contact8", "median"),
        norm_pos=("norm_pos", "median"),
        is_helix=("is_helix", "max"),
        is_sheet=("is_sheet", "max"),
        msa_low=("msa_low", "max"),
        msa_high=("msa_high", "max"),
        structure_cov=("structure_cov", "median")).dropna().reset_index()
    job.info(f"residue-level rows: {len(rl)}; proteins: {rl.UniProt_ID.nunique()}")

    # within-protein demeaning (protein fixed effects)
    for col in ["D_evo_struct", "plddt", "disorder", "contact8", "norm_pos",
                "is_helix", "is_sheet", "msa_low", "msa_high"]:
        rl[col + "_dm"] = rl[col] - rl.groupby("UniProt_ID")[col].transform("mean")

    X0 = sm.add_constant(rl["plddt_dm"])
    Xf = sm.add_constant(rl[["plddt_dm", "disorder_dm", "contact8_dm",
                             "norm_pos_dm", "is_helix_dm", "is_sheet_dm",
                             "msa_low_dm", "msa_high_dm"]])
    y = rl["D_evo_struct_dm"]

    m0 = sm.OLS(y, X0).fit(cov_type="cluster", cov_kwds={"groups": rl["UniProt_ID"]})
    mf = sm.OLS(y, Xf).fit(cov_type="cluster", cov_kwds={"groups": rl["UniProt_ID"]})

    ci0 = m0.conf_int().loc["plddt_dm"]
    cif = mf.conf_int().loc["plddt_dm"]
    rows = [{"model": "pLDDT only",
             "beta_plddt": m0.params["plddt_dm"],
             "ci_low": ci0[0], "ci_high": ci0[1],
             "se_plddt": m0.bse["plddt_dm"],
             "p_plddt": m0.pvalues["plddt_dm"],
             "n": int(m0.nobs),
             "n_proteins": rl.UniProt_ID.nunique(),
             "r2": m0.rsquared}]
    rows.append({"model": "full covariates",
                 "beta_plddt": mf.params["plddt_dm"],
                 "ci_low": cif[0], "ci_high": cif[1],
                 "se_plddt": mf.bse["plddt_dm"],
                 "p_plddt": mf.pvalues["plddt_dm"],
                 "beta_disorder": mf.params["disorder_dm"],
                 "p_disorder": mf.pvalues["disorder_dm"],
                 "beta_contact": mf.params["contact8_dm"],
                 "beta_helix": mf.params["is_helix_dm"],
                 "beta_sheet": mf.params["is_sheet_dm"],
                 "n": int(mf.nobs),
                 "n_proteins": rl.UniProt_ID.nunique(),
                 "r2": mf.rsquared})
    out = pd.DataFrame(rows)
    out.to_csv(STATISTICS / "multivariable_regression.csv", index=False)
    for _, r in out.iterrows():
        job.info(f"{r.model}: beta_pLDDT = {r.beta_plddt:+.5f} "
                 f"(SE {r.se_plddt:.5f}, p = {r.p_plddt:.2e}, r2 = {r.r2:.4f})")
        if "beta_disorder" in r and not np.isnan(r.beta_disorder):
            job.info(f"  beta_disorder = {r.beta_disorder:+.5f} "
                     f"(p = {r.p_disorder:.2e}); "
                     f"beta_contact = {r.beta_contact:+.5f}; "
                     f"beta_helix = {r.beta_helix:+.5f}; "
                     f"beta_sheet = {r.beta_sheet:+.5f}")
    job.close()


if __name__ == "__main__":
    main()