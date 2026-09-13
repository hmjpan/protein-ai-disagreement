"""31_xgboost_ablation.py

Review-response analysis (item 12): does context add information to
XGBoost, or is the earlier XGBoost advantage simply nonlinearity?

5-fold UniProt-grouped CV of:
  XGBoost(scores only)
  XGBoost(context only)
  XGBoost(scores + context)
against uniform ensemble. Median protein-level Spearman per fold.

Outputs:
  results/statistics/xgboost_ablation.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold
from xgboost import XGBRegressor

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job, load_core_panel  # noqa: E402

SEED = 2026
CONTEXT = ["plddt_z", "norm_pos", "log_len", "msa_shallow", "msa_deep",
           "active_site", "binding_site", "domain", "transmembrane"]
# percentile-space family centroids (same scale as Y), matching 13's experts
EXPERTS = ["U_single_seq", "U_evolution", "U_structure"]


def main():
    job = Job("31_xgboost_ablation")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    feats = pd.read_parquet(PROCESSED / "uniprot_features.parquet")
    pmap = pd.read_parquet(PROCESSED / "uniprot_position_map.parquet")
    ref = pd.read_csv(PROCESSED.parent / "raw" / "reference" / "DMS_substitutions.csv")
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    ref_map = ref.set_index("DMS_id")
    core = load_core_panel(job)

    pm = pmap[pmap["mapped"]]
    feats_t = feats.merge(pm, left_on=["UniProt_ID", "position"],
                          right_on=["UniProt_ID", "uniprot_pos"], how="inner")
    feats_t["position"] = feats_t["target_pos"].astype(int)
    d = pd.get_dummies(feats_t.set_index(["UniProt_ID", "position"])["feature"])
    fpiv = d.groupby(level=[0, 1]).max().astype(bool).reset_index()

    df = disc.merge(struct, on=["DMS_id", "mutant"], how="left")
    df = df.merge(fpiv, on=["UniProt_ID", "position"], how="left")
    # percentile-space family centroids
    norm_u = pd.read_parquet(PROCESSED / "normalized_scores.parquet",
                             columns=["DMS_id", "mutant"] +
                             [f"u_{m}" for m in core["model"]])
    df = df.merge(norm_u, on=["DMS_id", "mutant"], how="left")
    fam_u = {f: [f"u_{m}" for m in core["model"]
                 if core.loc[core.model == m, "family"].iloc[0] == f]
             for f in ["single_seq", "evolution", "structure"]}
    for f, cols in fam_u.items():
        df[f"U_{f}"] = df[cols].mean(axis=1, skipna=True)
    df["seq_len"] = df["DMS_id"].map(ref_map["seq_len"])
    df["msa_depth"] = df["DMS_id"].map(ref_map["MSA_Neff_L_category"])
    df["norm_pos"] = df["position"] / df["seq_len"]
    df["log_len"] = np.log(df["seq_len"])
    df["msa_shallow"] = df["msa_depth"].eq("Low").astype(int)
    df["msa_deep"] = df["msa_depth"].eq("High").astype(int)
    for f in ["active_site", "binding_site", "domain", "transmembrane"]:
        if f not in df:
            df[f] = False
        df[f] = df[f].fillna(False).astype(int)
    pmed = df["plddt"].median()
    df["plddt_z"] = (df["plddt"].fillna(pmed) - df["plddt"].mean()) / df["plddt"].std()
    df = df.dropna(subset=EXPERTS + ["Y_deleter", "norm_pos", "log_len"])
    job.info(f"rows: {len(df)}; proteins: {df.UniProt_ID.nunique()}")

    groups = df["UniProt_ID"].to_numpy()
    Xc = df[CONTEXT].to_numpy(dtype=np.float32)
    S = df[EXPERTS].to_numpy(dtype=np.float32)
    y = df["Y_deleter"].to_numpy(dtype=np.float32)
    uni = df["UniProt_ID"].to_numpy()

    gkf = GroupKFold(n_splits=5)
    configs = {"scores_only": S, "context_only": Xc,
               "scores_context": np.hstack([S, Xc])}
    results = []
    for name, X in configs.items():
        meds = []
        for fold, (tr, te) in enumerate(gkf.split(X, y, groups)):
            model = XGBRegressor(n_estimators=300, max_depth=5,
                                 learning_rate=0.1, subsample=0.8,
                                 colsample_bytree=0.8, n_jobs=4,
                                 random_state=SEED, eval_metric="mae")
            model.fit(X[tr], y[tr])
            p = model.predict(X[te])
            rhos = []
            for g in np.unique(uni[te]):
                m = uni[te] == g
                if m.sum() >= 20:
                    rhos.append(spearmanr(p[m], y[te][m])[0])
            meds.append(np.median(rhos))
        results.append({"config": name, "fold_medians": meds,
                        "median_rho": np.median(meds)})
        job.info(f"XGBoost {name}: median fold rho = {np.median(meds):.4f}")

    # uniform ensemble reference
    ens = S.mean(axis=1)
    rhos = []
    for g in np.unique(uni):
        m = uni == g
        if m.sum() >= 20:
            rhos.append(spearmanr(ens[m], y[m])[0])
    results.append({"config": "uniform_ensemble", "fold_medians": [],
                    "median_rho": np.median(rhos)})
    job.info(f"uniform ensemble: median rho = {np.median(rhos):.4f}")

    out = pd.DataFrame(results)
    out.to_csv(STATISTICS / "xgboost_ablation.csv", index=False)
    job.close()


if __name__ == "__main__":
    main()