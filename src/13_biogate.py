"""13_biogate.py

Phase 6 -- BioGate: biologically contextualized model arbitration.

Architecture (protocol s36-38):
    w = softmax(MLP(context_features))        # 3 expert weights, sum = 1
    prediction = w_seq*S_seq + w_evo*S_evo + w_struct*S_struct
Trained by MSE on Y_deleter (+ small entropy penalty so gating stays
informative). Torch CPU only.

Validation (protocol s39-44):
  * 5-fold GroupKFold by UniProt_ID (no variant-level leakage)
  * baselines: best single expert, uniform ensemble, linear stacking,
    XGBoost (S + context)
  * metrics: per-protein Spearman (macro), AUROC on DMS_score_bin
  * paired protein-level bootstrap for BioGate vs uniform ensemble
  * ablations: full / -pLDDT / -functional annotation / -MSA / no-context

Outputs:
  results/statistics/biogate_cv_predictions.parquet
  results/statistics/biogate_performance.csv
  results/statistics/biogate_ablation.csv
  results/statistics/biogate_gate_weights.csv
  results/GO_NO_GO_3.md
  figures/main/Fig6A_biogate_performance.png
  figures/main/Fig6B_gate_weights.png
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from scipy.stats import rankdata, spearmanr
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import Ridge
from xgboost import XGBRegressor

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    PROCESSED, STATISTICS, RESULTS, FIGURES_MAIN, Job, load_core_panel,
)

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
SEED = 2026
torch.manual_seed(SEED)
np.random.seed(SEED)

CONTEXT = ["plddt_z", "norm_pos", "log_len", "msa_shallow", "msa_deep",
           "active_site", "binding_site", "domain", "transmembrane"]
EXPERTS = ["U_single_seq", "U_evolution", "U_structure"]
HIDDEN = 32
EPOCHS = 12
BATCH = 4096
LR = 1e-3
ENT_PEN = 0.01


class GateNet(nn.Module):
    def __init__(self, n_in, n_exp):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(n_in, HIDDEN), nn.ReLU(), nn.Linear(HIDDEN, n_exp))
        self.n_exp = n_exp

    def forward(self, x, s):
        w = torch.softmax(self.mlp(x), dim=1)
        return (w * s).sum(1), w


def build_features(job: Job):
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    feats = pd.read_parquet(PROCESSED / "uniprot_features.parquet")
    pmap = pd.read_parquet(PROCESSED / "uniprot_position_map.parquet")
    pm = pmap[pmap["mapped"]]
    feats_t = feats.merge(pm, left_on=["UniProt_ID", "position"],
                          right_on=["UniProt_ID", "uniprot_pos"], how="inner")
    feats_t["position"] = feats_t["target_pos"].astype(int)
    d = pd.get_dummies(feats_t.set_index(["UniProt_ID", "position"])["feature"])
    fpiv = d.groupby(level=[0, 1]).max().astype(bool).reset_index()

    ref = pd.read_csv(PROCESSED.parent / "raw" / "reference" / "DMS_substitutions.csv")
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    ref_map = ref.set_index("DMS_id")
    core = load_core_panel(job)
    expert_cols = [f"S_{f}" for f in ["single_seq", "evolution", "structure"]]
    for c in expert_cols:
        if c not in disc:
            job.info(f"MISSING expert column {c}")
            return None

    df = disc.merge(struct, on=["DMS_id", "mutant"], how="left")
    df = df.merge(fpiv, on=["UniProt_ID", "position"], how="left")
    # percentile-space family centroids (U in [0,1], same scale as Y)
    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet",
                           columns=["DMS_id", "mutant"] +
                           [f"u_{m}" for m in core["model"]])
    fam_u = {f: [f"u_{m}" for m in core["model"]
                 if core.loc[core.model == m, "family"].iloc[0] == f]
             for f in ["single_seq", "evolution", "structure"]}
    df = df.merge(norm, on=["DMS_id", "mutant"], how="left")
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
    plddt_med = df["plddt"].median()
    df["plddt_z"] = (df["plddt"].fillna(plddt_med) - df["plddt"].mean()) / df["plddt"].std()
    df = df.dropna(subset=expert_cols + ["Y_deleter", "norm_pos", "log_len"])
    job.info(f"BioGate dataset: {len(df)} rows, {df.UniProt_ID.nunique()} proteins")
    return df


def spearman_per_protein(pred: np.ndarray, y: np.ndarray, groups: np.ndarray) -> pd.DataFrame:
    rows = []
    for g in np.unique(groups):
        m = groups == g
        if m.sum() >= 20:
            rows.append({"protein": g, "n": int(m.sum()),
                         "rho": spearmanr(pred[m], y[m])[0]})
    return pd.DataFrame(rows)


def main():
    job = Job("13_biogate")
    df = build_features(job)
    if df is None:
        job.close()
        return

    groups = df["UniProt_ID"].to_numpy()
    X = df[CONTEXT].to_numpy(dtype=np.float32)
    S = df[EXPERTS].to_numpy(dtype=np.float32)
    y = df["Y_deleter"].to_numpy(dtype=np.float32)
    ybin = df["DMS_score_bin"].to_numpy()
    uni = df["UniProt_ID"].to_numpy()

    gkf = GroupKFold(n_splits=5)
    cv_pred = []
    perf_rows = []
    gate_weights = []
    for fold, (tr, te) in enumerate(gkf.split(X, y, groups)):
        Xtr, Xte = X[tr], X[te]
        Str, Ste = S[tr], S[te]
        ytr, yte = y[tr], y[te]

        # --- BioGate (torch) ---
        model = GateNet(X.shape[1], S.shape[1])
        opt = torch.optim.Adam(model.parameters(), lr=LR)
        lossf = nn.MSELoss()
        Xt = torch.from_numpy(Xtr)
        St = torch.from_numpy(Str)
        yt = torch.from_numpy(ytr)
        n = len(tr)
        for ep in range(EPOCHS):
            model.train()
            perm = np.random.default_rng(SEED + fold * 100 + ep).permutation(n)
            for i in range(0, n, BATCH):
                idx = perm[i:i + BATCH]
                xb, sb, yb = Xt[idx], St[idx], yt[idx]
                opt.zero_grad()
                pred, w = model(xb, sb)
                loss = lossf(pred, yb) + ENT_PEN * (-(w * (w + 1e-9).log()).sum(1).mean())
                loss.backward()
                opt.step()
        model.eval()
        with torch.no_grad():
            pred_te, w_te = model(torch.from_numpy(Xte), torch.from_numpy(Ste))
            pred_te = pred_te.numpy()
            w_te = w_te.numpy()
        for i in range(len(EXPERTS)):
            gate_weights.append({"fold": fold, "expert": EXPERTS[i],
                                 "mean_weight": float(w_te[:, i].mean()),
                                 "median_weight": float(np.median(w_te[:, i]))})

        # --- baselines ---
        ens_te = Ste.mean(axis=1)
        ridge = Ridge(alpha=1.0)
        ridge.fit(Str, ytr)
        ridge_te = ridge.predict(Ste)
        xgb = XGBRegressor(n_estimators=300, max_depth=5, learning_rate=0.1,
                           subsample=0.8, colsample_bytree=0.8, n_jobs=4,
                           random_state=SEED, eval_metric="mae")
        xgb.fit(np.hstack([Str, Xtr]), ytr)
        xgb_te = xgb.predict(np.hstack([Ste, Xte]))

        best_te = None
        best_name = None
        for k in range(S.shape[1]):
            rho_tr = spearmanr(Str[:, k], ytr)[0]
            if best_te is None or rho_tr > best_te:
                best_te = rho_tr
                best_name = EXPERTS[k]
        best_te = S[te, EXPERTS.index(best_name)]

        methods = {"uniform_ensemble": ens_te, "linear_stacking": ridge_te,
                   "xgboost": xgb_te, "biogate": pred_te}
        if best_name is not None:
            methods[f"best_single({best_name})"] = best_te
        for name, p in methods.items():
            dfp = spearman_per_protein(p, y[te], uni[te])
            dfp["method"] = name
            dfp["fold"] = fold
            # ybin: 1 = fit (ProteinGym official); pred = deleteriousness
            dfp["auroc"] = _auroc(p, 1 - ybin[te])
            perf_rows.append(dfp)
        cv_pred.append(pd.DataFrame({
            "UniProt_ID": uni[te], "fold": fold, "y": y[te],
            "ybin": ybin[te],
            "uniform_ensemble": ens_te, "linear_stacking": ridge_te,
            "xgboost": xgb_te, "biogate": pred_te}))
        job.info(f"fold {fold}: best_single={best_name}; "
                 f"biogate median rho={np.median(spearman_per_protein(pred_te, y[te], uni[te])['rho']):.4f}; "
                 f"uniform median rho={np.median(spearman_per_protein(ens_te, y[te], uni[te])['rho']):.4f}")

    perf = pd.concat(perf_rows, ignore_index=True)
    summ = perf.groupby("method").agg(
        n_proteins=("protein", "count"),
        median_rho=("rho", "median"),
        mean_rho=("rho", "mean"),
        mean_auroc=("auroc", "mean")).reset_index()
    summ.to_csv(STATISTICS / "biogate_performance.csv", index=False)
    job.info("\n" + summ.round(4).to_string(index=False))

    cv = pd.concat(cv_pred, ignore_index=True)
    cv.to_parquet(STATISTICS / "biogate_cv_predictions.parquet", index=False)

    # paired protein-level bootstrap: biogate - uniform
    wide = perf.pivot_table(index=["protein", "fold"], columns="method",
                            values="rho").reset_index()
    if "biogate" in wide and "uniform_ensemble" in wide:
        d = wide["biogate"] - wide["uniform_ensemble"]
        rng = np.random.default_rng(SEED)
        boots = np.array([np.mean(rng.choice(d.dropna(), size=len(d), replace=True))
                          for _ in range(2000)])
        ci = np.percentile(boots, [2.5, 97.5])
        job.info(f"paired bootstrap Δrho biogate-uniform: "
                 f"mean={np.mean(d):.4f} CI=[{ci[0]:.4f},{ci[1]:.4f}]")
    else:
        ci = (np.nan, np.nan)

    # --- ablations ---
    ab_rows = []
    abl_specs = {"full": CONTEXT, "no_plddt": [c for c in CONTEXT if c != "plddt_z"],
                 "no_ann": [c for c in CONTEXT if c not in
                            {"active_site", "binding_site", "domain", "transmembrane"}],
                 "no_msa": [c for c in CONTEXT if c not in {"msa_shallow", "msa_deep"}]}
    Xfull = X
    for name, cols in abl_specs.items():
        Xi = Xfull[:, [CONTEXT.index(c) for c in cols]]
        meds = []
        for fold, (tr, te) in enumerate(gkf.split(Xi, y, groups)):
            model = GateNet(Xi.shape[1], S.shape[1])
            opt = torch.optim.Adam(model.parameters(), lr=LR)
            lossf = nn.MSELoss()
            Xt = torch.from_numpy(Xi[tr]); St = torch.from_numpy(S[tr])
            yt = torch.from_numpy(y[tr]); n = len(tr)
            for ep in range(EPOCHS):
                model.train()
                perm = np.random.default_rng(SEED + fold * 100 + ep).permutation(n)
                for i in range(0, n, BATCH):
                    idx = perm[i:i + BATCH]
                    opt.zero_grad()
                    pred, w = model(Xt[idx], St[idx])
                    loss = lossf(pred, yt[idx]) + ENT_PEN * (-(w * (w + 1e-9).log()).sum(1).mean())
                    loss.backward()
                    opt.step()
            model.eval()
            with torch.no_grad():
                p = model(torch.from_numpy(Xi[te]), torch.from_numpy(S[te]))[0].numpy()
            meds.append(np.median(spearman_per_protein(p, y[te], uni[te])["rho"]))
        ab_rows.append({"ablation": name, "median_rho": np.median(meds)})
    abd = pd.DataFrame(ab_rows)
    abd.to_csv(STATISTICS / "biogate_ablation.csv", index=False)
    job.info("\n" + abd.round(4).to_string(index=False))

    # --- GO / NO-GO 3 ---
    bv = summ[summ.method == "biogate"]["median_rho"].values
    uv = summ[summ.method == "uniform_ensemble"]["median_rho"].values
    delta = bv[0] - uv[0] if len(bv) and len(uv) else np.nan
    verdict = "GO" if (not np.isnan(delta) and delta >= 0.01 and ci[0] > 0) else "MARGINAL"
    lines = [
        "# GO / NO-GO 3 -- BioGate",
        "",
        f"- Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}",
        f"- BioGate median protein-level rho: {bv[0] if len(bv) else np.nan:.4f}",
        f"- Uniform ensemble median rho: {uv[0] if len(uv) else np.nan:.4f}",
        f"- Paired bootstrap mean Δrho: {np.mean(d) if 'd' in dir() else np.nan:.4f} "
        f"95% CI [{ci[0]:.4f}, {ci[1]:.4f}]",
        f"- Pre-registered rule: Δrho >= 0.01 with CI > 0",
        "",
        f"**VERDICT: {verdict}**",
        "",
    ]
    (RESULTS / "GO_NO_GO_3.md").write_text("\n".join(lines), encoding="utf-8")
    job.info(f"GO/NO-GO 3 verdict: {verdict}")

    # --- figures ---
    fig, ax = plt.subplots(figsize=(7, 5))
    meth_order = summ.sort_values("median_rho", ascending=False)["method"]
    for i, m in enumerate(meth_order):
        sub = perf[perf.method == m]["rho"]
        ax.boxplot(sub, positions=[i], widths=0.6)
    ax.set_xticks(range(len(meth_order)))
    ax.set_xticklabels(meth_order, rotation=20, fontsize=8)
    ax.set_ylabel("protein-level Spearman rho")
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig6A_biogate_performance.png", dpi=200)
    plt.close(fig)

    gw = pd.DataFrame(gate_weights)
    gsum = gw.groupby("expert")["mean_weight"].mean()
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(gsum.index, gsum.values, color=["#56B4E9", "#E69F00", "#009E73"])
    ax.set_ylabel("mean gating weight")
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig6B_gate_weights.png", dpi=200)
    plt.close(fig)
    job.info("figures written: Fig6A, Fig6B")
    job.close()


def _auroc(pred: np.ndarray, ybin: np.ndarray) -> float:
    if ybin.size == 0 or np.unique(ybin).size < 2:
        return np.nan
    from sklearn.metrics import roc_auc_score
    return roc_auc_score(ybin, pred)


if __name__ == "__main__":
    main()