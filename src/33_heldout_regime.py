"""33_heldout_regime.py

Reviewer-response (item 5): TRUE held-out protein validation of regimes.

Protocol:
  * split the 186 proteins into 5 folds (GroupKFold by UniProt)
  * on each training fold fit GMM (K=6) on z-space family centroids; freeze
  * assign held-out protein variants to frozen regimes (argmax responsibility)
  * validation metrics per held-out protein:
      - regime size distribution (min regime support, fraction assigned)
      - regime-level median experimental Y (held-out data)
      - phenotype replication: Spearman(train regime-Y profile vs held-out
        regime-Y profile over the 6 regimes)
      - structure-dissenting regime damaging fraction on held-out proteins
      - cross-assay replication on held-out proteins only

Outputs: results/statistics/heldout_regime_{per_protein,summary}.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job  # noqa: E402

SEED = 2026
K = 6
S_COLS = ["S_single_seq", "S_evolution", "S_structure"]


def main():
    job = Job("33_heldout_regime")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    keep = ~disc[S_COLS].isna().any(axis=1)
    d = disc[keep].reset_index(drop=True)
    X = d[S_COLS].to_numpy(dtype=float)
    y = d["Y_deleter"].to_numpy()
    prot = d["UniProt_ID"].to_numpy()
    job.info(f"variants: {len(d)}; proteins: {len(np.unique(prot))}")

    rows = []
    gkf = GroupKFold(n_splits=5)
    for fold, (tr, te) in enumerate(gkf.split(X, y, prot)):
        gmm = GaussianMixture(n_components=K, random_state=SEED + fold,
                              covariance_type="full", max_iter=300, n_init=3)
        gmm.fit(X[tr])
        # train regime-Y profile (global)
        lab_tr = gmm.predict(X[tr])
        prof_tr = pd.Series(y[tr]).groupby(lab_tr).median()
        # held-out assignment + profile
        lab_te = gmm.predict(X[te])
        dte = pd.DataFrame({"prot": prot[te], "reg": lab_te, "y": y[te],
                            "dms": d["DMS_id"].to_numpy()[te]})
        # structure-dissenting regime identified on TRAIN
        cent = pd.DataFrame(gmm.means_, columns=S_COLS)
        diss = cent[(cent["S_single_seq"] > 0.3) & (cent["S_evolution"] > 0.3)
                    & (cent["S_structure"] < 0.3)].index.tolist()
        for uni, g in dte.groupby("prot"):
            if len(g) < 100:
                continue
            prof_te = g.groupby("reg")["y"].median()
            common = prof_tr.index.intersection(prof_te.index)
            rho = spearmanr(prof_tr[common], prof_te[common])[0] if len(common) >= 4 else np.nan
            dis_rows = g[g.reg.isin(diss)]
            rows.append({"UniProt_ID": uni, "fold": fold, "n": len(g),
                         "n_regimes": g.reg.nunique(),
                         "rho_train_profile": rho,
                         "frac_assigned_dissenting": len(dis_rows) / len(g),
                         "dissenting_Y": dis_rows["y"].median() if len(dis_rows) > 20 else np.nan})
    rep = pd.DataFrame(rows)
    rep.to_csv(STATISTICS / "heldout_regime_per_protein.csv", index=False)

    # cross-assay replication restricted to HELD-OUT proteins
    lab_all = []
    for fold, (tr, te) in enumerate(gkf.split(X, y, prot)):
        gmm = GaussianMixture(n_components=K, random_state=SEED + fold,
                              covariance_type="full", max_iter=300, n_init=3)
        gmm.fit(X[tr])
        lab = gmm.predict(X[te])
        t = d.iloc[te].copy()
        t["reg"] = lab
        lab_all.append(t)
    allo = pd.concat(lab_all)
    rep2_rows = []
    for uni, g in allo.groupby("UniProt_ID"):
        assays = g["DMS_id"].unique()
        if len(assays) < 2:
            continue
        for a in range(len(assays) - 1):
            ga = g[g.DMS_id == assays[a]]; gb = g[g.DMS_id == assays[a + 1]]
            ya = ga.groupby("reg")["Y_deleter"].median()
            yb = gb.groupby("reg")["Y_deleter"].median()
            na = ga.groupby("reg").size(); nb = gb.groupby("reg").size()
            common = [r for r in ya.index if r in yb.index and na[r] >= 20 and nb[r] >= 20]
            if len(common) >= 3:
                rep2_rows.append(spearmanr(ya[common], yb[common])[0])
    summary = {
        "n_proteins": len(rep),
        "median_rho_train_profile": float(np.nanmedian(rep.rho_train_profile)),
        "q25_rho": float(np.nanpercentile(rep.rho_train_profile, 25)),
        "q75_rho": float(np.nanpercentile(rep.rho_train_profile, 75)),
        "frac_rho_positive": float((rep.rho_train_profile > 0).mean()),
        "dissenting_Y_median_heldout": float(np.nanmedian(rep.dissenting_Y)),
        "frac_dissenting_Y_gt05": float((rep.dissenting_Y > 0.5).mean()),
        "heldout_crossassay_rho_median": float(np.median(rep2_rows)) if rep2_rows else None,
        "heldout_crossassay_n_pairs": len(rep2_rows),
        "frac_pairs_above_chance": None,
    }
    # chance for cross-assay rho on shared regimes
    rng = np.random.default_rng(SEED)
    chance = []
    for _ in range(5000):
        vals = rng.normal(size=len(common))
        chance.append(np.abs(spearmanr(vals, rng.normal(size=len(common)))[0]))
    summary["frac_pairs_above_chance"] = float(np.mean(
        np.array(rep2_rows) > np.median(chance))) if rep2_rows else None
    pd.DataFrame([summary]).to_csv(STATISTICS / "heldout_regime_summary.csv", index=False)
    for k, v in summary.items():
        job.info(f"{k}: {v}")
    job.close()


if __name__ == "__main__":
    main()