"""32_identity_threshold_sensitivity.py

Review-response analysis (item 10): sensitivity of the strict non-overlap
clinical transfer to the sequence-identity threshold (30% / 50% / 70%).

For each threshold, report: number of excluded clinical proteins, retained
variants, and median per-protein uniform-ensemble AUROC.

Outputs:
  results/statistics/identity_threshold_sensitivity.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from difflib import SequenceMatcher
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import PROCESSED, STATISTICS, Job, read_reference_clinical  # noqa: E402

THRESHOLDS = [0.30, 0.50, 0.70]


def max_identity(target, dms_targets):
    tmers = {target[i:i + 8] for i in range(max(0, len(target) - 7))}
    best = 0.0
    for cand in dms_targets:
        cmers = set()
        for i in range(max(0, len(cand) - 7)):
            cmers.add(cand[i:i + 8])
        shared = len(tmers & cmers) / max(len(tmers), 1)
        if shared < 0.25:
            continue
        r = SequenceMatcher(None, target, cand, autojunk=False).ratio()
        if r > best:
            best = r
        if best >= 0.70:
            break
    return best


def main():
    job = Job("32_identity_threshold_sensitivity")
    ref = read_reference_clinical()
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    dref = pd.read_csv(PROCESSED.parent / "raw" / "reference" / "DMS_substitutions.csv")
    dms_targets = dref["target_seq"].astype(str).tolist()
    cli = pd.read_parquet(PROCESSED / "clinical_disagreement.parquet")
    score_dir = PROCESSED.parent / "raw" / "proteingym_clinical_scores"
    job.info(f"clinical proteins: {len(ref)}; DMS proteins: {len(dms_targets)}")

    # max identity per clinical protein
    ident = {}
    for row in ref.itertuples():
        ident[row.DMS_id] = max_identity(str(row.target_seq), dms_targets)
    ident_df = pd.DataFrame({"DMS_id": list(ident.keys()),
                             "max_identity": list(ident.values())})
    ident_df.to_csv(STATISTICS / "clinical_max_identity.csv", index=False)

    rows = []
    for thr in THRESHOLDS:
        excl = set(k for k, v in ident.items() if v >= thr)
        keep = ~cli["DMS_id"].isin(excl)
        job.info(f"threshold >= {thr:.0%}: excluded {len(excl)} proteins; "
                 f"retained {int(keep.sum())} variants")
        # uniform AUROC on retained proteins; clinical panel = core-model
        # intersection with released clinical scores
        from common import TABLES
        panel = pd.read_csv(TABLES / "model_panel.csv")
        core = panel[panel["panel"] == "core"]["model"].tolist()
        aucs = []
        for dms_id, g in cli[keep].groupby("DMS_id"):
            p = score_dir / f"{dms_id}.csv"
            if not p.exists():
                continue
            sc = pd.read_csv(p)
            lab = sc["DMS_bin_score"].astype(str).str.strip().str.lower()
            sc["label"] = lab.map({"pathogenic": 1, "benign": 0})
            sc = sc[sc["label"].notna()]
            yl = sc["label"].to_numpy().astype(int)
            if yl.size == 0 or np.unique(yl).size < 2:
                continue
            models = [m for m in core if m in sc.columns]
            if len(models) < 3:
                continue
            S = sc[models].to_numpy(dtype=float)
            ens = 1 - np.nanmean(S, axis=1)
            aucs.append(roc_auc_score(yl, ens))
        rows.append({"threshold": thr, "n_excluded": len(excl),
                     "n_variants": int(keep.sum()), "n_proteins": len(aucs),
                     "median_auroc": np.median(aucs)})
        job.info(f"  uniform AUROC (n={len(aucs)} proteins): {np.median(aucs):.4f}")

    out = pd.DataFrame(rows)
    out.to_csv(STATISTICS / "identity_threshold_sensitivity.csv", index=False)
    job.close()


if __name__ == "__main__":
    main()