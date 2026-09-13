"""42_homology_standard.py

Reviewer-response: replace the 8-mer heuristic similarity with STANDARD
protein homology search for the strict non-overlap clinical set.

Tool: HMMER3 `phmmer` (Pyhmmer bindings) -- query = each clinical benchmark
target protein (HMM built per query), database = the DMS assay target
sequences. E-value cutoff 1e-3.

Per clinical protein we record:
  * best-hit target, bit score, E-value
  * pairwise identity over aligned residues (from the alignment strings)
  * query coverage = aligned span / full query length

Exclusion rule (standard): identity >= T AND query coverage >= 70%,
T in {30, 50, 70}%.
Also reports per-protein uniform-ensemble AUROC (identical computation to
32_identity_threshold_sensitivity) and the median AUROC of retained
proteins at each threshold, plus overlap with the 35 proteins excluded by
the old heuristic.

Outputs:
  results/statistics/homology_standard_per_protein.csv
  results/statistics/homology_standard_sensitivity.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.hmmer import phmmer

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (PROCESSED, RAW, STATISTICS, TABLES, Job,
                    read_reference_clinical)  # noqa: E402

QCOV_MIN = 70.0
ZS = ["TranceptEVE_L", "GEMME", "EVE", "ESM1b", "PoET"]
THRESHOLDS = [30, 50, 70]


def _dec(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def main():
    job = Job("42_homology_standard")
    ref = read_reference_clinical()
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    dref = pd.read_csv(RAW / "reference" / "DMS_substitutions.csv")
    dref["DMS_id"] = dref["DMS_id"].astype(str)
    targets = dref[["DMS_id", "target_seq"]].drop_duplicates().astype(str)
    job.info(f"queries: {len(ref)} clinical; db: {len(targets)} DMS targets")

    abc = Alphabet.amino()
    db = [TextSequence(name=t, sequence=s).digitize(abc)
          for t, s in zip(targets["DMS_id"], targets["target_seq"])]
    queries = [TextSequence(name=i, sequence=str(s)).digitize(abc)
               for i, s in zip(ref["DMS_id"], ref["target_seq"])]

    rows = []
    for k, th in enumerate(phmmer(queries, db, cpus=8, E=1e-3,
                                  domE=1e-3), 1):
        q = _dec(th.query.name)
        qlen = int(th.query.L)
        best = None
        for hit in th.reported:
            for dom in hit.domains:
                if dom.i_evalue > 1e-3:
                    continue
                al = dom.alignment

                def _s(x):
                    return x.decode() if isinstance(x, bytes) else str(x)
                hs, ts = _s(al.hmm_sequence), _s(al.target_sequence)
                cols = sum(1 for a, b in zip(hs, ts)
                           if a not in "-." and b not in "-.")
                if not cols:
                    continue
                matches = sum(1 for a, b in zip(hs, ts)
                              if a not in "-." and b not in "-."
                              and a.upper() == b.upper())
                ident = 100.0 * matches / cols
                span = al.hmm_to - al.hmm_from + 1
                cov = 100.0 * min(1.0, span / max(qlen, 1))
                key = (ident * cov, dom.score)
                if best is None or key > best[0]:
                    best = (key, _s(hit.name), ident, cov,
                            float(dom.i_evalue), float(dom.score))
        if best is None:
            rows.append({"DMS_id": q, "best_target": "", "max_pident": 0.0,
                         "qcov_pct": 0.0, "evalue": np.nan, "bitscore": 0.0})
        else:
            rows.append({"DMS_id": q, "best_target": best[1],
                         "max_pident": best[2], "qcov_pct": best[3],
                         "evalue": best[4], "bitscore": best[5]})
        if k % 500 == 0:
            job.info(f"queries done: {k}")
    per = pd.DataFrame(rows)
    per["DMS_id"] = per["DMS_id"].astype(str)
    per.to_csv(STATISTICS / "homology_standard_per_protein.csv", index=False)
    n_hits = int((per["max_pident"] > 0).sum())
    job.info(f"queries with any E<=1e-3 hit: {n_hits}")

    old = pd.read_csv(STATISTICS / "clinical_max_identity.csv")
    old35 = set(old.loc[old["max_identity"] >= 0.70, "DMS_id"].astype(str))

    # uniform AUROC per protein (identical to 32)
    cli_dir = PROCESSED.parent / "raw" / "proteingym_clinical_scores"
    panel = pd.read_csv(TABLES / "model_panel.csv")
    core = panel[panel["panel"] == "core"]["model"].tolist()
    auroc = {}
    for dms_id in per["DMS_id"]:
        p = cli_dir / f"{dms_id}.csv"
        if not p.exists():
            continue
        sc = pd.read_csv(p)
        lab = sc["DMS_bin_score"].astype(str).str.strip().str.lower()
        sc["label"] = lab.map({"pathogenic": 1, "benign": 0})
        sc = sc[sc["label"].notna()]
        yl = sc["label"].to_numpy().astype(int)
        if yl.size < 10 or np.unique(yl).size < 2:
            continue
        models = [m for m in ZS if m in sc.columns]
        if len(models) < 3:
            continue
        M = sc[models].to_numpy(dtype=float)
        for j, m in enumerate(models):
            if m == "PoET":
                M[:, j] = -M[:, j]  # official clinical metadata dir=-1
        R = np.apply_along_axis(rankdata, 0, M)
        ens = 1 - np.nanmean(R, axis=1)
        auroc[dms_id] = roc_auc_score(yl, ens)
    per["uniform_auroc"] = per["DMS_id"].map(auroc)
    per.to_csv(STATISTICS / "homology_standard_per_protein.csv", index=False)

    rows2 = []
    for t in THRESHOLDS:
        excl = set(per.loc[(per["max_pident"] >= t)
                           & (per["qcov_pct"] >= QCOV_MIN), "DMS_id"].astype(str))
        kept = per[~per["DMS_id"].isin(excl)]["uniform_auroc"].dropna()
        rows2.append({"threshold_pident": t, "qcov_min_pct": QCOV_MIN,
                      "n_excluded": len(excl),
                      "overlap_with_old35": len(excl & old35),
                      "n_proteins_with_auroc": len(kept),
                      "median_uniform_auroc": float(np.median(kept))})
        job.info(f"pident>={t} & qcov>={QCOV_MIN:.0f}: excluded {len(excl)} "
                 f"(overlap old35 {len(excl & old35)}); "
                 f"median AUROC {np.median(kept):.4f}")
    pd.DataFrame(rows2).to_csv(STATISTICS / "homology_standard_sensitivity.csv",
                               index=False)
    job.close()


if __name__ == "__main__":
    main()
