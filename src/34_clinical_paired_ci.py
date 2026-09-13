"""34_clinical_paired_ci.py

Reviewer-response (items 6/8): paired protein-level bootstrap CI for
clinical BioGate vs uniform vs best-single differences.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import STATISTICS, Job  # noqa: E402

SEED = 2026
N = 10000


def paired(d, a, b):
    rng = np.random.default_rng(SEED)
    v = d[a].to_numpy() - d[b].to_numpy()
    boots = [np.median(rng.choice(v, size=len(v), replace=True)) for _ in range(N)]
    return v, float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def main():
    job = Job("34_clinical_paired_ci")
    rows = []
    for label, f in [("full", "clinical_transfer_performance.csv"),
                     ("strict", "clinical_transfer_nonoverlap.csv")]:
        d = pd.read_csv(STATISTICS / f)
        for a, b in [("auroc_biogate", "auroc_uniform"),
                     ("auroc_best", "auroc_uniform"),
                     ("auroc_best", "auroc_biogate")]:
            v, lo, hi = paired(d, a, b)
            rows.append({"set": label, "comparison": f"{a} - {b}",
                         "n_proteins": len(d), "median_delta": float(np.median(v)),
                         "ci_low": lo, "ci_high": hi, "frac_positive": float((v > 0).mean())})
            job.info(f"[{label}] {a} - {b}: median {np.median(v):+.4f} "
                     f"[{lo:+.4f}, {hi:+.4f}], {((v>0).mean()):.2f} positive, n={len(d)}")
    pd.DataFrame(rows).to_csv(STATISTICS / "clinical_paired_deltas.csv", index=False)
    job.close()


if __name__ == "__main__":
    main()