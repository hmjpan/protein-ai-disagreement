"""15_case_studies.py

Phase 7b -- biomedical case studies.

Candidates chosen from the data AFTER the main analyses (protocol s50; never
pre-selected around a narrative): TP53, BRCA1, PTEN -- each has dense DMS
coverage, high-quality structures, and clear functional biology.

Visualisation: AlphaFold CA-trace colored by residue-level median model
disagreement; alongside per-family median prediction scores. Outputs are
machine-readable tables + PNG renders (matplotlib 3D -> 2D projection).

Outputs:
  results/tables/case_study_<GENE>_residue.csv
  figures/main/Fig7_case_<GENE>_disagreement.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    PROCESSED, STRUCTURES_DIR, STATISTICS, TABLES, FIGURES_MAIN, Job,
    load_core_panel, read_reference_dms,
)

CASE_GENES = ["P53_HUMAN", "BRCA1_HUMAN", "PTEN_HUMAN"]
CMAP = plt.cm.viridis


def ca_trace(pdb_path: Path):
    xs, ys, zs, resnums = [], [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
                resnums.append(int(line[22:26]))
    return np.array(xs), np.array(ys), np.array(zs), np.array(resnums)


def main():
    job = Job("15_case_studies")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    df = disc.merge(struct, on=["DMS_id", "mutant"], how="left")
    core = load_core_panel(job)
    ref = read_reference_dms()
    ref["DMS_id"] = ref["DMS_id"].astype(str)
    uni_of = dict(zip(ref["DMS_id"], ref["UniProt_ID"]))

    CASE_LETTER = {"P53_HUMAN": "B", "BRCA1_HUMAN": "C", "PTEN_HUMAN": "D"}
    case_letter = CASE_LETTER
    for gene in CASE_GENES:
        sel = df[df["UniProt_ID"].eq(gene)]
        if sel.empty:
            job.info(f"{gene}: no DMS data")
            continue
        uni_id = gene
        pdb = STRUCTURES_DIR / f"{uni_id}.pdb"
        if not pdb.exists():
            job.info(f"{gene}: no AF2 structure file")
            continue
        xs, ys, zs, resnums = ca_trace(pdb)
        # residue-level medians
        rl = sel.groupby("position").agg(
            D=("D_std", "median"),
            Sseq=("S_single_seq", "median"),
            Sevo=("S_evolution", "median"),
            Sstruct=("S_structure", "median"),
            plddt=("plddt", "median"),
            n=("D_std", "size")).reset_index()
        rl.to_csv(TABLES / f"case_study_{gene}_residue.csv", index=False)
        job.info(f"{gene}: {len(rl)} residues with DMS; "
                 f"n_variants={len(sel)}; "
                 f"top-D residue: {rl.loc[rl.D.idxmax(), 'position']} "
                 f"(D={rl.D.max():.3f})")
        # map DMS positions to structure residue numbers via pLDDT map
        pos_to_d = dict(zip(rl["position"], rl["D"]))
        plddt_map = dict(zip(rl["position"], rl["plddt"]))
        D_on_struct = np.array([pos_to_d.get(r, np.nan) for r in resnums])
        p_on_struct = np.array([plddt_map.get(r, np.nan) for r in resnums])
        ok = ~np.isnan(D_on_struct)
        if ok.sum() < 10:
            job.info(f"{gene}: too few mapped residues for render")
            continue

        # 2D projection for rendering (PCA of CA coordinates)
        coords = np.vstack([xs, ys, zs]).T
        c = coords - coords.mean(axis=0)
        _, _, Vt = np.linalg.svd(c)
        proj = c @ Vt[:2].T

        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
        for ax, (col, title) in zip(axes, [
                ("D", "AI disagreement (median)"),
                ("Sstruct", "structure-family score"),
                ("plddt", "AlphaFold pLDDT")]):
            vals = np.array([(dict(zip(rl.position, rl[col])).get(r, np.nan))
                             for r in resnums])
            vmin, vmax = np.nanpercentile(vals[ok], [2, 98])
            seg = np.full(len(proj) - 1, np.nan)
            for i in range(len(proj) - 1):
                v = np.nanmean([vals[i], vals[i + 1]])
                if not np.isnan(v):
                    seg[i] = v
            sc = ax.scatter(proj[:, 0], proj[:, 1], c=vals, cmap=CMAP,
                            s=8, vmin=vmin, vmax=vmax)
            for i in range(len(proj) - 1):
                if not np.isnan(seg[i]):
                    ax.plot(proj[i:i + 2, 0], proj[i:i + 2, 1],
                            color=CMAP((seg[i] - vmin) / (vmax - vmin)),
                            lw=1.2, alpha=0.8)
            ax.set_title(f"{gene}: {title}")
            ax.axis("off")
            fig.colorbar(sc, ax=ax, fraction=0.04)
        fig.tight_layout()
        fig.savefig(FIGURES_MAIN / f"Fig4{case_letter[gene]}_case_{gene}.png", dpi=300)
        plt.close(fig)
        job.info(f"{gene}: figure written")
    job.close()


if __name__ == "__main__":
    main()