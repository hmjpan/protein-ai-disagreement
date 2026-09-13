"""17_structure_mechanics.py

Enhancement analysis (review v2, item 7): do disagreement and family
expertise track residue-level structural mechanics?

Per-residue features from the ProteinGym AlphaFold2 PDBs:
  * secondary structure (pydssp: H / E / C)
  * CA contact density (8 A and 10 A neighbor counts) -- burial proxy
  * burial class (core / surface from contact-density tertiles, within protein)

Outcomes (all aggregated per protein, then summarized -- no pooled P values):
  * D_evo_struct and D_seq_struct by SS class and burial class
  * structure-family advantage by burial class
  * Spearman(contact density, D_evo_struct) within protein

Outputs:
  data/processed/structure_mechanics.parquet
  results/statistics/mechanics_disagreement_effects.csv
  results/statistics/mechanics_family_advantage.csv
  figures/main/Fig6A_mechanics_disagreement.png
  figures/main/Fig6B_mechanics_advantage.png
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from scipy.stats import spearmanr

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    PROCESSED, STRUCTURES_DIR, STATISTICS, FIGURES_MAIN, Job, load_core_panel,
)

SS_MAP = {"H": "helix", "E": "sheet", "-": "loop"}


def backbone_features(pdb_path: Path):
    """Return per-residue dict of (resnum, ss3, contact8, contact10)."""
    import pydssp
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", str(pdb_path))
    coords, resnums = [], []
    for chain in next(iter(s)):
        for res in chain:
            if res.id[0] != " ":
                continue
            try:
                coords.append([res["N"].get_coord(), res["CA"].get_coord(),
                               res["C"].get_coord(), res["O"].get_coord()])
                resnums.append(res.id[1])
            except KeyError:
                pass
    coords = np.array(coords)
    ss = pydssp.assign(coords) if len(coords) else np.array([])
    ca = coords[:, 1, :]
    d = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=2)
    np.fill_diagonal(d, 999.0)
    n8 = (d < 8.0).sum(axis=1)
    n10 = (d < 10.0).sum(axis=1)
    return resnums, ss, n8, n10


def main():
    job = Job("17_structure_mechanics")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    df = disc.merge(struct, on=["DMS_id", "mutant"], how="left")

    feat_rows = []
    for pdb in sorted(STRUCTURES_DIR.glob("*.pdb")):
        uni = pdb.stem
        if uni not in set(df["UniProt_ID"]):
            continue
        try:
            resnums, ss, n8, n10 = backbone_features(pdb)
        except Exception as e:
            job.info(f"{uni}: backbone failed {e}")
            continue
        ss3 = np.array([SS_MAP.get(x, "loop") for x in ss])
        for r, s3, c8, c10 in zip(resnums, ss3, n8, n10):
            feat_rows.append({"UniProt_ID": uni, "position": r,
                              "ss3": s3, "contact8": int(c8), "contact10": int(c10)})
    mech = pd.DataFrame(feat_rows)
    job.info(f"structure mechanics residues: {len(mech)} from "
             f"{mech.UniProt_ID.nunique()} proteins")
    # burial tertiles within protein
    grp = mech.groupby("UniProt_ID")["contact8"]
    mech["c8_rank"] = grp.rank(pct=True)
    mech["burial"] = np.where(mech["c8_rank"] >= 0.67, "core",
                              np.where(mech["c8_rank"] <= 0.33, "surface", "medium"))
    mech.to_parquet(PROCESSED / "structure_mechanics.parquet", index=False)

    df = df.merge(mech.drop(columns=["contact10"]), on=["UniProt_ID", "position"],
                  how="left")
    job.info(f"merged mechanics onto variants: {df['ss3'].notna().mean():.3f}")

    # ---- disagreement by SS class / burial (within-protein) ----
    eff_rows = []
    for uni, g in df.groupby("UniProt_ID"):
        g = g.dropna(subset=["ss3"])
        if len(g) < 100:
            continue
        base = {"UniProt_ID": uni}
        for col, name in [("D_evo_struct", "D_evo_struct"),
                          ("D_seq_struct", "D_seq_struct"),
                          ("D_std", "D_std")]:
            if col not in g:
                continue
            for cls in ["helix", "sheet", "loop"]:
                med = g.loc[g.ss3 == cls, col].median()
                base[f"{name}_{cls}"] = med if len(g.loc[g.ss3 == cls]) else np.nan
            for cls in ["core", "surface"]:
                med = g.loc[g.burial == cls, col].median()
                base[f"{name}_{cls}"] = med if len(g.loc[g.burial == cls]) else np.nan
            ok = g[col].notna() & g["contact8"].notna()
            if ok.sum() > 50:
                base[f"rho_{name}_contact"] = spearmanr(
                    g.loc[ok, col], g.loc[ok, "contact8"])[0]
        eff_rows.append(base)
    eff = pd.DataFrame(eff_rows)
    eff.to_csv(STATISTICS / "mechanics_disagreement_effects.csv", index=False)
    for c in ["D_evo_struct_helix", "D_evo_struct_sheet", "D_evo_struct_loop",
              "D_evo_struct_core", "D_evo_struct_surface", "rho_D_evo_struct_contact"]:
        if c in eff:
            v = eff[c].dropna()
            job.info(f"{c}: median={np.median(v):.4f} (n={len(v)} proteins)")

    # ---- family advantage by burial / SS (percentile u scale) ----
    core = load_core_panel(job)
    norm = pd.read_parquet(PROCESSED / "normalized_scores.parquet",
                           columns=["DMS_id", "mutant"] +
                           [f"u_{m}" for m in core["model"]])
    fam_u = {f: [f"u_{m}" for m in core["model"]
                 if core.loc[core.model == m, "family"].iloc[0] == f]
             for f in ["single_seq", "evolution", "structure"]}
    for f, cols in fam_u.items():
        df[f"U_{f}"] = norm[cols].mean(axis=1, skipna=True)
    fam_labels = {"U_single_seq": "seq", "U_evolution": "evo", "U_structure": "struct"}
    df["error"] = (df[list(fam_labels)].mean(axis=1) - df["Y_deleter"]).abs()
    adv_rows = []
    for fam, lab in fam_labels.items():
        other = [c for c in fam_labels if c != fam]
        E_fam = (df[fam] - df["Y_deleter"]).abs()
        E_oth = df[other].sub(df["Y_deleter"], axis=0).abs().mean(axis=1)
        adv = E_oth - E_fam
        for ctx, col in [("core", "burial"), ("surface", "burial"),
                         ("helix", "ss3"), ("sheet", "ss3"), ("loop", "ss3")]:
            sel = df[df[col].eq(ctx)]
            if len(sel) < 500:
                continue
            meds = pd.DataFrame({"a": adv[sel.index],
                                 "u": sel["UniProt_ID"]}).groupby("u")["a"].median()
            adv_rows.append({"family": lab, "context": ctx, "context_type": col,
                             "n": len(sel), "median_adv": meds.median(),
                             "frac_pos": (meds > 0).mean()})
    adv_df = pd.DataFrame(adv_rows)
    adv_df.to_csv(STATISTICS / "mechanics_family_advantage.csv", index=False)
    for _, r in adv_df.iterrows():
        job.info(f"  {r.family} @ {r.context}: adv={r.median_adv:+.4f} "
                 f"(frac_pos {r.frac_pos:.2f}, n={r.n})")

    # ---- figures ----
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, c in enumerate(["D_evo_struct_core", "D_evo_struct_surface",
                           "D_evo_struct_helix", "D_evo_struct_sheet",
                           "D_evo_struct_loop"]):
        v = eff[c].dropna()
        ax.boxplot([v], positions=[i], widths=0.6)
    ax.set_xticks(range(5))
    ax.set_xticklabels(["core", "surface", "helix", "sheet", "loop"])
    ax.set_ylabel("median evolution-vs-structure disagreement (protein-level)")
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig6A_mechanics_disagreement.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    labs = sorted(adv_df["family"].unique())
    ctxs = ["core", "surface", "helix", "sheet", "loop"]
    width = 0.25
    for i, lab in enumerate(labs):
        sub = adv_df[adv_df.family == lab]
        vals = [sub.loc[sub.context == c, "median_adv"].iloc[0] for c in ctxs]
        ax.bar(np.arange(len(ctxs)) + i * width, vals, width, label=lab)
    ax.set_xticks(np.arange(len(ctxs)) + width)
    ax.set_xticklabels(ctxs)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("median family advantage")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES_MAIN / "Fig6B_mechanics_advantage.png", dpi=300)
    plt.close(fig)
    job.info("figures written: Fig6A, Fig6B")
    job.close()


if __name__ == "__main__":
    main()