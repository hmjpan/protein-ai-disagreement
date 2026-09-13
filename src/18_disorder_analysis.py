"""18_disorder_analysis.py

Enhancement analysis (review v2, item 6): does evolution-vs-structure
disagreement specifically accumulate in intrinsically disordered regions?

Two independent disorder annotations:
  * UniProt curated "Disordered" regions (features of type Region with
    disordered description), transferred to DMS coordinates via the
    alignment map produced in 10.
  * AlphaFold pLDDT < 50 as a disorder-correlated proxy (kept under the
    'low structural confidence' terminology; Piovesan et al. 2022 show
    pLDDT predicts experimental disorder -- cited in Discussion, not used
    to relabel our metric).

Outcomes (within-protein, then summarized):
  * median D_evo_struct / D_seq_struct / D_std in IDR vs ordered residues
  * family advantage in IDR vs ordered
  * overlap of UniProt IDR with pLDDT<50 (annotation agreement)

Outputs:
  data/processed/disorder_features.parquet
  results/statistics/disorder_disagreement_effects.csv
  results/statistics/disorder_family_advantage.csv
  figures/main/Fig6C_disorder_disagreement.png
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
    PROCESSED, STATISTICS, FIGURES_MAIN, Job,
)


def main():
    job = Job("18_disorder_analysis")
    disc = pd.read_parquet(PROCESSED / "disagreement_scores.parquet")
    struct = pd.read_parquet(PROCESSED / "structure_features.parquet",
                             columns=["DMS_id", "mutant", "plddt"])
    feats = pd.read_parquet(PROCESSED / "uniprot_features.parquet")
    pmap = pd.read_parquet(PROCESSED / "uniprot_position_map.parquet")
    df = disc.merge(struct, on=["DMS_id", "mutant"], how="left")

    # UniProt disordered regions (Region feature with disordered description)
    dis_rows = []
    for data in _iter_uniprot_json(job):
        pass  # placeholder replaced below
    # (feature extraction done from the stored parquet of features plus
    # descriptions: re-fetch descriptions from cached JSONs)
    from common import RAW
    import json
    dis_pos = {}
    for p in sorted((RAW / "uniprot").glob("*.json")):
        uni = p.stem
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        for feat in data.get("features", []):
            desc = (feat.get("description") or "").lower()
            if "disord" in desc or "flexib" in desc:
                loc = feat.get("location", {})
                s = loc.get("start", {}).get("value")
                e = loc.get("end", {}).get("value")
                if s:
                    dis_pos.setdefault(uni, []).append((int(s), int(e)))
    job.info(f"proteins with UniProt disordered regions: {len(dis_pos)}")

    # map UniProt positions -> DMS positions
    pm = pmap[pmap["mapped"]]
    dis_t = []
    for uni, ranges in dis_pos.items():
        sub = pm[pm["UniProt_ID"] == uni]
        upos = set(sub["uniprot_pos"].astype(int))
        for s, e in ranges:
            for pos in range(s, e + 1):
                t = sub.loc[sub["uniprot_pos"] == pos, "target_pos"]
                if len(t):
                    dis_t.append({"UniProt_ID": uni, "position": int(t.iloc[0])})
    dis_df = pd.DataFrame(dis_t)
    job.info(f"disordered positions mapped to DMS coordinates: {len(dis_df)}")

    df = df.merge(dis_df.assign(disorder=True), on=["UniProt_ID", "position"],
                  how="left")
    df["disorder"] = df["disorder"].fillna(False).astype(bool)
    df["lowconf"] = df["plddt"] < 50
    df.to_parquet(PROCESSED / "disorder_features.parquet", index=False)

    # annotation agreement
    ag = df[["UniProt_ID", "disorder", "lowconf", "plddt"]].dropna(subset=["plddt"])
    job.info(f"UniProt disorder rate: {ag.disorder.mean():.3f}; "
             f"pLDDT<50 rate: {ag.lowconf.mean():.3f}")
    both = ag[ag.disorder | ag.lowconf]
    job.info(f"overlap: P(pLDDT<50 | UniProt disorder)="
             f"{both.loc[both.disorder, 'lowconf'].mean():.3f}; "
             f"P(disorder | pLDDT<50)={both.loc[both.lowconf, 'disorder'].mean():.3f}")

    # within-protein effects
    eff_rows = []
    for uni, g in df.groupby("UniProt_ID"):
        if g["disorder"].sum() < 20:
            continue
        base = {"UniProt_ID": uni, "n_idr": int(g["disorder"].sum())}
        for col in ["D_evo_struct", "D_seq_struct", "D_std"]:
            if col not in g:
                continue
            base[f"{col}_idr"] = g.loc[g.disorder, col].median()
            base[f"{col}_ord"] = g.loc[~g.disorder, col].median()
        eff_rows.append(base)
    eff = pd.DataFrame(eff_rows)
    eff.to_csv(STATISTICS / "disorder_disagreement_effects.csv", index=False)
    if len(eff):
        for col in ["D_evo_struct", "D_seq_struct", "D_std"]:
            d = eff[f"{col}_idr"] - eff[f"{col}_ord"]
            job.info(f"{col}: median IDR-ordered diff = {np.median(d):+.4f} "
                     f"({len(eff)} proteins; frac positive {(d > 0).mean():.3f})")

    # family advantage in IDR vs ordered
    fam_labels = {"S_single_seq": "seq", "S_evolution": "evo", "S_structure": "struct"}
    df["error"] = (df[["S_single_seq", "S_evolution", "S_structure"]]
                   .mean(axis=1) - df["Y_deleter"]).abs()
    adv_rows = []
    for fam, lab in fam_labels.items():
        other = [c for c in fam_labels if c != fam]
        E_fam = (df[fam] - df["Y_deleter"]).abs()
        E_oth = df[other].sub(df["Y_deleter"], axis=0).abs().mean(axis=1)
        adv = E_oth - E_fam
        for ctx, col in [("IDR", "disorder"), ("ordered", "disorder"),
                         ("lowconf", "lowconf"), ("highconf", "lowconf")]:
            sel = df[df[col].eq(ctx in ("IDR", "lowconf"))]
            if len(sel) < 500:
                continue
            meds = pd.DataFrame({"a": adv[sel.index],
                                 "u": sel["UniProt_ID"]}).groupby("u")["a"].median()
            adv_rows.append({"family": lab, "context": ctx, "n": len(sel),
                             "median_adv": meds.median(), "frac_pos": (meds > 0).mean()})
    adv_df = pd.DataFrame(adv_rows)
    adv_df.to_csv(STATISTICS / "disorder_family_advantage.csv", index=False)
    for _, r in adv_df.iterrows():
        job.info(f"  {r.family} @ {r.context}: adv={r.median_adv:+.4f} "
                 f"(frac_pos {r.frac_pos:.2f}, n={r.n})")

    # figure
    if len(eff):
        fig, ax = plt.subplots(figsize=(6, 5))
        for i, col in enumerate(["D_evo_struct", "D_seq_struct", "D_std"]):
            ax.boxplot([eff[f"{col}_idr"] - eff[f"{col}_ord"]], positions=[i], widths=0.5)
        ax.axhline(0, color="red", ls="--", lw=1)
        ax.set_xticks(range(3))
        ax.set_xticklabels(["D_evo_struct", "D_seq_struct", "D_std"])
        ax.set_ylabel("median disagreement (IDR - ordered), per protein")
        fig.tight_layout()
        fig.savefig(FIGURES_MAIN / "Fig6C_disorder_disagreement.png", dpi=300)
        plt.close(fig)
    job.info("figures written: Fig6C")
    job.close()


def _iter_uniprot_json(job):
    return []


if __name__ == "__main__":
    main()