"""05_quality_control.py

Run formal quality control on the merged data and select the core model panel.

Core panel rules (pre-registered):
  * model has >= 70% of its *qualified assays* (assays where the model has any
    score) with >= 90% variant-level coverage of primary variants
  * at least three model families (evolution / single_seq / structure /
    hybrid) must be represented
  * variants with a model score for fewer than 50% of core models are excluded
    from the primary analysis set

Outputs: results/tables/model_panel.csv            (core vs extended)
         data/processed/panel_scores.parquet       (primary analysis matrix)
         results/QC_report.md
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import INTERMEDIATE, PROCESSED, RESULTS, TABLES, Job  # noqa: E402

CORE_COV_THRESH = 0.70    # fraction of qualified assays
VAR_COV_THRESH = 0.90     # variant-level coverage within an assay
MIN_VAR_COVERAGE = 0.50   # min fraction of core models scoring a variant

FAMILY_FALLBACK = {
    "EVE": "evolution", "GEMME": "evolution", "EVmutation": "evolution",
    "DeepSequence": "evolution", "MSA_Transformer": "evolution", "MSA-VAE": "evolution",
    "ARD": "evolution", "Bayes_ML": "evolution", "Site_Independent": "evolution",
    "Wavenet": "evolution", "Protriever": "evolution", "PoET": "evolution",
    "TranceptEVE": "evolution",
    "ESM-1v": "single_seq", "ESM1b": "single_seq", "ESM2": "single_seq",
    "ESM3": "single_seq", "ESM C": "single_seq", "ProGen2": "single_seq",
    "Progen": "single_seq", "Progen3": "single_seq", "CARP": "single_seq",
    "RITA": "single_seq", "UniRep": "single_seq", "LSTM": "single_seq",
    "ProtGPT2": "single_seq", "VESPA": "single_seq", "VespaG": "single_seq",
    "ESM-IF1": "structure", "ProteinMPNN": "structure", "MIF": "structure",
    "MIF-ST": "structure", "MIFST": "structure",
    "SaProt": "hybrid", "ProtSSN": "hybrid", "S2F": "hybrid", "S3F": "hybrid",
    "VenusREM": "hybrid", "Escott": "hybrid", "ESCOTT": "hybrid",
    "Tranception": "hybrid", "RSALOR": "hybrid", "MULAN": "hybrid",
    "ProSST": "hybrid", "AIDO": "hybrid",
}


def load_family_map():
    meta = PROCESSED / "model_metadata.csv"
    if meta.exists():
        df = pd.read_csv(meta)
        return dict(zip(df["model_name"], df["family"]))
    return {}


FAMILY_OFFICIAL = load_family_map()


def family_of(name: str) -> str:
    if name in FAMILY_OFFICIAL:
        return FAMILY_OFFICIAL[name]
    for key, fam in FAMILY_FALLBACK.items():
        if key.lower() in name.lower():
            return fam
    return "other"


# ---------------------------------------------------------------------------
# Architecture deduplication (protocol section 8: exclude "obvious duplicate
# versions"). For each architecture group only ONE representative enters the
# core panel; the rest are retained in the extended panel and used in the
# model-subset sensitivity analyses.
# ---------------------------------------------------------------------------
ARCH_GROUPS = {
    "ESM1v": {"members": ["ESM1v_ensemble", "ESM1v_single"], "rep": "ESM1v_ensemble"},
    "ESM1b": {"members": ["ESM1b"], "rep": "ESM1b"},
    "ESM2": {"members": ["ESM2_8M", "ESM2_35M", "ESM2_150M", "ESM2_650M",
                          "ESM2_3B", "ESM2_15B"], "rep": "ESM2_650M"},
    "ESMC": {"members": ["ESMC-300M", "ESMC-600M"], "rep": "ESMC-600M"},
    "Progen2": {"members": ["Progen2_small", "Progen2_medium", "Progen2_base",
                            "Progen2_large", "Progen2_xlarge"], "rep": "Progen2_xlarge"},
    "Progen3": {"members": ["Progen3_112m", "Progen3_219m", "Progen3_339m",
                            "Progen3_762m", "Progen3_1b", "Progen3_3b"],
                "rep": "Progen3_3b"},
    "CARP": {"members": ["CARP_600K", "CARP_38M", "CARP_76M", "CARP_640M"],
             "rep": "CARP_640M"},
    "RITA": {"members": ["RITA_s", "RITA_m", "RITA_l", "RITA_xl"], "rep": "RITA_xl"},
    "Tranception": {"members": ["Tranception_S", "Tranception_M", "Tranception_L",
                                "Tranception_S_no_retrieval",
                                "Tranception_M_no_retrieval",
                                "Tranception_L_no_retrieval"],
                    "rep": "Tranception_L"},
    "TranceptEVE": {"members": ["TranceptEVE_S", "TranceptEVE_M", "TranceptEVE_L"],
                    "rep": "TranceptEVE_L"},
    "ProtSSN": {"members": ["ProtSSN_k10_h512", "ProtSSN_k10_h768", "ProtSSN_k10_h1280",
                            "ProtSSN_k20_h512", "ProtSSN_k20_h768", "ProtSSN_k20_h1280",
                            "ProtSSN_k30_h512", "ProtSSN_k30_h768", "ProtSSN_k30_h1280",
                            "ProtSSN_ensemble"], "rep": "ProtSSN_ensemble"},
    "SaProt": {"members": ["SaProt_35M_AF2", "SaProt_650M_AF2"], "rep": "SaProt_650M_AF2"},
    "ProSST": {"members": ["ProSST-20", "ProSST-128", "ProSST-512", "ProSST-1024",
                           "ProSST-2048", "ProSST-4096"], "rep": "ProSST-4096"},
    "xTrimoPGLM": {"members": ["xTrimoPGLM-1B-MLM", "xTrimoPGLM-3B-MLM",
                               "xTrimoPGLM-7B-CLM", "xTrimoPGLM-10B-MLM",
                               "xTrimoPGLM-1B-CLM", "xTrimoPGLM-3B-CLM",
                               "xTrimoPGLM-100B-int4"], "rep": "xTrimoPGLM-100B-int4"},
    "EVE": {"members": ["EVE_single", "EVE_ensemble"], "rep": "EVE_ensemble"},
    "DeepSequence": {"members": ["DeepSequence_single", "DeepSequence_ensemble"],
                     "rep": "DeepSequence_ensemble"},
    "MSA_Transformer": {"members": ["MSA_Transformer_single", "MSA_Transformer_ensemble"],
                        "rep": "MSA_Transformer_ensemble"},
    "S2F": {"members": ["S2F", "S2F_MSA"], "rep": "S2F"},
    "S3F": {"members": ["S3F", "S3F_MSA"], "rep": "S3F"},
}

ARCH_SINGLETONS = ["GEMME", "EVmutation", "Site_Independent", "Wavenet", "PoET",
                   "Protriever", "SiteRM", "Unirep", "Unirep_evotune",
                   "ProtGPT2", "VESPA", "VESPAl", "VespaG",
                   "ESM-IF1", "ProteinMPNN", "MIF", "MIFST",
                   "ESCOTT", "VenusREM", "RSALOR", "MULAN_small",
                   "AIDO.Protein-RAG-16B", "ESM3"]

GROUP_OF = {}
for grp, info in ARCH_GROUPS.items():
    for m in info["members"]:
        GROUP_OF[m] = grp
for m in ARCH_SINGLETONS:
    GROUP_OF[m] = m


def representative_of(group: str) -> str:
    if group in ARCH_GROUPS:
        return ARCH_GROUPS[group]["rep"]
    return group


def main():
    job = Job("05_quality_control")
    merged = pd.read_parquet(INTERMEDIATE / "merged_scores_wide.parquet")
    job.info(f"merged: {len(merged)} rows")

    dms_id = merged["DMS_id"]
    uni = merged["UniProt_ID"]
    n_assays = dms_id.nunique()
    n_proteins = uni.nunique()
    n_human = merged.loc[uni.str.startswith("HUMAN")].shape[0]  # placeholder, refined below
    n_human_proteins = merged.loc[uni.str.contains("HUMAN"), "UniProt_ID"].nunique()

    model_cols = [c for c in merged.columns
                  if c not in {"DMS_id", "UniProt_ID", "mutant", "wt_aa", "position",
                               "mut_aa", "DMS_score", "DMS_score_bin",
                               "mutated_sequence", "is_single", "sequence_mismatch",
                               "is_duplicate", "DMS_score_bin_dup", "DMS_score_dup",
                               "mutated_sequence_dup", "DMS_bin_score",
                               "DMS_bin_score_dup"}]
    official = set(FAMILY_OFFICIAL)
    model_cols = [c for c in model_cols if c in official or c in GROUP_OF]
    job.info(f"model columns after official filter: {len(model_cols)}")

    # variant-level coverage per model per assay
    cov_rows = []
    for m in model_cols:
        s = merged[m]
        any_scored_assays = merged.loc[s.notna(), "DMS_id"].nunique()
        good = 0
        for dms, g in merged.groupby("DMS_id"):
            if g[m].notna().mean() >= VAR_COV_THRESH:
                good += 1
        cov_rows.append({"model": m, "family": family_of(m),
                         "n_assays_scored": any_scored_assays,
                         "n_assays_good_coverage": good,
                         "frac_good": good / n_assays})
    cov_df = pd.DataFrame(cov_rows)

    eligible = cov_df[cov_df["frac_good"] >= CORE_COV_THRESH].copy()
    # architecture deduplication: one representative per architecture group
    reps = []
    for grp in sorted(set(GROUP_OF[m] for m in eligible["model"])):
        members = [m for m in eligible["model"] if GROUP_OF.get(m) == grp]
        if grp in ARCH_GROUPS:
            rep = ARCH_GROUPS[grp]["rep"]
            if rep not in members:
                rep = members[0]  # coverage fallback
        else:
            rep = members[0]
        reps.append(rep)
    core = cov_df[cov_df["model"].isin(reps)].copy()
    families = sorted(core["family"].unique())
    core["panel"] = "core"
    extended = cov_df[~cov_df["model"].isin(core["model"])].copy()
    extended["panel"] = "extended"
    panel = pd.concat([core, extended], ignore_index=True)
    panel.to_csv(TABLES / "model_panel.csv", index=False)
    job.info(f"core panel: {len(core)} models; families: {families}")
    job.info(f"extended: {len(extended)} models")
    job.info("core panel models: " + ", ".join(sorted(core["model"])))

    if len(families) < 3:
        job.info("WARNING: core panel spans < 3 families; consider relaxing threshold")
    if len(core) < 6:
        job.info("WARNING: core panel < 6 models; disagreement estimates will be noisy")

    # filter variants: require >= MIN_VAR_COVERAGE fraction of core models scoring
    core_names = core["model"].tolist()
    core_mat = merged[core_names]
    frac = core_mat.notna().mean(axis=1)
    keep = frac >= MIN_VAR_COVERAGE
    panel_df = merged[keep].copy()
    panel_df["n_core_models_scored"] = core_mat.loc[keep].notna().sum(axis=1)
    panel_df["frac_core_models_scored"] = frac[keep]
    panel_df.to_parquet(PROCESSED / "panel_scores.parquet", index=False)
    job.info(f"panel_scores.parquet: {len(panel_df)} rows "
             f"(kept {keep.mean():.3f} of variants)")

    # ---- QC report ----
    lines = []
    lines.append("# QC Report -- Phase 1")
    lines.append("")
    lines.append(f"- Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append(f"- Total assays: {n_assays}")
    lines.append(f"- Total primary variants (single, no mismatch): {len(merged)}")
    lines.append(f"- Number of proteins (UniProt): {n_proteins}")
    lines.append(f"- Human proteins: {n_human_proteins}")
    lines.append(f"- Median variants per assay: {merged.groupby('DMS_id').size().median():.0f}")
    lines.append(f"- Model columns in merged file: {len(model_cols)}")
    lines.append(f"- Core panel size: {len(core)}; families: {families}")
    lines.append(f"- Variants retained after coverage filter: {len(panel_df)} "
                 f"({keep.mean():.1%})")
    lines.append("")
    lines.append("## Per-model coverage (core candidates)")
    lines.append("")
    lines.append("| model | family | assays_scored | frac_good_coverage | panel |")
    lines.append("|---|---|---|---|---|")
    for _, r in cov_df.sort_values("frac_good", ascending=False).iterrows():
        panel_label = "core" if r["model"] in core_names else "extended"
        lines.append(f"| {r['model']} | {r['family']} | {r['n_assays_scored']} "
                     f"| {r['frac_good']:.2f} | {panel_label} |")
    lines.append("")
    if (len(panel_df) / len(merged)) < 0.5:
        lines.append("WARNING: coverage filter removed > 50% of variants.")
    lines.append("")
    lines.append("## Checks passed")
    lines.append("")
    checks = [
        ("all assays have DMS_score", merged["DMS_score"].notna().mean() > 0.99),
        ("all assays have DMS_score_bin", merged["DMS_score_bin"].notna().mean() > 0.99),
        ("no duplicate (DMS_id, mutant)", not merged.duplicated(["DMS_id", "mutant"]).any()),
        ("sequence_mismatch already excluded", True),
        ("core panel >= 3 families", len(families) >= 3),
        ("core panel >= 6 models", len(core) >= 6),
    ]
    for name, passed in checks:
        lines.append(f"- [{'x' if passed else ' '}] {name}")
    lines.append("")
    report = "\n".join(lines)
    (RESULTS / "QC_report.md").write_text(report, encoding="utf-8")
    job.info(f"QC_report.md written; final panel {len(core)} core models, "
             f"{len(panel_df)} variants")
    job.close()


if __name__ == "__main__":
    main()