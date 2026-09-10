# -*- coding: utf-8 -*-
"""Build supplementary.docx: one page per supplementary figure (S1-S11) and
real Word tables for Tables S1-S8, from the pipeline result CSVs.

Captions for the tables are parsed live from manuscript/submission/
manuscript.md so there is a single source of truth.
"""
import csv
import io
import re
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
STAT = ROOT / "results" / "statistics"
OUT = ROOT / "manuscript" / "submission"

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt


def fmt(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return str(v)
    if abs(x - round(x)) < 1e-12 and abs(x) < 1e6:
        return f"{int(round(x)):,}" if abs(x) >= 1000 else str(int(round(x)))
    return f"{x:.3f}"


def rows_of(fname, cols=None, keep=None):
    with open(STAT / fname, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if keep:
        rows = [r for r in rows if keep(r)]
    if cols:
        rows = [{c: r[c] for c in cols} for r in rows]
    return rows


FIG_S = {
    1: "Figure S1. Assay-level distribution of the Spearman correlation between "
       "disagreement and ensemble prediction error (median 0.032; Section 2.6).",
    2: "Figure S2. Mean ensemble prediction error by disagreement decile "
       "(monotonic: top decile 0.214 vs bottom decile 0.201).",
    3: "Figure S3. Disagreement versus distance to the nearest UniProt domain "
       "boundary (no boundary effect).",
    4: "Figure S4. BioGate 5-fold UniProt-grouped cross-validation: protein-level "
       "Spearman distributions against all baselines.",
    5: "Figure S5. BioGate and XGBoost ablations (scores-only, context-only, "
       "scores+context; Table S5).",
    6: "Figure S6. Leave-one-model-out sensitivity of the disagreement-error "
       "association (range 0.027-0.039).",
    7: "Figure S7. Percentile (u) versus inverse-normal (z) disagreement scale: "
       "median assay-level rho 0.145 vs 0.032 at fixed error definition.",
    8: "Figure S8. Disagreement-error association by assay selection type "
       "(stability 0.076, binding 0.057; activity 0.003, organismal fitness -0.004).",
    9: "Figure S9. Residue mechanics (secondary structure, burial, contact "
       "density) versus family advantage or disagreement: null result.",
    10: "Figure S10. ESM-IF1 gain from experimental structures over AlphaFold2 "
        "by assay mean pLDDT (64 assays; low-confidence assays show MORE negative "
        "gain: -0.081 vs -0.023, difference -0.058, 95% CI -0.099 to -0.013).",
    11: "Figure S11. SSEmb distinct-architecture replication: residue-level "
        "Spearman between pLDDT and evolution-vs-SSEmb disagreement (median -0.067; "
        "65.6% of 151 proteins negative).",
}

md = (OUT / "manuscript.md").read_text(encoding="utf-8")
seg = md[md.find("**Table S1."):]
seg = seg[:seg.find("## Figure captions") if "## Figure captions" in seg[:200] else len(seg)]
TBL_CAP = {}
parts = re.split(r"(?=\*\*Table S\d)", seg)
for pt in parts:
    m = re.match(r"\*\*Table (S\d)\.(.*?)\n\n", pt, re.S)
    if not m:
        m = re.match(r"\*\*Table (S\d)\.(.*)$", pt, re.S)
    if m:
        TBL_CAP[m.group(1)] = "Table " + m.group(1) + "." + m.group(2)
assert len(TBL_CAP) == 8, sorted(TBL_CAP)

doc = Document()
st = doc.styles["Normal"]
st.font.name = "Times New Roman"
st.font.size = Pt(10.5)


def para(text, bold=False, size=None, center=False):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.bold = bold
    if size:
        r.font.size = Pt(size)
    if center:
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    return p


def add_table(headers, data, widths=None):
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = "Table Grid"
    for c, h in zip(t.rows[0].cells, headers):
        run = c.paragraphs[0].add_run(h)
        run.bold = True
        run.font.size = Pt(9)
    for row in data:
        cells = t.add_row().cells
        for c, v in zip(cells, row):
            r = c.paragraphs[0].add_run(str(v))
            r.font.size = Pt(9)
    return t


# ---- cover page ----
para("Supplementary Information", bold=True, size=16, center=True)
para("Protein AI Model Disagreement Tracks AlphaFold Structural Confidence",
     bold=True, size=13, center=True)
para("")
para("Contents", bold=True)
para("Supplementary Figures S1-S11 (one figure per page).", size=10.5)
para("Supplementary Tables S1-S8.", size=10.5)
para("")
para("All underlying machine-readable tables and pipeline code (scripts 01-37) "
     "are archived with the manuscript at "
     "https://github.com/hmjpan/protein-ai-disagreement (results/ and figures/ "
     "directories). Statistical aggregation is protein-level with bootstrap "
     "inference throughout; the analysis plan and GO/NO-GO thresholds were fixed "
     "in config.yaml before the final analyses.", size=10.5)

# ---- supplementary figures ----
doc.add_page_break()
para("Supplementary Figures", bold=True, size=14)
figdir = ROOT / "figures" / "supplementary"
for n in range(1, 12):
    files = sorted(f for f in figdir.iterdir()
                   if re.match(rf"FigS{n}[_\.]", f.name))
    doc.add_page_break()
    para(FIG_S[n], bold=True, size=10.5)
    for f in files:
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.add_run().add_picture(str(f), width=Inches(6.2))

# ---- supplementary tables ----
doc.add_page_break()
para("Supplementary Tables", bold=True, size=14)


def cap(key):
    doc.add_page_break()
    para(TBL_CAP[key].replace("\n", " "), bold=True, size=10.5)


# S1: disagreement definitions
cap("S1")
add_table(["Definition", "n proteins", "median rho", "95% CI",
           "frac proteins negative"],
          [[r["definition"], r["n_proteins"], fmt(r["median_rho"]),
            f"{fmt(r['ci_low'])} to {fmt(r['ci_high'])}", fmt(r["frac_negative"])]
           for r in rows_of("disagreement_definition_plddt_summary.csv")])

# S2: confound controls
cap("S2")
asym = rows_of("confound_control_asymmetry.csv")
mr_e = sorted(float(r["rho_evo_struct"]) for r in asym)
mr_s = sorted(float(r["rho_seq_struct"]) for r in asym)
med = lambda v: v[len(v) // 2]
para("(i) Asymmetry of the residue-level association:", bold=True, size=10)
add_table(["Pair", "median rho vs pLDDT", "n proteins"],
          [["evolution vs structure", f"{med(mr_e):.3f}", len(asym)],
           ["sequence vs structure", f"{med(mr_s):.3f}", len(asym)],
           ["% proteins with more negative evo-vs-struct rho",
            f"{sum(1 for r in asym if float(r['rho_evo_struct']) < float(r['rho_seq_struct'])) / len(asym) * 100:.0f}%", "-"]])
para("")
para("(ii) Structure-model DMS correlation by pLDDT bin:", bold=True, size=10)
add_table(["pLDDT bin", "model", "n variants", "n assays", "median rho DMS"],
          [[r["plddt_bin"], r["model"], fmt(r["n_variants"]), r["n_assays"],
            fmt(r["median_rho_dms"])]
           for r in rows_of("confound_control_accuracy.csv")])
para("")
para("(iii) Structure-model score coverage:", bold=True, size=10)
add_table(["pLDDT bin", "model", "fraction scored"],
          [[r["plddt_bin"], r["model"], fmt(r["frac_scored"])]
           for r in rows_of("confound_control_coverage.csv")])

# S3: regime robustness across K
cap("S3")
ari = {(r["K_other"], r["K_ref"]): float(r["ARI"])
       for r in rows_of("regime_k_ari.csv")}
data3 = []
for r in rows_of("regime_k_robustness.csv"):
    k = r["K"]
    a = "1.000" if k == "6" else fmt(ari.get((k, "6")))
    data3.append([k, r["has_structure_dissenting"], fmt(r["consensus_tolerant_Y"]),
                  fmt(r["consensus_damaging_Y"]), fmt(r["mean_crossassay_agreement"]),
                  r["n_crossassay_pairs"], a])
add_table(["K", "structure-dissenting regime present", "Y consensus-tolerant",
           "Y consensus-damaging", "cross-assay agreement", "cross-assay pairs",
           "ARI vs K=6"], data3)

# S4: identity thresholds
cap("S4")
add_table(["Identity threshold", "proteins excluded", "variants", "proteins",
           "median uniform AUROC"],
          [[f"{int(float(r['threshold']) * 100)}%", r["n_excluded"],
            fmt(r["n_variants"]), fmt(r["n_proteins"]), fmt(r["median_auroc"])]
           for r in rows_of("identity_threshold_sensitivity.csv")])

# S5: XGBoost ablations
cap("S5")
add_table(["Configuration", "median protein-level Spearman"],
          [[r["config"], fmt(r["median_rho"])]
           for r in rows_of("xgboost_ablation.csv")])

# S6: model panel composition
cap("S6")
add_table(["Model", "Official family", "Analysis family", "Uses MSA",
           "Uses structure", "Structure input AlphaFold", "In mechanistic panel"],
          [[r["model"], r["official_family"], r["family_in_analysis"],
            r["uses_MSA"], r["uses_structure"], r["structure_input_is_AlphaFold"],
            r["in_mechanistic_panel"]]
           for r in rows_of("model_panel_table.csv")])

# S7: clinical paired deltas
cap("S7")
add_table(["Set", "Comparison", "n proteins", "median delta AUROC",
           "95% CI", "frac positive"],
          [[r["set"], r["comparison"], fmt(r["n_proteins"]), fmt(r["median_delta"]),
            f"[{fmt(r['ci_low'])}, {fmt(r['ci_high'])}]", fmt(r["frac_positive"])]
           for r in rows_of("clinical_paired_deltas.csv")])

# S8: held-out regime generalization
cap("S8")
h = rows_of("heldout_regime_summary.csv")[0]
add_table(["Statistic", "Value"],
          [["Proteins with held-out support", h["n_proteins"]],
           ["Median Spearman, held-out regime profile vs training profile",
            f"{fmt(h['median_rho_train_profile'])} (IQR {fmt(h['q25_rho'])}-{fmt(h['q75_rho'])})"],
           ["Proteins with positive correlation", fmt(h["frac_rho_positive"])],
           ["Dissenting regime held-out Y (median)", fmt(h["dissenting_Y_median_heldout"])],
           ["Proteins where dissenting regime damaging (Y > 0.5)",
            fmt(h["frac_dissenting_Y_gt05"])],
           ["Held-out cross-assay replication, median rho",
            f"{fmt(h['heldout_crossassay_rho_median'])} ({h['heldout_crossassay_n_pairs']} pairs)"],
           ["Cross-assay pairs above chance", fmt(h["frac_pairs_above_chance"])]])
para("Per-protein values: results/statistics/heldout_regime_per_protein.csv.",
     size=9)

doc.save(str(OUT / "supplementary.docx"))
print("supplementary.docx:", (OUT / "supplementary.docx").stat().st_size,
      "bytes | tables:", len(TBL_CAP), "| figs:", len(FIG_S))
