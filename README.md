# Protein AI Disagreement: Revealing Hidden Biological Regimes of Missense Variation

Scientific question: **Why do different protein AI models disagree on the same
amino-acid substitution, and does this disagreement carry biological
structure?**

This is a *discovery-driven* (not benchmark) study. The pipeline is designed so
that hypotheses are tested, not retrofitted. See `config.yaml` for all
pre-registered thresholds.

## Pipeline (Phase 1)

| Script | Purpose | Key outputs |
|---|---|---|
| `01_download_data.py` | verify data integrity, fetch official model score directions | `MANIFEST.md`, `model_direction.csv` |
| `02_inventory_data.py` | assay / model / coverage inventories | `TableS1/S2/S3` |
| `03_parse_dms.py` | parse assays, flag single substitutions & sequence mismatches | `processed_single_mutations.parquet` |
| `04_merge_model_scores.py` | merge zero-shot model scores | `merged_scores_wide.parquet` |
| `05_quality_control.py` | QC + core model panel selection | `model_panel.csv`, `panel_scores.parquet`, `QC_report.md` |
| `06_normalize_scores.py` | per-assay rank normalization, direction alignment | `normalized_scores.parquet` |
| `07_compute_disagreement.py` | D_std / D_MAD / family disagreement | `disagreement_scores.parquet`, `TableS4` |

## Hard rules (non-negotiable)

1. No fabricated numbers. Every figure must have a machine-readable source table.
2. No direction inference from DMS labels (official configs only).
3. Protein-level (UniProt) splits only; no random mutation splits.
4. AlphaFold low-pLDDT != intrinsically disordered region (wording rule).
5. Association != causation (wording rule).
6. Negative results are reported, never deleted.
7. GO/NO-GO checkpoints are evaluated against `config.yaml` thresholds.

## Data provenance

- ProteinGym v1.3 (OATML-Markslab), DOI: 10.5281/zenodo.15293562
- DMS substitutions: 217 assays / ~2.7M variants
- Clinical substitutions: 2,525 proteins / ~63K variants
- All downloads from `https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/`

## Reproduction

```bash
python src/01_download_data.py
python src/02_inventory_data.py
python src/03_parse_dms.py
python src/04_merge_model_scores.py
python src/05_quality_control.py
python src/06_normalize_scores.py
python src/07_compute_disagreement.py
```

All scripts log inputs/outputs/counts/runtime to `logs/`.
## Current manuscript status
- Analysis scripts: src/01-src/27 (01-07 Phase 1; 08-16 main analyses; 17-26 review-response enhancements; 27 panel residue-level verification)
- Submission package (v1.2): manuscript/submission/ (manuscript.md, manuscript.docx, cover_letter.md)
- Figures: 20 main + 9 supplementary, all regenerable via src/16_make_figures.py
