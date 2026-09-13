# Model Disagreement as a Scientific Observable: Structure and Boundaries of Protein AI Divergence

Analysis code, results tables and figures for the manuscript
*"Model Disagreement as a Scientific Observable: Structure and Boundaries of Protein AI Divergence"*.

Every figure in the manuscript has a machine-readable source table under
`results/`, and all analyses run from public data with fixed settings
recorded in `config.yaml`.

## Data dependencies (all public)

1. **ProteinGym v1.3** -- https://doi.org/10.5281/zenodo.15293562
   DMS substitution assays, precomputed zero-shot model scores, the clinical
   benchmark, and AlphaFold2 structure files. Downloaded by
   `src/01_download_data.py`.
2. **UniProt REST API** -- features and curated disordered regions
   (`src/09_structure_annotations.py`, `src/36_aligner_mapping.py`).
3. **Gitter-lab experimental-structure benchmark** --
   https://doi.org/10.5281/zenodo.13819824 and
   https://github.com/gitter-lab/benchmarking-structure-based-models
   (ESM-IF1 scores on experimental structures; SSEmb scores) used by
   `src/35_expstruct_control.py`.

## Layout

```
src/            numbered pipeline scripts (01-44); run in numeric order
                41 matched two-family DMS/clinical control
                42 standard-homology exclusion via HMMER3 phmmer (pyhmmer)
                43 structure-subset sensitivity; 44 held-out consensus removal
config.yaml     fixed analysis settings (seed, panels, thresholds, QC rules)
results/        every table/CSV/parquet behind each figure and statistic
figures/        main panels, composed Figure 1-6 files (figures/composed/)
                and supplementary S1-S11
manuscript/     submission manuscript, Supplementary Information, cover letter
logs/           run logs of the pipeline
requirements.txt
```

## Reproducing the analysis

```bash
pip install -r requirements.txt
python src/01_download_data.py      # downloads public data (large)
for s in src/0[2-9]_*.py src/1?_*.py src/2?_*.py src/3?_*.py src/4?_*.py; do python "$s"; done
```

`src/42_homology_standard.py` requires `pyhmmer` (HMMER3). The confirmatory
analysis specification and GO/NO-GO thresholds are recorded in `config.yaml`,
frozen before the reported runs (git history timestamps each version).

## License

MIT (see `LICENSE`). Third-party data remain subject to their original terms
(ProteinGym Zenodo record; UniProt and ClinVar terms of use).
