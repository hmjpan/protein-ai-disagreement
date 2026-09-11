# Cover letter -- AI for Science (AI4S)

Dear Editors,

We are pleased to submit our manuscript entitled:

**"Protein AI Model Disagreement Tracks AlphaFold Structural Confidence"**

for consideration for publication in *AI for Science*.

**Why this fits AI for Science.** This is an AI-driven biological discovery
study, not a benchmark. We use 40 architecturally diverse protein AI models
as a scientific instrument and show that their *disagreement* -- usually
treated as noise -- is a structured observable that tracks the biological
organization of the mutational landscape: evolution-vs-structure disagreement accumulates where AlphaFold
structural confidence is low, total disagreement is elevated in intrinsically disordered
regions, resolves into reproducible biological regimes, and reveals that
relative model-family expertise is context-specific only for curated
disorder. The framing -- "disagreement as a biological probe" -- is, to our
knowledge, new, and it falls squarely within the journal's mission of
transformative AI applications in the life sciences.

**Key results (all with protein-level statistics; analysis plan
fixed in the versioned project configuration before analysis; all
accuracy comparisons in percentile space).**
(1) Residue-level Spearman between structural confidence and
evolution-vs-structure disagreement is -0.04 (95% CI -0.06 to -0.03,
p < 0.001), robust to multivariable adjustment and to four complementary
disagreement definitions, with residues in the lowest pLDDT decile
~1.5-1.8x more likely to fall in the top disagreement decile (pooled OR
1.49, 95% CI 1.35-1.65); (2) regime-level experimental phenotypes replicate
across independent DMS assays (median Spearman 0.94), with consensus
regimes aligned with experiment in 100% of pairs; (3) relative family
expertise is context-specific only for curated disorder (single-sequence
models lose relative advantage in IDR, -0.043 vs +0.038), not along the
confidence axis; (4) a DMS-trained biologically gated ensemble transfers to
clinical variants without clinical-label retraining (median protein-level
AUROC 0.899 vs 0.867 uniform on a strictly non-overlapping 788-protein
set), while disagreement-as-uncertainty fails in clinical annotation -- a
boundary condition we report explicitly. All negative results are
reported, and the shared-input confound between structural confidence and
structure-conditioned model scores is analyzed with three partial controls
and stated as the primary limitation.

**Rigour.** Every figure has a machine-readable source table, statistical
aggregation is protein-level with bootstrap inference, and all numbers are
reproducible from public data and code (ProteinGym v1.3). A separate
Supplementary Information file contains Figures S1-S11 (one per page) and
Tables S1-S9; main-text Tables 1-4 are typeset in the manuscript and main
Figures 1-6 appear on their own pages at the end of the manuscript file.

The manuscript is original, has not been published or submitted elsewhere,
and all authors approve its submission. There are no competing interests.

**Suggested reviewers** (none has authored a model used in our panel; no
co-authorship or institutional overlap; contact details to be matched by the
editorial office):

1. Jesse C. Bloom (Fred Hutch Cancer Center, Seattle) -- deep mutational
   scanning and quantitative evaluation of protein variant effect predictors.
2. Rohit Singh (Broad Institute of MIT and Harvard) -- computational
   interpretation of missense variation and clinical variant classification.
3. Susan Marqusee (University of California, Berkeley) -- high-throughput
   protein stability assays and biophysics of intrinsically disordered regions.
4. David T. Jones (University College Dublin) -- large-scale assessment of
   protein function and fitness prediction methods.

Sincerely,

[Author names and affiliations]
Corresponding author: [name, email]
