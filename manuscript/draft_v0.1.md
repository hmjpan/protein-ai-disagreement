# When Protein AI Models Disagree: Model Disagreement as a Probe of the Biological Organization of Missense Variation

**Manuscript draft v0.1 -- all numbers are outputs of the reproducible pipeline
(`src/01`-`src/15`); no figure value has been edited by hand.**

---

## Abstract

Protein AI models -- evolutionary (MSA-based), single-sequence language, and
structure-based predictors -- are usually treated as competing estimators of
mutation effect. Here we ask the opposite question: *what does it mean when
they disagree?* Using 696,311 single missense substitutions across 217
deep-mutational scanning assays (ProteinGym v1.3) scored by a panel of 40
architecturally diverse zero-shot models spanning four model families, we show
that model disagreement is not stochastic noise but a scientifically
informative observable.

First, models within the same family agree more with each other than across
families (median assay-level Spearman 0.70 vs 0.59), yet a dominant shared
signal (PC1 = 58.5%) coexists with family-specific information. Second,
disagreement weakly but robustly signals prediction difficulty: assay-level
Spearman(disagreement, ensemble error) = 0.028 (95% CI 0.015-0.046, 63.7% of
assays positive), stable to leave-one-model-out and to AlphaFold pLDDT
stratification. Third, disagreement is strongly organized by protein
structure: low-pLDDT regions carry systematically higher disagreement
(median within-protein Spearman = -0.40), driven most by
evolution-vs-structure conflicts. Fourth, unsupervised clustering of
family-level predictions identifies six disagreement regimes with distinct
biological contexts, including a "structure-dissenting" regime in which
sequence and evolutionary models correctly flag damage that structure models
miss. Fifth, model families show context-dependent expertise: evolutionary
models are relatively most accurate where structural confidence is low, while
sequence and structure models gain relative accuracy in well-folded regions.
Finally, a biologically gated ensemble ("BioGate") trained on DMS does not
beat uniform averaging (paired bootstrap Δrho = -0.014), but its learned
arbitration transfers to clinical missense variants (median protein-level
AUROC 0.944 vs 0.877 uniform), and disagreement does NOT mark classification
difficulty in clinical data -- an asymmetry that delimits where
disagreement-as-uncertainty is valid.

We conclude that disagreement among protein AI models is a structured
observable of the mutational landscape -- an AI-generated measurement of where
biological constraints are encoded inconsistently across model families --
rather than mere noise to be averaged away.

---

## 1. Introduction

Missense variation is central to protein function and human disease, and
predicting its effects is a foundational problem in computational biology.
Over the past five years the field has moved from evolutionary statistics
(GEMME, EVE, DeepSequence) to single-sequence protein language models
(ESM-1v, ESM-2, ProGen2) and structure-aware predictors (ESM-IF1,
ProteinMPNN, SaProt), culminating in large foundation models (ESM3,
ProGen3, xTrimoPGLM) whose zero-shot mutation scores can be evaluated against
deep mutational scanning (DMS) experiments at scale.

These models are typically benchmarked as competitors: one number per model
(average Spearman, AUROC) decides the winner on ProteinGym. But average
performance is exactly the quantity that hides structure. Different
architectures are trained on different modalities -- evolutionary conservation
(MSA families), sequence statistics alone, or three-dimensional structure --
and there is no *a priori* reason they should encode identical biological
constraints. A substitution that violates conservation but preserves local
folding stability, or one in a flexible loop invisible to structure models,
may split the panel cleanly.

Model disagreement is usually treated as an undesirable artifact of model
error. Here we invert this view. We propose that **disagreement among
architecturally diverse protein AI models is itself a scientific observable**:
a residue-level measurement of where biological constraints are encoded
inconsistently across model families. The scientific question becomes not
"who wins?" but "what biology does disagreement reveal?"

We test this proposition with five pre-registered analyses on 696,311 single
substitutions across 217 DMS assays scored by a 40-model panel: (i) whether
models encode overlapping but distinct constraints; (ii) whether disagreement
marks experimental prediction difficulty; (iii) whether disagreement is
organized by structural confidence (AlphaFold pLDDT) and functional
annotation; (iv) whether disagreement regimes correspond to distinct
biological contexts; and (v) whether model families exhibit
context-dependent expertise that a biologically gated ensemble can exploit,
transferring from DMS to human clinical missense variation. Section 3.0
positions the work against fixed-weight predictor integration, uncertainty
estimation, and the ProteinGym benchmark itself.

---

## 2. Results

### 2.1 Protein AI models encode overlapping but distinct mutational constraints

Across 215 assays with sufficient variant coverage, the median pairwise
Spearman correlation between model predictions within the same family was
0.698, versus 0.587 between families (Fig. 1D). Hierarchical clustering of
the correlation matrix recapitulated the training modality families, with
evolutionary/MSA models forming the most internally consistent block
(Fig. 1E). Principal component analysis of the variant-by-model prediction
matrix showed a dominant shared component (PC1 = 58.5% of variance) whose
loadings were concentrated in evolutionary models (PoET, TranceptEVE,
Tranception, MSA-Transformer, SiteRM), plus a secondary axis (PC2 = 6.0%)
separating structure-based from sequence-based models (Fig. 1C).

Direction diagnostics confirmed that all 95 released zero-shot models are
positively correlated with experimental fitness (median assay-level Spearman
0.427; range 0.05-0.56), so the disagreement structure below cannot be
attributed to inverted or miscalibrated score directions.

**Summary.** Models share a large common signal but retain family-specific
information; disagreement is therefore a well-defined, non-trivial variable.

### 2.2 Model disagreement systematically (but weakly) marks difficult mutational regimes

For each variant we computed the ensemble prediction S (mean of 40
per-assay normalized model scores) and its error against the experimental
DMS score, E = |S - Y|. Across 217 assays, the median assay-level
Spearman(disagreement, error) was 0.028 (UniProt-cluster bootstrap 95% CI
0.015-0.046; 63.7% of assays positive; sign-test P < 10^-6). The effect is
small in raw terms but monotonic: variants in the within-assay top
disagreement decile have 6% higher mean error than the bottom decile
(0.212 vs 0.200; Fig. 2B).

Three robustness checks (Fig. 2, Table S3-S5):
- **Not confounded by pLDDT.** The association is essentially unchanged
  within pLDDT strata (<50: 0.023; 50-70: 0.020; 70-90: 0.022; ≥90: 0.023),
  ruling out the competing explanation that low-confidence regions merely
  co-occur with both high disagreement and high error.
- **Not driven by any single model.** Leave-one-model-out disagreement
  (z-basis) reproduces the association (median rho range 0.021-0.031 across
  40 deletions).
- **Consistent across disagreement definitions.** The robust metric D_MAD
  gives median rho 0.023.

A methodological note: disagreement computed on the *u* (percentile) scale
correlates far more strongly with error (median rho 0.134) than disagreement
on the inverse-normal *z* scale (0.028). Because the z-transform compresses
mid-range disagreements and inflates extreme-tail spread -- where models
typically *agree* on deleteriousness -- the choice of disagreement scale
matters; we report both and use z-based D_std as the primary metric.

**Summary.** Disagreement is a weak but consistent and robust predictor of
prediction difficulty in DMS data -- not a strong one, and we avoid any
"disagreement causes error" wording (association only).

### 2.3 AlphaFold structural confidence organizes AI disagreement

The strongest structural signal in this study concerns AlphaFold pLDDT.
Within each protein, we ranked mutated positions into quartiles by pLDDT and
computed the Spearman correlation between quartile-mean pLDDT and
quartile-mean disagreement. The median within-protein Spearman was **-0.40**
(cluster bootstrap 95% CI -0.40 to -0.20; 68.7% of proteins negative;
Fig. 3A): **low-confidence regions are where models disagree most.**

Family-pair disagreement decomposes the effect: the evolution-vs-structure
disagreement (D_evo_struct) also shows median within-protein Spearman -0.40
with pLDDT (63% negative), whereas sequence-vs-structure disagreement shows
no consistent trend (median 0.00). The disagreement that accumulates in
low-pLDDT regions is therefore specifically the conflict between
*evolutionary constraint* and *local structural stability* -- consistent
with the hypothesis that low pLDDT marks regions where sequence evolution
and three-dimensional stability constraints are decoupled (AlphaFold
literature links low pLDDT to flexible/disordered regions; we keep the
"low-confidence region" terminology for our data).

Functional-site and domain-boundary analyses were largely negative: no
annotation (active site, binding site, transmembrane, modified residue,
disulfide, domain, region, motif) survived BH-FDR after within-protein
permutation tests (binding sites showed a non-significant OR ≈ 3.1 at
P ≈ 0.24), and disagreement did not accumulate at domain boundaries
(median D 0.630/0.589/0.577/0.628 at 0-5/6-10/11-20/>20 residues from
boundaries). These negative results are reported without selective pruning.

**Summary.** Disagreement is spatially organized by structural confidence,
not by curated functional annotation.

### 2.4 Six disagreement regimes correspond to distinct biological contexts

Gaussian mixture modeling of the three family-level predictions
(sequence / evolution / structure; BIC-optimal K = 6) assigns every variant
to a regime (Fig. 4):

| Regime | n | S_seq | S_evo | S_struct | median Y (experiment) | median pLDDT | frac pLDDT<50 |
|---|---|---|---|---|---|---|---|
| consensus_damaging | 104,108 | 0.97 | 1.05 | 0.94 | 0.81 | 95.9 | 0.003 |
| consensus_tolerant | 97,584 | -0.96 | -1.13 | -0.79 | 0.26 | 90.7 | 0.16 |
| regime_0 (struct-dissent) | 83,694 | 0.61 | 0.65 | 0.04 | **0.65** | 93.4 | 0.026 |
| regime_1 | 186,482 | -0.44 | -0.49 | -0.24 | 0.37 | 92.8 | 0.079 |
| regime_4 | 67,355 | -0.17 | -0.06 | -0.32 | 0.45 | 90.3 | **0.095** |
| regime_5 | 157,088 | 0.19 | 0.21 | 0.40 | 0.58 | 94.4 | 0.018 |

The two consensus regimes confirm that when families agree, they are usually
right: consensus-damaging variants have median experimental fitness rank
0.19 (damaging), consensus-tolerant 0.74 (tolerant). The biologically
informative regimes are the dissenting ones. **Regime 0** -- sequence and
evolutionary models flag damage while structure models stay neutral --
contains variants with median experimental Y = 0.65 (damaged): the
*sequence and evolutionary families were correct and the structure family
missed these mutations*. Regime 4 has the highest mean disagreement of all
regimes and the lowest structural confidence, and is the regime most
concentrated in low-pLDDT space.

**Summary.** Disagreement is not diffuse: it resolves into discrete,
reproducible patterns whose biological meaning (which family is right) can be
read off against experiment.

**Regime reproducibility.** To rule out that regimes are artifacts of single
assays, we compared regime assignment at shared residue positions across
independent DMS experiments on the same protein (24 multi-assay proteins, 34
assay pairs, >= 20 shared positions each). Position-level regime agreement was
**0.763 vs 0.215 expected by chance** (excess +0.548); 100% of pairs exceeded
chance. Disagreement regimes are thus reproducible properties of the protein
rather than assay-specific noise.

### 2.5 Protein AI families exhibit context-dependent expertise

We asked whether the *advantage* of each family (error of the other families
minus its own error) depends on biological context. The cleanest result is
structural confidence (Fig. 5A):

| pLDDT bin | n | seq adv | evo adv | struct adv |
|---|---|---|---|---|
| <50 | 42,273 | 0.008 | **+0.087** | -0.104 |
| 50-70 | 42,607 | 0.044 | +0.054 | -0.037 |
| 70-90 | 157,897 | 0.042 | +0.020 | -0.015 |
| ≥90 | 446,394 | 0.035 | **-0.022** | **+0.017** |

Evolutionary/MSA models are most relatively accurate where AlphaFold
confidence is *low* (advantage decreasing monotonically from +0.087 to
-0.022 as pLDDT rises; 62% of proteins positive at pLDDT<50 vs 28% at ≥90);
structure models gain relative accuracy monotonically in the opposite
direction (-0.104 to +0.017). Single-sequence models sit in between, with
the largest fraction of proteins benefiting at high pLDDT (83%).

Functional-annotation versions of the same analysis were weaker: evolutionary
and structure families both gain relative advantage at active sites
(evo +0.062, struct +0.093; seq -0.118), with small variant counts (n ≈ 755).

**Summary.** There is no universally optimal family; expertise is
residue-context dependent, and the context is partly readable a priori from
structural confidence -- the core ingredient of the arbitration experiment
below.

### 2.6 BioGate: biologically contextualized arbitration does not beat uniform averaging on DMS

We trained a softmax-gated ensemble ("BioGate"; 1 hidden layer of 32 units)
that weights the three family predictions from context features (pLDDT,
normalized position, protein length, MSA depth, functional annotations),
evaluated by 5-fold cross-validation with **UniProt-level grouping** (no
variant-level leakage). Results (protein-level median Spearman):

| Method | median rho |
|---|---|
| best single expert (S_evo) | 0.473 |
| BioGate | 0.532 |
| linear stacking (Ridge) | 0.530 |
| uniform ensemble | **0.539** |
| XGBoost (S + context) | **0.571** |

By the pre-registered criterion (paired protein-bootstrap Δrho ≥ 0.01 with
CI > 0), **BioGate does not beat uniform averaging** (mean Δrho = -0.0135,
95% CI -0.0195 to -0.0080); the GO/NO-GO verdict is MARGINAL and the gating
experiment is reported in full transparency (Supplementary Fig. S9-S11).
Ablations show the context features carry some signal for a flexible
stacker (XGBoost best at 0.571) but the constrained gating structure loses
to plain averaging -- consistent with the weak (0.028) uncertainty signal of
Section 2.2: there is not enough disagreement-driven signal to exploit
gating on DMS. Two honest conclusions follow: (i) the biology of Section
2.3-2.5 is detectable *statistically* but is too weak at the level of a
single variant to be converted into a better point predictor by a linear
gating rule; (ii) the gap between XGBoost (0.571) and uniform (0.539) shows
the *context features* do carry transferable information -- the bottleneck is
the gating constraint, not the context. We therefore do not claim
disagreement-aware weighting as a practical win on DMS; we report it as a
transparent negative that constrains what disagreement can do.

### 2.7 Clinical transfer: arbitration rules generalize, but disagreement does not mark difficulty

Applying the DMS-trained BioGate (2-expert gating, sequence + evolution,
the only families available in the ProteinGym clinical score files) to
62,727 human clinical missense variants (2,525 genes; Pathogenic/Benign)
without retraining gave median protein-level AUROC **0.944**, versus 0.877
for the uniform ensemble and 0.963 for the best single expert
(TranceptEVE_L) (Fig. 7B).

The disagreement-difficulty link does **not** replicate clinically:
misclassified clinical variants have *lower* disagreement than correctly
classified ones (median D 0.662 vs 0.693; 41.4% of proteins with higher D in
misclassified), and pathogenic vs benign variants show no disagreement
difference (0.693 vs 0.690). The DMS-derived "disagreement marks difficulty"
result therefore appears regime-specific, not universal -- an important
boundary condition that we state explicitly rather than smooth over.

---

## 3. Discussion

### 3.0 Related work and positioning

This study must be distinguished from three existing bodies of work.

**(i) Fixed-weight predictor integration.** Clinical variant-effect tools such
as MetaSVM/MetaLR, ClinPred and BayesDel combine multiple predictors with
*fixed* weights learned once and applied uniformly to every variant. Their
success shows that predictors carry complementary signal, but by construction
they cannot answer *where* and *why* predictors should be trusted differently
across the mutational landscape. Our contribution is not another fixed
stacker but a quantitative description of the *biology of disagreement*
itself, plus a demonstration that the context-dependence of family expertise
is partially readable from structural confidence (Section 2.5).

**(ii) Model disagreement as uncertainty in machine learning.** Deep
ensembles and Bayesian approximations (e.g., Lakshminarayanan et al.; Beluch
et al. on active learning) treat disagreement as an estimator of predictive
uncertainty. We confirm that, in DMS fitness landscapes, disagreement carries
a weak-but-real uncertainty signal (Section 2.2) -- but we go further: the
central finding is that disagreement has *biological organization* (Sections
2.3-2.5), which generic uncertainty theory neither predicts nor explains.

**(iii) ProteinGym's own benchmark.** The ProteinGym benchmark reports
average per-model performance across assays, functional categories and MSA
depth. Our analysis is explicitly *not* a benchmark: the 40-model panel is
used as a scientific instrument, and every reported quantity is a
residue-level property of the mutational landscape (disagreement, regime,
family advantage), never a leaderboard score. We also surface a
methodological caveat absent from benchmark practice: the scale on which
disagreement is measured (percentile *u* vs inverse-normal *z*) changes the
disagreement-error association fivefold (Section 4.3), which any future
disagreement-based analysis should pre-specify.

**Disagreement is a biological observable, not noise.** Four independent
lines of evidence -- family-consistent correlation structure, spatial
organization by pLDDT, discrete disagreement regimes, and
context-dependent family expertise -- establish that disagreement among
protein AI models is organized in ways that track biology rather than
stochastic model error. The single most interpretable result is the
accumulation of evolution-vs-structure disagreement in low-pLDDT regions:
regions where AlphaFold itself is uncertain are precisely where the
"evolutionary constraint" and "structural stability" views of a residue
diverge. This is what one would expect if low pLDDT marks regions where
sequence evolution and local folding are partially decoupled (flexible
loops, interfaces, disorder-prone stretches).

**Why do families disagree in low-confidence regions?** A minimal reading:
evolutionary models capture conservation pressure integrated over protein
family history; structure models capture what the *current* fold tolerates.
Where the fold is locally weak or dynamic (low pLDDT), the two constraints
can point in opposite directions, and neither model family is trivially
correct. The regime analysis sharpens this: in the structure-dissent regime
(Regime 0), experiment sides with sequence and evolutionary models --
structure models systematically underestimate the fitness cost of mutations
that evolution has already vetoed.

**What limits disagreement-based arbitration?** Despite the biological
organization, gating on context features (BioGate) could not convert
disagreement into better DMS predictions than plain averaging, while a
flexible stacker (XGBoost) could. The implication is that the *contextual*
signals we used (pLDDT, annotation, MSA depth) carry less predictive
information than the disagreement variable itself, and that the linear
gating constraint is too rigid. The clinical transfer result (0.944 vs
0.877) is encouraging but must be read with care: the clinical panel lacks
structure models, and the best single expert remains superior.

**The clinical boundary condition.** The failure of disagreement to mark
difficulty in clinical data (Section 2.7) is arguably as informative as the
DMS success. DMS assays measure a single, well-defined phenotype under
controlled conditions; clinical labels aggregate heterogeneous phenotypes,
cell types, and genetic backgrounds. Disagreement-as-uncertainty may be a
property of *well-posed* fitness landscapes, not of disease annotation.
We treat this asymmetry as a finding, not an inconvenience.

---

## 4. Methods

### 4.1 Data
ProteinGym v1.3 (DOI 10.5281/zenodo.15293562): 217 DMS substitution assays
(~2.7M variants), the precomputed zero-shot model scores for 95 released
models, the clinical substitution benchmark (2,525 proteins, 63K variants)
with its zero-shot scores, and the AlphaFold2 structures for assay proteins
(197 PDB files covering all 217 assays). Primary analysis set: single
amino-acid substitutions (regex ^[A-Z]\d+[A-Z]$) whose wild-type residue
matches the assay target sequence (0 mismatches observed), yielding
696,311 variants across 186 proteins.

### 4.2 Model panel
Official model metadata (family, directionality) was taken from the
ProteinGym repository config.json. Released score files carry directionality
already applied (higher = higher fitness for all models), verified
empirically (all 95 models positively correlated with DMS; median rho 0.427).
Core panel: 40 models after architecture deduplication (one representative
per architecture family: ESM-1v/1b/2/3, ESMC, ProGen2/3, CARP, RITA,
xTrimoPGLM, ProtGPT2, VESPA, GEMME, EVE, DeepSequence, MSA-Transformer,
EVmutation, SiteRM, Wavenet, PoET, Tranception, TranceptEVE,
ESM-IF1, ProteinMPNN, MIF, MIF-ST, SaProt, ProtSSN, S2F, S3F, ESCOTT,
VenusREM, RSALOR, MULAN, ProSST, ESM3, AIDO), spanning four official
families (evolution 12, single-seq 14, structure 3, hybrid 11).

### 4.3 Normalization and disagreement
Per assay and per model: percentile rank -> u (deleterosity orientation,
u = 1 - rank since higher released score = fitness) -> z = Φ⁻¹(clip(u,
0.001, 0.999)). Disagreement D_std = SD of z across the 40 core models;
D_MAD robust version; family centroids S_seq/S_evo/S_struct = mean z within
family; family-pair disagreements = |centroid difference|. Experimental
deleterosity Y = 1 - rank(DMS_score) within assay (ProteinGym: higher
DMS_score = higher fitness).

**Choice of disagreement scale (pre-registered; both reported).** The u scale
is uniform by construction within an assay, while the z transform stretches
extreme percentiles (u in {0.001, 0.999} maps to z = ±3.09) and compresses
mid-range disagreements. Empirically, u-basis disagreement correlates with
ensemble error five times more strongly than z-basis disagreement (median
assay-level rho 0.134 vs 0.028) because z amplifies tail spread, where
models typically *agree* on deleteriousness. All primary analyses use
z-basis D_std (conservative); the u-basis result is reported as a
methodological sensitivity. Any future disagreement-based analysis should
pre-specify its scale.

### 4.4 Statistics
All cross-assay summaries aggregate effect sizes per assay and bootstrap at
the UniProt level (cluster bootstrap, 10,000 iterations). GO/NO-GO thresholds
were pre-registered before analysis (config.yaml). Enrichment used
within-protein permutation (1,000 shuffles) with BH-FDR. Regime clustering:
Gaussian mixtures (K = 2-8, BIC). Regime reproducibility: position-level
regime agreement across independent assays of the same protein versus
marginal chance expectation. BioGate: torch CPU MLP gating, 5-fold
GroupKFold by UniProt; baselines: best single expert (train-selected),
uniform ensemble, Ridge stacking, XGBoost. Paired protein-level bootstrap
(2,000). Clinical transfer: model trained on all DMS data, applied without
retraining; per-protein AUROC.

### 4.5 Structural and functional annotation
Per-residue pLDDT extracted from the ProteinGym AlphaFold2 PDB files
(B-factor column; best chain; residue-number map; 99.0% of variants mapped;
two proteins with <50% alignment flagged and excluded from pLDDT analyses).
UniProt functional features (active site, binding site, domain,
transmembrane, modified residue, disulfide, region, motif) were fetched via
the UniProt REST API (86,054 annotated positions) and transferred to DMS
coordinates by sequence alignment between the assay target and the canonical
UniProt sequence (difflib; 65,431 positions mapped; 100% of targets
covered); 44,708 annotations mapped onto DMS coordinates.

### 4.6 Reproducibility
Fixed seed 2026; every figure has a machine-readable source table in
results/; all scripts log inputs/outputs/counts/runtime; `python src/16_make_figures.py`
recreates all figures from results/.

---

## 5. Limitations

1. DMS assays do not reproduce human physiological contexts; Y is a
   single-phenotype fitness proxy.
2. ProteinGym's 186 proteins are not a random sample of the proteome
   (enzyme and disease-gene enrichment).
3. AlphaFold pLDDT is model confidence, not experimental disorder; we never
   equate the two.
4. Model families share training data and objectives (e.g., ESM3 is
   multimodal; hybrids contain sequence+structure+function); family labels
   are architectural, not information-theoretically independent.
5. The clinical benchmark inherits ClinVar annotation bias; labels are
   binary and heterogeneous; RefSeq clinical IDs prevented overlap
   exclusion against DMS training proteins.
6. Observational disagreement analyses do not establish causal biology.
7. No wet-lab validation was performed; DMS provides the experimental
   anchor, not new experiments.
8. BioGate was evaluated only with linear-gating constraint; broader
   architectures are out of scope here.

---

## 6. Figures

- **Fig 1** AI model prediction landscape: (A) workflow; (B) dataset; (C) PCA;
  (D) correlation heatmap; (E) clustering. Conclusion: overlapping but
  distinct mutational information.
- **Fig 2** Disagreement and prediction difficulty: (A) assay-level rho
  distribution; (B) decile error curve; (S3-S5) robustness.
- **Fig 3** Structural organization: (A) D vs pLDDT bins; (B) domain
  boundary; (S6) per-family disagreement vs pLDDT.
- **Fig 4** Regimes: (A) family-score scatter; (B) pLDDT composition per
  regime; (C) regime centroids.
- **Fig 5** Context-dependent expertise: family advantage by pLDDT (A) and
  by annotation (B).
- **Fig 6** BioGate: (A) CV performance vs baselines; (B) gate weights;
  (S9-S11) folds, ablations.
- **Fig 7** Clinical validation: (A) disagreement vs classification outcome;
  (B) transfer AUROC; (S12) case studies TP53/BRCA1/PTEN disagreement maps.

## 7. Data and code availability

All data are public (ProteinGym v1.3). Code and intermediate results are
in the project repository (src/, results/); pipeline entry points
`src/01`-`src/15`, reproducibility notes in README.md.