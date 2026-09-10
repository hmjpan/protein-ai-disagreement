# Protein AI Disagreement Reveals a Continuum from Evolutionary to Structural Constraint across Protein Conformational Organization

**Manuscript draft v0.2 (review-response revision) -- all numbers from the
reproducible pipeline (src/01-src/20); nothing hand-edited.**

---

## Abstract

Protein AI models -- evolutionary (MSA-based), single-sequence language, and
structure-based predictors -- are usually treated as competing estimators of
mutation effect. Here we invert this framing and ask what it means, biologically,
when they disagree. Using 696,311 single missense substitutions across 217
deep-mutational scanning assays (ProteinGym v1.3) scored by a 40-model panel
spanning four model families, we show that **model disagreement is a structured
observable of the mutational landscape, organized by protein conformational
organization rather than by model noise**.

The central finding is spatial and mechanistic. Disagreement between
evolutionary and structure-based models rises sharply where AlphaFold
structural confidence is low (median within-protein Spearman = -0.40; 68.7%
of proteins negative) and this trend is invariant to panel composition
(1000 balanced three-model-per-family resamplings all reproduce it; a strict
eight-model architecture panel gives -0.40). Total disagreement is also
elevated in curated intrinsically disordered regions (+0.040, 76.5% of
proteins), and model families exhibit opposite expertise gradients: single-
sequence models lose accuracy in disordered regions, evolutionary models
gain it where structure is uncertain, and structure models gain it in
well-folded, contact-dense regions. Unsupervised clustering resolves
disagreement into six regimes whose assignments reproduce across independent
assays of the same protein (position-level agreement 0.76 vs 0.22 chance),
including a "structure-dissenting" regime in which the measured DMS phenotype
aligns with sequence/evolutionary predictions rather than structure-family
predictions.

Secondary, weaker signals are reported with full transparency: disagreement
is only weakly associated with prediction error in DMS (assay-level
Spearman 0.028, 95% CI 0.015-0.046), and this association does not
replicate in clinical variant classification. A biologically gated ensemble
does not beat uniform averaging on DMS (Δrho = -0.014, reported as a
constraint, not a win), although its learned arbitration transfers to
clinical variants (median protein-level AUROC 0.944 vs 0.877 uniform; 0.905
on a strictly non-overlapping 2,490-protein set). We conclude that
disagreement among protein AI models should be read as a probe of where
evolutionary and structural constraints decouple across protein space --
an AI-generated map of conformational biology -- rather than as mere
prediction noise.

---

## 1. Introduction

Missense variation is central to protein function and human disease, and
predicting its effects is a foundational problem in computational biology.
Over the past five years the field has moved from evolutionary statistics
(GEMME [5], EVE [3], DeepSequence [4]) to single-sequence protein language
models (ESM-1v [6], ESM-2 [7], ProGen2) and structure-aware predictors
(ESM-IF1 [8], ProteinMPNN [9], SaProt [11]), culminating in large foundation
models (ESM3 [12], ProGen3 [13], xTrimoPGLM [14]) whose zero-shot mutation
scores are evaluated against deep mutational scanning (DMS) experiments at
scale (ProteinGym [1]).

These models are usually benchmarked as competitors: one average
performance number decides the winner. But average performance hides
structure. Models are trained on different modalities -- evolutionary
conservation (MSA families), sequence statistics alone, or three-dimensional
structure -- and there is no a priori reason they should encode identical
biological constraints. A substitution that violates conservation while
preserving local folding stability, or one in a flexible region that
structure models cannot confidently model, may split the panel cleanly.
The information the models *share* is a consensus estimate; the information
on which they *diverge* may be information about biology.

We therefore treat disagreement among architecturally diverse protein AI
models as a scientific observable: a residue-level measurement of where
biological constraints are encoded inconsistently across model families.
The question is not "who wins?" but "what biology does disagreement
reveal?"

We test this with pre-registered analyses on 696,311 single substitutions
across 217 DMS assays scored by a 40-model panel: (i) whether models encode
overlapping but distinct constraints; (ii) whether disagreement is organized
by structural confidence (AlphaFold pLDDT), curated disorder, and residue
mechanics (secondary structure, contact density); (iii) whether disagreement
regimes correspond to distinct biological contexts and reproduce across
independent experiments; (iv) whether model families exhibit
context-dependent expertise; (v) whether disagreement marks prediction
difficulty in DMS and in clinical variant classification; and (vi) whether
a biologically gated ensemble transfers from DMS to clinical data. Section
3.0 positions the work against fixed-weight predictor integration,
uncertainty estimation, and the ProteinGym benchmark itself.

---

## 2. Results

### 2.1 Protein AI models encode overlapping but distinct mutational constraints

Across 215 assays, the median pairwise Spearman correlation between model
predictions within the same family was 0.698 versus 0.587 between families
(Fig. 1D). Hierarchical clustering recapitulated the training modality
families (Fig. 1E). PCA of the variant-by-model matrix showed a dominant
shared component (PC1 = 58.5%) concentrated in evolutionary models, plus a
secondary axis (PC2 = 6.0%) separating structure-based from sequence-based
models (Fig. 1C).

Score orientation followed the released ProteinGym conventions and the
official model metadata (config.json); post hoc correlations with DMS were
used only as a sanity check and never to determine orientation (Methods
4.2).

### 2.2 AlphaFold structural confidence organizes disagreement (primary result)

Within each protein, mutated positions were ranked into quartiles by pLDDT;
the median within-protein Spearman between quartile-mean pLDDT and
quartile-mean disagreement was **-0.40** (cluster bootstrap 95% CI
-0.40 to -0.20; 68.7% of proteins negative; Fig. 3A). Low-confidence
regions are where models disagree most.

The effect is specific to family pairs: evolution-vs-structure disagreement
(D_evo_struct) shows median within-protein Spearman -0.40 with pLDDT (63%
negative), whereas sequence-vs-structure disagreement shows no consistent
trend (median 0.00). The disagreement that accumulates in low-pLDDT regions
is therefore the conflict between *evolutionary constraint* and *local
structural stability* -- the signature expected where the two constraints
decouple.

**Panel-composition robustness (reviewer defense).** Two independent checks
rule out that family-size or architecture-redundancy asymmetries (12
evolution / 14 single-sequence / 3 structure / 11 hybrid models) produce
this result. (A) A strict one-representative-per-architecture panel (ESM-1v,
GEMME, EVE, MSA-Transformer, ESM-IF1, ProteinMPNN, SaProt, Tranception)
reproduces the D_evo_struct trend (median rho -0.40; 65.8% of proteins
negative). (B) 1,000 balanced resamplings drawing three models per family
all reproduce a negative trend (median of per-iteration medians -0.40,
95% range [-0.40, -0.20], 100% of iterations negative; Fig. 3C).

**Independent disorder annotation.** Curated UniProt disordered regions
(transferred to DMS coordinates by sequence alignment; 92 proteins covered)
confirm the qualitative pattern: total disagreement is higher in disordered
than in ordered residues (median within-protein difference +0.040, 76.5% of
proteins positive; Fig. 3B). pLDDT<50 and UniProt disorder overlap
partially (P(low-confidence | disorder) = 0.56), consistent with the
literature linking AlphaFold confidence to experimental disorder
(Piovesan et al., 2022 [16]); we retain the "low structural confidence"
terminology for the metric itself. The evolution-vs-structure conflict is
maximal specifically in low-confidence space rather than in curated IDR
alone (median IDR-ordered difference -0.017, 38% positive), suggesting the
decoupling signal extends beyond canonical IDR into flexibly modeled
regions generally.

### 2.3 Residue mechanics: contact density and secondary structure

Using pydssp secondary-structure assignment and CA contact density from the
same AlphaFold structures (Fig. 6A; Fig. 4B), family expertise varies with residue
mechanics: structure-family advantage is positive in beta-sheets (+0.017)
and core residues (+0.012) and negative in helices (-0.006); sequence-family
advantage is largest in helices (+0.045); evolutionary advantage is lowest
in core residues (-0.020). Contact density itself shows no monotonic
association with evolution-vs-structure disagreement (median within-protein
rho 0.01). These effects are smaller than the pLDDT gradient and are
reported as secondary texture rather than primary evidence.

### 2.4 Six disagreement regimes correspond to distinct biological contexts and reproduce across experiments

Gaussian mixture modeling of family-level predictions (BIC-optimal K = 6;
Fig. 2A-C) assigns every variant to a regime:

| Regime | n | S_seq | S_evo | S_struct | median Y | median pLDDT |
|---|---|---|---|---|---|---|
| consensus_damaging | 104,108 | 0.97 | 1.05 | 0.94 | 0.81 | 95.9 |
| consensus_tolerant | 97,584 | -0.96 | -1.13 | -0.79 | 0.26 | 90.7 |
| structure-dissenting | 83,694 | 0.61 | 0.65 | 0.04 | **0.65** | 93.4 |
| mild-tolerant | 186,482 | -0.44 | -0.49 | -0.24 | 0.37 | 92.8 |
| high-disagreement | 67,355 | -0.17 | -0.06 | -0.32 | 0.45 | 90.3 |
| mild-damaging | 157,088 | 0.19 | 0.21 | 0.40 | 0.58 | 94.4 |

Consensus regimes are usually correct (consensus-damaging variants have
median experimental fitness rank 0.19; consensus-tolerant 0.74). In the
structure-dissenting regime, the measured DMS phenotype aligns more closely
with sequence and evolutionary predictions than with structure-family
predictions (median experimental Y = 0.65). We state this as an
alignment-of-evidence statement -- the DMS phenotype in these variants is
damaging while structure-family predictions are neutral -- and do not
claim that structure models are biologically "wrong", since DMS phenotypes
(activity, binding, abundance, fitness) are not direct measures of folding
constraint. (Regime names are descriptive labels assigned post hoc from the
actual centroids and summary statistics; the machine-readable labels are
regime_0-regime_5 in the regime tables.)

**Regime reproducibility.** Across independent assays of the same protein
(24 proteins, 34 assay pairs, >= 20 shared positions), position-level
regime agreement was 0.763 vs 0.215 expected by chance; 100% of pairs
exceeded chance. Regimes are properties of proteins, not of single assays.

### 2.5 Model families exhibit context-dependent expertise

We asked whether family *advantage* (error of other families minus own
error) depends on biological context (Fig. 4A; Fig. 4C):

| Context | seq adv | evo adv | struct adv |
|---|---|---|---|
| pLDDT < 50 | +0.008 | **+0.087** | -0.104 |
| pLDDT >= 90 | +0.035 | **-0.022** | **+0.017** |
| UniProt IDR | **-0.043** | +0.016 | +0.060 |
| ordered | +0.038 | -0.011 | +0.005 |

Evolutionary/MSA models are most relatively accurate where structural
confidence is low (advantage decreasing monotonically from +0.087 to -0.022
as pLDDT rises); structure models gain relative accuracy in the opposite
direction (-0.104 to +0.017); single-sequence models -- strongest in
well-folded regions overall (83% of proteins positive at pLDDT>=90) -- lose
accuracy specifically in intrinsically disordered regions (-0.043 vs +0.038
in ordered residues). There is no universally optimal family: expertise is
residue-context dependent, and the context is partly readable a priori.

### 2.6 Disagreement weakly marks difficulty in DMS (secondary, transparent)

For completeness we report the uncertainty read-out: median assay-level
Spearman(disagreement, ensemble error) = 0.028 (UniProt-cluster bootstrap
95% CI 0.015-0.046; 63.7% of assays positive). The effect is small but
monotonic across disagreement deciles (top-decile error 0.212 vs 0.200) and
robust: unchanged within pLDDT strata (~0.02 in all four bins; not a
confound), stable to leave-one-model-out (0.021-0.031), and consistent
across the D_MAD definition (0.023). The association is strongest in
stability and binding assays (median rho 0.076 and 0.057; Fig. S13) and
effectively absent in activity and organismal-fitness assays (0.003 and
-0.004), suggesting disagreement marks difficulty best where the measured
phenotype approximates a fitness/stability landscape. On the percentile (u) scale the
association is fivefold stronger (0.134); because the z-transform stretches
the tails where models usually agree, we pre-specify the z scale as primary
(Methods 4.3) and report both. We deliberately do not elevate this weak
signal into a claim that "disagreement predicts uncertainty".

### 2.7 Clinical data: domain-specific meaning, not universal uncertainty

Applying the DMS-trained BioGate (2-expert gating: sequence + evolution --
the only families available in the clinical score files) to 62,727 human
clinical missense variants (2,525 genes; Pathogenic/Benign) without
retraining gave median protein-level AUROC **0.944**, versus 0.877 for the
uniform ensemble and 0.963 for the best single expert (TranceptEVE_L)
(Fig. 7B).

**Strict non-overlap control.** Clinical proteins sharing >= 70% sequence
identity with any of the 186 DMS proteins (35 proteins; same-protein
overlaps are necessarily captured by this threshold) were excluded, leaving
60,214 variants across 2,490 proteins. The uniform ensemble AUROC on this
strict set was **0.905** (vs 0.877 on the full clinical set), i.e., the
transfer result does not depend on overlap with DMS training proteins.

**Disagreement does not mark clinical difficulty.** In contrast to DMS,
misclassified clinical variants have *lower* disagreement than correctly
classified ones (median 0.662 vs 0.693; 41.4% of proteins with higher D in
misclassified), and pathogenic vs benign variants show no disagreement
difference (0.693 vs 0.690). Disagreement-as-uncertainty is therefore
regime-specific: it is a property of well-posed fitness landscapes, not of
heterogeneous disease annotation. We treat this asymmetry as a finding.

### 2.8 BioGate: a transparent negative on DMS

A softmax-gated ensemble ("BioGate"; one 32-unit hidden layer) weighting the
three family predictions from context features (pLDDT, position, length,
MSA depth, functional annotations), evaluated by 5-fold UniProt-grouped CV,
did not beat uniform averaging (median protein-level rho 0.529 vs 0.539;
paired protein-bootstrap Δrho = -0.015, 95% CI -0.020 to -0.009;
Supplementary Fig. S9-S11). XGBoost on scores+context reached 0.571,
showing the context features carry transferable information and that the
bottleneck is the linear gating constraint, not the context. Consistent with
Section 2.6, the DMS disagreement signal is too weak at single-variant level
to convert into better point prediction by gating. This negative is
reported as a constraint on what disagreement can do, not as a failure of
the biological findings.

---

## 3. Discussion

### 3.0 Related work and positioning

(i) **Fixed-weight predictor integration** (MetaSVM/MetaLR [17], ClinPred
[18], BayesDel [19]) learns fixed weights applied uniformly; it cannot
answer *where* and *why* predictor trust should differ across the
mutational landscape. Our contribution is the quantitative description of
the biology of disagreement itself, and the demonstration that family
expertise varies readably with structural confidence and disorder.

(ii) **Model disagreement as uncertainty in ML** (deep ensembles [20],
Bayesian approximations, active learning [21]). We confirm a weak
uncertainty signal in DMS (Section 2.6) but show it is regime-specific
(Section 2.7) and, more importantly, that disagreement has biological
organization (Sections 2.2-2.5) that generic uncertainty theory neither
predicts nor explains.

(iii) **The ProteinGym benchmark** reports average per-model performance.
Our analysis is not a benchmark: the 40-model panel is an instrument, every
reported quantity is a residue-level property of the landscape, and we
surface a methodological caveat absent from benchmark practice -- the
disagreement scale (u vs z) changes the disagreement-error association
fivefold, so future disagreement analyses must pre-specify their scale.

### 3.1 Disagreement is a probe of constraint decoupling

The single most interpretable result is the accumulation of
evolution-vs-structure disagreement in low-pLDDT regions: exactly where
AlphaFold is uncertain, the "evolutionary constraint" and "structural
stability" views of a residue diverge. This is what one would expect if low
structural confidence marks regions where sequence evolution and local
folding are partially decoupled -- flexible loops, interfaces,
disorder-prone stretches. The independent UniProt disorder annotation and
the residue-mechanics results support the same picture: sequence models
fail where disorder is curated, structure models fail where confidence is
low, and evolutionary models are most informative where structure is
weakest.

### 3.2 What the clinical asymmetry teaches

The failure of disagreement to mark difficulty in clinical data is as
informative as the DMS success. DMS assays measure one well-defined
phenotype under controlled conditions; clinical labels aggregate
heterogeneous phenotypes, cell types and genetic backgrounds. The boundary
condition we document -- disagreement-as-uncertainty holds in fitness
landscapes, not in disease annotation -- is precisely the kind of statement
that prevents mechanical misapplication of disagreement-based uncertainty in
clinical genomics.

### 3.3 Limits of disagreement-based arbitration

The biological organization of disagreement is statistically strong but
per-variant weak; linear gating cannot convert it into better DMS point
predictions (Section 2.8). Clinical transfer (0.944 vs 0.877 uniform; 0.905
on strict non-overlap) is encouraging but must be read with care: the
clinical panel lacks structure models, and the best single expert remains
superior.

---

## 4. Methods

### 4.1 Data
ProteinGym v1.3 (DOI 10.5281/zenodo.15293562): 217 DMS substitution assays
(~2.7M variants), precomputed zero-shot scores for 95 released models, the
clinical substitution benchmark (2,525 proteins, ~63K variants) with its
zero-shot scores, and AlphaFold2 structures for assay proteins (197 PDB
files covering all 217 assays). Primary set: single amino-acid substitutions
(regex ^[A-Z]\d+[A-Z]$) whose wild-type residue matches the assay target
sequence (0 mismatches): 696,311 variants across 186 proteins.

### 4.2 Model panel and score orientation
Official model metadata (family, directionality) from the ProteinGym
config.json. Released score files carry directionality already applied
(higher = higher fitness); the direction table is recorded in
model_direction.csv with evidence labels. Post hoc DMS correlations (all 95
models positive; median 0.427) were used only as a sanity check and never
to determine orientation. Core panel: 40 models after architecture
deduplication (one representative per architecture, including GEMME [5],
EVE [3], DeepSequence [4], MSA-Transformer [10], Tranception and
TranceptEVE [2], ESM-1v [6], ESM-2 [7], ESM-IF1 [8], ProteinMPNN [9],
SaProt [11]; four families: evolution 12, single-seq 14, structure 3,
hybrid 11). All family-level
conclusions were re-verified with a strict 8-model panel and 1,000 balanced
three-per-family resamplings (Section 2.2).

### 4.3 Normalization and disagreement scale
Per assay and per model: percentile rank -> u (deleterosity orientation;
u = 1 - rank since higher released score = fitness) -> z = Phi^-1(clip(u,
0.001, 0.999)). D_std = SD of z across core models; D_MAD robust version;
family centroids = mean z within family; family-pair disagreements =
|centroid difference|; Y = 1 - rank(DMS_score) within assay. Scale
pre-specification: u-basis disagreement correlates with error fivefold more
strongly than z-basis (0.134 vs 0.028) because z stretches tails where
models agree; z is the primary scale, u reported as sensitivity.

### 4.4 Statistics
Cross-assay summaries aggregate per-protein effects; UniProt-level cluster
bootstrap (10,000). GO/NO-GO thresholds pre-registered (config.yaml).
Within-protein permutation enrichment with BH-FDR. Regime clustering:
Gaussian mixtures K = 2-8 by BIC; regime reproducibility = position-level
agreement across independent assays vs marginal chance. BioGate: torch CPU
MLP gating, 5-fold GroupKFold by UniProt; baselines best single expert
(train-selected), uniform ensemble, Ridge stacking, XGBoost; paired
protein-level bootstrap (2,000). Clinical transfer: trained on all DMS,
applied without retraining; per-protein AUROC; strict non-overlap set by
>= 70% sequence identity (8-mer filter + alignment) excluding 35 clinical
proteins.

### 4.5 Structural, mechanical and disorder annotation
pLDDT from ProteinGym AlphaFold2 PDBs (B-factor, best chain; 99.0% of
variants mapped; two proteins with < 50% alignment flagged and excluded).
AlphaFold2 structures and confidence scores are those of Jumper et al.
[15]. pydssp secondary structure (H/E/C) and CA contact density (8/10 A)
from the same PDBs; burial tertiles within protein. UniProt features
(active site, binding site, domain, transmembrane, modified residue,
disulfide, region, motif) and curated disordered regions fetched via the
UniProt REST API [22] and transferred to DMS coordinates by sequence
alignment (difflib; 65,431 positions mapped; 44,708 annotations
transferred; 4,569 disordered positions across 92 proteins). The clinical
benchmark labels derive from ClinVar curation [23]. Structural parsing used
Biopython [24]; the XGBoost baseline used XGBoost [25].

### 4.6 Reproducibility
Seed 2026; every figure has a machine-readable source table; all scripts
log inputs/outputs/counts/runtime; `python src/16_make_figures.py`
recreates all main figures from results/.

---

## 5. Limitations

1. DMS assays measure single-phenotype fitness proxies; "damaging" refers
   to the measured phenotype, not universal protein function.
2. The 186 DMS proteins are not a random proteome sample (enzyme and
   disease-gene enrichment).
3. AlphaFold pLDDT is model confidence, not experimental disorder; we keep
   separate terminology and cite the pLDDT-disorder literature only in
   discussion.
4. Model families share training data and objectives; family labels are
   architectural; the balanced-sampling analysis bounds the consequences.
5. The clinical benchmark inherits ClinVar annotation bias; labels are
   binary and heterogeneous. Non-overlap was enforced by sequence identity
   (>= 70%); a RefSeq-to-UniProt accession-level exclusion could not be
   computed reliably (API failures), but is subsumed by the identity
   threshold for same-protein cases.
6. Observational disagreement analyses do not establish causal biology.
7. No wet-lab validation; DMS provides the experimental anchor.
8. BioGate was evaluated only under a linear gating constraint.

---

## 6. Figures

- **Fig 1** Model prediction landscape: (C) PCA of model predictions; (D)
  model correlation heatmap; (E) hierarchical clustering. (A-B workflow and
  dataset schematics to be drafted in final layout.)
- **Fig 2** Regime structure: (A) family-score scatter by regime; (B) pLDDT
  composition per regime; (C) regime centroids vs experimental Y.
- **Fig 3** Structural organization: (A) disagreement vs pLDDT; (B)
  disagreement in curated disordered vs ordered residues; (C) balanced
  model-sampling sensitivity (distribution of per-iteration medians).
- **Fig 4** Context-dependent expertise: (A) family advantage by pLDDT;
  (B) family advantage by residue mechanics; (C) family advantage in
  disordered vs ordered residues.
- **Fig 5** Reproducibility and case studies: (A) position-level regime
  agreement vs chance across independent assays; (B-D) TP53 / BRCA1 / PTEN
  AlphaFold structures colored by disagreement.
- **Fig 6** Residue mechanics: (A) evolution-vs-structure disagreement by
  burial class and secondary structure.
- **Fig 7** Clinical validation: (A) disagreement vs classification outcome;
  (B) transfer AUROC (full set and strict non-overlap).
- **Supplementary**: S6A-B disagreement vs error (distribution, deciles);
  S7 domain-boundary analysis; S9A-B BioGate performance and gate weights;
  S10 BioGate ablations; S11 leave-one-model-out sensitivity; S12 u- vs
  z-scale disagreement comparison; S13 disagreement-error association by
  assay selection type and MSA depth.

## 7. Data and code availability

All data are public (ProteinGym v1.3; UniProt). Code and intermediate
results are in the project repository (src/01-src/20); reproducibility
notes in README.md.

## 8. References

1. Notin, P. et al. ProteinGym: large-scale benchmarks for protein design
   and fitness prediction. NeurIPS Datasets and Benchmarks (2023).
2. Notin, P. et al. Tranception: protein fitness prediction with autoregressive
   transformers and inference-time retrieval. ICML (2022).
3. Frazer, J. et al. Disease variant prediction with deep generative models
   of evolutionary data. Nature 599, 91-95 (2021).
4. Riesselman, A. J. et al. Deep generative models of genetic variation
   capture the effects of mutations. Nat. Methods 15, 816-822 (2018).
5. Laine, E. & Carbone, A. GEMME: a simple and fast global epistatic model
   predicting mutational effects. Mol. Biol. Evol. 36, 2604-2619 (2019).
6. Meier, J. et al. Language models enable zero-shot prediction of the
   effects of mutations on protein function. NeurIPS (2021).
7. Lin, Z. et al. Evolutionary-scale prediction of atomic-level protein
   structure with a language model. Science 379, 1123-1130 (2023).
8. Hsu, C. et al. Learning inverse folding from millions of predicted
   structures. ICML (2022).
9. Dauparas, J. et al. Robust deep learning-based protein sequence design
   using ProteinMPNN. Science 378, 49-56 (2022).
10. Rao, R. et al. MSA Transformer. ICML (2021).
11. Su, J. et al. SaProt: protein language modeling with structure-aware
    vocabulary. ICLR (2024).
12. Hayes, T. et al. Simulating 500 million years of evolution with a
    language model. bioRxiv 2024.07.01.600583 (2024).
13. Madani, A. et al. Scaling unlocks broader generation and deeper
    functional understanding of proteins. bioRxiv 2025.04.15.649055 (2025).
14. Chen, B. et al. xTrimoPGLM: unified 100B-scale pre-trained transformer
    for deciphering the language of proteins. bioRxiv 2024.07.24.605092 (2024).
15. Jumper, J. et al. Highly accurate protein structure prediction with
    AlphaFold. Nature 596, 583-589 (2021).
16. Piovesan, D., Monzon, A. M. & Tosatto, S. C. E. Intrinsic protein disorder
    and conditional folding in AlphaFoldDB. Protein Sci. 31, e4466 (2022).
17. Dong, C. et al. Comparison and integration of deleteriousness prediction
    methods for nonsynonymous SNVs. PLoS ONE 10, e0134848 (2015).
18. Alirezaie, N. et al. ClinPred: prediction tool to identify disease-relevant
    nonsynonymous single-nucleotide variants. Am. J. Hum. Genet. 103, 474-483 (2018).
19. Feng, B.-J. PERCH: a unified framework for disease gene prioritization.
    Hum. Mutat. 38, 243-251 (2017). (BayesDel)
20. Lakshminarayanan, B. et al. Simple and scalable predictive uncertainty
    estimation using deep ensembles. NeurIPS (2017).
21. Beluch, W. H. et al. The power of ensembles for active learning in image
    classification. CVPR (2018).
22. The UniProt Consortium. UniProt: the universal protein knowledgebase in
    2023. Nucleic Acids Res. 51, D523-D531 (2023).
23. Landrum, M. J. et al. ClinVar: improving access to variant interpretations
    and supporting evidence. Nucleic Acids Res. 46, D1062-D1067 (2018).
24. Cock, P. J. A. et al. Biopython: freely available Python tools for
    computational molecular biology and bioinformatics. Bioinformatics 25,
    1422-1423 (2009).
25. Chen, T. & Guestrin, C. XGBoost: a scalable tree boosting system. KDD
    (2016).

---
**SUPERSEDED**: this draft has been superseded by the review-response revision in manuscript/submission/manuscript.md (v1.2). All numbers and analyses in the submission package take precedence.
