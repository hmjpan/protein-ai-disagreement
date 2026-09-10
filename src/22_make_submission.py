"""22_make_submission.py  (v1.1 -- review-response revision)

Assemble the submission package for AI for Science (IOP Publishing):
  manuscript/submission/manuscript.docx   (IOP-style, Word)
  manuscript/submission/manuscript.md
  manuscript/submission/cover_letter.md
  manuscript/submission/figures_list.md

Revision log (review response):
  v1.1  residue-level pLDDT-disagreement is the primary statistic
        (quartile curve demoted to visualization); five disagreement
        definitions cross-checked; GMM K-robustness reported; panel
        composition caveats stated; "structural stability" replaced by
        "structure-conditioned constraint"; clinical transfer framed as
        cross-domain transfer without clinical-label retraining; BioGate
        demoted to supporting experiment; VESPA reclassified as a
        sequence-embedding predictor; 217/215/186 clarified.
"""
import os
import re
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT  # noqa: E402

OUT = ROOT / "manuscript" / "submission"
OUT.mkdir(parents=True, exist_ok=True)

TITLE = ("Protein AI Model Disagreement Tracks AlphaFold Structural Confidence")

ABSTRACT = """Protein AI models -- evolutionary (MSA), single-sequence language, and structure-conditioned predictors -- are usually treated as competing estimators of mutation effect. Here we ask what it means, biologically, when they disagree. Using 696,311 single missense substitutions across 217 deep-mutational scanning assays (ProteinGym v1.3) scored by 40 architecturally diverse zero-shot models, we show that disagreement contains reproducible biological structure: it tracks AlphaFold structural confidence.

Within proteins, evolution-vs-structure disagreement decreases monotonically with pLDDT (median residue-level Spearman -0.04, 95% CI -0.06 to -0.03; 67.5% of 154 proteins negative; p < 0.001), robust to multivariable adjustment, panel resampling and four definitions of total disagreement; low-confidence residues are ~1.5-1.8x more likely to be maximally disagreeable. Two external controls argue against an input-quality artifact: ESM-IF1 rescoring on experimental structures (64 assays) runs opposite to that explanation, and SSEmb, a distinct joint sequence-structure architecture, reproduces the trend. Curated disorder marks a related but distinct regime: total, not pair-specific, disagreement is elevated there and single-sequence models lose relative advantage. Disagreement regimes, discrete summaries of a continuous geometry, generalize to held-out proteins (frozen-centroid protein-level cross-validation: median Spearman 0.94; 100% of 185 positive) and replicate across independent assays (median 0.94; 37 pairs overall, 28 on held-out proteins).

Disagreement is only weakly associated with prediction error in DMS (assay-level Spearman 0.032) and not at all in clinical variant classification. A biologically gated ensemble does not beat uniform averaging on DMS (paired delta +0.000; 95% CI -0.008 to +0.008), yet a DMS-trained two-family gate transfers to clinical data without retraining on clinical labels (AUROC +0.021; 95% CI +0.011 to +0.026 on a strictly non-overlapping set), though the best single expert remains superior. Disagreement among protein AI models is an AI-generated map of where evolutionarily informed and structure-conditioned predictions diverge across protein space; its primary value is explanatory, not predictive."""

KEYWORDS = ["variant effect prediction", "protein language models",
            "deep mutational scanning", "model disagreement",
            "AlphaFold structural confidence"]

BODY = r"""## 1. Introduction

Missense variation is central to protein function and human disease, and predicting its effects is a foundational problem in computational biology. Over the past five years the field has moved from evolutionary statistics (GEMME [5], EVE [3], DeepSequence [4], EVmutation [33]) to single-sequence protein language models (ESM-1v [6], ESM-1b [29], ESM-2 [7], ProtTrans [30], ProGen [32], ProGen2 [31], VESPA [27]) and structure-conditioned predictors (ESM-IF1 [8], ProteinMPNN [9], SaProt [11], MIF/MIF-ST [26]), culminating in large foundation models (ESM3 [12], ProGen3 [13], xTrimoPGLM [14], AlphaMissense [34]) whose zero-shot mutation scores are evaluated against deep mutational scanning (DMS) experiments at scale (ProteinGym [1]; DMS methodology [28, 41]).

These models are usually benchmarked as competitors: one average performance number decides the winner. But average performance hides structure. Models are trained on different modalities -- evolutionary conservation (MSA families), sequence statistics alone, or three-dimensional structure -- and there is no a priori reason they should encode identical biological constraints. A substitution that violates conservation while preserving local folding compatibility, or one in a flexible region for which the available structural representation is uncertain, may split the panel cleanly. The information the models share is a consensus estimate; the information on which they diverge may be information about biology.

We therefore treat disagreement among architecturally diverse protein AI models as a scientific observable: a residue-level measurement of where biological constraints are encoded inconsistently across model families. The question is not "who wins?" but "what biology does disagreement reveal?"

We test this with a confirmatory analysis, whose plan and GO/NO-GO thresholds were fixed in a versioned project configuration before the final analyses, on 696,311 single substitutions across 217 DMS assays scored by a 40-model panel: (i) whether models encode overlapping but distinct constraints; (ii) whether disagreement tracks structural confidence (AlphaFold pLDDT), curated disorder, and residue mechanics; (iii) whether disagreement regimes correspond to distinct biological contexts and reproduce across independent experiments; (iv) whether model families exhibit context-dependent expertise; (v) whether disagreement marks prediction difficulty in DMS and in clinical variant classification; and (vi) whether a biologically gated ensemble transfers across domains. Section 3.0 positions the work against fixed-weight predictor integration, uncertainty estimation, and the ProteinGym benchmark itself.

## 2. Results

### 2.1 Protein AI models encode overlapping but distinct mutational constraints

Across the 215 assays with sufficient variant coverage for pairwise correlation (2 of 217 assays contained fewer than 50 qualifying mutations and were excluded from this specific analysis only; all 217 assays enter the primary disagreement analyses), the median pairwise Spearman correlation between model predictions within the same family was 0.698 versus 0.587 between families. PCA of the variant-by-model matrix showed a dominant shared component (PC1 = 58.5%) concentrated in evolutionary models, plus a secondary axis (PC2 = 6.0%) separating structure-conditioned from sequence-based models (figure 1(c)); pairwise correlations (figure 1(d)) and hierarchical clustering (figure 1(e)) recapitulated the training modality families.

Score orientation followed the released ProteinGym conventions and the official model metadata (config.json); post hoc correlations with DMS were used only as a sanity check and never to determine orientation (Methods 4.2).

### 2.2 AlphaFold structural confidence tracks disagreement (primary result)

The primary statistic is protein-wise residue-level correlation. Within each protein, mutations at the same residue were aggregated to the residue median of disagreement, and Spearman(pLDDT, disagreement) was computed directly over residues (154 proteins with >= 50 mapped residues). The median within-protein Spearman between pLDDT and evolution-vs-structure disagreement was **-0.04** (protein-bootstrap 95% CI -0.06 to -0.03; 67.5% of proteins negative; one-sample Wilcoxon p < 0.001), and -0.08 for total disagreement (71.4% negative; the distribution of within-protein residue-level correlations is shown in figure 2(a)). The effect is small in absolute terms but highly consistent across proteins.

The effect is also interpretable in intuitive terms (figure 2(b); Table S1): within proteins, residues in the lowest pLDDT decile have higher evolution-vs-structure disagreement than residues in the highest decile (median difference +0.057, bootstrap 95% CI +0.013 to +0.084; Cohen's d = 0.15, CI 0.04-0.23), and low-pLDDT-decile residues are substantially more likely to fall in the protein-internal top disagreement decile (median protein-level odds ratio 1.84, 95% CI 1.62-2.04; pooled OR 1.49, 95% CI 1.35-1.65). That is, residues in the least confident structural regions are roughly 1.5-1.8x more likely to be among the most disagreeable in the protein.

The trend is robust to the total-disagreement definition (Methods 4.3): z-based SD (-0.078), percentile SD (-0.115), MAD (-0.076) and mean pairwise distance (-0.089) all give negative residue-level associations with confidence intervals excluding zero and 71-77% of proteins negative (figure 2(c); Table S1). A fifth candidate metric (variance of per-variant ranks) is mathematically redundant with the percentile SD and is not used (Methods 4.3).

The effect is specific to family pairs: evolution-vs-structure disagreement shows the negative trend, whereas sequence-vs-structure disagreement does not (median residue-level rho 0.01, CI overlapping zero; Table S2). Disagreement therefore accumulates in low-confidence regions specifically as a conflict between evolutionary-model and structure-conditioned model predictions.

Confound controls (shared AlphaFold input). The structure-conditioned models and the pLDDT metric derive from the same AlphaFold structures, so low pLDDT might simply degrade structure-model inputs. We report three partial controls (Table S2). (i) Asymmetry: if disagreement were a pure input-quality artifact, sequence-vs-structure disagreement should rise at low pLDDT exactly like evolution-vs-structure disagreement, since the same structure models are involved; it does not (median residue-level rho +0.01 vs -0.04; asymmetry -0.03, 65% of proteins in the evo-struct direction). (ii) Structure-model accuracy: structure-conditioned model scores do lose DMS-correlation at low pLDDT (e.g., ESM-IF1 median per-assay rho 0.04 at pLDDT<50 vs 0.47 at >=90), confirming that their scores carry less information there -- the input-quality pathway is real and is a documented limitation, not fully resolvable without experimental structures. (iii) Coverage: structure-model score availability is complete at all pLDDT levels (100%), ruling out missingness-driven artifacts. Beyond these internal controls we use two external datasets. (iv) Experimental-structure gain: for 64 assays that have ESM-IF1 scores computed on BOTH experimental PDB structures and AlphaFold2 structures [56], the assay-level gain (experimental minus AlphaFold2) is negative on average (-0.057) and significantly MORE negative for low-confidence assays (median gain -0.081 vs -0.023; difference -0.058, bootstrap 95% CI -0.099 to -0.013, Mann-Whitney p = 0.011; Spearman(mean pLDDT, gain) = +0.38, p = 0.002; figure S10): replacing AlphaFold inputs with experimental structures does not rescue low-confidence assays, contrary to what a pure input-quality explanation predicts. (v) Distinct-architecture replication: SSEmb -- which jointly embeds MSA information and structure in one architecture (an MSA-Transformer language model combined with a structure-conditioned graph network) and was trained under a different pipeline [57] -- shows the same residue-level trend for evolution-vs-SSEmb disagreement (median rho -0.067, 65.6% of 151 proteins negative; figure S11), matching the inverse-folding result (-0.042). A variant-level experimental-structure re-scoring of all models remains the decisive future control and is stated as a limitation.

Panel-composition checks. (a) A strict one-representative-per-architecture panel (ESM-1v, GEMME, EVE, MSA-Transformer, ESM-IF1, ProteinMPNN, SaProt, Tranception) reproduces the residue-level trend (median rho -0.06; 68% of proteins negative). (b) 1,000 resamplings drawing three models from the evolutionary and sequence families reproduce the trend at residue level (median of per-iteration medians -0.05, 95% range [-0.08, -0.02]; 100% of iterations negative; figure 2(c), right). Because only three structure-family models were available in the released scores, the structure panel was fixed in these resamplings; composition robustness is therefore established for the evolutionary and sequence families, and the structure family is checked by the architecture panel in (a).

Independent disorder annotation. Curated UniProt disordered regions (transferred to DMS coordinates by sequence alignment; 92 proteins covered) confirm that total disagreement is higher in disordered than in ordered residues (median within-protein difference +0.040, 76.5% of proteins positive; figure 2(d)). However, the evolution-vs-structure conflict does not concentrate in curated IDR (median IDR-ordered difference -0.017, 38% positive). The primary signal is therefore located in the broader low-confidence continuum -- regions consistent with flexibility, conformational heterogeneity, or weak structural constraint -- rather than in canonical disorder alone. pLDDT<50 and UniProt disorder overlap partially (P(low-confidence | disorder) = 0.56), consistent with the literature linking AlphaFold confidence to experimental disorder [16, 39]; we retain the "low structural confidence" terminology for the metric itself.

### 2.3 Residue mechanics: contact density and secondary structure

Using pydssp secondary-structure assignment and CA contact density from the same AlphaFold structures, we tested whether family relative advantage varies with residue mechanics. In the percentile space used throughout for accuracy comparisons, no mechanical context showed a consistent family advantage (all median advantages within +/-0.013; helix/sheet/core/surface/loop), and contact density showed no monotonic association with evolution-vs-structure disagreement (median within-protein rho 0.01). These analyses are reported as a null result: residue mechanics, as captured by secondary structure and contact density, does not shift relative family expertise beyond the confidence axis (figure S9). Disagreement likewise does not concentrate at UniProt domain boundaries (median disagreement is flat across distance-to-boundary bins: 0.63, 0.59, 0.58 and 0.63 for 0-5, 6-10, 11-20 and >20 residues; figure S3).

### 2.4 Disagreement regimes: discrete summaries of a continuous geometry

Gaussian mixture modeling of the three mechanistic family predictions (BIC-optimal [53] K = 6; figure 3(a)-(c)) assigns every variant to a regime (Table 1). We treat regimes as discrete summaries of an underlying continuous disagreement geometry, not as evidence of discrete natural categories.

Consensus regimes align with the measured DMS phenotype: consensus-damaging variants have median experimental fitness rank 0.19 (damaged), consensus-tolerant 0.74 (tolerant). In the structure-dissenting regime, the measured DMS phenotype aligns more closely with sequence and evolutionary predictions than with structure-conditioned predictions (median experimental Y = 0.65). We state this as an alignment-of-evidence statement and do not claim that structure-conditioned models are "wrong": DMS phenotypes (activity, binding, abundance, fitness) are not direct measures of folding compatibility. (Regime names are descriptive labels assigned post hoc from the actual centroids; machine-readable labels are regime_0-regime_5.)

Robustness of the regime structure. Consensus regimes and their experimental alignment are stable across the number of clusters (K = 4-7; Table S3). The structure-dissenting pattern emerges at K >= 6; at K = 4-5 it is absorbed into broader regimes, so this specific pattern is a K-dependent characterization rather than a universal feature.

Experimental replication of regime-phenotype relationships. Because regime labels derive from model predictions, agreement of assignments across independent assays of the same protein partly reflects deterministic model behaviour on the same sequence (a form of circularity). The scientifically meaningful replication is whether regime-phenotype relationships hold across assays: for 37 assay pairs sharing >= 3 regimes with >= 20 variants each, the median Spearman between regime-level experimental Y across assays was **0.94**; consensus-damaging regimes were damaging (Y > 0.5) in both assays of every pair (100%), consensus-tolerant regimes tolerated in every pair (100%), and the structure-dissenting regime was damaging in both assays in 92% of pairs where it was present (figure 4(a), inset). Separately, 34 of these 37 pairs additionally had >= 20 shared residue positions, allowing position-level assignment agreement (0.763 vs 0.215 chance; 100% of pairs above chance). The two analyses use different inclusion criteria and different pair counts; it is the phenotype replication that establishes experimental reproducibility. Held-out protein generalization. With a stricter test - protein-level cross-validation with frozen GMM centroids (models fitted on four training folds, assignments made to completely unseen held-out proteins; 185 proteins) - the regime-level phenotype profile of held-out proteins reproduced the training profile with median Spearman 0.94 (interquartile 0.89-1.00; 100% of proteins positive), and the structure-dissenting regime was damaging on held-out data (median Y 0.61; damaging in 72% of proteins with sufficient held-out support; Table S8). Regimes therefore generalize beyond the proteins used to define them. Adjusted Rand index between K = 6 and K = 5 is 0.58 (K = 4: 0.35, K = 7: 0.46). Illustrative per-protein disagreement maps are shown in figure 4(b)-(d).

### 2.5 Relative family expertise: an IDR-specific effect, no confidence gradient

We asked whether family advantage (error of other families minus own error, all in percentile space) depends on biological context (figure 5(a)-(c); Table 2). Two results follow. First, relative expertise does NOT vary along the structural-confidence axis: family advantages are small and flat across pLDDT bins (single-sequence +0.008 to +0.022; evolution +0.002 to +0.005; structure -0.002 to -0.007). Second, curated disorder shows a specific effect: single-sequence models lose relative predictive advantage in intrinsically disordered regions (-0.043) versus ordered residues (+0.038; 88% of proteins positive in ordered), while structure-conditioned models gain relative advantage in IDR (+0.060 vs +0.005). Evolutionary models show no consistent effect. We conclude that the readable context-dependence of family expertise is confined to disorder annotation, not to the confidence continuum, and is modest in size.

### 2.6 Disagreement is only weakly predictive (secondary, transparent)

The uncertainty read-out is weak and we keep it secondary: median assay-level Spearman(disagreement, ensemble error) = 0.032 (UniProt-cluster bootstrap 95% CI 0.023-0.047; 70.7% of assays positive; distribution in figure S1), with error defined in percentile space (|ensemble prediction - Y|, both on the 0-1 scale). The effect is monotonic across disagreement deciles (top-decile error 0.214 vs 0.201; figure S2), unchanged within pLDDT strata (0.027-0.040; not a confound), stable to leave-one-model-out (0.027-0.039; figure S6), and consistent across the D_MAD definition (0.026). The association is strongest in stability and binding assays (median rho 0.076 and 0.057; figure S8) and effectively absent in activity and organismal-fitness assays (0.003 and -0.004). The primary value of disagreement is explanatory rather than predictive.

### 2.7 Clinical data: cross-domain transfer of the gating function, not better prediction

The clinical benchmark contains 2,525 genes / ~63K variants; of these, 1,530 proteins had complete two-family gate inputs and a computable per-protein AUROC, and 788 remained after the strict non-overlap filter (below). Applying the DMS-trained BioGate (a separate two-expert gate -- sequence + evolution, the only families available in the clinical score files -- trained exclusively on DMS, frozen, and applied once to the clinical benchmark without clinical-label retraining) gave median protein-level AUROC 0.944 on the full clinical set (n = 1,530), versus 0.877 for the uniform ensemble and 0.963 for the best single expert (TranceptEVE_L) figure 6(a). BioGate is therefore a supporting experiment: it shows that the arbitration function learned from DMS transfers across domains, not that gating outperforms the best base predictor. All base models are pre-trained on large sequence databases and have seen sequence space broadly; "transfer" here strictly means that the DMS-trained gating function was applied without clinical-label retraining, and the two-expert gate is a distinct model from the three-expert DMS BioGate (Methods 4.4).

Strict non-overlap control. Clinical proteins sharing >= 70% sequence identity with any of the 186 DMS proteins used to train the gate (35 proteins; same-protein overlaps are necessarily captured by this threshold) were excluded, leaving 60,214 variants across 2,490 proteins. On this strict set the uniform ensemble AUROC was 0.905 (all 2,490 proteins), and restricting to the 788 proteins with complete two-family gate inputs, AUROC was 0.867 (uniform), 0.899 (BioGate) and 0.974 (best single, reselected within the strict set) (Table 3). The transfer is insensitive to the identity threshold (30% / 50% / 70% all exclude the same 35 proteins and give uniform AUROC 0.905; Table S4). The DMS-trained gating function therefore retains a modest relative advantage over uniform averaging on a strictly non-overlapping set (paired median delta AUROC +0.021, 95% CI +0.011 to +0.026, 59% of proteins positive; Table S7); it does not, however, exceed the best base expert (paired delta -0.061, 93% of proteins where best > BioGate). The gate was a separately trained two-family model (sequence + evolution) learned exclusively on DMS, frozen, and applied once - a gating STRATEGY transferred, not the three-expert DMS gate itself. The gate-complete subsets (1,530 full; 788 strict) are defined by input completeness and are not representative of the full 2,525-gene benchmark.

Disagreement does not mark clinical difficulty. In contrast to DMS, misclassified clinical variants have lower disagreement than correctly classified ones (median 0.662 vs 0.693; 41.4% of proteins with higher D in misclassified; (figure 6(b)), and pathogenic vs benign variants show no disagreement difference (0.693 vs 0.690). Disagreement-as-uncertainty is therefore domain-specific: a property of well-posed fitness landscapes, not of heterogeneous disease annotation.

### 2.8 BioGate on DMS: a transparent negative

A softmax-gated ensemble ("BioGate"; one 32-unit hidden layer) weighting the three percentile-space family centroids from context features (pLDDT, position, length, MSA depth, functional annotations), evaluated by 5-fold UniProt-grouped CV, did not beat uniform averaging: median protein-level rho 0.551 vs 0.542, while the median within-protein paired difference was 0.000 (bootstrap 95% CI -0.008 to +0.008); the two summaries differ because the median of a difference is not the difference of medians (Supplementary figure S4) (Table 4). XGBoost ablations show context does add signal (scores-only 0.509, with context 0.551 under the per-fold-median aggregation of Table S5; 0.573 under the pooled aggregation of Table 4), but the constrained convex-mixture gate cannot exploit it. BioGate is a constraint result: with a per-variant uncertainty signal this weak, gating cannot improve point prediction (details in Supplementary figures S4-S6).

### 2.9 The uncertainty read-out is scale-dependent (methodological note)

Disagreement-based uncertainty is not invariant to monotonic score transformation: holding the prediction error fixed (|ensemble prediction - Y| in percentile space), the disagreement-error association is ~4.5x stronger when disagreement is measured on the percentile (u) scale (median assay-level rho 0.145) than on the inverse-normal (z) scale (0.032), because the z-transform stretches the tails, where models typically agree on deleteriousness (figure S7). We fixed the z scale as primary in the versioned analysis configuration; the substantive conclusion of this paper -- that disagreement is a weak uncertainty signal in DMS and none in clinical annotation -- holds on either scale, but the effect magnitude is scale-dependent. Any future use of disagreement as confidence should therefore pre-specify the scale, and our primary value claim for disagreement is explanatory rather than predictive.

## 3. Discussion

### 3.0 Related work and positioning

(i) Meta-predictors and integrative scores (MetaSVM/MetaLR [17], ClinPred [18], BayesDel (distributed within the PERCH framework) [19], REVEL [35], CADD [36]) combine evolutionary, conservation and structural features through fitted or heuristic weights applied uniformly across variants; they cannot answer where and why predictor trust should differ across the mutational landscape. Our contribution is the quantitative description of the biology of disagreement itself, and the demonstration that relative family expertise is context-dependent only in a limited sense: single-sequence models lose relative advantage in curated disorder, while no expertise gradient exists along the structural-confidence axis (Section 2.5).

(ii) Model disagreement as uncertainty in ML (deep ensembles [20], Bayesian approximations, active learning [21]). We confirm a weak uncertainty signal in DMS (Section 2.6) but show it is domain-specific (Section 2.7) and, more importantly, that disagreement carries a reproducible biological structure (Sections 2.2-2.5) that generic uncertainty theory neither predicts nor explains.

(iii) The ProteinGym benchmark [1] reports average per-model performance. Our analysis is not a benchmark: the 40-model panel is an instrument, every reported quantity is a residue-level property of the landscape, and we surface a methodological caveat absent from benchmark practice -- the disagreement scale (u vs z) changes the disagreement-error association ~4.5x, so future disagreement analyses must pre-specify their scale.

### 3.1 Disagreement as a probe of prediction-space divergence

The central observation is that evolution-vs-structure disagreement accumulates where AlphaFold confidence is low. Mechanistically, this is what one would expect if the constraints encoded by confidently folded structures and the constraints encoded by evolutionary sequence statistics become partially decoupled in flexible, poorly constrained and conformationally heterogeneous regions [16, 39, 40]. We deliberately state this as a plausible biological interpretation rather than a demonstrated mechanism: pLDDT is a model confidence, not an experimental stability measurement, and structure-conditioned model scores reflect residue/sequence compatibility with a given backbone, not free energy of folding.

Two caveats bound the interpretation. First, the structure-conditioned models and the confidence metric share AlphaFold structures as input. The external controls of Section 2.2 argue against a simple input-quality explanation - experimental structures do not rescue (they hurt) low-confidence assays for ESM-IF1, and a structure-conditioned model with a distinct joint sequence-structure architecture reproduces the trend - but these controls are assay-level or single-model, and a full variant-level experimental-structure re-scoring [54, 56] remains the decisive test and the primary limitation. Throughout, "structural constraints" should be read as "constraints encoded by the evaluated structure-conditioned models".

Second, the interpretation is supported by three complementary lines of evidence rather than by any single statistic: (i) the residue-level disagreement-pLDDT gradient, small but consistent, robust across four disagreement definitions and to multivariable adjustment; (ii) the IDR-specific loss of relative predictive advantage of single-sequence models, replicated against curated disorder; and (iii) the reproducible regime-phenotype structure, whose consensus components replicate experimentally across independent assays. Each line is individually modest; together they indicate that the divergence between evolutionarily informed and structure-conditioned model predictions tracks the AlphaFold confidence landscape of proteins.

### 3.2 What the clinical asymmetry teaches

The failure of disagreement to mark difficulty in clinical data is as informative as the weak but reproducible DMS association. DMS assays measure one well-defined phenotype under controlled conditions [41]; clinical labels aggregate heterogeneous phenotypes, cell types and genetic backgrounds curated across submitters [23]. The boundary condition we document -- a weak disagreement-error association is detectable in fitness landscapes but does not generalize to clinical annotation -- prevents mechanical misapplication of disagreement-based confidence in clinical genomics.

### 3.3 Limits of disagreement-based arbitration

The reproducible structure of disagreement is statistically consistent but per-variant weak; the constrained convex-mixture gate cannot convert it into better DMS point predictions (Section 2.8). Clinical transfer is modest: BioGate achieved 0.944 versus 0.877 for uniform averaging on the full gate-complete set, and 0.899 versus 0.867 on the strict non-overlapping gate-complete subset; the best single expert remained superior (0.963 full; 0.974 strict). These results position the gated ensemble as a supporting, exploratory component rather than a predictive contribution.

## 4. Methods

### 4.1 Data
ProteinGym v1.3 [1] (DOI 10.5281/zenodo.15293562): 217 DMS substitution assays (~2.7M variants), precomputed zero-shot scores for 95 released models, the clinical substitution benchmark (2,525 proteins, ~63K variants) with its zero-shot scores, and AlphaFold2 structures [15, 38] for assay proteins (197 PDB files covering all 217 assays). Primary set: single amino-acid substitutions (regex `^[A-Z][0-9]+[A-Z]$`) whose wild-type residue matches the assay target sequence (0 mismatches): 696,311 variants across 186 proteins (multiple assays can share one protein). Two of the 217 assays contained fewer than 50 qualifying mutations and were excluded only from the pairwise model-correlation analysis (Section 2.1); all 217 assays enter the primary disagreement analyses.

### 4.2 Model panels and score orientation
Official model metadata (family, directionality) from the ProteinGym config.json. Released score files carry directionality already applied (higher = higher fitness); the direction table is recorded with evidence labels. Post hoc DMS correlations (all 95 models positive; median 0.427) were used only as a sanity check and never to determine orientation.

Two panels serve different purposes. (a) The 40-model core panel (one representative per architecture, including GEMME [5], EVE [3], DeepSequence [4], MSA-Transformer [10], Tranception and TranceptEVE [2], ESM-1v [6], ESM-2 [7], ESM-IF1 [8], ProteinMPNN [9], SaProt [11], MIF/MIF-ST [26]) quantifies total disagreement, PCA and pairwise model analysis; it spans four official families (evolution 12, single-seq 14, structure 3, hybrid 11). (b) The three-family mechanistic panel (evolution, single-sequence, structure) is used for all modality-attributable analyses (family centroids, regimes, advantage); hybrid models were included in (a) but excluded from (b) because their scores cannot be uniquely attributed to sequence-, evolutionary- or structure-derived constraints. All family-level conclusions were re-verified with a strict 8-model architecture panel and 1,000 resamplings of the evolutionary and sequence families (structure family fixed by availability; Section 2.2).

Score coverage and missingness. All 696,311 retained variants had complete scores from all 40 core models and all three mechanistic families after the panel coverage filter; no imputation was required. Table S6 lists, for every core model, its architecture, input modalities (MSA/sequence/structure), whether its structure input is AlphaFold-derived, and mechanistic-panel membership. Held-out regime protocol: 5-fold GroupKFold by UniProt; GMMs fitted on training proteins only; held-out variants assigned to frozen centroids; phenotype-profile agreement scored per held-out protein (Table S8).

### 4.3 Normalization, disagreement scale and definitions
Per assay and per model: percentile rank -> u (deleterosity orientation; u = 1 - rank since higher released score = fitness) -> z = Phi^-1(clip(u, 0.001, 0.999)). Y = 1 - rank(DMS_score) within assay (higher Y = more damaged; Y on the 0-1 scale).

Two family-centroid spaces are used for distinct purposes:
$$Z_{{\mathrm{{fam}},i}}=\frac{{1}}{{|F|}}\sum_{{m\in F}}z_{{im}}\quad\text{{(z-space; regime clustering and Table 1)}}$$

$$U_{{\mathrm{{fam}},i}}=\frac{{1}}{{|F|}}\sum_{{m\in F}}u_{{im}}\quad\text{{(percentile space; accuracy, advantage, BioGate)}}$$

All accuracy comparisons and predictions (family ensembles, errors, gating)
are therefore computed in the percentile space, so that predictions and
outcomes share the same 0-1 scale.

Primary quantities (explicit definitions):
$$D_\mathrm{std},i=\sqrt{{\frac{{1}}{{M}}\sum_m (z_{{im}}-\bar z_i)^2}}$$

$$D_\mathrm{{MAD}},i=\mathrm{{median}}_m|z_{{im}}-\mathrm{{median}}_m z_i|$$

$$D_u,i=\mathrm{{SD}}_m(u_{{im}})$$

$$D_\mathrm{{pair}},i=\frac{{2}}{{M(M-1)}}\sum_{{m<k}}|z_{{im}}-z_{{ik}}|$$

$$D_{{\mathrm{{evo}}\text{{-}}\mathrm{{struct}},i}}=|U_{{\mathrm{{evolution}},i}}-U_{{\mathrm{{structure}},i}}|$$

$$S_{{\mathrm{{ens}},i}}=\frac{{1}}{{3}}\sum_f U_{{f,i}},\qquad \mathrm{{error}}_i=|S_{{\mathrm{{ens}},i}}-Y_i|$$

Four disagreement definitions are used (D_std, D_u, D_MAD, D_pair); a fifth
candidate (variance of per-variant ranks) is mathematically redundant with
D_u (var(u) = SD(u)^2) and is not used. The pLDDT organization result holds
under all four (Section 2.2). Scale pre-specification for the uncertainty
analysis: holding error fixed in percentile space, u-basis disagreement
correlates with error ~4.5x more strongly than z-basis (0.145 vs 0.032)
because z stretches tails where models agree; the z scale is primary and u
reported as sensitivity (Section 2.9).

### 4.4 Statistics
Primary pLDDT analysis: residue-level Spearman per protein, aggregated by median; protein-level bootstrap (2,000) confidence intervals; one-sample Wilcoxon test; four disagreement definitions; protein-fixed-effects regression (within-protein demeaning, protein-clustered standard errors) adjusting for disorder, contact density, secondary structure, position and MSA depth (Table S2). Cross-assay summaries aggregate per-protein effects; UniProt-level cluster bootstrap (10,000) [51]. GO/NO-GO thresholds were fixed in the versioned project configuration before the confirmatory analyses (config.yaml; not a public registry). Within-protein permutation enrichment with Benjamini-Hochberg FDR control [52]. Regime clustering: Gaussian mixtures K = 2-8 by BIC [53]; K-robustness reported (K = 4-7, adjusted Rand index); experimental replication evaluated as regime-level Y correlation across independent assays of the same protein. BioGate: torch [49] CPU MLP gating, 5-fold GroupKFold by UniProt; baselines best single expert (train-selected), uniform ensemble, Ridge stacking, XGBoost [25] (scores-only / context-only / scores+context ablations); paired protein-level bootstrap (2,000) [51]. Clinical transfer: a separate two-expert gate (sequence + evolution, the only families in the clinical score files) was trained exclusively on DMS, frozen, and applied once to the clinical benchmark; no clinical label was used in any model selection; per-protein AUROC; strict non-overlap set by >= 70% sequence identity (8-mer filter + alignment; 30% / 50% / 70% sensitivity in Table S4) excluding 35 clinical proteins. Software: Python ecosystem (numpy [44], scipy [45], pandas [47], scikit-learn [43], statsmodels [48], matplotlib [46]), Biopython [24], PyTorch [49], pydssp [50]; transformer architectures follow Vaswani et al. [42].

### 4.5 Structural, mechanical and disorder annotation
pLDDT from ProteinGym AlphaFold2 PDBs (B-factor, best chain; 99.0% of variants mapped; two proteins with < 50% alignment flagged and excluded from pLDDT analyses). AlphaFold2 structures and confidence scores are those of Jumper et al. [15]. pydssp secondary structure (H/E/C) and CA contact density (8/10 A) from the same PDBs; burial tertiles within protein. UniProt features (active site, binding site, domain, transmembrane, modified residue, disulfide, region, motif) and curated disordered regions fetched via the UniProt REST API [22] and transferred to DMS coordinates by sequence alignment (Needleman-Wunsch global alignment, Biopython PairwiseAligner; mean per-protein coverage 0.999 and identity 0.993; 179 proteins mapped; 65,420 positions agree with the original heuristic matcher at 99.92%, and IDR mapping quality is reported per protein in mapping_quality.csv; 44,708 annotations transferred; 4,569 disordered positions across 92 proteins). The clinical benchmark labels derive from ClinVar curation [23]. Case-study proteins (TP53, BRCA1, PTEN; figure 4(b)-(d)) were selected a priori as clinically prominent proteins with dense DMS coverage (>1,000 mutations each), high-quality structures, and broad pLDDT coverage; they are illustrative rather than exhaustive.

### 4.6 Reproducibility
Seed 2026; every figure has a machine-readable source table; all scripts log inputs/outputs/counts/runtime; `python src/16_make_figures.py` recreates all main figures from results/.

## 5. Conclusions

Protein AI disagreement is not merely model noise: it contains reproducible biological structure. It tracks the protein's structural-confidence landscape: evolution-vs-structure disagreement accumulates where AlphaFold confidence is low, total disagreement is elevated in intrinsically disordered regions, and the disagreement geometry decomposes into regimes whose phenotypes replicate across independent assays and generalize to held-out proteins. Relative model-family expertise shows limited context dependence: single-sequence models lose relative predictive advantage in curated disordered regions, whereas no consistent expertise gradient is observed across AlphaFold confidence -- disagreement and relative accuracy are distinct quantities. Disagreement is a weak uncertainty signal in DMS and none in clinical annotation -- a boundary condition that should discipline future use of disagreement-based confidence. Its primary value is explanatory: disagreement among protein AI models is an AI-generated map of where evolutionarily informed and structure-conditioned model predictions diverge across protein space.

## Data availability
All data are public (ProteinGym v1.3, DOI 10.5281/zenodo.15293562; UniProt). All analysis code (scripts 01-37), machine-readable source tables for every figure, and figure files are openly available at https://github.com/hmjpan/protein-ai-disagreement ; a versioned archive with a citable DOI will be deposited upon acceptance.

## References
[1] Notin P, Kollasch A, Ritter D, van Niekerk L, Paul S, Spinner H, Rollins N, Shaw A, Orenbuch R, Weitzman R, Frazer J, Dias M, Franceschi D, Gal Y and Marks D S 2023 ProteinGym: large-scale benchmarks for protein fitness prediction and design Proc. NeurIPS 2023 Datasets and Benchmarks Track
[2] Notin P, Dias M, Frazer J, Marchena-Hurtado J, Gomez A N, Marks D S and Gal Y 2022 Tranception: protein fitness prediction with autoregressive transformers and inference-time retrieval Proc. 39th Int. Conf. on Machine Learning (ICML)
[3] Frazer J, Notin P, Dias M, Gomez A, Min J K, Brock K, Gal Y and Marks D S 2021 Disease variant prediction with deep generative models of evolutionary data Nature 599 91-5
[4] Riesselman A J, Ingraham J B and Marks D S 2018 Deep generative models of genetic variation capture the effects of mutations Nat. Methods 15 816-22
[5] Laine E, Karami Y and Carbone A 2019 GEMME: a simple and fast global epistatic model predicting mutational effects Mol. Biol. Evol. 36 2604-19
[6] Meier J, Rao R, Verkuil R, Liu J, Sercu T and Rives A 2021 Language models enable zero-shot prediction of the effects of mutations on protein function Adv. Neural Inf. Process. Syst. 34
[7] Lin Z et al 2023 Evolutionary-scale prediction of atomic-level protein structure with a language model Science 379 1123-30
[8] Hsu C, Verkuil R, Liu J, Lin Z, Hie B, Sercu T, Lerer A and Rives A 2022 Learning inverse folding from millions of predicted structures Proc. 39th Int. Conf. on Machine Learning (ICML)
[9] Dauparas J et al 2022 Robust deep learning-based protein sequence design using ProteinMPNN Science 378 49-56
[10] Rao R, Liu J, Verkuil R, Meier J, Canny J, Abbeel P, Sercu T and Rives A 2021 MSA Transformer Proc. 38th Int. Conf. on Machine Learning (ICML)
[11] Su J, Han C, Zhou Y, Shan J, Zhou X and Yuan F 2024 SaProt: protein language modeling with structure-aware vocabulary Int. Conf. on Learning Representations (ICLR)
[12] Hayes T et al 2025 Simulating 500 million years of evolution with a language model Science 387 850-8
[13] Bhatnagar A, Jain S, Beazer J et al 2025 Scaling unlocks broader generation and deeper functional understanding of proteins bioRxiv 2025.04.15.649055
[14] Chen B et al 2023 xTrimoPGLM: unified 100B-scale pre-trained transformer for deciphering the language of proteins bioRxiv 2023.07.05.547496
[15] Jumper J et al 2021 Highly accurate protein structure prediction with AlphaFold Nature 596 583-9
[16] Piovesan D, Monzon A M and Tosatto S C E 2022 Intrinsic protein disorder and conditional folding in AlphaFoldDB Protein Sci. 31 e4466
[17] Dong C, Wei P, Jian X, Gibbs R, Boerwinkle E, Wang K and Liu X 2015 Comparison and integration of deleteriousness prediction methods for nonsynonymous SNVs in whole exome sequencing studies PLoS ONE 10 e0134848
[18] Alirezaie N, Kernohan K D, Hartley T, Majewski J and Richer T D 2018 ClinPred: prediction tool to identify disease-relevant nonsynonymous single-nucleotide variants Am. J. Hum. Genet. 103 474-83
[19] Feng B-J 2017 PERCH: a unified framework for disease gene prioritization Hum. Mutat. 38 243-51
[20] Lakshminarayanan B, Pritzel A and Blundell C 2017 Simple and scalable predictive uncertainty estimation using deep ensembles Adv. Neural Inf. Process. Syst. 30
[21] Beluch W H, Genewein T, Nurnberger A and Kohler J M 2018 The power of ensembles for active learning in image classification Proc. IEEE/CVF Conf. on Computer Vision and Pattern Recognition (CVPR)
[22] UniProt Consortium 2023 UniProt: the universal protein knowledgebase in 2023 Nucleic Acids Res. 51 D523-31
[23] Landrum M J et al 2018 ClinVar: improving access to variant interpretations and supporting evidence Nucleic Acids Res. 46 D1062-7
[24] Cock P J A et al 2009 Biopython: freely available Python tools for computational molecular biology and bioinformatics Bioinformatics 25 1422-3
[25] Chen T and Guestrin C 2016 XGBoost: a scalable tree boosting system Proc. 22nd ACM SIGKDD Int. Conf. on Knowledge Discovery and Data Mining
[26] Yang K K, Zanichelli N and Yeh H 2023 Masked inverse folding with sequence transfer for protein representation learning Protein Eng. Des. Sel. 36 gzad015
[27] Marquet C, Heinzinger M, Olenyi T, Dallago C, Erckert K, Bernhofer M, Nechaev D and Rost B 2022 Embeddings from protein language models predict conservation and variant effects Hum. Genet. 141 1629-47
[28] Tsuboyama K et al 2023 Mega-scale experimental analysis of protein folding stability in biology and design Nature 620 434-44
[29] Rives A et al 2021 Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences Proc. Natl Acad. Sci. USA 118 e2016239118
[30] Elnaggar A et al 2022 ProtTrans: toward understanding the language of life through self-supervised learning IEEE Trans. Pattern Anal. Mach. Intell. 44 7112-27
[31] Nijkamp E, Ruffolo J A, Weinstein E N, Naik N and Madani A 2023 ProGen2: exploring the boundaries of protein language models Cell Systems 14 968-978
[32] Madani A et al 2023 Large language models generate functional protein sequences across diverse families Nat. Biotechnol. 41 1099-106
[33] Hopf T A et al 2017 Mutation effects predicted from sequence co-variation Nat. Biotechnol. 35 128-35
[34] Cheng J et al 2023 Accurate proteome-wide missense variant effect prediction with AlphaMissense Science 381 eadg7492
[35] Ioannidis N M et al 2016 REVEL: an ensemble method for predicting the pathogenicity of rare missense variants Am. J. Hum. Genet. 99 877-85
[36] Rentzsch P, Witten D, Cooper G M, Shendure J and Kircher M 2019 CADD: predicting the deleteriousness of variants throughout the human genome Nucleic Acids Res. 47 D886-94
[38] Varadi M et al 2022 AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models Nucleic Acids Res. 50 D439-44
[39] Ruff K M and Pappu R V 2021 AlphaFold and implications for intrinsically disordered proteins J. Mol. Biol. 433 167208
[40] Necci M, Piovesan D, CAID Predictors, DisProt Curators and Tosatto S C E 2021 Critical assessment of protein intrinsic disorder prediction Nat. Methods 18 472-81
[41] Fowler D M and Fields S 2014 Deep mutational scanning: a new style of protein science Nat. Methods 11 801-7
[42] Vaswani A et al 2017 Attention is all you need Adv. Neural Inf. Process. Syst. 30
[43] Pedregosa F et al 2011 Scikit-learn: machine learning in Python J. Mach. Learn. Res. 12 2825-30
[44] Harris C R et al 2020 Array programming with NumPy Nature 585 357-62
[45] Virtanen P et al 2020 SciPy 1.0: fundamental algorithms for scientific computing in Python Nat. Methods 17 261-72
[46] Hunter J D 2007 Matplotlib: a 2D graphics environment Comput. Sci. Eng. 9 90-5
[47] McKinney W 2010 Data structures for statistical computing in Python Proc. 9th Python in Science Conf.
[48] Seabold S and Perktold J 2010 statsmodels: econometric and statistical modeling with python Proc. 9th Python in Science Conf.
[49] Paszke A et al 2019 PyTorch: an imperative style, high-performance deep learning library Adv. Neural Inf. Process. Syst. 32
[50] Thomas D pydssp: a fast secondary structure assigner for protein structures https://github.com/continuum-continuum/pydssp (accessed 2026)
[51] Efron B and Tibshirani R J 1993 An Introduction to the Bootstrap (New York: Chapman and Hall)
[52] Benjamini Y and Hochberg Y 1995 Controlling the false discovery rate: a practical and powerful approach to multiple testing J. R. Stat. Soc. B 57 289-300
[53] Schwarz G 1978 Estimating the dimension of a model Ann. Stat. 6 461-4
[54] Gitter laboratory benchmarking structure-based models on ProteinGym https://github.com/gitter-lab/benchmarking-structure-based-models (accessed 2026)
[56] Sharma A and Gitter A 2025 Exploring zero-shot structure-based protein fitness prediction Proc. ICLR Workshop on Generative and Experimental Perspectives for Biomolecular Design; experimental-structure scores at Zenodo 10.5281/zenodo.13821399
[57] Blaabjerg L M, Jonsson N, Boomsma W, Stein A and Lindorff-Larsen K 2024 SSEmb: a joint embedding of protein sequence and structure enables robust variant effect predictions Nat. Commun. 15 9646

## Figure captions
Figure 1. Protein AI model prediction landscape. (a) workflow and (b) dataset composition; (c) PCA of the variant-by-model prediction matrix; (d) median assay-level Spearman correlation heatmap; (e) hierarchical clustering of models by prediction correlation.
Figure 2. Structural-confidence gradient in disagreement (all primary statistics residue-level). (a) distribution of within-protein residue-level Spearman between pLDDT and evolution-vs-structure disagreement (median -0.04, red dashed; 67.5% negative); (b) effect sizes: lowest-vs-highest pLDDT decile difference and enrichment odds ratio; (c) four disagreement definitions (left) and 1,000 balanced family resamplings (right; red line = median -0.05); (d) disagreement in curated disordered vs ordered residues.
Figure 3. Disagreement regimes (K = 6). (a) family-prediction scatter colored by regime; (b) pLDDT composition per regime; (c) regime centroids vs experimental Y.
Figure 4. Reproducibility and case studies. (a) position-level regime agreement vs chance across independent assay pairs, with inset: regime-level experimental Y correlation; (b) TP53, (c) BRCA1, (d) PTEN AlphaFold structures colored by median residue disagreement.
Figure 5. Relative family expertise (percentile space). (a) family advantage across pLDDT bins (flat); (b) family advantage by residue mechanics (flat); (c) family advantage in disordered vs ordered residues (IDR-specific).
Figure 6. Clinical validation. (a) transfer AUROC on the full clinical set and the strict non-overlap set; (b) disagreement vs classification outcome.
Supplementary figures: S1 disagreement-error distribution; S2 disagreement deciles; S3 domain-boundary analysis; S4 BioGate CV performance; S5 BioGate ablations; S6 leave-one-model-out sensitivity; S7 u- vs z-scale comparison; S8 disagreement-error association by assay selection type; S9 residue mechanics (null result); S10 experimental-structure gain by assay pLDDT; S11 SSEmb distinct-architecture replication.

## Table captions

Tables 1-4 are main-text tables; Tables S1-S8 are supplementary.

**Table 1. Disagreement regimes (GMM, K = 6).** Family centroids are z-scale
means; Y = experimental deleteriousness (1 - DMS rank; 0 = tolerated,
1 = damaged); regime names are descriptive labels from the actual centroids
(machine labels regime_0-regime_5).

| regime | n | S_seq | S_evo | S_struct | median Y | median pLDDT |
|---|---|---|---|---|---|---|
| consensus_damaging | 104,108 | 0.97 | 1.05 | 0.94 | 0.81 | 95.9 |
| consensus_tolerant | 97,584 | -0.96 | -1.13 | -0.79 | 0.26 | 90.7 |
| structure-dissenting | 83,694 | 0.61 | 0.65 | 0.04 | 0.65 | 93.4 |
| mild-tolerant | 186,482 | -0.44 | -0.49 | -0.24 | 0.37 | 92.8 |
| high-disagreement | 67,355 | -0.17 | -0.06 | -0.32 | 0.45 | 90.3 |
| mild-damaging | 157,088 | 0.19 | 0.21 | 0.40 | 0.58 | 94.4 |

**Table 2. Relative family advantage by biological context (percentile
space).** Advantage = median error of the other two families minus the
family's own error; positive = family more accurate; protein-level
aggregation. No context shows a consistent advantage except IDR: seq
-0.043 (IDR) vs +0.038 (ordered); struct +0.060 vs +0.005; evo +0.016 vs
-0.011. pLDDT bins are flat (seq +0.008 to +0.022; evo +0.002 to +0.005;
struct -0.002 to -0.007). Mechanics (core/surface/helix/sheet/loop): all
within +/-0.013.

**Table 3. Clinical transfer AUROC (median per protein).** Full set
(n = 1,530 of 2,525 benchmark genes with complete gate inputs): uniform
0.877, BioGate 0.944, best single 0.963. Strict non-overlap set,
gate-complete subset (n = 788): uniform 0.867, BioGate 0.899, best single
0.974 (reselected within the strict set). Uniform AUROC over all 2,490
strict-set proteins: 0.905 (identical at 30% / 50% / 70% identity
thresholds; Table S4).

**Table 4. BioGate 5-fold UniProt-grouped cross-validation (median
protein-level Spearman, percentile space).** best single 0.473; linear
stacking 0.528; uniform ensemble 0.542; BioGate 0.551; XGBoost
(scores + context) 0.573 (pooled out-of-fold protein median, the
aggregation used in this table; the per-fold-median aggregation of
Table S5 yields 0.551 for the same configuration). The median
within-protein paired difference
(BioGate - uniform) was 0.000 (bootstrap 95% CI -0.008 to +0.008).
XGBoost ablations (Table S5): scores only 0.509, context only 0.167,
scores + context 0.551.

**Table S1. Residue-level pLDDT association across disagreement
definitions.** D_std -0.078, D_MAD -0.076, D_pair -0.089, D_u -0.115
(median rho; all CIs exclude zero; 71-77% of proteins negative; four
complementary definitions). Effect sizes (D_evo_struct): lowest-vs-highest
pLDDT decile difference +0.057 (CI +0.013 to +0.084); Cohen's d 0.15
(CI 0.04-0.23); enrichment OR (low-pLDDT decile vs top-disagreement decile)
median 1.84 (CI 1.62-2.04), pooled 1.49 (CI 1.35-1.65).

**Table S2. Confound controls for the shared AlphaFold input.** (i)
Asymmetry: residue-level rho(seq-vs-struct, pLDDT) = +0.01 vs
rho(evo-vs-struct, pLDDT) = -0.04 (65% of proteins in the evo-struct
direction); (ii) structure-model DMS correlation by pLDDT bin (ESM-IF1:
0.04 at <50, 0.47 at >=90; MIF: 0.03 vs 0.43; ProteinMPNN: 0.01 vs 0.22);
(iii) structure-model score coverage = 100% at all pLDDT levels.

**Table S3. Regime robustness across K = 4-7.** Consensus regimes track
experimental Y at every K (Y_tolerant 0.26-0.30; Y_damaging 0.73-0.79);
cross-assay position-level agreement 0.75-0.84 at every K; the
structure-dissenting pattern emerges at K >= 6; ARI vs K = 6: 0.35 (K = 4),
0.58 (K = 5), 0.46 (K = 7).

**Table S4. Identity-threshold sensitivity.** 30% / 50% / 70% each exclude
the same 35 clinical proteins; uniform AUROC 0.905 (n = 2,490) in all
three settings.

**Table S5. XGBoost ablations (5-fold UniProt-grouped CV, percentile
space; median of per-fold medians).** scores only 0.509; context only 0.167; scores + context 0.551;
uniform ensemble 0.542.

**Table S6. Model panel composition.** Per core model: architecture, input
modalities, AlphaFold-structure input, mechanistic-panel membership
(model_panel_table.csv).

**Table S7. Clinical paired deltas (10,000 protein bootstrap).** Full set:
BioGate - uniform +0.039 [0.033, 0.045], best - uniform +0.056 [0.048,
0.064]. Strict set: BioGate - uniform +0.021 [0.011, 0.026], best -
uniform +0.085 [0.077, 0.091], best - BioGate +0.061 [0.056, 0.066].

**Table S8. Held-out protein regime generalization.** Median regime-profile
Spearman 0.94 (IQR 0.89-1.00, 100% positive, 185 proteins); dissenting
regime held-out Y median 0.61, damaging in 72% of supported proteins;
cross-assay replication on held-out proteins median 0.94 (28 pairs,
100% above chance)."""

# write markdown version
md = (f"# {TITLE}\n\n"
      f"## Abstract\n\n{ABSTRACT}\n\n"
      f"**Keywords:** {', '.join(KEYWORDS)}\n\n"
      f"{BODY}\n")

# --- IOP sequential citation renumbering (numbers follow first citation) ---
_body, _rest = md.split("## References", 1)
_refsec, _tail = _rest.split("## Figure captions", 1)
_order = []
for _m in __import__("re").finditer(r"\[(\d+(?:,\s*\d+)*)\]", _body):
    for _x in __import__("re").split(r",\s*", _m.group(1)):
        _x = int(_x)
        if _x > 0 and _x not in _order:
            _order.append(_x)
_old_to_new = {_old: _i + 1 for _i, _old in enumerate(_order)}
_refs = {}
for _m in __import__("re").finditer(r"^\[(\d+)\] (.*)$", _refsec, __import__("re").M):
    _refs[int(_m.group(1))] = _m.group(2)
_new_refsec = "\n".join(f"[{_new}] {_refs[_old]}" for _new, _old in enumerate(_order, 1))
def _repl(_m):
    _nums = [int(_x) for _x in __import__("re").split(r",\s*", _m.group(1)) if int(_x) > 0]
    return "[" + ",".join(str(_old_to_new[_n]) for _n in _nums) + "]"
_new_body = __import__("re").sub(r"\[(\d+(?:,\s*\d+)*)\]", _repl, _body)
md = _new_body + "## References\n" + _new_refsec + "\n\n## Figure captions" + _tail

# --- order Table captions by number ---
_cap_split = md.split("## Table captions")
if len(_cap_split) == 2:
    _cap_blocks = re.split(r"(?=\*\*Table )", _cap_split[1])
    _head = _cap_blocks[0]
    _ok = {"Table 1": 1, "Table 2": 2, "Table 3": 3, "Table 4": 4,
           "Table S1": 5, "Table S2": 6, "Table S3": 7, "Table S4": 8,
           "Table S5": 9}
    _sorted = sorted(_cap_blocks[1:], key=lambda _b: _ok.get(_b.split(".")[0].strip("* "), 99))
    md = _cap_split[0] + "## Table captions" + _head + "".join(_sorted)
(OUT / "manuscript.md").write_text(md, encoding="utf-8")

cover = f"""# Cover letter -- AI for Science (AI4S)

Dear Editors,

We are pleased to submit our manuscript entitled:

**"{TITLE}"**

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
Tables S1-S8; main-text Tables 1-4 are typeset in the manuscript and main
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
"""
(OUT / "cover_letter.md").write_text(cover, encoding="utf-8")

(OUT / "figures_list.md").write_text(
    "Figure/Table caption list: see the captions sections in manuscript.md. Main figures 1-6: figures/main/. Supplementary S1-S11: figures/supplementary/.\n",
    encoding="utf-8")
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

doc = Document()
style = doc.styles["Normal"]
style.font.name = "Times New Roman"
style.font.size = Pt(11)

_p = doc.add_paragraph(TITLE)
_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
for _r in _p.runs:
    _r.font.size = Pt(14)
    _r.font.bold = True

_p = doc.add_paragraph("[Author names and affiliations to be added]")
_p.alignment = WD_ALIGN_PARAGRAPH.CENTER

doc.add_paragraph("Abstract").bold = True
for _para in ABSTRACT.split("\n\n"):
    doc.add_paragraph(_para)
doc.add_paragraph("Keywords: " + ", ".join(KEYWORDS)).bold = True


def _beautify_equation(text):
    return (text
            .replace("sqrt(", "\u221a(")
            .replace("*", "\u00b7")
            .replace("Phi^-1", "\u03a6\u207b\u00b9")
            .replace("sum_m", "\u03a3_m")
            .replace("sum_{m<k}", "\u03a3_{m<k}"))


def add_rich_paragraph(doc, text, bold=False, center=False):
    _p = doc.add_paragraph()
    if center:
        _p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _parts = text.split("**")
    for _i, _part in enumerate(_parts):
        if not _part:
            continue
        _run = _p.add_run(_part)
        _run.font.bold = bold or (_i % 2 == 1)
    return _p


def add_markdown_table(doc, lines):
    _rows = []
    for _line in lines:
        _cells = [_c.strip() for _c in _line.strip().strip("|").split("|")]
        if set(_cells) <= {"---", ":", ":--", "--:", ":-:", ""}:
            continue
        _rows.append(_cells)
    if not _rows:
        return
    _ncols = max(len(_r) for _r in _rows)
    _table = doc.add_table(rows=len(_rows), cols=_ncols)
    _table.style = "Table Grid"
    for _i, _r in enumerate(_rows):
        for _j in range(_ncols):
            _cell = _table.cell(_i, _j)
            _cell.text = _r[_j] if _j < len(_r) else ""
            if _i == 0:
                for _para in _cell.paragraphs:
                    for _run in _para.runs:
                        _run.font.bold = True
    doc.add_paragraph("")


_lines = BODY.split("\n")
_i = 0
while _i < len(_lines):
    _line = _lines[_i].strip()
    if not _line:
        _i += 1
        continue
    if _line.startswith("#"):
        doc.add_heading(_line.lstrip("# ").strip(),
                        level=1 if _line.startswith("## ") else 2)
        _i += 1
        continue
    if re.match(r"^Figure [1-6]\. ", _line):
        n = _line.split(".")[0].strip().split(" ")[1]
        cdir = ROOT / "figures" / "composed"
        cf = cdir / f"Fig{n}_full.png"
        if cf.exists():
            doc.add_picture(os.fspath(cf), width=Inches(6.0))
        add_rich_paragraph(doc, _line)
        _i += 1
        continue
    if _line.startswith("Figure 1.") or _line.startswith("Figure 2.") \
       or _line.startswith("Figure 3.") or _line.startswith("Figure 4.") \
       or _line.startswith("Figure 5.") or _line.startswith("Figure 6."):
        n = _line.split(".")[0].strip().split(" ")[1]
        cdir = ROOT / "figures" / "composed"
        cf = cdir / f"Fig{n}_full.png"
        if cf.exists():
            try:
                doc.add_picture(os.fspath(cf), width=Inches(6.0))
            except Exception:
                pass
        add_rich_paragraph(doc, _line)
        _i += 1
        continue
    if _line.startswith("|"):
        _block = [_lines[_i]]
        _j = _i + 1
        while _j < len(_lines) and _lines[_j].strip().startswith("|"):
            _block.append(_lines[_j])
            _j += 1
        add_markdown_table(doc, _block)
        _i = _j
        continue
    _block = [_line]
    _j = _i + 1
    if not _lines[_i].startswith("  "):
        while _j < len(_lines):
            _nxt = _lines[_j].strip()
            if not _nxt or _nxt.startswith("#") or _nxt.startswith("|") \
               or _lines[_j].startswith("  "):
                break
            _block.append(_nxt)
            _j += 1
    _para_text = " ".join(_block)
    if _lines[_i].startswith("  "):
        add_rich_paragraph(doc, _beautify_equation(_para_text), center=True)
    else:
        add_rich_paragraph(doc, _para_text)
    _i = _j


doc.save(str(OUT / "manuscript.docx"))

# --- pandoc pass: LaTeX math -> native Word OMML, embed figures ---
try:
    import pypandoc
    import csv as _csv
    def _csv_rows(_f):
        with open(ROOT / "results" / "statistics" / _f, encoding="utf-8") as _fh:
            return list(_csv.DictReader(_fh))
    def _f3(_v):
        try:
            _x = float(_v)
        except (TypeError, ValueError):
            return str(_v)
        if abs(_x - round(_x)) < 1e-12 and abs(_x) < 1e6:
            _i2 = int(round(_x))
            return f"{_i2:,}" if abs(_i2) >= 1000 else str(_i2)
        return f"{_x:.3f}"

    _fam2 = {"seq": "seq", "evo": "evolution", "evolution": "evolution",
             "struct": "structure", "structure": "structure"}
    _ctx_order = ([("pLDDT " + b) for b in ("bin",)] * 0 +
                  ["pLDDT<50", "50-70", "70-90", ">=90",
                   "IDR", "ordered", "core", "surface", "helix", "sheet", "loop"])
    _adv = {}
    for _r in _csv_rows("family_advantage_plddt.csv"):
        _adv[(_r["plddt_bin"], _fam2[_r["family"]])] = float(_r["median_adv"])
    for _r in _csv_rows("disorder_family_advantage.csv"):
        _adv[(_r["context"], _fam2[_r["family"]])] = float(_r["median_adv"])
    for _r in _csv_rows("mechanics_family_advantage.csv"):
        _adv[(_r["context"], _fam2[_r["family"]])] = float(_r["median_adv"])
    _t2 = ["**Table 2 (data).** Median percentile-space family advantage by context:", "",
           "| Context | seq | evolution | structure |", "|---|---|---|---|"]
    for _c in _ctx_order:
        _vs = [_adv.get((_c, _fa)) for _fa in ("seq", "evolution", "structure")]
        if any(v is None for v in _vs):
            continue
        _t2.append("| " + _c + " | " + " | ".join(f"{v:.3f}" for v in _vs) + " |")
    _t2.append("")
    _t2.append("*Positive = family more accurate than the other two (median of "
               "protein-level advantages). Contexts not captured by a CSV are omitted.*")
    _t2 = "\n".join(_t2)

    import statistics as _st2
    def _med(_rows, _k):
        return _st2.median(float(r[_k]) for r in _rows)
    _cf = _csv_rows("clinical_transfer_performance.csv")
    _cs = _csv_rows("clinical_transfer_nonoverlap.csv")
    _t3 = "\n".join([
        "**Table 3 (data).** Median per-protein clinical-transfer AUROC:", "",
        "| Set | n proteins | uniform | BioGate | best single |",
        "|---|---|---|---|---|",
        f"| Full clinical benchmark, gate-complete | {len(_cf)} | "
        f"{_med(_cf, 'auroc_uniform'):.3f} | {_med(_cf, 'auroc_biogate'):.3f} | "
        f"{_med(_cf, 'auroc_best'):.3f} |",
        f"| Strict non-overlap, gate-complete | {len(_cs)} | "
        f"{_med(_cs, 'auroc_uniform'):.3f} | {_med(_cs, 'auroc_biogate'):.3f} | "
        f"{_med(_cs, 'auroc_best'):.3f} |",
        "| Strict non-overlap, all proteins (uniform only) | 2,490 | 0.905 | - | - |"])

    _t4rows = _csv_rows("biogate_performance.csv")
    _t4 = ["**Table 4 (data).** BioGate 5-fold UniProt-grouped cross-validation "
           "(protein-level, percentile space):", "",
           "| Method | proteins | median Spearman | mean Spearman | mean AUROC |",
           "|---|---|---|---|---|"]
    for _r in _t4rows:
        _t4.append(f"| {_r['method']} | {_r['n_proteins']} | "
                   f"{float(_r['median_rho']):.3f} | {float(_r['mean_rho']):.3f} | "
                   f"{float(_r['mean_auroc']):.3f} |")
    _t4 = "\n".join(_t4)
    _chunks = md.split("\n\n")
    for _tn, _tmd in (("**Table 2.", _t2), ("**Table 3.", _t3), ("**Table 4.", _t4)):
        for _ci, _ck in enumerate(_chunks):
            if _ck.startswith(_tn):
                _chunks.insert(_ci + 1, _tmd)
                break
    md = "\n\n".join(_chunks)

    t = md
    cdir = ROOT / "figures" / "composed"
    _pbr = ('```{=openxml}\n'
            '<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n'
            '```')
    _blocks = ["## Figures (embedded)"]
    for _gn in sorted(int(_m.group(1)) for _fn in os.listdir(cdir)
                      for _m in [re.match(r"Fig(\d+)_full\.png$", _fn)] if _m):
        _blocks.append(_pbr)
        _blocks.append(f"### Figure {_gn}")
        _blocks.append(f"![]({(cdir / f'Fig{_gn}_full.png').as_posix()})")
    img_lines = "\n\n".join(_blocks)
    t2 = t.replace("## Figure captions",
                   img_lines + "\n\n## Figure captions")
    t2 = re.sub(r"\n(?=\*\*Table )", "\n\n", t2)
    t2 = re.sub(r"\n(?=\[\d+\] )", "\n\n", t2)          # reference items
    t2 = re.sub(r"\n(?=Figure \d+\.)", "\n\n", t2)       # figure captions
    t2 = re.sub(r"\n(?=\*\*Table )", "\n\n", t2)          # table captions
    t2 = re.sub(r"\n(?=Supplementary figures:)", "\n\n", t2)
    tmp = OUT / "_pandoc_src.md"
    tmp.write_text(t2, encoding="utf-8")
    (OUT / "manuscript.md").write_text(md, encoding="utf-8")
    pypandoc.convert_file(str(tmp), "docx", format="markdown+tex_math_dollars",
                          outputfile=str(OUT / "manuscript.docx"),
                          extra_args=["--resource-path=" + str(ROOT)])
    tmp.unlink()
    print("docx re-rendered via pandoc (OMML math, embedded figures)")
    from docx import Document as _Doc
    from docx.shared import Pt as _Pt
    from docx.oxml.ns import qn as _qn
    _d = _Doc(str(OUT / "manuscript.docx"))
    _n = _d.styles["Normal"]
    _n.font.name = "Times New Roman"
    _n.font.size = _Pt(11)
    _n.element.rPr.rFonts.set(_qn("w:eastAsia"), "SimSun")
    for _sn in ("Heading 1", "Heading 2", "Heading 3"):
        try:
            _st = _d.styles[_sn]
            _st.font.name = "Times New Roman"
            _st.element.get_or_add_rPr().get_or_add_rFonts().set(_qn("w:eastAsia"), "SimHei")
        except KeyError:
            pass
    _d.save(str(OUT / "manuscript.docx"))
    print("docx font set to Times New Roman")
except Exception as _e:
    print("pandoc pass skipped:", _e)
print("submission package v1.3 written to", OUT)
for _f in sorted(OUT.iterdir()):
    print(" ", _f.name, _f.stat().st_size)
