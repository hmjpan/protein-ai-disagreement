# GO / NO-GO 1 -- Disagreement vs prediction error

- Date: 2026-09-06 20:49
- Assays analysed: 215
- Median assay-level Spearman(D_std, error): 0.0315
- UniProt-cluster bootstrap 95% CI: [0.0228, 0.0465]
- Fraction of assays with rho > 0: 0.707
- Decile 10 vs 1 mean error: 0.2144 vs 0.2006
- pLDDT-stratified median rho: <50=0.040; 50-70=0.027; 70-90=0.036; >=90=0.029
- Leave-one-model-out median rho range: [0.0271, 0.0391]
- D_MAD sensitivity median rho: 0.0255
- u-basis disagreement median rho: 0.1447 (methodological note)
- Pre-registered rule (protocol s60): PIVOT only if median rho < 0.05 AND bootstrap CI crosses 0

**VERDICT: GO**

NOTE: median rho 0.0315 is below the 0.05 effect-size target even though the CI excludes 0; disagreement-as-uncertainty is reported as a weak, consistent signal (no 'strongly predicts' wording).