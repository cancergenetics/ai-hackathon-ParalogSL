# Adversarial review: Tier B dev gain (GBDT A+B, T_pf 0.143)

Reviewer: an independent agent that did not build the harness or the features. Only dev labels were read. Nothing was committed. Every model was refit with the harness's own `GBDTModel`, using the selected parameters (GBDT_17 = GBDT_AB_17), species-balanced weights and the stored `fold_s*` columns, with all 5 seeds from `config/seeds.yaml`.

- Scripts: `scripts/review/tierB_{common,core,homology,ascertainment,concentration,importance}.py`
- Results: `results/review/tierB_gain/{core,homology,ascertainment,concentration,importance}.json`
- Review bootstraps are stratified by group, with B=2000 (B=1000 for the cell analyses), seed 11.

## VERDICT: PARTLY ARTEFACT

The gain contains real signal: I found no leakage, it survives every control, and it is not concentrated in a few families. But most of the human gain (0.247) is ranking *between* screening-ascertainment cells. When ascertainment is held fixed, the gain over identity shrinks to about 0.01–0.03. Yeast is the species without library ascertainment, and its gain is 0.040 [−0.024, 0.097]. T_pf = 0.143 therefore overstates the transferable gain severalfold.

## 1. Reproduction and controls (core.json)

**Reproduction:** exact. T_pf for A+B is 0.1432 (human 0.2468, yeast 0.0396). For A alone it is 0.0462.

**Tier B rows shuffled within species:** T_pf falls to 0.040 (human 0.082), which is the A level. The gain is carried by Tier B values that are attached to the right pairs.

**Labels shuffled:**

| species | A+B AP over 5 permutations | prevalence |
|---|---|---|
| human | 0.025–0.030 | 0.028 |
| yeast | 0.055–0.084 | 0.060 |

**Tier B alone:** T_pf 0.070 (human 0.147, yeast −0.007).

**Seeds:**
- T_pf for A+B per seed: 0.137, 0.137, 0.138, 0.160, 0.143.
- ΔT (A+B − A): 0.092–0.107. Seed spread is negligible.

## 2. Leakage (homology.json; tierB.parquet.meta.json)

**Self-exclusion works.** It removes HUMAN, PANTR and GORGO for human, and YEAST for yeast. No feature path uses labels.

**Macaque (d = 0.011) is kept for human.** I did not recompute Tier B with `--panel-dmin 0.03`, so this is unverified.

**Homology between top contributors and training positives:**
- No top-20 human contributor to ΔAP (A+B − A) has a homologue with bsr ≥ 0.25 among the same fold's positive training pairs in either species. This used my own MMseqs2 all-vs-all search over the dev genes.
- 25% of the top 20 have any homologous hit at all, against 31% of all positives and 51% of negatives.
- Contribution does not correlate with homology to training positives (Spearman 0.03).

**Stricter regrouping.** I merged groups by any hit (giant group of 3,432 pairs) or by bsr ≥ 0.15:

| regrouping | T_pf A+B | T_pf A |
|---|---|---|
| any hit | 0.150 | 0.048 |
| bsr ≥ 0.15 | 0.151 | 0.035 |

The gain does not depend on fold-level homology, so family memorisation is ruled out.

## 3. Concentration (concentration.json)

Dropping each of the 5 largest positive-holding groups in turn (train and eval) gives human per-fold ΔAP 0.241–0.262. Dropping all 5 together gives 0.262 (baseline 0.247).

## 4. Ascertainment — the decisive part (ascertainment.json, concentration.json, importance.json)

**Within each source library**, primary label:
- A+B − A is positive in 9/9 libraries; the CI excludes 0 in 8/9.
- A+B − B1 is positive in 9/9; the CI excludes 0 in 6/9.

**Within each library, scored against the library's own SL call:**
- A+B − A is positive in 9/9; the CI excludes 0 in 6/9.
- A+B − B1 is positive in 8/9; the CI excludes 0 in 3/9.

**By testable-studies stratum**, A+B − B1 (pooled):

| stratum | ΔAP [CI] |
|---|---|
| 1 | +0.017 [−0.001, 0.093] |
| 2 | +0.010 [−0.096, 0.112] |
| ≥3 | +0.151 [0.043, 0.258] |

In the ≥3 stratum the per-fold Δ is +0.150 for A+B and −0.011 for A.

**Within exact study-set cells.** These hold identical ascertainment: 43 cells, 105 positives.
- A+B − B1 = **+0.010 [−0.057, 0.063]**.
- A+B − A = +0.096 [0.039, 0.135].
- A − B1 = −0.086.
- In the 9 large cells (≥20 pairs, 43 positives), positive-weighted AP is B1 0.155, A 0.118, A+B 0.172.

**Ascertainment covariates regressed out, at score level.** The covariates are n_testable_studies plus 9 library indicators (human, grouped-CV LR):

| model | per-fold ΔAP vs B1 |
|---|---|
| covariates alone | +0.291 |
| covariates + B1 | 0.342 |
| covariates + B1 + A+B | 0.357 |
| covariates + B1 + A + A+B | 0.367 |

- Adding the A+B score therefore adds +0.015 over covariates + B1, and +0.030 over covariates + B1 + A.
- In-sample LR test for the A+B score: χ² = 37, p = 1e-9. Within the ≥3 stratum it adds +0.017.

**Regressed out at feature level, or by training.**
- Tier B residualised on the same covariates: T_pf 0.172.
- D12(c) deconfounded training: 0.142.
- Both leave the gain intact, so the dependence is not linear in the covariates.

**Label-free check: does Tier B predict screening?** Target: whether a pair was screened, among Compara paralog pairs.

| feature set | AUROC, all pairs | AUROC, pairs touching a dev gene |
|---|---|---|
| Tier B | 0.84 | 0.81 |
| A | 0.90 | 0.90 |
| A+B | 0.91 | 0.91 |
| identity | 0.69 | 0.80 |

Tier B is strongly predictive of screening, but it adds almost nothing to A for that task. Used as an SL score, a screening propensity built from Tier B gives human ΔAP +0.107. However, the A+B model adds +0.114 on top of B1 + propensity, which suggests the gain is not simply propensity.

**Leave-one-species-out (yeast → human), per-fold ΔAP vs B1:** A 0.083, A+B 0.144. The yeast-trained A+B score adds +0.007 once the human ascertainment covariates are included.

## 5. Biology (importance.json)

**Held-out permutation importance**, as the drop in AP:

| species | Tier B block | Tier A block |
|---|---|---|
| human | 0.239 | 0.286 |
| yeast | 0.271 | 0.463 |

The top Tier B features are `tb_n_post_single`, `tb_rel_rate` and `tb_cons_min_pct`.

**Effect directions agree across species** (standardised difference, positives vs negatives):
- More 1:2 single-ortholog patterns: `tb_n_post_single` +1.08 human, +1.60 yeast.
- Slower relative rate: `tb_rel_rate` −0.70 human, −0.63 yeast.
- Higher conservation: `tb_cons_min_pct` +0.38 human, +0.28 yeast.

This fits the LITERATURE.md priors: conservation and age as predictors (De Kegel 2021), and asymmetric divergence lowering redundancy.

**"Older duplicates are more often SL" is not supported.** The direction of `tb_depth_pct` reverses between species (human −0.01, yeast −1.36) and between human strata (stratum 1 +0.51, ≥3 −0.22).

**The human signal sits in heavily screened pairs.** `tb_n_post_single` has a standardised difference of 0.93 in the ≥3 stratum but only 0.14 in stratum 1.

## Required before adoption

1. Report the within-cell and covariate-adjusted ΔAP alongside T_pf, and consider a within-cell guard.
2. Treat about 0.01–0.05 as the realistic transferable gain, not 0.14.
3. Rerun Tier B with `--panel-dmin 0.03` (macaque removed).
4. Judge Tier B on the unseen-species test, where library ascertainment is absent.
