# Independent review: M2 development evaluation harness and first Tier A dev report

Reviewer: independent agent (Opus 5.5). I did not write the harness. Date 2026-10-28.
Reviewed: harness at commit 17b69c5 (`src/slparalog/eval/*`, `scripts/eval/dev_report.py`, `tests/test_eval_harness.py`,
`config/models.yaml`); report and log at cf3509a (`results/dev/dev_primary_tierAprev_17b69c5.json/.md`,
`results/experiments.tsv`). Spec: DESIGN v3 §5–§6, RESEARCH_TREE D8–D10, `config/success.yaml`,
`reviews/round2/statistics.md` (items 6, 11, 13, 15).
Data read: `data/labels/dev/{human,yeast}_dev.tsv` and `data/interim/features/{human,yeast}/tierA_prev_prefiltercap.parquet`
(the inputs recorded in the report JSON). I did not read `data/labels/sealed/`, and I computed no score on any non-dev pair.
Reviewer scripts and outputs: `data/interim/review_eval/work/` (gitignored). They were run from the session scratchpad,
under a load average of about 147. Every number below comes from those outputs or from the committed report JSON
(`results/dev/dev_primary_tierAprev_17b69c5.json`, abbreviated **RJ**).

## Verdict: ACCEPT WITH FIXES

I found no implementation bug that would bias a decision. The harness computes what the spec says. I reproduced the
key numbers with my own code: the B1 variants, LR dev AP, T, its bootstrap CI and the ≥3-studies stratum.

Two findings matter for how the results are read. Neither is a code bug.
- **(M1)** The pre-registered "pooled OOF" T is about half the paired per-fold ΔAP: LR T is 0.026 vs 0.058. The cause is
  not the weighting scheme, and is otherwise unresolved.
- **(M2)** The human gain over B1 disappears once ascertainment strata are respected.

The remaining fixes are minor: a degenerate B3, G5 not exercised, binning documentation, and missing R2 item-13 strata.
The report was also run on the pre-fix caller tables and with a dirty `registry.py`.

## 1. Independent re-computation (key numbers)

These numbers come from my own pipeline: own label parsing, unordered-pair join, weights, LR pipeline, fold loop and
bootstrap. Only `species_balanced_weights` and `GBDTModel` were imported, for direct comparison.

| quantity | RJ | reviewer |
|---|---|---|
| primary population human / yeast (pos/neg) | 4415 (123/4292) / 996 (60/936) | same |
| B1 `ident_global` AP human / yeast | 0.118 / 0.567 | 0.1183 / 0.5671 |
| all 7 B1 variants (both species) | table in RJ `b1_variants` | identical to 3 dp; `ident_global` max mean AP |
| LR_00 (C=0.01) AP human / yeast | 0.179 / 0.558 | 0.1793 / 0.5582 |
| T (LR_sel vs B1), point | 0.0261 | 0.0261 |
| T 95% CI, SE (B=10,000 in RJ; B=2,000 own RNG) | [−0.011, 0.072], 0.0211 | [−0.0106, 0.0709], 0.0206 |
| human ΔAP CI (LR_sel) | [0.025, 0.108] | [0.026, 0.110] |
| LR_01 (C=0.1) T (the family's best) | 0.033 | 0.0335 |
| GBDT_17 / GBDT_09 per-seed AP, seeds 11 and 23 (`per_seed_ap`) | e.g. 0.151/0.5398, 0.1353/0.5641 (s11) | identical to 4 dp |
| human n_studies ≥3 stratum: n, pos, B1, LR | 356, 75, 0.451, 0.446 | 356, 75, 0.451, 0.446 |
| human n_studies = 2 stratum: B1, LR | 0.143, 0.084 | 0.143, 0.084 |
| bootstrap stratum group counts | h0 820, h1 88, m0 26, m1 24, y0 202, y1 40 | same |

The bootstrap CI matches to within Monte-Carlo error, from an independent draw scheme and sklearn's weighted AP.

## 2. Checklist

1. **OOF.** `cv.oof_scores` trains on `folds != k` and predicts `folds == k` for each seed's stored `fold_s<seed>` column.
   I checked that every group sits in exactly one fold for all 5 seeds; `check_folds` also asserts this. Each pair is
   predicted once per seed and the 5 predictions are averaged. Seed-averaging gives no material ensemble gain: LR_00
   mean per-seed pooled T is 0.0288 vs 0.0261 seed-averaged. Pooled vs per-fold: see **M1**.
2. **Weights.** Verified numerically on seed 11, fold 0 training rows (N = 3,487):
   - each species has total weight 1,743.5;
   - weighted prevalence is 0.048576 in both species (raw prevalence 0.0368 human, 0.0604 yeast);
   - the weights sum to N;
   - an independent formula agrees to 0.0.

   Single-species (LOSO) training gives uniform weights.
3. **Group bootstrap.**
   - Strata = composition × has_pos over the primary population.
   - Within each stratum, groups are drawn with replacement, as many as the stratum holds.
   - Every statistic uses the same multiplicity matrix, so Δ and T are paired.
   - CIs are percentile.
   - Seed 11 and B = 10,000 are read from `config/success.yaml`.
   - The fast AP equals sklearn AP with integer weights. It is re-checked at run time on 2 resamples per call.
4. **B1.**
   - Selected on dev data only (max species-mean AP, ties in config order).
   - Scored on the same rows as the model: the same `df` rows are used, and there are 0 NaN in `ident_global` and in all
     50 features. B1 variants are unfitted.
   - The explicit G5 guard is not called (minor issue m2).
5. **B2/B3/B4.**
   - B2 = LR on (B1 variant, log1p fam_mc_size, rbh), OOF on the same folds. Correct.
   - B4 = `n_testable` as a raw score. Correct.
   - B3 is degenerate (m1).
6. **GBDT binning.** Correct as described. Measured effect vs sklearn 1.9.1 native weighted binning is small (m3).
7. **One-SE rule.**
   - Implemented as declared: threshold T_best − bootstrap SD(T_best), then the simplest by the declared order, with
     −T as tie-break.
   - GBDT_17 is indeed the simplest of the 29 eligible configs.
   - See o1.
8. **LOSO.**
   - Trains on species X with `folds != k` and tests species Y with `fold == k` (R2 item 14), so mixed groups (50)
     never cross.
   - A test covers it.
9. **Controls.**
   - **Shuffled labels:** permuted within species and fully refit. AP is close to prevalence: human 0.027 vs 0.028.
     Yeast GBDT is 0.071 vs 0.060, inside the B0 per-seed range 0.058–0.074 (RJ `b0_per_seed_ap`).
   - **Species-ID and source-library:** both are grouped CV with unweighted GBDT and computed sensibly. Results
     (interpretation in M2):
     - species-ID AUROC 0.957 [0.879, 1.000];
     - source-library AUROC 0.928 [0.885, 0.955]: 61% of positives vs 6.5% of negatives are in ≥3 studies.
10. **≥3-studies stratum.** Reproduced (table above). Subset AP on pooled OOF scores. The ΔAP CI uses the same
    bootstrap, restricted to the subset's groups.

## 3. Issues

### M1 (moderate; specification, cause unresolved): pooled-OOF T is about half the paired per-fold ΔAP
`config/models.yaml` and R2 item 15 define T on OOF scores pooled across the 5 folds of 5 separately fitted models. B1
needs no fit, so its pooled ranking is internally consistent. The model's pooled ranking also mixes between-fold
differences in score level. LR_00 on the same pairs and folds:

| | human ΔAP | yeast ΔAP | T |
|---|---|---|---|
| pooled OOF (as reported) | 0.061 | −0.009 | 0.026 |
| mean over 25 seed×fold test sets of ΔAP (paired, same fold) | 0.095 (0.235 − 0.140) | 0.022 (0.608 − 0.586) | 0.058 |
| LR_01 (C=0.1) pooled / per-fold | 0.060 / 0.108 | 0.007 / 0.023 | 0.033 / 0.066 |

The per-fold mean is itself a different estimand. AP on small folds (about 12 yeast and 24 human positives) is biased,
and folds are equally weighted. That bias largely cancels in a paired Δ, but not exactly.

**What it is not.** In every seed, human fold 0 holds the forced-to-dev 0-positive family, so fold-0 test prevalence
is 0.014 against 0.037 in the other folds. I suspected the per-fold weighted target π* as the cause, since the training
prevalence is 0.0368 when fold 0 is held out vs 0.0262 otherwise, and π* sets the model's intercept. Refitting LR_00
with a fixed π* (0.04405, the dev mean) for every fold leaves the gap unchanged:
- pooled T = 0.0248 (human AP 0.1761, yeast 0.5588);
- per-fold human 0.234 − 0.140, yeast 0.608 − 0.586, so per-fold T is about 0.058 as before.

So the gap is not caused by the weighting scheme. It is either other fold-to-fold score-level shifts, or simply that the
two estimands differ: AP is not additive over folds, and the model may gain more in some folds than others.
**I have not resolved which.**

Consequences, if the pooled figure is the biased one:
- The dev T understates what a single model scoring one held-out set would show.
- The one-SE rule and the RESEARCH_TREE stop rule both use this T, so they are biased toward rejecting ideas and toward
  configs whose score level is stable across folds.
- This does not inflate any claim, and the test-time T (one model, one test set) is unaffected.

Recommendation for the PI:
- Report the paired per-fold ΔAP, and its T, next to pooled T in every dev report.
- Consider whether selection should use it. That needs a protected-config change, so it is the PI's call.
- Do not change anything silently.
- A cheap diagnostic to settle the cause: a simulation with a known-good fixed score plus fold-specific monotone
  distortions, or pooled AP after within-fold quantile normalisation of the model scores.

### M2 (major; interpretation, harness correct): human ΔAP vs B1 is not robust to ascertainment stratification
From RJ `strata.D8_human`, AP by number of testable studies:

| testable studies | B1 | LR_sel | GBDT_sel | B2 |
|---|---|---|---|---|
| 1 | 0.009 | 0.009 | 0.011 | |
| 2 | 0.143 | 0.084 | 0.076 | |
| ≥3 | 0.451 | 0.446 | 0.384 | |
| stratum-averaged | 0.201 | 0.180 | 0.157 | 0.203 |

- In the 2-studies stratum, LR_sel ΔAP is −0.059 [−0.162, −0.001] and GBDT_sel ΔAP is −0.066 [−0.174, −0.006].
- Together with source-library AUROC 0.928, the human pooled gain (+0.061 LR) looks like learned ascertainment (which
  pairs were heavily screened), not SL signal beyond identity.
- Yeast shows no gain.

This is what D8 was designed to catch. The harness surfaces it correctly. It should be stated explicitly in the dev
summary, not left to the strata list.

### m1 (minor): B3 is degenerate by construction, and its numbers are a fold artifact
`b3_diagnostics.frac_pairs_with_gene_seen_in_training = 0.0`. Split groups are gene-connected components, so a test
gene never has training pairs. Every B3 score is therefore the species' training-fold prevalence, a per-fold constant.
Its AP and AUROC (human AUROC 0.457 < 0.5) reflect fold-prevalence anti-correlation only. **Fix:** report B3 as "not
informative: 0% gene overlap by construction", or drop its AP from the tables.

### m2 (minor): the G5 pair-set guard is not exercised
`metrics.ap(..., reference_pair_ids=...)` and `ap_pair` are never called by `dev_report.py`. Identical pair sets hold
here by construction: same rows, and B1 has no NaN. Otherwise `method_stats` would raise. **Fix:**
- assert `np.isfinite(S["B1"]).all()`, or route the model-vs-B1 point estimates through `ap_pair`;
- `final_test.py` must use the guard.

### m3 (minor): GBDT unweighted-quantile pre-binning
sklearn 1.9.1 uses weighted percentiles for bin edges when `sample_weight` is given (`binning.py`
`_find_binning_thresholds`). The wrapper's bin edges are the unweighted ones. The test checks only the unweighted case,
which is exact.

Measured against native weighted binning, same folds, weights and seed:

| config, seed | T (pre-binned) | T (native) | Spearman of predictions |
|---|---|---|---|
| GBDT_17, s11 | 0.0027 | −0.0007 | 0.989 |
| GBDT_17, s23 | 0.0041 | 0.0056 | 0.990 |
| GBDT_09, s11 | 0.0070 | 0.0030 | 0.958 |
| GBDT_09, s23 | 0.0079 | 0.0059 | 0.945 |

The largest per-species AP difference is 0.014, and single predictions differ by up to 0.35.
- |ΔT| ≤ 0.004 is far below SE(T) ≈ 0.021, so this does not matter for decisions now.
- It is a deviation from "sklearn HGB with sample weights". The docstring should say so; currently it says only
  "exactly the binning an unweighted fit uses".
- The size of the prediction differences shows that GBDT seeds and bins carry noise comparable to inter-config gaps.

### m4 (minor): R2 item 13 strata missing; R2 item 6 BCa check missing; ribosomal strata unavailable
- The source-library control is positive (CI excludes 0.5). R2 item 13 then requires T within n_lines_tested strata
  (2–3, 4–10, >10). Only the D8 n_studies strata (1, 2, ≥3) are reported.
- The BCa check (R2 item 6) is not implemented.
- The ribosomal strata are "unavailable" because the UniProt KW-0689 files are missing. The D10 pre-declaration needs
  them.

### m5 (minor; provenance)
- RJ `git_dirty_paths` = `src/slparalog/features/registry.py`. The uncommitted diff only adds Tier B entries; the Tier A
  feature list is unaffected, and `registry_hash` is recorded.
- The inputs are the **pre-fix caller tables** (`tierA_prev_prefiltercap.parquet`). None of these numbers should drive
  adoption until the report is re-run on the fixed caller.
- `config/models.yaml` first appears in the harness commit 17b69c5 (04:38 UTC); the report ran after it. I cannot verify
  from git that no real-data score was computed during development before that commit.

### o1 (observation; per spec): the one-SE rule is nearly non-selective here
SE(T_best) ≈ 0.021 is the marginal bootstrap SD of T. Configs are scored on the same pairs, so their T differences are
much less noisy than that. Result: 29/32 GBDT configs are eligible, and the simplest (T = 0.010) is chosen over the best
(0.026). This follows R2 item 15, but in practice the rule means "pick the simplest". An SE of the paired difference to
the best would be more discriminating; that is the PI's call.

## 4. Not verified
- I did not re-run the full grid end-to-end. I re-ran LR C ∈ {0.01, 0.1} on all seeds and GBDT_17 and GBDT_09 on 2
  seeds. I did not re-run B2, LOSO numbers or the controls; I checked their code paths only.
- The end-to-end block was not re-computed.
- `pytest tests/test_eval_harness.py` was not run, because of machine load.
