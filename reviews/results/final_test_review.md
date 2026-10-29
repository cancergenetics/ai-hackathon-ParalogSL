# Independent review: scripts/eval/final_test.py (commit a388a30) and tests/test_final_test_guard.py

Reviewer: an independent agent that did not write this code. Nothing under `data/labels/sealed/` was read, listed or
opened. Real mode was not run. Code was not modified.
Checked against: DESIGN.md v3 §5, §6.1–6.4 and §8; config/success.yaml; config/labels_unseen.yaml; config/models.yaml;
reviews/round2/statistics.md (items 5, 6, 11, 12) and agent_failure_modes.md (G15); RESEARCH_TREE D8, D10, D12, D17, D18;
src/slparalog/eval/*, scripts/eval/dev_report.py, src/slparalog/train.py and model.py; and
results/dev/dev_primary_AB_ABC_noasym_6bb9438.json, which holds the GBDT_AB_17 selection.

## VERDICT: READY WITH FIXES

The primary statistic is computed as pre-registered and as on dev, and the single-access machinery is sound. Before the
access, three gaps need fixing (R1–R3). Each is cheap to fix now, and none can be fixed after the access:
- R1: nothing checks that the frozen model *is* GBDT_AB_17 on A+B.
- R2: several pre-registered secondaries that need test labels are missing.
- R3: a failure in a secondary analysis discards the primary result, yet still uses up the access.

## What was run

- `pytest tests/test_final_test_guard.py`: **46 passed** (23.6 s).
- `final_test.py --dry-run-toy --out-dir <scratchpad>`: the toy run completed. T = 0.306 [0.050, 0.417] with B = 10,000 and
  seed 11, and the verdict was "CLEARLY BEATS" on the stand-ins. All secondary rows rendered, and the run never touched
  sealed paths or the access log.
- I did not verify anything that needs the real frozen model, tag or sealed files. At review time `results/frozen/` and
  tag `final-model-v1` do not exist.

## Q1 Primary T

| item | status | evidence |
|---|---|---|
| B1 variant | OK in practice, weakly checked | `config/models.yaml b1.selected: ident_global` = dev `b1_selected` (ident_global). `train` writes it to model.json. But `check_model` only tests membership: `if v not in baselines.B1_VARIANTS:` (l.220). |
| Same-pair-set guard | OK | Per species, `ap_pair(mf, ...)` (l.677–678) runs on frames built from one `d`. Model, B1 and B2 must be finite on the whole population (l.817–819). |
| Caller-missed labelled pairs | **Primary is correct; the end-to-end analysis is missing.** | The brief's premise is not what is pre-registered. DESIGN §6.2, statistics R2 item 5 and the dev harness all define the primary population as labelled ∩ called. dataset.py says: "labelled pairs that the caller misses are excluded from it … (they enter only the end-to-end analysis, §6.4)". final_test does the same: `pops[sp]["primary"]` (l.793) and `labelled_pairs_not_called_excluded_from_primary` (l.853). Dev's minimum-score treatment exists only in its `end_to_end` block (dev_report.py l.638–651). final_test has **no** end-to-end block (see R2). |
| Bootstrap | OK | `make_bootstrap(S, B, seed)` with B and seed from success.yaml (l.785, 823). `check_config` pins strata, `ci_type: percentile` and `relative_gain_on`. `species_mean_delta` and `ci` are the dev functions. Resamples where a species is NaN are dropped and counted (`n_dropped`). |
| "Clearly beats" | OK | `rel = mean(AP_m/AP_B1 − 1)` on point estimates (l.693). `clearly = passed and rel_ok`, with the threshold read from config (l.671). |
| Per-species ΔAP > 0 | OK | `each_pos = all(... ["delta_vs_B1"]["point"] > 0 ...)` (l.695). |
| Differences from dev | Minor | (a) `delta_block` uses `min_pos_for_ci=1` (l.663). Dev and statistics R2 item 11 use `controls.min_pos_for_ci: 10`, giving counts only below 10 positives (R5). (b) The study-set-cell CI reuses the full test bootstrap, whereas the D17 review re-bootstrapped the subset. This is acceptable, but state it. (c) `excl0` also accepts `hi < 0` (l.694). This has no effect while `each_pos` is required. |

## Q2 Model and training data

- **Training labels:** `check_model` requires each training label file to equal `data/labels/dev/<name>` by sha256
  (l.228–231), and `train` reads labels through `load_labels`, which refuses sealed paths. Test labels therefore cannot
  enter the frozen model.
- **Yeast-only refit:** uses `dev_dir/"yeast_dev.tsv"` and the model's own params, features and seed, and checks that it
  reproduces them (l.598–602).
- **B2:** fitted on dev pairs only (l.578).
- **Timing:** both refits run inside preflight, before the access record, and are cached (l.531–534).
- **Not verified, R1:** that model.json has `feature_set == "A+B"`, `model_type == "gbdt"`, the GBDT_AB_17 params
  `{l2 1.0, lr 0.03, max_iter 100, leaves 7, min_leaf 50}`, the 73 dev A+B features with no `tb_asym_*`, and
  `training.tier == "primary"`. `train`'s default params happen to equal GBDT_AB_17, but if the PI passes the wrong
  `--feature-set`, every precondition still passes.

## Q3 Secondaries

- **No post-hoc analyses.** Everything present is pre-declared:
  - D17(c), all four analyses;
  - Δ vs B2 (DESIGN §6.2);
  - human strata 1/2/≥3 (D8(b));
  - the Arabidopsis broad tier (labels_unseen.yaml).
- **Unseen species:** scored with the same frozen model, B1 variant and B2 (l.730–732), and the Arabidopsis strict tier
  is listed first. However, `render_md` does not mark strict as the headline (R6).
- **Missing.** These analyses are pre-registered and need sealed labels, so they can never be computed after this access:
  - end-to-end (DESIGN §6.4; statistics item 5 says "co-reported");
  - lenient, strict and ≥2-studies tiers (§6.2; statistics item 11, level 2);
  - leakage tertiles with a trend CI (§6.1; the test files carry `leakage_bsr`);
  - human→yeast LOSO (§6.2 "LOSO ΔAP"; statistics item 11 says "both directions");
  - B4 next to every human AP (D8(a)), and the stratum-averaged AP (D8(b));
  - results without Harle/Chan-dependent labels (D8(c));
  - yeast results without orientation-conflict positives and with/without ribosomal pairs (D10);
  - Arabidopsis results without lethal-single positives (D4).

  See R2.

## Q4 Single access

- **Order of operations is correct.** The order is: preconditions; then preflight (all four species, model, yeast-only
  refit, B1, B2, NaN checks); then an fsync'd `access_started` (l.939); then the sealed reads.
- **Preflight failures cannot use up the access.** Any exception in `RealBackend` or preflight becomes a Refusal before the
  record is written (l.929–938).
- **Crashes count as an access.** `except BaseException` logs `access_failed`. A hard kill leaves `access_started`, and
  `check_no_prior_access` blocks the next run. Both cases are tested.
- **Risks after the access record:**
  - Human and yeast tables are protected before the access by `test_splits.py::test_sealed_manifest_and_protected_hashes`.
    That test strictly checks the manifest against protected.sha256 and the tables against the manifest, and it runs in
    the full suite recorded by the marker.
  - **The unseen tables are verified only inside `read_sealed`.** They are read in `sec["unseen_species"]` (l.848), and
    outputs are written only at the end (l.948). A mismatch, or any secondary exception, therefore logs `access_failed`
    and throws away an already computed primary T (R3).
- **Deviations from DESIGN §6.3 and G15:**
  - The log lives at `results/final/access_log.tsv`, not `results/test_access_log.tsv`.
  - Refusals are logged untracked in the git dir.
  - Log rows lack `sealed_manifest_sha256` (R6).
- **Before the access:** protected.sha256 does not yet list `scripts/eval/final_test.py` or
  `tests/test_final_test_guard.py`, so the script refuses at present. This is correct behaviour; the PI must add both
  entries with a D-note.

## Q5 Sealed path and outputs

- **One sealed reader.** `read_sealed` (l.403–434) is the only function that opens sealed files. Labels reach
  `load_dev_labels` as `raw=`. Other sealed-path uses are git metadata (`check_sealed_committed`) or skips (`check_protected`).
- **No pair-level labels are written.** Outputs are aggregates only (AP, counts, CIs, hashes), so nothing written to
  results/ exposes test labels to later dev work.
- **Unprotected unseen feature tables.** The unseen-species feature tables are gitignored and are not checksummed against
  the tag (R6).

## Required fixes

- **R1. Pin the model in `check_model`.** Assert the following against `results/dev/dev_primary_AB_ABC_noasym_6bb9438.json`
  `selection.gbdt_ab.selected` and the A+B feature list:
  - `meta["feature_set"] == "A+B"` and `meta["model_type"] == "gbdt"`;
  - `meta["params"]` equals the GBDT_AB_17 params;
  - `meta["features"]` equals dev `feature_sets["A+B"]`, with no `tb_asym_*`;
  - `meta["training"]["tier"] == "primary"`;
  - `meta["b1_variant"] == mcfg["b1"]["selected"] == "ident_global"`.
- **R2. Settle the missing pre-registered secondaries before the access.** For each analysis listed under Q3, either add
  it, reusing the dev code (dev_report end_to_end, strata, `strat_avg`, `TIERS`), or have the PI drop it explicitly in a
  D-note. At minimum, add end-to-end, the tiers and the leakage tertiles. They are normative, and after the access they
  are unrecoverable.
- **R3. Protect the primary result from secondary failures.**
  1. Hash-verify every sealed file named in `sealed_names` immediately after `access_started`.
  2. Compute and write the primary block first (JSON, fsync).
  3. Wrap each secondary block in try/except and record the error in the JSON instead of aborting.

## Should

- **R4.** `read_sealed` compares the manifest strictly to protected.sha256 (l.418), but `check_protected` and G3 accept a
  D-note override. Today the strict `test_splits` check closes this gap. Document that dependence, or check that the
  manifest entry exists during preconditions.
- **R5.** Use `controls.min_pos_for_ci` (10), as dev does.
- **R6.** Several smaller fixes:
  - Mark Arabidopsis strict as the headline in the md report.
  - Put the log at the DESIGN path, or record the change in a D-note.
  - Add `sealed_manifest_sha256` to log rows.
  - Record unseen-species feature-table sha256 in model.json at the tag.
  - Report the CI of the relative gain and AP above chance (statistics item 12).
  - Use `T_ci["lo"] > 0`.
- **R7.** Uncommitted working-tree changes (registry `columns_hash`, model.py, train.py) are in flight. `table_failures`
  still compares the global `registry_hash`, so once that lands preflight will refuse until both are aligned. This is safe,
  because it blocks before any access, but the two must be coordinated.

## Re-review of 563b15d

Verdict: **READY**. Checked against `git diff 8ecda54 563b15d` (final_test.py, controls.py, dev_report.py). Sealed data was not
read, and `final_test.py` was run only with `--dry-run-toy`.

- **R1 fixed.** `frozen_spec_failures` (called from `check_model`) compares `meta` against values derived from config and the
  registry:
  - `feature_set == "A+B"` and `model_type == "gbdt"`;
  - `params == models.yaml anchoring.ensemble.gbdt_params`, which I checked equals dev `selection.gbdt_ab.selected.params`
    (GBDT_AB_17);
  - `features == registry A/B model features minus models.yaml excluded_features` (the two tb_asym_* are in that list). I checked
    that this equals dev `feature_sets["A+B"]` (73 features, no `tb_asym_*`);
  - `training.tier == "primary"`;
  - `seeds` equals `config/seeds.yaml`;
  - `b1_variant == models.yaml b1.selected` (ident_global).
  Tests `test_refuses_model_that_is_not_the_frozen_spec` and `test_frozen_spec_follows_config_not_literals` cover these.
- **R2 fixed.**
  - `sec_end_to_end` calls `controls.end_to_end_ap`, which is also used by dev_report.
  - `sec_tiers` loops `dataset.TIERS` (lenient, two_studies, strict, biogrid), each with `load_dev_labels`, `population` and `primary_block`.
  - `sec_leakage_tertiles`, `sec_human_strata_D8` (strat_avg_ap, B4, Harle/Chan, n-lines), `sec_yeast_and_ribosomal_D10`,
    both only-model analyses and the D4 lethal-single block are present.
  - The primary population is built before any secondary and is unchanged. Uncalled pairs are excluded from primary and appear
    only in end-to-end and in `labelled_pairs_not_called_excluded_from_primary`.
  - Nothing outside the pre-registered list was added. The one nuance is the `biogrid` tier, which is in `TIERS` and in dev but
    was not named in the review. Minor.
- **R3 fixed.**
  - `verify_sealed` hashes every `sealed_names` file against the manifest (`read_sealed(verify_only=True)`). It runs before
    `access_started` and refuses without consuming the access.
  - `compute` passes the primary block to `on_primary`, which calls `write_outputs` with `_write_durable` (fsync and atomic
    replace), before any secondary runs.
  - Each secondary runs in its own try/except and records the error in `secondary_errors`. Status becomes
    `complete_with_secondary_errors`, and the dry-run shows this: the toy fixture has no biogrid column, so the toy run records
    `KeyError: 'biogrid'` and finishes. Test `test_sealed_hashes_are_verified_before_the_access_record` covers the hash order.
- **Refactor.**
  - The old dev_report `strat_avg` body and the new `controls.strat_avg_ap` are logically identical. The only change is
    `S.loc[mm, col]` -> `s[mm]` (a numpy array in place of a Series), which does not change `summary`'s AP.
  - The old end-to-end loop body is identical to `controls.end_to_end_ap`.
  - I did not re-run dev_report to compare numbers, only the code.
- **Tests.**
  - `tests/test_final_test_guard.py`: 69 passed.
  - `--dry-run-toy`: completes.
- **Caveat (not blocking).** The real yeast test table's `biogrid` column cannot be confirmed without reading sealed data. If it is
  absent, only the biogrid tier fails, and that failure is recorded.
