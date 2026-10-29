# Round 2 review: statistics and ML evaluation

Reviewer role: statistician (ML evaluation). Files read: REVIEW_BRIEF.md, PROMPT.md, DESIGN.md v2,
reviews/round1/statistics.md, reviews/round1/PI_response.md, RESEARCH_TREE.md (D2), config/seeds.yaml, and the yeast label
survey (DMF section).

I reran the round-1 power simulation on the v2 decision rule. It is small: binormal scores, family clusters, a stratified
paired cluster bootstrap, 120–200 simulated datasets and B = 300.
- Script (scratch, not part of the project):
  `<scratchpad>`
- Assumed composition: 30% test. Human ≈ 1,300 pairs with ≈ 91 positives. Yeast ≈ 600 pairs with ≈ 36 positives (and a
  variant with 25). B1 AUROC 0.65 for human and 0.70 for yeast.

| Scenario | Mean T (species-mean ΔAP) | SD(T) | P(CI of T > 0) | P(full v2 rule) | P(rule and relative gain ≥ 0.20) |
|---|---|---|---|---|---|
| No true effect, ρ = 0.7 | −0.002 | 0.016 | 0.01 | 0.00 | 0.00 |
| ΔAUROC +0.05, ρ = 0.7 | +0.033 | 0.019 | 0.37 | 0.37 | 0.37 |
| ΔAUROC +0.05, ρ = 0.9 | +0.034 | 0.013 | 0.75 | 0.74 | 0.62 |
| ΔAUROC +0.10, ρ = 0.7 | +0.082 | 0.025 | 0.95 | 0.95 | 0.95 |
| ΔAUROC +0.05, ρ = 0.7, yeast 25 positives | +0.033 | 0.024 | 0.32 | 0.32 | 0.32 |

What the simulation shows:
- **Power.** The 80%-power MDE of the v2 rule is about 0.035–0.055 AP (2.8 × SD(T)).
- **The "MDE > 0.10" fallback trigger would never fire** at any plausible composition. The test can therefore stay
  underpowered for the realistic effect (+0.05 AUROC ≈ +0.033 AP; power 0.37 at ρ = 0.7) and the fallback still never
  activates.
- **The per-species "ΔAP_s > 0" conjunct costs almost no power.** That is good.
- **The stratified bootstrap is, if anything, slightly conservative.** Under the null, 1% of CIs excluded 0 against a
  nominal 2.5%, with Monte-Carlo SE ≈ 0.01.

## 1. Overall verdict

v2 fixes most of my round-1 must-items in substance:
- a species-mean primary criterion with a stratified paired bootstrap;
- a larger test set;
- a reported CI for Arabidopsis instead of a coin-flip clause;
- a pre-registered set of B1 variants;
- an experiment registry with the one-SE rule;
- the human penetrance rule;
- measured leakage.

Four new elements are under-specified in ways that reopen forking paths:
- the MDE trigger is set where it can never fire, and its simulation inputs are free;
- the "co-primary nested CV" has no decision rule or α allocation, and is partly selection-contaminated;
- Leiden partitioning can cut labelled-pair edges, and its resolution and cap are unset;
- the panel-masking augmentation is statistically harmless only if the 4 copies of a pair are kept together everywhere
  (folds, early stopping, bootstrap, leaf-size constraints) and evaluation uses one row per pair.

None of these needs a redesign. Each is a few lines in DESIGN.md or `config/`, and all must be fixed before M1 freezes
the splits.

### Status of round-1 must-items

| R1 item | Status in v2 | Remaining gap |
|---|---|---|
| 1 Powered primary + MDE | Mostly fixed | Trigger level and simulation inputs (items 1–2 below) |
| 2 Larger yeast test | Fixed (~30%, ≥30 positives target) | The yeast count under the new DMF rule is unverified (item 9) |
| 3 Relative-gain clause | Partly fixed (mean of per-species ratios, raw AP) | Say "point estimate"; report the above-chance scale (item 12) |
| 4 Arabidopsis | Fixed | — |
| 5 Cluster construction | Partly fixed | Leiden can cut labelled-pair edges; resolution and cap are free (item 3) |
| 6 What the large-group rule removes | Partly fixed | Leakage dose–response instead of only a no-hit subset (item 4) |
| 7 Human ascertainment | Fixed (penetrance rule, source-library control) | Stratify by n_lines if the control is positive (item 13) |
| 8 Fair B1 | Mostly fixed | The primary pair population and B1's scores on caller-missed pairs are unstated (item 5) |
| 9 Bootstrap specification | Mostly fixed | "Stratified by label" is undefined at group level; mixed-species groups (item 6) |
| 10 Selection discipline | Fixed | The Holm/exploratory labelling of ablations is not stated (item 14) |

## 2. Concrete changes

### Must

1. **Lower the MDE trigger, and apply it to the whole decision rule.**
   - Replace "if the 80%-power MDE of the primary statistic exceeds 0.10 AP" with "if the 80% power of the full primary rule
     (CI of T excludes 0 **and** ΔAP_s > 0 in both species) at **T = 0.04 AP** is below 0.8".
   - 0.04 AP ≈ ΔAUROC +0.05–0.06 here. That is the realistic effect, given that Dennler & Ryan found identity hard to beat.
     It is also roughly a 20% relative gain at the expected B1 AP.

   *Why:* the simulation gives MDE ≈ 0.035–0.055 at the planned composition, so a 0.10 trigger is inert. The fallback would
   never activate in the regime it was written for.

2. **Pre-register the inputs of the MDE simulation** in `config/mde.yaml`, committed and hashed before it runs.
   - **Score model:** binormal.
   - **B1 AUROC:** human 0.65, yeast 0.70, from the literature, not from our data.
   - **Within-class model–B1 correlation:** ρ = 0.8 is decisive; 0.7 and 0.9 are reported.
   - **Group structure:** cluster-level random effect, with the real group IDs and group sizes of the frozen test pairs.
   - **Simulation size:** B ≥ 1,000 and ≥ 500 simulated datasets.
   - **No scores of any kind are computed on test pairs.** That includes identity, which is a feature.

   *Why:* ρ alone moves power from 0.37 to 0.75. Leaving it free lets the trigger outcome be chosen. Computing B1 on test
   pairs to "calibrate" the simulation would be a test access.

3. **Make labelled-pair edges uncuttable, and fix Leiden's free parameters.**
   1. Contract every connected component of the *labelled-pair* graph (both species) into a supernode before Leiden.
      Otherwise Leiden can place gene A and gene B of a labelled pair in different groups. The pair then belongs to no
      split, or to a split with only one of its genes.
   2. Pre-register the cap and the resolution rule:
      - cap = max(5% of a species' labelled positives, 5% of its labelled pairs) per group;
      - the resolution γ is the smallest value on a fixed grid at which every group is ≤ cap;
      - the seed is the first seed in `config/seeds.yaml`;
      - edge weight = bitscore / min(self-bitscore), not raw bitscore. Raw bitscore scales with length, so long
        multidomain proteins would dominate the cut.
   3. Report the cut-weight fraction, the number of groups, and the size distribution.

   *Why:* round-1 item 5's premise (no pair spans groups) is again not guaranteed. γ and the cap are tunable knobs that
   change test composition.

4. **Report leakage as a dose–response, not only as a strict no-hit subset.**
   - Define leakage(pair) = max over its two genes of the bitscore ratio (bitscore / self-bitscore of the test gene) to any
     gene in a *development labelled pair* of either species.
   - Report T and per-species ΔAP within leakage tertiles, with a bootstrap CI for the trend (highest minus lowest tertile).
   - Keep the no-hit subset, but describe its composition (size, identity, family size, species), and do not read its ΔAP
     as "the leakage-free effect". It is a different population: isolated, mostly small families.

   *Why:* B1 is not trained, so homology leakage can only inflate the model side of ΔAP. A trend of ΔAP rising with
   leakage is the direct signature of memorisation.

5. **State the primary pair population and how B1 scores caller-missed pairs.**
   - Primary population = test labelled pairs called by the caller. The model and B1 are scored on exactly these pairs.
   - The end-to-end analysis (all labelled test pairs, missed pairs at the minimum model score, ties treated as one block)
     is co-reported. In it, B1 comes from a direct pairwise alignment independent of the caller, with no alignment = 0.
   - Report caller recall per species and tier on test pairs next to T.

   *Why:* if the two populations are not fixed in advance, one could pick whichever gives the better ΔAP after seeing the
   result. Round 1 raised this risk and v2 does not resolve it.

6. **Make the bootstrap well defined on joint human–yeast groups.**
   - Groups now contain pairs of both species, so "stratified by species and label" is undefined at the resampling unit.
   - Resampling unit: the Leiden group (supernode-respecting).
   - Strata: {human-only, yeast-only, mixed} × {contains ≥ 1 positive, none}.
   - In each resample, T = mean of per-species ΔAP computed from the pairs drawn.
   - Use 10,000 resamples with percentile CIs and BCa as a check. Log and drop resamples with no positives in a species,
     and report how many were dropped.

   *Why:* this is what the simulation shows to be valid. The current wording admits several different procedures.

7. **Keep panel-masking copies together everywhere, and evaluate on one row per pair.** Duplicating each pair under k ∈
   {0,1,2,3} masks with weight 1/4 is statistically neutral for estimation **only if**:
   1. **Fold and split assignment** is by group ID, never by row. All 4 copies share a fold, including in LOSO and in the
      nested CV.
   2. **Class and species weights** are computed on pair counts and then multiplied by 1/4, so every pair still sums to the
      intended weight.
   3. **Row-count hyper-parameters** are rescaled, or replaced by Hessian-based ones. LightGBM's `min_data_in_leaf` and
      HistGB's `min_samples_leaf` count rows, so a leaf of 20 rows may hold 5 pairs. Use `min_sum_hessian_in_leaf` (HistGB
      has no direct equivalent), or multiply the count by 4. The same applies to row bagging (`bagging_fraction`).
   4. **Internal early-stopping validation** is a group-aware split. A random validation fraction leaks copies across
      train and validation, which ends training late and optimistically.
   5. **All reported metrics, SEs for the one-SE rule, and test, LOSO and unseen-species scores** use **k = 0 rows only**,
      one per pair. These are the features the deployed tool computes. Masked rows enter evaluation only in the
      panel-thinning test, one mask at a time.

   Under those conditions the augmentation does not distort CV or the bootstrap. The bootstrap resamples test predictions,
   not training rows, and it never sees copies. Add a unit test asserting that every pair ID has its copies in a single
   fold and that evaluation frames contain unique pair IDs.

8. **Give the co-primary nested CV a decision rule, an α split and a contamination label.** If item 1 triggers:
   - **α split.** Each estimate gets α = 0.025 (97.5% CIs). Pass if either CI excludes 0 **and** the other estimate's point
     value is > 0, plus ΔAP_s > 0 in both species for the estimate that passes. Any other "either one" rule inflates
     type-I error.
   - **Order of operations.** The nested CV necessarily reads test labels, so it runs only **after** the single sealed test
     access is logged, as a second logged access of the same frozen config.
   - **Inner loop.** It re-runs every *automated* choice in each outer fold: grid + one-SE rule, B1 variant, caller
     thresholds, and the augmentation on/off decision if that was made by a rule. Manual feature-set choices stay fixed at
     the frozen config and are disclosed as contaminating.
   - **Reporting.** Report the nested-CV ΔAP separately for held-out groups that came from the original test split.
     Those are uncontaminated by development. Report the remaining groups separately as well.

   *Why:* as written, "supplemented by nested CV as a co-primary estimate" lets the more favourable of two estimates
   decide. Agent-driven feature engineering on the development set also biases a CV over the development data upward.

### Should

9. **Recompute the expected yeast positive count under the new DMF rule before relying on "≥30 test positives".**
   - With healthy single mutants (SMF ≈ 0.9–1.0), DMF ≤ 0.5 already requires ε ≲ −0.3 to −0.4. So the v2 primary tier is
     close to the round-1 strict tier, plus pairs with sick single mutants. The PI's claim that it "raises yeast positive
     counts" is unverified.
   - Report the count at M1. If yeast test positives fall below 30 at a 30% test share, apply the pre-stated rule (raise
     the yeast test share up to 40%) before any feature exists.
   - Store min(SMF) per pair and report AP within SMF strata. The DMF rule makes positivity depend on single-mutant
     fitness, which is an ascertainment axis like n_lines_tested in human. SMF is fine as an evaluation stratifier and
     forbidden as a feature.

10. **Make the yeast sensitivity tiers nested.** The strict tier "DMF ≤ 0.2 **or** ε < −0.35" is not a subset of the
    primary tier:
    - a pair with ε < −0.35 and DMF 0.6 is strict-positive but primary-ambiguous;
    - a pair with DMF ≤ 0.2 and no significant ε is strict-positive with no interaction at all.

    Define strict = primary ∧ (DMF ≤ 0.2 ∨ ε < −0.35). Lenient should be a superset of primary, or be documented as an
    alternative definition rather than a looser one.

    *Why:* a sensitivity analysis is interpretable only when moving along the tiers changes stringency alone.

11. **Fix the reporting hierarchy for tiers, strata and secondary comparisons in `config/report_plan.yaml`** before
    features exist.
    - **Level 1:** T and the rule (one confirmatory test).
    - **Level 2:** pre-listed secondary estimates (Δ vs B2, LOSO both directions, unseen species, end-to-end, tiers). Each
      is shown with a 95% CI and no pass/fail.
    - **Level 3:** strata (WGD/SSD, ribosomal, identity tertiles, leakage tertiles, SMF, n_lines, source study), shown as
      a forest plot.
    - Every listed estimate is reported whatever its result. A stratum is reported with a CI only if it holds ≥ 10
      positives; otherwise counts are given.
    - Statements that an effect differs between strata need a CI on the *difference* (interaction), not two
      non-overlapping or overlapping CIs.
    - Nothing is added to the plan after the test access without an "added post hoc" label.

    *Why:* v2 now has about 25 tier × stratum × comparator cells on a test set with 25–40 yeast positives. Without a fixed
    plan, selective emphasis is the main residual risk.

12. **Pin the relative-gain clause completely.** It is judged on the point estimate of mean_s(AP_model,s / AP_B1,s − 1),
    with a bootstrap CI reported. Also report the same quantity on the above-chance scale, (AP − prev)/(1 − prev).
    *Why:* in the simulation the 20% clause is binding at ρ = 0.9 (0.74 → 0.62). Whether it is judged on the point
    estimate or the CI has to be fixed in advance.

13. **Human ascertainment follow-through.** If the source-library or n_lines_tested control predicts above chance
    (grouped-CV AUROC CI excludes 0.5), report T within n_lines_tested strata (2–3, 4–10, > 10) as a pre-registered
    sensitivity analysis. Otherwise state that the control was negative.

14. **Restore the remaining round-1 should-items that v2 does not mention:**
    - ablations are exploratory, and a Holm correction applies before any "tier X helps" claim (R1 10);
    - deterministic test-split construction with logged redraws and acceptance criteria on label counts only (R1 13);
    - for LOSO with joint groups, train on species X in folds ≠ k and test on species Y in fold k. Testing on all Y would
      reuse orthologous groups (R1 11);
    - per-pair leakage (item 4) also computed for unseen-species pairs, with AP by leakage tertile (R1 14);
    - covariate-shift weights for the genome-wide-prevalence precision (R1 15);
    - a label-noise ceiling from overlapping screens (R1 19).

15. **Define the one-SE rule's SE.** Use the cluster-bootstrap SD of pooled out-of-fold T over k = 0 rows, with OOF scores
    averaged over the 5 seed repeats. Do not use the SD across folds, which is anti-conservative.

16. **Pre-declare the adoption rule for the masking augmentation itself.** Compare models with and without augmentation on
    development CV. Adopt augmentation if:
    - T at k = 0 is not worse by more than 1 SE; **and**
    - mean T over masks k = 2, 3 improves by more than 1 SE.

    The panel-thinning test then reports both models under every mask.

    *Why:* the 1/4 weighting puts 75% of the training weight on masked distributions while human and yeast are evaluated
    at k = 0. Whether that trade is worth it is an empirical question and should be settled by a fixed rule.

### Could

17. Bootstrap sensitivity with a coarser resampling unit (Leiden at a lower resolution, e.g. γ/4). Groups joined by cut
    homology edges are not independent, so group-level resampling can slightly under-state variance.
18. Report the permutation-null AP per test set (R1 21). With 25–40 yeast positives, the small-sample chance level is above
    prevalence.
19. Report how sensitive the LightGBM fit is to duplicating copies: fit once with the 4 copies and once with the k = 0 rows
    only at weight 1, and confirm that OOF T at k = 0 differs by less than 1 SE when masks are uninformative (for example
    on a panel with no close bins).

## 3. Risks the design does not address

- **The test can stay underpowered for the realistic effect even with every fix.** At ρ = 0.7 and ΔAUROC +0.05, power is
  about 0.37. A non-significant result must be reported as "not distinguishable at this n" with the CI, not as "does not
  beat identity". The MDE should be printed next to T in REPORT.md.
- **Leiden residual homology is concentrated in the giant component.** That is where large families sit (kinases, GTPases),
  so leakage and family size are confounded. Item 4's tertile analysis should be read together with family-size strata.
- **The DMF rule ties yeast positivity to single-mutant fitness.** Conservation and family features correlate with that
  fitness, which reproduces the essentiality signal the features are forbidden to use (item 9 stratifier).
- **Masked-training mismatch.** If the augmentation helps unseen species but costs k = 0 performance, the primary T is
  computed on a model tuned for a different target. Item 16's rule decides this explicitly.
- **Multiplicity across many reported cells.** About 25 pre-listed cells mean that 1–2 "significant" secondary CIs are
  expected under the null. Item 11's hierarchy is the protection.

## 4. Fine as it is

- T = species-mean ΔAP with ΔAP_s > 0 in both species: valid, nearly free in power, and of correct size in simulation.
- 10,000-resample paired, stratified cluster bootstrap and the leave-one-group-out influence check.
- ~30% test share with a minimum on yeast positives, MDE computed at M1 from the real composition before any feature
  exists (once items 1–2 are applied).
- The Arabidopsis and other unseen-species results reported as ΔAP with a CI, pre-registered as underpowered.
- The B1 variant set chosen on development folds and scored on the same pairs.
- B2 as a pre-registered secondary comparator.
- The experiment registry, the one-SE rule and a single frozen model.
- The sealed test directory, the tag and the access log, and the re-run rule that reports both numbers.
- The human robust tier, with ≥ 2 lines and (≥ 2 studies or penetrance ≥ 0.25), plus a borderline margin for negatives.
- Training features produced by the shipped CLI, which removes train/inference skew.
- Pre-registering the panel-thinning test. Weighting each pair at 1/4 per mask, so the total weight per pair is 1, is the
  right construction once item 7's conditions hold.
