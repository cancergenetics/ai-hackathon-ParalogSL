# Round 1 review — Statistics / ML evaluation (pair-input biological prediction)

Reviewer role: statistician, ML evaluation. Files read: PROMPT.md, DESIGN.md v1, LITERATURE.md, docs/lit/human_labels.md,
docs/lit/yeast_and_other_species_labels.md, docs/lit/methods.md (evaluation sections), config/seeds.yaml.
A small power simulation (binormal scores, family clusters, paired cluster bootstrap; ~2 CPU-min) was run to put numbers
on §2. Script: `<scratchpad>`
(scratch, not part of the project; the PI may copy it into `scripts/eval/` if useful).

## 1. Overall verdict

The evaluation skeleton is sound: splits grouped by homology cluster across both species, a frozen and logged single-access
test set, paired cluster bootstrap, prevalence and lift reported, leakage baselines, and end-to-end scoring. The weak points
are quantitative, not conceptual. (a) The pre-registered success criterion is badly underpowered. Under realistic effect
sizes, the chance that both the human and the yeast CI exclude 0 is about 5–45%. The Arabidopsis clause is a coin flip under
the null. (b) The cluster construction does not guarantee its own premise. (c) The human and yeast label rules depend on
how often a pair was measured, which is confounded with library-design criteria that correlate with the features.
These must be fixed before labels and splits are frozen at M1.

## 2. Concrete changes

### Expected test-set sizes (inputs to everything below; to be recomputed by script at M1)

| Set | Labelled pairs (est.) | Positives (est.) | Prevalence | In frozen test (20% of clusters, large clusters excluded) |
|---|---|---|---|---|
| Human, primary (≥2 lines) | 3–5k confident | 300–600 | 5–10% | **≈50–120 pos**, 600–1,200 neg; perhaps 30–80 distinct positive-containing clusters |
| Yeast, primary (ε < −0.35 ∪ Dean lethal ∪ BioGRID-LT) | 0.8–1.8k | ≈60–150 (ε ≤ −0.35 alone gives ≈4.8% prevalence per Dennler; Dean adds ≤48, overlapping) | ≈5–10% | **≈12–30 pos** |
| Yeast, lenient (ε < −0.12) | same | ≈200–500 | 20–35% | ≈40–100 pos |
| Arabidopsis RD4 strict / RD9 broad | 161 / 300 | 51 / 190 | 32% / 63% | all (held out) |
| C. elegans (Tischler) | 143 | 16 | 11% | all |
| S. pombe | ≈50–200 [U] | unknown, probably <30 | — | all |

Each yeast development fold (5-fold) will hold only about 10–25 positives.

### Power of the current success criterion (simulation; identity AUROC 0.65 human / 0.70 yeast; model–identity score correlation ρ within class)

| Scenario (test set) | ΔAUROC | true ΔAP | SD(ΔAP) | P(95% CI > 0) |
|---|---|---|---|---|
| Human, 70 pos, prev 7%, ρ = 0.7 | +0.05 | +0.029 | 0.023 | 0.35 |
| Human, same, ρ = 0.9 | +0.05 | +0.031 | 0.017 | 0.65 |
| Human, 70 pos, +0.10 | +0.10 | +0.070 | 0.028 | 0.82 |
| Human, 42 pos | +0.05 | +0.030 | 0.031 | 0.22 |
| Yeast, 21 pos, prev 6%, ρ = 0.7 | +0.05 | +0.040 | 0.049 | 0.20 |
| Yeast, same, ρ = 0.9 | +0.05 | +0.041 | 0.037 | 0.31 |
| Yeast, 31 pos (30% test fraction), ρ = 0.9 | +0.05 | +0.037 | 0.025 | 0.43 |
| Yeast, 21 pos, +0.10 | +0.10 | +0.094 | 0.059 | 0.56 |
| Yeast, 52 pos (prev 15%) | +0.05 | +0.055 | 0.044 | 0.28 |
| Arabidopsis RD4, 51 pos | +0.05 / +0.10 | +0.055 / +0.120 | 0.042 | 0.20 / 0.68 |
| Arabidopsis, **no true effect** | 0 | 0 | 0.041 | P(point Δ > 0) = **0.50** |
| C. elegans, 16 pos | +0.10 | +0.115 | 0.073 | 0.36 |

Joint power for "CI > 0 in human **and** in yeast" is the product: about 0.07 at ΔAUROC +0.05 (ρ 0.7), about 0.20 at ρ 0.9,
and about 0.46 at ΔAUROC +0.10. Dennler & Ryan found identity hard to beat on screened human pairs, so a gain of about +0.05
AUROC is the realistic case. Two further findings from the simulation:
- The cluster bootstrap SD matched the true sampling SD here. Its calibration will be worse with heavy-tailed real cluster sizes.
- AP has a small-sample upward bias: with about 21 positives, mean sample AP was 0.172 against a population AP of 0.142, i.e. +20%.

### Must

1. **Re-specify the confirmatory test so that it has usable power, and compute the minimum detectable effect at M1.**
   After labels and clusters exist, and before any feature is computed, rerun the power simulation on the actual test
   composition. Record the minimum detectable ΔAP per species in RESEARCH_TREE.md. Then pre-register:
   - **Primary:** species-averaged ΔAP (equal species weights; normalised AP, see item 3), model − B1. Use a 95% paired
     cluster-bootstrap CI, stratified by species, that must exclude 0. In addition, the point estimate must be > 0 in each of
     human and yeast.
   - **Secondary:** per-species 95% CIs, reported as such.
   *Why:* the current conjunctive test (CI > 0 in both species) has joint power of about 0.1–0.5. A failure would then be
   uninformative, and that creates pressure to re-access the test set.

2. **Enlarge the yeast test set.** Either allocate about 30–35% of the clusters that contain yeast positives to the test set
   (the human fraction stays at 20%), or reconsider the ε < −0.35 primary positive threshold. The survey recommended
   ε < −0.12 for primary and the stricter threshold as a sensitivity analysis.
   *Why:* yeast is the variance bottleneck, with about 12–30 test positives. The test fraction is cheap to change before freeze.
   Choose by the rule in item 1, blind to features.

3. **Make the "≥20% relative gain" clause unambiguous and prevalence-robust.** Pre-specify:
   - the formula: ratio of species-mean AP versus mean of per-species ratios;
   - the scale: raw AP versus AP above chance, (AP − prev)/(1 − prev);
   - that the clause is judged on point estimates, with its bootstrap CI reported.
   *Why:* the same result can show a 20% raw gain and a 50% gain above chance. AP is biased upward at small n_pos. Leaving
   this open is a forking path.

4. **Replace "Arabidopsis Δ > 0 point estimate" with a reported ΔAP and 95% CI.** State in advance that an unseen-species
   generalisation claim requires the CI to exclude 0. Otherwise, report the result as "consistent with / not distinguishable
   from identity".
   *Why:* under no true effect the point estimate is > 0 with probability 0.50. Power is only about 0.2 at ΔAUROC +0.05.

5. **Fix the cluster construction, which is not "both genes in the same cluster by construction".**
   - The caller has no coverage floor (§3), while clustering requires coverage ≥ 0.3 (§6.1). A called low-coverage pair can
     therefore straddle two clusters. Run union–find over the clustering edges **plus every labelled-pair edge** (and every
     called-pair edge) so that pairs never span clusters.
   - The leakage audit must count genes, clusters and **any MMseqs2 hit between a test gene and a train gene at the
     clustering threshold**. All three counts must be 0.
   - Report the cluster size distribution (largest component as a share of labelled pairs and of positives) **before**
     deciding the >5% rule. Single-linkage at E ≤ 1e-3 tends to chain multidomain proteins into a giant component (kinases,
     small GTPases, WD40, and others).
   *Why:* leakage-free splitting is the premise of every number.

6. **Quantify what the ">5% cluster → development" rule removes from the test set.** Report the share of pairs and positives
   excluded from test eligibility, and test-versus-development distributions of identity, family size, closest-pair flag and
   source study. For every test pair, report its maximum identity to any training gene, and give AP stratified by that value.
   *Why:* the rule systematically drops large families from the test set. That is the regime where family-context features
   should help and identity is weakest, so it biases model − B1 in an unknown direction.

7. **Decouple human positive/negative calls from testing intensity.** DESIGN dropped the survey's penetrance condition. Under
   "SL in ≥ 2 lines" with 2–27 lines tested per pair, positives concentrate in heavily screened libraries. Those libraries
   were selected on criteria that correlate with our features: Harle, ≥45% identity plus a lethal fly/worm ortholog
   (conservation); Klingbeil, domain-focused signalling and epigenetic genes; Ito, 11 lines.
   - Restore the survey rule: ≥ 2 lines **and** (≥ 2 studies or penetrance ≥ 0.25).
   - Store n_lines_tested and source study per pair.
   - Add a control: predict "source library" and n_lines_tested from the sequence features (grouped CV). If this is
     predictable, report primary AP within n_lines_tested strata as well.
   *Why:* otherwise conservation and identity features can earn credit for ascertainment rather than biology.

8. **Make the identity baseline (B1) provably fair.**
   - Pre-register the finite candidate set: local fident, either direction, min/max/mean; global (Needleman–Wunsch or
     MMseqs2 `--alignment-mode 3`) identity; identity × coverage; normalised bitscore.
   - Choose the variant on development folds with the same CV as the model.
   - Evaluate B1 on exactly the same pairs as the model.
   - For the end-to-end analysis (§6.4), compute B1 by direct pairwise alignment of every labelled pair, independent of the
     caller.
   *Why:* local identity with no coverage floor is a degraded "identity". Beating a weak variant would not support "clearly
   beats sequence identity". B2 must also be reported with a CI next to B1 (it already is planned; keep it).

9. **Fully specify the bootstrap.**
   - Paired, cluster-level resampling, stratified by species and by whether a cluster contains a positive. Stratification
     keeps the prevalence roughly stable and avoids resamples with no positives.
   - 10,000 resamples for confirmatory CIs (1,000 for exploratory). Use percentile CIs, with BCa as a check.
   - Log and drop degenerate resamples.
   - Report the number of distinct positive clusters and the largest cluster's share of positives.
   - Add a leave-one-cluster-out ΔAP influence plot, showing that no single family carries the result.
   *Why:* with 12–120 clustered positives, an unspecified bootstrap is a forking path, and percentile tails at B = 1,000
   are noisy.

10. **Constrain model development to limit selection bias, and log it.**
    - Every evaluated configuration (feature set × model × hyper-parameters × caller thresholds) is appended to
      `results/dev_runs.tsv`. The count is reported in REPORT.md.
    - Selection uses pooled out-of-fold AP over the 5 repeats, compared with a paired cluster bootstrap.
    - Apply a one-SE rule: choose the simplest configuration within 1 SE of the best.
    - Freeze exactly one final model, plus B1/B2, before test access.
    - Ablation results are labelled exploratory, with CIs, no significance claims and no "Tier C helps" claims, unless the
      paired CI excludes 0 after Holm correction across the pre-listed ablation families.
    *Why:* with about 10–25 yeast positives per fold, the CV noise floor in ΔAP is about 0.04–0.05. Picking the best of many
    noisy ablations reliably overfits the development data.

### Should

11. **Redesign LOSO as cross-species grouped CV.** Train on species X using clusters in folds ≠ k, and test on species Y
    using clusters in fold k. Do not "cluster-filter" Y down to families with no X homolog: most yeast paralog families have
    human homologs, so the filtered set would be small and atypical. Report both directions.
    - Yeast → human trains on only about 50–120 positives.
    - Human → yeast is evaluated on only about 60–150 positives.
    *Why:* this keeps homology separation and uses all pairs. Interpret the asymmetric power accordingly.

12. **Use StratifiedGroupKFold** (strata = species × has-positive; groups = cluster) for the 5 × 5 CV.
    - CV uncertainty comes from a cluster bootstrap of pooled out-of-fold predictions, not from the SD across folds, which
      is anti-conservative with overlapping training sets. Use the Nadeau–Bengio correction if fold-level tests are used.
    - Rank-normalise scores within each fold before pooling, or check that pooled AP and mean per-fold AP agree.
    *Why:* unbalanced folds, and miscalibration across fold models.

13. **Build the test split deterministically, with one fixed seed and stated acceptance criteria.**
    - Acceptance criteria: per-species positive share 20–35%; identity tertiles and source studies represented.
    - At most N redraws, judged only on label counts. Log every draw.
    *Why:* re-drawing until the split "looks balanced" is a hidden selection step.

14. **Stratify the unseen-species results by homology to the training families.** For each test pair, record whether either
    gene hits a training-labelled cluster at the clustering threshold, and report AP for both strata.
    *Why:* "unseen species" is not the same as "unseen family". S. pombe especially is close to S. cerevisiae.

15. **Correct the genome-wide-prevalence reweighting.** Scaling the screened negatives assumes they represent all paralog
    pairs; they do not (median identity 54% vs 32%; 89% vs 13% closest pairs).
    - Use covariate-shift weights: bins of identity × family size × closest flag, matched to the tool's full de novo pair set.
    - Show a sensitivity band for prevalence 0.2–0.5%.
    - Label the result model-based.

16. **Define the B3 gene-prior baseline where it is computable.** Use a random pair split, and a gene-disjoint split that
    shares families. Under the cluster split, every test gene is unseen, so B3 equals the training mean.
    *Why:* B3 measures leakage only in splits where leakage is possible.

17. **Weight large families in training and as an evaluation sensitivity.** Family pair counts grow as O(n²).
    - Training: cap the total sample weight per cluster.
    - Evaluation: report AP with cluster-capped weights alongside unweighted AP.
    *Why:* a few big families should not dominate the fit or the AP.

18. **Remove the yeast label dependence on the number of measurements.** "ε < −0.35, p < 0.05 in **any** orientation or
    allele" (positive) and "**never** ε < −0.08, p < 0.05" (negative) both depend on how many measurements a pair has. TS
    alleles of essential genes add measurements, and essentiality is a forbidden signal that correlates with SL.
    - Aggregate to one estimate per pair from deletion alleles: an inverse-variance mean of ε across orientations, with a
      combined p-value.
    - Store n_measurements, and show the label rate as a function of n_measurements.
    - Keep pairs with an essential member as a separate reported stratum, as the survey suggested.

19. **Report a label-noise ceiling.** Where screens overlap, compute the AP that one study's calls achieve in predicting
    another study's calls on shared pairs and lines. Report the lenient-tier AP alongside.
    *Why:* replicate precision/recall is about 0.40/0.32. Label noise caps the achievable AP and shrinks ΔAP, which is part
    of the power problem. Readers need the ceiling to judge the numbers.

20. **Pre-specify the estimator.**
    - Step-wise average precision (`sklearn.metrics.average_precision_score`), never trapezoidal AUPRC, which is optimistic.
    - Tie handling: tied blocks are treated as one threshold. This matters for missed pairs set to the minimum in §6.4, and
      for identity ties.
    - Hyper-parameter selection uses the mean of normalised AP across species, not raw AP. Raw AP is prevalence-dependent.

### Could

21. Report a permutation-null distribution of AP for each test set (labels shuffled across clusters). This gives the
    small-sample-biased chance level, which is above the prevalence. Report the observed lift against it.
22. As a second view of the primary comparison, run a cluster-level paired permutation test (swap the two models' within-cluster ranks).
23. Add a secondary, higher-powered estimate: repeated grouped CV over development + test with the frozen pipeline. Label it
    clearly as partly selection-contaminated. Use it only to narrow uncertainty, never to replace the test result.
24. Add a study-held-out secondary test (Harle + Chan) to measure robustness to platform bias, as in human_labels.md §5.
25. Calibrate the output `sl_score` (it is class- and species-weighted, so it is not a probability), or document it as a
    rank score.

## 3. Risks the design does not address

- **Abuse of the "re-run on bug" clause.** A low-powered test that "fails" invites a second access. Define a bug narrowly:
  a failure of a verification test that existed before the first access, or a mismatch between documented and executed code.
  Pre-commit the analysis script's hash in the access log. Any second run reports both results.
- **Label rules tuned after seeing features.** The 13-pair check and the Dennler agreement check are fine. Any threshold
  change after Tier A features exist must be logged as a decision note, with the reason, and rerun blind to model scores.
- **Different SL constructs across test species.** Arabidopsis uses developmental lethality, C. elegans an RNAi synthetic
  phenotype, yeast colony-size ε, human CRISPR in cancer lines. AP values are not comparable across species, and pooled
  "unseen-species AP" should not be reported.
- **Heterogeneous meaning of human labels.** Per-study thresholds produce different false-positive rates, and study is
  confounded with library design (item 7).
- **The caller's recall interacts with the comparison.** If recall of labelled pairs is low or depends on identity, the
  end-to-end and called-pairs-only analyses can disagree in sign. Pre-specify that called-pairs-only is the primary
  comparison (same pairs for model and B1) and that end-to-end is co-reported.
- **PLM memorisation of families.** Pretraining on UniRef is not label leakage. But per-gene embedding features (Tier C8)
  can identify families, and only the homology split protects against that. Item 6 (nearest-training-identity
  stratification) is the check.
- **Small n per stratum.** Identity-stratified analyses of Tier C/D with about 3–5 bins will have only a few yeast positives
  per bin. Pool development folds and report CIs, or the "conditional on identity" check will be noise.

## 4. Fine as it is

- Splitting by homology cluster jointly over human and yeast proteomes (subject to item 5), with test families never seen in
  training in either species.
- A single frozen test split with an access log, and one planned access; refitting on all data only after the test result
  is recorded.
- Per-species AUPRC with prevalence and lift, AUROC as secondary, the lenient tier as sensitivity, and no pooled
  cross-species AP.
- Ambiguous pairs excluded from the primary metric; tested negatives only, with no "unlabelled = negative".
- Paired bootstrap for Δ (model − B1), and the B2 trivial-model and B3 gene-prior baselines (subject to item 16).
- The random-split leakage-inflation estimate, the shuffled-label control and the leakage-audit script.
- Unseen species never used for training, selection or thresholds.
- Fixed seed list; numbers produced by scripts into `results/*.json`.
