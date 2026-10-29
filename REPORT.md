# REPORT — sequence-only prediction of synthetic-lethal paralog pairs

Final report, 2026-10-28. Goal: `PROMPT.md`. Every number below is copied from the committed file cited next to it; the commit is given where it helps.

## 1. Summary

- `slparalog` is an installable CLI. It takes one protein FASTA, calls all paralog pairs with MMseqs2 and scores each pair with a gradient-boosted tree model. The model uses 73 sequence-only features: 50 from a self-search of the proteome and 23 comparative features from an 80-proteome reference panel.
- On the held-out test set (single access, `results/final/final_test.md`, a88cff1) the tool **clearly beats sequence identity**: T = 0.297 [0.131, 0.395].
  - Human AP is 0.349 against 0.106 for identity (prevalence 0.029).
  - Yeast AP is 0.541 against 0.190 (prevalence 0.054).
  - The model also beats B2 (identity + family size + closest-pair flag) in both species.
- **Caveat 1, ascertainment.** Most of the human headline ranks pairs *between* screening-ascertainment cells. Within identical study-set cells the gain over identity is −0.020 [−0.098, 0.060], as the adversarial dev review predicted (D17). A model trained on yeast alone still lifts human AP to 0.259, so part of the signal is transferable biology.
- **Caveat 2, unseen species.** The gain transfers to *C. elegans* (ΔAP +0.178 [0.064, 0.349]) but not to *Arabidopsis* (strict ΔAP −0.005 [−0.111, 0.097]), where B2 also scores higher than the model. Plants are therefore not supported.
- **Runtime on a full human proteome (20,770 proteins, 686,225 scored pairs):** 99.8 min wall time with 10 threads, 2.52 GB peak memory (process tree), from an empty cache (`results/bench/e2e_human.json`). The shipped pipeline reproduces the evaluated model's scores exactly.
- **Robustness (simulated draft proteomes, development labels).**
  - The model tolerates fragments, split genes and fusions better than identity does.
  - Missing genes cost both methods equally.
  - With all defects combined, the human margin over identity roughly halves (§6a).
- **Other caveats:**
  - The primary comparison was pre-declared underpowered.
  - The release model (refit on all labels) has no independent evaluation.

## 2. Problem and constraints

The input is one protein FASTA with one protein per gene. The output is every paralog pair the tool calls, each with an SL score. Every feature must be computed from sequences alone, because the tool must work on newly sequenced, unannotated proteomes, plants and insects included (`PROMPT.md`).

Allowed inputs:
- the query sequences;
- a fixed, versioned panel of other species' proteomes, with headers stripped;
- pretrained sequence models.

Forbidden inputs, both as features and as auxiliary training targets (D1):
- expression, PPI, essentiality and genetic interactions;
- GO, Pfam or any other annotation;
- labels of orthologous pairs;
- hand-entered knowledge of specific genes.

Annotations such as YGOB, OHNOLOGS or ribosomal keywords appear only in stratified evaluation (DESIGN §1). Supervision comes from yeast and human SL labels only. The metric is AUPRC (`sklearn.metrics.average_precision_score`). Hardware: one Apple M4 laptop, 10 cores, 24 GB, no GPU.

## 3. Labels delivered

All labels were built by scripts from primary tables, under rules committed before any label code ran (`config/labels_human.yaml`, `config/labels_yeast.yaml`, `config/labels_unseen.yaml`). No label came from SLKB or SynLethDB.

**SL definition, the same in both species (DESIGN §2.3):** a significant negative genetic interaction plus a strong absolute double-mutant defect, in a context where both single mutants are testable.

**Label universe:** tested pairs ∩ Ensembl Compara 116 paralogs. Compara defines the universe only; it is never a feature.

- **Human:** 9 combinatorial screens; per cell line, calls come from harmonised LFCs.
  - Primary positive: SL in ≥2 cell lines, plus either ≥2 studies with penetrance ≥0.10 or penetrance ≥0.25.
  - Primary negative: testable in ≥2 cell lines and SL in none.
- **Yeast:** Costanzo 2016 SGA, null alleles only.
  - Positive: ε < −0.08, p < 0.05 and DMF ≤ 0.5.
  - Negative: DMF > 0.5 with ε ≥ −0.04, or ε ≥ −0.08 with p ≥ 0.05.
  - Dean 2008 contributes under its own rule. BioGRID low-throughput SL forms a separate tier only.

| species / set | tier | positives | negatives | source file |
|---|---|---|---|---|
| human, in universe | primary | 205 | 7,026 | `results/splits/split_summary.json` |
| human, all tested pairs (11,125) | primary / lenient / ≥2-studies | 222 / 907 / 152 | 8,040 / 9,401 / 2,196 | `results/labels/human_label_qc.json` |
| yeast, in universe | primary | 104 | 2,049 | `results/labels/yeast_label_qc.json` |
| yeast, in universe | lenient / strict / BioGRID tier | 183 / 68 / 227 | 2,047 / 2,049 / 1,998 | same |
| development split | primary, human / yeast | 123 / 62 | 4,297 / 1,143 | RESEARCH_TREE D9 |
| test split | primary, human / yeast | 82 / 42 | 2,729 / 906 | D9; `results/final/final_test.json` |
| *Arabidopsis* (unseen; Cusack 2021) | strict / broad | 51 / 190 | 107 / 107 | `results/labels/unseen_label_counts.json` |
| *C. elegans* (unseen; Tischler 2006) | primary | 16 | 120 | same |
| *S. pombe* | – | unavailable (the E-MAP matrix needs a manual browser download) | | same |
| insects | – | no test set with tested negatives exists | | same |

**QC:**
- Human: 12 of 12 published counts reproduced (D8).
- Yeast: an exact independent re-implementation agrees on all 9.1 M pairs (D10).
- Both label sets passed independent review (`reviews/results/human_labels_review.md`, `yeast_labels_review.md`).
- Ascertainment is strong: the number of testable lines alone gives AUPRC 0.127 against a prevalence of 0.027 (D8).

**Files:**

| content | path |
|---|---|
| all labels | `data/labels/human_labels_all.tsv`, `data/labels/yeast_labels_all.tsv` |
| development split | `data/labels/dev/` |
| release training labels (development + test) | `data/labels/release/` (`results/release/release_model.json`) |
| test and unseen labels | `data/labels/sealed/`, with `MANIFEST.sha256` |

**Licences.** Sources and checksums are listed in `DATA.md`. Committed rows are derived binary calls with per-source evidence codes; raw tables are never committed.

| sources | licence |
|---|---|
| Human: Dede, Thompson, in4mer/Chou, Harle, Chan | CC BY 4.0 |
| Human: Parrish, Klingbeil | CC BY-NC-ND 4.0 |
| Human: Ito, CHyMErA | no open licence on the supplementary tables |
| Yeast: Costanzo 2016 (TheCellMap) | none given (treated as research use) |
| Yeast: Dean 2008 | CC BY 4.0 |
| Yeast: BioGRID | MIT |
| Unseen: Cusack 2021 | CC BY 4.0 |
| Unseen: Tischler 2006 | CC BY |

- Human rows carry `license_classes` and `only_nonopen_evidence`.
- 44 of the 222 human primary positives rest only on non-open evidence (`human_label_qc.json`).
- Before any public release of the label files, the derived calls from non-open sources need a licence check.

## 4. Approach

**Paralog caller (DESIGN §3, `config/caller.yaml`).**
- Search: MMseqs2 all-vs-all with `-s 7.5`, `--max-seqs 1000` and `-a`.
- A pair is called if the bitscore-based E′ ≤ 1e-3 in either direction. E′ uses a fixed reference size, so calling does not depend on proteome size. There is no identity or coverage floor.
- Capped queries are re-searched only if their *called* hits reach 95% of the cap (D13).
- Called pairs: human 686,225, yeast 18,014.
- Recall of Compara pairs: human 89%, yeast 82%; 98–100% above 30% identity (D11).

**Feature tiers.** All features are symmetric in A and B and listed in a registry. A test checks that no registered feature comes from a forbidden source.

| tier | content | in frozen model | status |
|---|---|---|---|
| A, self-search | identities (local, Ensembl-style, global NW); coverage; bitscore ratios; length ratio; family context (hit counts at several bitscore-ratio levels, partner rank, RBH, shared third paralogs, best other paralog, mutual-coverage family size); per-gene length and low complexity; proteome percentiles | 50 | adopted |
| B, reference panel (80 proteomes, ≈1.55 M proteins) | one MMseqs2 search per panel species; divergence calibration d_s (Kimura distance over a core set, no species tree); pre/post duplication classification and depth bounds; relative rate; conservation; post-duplication single-ortholog count; species rank; outgroup identity | 23 | adopted; `tb_asym_norm` and `tb_asym_abs` excluded (94–98% NaN; Tier B review) |
| C, protein language model | ESM-2 150M mean-pooled embedding similarity and neighbourhood | 14 built | not adopted (D17b) |
| D, predicted structure | ProstT5 3Di + Foldseek | – | not built |

**Model.**
- scikit-learn HistGradientBoosting, configuration `GBDT_AB_17`: learning rate 0.03, 100 iterations, 7 leaves, min leaf 50, L2 1.0.
- Each species gets equal total weight and equal weighted prevalence.
- Frozen model `AB-gbdt-7ff069ac`, seed 11. It was trained on all 5,411 primary development pairs (human 4,415, yeast 996) at clean commit 563b15d (D19; `results/frozen/model/model.json`).
- The score ranks pairs. It is not a calibrated probability (README).

**Baselines (DESIGN §5, D8).**

| baseline | definition |
|---|---|
| B0 | random |
| B1 | identity: the best of 7 pre-registered identity variants on development folds, which was `ident_global` (Needleman–Wunsch) |
| B2 | logistic regression on identity + family size + RBH flag |
| B3 | gene-prior leakage baseline; uninformative by construction, because the splits keep homologous genes together |
| B4 | number of testable lines; a human ascertainment diagnostic, never a tool input |

**Splits (DESIGN §6.1, D9).**
- Construction: a homology graph over all labelled genes of both species, with each labelled pair contracted to a supernode, then connected components.
- Groups holding more than 10% of a species' labelled pairs go to development. This moved one 238-gene human family (1,611 pairs, 0 primary positives).
- Test fraction 0.40, one draw.
- Power at T = 0.04 AP was 0.588 (`results/mde.json`), so the primary comparison was pre-declared underpowered.
- Development: 5-fold grouped CV × 5 seeds {11, 23, 37, 41, 59}.

**Statistics.**
- Primary statistic T = mean over {human, yeast} of ΔAP(model − B1) on test.
- Pass: the 95% CI excludes 0 and ΔAP > 0 in each species. "Clearly beats" also needs a mean relative gain ≥ 0.20.
- CIs: 10,000 group-bootstrap resamples, stratified by species composition × positive-holding, percentile intervals.
- On development, selection uses T_pf, the mean per-fold paired ΔAP vs B1 (D12a).

**Adoption rule (`config/models.yaml`; D12 as restated in 95ac397).** The one-SE-selected configuration of each family is adopted only if all three hold:
1. the T_pf CI excludes 0;
2. per-fold ΔAP > 0 in both species;
3. the human stratum-averaged AP (strata: 1, 2, ≥3 testing studies) ≥ that of B1.

**Guards.**
- `load_labels()` refuses `data/labels/sealed/`.
- Tests check protected-file hashes.
- Feature builders never see labels.
- `scripts/eval/final_test.py` checks its own preconditions (tag `final-model-v1`, clean tree, model hash and frozen spec, sealed manifest, recorded test pass, exactly one primary access) and appends to `results/final/access_log.tsv`.
- Every key result had an independent review (`reviews/results/`).

## 5. Development results

Population for every development run: human 4,415 called pairs (123 positive, prevalence 0.0279) and yeast 996 (60 positive, prevalence 0.0602). All runs and their pre-declared rules are logged in `results/experiments.tsv`.

| method (selected configuration) | human AP | yeast AP | T_pf [95% CI] | human stratum-avg AP | outcome | file |
|---|---|---|---|---|---|---|
| B1 identity | 0.118 | 0.567 | – | 0.201 | guard reference | noasym |
| B2 | 0.171 | 0.568 | 0.068 [0.036, 0.089] | 0.203 | comparator | noasym |
| B4 ascertainment (human only) | 0.153 | – | – | 0.111 | diagnostic | noasym |
| LR, Tier A | 0.179 | 0.558 | 0.058 [0.026, 0.099] | 0.180 | fails guard | noasym |
| GBDT, Tier A | 0.165 | 0.548 | 0.046 [0.011, 0.084] | 0.157 | fails guard | noasym |
| LR / GBDT, Tier A deconfounded | 0.174 / 0.158 | 0.559 / 0.568 | 0.052 / 0.041 | 0.179 / 0.168 | fails guard | dc |
| GBDT monotone in identity | 0.166 | 0.545 | 0.043 [0.008, 0.082] | 0.159 | fails guard | ablation |
| LR / GBDT rank-ensemble with B1 | 0.180 / 0.156 | 0.603 / 0.605 | 0.052 / 0.040 | 0.200 / 0.194 | fails guard | ablation |
| LR / GBDT, A+C | 0.184 / 0.193 | 0.567 / 0.543 | 0.064 / 0.059 | 0.184 / 0.185 | fails guard | AC |
| LR, A+B | 0.265 | 0.581 | 0.099 [0.062, 0.150] | 0.233 | passes | noasym |
| **GBDT, A+B** | **0.361** | **0.591** | **0.143 [0.095, 0.187]** | **0.260** | **adopted, frozen** | noasym |
| GBDT, A+B+C | 0.374 | 0.591 | 0.154 [0.103, 0.198] | 0.266 | passes; not adopted (D17b) | noasym |

File key (`results/dev/`):

| key | file | commit |
|---|---|---|
| noasym | `dev_primary_AB_ABC_noasym_6bb9438.json` | e6849e8 |
| ablation | `dev_primary_ablation_anchor_95ac397.json` | 6574597 |
| dc | `dev_primary_dc_tierA_68ca19a.json` | 157aa80 |
| AC | `dev_primary_AC_tierA_tierC_a1917e0.json` | 8a321fe |

**Controls (noasym):**
- Shuffled labels give AP 0.027 (human) and 0.069 (yeast), against prevalences of 0.0279 and 0.0602.
- Species can be identified from the features with AUROC 0.955 (Tier A) and 0.998 (A+B).
- The human source library can be identified with AUROC 0.940 (A+B).

**Tier A is mostly ascertainment, not biology (D12; ablation file).**
- Tier A's human gains disappear within testing-study strata.
- Similarity features alone give a negative T_pf.
- Adding the family-context group gives LR a T_pf of 0.084, but its stratum-averaged AP is 0.182, below 0.201. Tier A's apparent signal therefore lives in the family-context features and tracks screening design.
- Deconfounded training, with n_testable_studies as a training-only covariate (D12c), did not fix this.
- Identity-anchored candidates also fail the guard. The best is the LR rank-ensemble at 0.200 against 0.201.

**Tier B: adopted by rule, then reviewed adversarially (D15, D17; `reviews/results/tierB_gain_review.md`, 457ed8c).** The verdict was PARTLY ARTEFACT.
- **No leakage.**
  - Merged-homology regrouping gives T_pf 0.150.
  - Dropping the top 5 families keeps the human ΔAP at 0.241–0.262.
  - Shuffling Tier B rows within species brings T_pf back to 0.040.
- **Real, transferable signal.**
  - A+B beats A in 9 of 9 human source libraries; the CI excludes 0 in 8.
  - In the ≥3-studies stratum it beats B1 by +0.151 [0.043, 0.258].
  - A model trained on yeast alone lifts human AP from 0.083 (A) to 0.144 (A+B).
  - The driving features have the same direction of effect in both species: `tb_n_post_single`, the number of panel species with a single post-duplication ortholog (+); `tb_rel_rate`, relative evolutionary rate (−); conservation (+).
- **But ascertainment carries much of T_pf.**
  - Within exact study-set cells, A+B vs B1 is only +0.010 [−0.057, 0.063].
  - Given the ascertainment covariates plus B1, the A+B score adds +0.015.
- The model was adopted by the pre-declared rule, the caveat was carried forward, and the within-cell analysis was made a pre-declared test secondary (D17c).

**Tier C not adopted (D17b).**
- A+C vs A: ΔT_pf 0.006 [−0.000, 0.014] for LR and 0.013 [−0.010, 0.035] for GBDT.
- Over A+B, Tier C adds about +0.014 in human and nothing in yeast.
- The embedding step costs 5,982 s for the human proteome (Low Power Mode; `results/bench/plm_esm2_150m_human.json`).
- It remains an optional `--tier-c` path with a `plm` install extra (D16).

**Masking not adopted (D18; `results/dev/dev_thinning_GBDT_AB_17.json`, eea23bc).**

| condition (train → predict) | T_pf [CI] | human ΔAP |
|---|---|---|
| natural → natural | 0.143 [0.095, 0.187] | 0.247 |
| natural → thinned, c = 0.03 (removes macaque) | 0.143 | 0.247 |
| natural → thinned, c = 0.35 (human panel 77 → 65) | 0.127 [0.083, 0.168] | 0.219 |
| thinned → thinned, c = 0.35 | 0.145 | 0.245 |
| masking-augmented → thinned, c = 0.35 | 0.132 | 0.232 |

- At c = 0.35 the drop is 0.016 [0.003, 0.030]. It comes from a mismatch between the training and prediction panels, not from missing information.
- Augmentation recovered 34% of the drop, short of the 50% required, so it was not adopted.

**Tier B validation, label-free (`results/features/tierB_validation.json`, 3d5209d).**
- Spearman correlation of d_s with TimeTree divergence: human 0.78, yeast 0.87, worm 0.77, *Arabidopsis* 0.59.
- Human 2R ohnologs are recovered (WGD-signature AUROC 0.92).
- The yeast WGD depth signature is not recovered: pre-WGD *Torulaspora* sits closer than the post-WGD species (a long-branch effect). Only `tb_n_post_single` separates yeast ohnologs (AUROC 0.91).

## 6. Test-set results

Single access: 2026-10-28, 11:59:51–12:00:18 UTC, tag `final-model-v1` = 29e8fc8 (`results/final/access_log.tsv`). No secondary failed. All numbers come from `results/final/final_test.json` and `final_test.md` (a88cff1).

**Primary.** Called, labelled, primary-tier pairs; the model and B1 are scored on the same pairs.

| species | n | pos | neg | prevalence | AP model | AP B1 | AP B2 | ΔAP vs B1 [95% CI] | ΔAP vs B2 [95% CI] |
|---|---|---|---|---|---|---|---|---|---|
| human | 2,805 | 82 | 2,723 | 0.029 | 0.349 | 0.106 | 0.148 | 0.244 [0.156, 0.335] | 0.201 [0.122, 0.285] |
| yeast | 743 | 40 | 703 | 0.054 | 0.541 | 0.190 | 0.434 | 0.351 [0.050, 0.523] | 0.108 [0.010, 0.187] |

- T = **0.297 [0.131, 0.395]**, with relative gain 2.08 against the 0.20 threshold. T vs B2 = 0.154 [0.087, 0.210].
- Verdict: **CLEARLY BEATS identity**. Pre-declared underpowered (power 0.59 at T = 0.04); the observed effect is far larger than that design target.
- Labelled pairs the caller missed are excluded from the primary population (DESIGN §6.2): human 6 (0 positives), yeast 205 (2 positives).

**Pre-declared secondaries.**

| analysis | n | pos | prevalence | AP model | AP B1 | ΔAP vs B1 [95% CI] |
|---|---|---|---|---|---|---|
| human, 1 testing study | 1,885 | 10 | 0.005 | 0.010 | 0.005 | 0.005 [0.001, 0.022] |
| human, 2 testing studies | 659 | 23 | 0.035 | 0.333 | 0.280 | 0.054 [−0.146, 0.225] |
| human, ≥3 testing studies | 261 | 49 | 0.188 | 0.569 | 0.492 | 0.077 [−0.046, 0.192] |
| **human, within exact study-set cells (33 cells)** | 2,054 | 67 | 0.033 | 0.566 | 0.585 | **−0.020 [−0.098, 0.060]** |
| human, yeast-only-trained model | 2,805 | 82 | 0.029 | 0.259 | 0.106 | 0.153 [0.085, 0.244] |
| yeast, human-only-trained model | 743 | 40 | 0.054 | 0.497 | 0.190 | 0.307 [−0.031, 0.498] |
| human, without Harle/Chan-dependent labels | 2,773 | 70 | 0.025 | 0.287 | 0.092 | 0.195 [0.109, 0.289] |
| yeast, without orientation-conflict positives | 742 | 39 | 0.053 | 0.529 | 0.181 | 0.348 [0.054, 0.517] |
| human, lenient tier | | 337 | | 0.223 | 0.134 | 0.088 [0.049, 0.136] |
| human, ≥2-studies tier | | 61 | | 0.406 | 0.210 | 0.196 [0.094, 0.286] |
| yeast, lenient / strict / BioGRID tiers | | 58 / 26 / 86 | | 0.507 / 0.574 / 0.643 | 0.195 / 0.210 / 0.261 | 0.312 / 0.363 / 0.383 (all CIs exclude 0) |
| end-to-end (missed pairs at minimum score), human / yeast | | 82 / 42 | | 0.349 / 0.518 | 0.106 / 0.183 | 0.244 / 0.335 |

Leakage tertiles split test pairs by their maximum bitscore ratio to any development gene:

| leakage tertile | ΔAP [95% CI] |
|---|---|
| T1 | 0.218 [0.057, 0.342] |
| T2 | 0.194 [−0.030, 0.465] |
| T3 | 0.134 [−0.000, 0.270] |

The trend slope is −0.042 [−0.130, 0.066], so there is no evidence that the gain depends on homology to training genes.

**Unseen species.** Pre-declared as descriptive only, with no pass/fail.

| species / tier | n | pos | prevalence | AP model | AP B1 | AP B2 | ΔAP vs B1 [95% CI] | ΔAP vs B2 [95% CI] |
|---|---|---|---|---|---|---|---|---|
| *Arabidopsis*, strict (headline) | 158 | 51 | 0.323 | 0.593 | 0.598 | 0.692 | −0.005 [−0.111, 0.097] | −0.098 [−0.183, 0.003] |
| *Arabidopsis*, strict, without lethal-single positives | 155 | 48 | 0.310 | 0.597 | 0.582 | 0.687 | 0.015 [−0.101, 0.115] | – |
| *Arabidopsis*, broad | 297 | 190 | 0.640 | 0.778 | 0.757 | 0.810 | 0.021 [−0.023, 0.062] | −0.032 [−0.068, 0.006] |
| *C. elegans* | 135 | 16 | 0.119 | 0.358 | 0.180 | 0.187 | 0.178 [0.064, 0.349] | 0.171 [0.061, 0.328] |

**What the test means.**
1. **The pre-registered criterion is met with a wide margin in both species.**
2. **The human headline (ΔAP 0.244) is mostly a between-cell effect.**
   - Human positives concentrate in pairs screened in many studies and cell lines.
   - The model's features (family context, conservation, single-ortholog patterns) also predict *which pairs get screened*: Tier B predicts screening with AUROC 0.84 (review).
   - Within identical study-set cells the model does not beat identity (−0.020 [−0.098, 0.060]); in the ≥3-studies stratum the gain is +0.077 [−0.046, 0.192].
   - The honest human claim: a large improvement in ranking pairs as they are screened today, but no demonstrated gain among pairs screened under identical designs.
   - Two results argue that part of the signal is biology: the yeast-only model lifts human AP to 0.259 without seeing any human ascertainment, and the worm gain is +0.178.
3. **Yeast has no library ascertainment.** Its gain (0.351 [0.050, 0.523]) is large, but the CI is wide. Test B1 is 0.190 against 0.567 on development. The cause is unexplained: no yeast group was forced into development (`results/splits/split_summary.json`), so it is presumably group composition. It was not investigated after the access.
4. **Unseen species split.** *Arabidopsis* shows no gain over identity and sits below B2 on both tiers. *C. elegans* (16 positives) shows a clear gain.

### 6a. Robustness to degraded proteomes (pre-registered, exploratory, development labels only)

**Design** (DESIGN §6.5; declared in b58dbc2, results in a84a4f1, `results/dev/dev_robustness_GBDT_AB_17.json`):
- Each proteome was degraded with seed 11 in one of these ways:
  - fragments: 20% of proteins truncated to 50–80% of their length;
  - split: 5% of genes split in two;
  - fused: 2% of genes fused with a donor gene;
  - missing: 15% of genes removed;
  - haplotig: near-identical copies (1% substitutions) of 5% of genes;
  - combined: all of the above.
- Tier A and B features were then recomputed with the shipped `predict`.
- Evaluation followed the thinning design: the model (GBDT_AB_17, A+B) is trained on natural features in each development fold and evaluated on the degraded features of the held-out fold.
- Labelled pairs that lose a gene, or are no longer called, are scored at the minimum for both the model and B1.
- The run took place after the test access and does not touch it. Yeast was run for every kind; human only for combined, for compute reasons.

| species | degradation | model AP | drop [95% CI] | B1 AP | drop | pairs at min (pos) |
|---|---|---|---|---|---|---|
| yeast | none (identity control) | 0.626 → 0.626 | 0 | 0.586 → 0.586 | 0 | – |
| yeast | fragments | 0.626 → 0.621 | 0.005 [−0.017, 0.022] | 0.586 → 0.538 | 0.048 [0.001, 0.088] | |
| yeast | split | 0.626 → 0.609 | 0.017 [0, 0.037] | 0.586 → 0.546 | 0.041 [0.008, 0.079] | |
| yeast | fused | 0.626 → 0.628 | −0.002 [−0.008, 0.004] | 0.586 → 0.591 | −0.005 | |
| yeast | missing (15%) | 0.626 → 0.442 | 0.184 [0.091, 0.276] | 0.586 → 0.422 | 0.164 | |
| yeast | haplotig | 0.626 → 0.604 | 0.022 [−0.004, 0.048] | 0.586 → 0.586 | 0 | |
| yeast | combined | 0.626 → 0.393 | 0.233 [0.126, 0.327] | 0.586 → 0.329 | 0.257 | 343 (16) |
| human | combined | 0.387 → 0.228 | 0.159 [0.098, 0.225] | 0.140 → 0.111 | 0.030 [0, 0.071] | 1,130 (35) |

AP here is the per-fold mean, so it differs from the pooled development APs in §5. With both species combined, T_pf goes from 0.143 to 0.090, a drop of 0.053 [0.008, 0.096].

**Reading.**
- Fragments and split genes hurt identity more than the model, so the model tolerates fragmented gene models better than identity does.
- Fusions do nothing.
- Haplotig copies cost the model a little (0.022) and identity nothing. Plausibly, the extra near-identical partners disturb the family-context and rank features.
- Missing genes dominate the loss for both methods. This is unavoidable: a pair whose gene is absent cannot be scored.
- In human (combined), the model loses more AP than identity (0.159 vs 0.030), but still ranks above it after degradation (0.228 vs 0.111). Its margin over identity roughly halves, from 0.247 to 0.117.
- The degradation levels were chosen in the design, not estimated from real assemblies. The size of the effect on real draft proteomes is unknown.

## 7. Runtime and memory

**Full human proteome, end to end (the deliverable; `results/bench/e2e_human.json`, `scripts/bench/e2e_bench.py`).**

| item | value |
|---|---|
| command | `slparalog predict human.fasta --model results/frozen/model -t 10` (Tier A self-search, Tier B panel search against 80 proteomes, scoring) |
| input | 20,770 proteins; 686,225 pairs scored |
| start state | empty workdir and empty panel-hit cache |
| wall time | **5,990 s (99.8 min)**, of which the Tier B panel search ≈5,277 s |
| CPU time | 47,240 s user + 535 s system |
| peak memory | **2.52 GB** process tree, 2.39 GB largest single process |
| peak scratch disk | 1.52 GB (deleted afterwards) |
| machine | Apple M4, 10 cores (4P + 6E), Low Power Mode off |
| load | 1-min load average 4.6 at start, climbing to about 60 during the run. Mostly the run's own MMseqs2 threads, plus background system daemons and one unrelated user process |
| commit | 211c264 |

- **One-off setup.** `slparalog setup --panel` built the 80 per-species MMseqs2 databases in 16 s, with 0.16 GB peak memory and 2.2 GB on disk (`results/bench/e2e_setup.json`). Network download time is not included, because the checksummed panel files were already local.
- **Yeast end to end:** 1,115 s and 1.46 GB (`results/bench/e2e_yeast.json`).
- **Consistency.** For both species, the shipped `predict` output has the same pair set as the evaluated feature tables. Its scores are identical to the frozen model's scores computed from those tables: maximum difference 0.0 (`results/bench/e2e_consistency_*.json`).
- **Against the design target:** the target was ≤2 h and ≤8 GB for human, and both are met.

**Per-tier timings from development**, under heavy parallel load and partly in Low Power Mode (LPM), for reference only:

| step | species | wall | peak RSS | LPM | file |
|---|---|---|---|---|---|
| Tier A | human | 7,057 s | 1.66 GB | 1→0 | `results/bench/tierA_human.json` |
| Tier B panel search | human | 12,832 s | 4.76 GB | 1→0 | `tierB_human_search.json` |
| ESM-2 150M embeddings (Tier C, not shipped by default) | human | 5,982 s | 1.29 GB | 1 | `plm_esm2_150m_human.json` |

## 8. Limitations

1. **Ascertainment of the human labels.**
   - The human test gain is largely between screening-ascertainment cells; within cells it is −0.020 [−0.098, 0.060].
   - Screens chose pairs non-randomly, and sequence features predict that choice.
   - Human AP should not be read as precision on unscreened pairs.
2. ***Arabidopsis.***
   - There is no gain over identity, and B2 is higher than the model.
   - The plant panel is sparse, and d_s agrees less with TimeTree there (0.59).
   - The Cusack labels measure a different endpoint (genetic redundancy in a plant) and have high prevalence (0.32 strict).
   - The non-gain was not used to retune anything (D20).
3. **Yeast long-branch WGD.** Tier B depth does not recover the yeast WGD; only `tb_n_post_single` carries that signal (D14).
4. **Statistical power.**
   - The comparison was pre-declared underpowered.
   - The yeast test has 40 called positives, and its ΔAP lower bound is 0.050.
   - The unseen sets are small (worm: 16 positives).
5. **Panel density.**
   - The model loses a little (0.016) when the nearest panel species are absent.
   - Insects are untested.
   - The tool emits QC warnings when the nearest panel species is distant or the core set is small.
6. **Timings.** The end-to-end benchmark ran on a laptop that also carried system daemons and one unrelated user process. Development per-tier timings were taken under heavy load.
7. **The release model has no evaluation.** `AB-gbdt-332d4d20` (`results/release/model`, f8cf438) refits the frozen spec on all 8,959 called primary pairs from development + test: human 7,220 (205 positives), yeast 1,739 (100 positives). No held-out evaluation of it exists (`results/release/release_model.json`).
8. **Draft-quality proteomes.**
   - With the combined simulated degradation, the margin over identity roughly halves in human (§6a).
   - Haplotig-rich or incomplete assemblies will score worse.
   - The tool's QC flags (`near_identical`, `fragment_suspect`) mark some of these cases.
9. **Uncalled pairs.**
   - The caller misses 22% of labelled yeast test pairs (205 of 948, including 2 positives), but only 6 human pairs.
   - End to end, yeast AP is 0.518 instead of 0.541.
10. **Score semantics.** `sl_score` is a ranking score, not a calibrated probability.
11. **Pre-registered items not done.**
    - WGD vs SSD stratified test AP and expected precision at genome-wide prevalence (DESIGN §6.2) were not computed in the single access. They cannot be added now without a second access.
    - Tier D was not built.

## 9. Process and integrity

**Pre-registration.**
- Three design review rounds (`reviews/round1–3/`).
- Label rules were committed before any label code ran.
- Splits and power were frozen before any score existed (D9; `results/mde.json`).
- Each development batch had a declare-only commit before it ran (e.g. b2cf739, 2ab485d, 6e3c00a, 048c11a).
- The test secondaries were fixed in D17c before freezing.
- `final_test.py` passed independent review: first READY WITH FIXES (R1 model pinning, R2 missing secondaries, R3 a failing secondary could lose the primary), then READY (`reviews/results/final_test_review.md`).

**Single access.**
- There was exactly one primary access, logged with the git head, model hash, config hash and sealed-manifest hash (`results/final/access_log.tsv`).
- No model or analysis changed after it (D20).
- The access log lives at `results/final/access_log.tsv`, not at the `results/test_access_log.tsv` path named in DESIGN §6.3.

**Deviations and incidents.**

| note | what happened | effect |
|---|---|---|
| D7 | Label-blind feature code started before the splits were frozen | The power simulation uses no features |
| D10 | A label reviewer computed one identity-only AP on unsplit labels | The value was never recorded; briefs now forbid scoring labelled pairs |
| D13, D18 | Production runs started from uncommitted code | Later committed and verified (D18: natural → natural reproduced 0.1432 exactly); Tier C was re-run from a clean tree with identical values (00991c6) |
| 239c3fd | The whole-registry hash invalidated unchanged feature tables | Replaced by per-table column-definition hashes; the existing tables were re-stamped after a check against git history |
| D19 | D17 had misstated the primary population | The final-test review caught it before the access |
| missing source | Kuzmin 2020 was not obtained, and *S. pombe* and the Parrish Cell Reports tables need manual browser downloads | Parrish calls come from the in4mer harmonised LFCs |

**The D12 rule refinement, made after Tier A development results.**
- After the first Tier A development report showed that the human "gains" were learned ascertainment, the PI added:
  - the human stratum-averaged-AP guard;
  - per-fold paired T_pf as the selection statistic;
  - the deconfounding option.
- The adoption condition was later restated (95ac397).
- Only development data had been seen, and both changes made adoption stricter: they rejected every Tier A configuration. Every adoption decision is the same under both wordings (GBDT A+B: T_pf − SE = 0.119; CI lower bound 0.095).
- It is nonetheless data-informed rule-making. The retro lists it as a pattern to avoid (`reviews/retro/M3_retro.md`).

## 10. Follow-ups and remaining ideas

1. **Ascertainment-robust human evaluation and training.**
   - A within-cell guard during development.
   - Unbiased human labels, e.g. all-paralog screens with full non-hit tables.
   - Objectives that compare pairs only within study-set cells.
2. **Plants.**
   - Add plant proteomes to the panel by a written rule (gymnosperm, lycophyte, more angiosperm families).
   - Re-check d_s calibration.
   - Analyse the *Arabidopsis* non-gain label-free first.
3. **Tree-based duplication depth.** Use gene-tree or species-overlap reconciliation over the panel hits instead of bitscore depth bounds, to fix the yeast long-branch WGD failure.
4. **Calibrated probabilities.** Per-species recalibration to a stated prevalence, and expected precision at genome-wide prevalence on a future independent test set.
5. **More unseen species.**
   - *S. pombe* E-MAP (manual download pending).
   - Insects, if a screen with tested negatives appears.
   - A second plant set.
6. **Other representations.**
   - Tier D (ProstT5 3Di + Foldseek).
   - ProtT5 or ESM C embeddings. Low priority: Tier C added about 0.01 here.
7. **Runtime.**
   - Tier B dominates (≈88% of human wall time). Panel search could be restricted to Tier A genes, or species could be searched in parallel with fewer threads each.

## 11. Reproduction

Environment: `.venv` (Python 3.12), MMseqs2 at `bin/mmseqs` (pinned in `config/tools.yaml`), seeds in `config/seeds.yaml`.

```bash
# install and use the tool (README)
python -m venv .venv && . .venv/bin/activate && pip install -e ".[dev]"   # add .[plm] for Tier C
slparalog setup --panel                        # 80 panel proteomes and per-species MMseqs2 DBs, sha256-verified
slparalog setup --model results/release/model  # install the release model as the default
slparalog predict proteome.faa -o pairs.tsv -t 10
pytest

# data and labels
.venv/bin/python scripts/data/download.py
.venv/bin/python scripts/data/build_gene_proteome.py
.venv/bin/python scripts/data/build_panel.py all
.venv/bin/python scripts/labels/build_human_labels.py
.venv/bin/python scripts/labels/build_yeast_labels.py
.venv/bin/python scripts/labels/apply_universe.py
.venv/bin/python scripts/eval/power_sim.py
.venv/bin/python scripts/labels/make_splits.py all

# features
slparalog predict human.faa -o human_feats.parquet --features-only --tier-b -t 10   # likewise yeast, arabidopsis, worm

# development runs (see results/experiments.tsv for every declared configuration)
.venv/bin/python scripts/eval/dev_report.py --help

# frozen model (D19) and the single test access (already done; must not be repeated)
slparalog train --features human=tierA.parquet --features human=tierB.parquet --features yeast=... \
    --labels human=data/labels/dev/human_dev.tsv --labels yeast=data/labels/dev/yeast_dev.tsv \
    --feature-set A+B --model-type gbdt --seed 11 -o results/frozen/model   # GBDT_AB_17 params
# scripts/eval/final_test.py

# release model
.venv/bin/python scripts/eval/build_release_labels.py
.venv/bin/python scripts/eval/write_release_model_record.py

# benchmarks
.venv/bin/python scripts/bench/e2e_bench.py setup
.venv/bin/python scripts/bench/e2e_bench.py predict human
```

Every development and test JSON records its git hash, input hashes, seeds and bootstrap settings. The frozen and release model directories record their training files, label checksums and parameters.
