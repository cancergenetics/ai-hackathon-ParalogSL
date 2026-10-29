# Independent review: human paralog SL labels

- reviewer: independent review agent (Claude Opus 5.5). It did not write the label code, the config or the labels.
- producer: label-builder agent (human), commit `1a3c795`
- commit: `1a3c795` (labels, code, QC). Config history: `513c1e5` pre-registration, then `c4e5039` D4 amendment, then `1a3c795` §7 post-run amendments. Reviewed at HEAD `bfae7f8`. No human-label file changed between `1a3c795` and `bfae7f8`.
- rerun_command: the two scripts and `config/labels_human.yaml` were copied into a clean scratch root `data/interim/review_rerun_human/`. In that root, `data/raw` and `data/interim/universe` are symlinked to the project copies. Then `.venv/bin/python data/interim/review_rerun_human/scripts/labels/build_human_labels.py` was run (121 s). The script has no output-path option, and `ROOT` comes from `__file__`, so this writes nothing into the tracked tree.
- inputs_sha256 (16-character prefixes):
  - config `264f26474b998895`; build script `aee6eafe7e4aaa30`; parsers `cc680696bf93105a`
  - Supp_table1 `7058f8d41a06a308`; Dede S2 `122d0d7271d1dbc8`; Thompson SD1 `6ab438558c086c39`; Ito `c3dd8427d78b0fca`; CHyMErA S8 `af0439c76a37983b`
  - Chou T1 `8ced00d405f65a5d`; Chou gene LFC `409ae4e4703582b0`; Klingbeil S6 `ff5c25c4e4575e31` / S2 `897016827bef31ab`
  - Harle `297ed1ee52a753ea`; Chan T6 `6fe74cc8958e899a`
  - UniProt idmapping `214ef098abf898d2`, gene2acc `0e091cecfad9e201`, fasta `64299d01e8f41fb6`
  - Compara universe (git-ignored) `68e7e62ab17922f1`
- match: **yes, exact (tolerance 0).**
  - `data/labels/human_labels_all.tsv` sha256 `84e2d7e5…`, identical to the committed file.
  - `human_calls.tsv` `396eea18…` and `human_mapping_report.tsv` `32e36680…`, both identical.
  - `human_label_qc.json` differs only in the key order inside `ensembl_multi_xref_resolution` (issue 7). The values are identical, apart from `git_commit_at_build`.
- Known-answer and rule tests: `tests/test_label_known_answers_human.py` and `tests/test_human_label_rules.py` give 20 passed.

## Verdict: **ACCEPT WITH FIXES**

The code implements the pre-registered rules (R2-5/6/7/15 as amended by D4). An independent reimplementation reproduces every label exactly. Re-running the builder gives byte-identical labels.

There is one real bug: one cell line is not merged across studies. It changes no tier. The other open points are not bugs:
- a strong ascertainment structure that is inherent to the sources and the rule. It must be carried into evaluation.
- the sanctioned Harle/Chan fallback. It affects 22% of positives and needs a flag and a sensitivity analysis.

Tier counts (from `results/labels/human_label_qc.json`): 11,125 pairs.

| Tier | Positive | Negative | Ambiguous | Untestable |
|---|---|---|---|---|
| Primary | 222 | 8,040 | 2,046 | 817 |
| Primary, inside the Compara universe | 205 | 7,026 | | |

Prevalence among labelled primary pairs is 0.027.

## Checks performed

### 1. Code vs pre-registered config and normative rules

**R2-5 line call** (`call_rows`), all correct:
- m = median and s = 1.4826·MAD of per-gene-deduplicated singles; T = max(2s, 0.5).
- testable = (a+b) > −T and no study flag.
- SL = testable ∧ GI_strict ∧ ab ≤ −T.
- clear negative = testable ∧ ¬SL ∧ (¬GI_relaxed ∨ ab > −T/2).
- Missing LFC or significance gives *untested* (X), not negative.
- The no-LFC rule for Harle and Chan follows config §2.

**R2-6 per-study rules**, all match config:
- Every study's GI_strict and GI_relaxed thresholds.
- Ito: the 90th-percentile rule is dropped; pairs with no GEMINI score give X; anchor-screen tables are not read.
- CHyMErA: GI at either timepoint with effect < 0, DKO from the late sheet, singles from `Paralog*.Mean`, not `expected_LFC`. Torin sheets are not read.
- in4mer: `MelJuso_sel_T21` is not read.
- Klingbeil domain collapse: SL if more than 50% of testable domain rows are SL, done per symbol pair × line.
- Thompson authors' rule plus filter flag.
- Parrish fallback ZdLFC: a Chou re-implementation, checked against Chou's notebook on disk.

**R2-7 aggregation and R2-15:**
- Merge rule: S + B gives S; S + N gives D; all N gives N; U and X are ignored.
- n_testable counts merged S/N/B/D lines; penetrance = n_sl / n_testable, over **testable lines only**.
- n_studies_sl counts only studies whose S call sits on a merged-S line.
- The primary, lenient and two-studies tiers match the config, including the D4 floor of 0.10 on the two-studies tier.

**Independent re-aggregation.** I re-wrote the line-merge and tier rules from the config text and ran them on `human_calls.tsv`. Result: 0 mismatches over 11,125 pairs in `label_primary`, `label_lenient`, `label_two_studies`, `n_testable` and `n_sl`.

**End-to-end spot check from raw tables.** I wrote fresh code: no `human_parsers` import, my own m/T computation, and my own Klingbeil domain collapse. It covers 12 pairs from all 9 studies. Every value matches the built label file: label, n_testable, n_sl, n_borderline, penetrance and studies_sl.

| Pair (studies) | Primary | n_test / n_sl / n_bord |
|---|---|---|
| GSK3A–GSK3B (Dede, Parrish, Ito, in4mer, Klingbeil) | 1 | 33/19/5 |
| SMARCC1–SMARCC2 (Harle, Chan, Klingbeil + 4 more) | 1 | 49/31/10 |
| FAM50A–FAM50B (Harle, Chan, Parrish, Thompson, in4mer) | 1 | 28/9/12 |
| SPATA2–SPATA2L (CHyMErA) | 1 | 2/2/0 |
| YWHAG–YWHAZ (Ito, in4mer) | 1 | 8/2/3 |
| AGPAT3–AGPAT4 (CHyMErA, in4mer, Parrish) | 0 | 7/0/0 |
| AOX1–XDH (Harle, Chan, Thompson) | 0 | 13/0/0 |
| RAPGEF2–RAPGEF6 (Dede, in4mer, Parrish) | 0 | 6/0/0 |
| CYTIP–GRASP/TAMALIN (Thompson) | −1, borderline_excess | 3/0/1 |
| MAP4K2–MAP4K5 (Ito, in4mer, Klingbeil) | −1, single_line_sl | 26/1/1 |
| ATP1A3–ATP2C2 (Klingbeil) | −2 | all X |
| TRIM11–TRIM17 (Ito) | −2 | all X |

**Duplicates and canonicalisation:**
- No duplicated (gene_a, gene_b) pair, and gene_a < gene_b in every row.
- No symbol pair occurs twice.
- No unordered duplicates in any Supp_table1 sheet.
- No Harle pair × line duplicates.

**Gene mapping:**
- Cross-check against `data/raw/hgnc/hgnc_complete_set_2026-10-06.tsv` (on disk, but not used by the builder): 22,044 symbol→ENSG uses agree, 12 disagree (7 symbols), 226 symbols are not HGNC-approved.
- All 8,068 label gene IDs exist in the Ensembl 116 peptide FASTA.
- Thompson SD3 Ensembl IDs: 2,567 same, 8 different. Six of the 8 differ because Thompson's older olfactory-receptor IDs are absent from Compara 116, while ours are present.

### 2. Clean re-run
Exact match, as stated in the header.

### 3. Post-run amendments
- `git diff c4e5039 1a3c795 -- config/labels_human.yaml` only appends §7. No threshold, tier, study rule or line map changed.
- gene2acc layout fix and OpenTargets/reference-proteome haplotype resolution: these touch ID mapping only. They were triggered by mapping-rate failures, not by label outcomes. The HGNC agreement of 99.95% above supports them.
- Ito header fix (header row 3, 0-based): a parser bug. The Table 10 size test (4,859) pins it.
- Chou test skip: `test_chou_zscore_reimplementation_matches_t1` restricts the library to pairs whose single-gene rows exist. This mirrors `calculate_fcOBS_fcEXP_dLFC` in Chou's notebook on disk, which skips them on KeyError. The tolerance stays 0.01, and no skip/xfail marker is used. Justified.
  - It is not recorded in config §7, which is fine because it is a test, not a rule.
  - The tests were first committed together with the code, so git cannot show *when* this restriction was added.

### 4. Bias and bug audit
- The Klingbeil domain collapse is per gene pair × line, over testable rows only. My independent collapse agrees.
- Collapse audit of duplicated rows: no non-Klingbeil group mixes U with a testable call. CHyMErA duplicates are UU 3, SS 1, BBB 1, BB 1.
- The same physical line counted twice: see issue 1.
- Ascertainment: see issue 2.

### 5. Harle/Chan DKO-depth gap
See issue 3.

## Issues

| # | Severity | Issue | Evidence | Fix |
|---|---|---|---|---|
| 1 | **Minor (bug)** | Same cell line counted twice. Harle S2 gives the Cellosaurus IDs of HPAF-II and CFPAC-1 with a leading **non-breaking space** (`'\xa0CVCL_0313'`). Harle HPAF-II therefore becomes line `name:\xa0CVCL_0313` and is not merged with Klingbeil HPAFII (`CVCL_0313`). | 46 pairs are in both studies. After stripping the space and re-aggregating, n_testable falls by 1 for 26 pairs and n_sl changes for 1 pair. **0 tier changes** in primary, lenient or two-studies. All other overlaps merge correctly: A549, A375, MeWo, NCI-H1299, AsPC-1, RPE-1, MelJuso, MEL202, K562. | Strip `\xa0` and whitespace in `read_harle_s2`. Add a test that every `line_id` matches `^CVCL_\w+$` or `^name:[A-Z0-9]+$`. Rebuild. Counts change, tiers do not. |
| 2 | **Major (interpretation; inherent, not a bug)** | Label status is strongly predicted by how much a pair was screened (see the breakdown below this table). | Reviewer audit `data/interim/review_rerun_human/rev/gap_audit.json` (uncommitted). | No code fix. Before any model result, pre-declare: (a) an "ascertainment baseline" scored like B1 (n_testable and n_studies_testable as the score), reported next to B1; (b) primary AUPRC stratified by n_studies_testable (1 vs ≥2); (c) the source-library control from DESIGN §2.3; (d) a sensitivity analysis without negatives that have n_testable = 2. Add the audit tables to the QC script. |
| 3 | **Moderate (sanctioned gap)** | Harle and Chan calls have no DKO-depth term (ab ≤ −T) and no (a+b) > −T testability. Only the authors' single-depletion and expression flags are used. In the LFC studies the DKO-depth term removes a large share of testable GI-strict calls (see the breakdown below this table). | `data/interim/review_rerun_human/rev/harle_reliance.json`, `gap_audit.json` | Acceptable for development labels: R2-5 explicitly allows "authors' calls only", and config §2 pre-registers it. Required: (a) a label column `positive_depends_on_harle_chan`, with primary AUPRC also reported without those 48 pairs; (b) PI decision on whether to implement the R2-5-preferred recount from S3/T5 counts (Harle repo, MIT) **before** splits and labels are frozen. That is the pre-registered preferred route, not tuning, as long as no model result exists. |
| 4 | Minor | Thresholds are hard-coded in `human_parsers.py`: Dede −3/−2, Ito 0.05/0.25, Harle −0.5/0.01/−0.3/0.1, and others. The module docstring says thresholds are only read from config. | All values match config today (checked line by line). | Read them from config, or add a test that asserts the parser constants equal the config. |
| 5 | Minor (latent) | CHyMErA `score=min(es, ls)` uses Python `min`, which returns NaN or a value depending on argument order. A NaN early effect would give X. | No NaN effect sizes in S8 (checked), so there is no effect today. | Use `np.nanmin`. |
| 6 | Minor | `collapse_dups(['S','U'])` returns S. The config says that disagreeing duplicate rows give B "if any row is testable". | 0 occurrences in the data (collapse audit). | Align the code with the config, or document the reading. |
| 7 | Minor | The QC JSON is not byte-reproducible. The insertion order of the `ens_resolution` Counter follows set iteration under hash randomisation. | The re-run diff is one moved line with the same values. | `sort_keys=True`, or sort that dict. |
| 8 | Minor | 7 symbols map to a different ENSG than the HGNC complete set gives: CT47A5, CTAG1B, PNLIPRP2, PRICKLE4, PRSS50, UBE2NL, ZBTB8B. Config §1 says the HGNC set is "not on disk", which is now stale: it was downloaded 00:26, after pre-registration. | HGNC cross-check above. | Report as a known mapping discrepancy. Optionally fix these 7 as a documented ID-only amendment. |
| 9 | Minor (reproducibility) | `in_compara_universe` comes from a git-ignored file that the script picks with `sorted(glob('*human*.tsv*'))[0]`. A second matching file would silently change the universe. | `data/interim/universe/human_paralog_pairs.tsv` (sha `68e7e62a…`). | Pin the path and sha256 in config. |

### Issue 2 in numbers (labelled primary pairs)

**By number of studies in which the pair is testable** (positive rate among labelled primary pairs):
- 1 study: 0.67% (n = 6,257)
- 2 studies: 4.0% (1,385)
- 3 studies: 10.3% (331)
- 5 studies: 39% (77)
- 6–8 studies: 56–58%

**By number of testable lines.** Used as the only score, n_testable gives:

| Pair set | AUROC | AUPRC | Prevalence |
|---|---|---|---|
| All labelled | 0.70 | 0.127 | 0.027 |
| Compara universe | 0.72 | 0.142 | 0.028 |

By line count:
- 2 testable lines: 7 positives vs 2,209 negatives. These 2,209 are 27% of all negatives, mostly 2-line in4mer sub-library pairs.
- Single-study in4mer pairs: positive rate 0.16% (n = 2,447).

**By library membership** (Harle or Chan): positive rate 26% vs 1.7% for non-members (AUROC 0.68).

**Why this is inherent and not a bug.** The rule needs at least 2 SL lines. For single-study pairs it also needs at least ceil(0.25·n) SL lines, which makes 2-line pairs need 2/2. Pairs tested by many labs are the ones that were prioritised as likely SL. A sequence model can exploit "is a much-screened, famous pair" through features such as identity or closest-pair status.

### Issue 3 in numbers

**Share of testable GI-strict calls removed by the DKO-depth term, by LFC study:**
- Ito 76%, CHyMErA 78%, in4mer 48%, Parrish 34%, Klingbeil 16%
- Thompson 5%, Dede 0%

**How much the positives depend on Harle and Chan:**
- **48 of 222 primary positives (22%)** are positive only because of Harle/Chan calls. Without them, 32 are ambiguous, 9 are absent, 5 untestable and 2 negative.
- 29 positives have SL support only from Harle/Chan.
- Within-study rate of pairs SL in at least one line: Harle 26%, Chan 21%, other studies 3–12%.

**Same-line agreement:** where Harle calls S and another study has a testable call on the same line and pair, the other study calls N in 34 of 86 cases (40%).

## Not verified
- Whether the label code ran before the config commits. Commit times are 00:23:50 and 00:24:06 for the config, and 00:39 for the code. I accept the provenance as plausible but cannot prove it.
- Thompson Supp_table1 timepoint (R2 risk 1, already flagged in config).
- Harle/Chan counts-based recomputation (not implemented).
- Whether the Compara universe file itself is correct; that belongs to a separate review.
