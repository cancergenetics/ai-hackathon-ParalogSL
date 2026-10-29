# Independent review: yeast paralog SL labels

Reviewer: independent agent (Opus 5.5). I did not produce these labels. Date 2026-10-28.
Reviewed at commit 9b4614a (labels and QC) and a79394b (code). `config/labels_yeast.yaml` v1 is unchanged since it
landed in 2c7ab4e. The config sha256 in the QC JSON is a6132bf4…, and 415356b records the shared-index race.
Inputs: DESIGN §2.3, the D4 note, reviews/round2/data.md R2-1/2/3/9/10/11, the builder, the parsers, `apply_universe.py`,
the tests, `data/labels/yeast_labels_all.tsv`, `results/labels/yeast_label_qc.json` and the yeast part of
`results/labels/known_answer_gaps.json`. I did not read `data/labels/sealed/` or the Dennler quarantine files; the builder's
whitelisted read of the quarantine files ran only inside its own re-run.
Reviewer scripts and outputs are in `data/interim/review_yeast/work/` (gitignored). Every number below comes from those
outputs or from the committed QC JSON.

## Verdict: ACCEPT WITH FIXES (all fixes minor; none changes a label)

The labels follow the pre-registered config and the normative R2 rules:
- An independent re-implementation, written from the config text, gives the same tier labels as the builder for every one of the 9.1 M pairs (§1).
- 15 pairs traced by hand from the raw rows all agree.
- A re-run of the builder gives byte-identical outputs.
- The +10.9 % stringent-count gap is explained by which raw files Costanzo's count covers, not by our rules.

The issues below are about reporting and interpretation. Two of them need pre-declared sensitivity analyses downstream: 7 tie-broken orientation-conflict positives, and ribosomal pairs making up 23 % of positives.

## 1. Code vs pre-registered rules (independent re-implementation)

`work/reimpl.py` streams `SGA_NxN.txt` and `SGA_ExN_NxE.txt`. It applies the allele whitelist by strain-ID suffix (sn query × dma
array, DMA26/DMA30), maps genes through the strain table, and computes E from the row SMFs with the strain-table fallback. It then applies the
POS/NEG/MID/testability rules and R2-3 aggregation (overall and per orientation), the lenient, strict and Dennler-comparable tiers,
Dean 2008 (its own pdfplumber regex parse), the BioGRID LT SL tier and the 50 kb linked rule (Ensembl 116 GTF), and finally compares all pairs with
`data/interim/labels/yeast_pair_labels_unrestricted.tsv`.

**Result (`work/reimpl_out/reimpl_report.json`, peak RSS 6.35 GB):**
- My pair set (9,114,384) is identical to the builder's once 4,829 builder-only rows are set aside. Those rows carry only BioGRID LT Synthetic-Growth-Defect evidence or essential-member SL evidence, and all of them are `untested` in every tier.
- On the shared pairs there are **0 differences** in `primary`, `lenient`, `strict`, `biogrid`, `dennler_comparable`, `orientation_conflict` and `linked`.
- My `n_T` differs on 602 rows. Those are the non-SGA pairs, which I fill with −1 and the builder with 0, so this is not a label difference.
- My independent Dean parse gives 289 / 201 / 49, the strain classes are exact, and my BioGRID filter gives 3,509 qualifying records and 2,840 pairs, the same as the QC.
- All three nesting relations hold (strict ⊂ primary ⊂ lenient positives; primary ⊂ biogrid positives; lenient negatives ⊂ primary negatives = strict negatives).

Measurement level (QC JSON `measurements`): 12,434,819 null × null rows, all from NxN; ExN_NxE has 0 sn×dma rows, as R2-2 says.
- 36,556 POS, 10,580,008 NEG, 610,368 MID, 1,207,887 UNTESTABLE, 0 MISSING.
- The rows with E > 0.5 that pass the v2 POS rule are exactly the 36,556 of data.md R2-1, and 17,956 more fall at E ≤ 0.5.
- Known-answer tests: `pytest tests/test_label_known_answers_yeast.py` gives 10 passed and 3 xfailed (strict), in 760 s. `tests/test_yeast_label_rules.py` gives 8 passed.

### Hand spot-checks from raw rows (`work/spot_rows.tsv`)
| pair (gene names from strain table) | raw evidence | expected by rules | file label |
|---|---|---|---|
| *rows redacted in the public export: per-pair calls and interaction scores from a source without a redistribution licence (see DATA.md)* | | | |

## 2. Reproducibility
I re-ran `build_yeast_labels.py` (which calls `apply_universe.py`) from HEAD code in a scratch copy, `data/interim/review_yeast/rerun/`.
- `yeast_labels_all.tsv` (sha256 1483b2d2…), the unrestricted pair table (e9908a8e…) and `yeast_calls.tsv` (cc198668…) are byte-identical to the committed run.
- The QC JSON differs only in `git_head_at_run` and `runtime_s`.
- Peak RSS was 7.77 GB and wall time 2,711 s (Low Power Mode, run concurrently with other reviewer jobs).

## 3. Known-answer gaps (xfail): the stringent count, +10.9 %
Diagnostics are in `work/stringent.json` and `work/stringent2.json`.

**Ruled out:**
- duplicate rows: 0 duplicated (query strain, array strain, Arraytype/Temp) rows in any of the four files;
- NaN p or ε: 0;
- `p ≤ / ε ≤` instead of `<`: changes the total by under 0.1 %;
- allele-level collapsing: 421,499 unique ordered allele-name pairs vs 421,628 rows.

**Cause: the readme counts do not include `SGA_DAmP.txt`.**
- BioGRID 5.0.262 holds 279,599 negative and 47,191 positive records for PUBMED:27708008. That is 326,790 records, exactly the readme's collapsed gene-level count, so BioGRID is the published set.
- Share of stringent rows whose (gene pair, sign) appears in that BioGRID set: ExE 99.1 %, ExN_NxE 93.9 %, NxN 90.0 %, DAmP 25.1 %.
- Only 6 BioGRID (pair, sign) combinations occur only in DAmP.
- Without DAmP, the stringent count is 156,162 + 165,772 + 62,608 = **384,542 (+1.18 % vs 380,059)**.
- The intermediate counts without DAmP are 545,600 negative (−0.8 % vs ~550,000) and 344,507 positive (−1.6 % vs ~350,000), both inside the pre-registered ±5 %.
- The remaining ~1.2 % in the stringent count is unexplained, for example by unstated post-hoc filtering.

**Effect on our rules: none.** Labels use only sn × dma rows of NxN, which are inside the published set, and no label
threshold depends on these counts. The `nxn_dmf_min` xfail is also benign: −0.1525 is reproduced on sn × dma rows.
Process note: according to its own comment, `test_nxn_dmf_min` was split out with an xfail marker after the first run. The assertion and tolerance are unchanged and G4 sanctions this, but it is the kind of edit AGENT_RULES §10 asks to keep visible.

## 4. Universe application
- **Universe file:** 5,650 Compara-116 pairs (2,701 within_species_paralog, 2,949 other_paralog), canonical, with no duplicates. All 2,585 genes are ordered-locus-name IDs. Matching the strain table needed no case or format fixes: 0 IDs match only case-insensitively and 0 have a non-ORF format.
- **Labelled table:** 2,875 rows. 2,715 have SGA measurements. The other 160 are `untested` rows carried only by Dean (69) or BioGRID-only evidence (91).
- **Universe pairs with no row (2,775):**
  - 1,173 have a gene that is not in the Costanzo strain table;
  - 1,238 have no possible sn × dma combination (TS/DAmP-only, i.e. essential genes);
  - 364 had an sn × dma combination that was not measured. Of these, 76 are ≤ 50 kb apart, 53 are 50–200 kb, 39 are 200–500 kb, 10 are > 500 kb on the same chromosome, and 186 are on different chromosomes. Some of them may be reciprocal pairs that the source silently removed because of opposite-sign conflicts (Costanzo readme).
- None of these 2,775 is a negative, as R2-4 requires. Every primary pos or neg with n_T = 0 comes from Dean (14 pos, 17 neg).
- **Linked rule:** SGA measurement coverage is depleted out to ~200 kb (fraction measured 0.018 at 100–200 kb vs 0.40 at 200–500 kb). The pre-registered window is 50 kb, but in the universe the only SGA-measured pair at 50–200 kb is YGR132C–YGR231C (196 kb, neg, Dean also neg). The window choice therefore changes nothing here. The universe has 14 linked pairs: 13 untested and 1 Dean negative.

## 5. BioGRID tier
- Filtering follows the config: genetic, Synthetic Lethality, Throughput contains "Low Throughput", both organism 559292, systematic names present, not self, both genes with an sn/dma strain.
- 3,509 qualifying records collapse to 2,840 pairs; a pair counts once however many records it has. 4,774 SL records are excluded because a member is essential.
- On universe pairs, 14 of the 432 LT SL records are "High Throughput|Low Throughput". This is allowed by the "contains" wording.
- None of the records on universe pairs comes from Dean 2008, so Dean and BioGRID evidence are not double-counted.
- The largest single source is Steffen 2012, with 49 records on universe pairs.
- `biogrid_only` = 123: 51 primary neg, 31 amb, 24 untested, 17 untestable.
- The 51 `biogrid_conflict` pairs are 2.5 % of primary negatives with literature SL evidence. They are a measurable floor on false negatives and are flagged as pre-registered.

## 6. Bias and composition (universe, primary tier: 104 pos / 2,049 neg, prevalence 0.0483)
- **No hub domination:** 104 positives involve 204 genes, and no gene has more than 2 positive pairs. The pairs touching the top-5 genes are 7.7 % of positives (strict 7.4 %, lenient 8.2 %).
- **Ribosomal pairs:** 24 of 104 primary positives (23 %) are ribosomal-protein pairs (RPL/RPS, from strain-table allele names; review only), against 7 of 2,049 negatives. In the biogrid tier the figure is 42 of 227. The DESIGN's pre-declared "with/without ribosomal" analysis is essential for yeast.
- **Fitness confound remains despite the floor:** median min-SMF is 0.877 for positives vs 0.990 for negatives. At measurement level, the POS-rule rate is 5.2 % at E 0.5–0.6 vs 0.03 % at E > 0.85.
- **Composition:**
  - 89 of 104 positives are `within_species_paralog` (WGD-type), against 907 of 2,049 negatives.
  - Median min Ensembl identity is 68.5 % for positives vs 17.2 % for negatives.
  - Negatives cluster in large families. Seripauperin-family genes have up to 23 labelled pairs each, so splits must group families, as §6.1 already plans.
- **Ascertainment:** 68 of 104 positives also have BioGRID LT SL records, and only 14 positives come from SGA alone.
- **Disclosure:** while checking composition I computed an identity-only AP once on the unsplit label set, which includes future test pairs. I did not record the value and do not report it. The PI should know that it was seen.

## Issues
| # | Severity | Issue | Evidence | Fix |
|---|---|---|---|---|
| 1 | minor (interpretation) | All SGA pairs have n_T ≤ 2, so the R2-3 rule reduces to "one POS with at most one NEG → pos", and `discordant` can never occur. The tie goes to positive: 7 of 104 universe positives (11,164 of 34,167 genome-wide) are one-POS-one-NEG orientation conflicts. | QC `sga_pairs.max_n_T` = 2, `discordant` = 0, `orientation_conflict_by_primary` | None to labels (pre-registered, D4). Downstream must report the pre-declared metrics without flagged pairs. |
| 2 | minor (QC labelling) | `untestable_sick_singles` (1,152,957 pairs) also contains pairs with no SMF. 679,131 measurements have E undefined, and the strain-table fallback never filled one (`smf_source.strain_table` = 0). | QC `measurements` | Report untestable pairs split by reason (no_smf vs E ≤ 0.5). |
| 3 | minor (provenance) | `git_head_at_run` = ea97950 predates the code commit a79394b, so the labels were produced from uncommitted code. | QC JSON. My HEAD re-run is byte-identical. | Note it in the tree, or re-run from HEAD so the recorded hash points at the code. |
| 4 | minor (cosmetic) | 59 universe rows have an empty `source`: their only evidence is BioGRID SL records with an essential member. | `primary_by_source['']` | Add a source code for them, or drop rows that carry no usable evidence. |
| 5 | info | Stringent and intermediate gaps: the readme excludes `SGA_DAmP.txt` (§3). | `work/stringent*.json` | Add this explanation to the `known_answer_gaps.json` note. Keep the xfails as they are. |
| 6 | info (downstream) | 23 % of positives are ribosomal; positives have lower min-SMF; 160 `untested` rows sit in the labels file. | §4, §6 | Loaders must treat `untested`, `untestable` and `amb` as unlabelled. Report results with/without ribosomal pairs and stratified by min_smf. |
