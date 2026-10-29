# Round 2 review: data curation (genetic-interaction labels)

**Reviewer role:** data curator for combinatorial CRISPR screens, yeast SGA, and other-species genetic data.

**Read:** REVIEW_BRIEF.md, PROMPT.md, DESIGN.md v2 (§1–§3, §6), reviews/round1/data.md, reviews/round1/PI_response.md, docs/lit/human_labels.md and docs/lit/yeast_and_other_species_labels.md.

**Data checked:** the files that the acquisition agent had placed in `data/raw/` by 2026-10-28. There is no `README_ACQUISITION.md`, so I worked from the per-source `SOURCE.json` files.

**Scope of my own checks:** everything below came from local files with small read-only scripts. I downloaded nothing. All numbers are measured facts about the source tables, used only to test whether the rules can be implemented; they are not proposals to tune against.
- **Ito and CHyMErA trial counts.** I applied the rules proposed in §2 once, as fixed in advance, and changed nothing after seeing the counts. The PI should do the same: freeze the thresholds, then run.
- **Scratch copy.** My only scratch file is an extracted copy of the Costanzo strain/SMF sheet at `<scratchpad>`.

---

## 1. Overall verdict

v2 adopts every round-1 must-item in principle: one SL definition, essential-member pairs untestable, the human robust tier, Flister kept out, gene as the unit, a caller-independent universe, and per-source rules. Some of these are now implementable from the files on disk. Others are not yet: "DKO criterion" and "borderline margin" are still unspecified per study, and two rules as written would bake gene-level single-mutant fitness into the labels.

The main problems:
1. **The yeast negative rule (DMF > 0.5 in all measurements) has no testability floor.** It makes pairs with sick singles positive-or-ambiguous by construction. Measured in the Costanzo NxN file: **33% of the measurements that pass the v2 positive rule come from the stratum where a negative call is effectively impossible.**
2. **The human absolute-DKO criterion cannot be computed from most authors' tables.** Dede S2, Thompson SD1 and Chou T1 hold no DKO or single-gene LFC. A harmonised single + DKO LFC table for Dede, Thompson, Parrish, Ito and CHyMErA now exists on disk (in4mer `Supp_table1.xlsx`), plus Chou's gene-level LFCs for in4mer. One uniform DKO + testability rule can therefore be applied to every study except Harle/Chan/Klingbeil.
3. **"All tested pairs from the sources" is not a usable label universe for yeast.** SGA tested about 23 M pairs.

All of this is fixable at M1 with no new sources. §2 gives exact rules.

### Status of round-1 must-items

| Round-1 item | v2 status | Remaining gap |
|---|---|---|
| D1 null alleles only (yeast) | Adopted | `_S` (suppressor-carrying) query strains are neither TS nor DAmP but are not clean nulls either. They must be excluded too (R2-2). |
| D2 one SL definition | Adopted in words | Yeast: no testability floor (R2-1). Human: DKO criterion unimplementable from the authors' tables (R2-5). |
| D3 aggregation | Adopted | The borderline margin is unspecified, and "no borderline line" penalises many-line libraries (R2-7). Same-line merging is unspecified. |
| D4 Flister out | Adopted | — |
| D5 gene unit | Adopted | Implementable. Exact mapping rules needed (R2-8). |
| D6 caller-independent universe | Adopted in words | "All tested pairs" is wrong for yeast and includes human non-paralog pairs (R2-4). |
| D7 viable ≠ negative | Adopted | Dean has no p-values and uses a growth-rate scale. Needs its own rule (R2-9). |
| D8 BioGRID-LT | Adopted | — |
| D9 per-source rules | Deferred to `config/labels.yaml` | Concrete per-study rules in R2-6, partly verified (Thompson reproduces SD5 exactly). |
| D12 unseen species | Adopted | Lloyd & Meinke single-mutant classes are not on disk, and no fallback is pre-registered (R2-13). |

The round-1 should-items D10 (measurement multiplicity), D11 (linked yeast pairs) and D16 (Dean as independent check) do not appear in v2 or in the PI response. I re-raise them as R2-3, R2-11 and R2-9.

---

## 2. Concrete changes

### R2-1 (must). Yeast: add a testability floor on the expected double-mutant fitness, and give the negative rule a margin.

**Answer to the PI's question.** Requiring DMF > 0.5 for negatives is sensible **only if** pairs whose *expected* DMF is already ≤ 0.5 are removed as untestable. Without that, the label becomes a function of single-mutant fitness.

Under the multiplicative model, DMF ≈ E + ε with E = SMF_query × SMF_array.
- **E ≤ 0.5.** The negative rule (DMF > 0.5 everywhere) can only pass through a positive interaction or noise. Any ε < −0.08 with p < 0.05 is automatically a positive. In this stratum the label is "positive or ambiguous" by construction.
- **E ≈ 1.** DMF ≤ 0.5 requires ε ≲ −0.5. The effective GI threshold therefore depends on the single-mutant fitness (SMF).

Measured in `SGA_NxN.txt` (deletion query × DMA rows only; 12,698,939 rows):
- 943,251 rows have a query SMF of NaN, so E is undefined.
- Of the 11,755,688 rows with E defined, 528,756 (4.5%) have E ≤ 0.5. A further 326,057 have 0.5 < E < 0.6 and 528,122 have 0.6 ≤ E < 0.7.
- **Measurements passing the v2 positive rule:** 36,556 with E > 0.5, plus 17,956 with E ≤ 0.5 (where ε < −0.08 already implies DMF ≲ 0.42). So **one third** of v2-positive measurements come from a stratum where a v2 negative is essentially impossible.
- **Consequence for the model:** features that track "this gene is sick when deleted" (conservation, family size) would be rewarded for predicting fitness rather than SL.

**Rule** (per measurement *m*, for null × null alleles only):

```yaml
yeast_sga:
  E_m: query_SMF * array_SMF            # use the row's own SMF columns (temperature-matched); if NaN, take the
                                        # strain's SMF at the row's temperature from strain_ids_and_single_mutant_fitness.xlsx;
                                        # if still NaN -> measurement untestable
  testable_m: E_m > 0.5                 # logical minimum: a non-interacting pair must not already meet the DMF cut
  POS_m: eps < -0.08 and p < 0.05 and DMF <= 0.5
  NEG_m: (eps >= -0.04 or (p >= 0.05 and eps >= -0.08)) and DMF > 0.5
  MID_m: otherwise (testable but neither)   # gives an epsilon margin of 0.04 between NEG and POS
  untestable_sick_singles: pair has >=1 null x null measurement but none with E_m > 0.5  -> own stratum, counted and reported
```

- Report positive rate by E bin (≤0.5 / 0.5–0.6 / 0.6–0.7 / 0.7–0.85 / >0.85) in `results/labels/`.
- This keeps v2's thresholds. It only adds the floor and the ε margin; the margin is the yeast analogue of the human "borderline" band.
- **A stricter alternative,** if the PI prefers full margin on both axes: testable iff E ≥ 0.7, and NEG additionally needs DMF > 0.6. That is log-symmetric: −0.5 log2 is half of the −1 log2 SL depth. It costs about another 7% of measurements. Pre-register one alternative; do not try both.

### R2-2 (must). Yeast allele whitelist: decide by strain-ID suffix, and drop suppressor strains.

- **Checked suffix classes.** In the Costanzo strain table (10,859 strains) the classes are `sn` 4,071, `dma` 3,893, `tsq` 1,167, `damp` 819, `tsa` 792 and **`S` 117**.
- **`_S` strains carry suppressors.** For example `YBR126C_S26`, allele `tps1-supp2`. They contribute 264,120 rows to NxN and also appear in ExN_NxE. They are not clean null alleles: exclude them from every tier.
- **`Arraytype/Temp` does not identify the allele.** The DAmP file contains DMA strains on `TSA26` arrays, e.g. `YBL007C_dma91 sla1 TSA26`.

**Primary-tier whitelist:**
- query suffix ∈ {`sn`} × array suffix ∈ {`dma`}, plus the reverse orientation;
- the files are `SGA_NxN.txt` and, for the null × null rows only, `SGA_ExN_NxE.txt` (`S|dma`/`sn|dma` rows; by my count the latter has no `sn|dma` rows);
- both temperatures, DMA26 and DMA30, kept as separate measurements.

**Gene mapping:** take the gene from the strain table's `Systematic gene name` column, not from string parsing. Assert the class counts above in a parser test.

**Other checks:** DMF reaches negative values (min −0.1525), so the floor is ≈0. DMF ≤ 0.5 is far from the floor, which is good. The strict DMF ≤ 0.2 tier is reachable (29,174 NxN null × null measurements).

### R2-3 (must). Yeast pair-level aggregation that does not reward measurement count (round-1 D10).

"Positive in some orientation" combined with "negative only if DMF > 0.5 in *all* measurements" pushes heavily measured pairs toward positive. Measurements multiply because of two orientations, two temperatures, multiple `sn`/`dma` strains per gene, and 2010/2016 strains.

```yaml
yeast_pair:
  n_T: number of testable null x null measurements (both orientations, both temperatures)
  positive: n_POS >= 1 and n_POS >= n_NEG
  discordant: n_POS >= 1 and n_NEG > n_POS          # ambiguous, flagged, reported
  negative: n_POS == 0 and n_NEG >= 1 and n_MID <= floor(0.25 * n_T)
  ambiguous: otherwise
```

- Report label rate against n_T (round-1 D10 audit).
- State in DATA.md that the raw file already silently dropped reciprocal pairs with opposite-sign significant scores (Costanzo readme).

### R2-4 (must). Define the label universe explicitly; "all tested pairs" is wrong for yeast.

DESIGN §3 says "Label universe is independent of the caller (all tested pairs from the sources)".
- **Yeast:** SGA tested about 23 M pairs, almost all non-paralogs. Taken literally, this makes the yeast label set a genome-wide GI set.
- **Human:** it brings in non-paralog library pairs: Thompson's 447 "Bioinformatically_derived_pair", 95 "SynLethDB_gene_pair" and 4 "Miscellaneous" pairs (verified in SD1), plus Harle/Chan non-paralog pairs (`sgrna_group` ≠ "Paralog" in Harle S4).

In the end-to-end metric (§6.4) these would count as caller-missed positives at minimum score.

**Rule:**
- `U_species` = within-species paralog pairs from **one pinned Ensembl Compara release**: the homology TSV for *S. cerevisiae* and for *H. sapiens*, mapped to gene IDs via R2-8.
- Compara is annotation used only to *select which pairs carry labels*, so this respects D1. It is the same choice Dennler & Ryan made with Ensembl 111, which keeps the cross-check comparable.
- Labels are attached only to pairs in U.
- Keep the authors' class columns (Thompson `Gene pair class`, Harle `sgrna_group`) as a logged cross-check, not as the rule.
- Add the Compara file to DATA.md. It is not on disk yet.
- Caller recall (§3) is then reported against U.
- Pairs in U that SGA or the screens never tested are **unlabelled**, not negatives.

### R2-5 (must). Human: one uniform, implementable DKO criterion and testability rule from harmonised LFCs.

What the tables actually contain (checked):

| Study | Authors' table on disk | Has single-gene LFC? | Has DKO LFC? | Has a significance call? |
|---|---|---|---|---|
| Dede | S2 = zdLFC only | no | no | zdLFC < −3 (reproduced: 17/21/19 per line, 24 pairs) |
| Thompson | SD1 = FDRs, per-gene MAGeCK/BAGEL calls, "Passes filter?" | no | no | yes (R2-6) |
| Parrish | Table S4 not obtained (PMC bot wall); GEO counts on disk | — | — | GI/FDR only in S4 |
| Ito | Tables 4/5 (raw/inferred LFC; singles = AAVS1-paired rows), 8–10 (GEMINI score/p/FDR) | yes | yes | yes |
| CHyMErA | S8 per line × timepoint: `Paralog1/2.Mean`, `Pair.Mean`, per-orientation FDR and `gi` flags | yes | yes | yes |
| in4mer | Chou T1 = dLFC/ZdLFC/RdLFC; **Chou `In4mer_4screens_gene_level_log_fold_change_cleandup.txt` (19,490 singles + 5,120 combos)** | yes | yes | ZdLFC |
| **in4mer `Supp_table1.xlsx` (Hart lab, CC BY 4.0)** | 21 sheets: Dede ×3, CHyMErA ×2, Parrish ×2, Thompson ×3, Ito ×11. Columns `gene1_mean, gene2_mean, exp_mean, obs_mean, CohenD, deltaMean` | **yes** | **yes** | effect size only |
| Harle / Chan | S4 / Suppl. T6: `mean_norm_gi`, `fdr`, `is_bassik_hit`, expression flags; raw guide counts incl. safe-targeting singles (S3 / Suppl. T5) | from counts | from counts | yes |

**Rule** (per study × line × the timepoint of R2-6):

```yaml
human_line_call:
  lfc_source: in4mer Supp_table1 (Dede, Thompson, Parrish, Ito, CHyMErA); Chou gene-level LFC (in4mer);
              Harle/Chan: recomputed from S3/T5 counts with the authors' code (Harle repo is MIT) or, if not done, authors' calls only
  m: median of all single-gene LFCs in that screen (gene1_mean/gene2_mean, deduplicated per gene)
  s: 1.4826 * MAD of the same
  T: max(2 * s, 0.5)                              # SL depth in log2 units below m; 0.5 floor protects tight screens
  a, b, ab: gene1_mean - m, gene2_mean - m, obs_mean - m
  testable: (a + b) > -T and not study_flag_untestable    # expected DKO must not already reach SL depth (same logic as yeast E > 0.5)
  SL: testable and GI_strict and ab <= -T
  clear_negative: testable and not SL and (not GI_relaxed or ab > -T/2)
  borderline: testable and not SL and not clear_negative
```

**Why this rule:**
- It is label-free and screen-scaled.
- It needs no reference essential gene list.
- It is the human counterpart of R2-1.
- It replaces "midpoint of essential/non-essential distributions", which no paralog-library table supports.

**Ito trial** (raw Table 4 LFC, rule fixed in advance, run once):
- T is 0.50–0.52 in every line, and 201–508 of 4,859 pairs per line are untestable.
- SL calls: 1,714 pair × line.
- Pair level, Ito alone: 210 positives (≥2 SL lines and penetrance ≥ 0.25).

Thompson illustrates why testability matters. `AARS2:AARS` in A375 has `gene2_mean` −2.28 and `exp_mean` −3.71. A pair like this has a positive dLFC and would otherwise become a clean "negative".

Studies whose testability cannot be determined (no LFCs, no flag) **do not enter the primary tier**; they go to the lenient tier only. Otherwise their pairs would be positive-or-ambiguous by construction, which is the library-membership confound from round-1 D3.

### R2-6 (must). Per-study rules for `config/labels.yaml`

Proposed values; the verified parts are marked [V].

| Study | Lines / timepoint | Exclusions | GI_strict | GI_relaxed (for borderline) | Untestable flag | Notes |
|---|---|---|---|---|---|---|
| Dede 2020 | A549, HT29, OVCAR8 | — | zdLFC < −3 [V: 17/21/19] | zdLFC < −2 | none; R2-5 testability from Supp_table1 | |
| Thompson 2021 | A375, MeWo, RPE-1; **D28** | non-paralog classes (R2-4) | `{L}_D28_rra_fdr_low` < 0.1 **and** `{L}_D28_t_fdr_low` < 0.1 **and** filter == "Yes" [V: reproduces SD5 exactly, 44/40/57] | D28 either FDR < 0.25, or D14 rra & t < 0.1 | `Passes {L} filter?` == "No" [V: 549/431/388 pairs; not a pure function of the BAGEL/MAGeCK columns, so use it as given] | Without the filter, rra&t gives 199/177/201 calls, about 4× the hits |
| Parrish 2021 | PC9, HeLa | — | GI < −0.5 and FDR < 0.1 (needs S4) | GI < −0.3 or FDR < 0.25 | R2-5 | **If S4 is not obtained by M1:** pre-register a fallback of Chou-style ZdLFC < −2 from Supp_table1 `deltaMean` (Z per screen as in Chou's `Figure02f-i_Step01` notebook) |
| Ito 2021 | 11 lines (Tables 4/8–10) | **Tables 6/7/11–13 (MELJUSO DMSO/ERKi anchor screen)** | GEMINI FDR (Table 10) < 0.05 | FDR < 0.25 | R2-5 | **Drop the "GEMINI ≥ 90th percentile within screen" rule.** It fixes about 10% hits per line whatever the screen quality; the absolute R2-5 criterion already removes Ito's inflation (FDR < 0.05 alone gives 325–610 per line; with the DKO criterion 48–287). Use all 11 lines unless a pre-stated QC reason excludes one (Kebabci used 10; say which and why) |
| CHyMErA 2020 | HAP1, RPE-1 | **Torin sheets and `HAP1 +/- Torin` column**; the 3–4 duplicated gene-pair rows per sheet (if duplicates disagree → line borderline) | `early_both_gi` or `late_both_gi` **and** the corresponding `*_effect_size` < 0 (the `gi` flags also mark *positive* GIs: 1–11 `both_gi` rows per sheet have effect > 0 [V]) | one orientation significant with effect < 0 | R2-5 | DKO at the **late** sheet (`HAP1.T18`, `RPE1.T24`). Do **not** use the `*_expected_LFC*` columns as single-gene LFCs: they are not on the single-LFC scale. Trial: HAP1 late 41 SL, RPE1 late 8 |
| in4mer 2024 / Chou 2025 | K562, A549 (prototype); MelJuso, A375 (Inzolia) | **`MelJuso_sel_T21` (selumetinib)** | ZdLFC < −2 (Chou T1) | ZdLFC < −1 | R2-5 from the gene-level LFC file | Prototype and Inzolia are separate screens; a pair in both counts once per line |
| Harle 2025 | 27 lines (S4) | `sgrna_group` ≠ Paralog only via R2-4 | `mean_norm_gi` < −0.5 and `fdr` < 0.01 and `is_bassik_hit` == 1 | `mean_norm_gi` < −0.3 or `fdr` < 0.1 | authors' single-depletion columns in S4 (confirm column names; also `target*__is_expressed`) | Known answer 117 SL / 355 never |
| Chan 2025 | 10 uveal lines (Suppl. T6) | as Harle | as Harle | as Harle | as Harle | Known answer 105 pairs SL in ≥1 line |
| Klingbeil 2024 | 22 lines | — | GEMINI FDR < 0.05 + R2-5 | FDR < 0.25 | R2-5 | **Not on disk** (`klingbeil2024/` holds only SOURCE.json). Domain collapse: a line is SL iff >50% of the pair's domain-combination rows are SL, borderline if 1–50%. Exclude from primary if not obtained by M1 |
| Chou rescoring of Dede/Parrish/Thompson | — | primary excluded (v2) | — | — | — | Sensitivity only |

### R2-7 (must). Pair aggregation: a borderline allowance, same-line merging, and line conflicts.

**(a) Borderline allowance.** v2's negative rule is "no line within the borderline margin". With many lines, almost every pair has some borderline line, so negatives are drawn mostly from few-line libraries. That is the round-1 D3 confound in reverse.
- In the Ito trial, borderline line calls are 384–831 per line.
- Pair-level Ito negatives: 1,796 with no borderline allowed vs 2,534 with ≤ floor(0.1·n_testable) borderline lines.
- **Rule:** negative = n_testable ≥ 2, n_SL = 0, and n_borderline ≤ floor(0.1·n_testable). This means 0 borderline lines allowed for <10 lines, 1 for Ito, 2 for Harle.

**(b) Same line across studies.** Normalise cell-line names to Cellosaurus IDs (Harle S2 already gives them). Overlaps: A549 (Dede, Ito, in4mer, Harle), MeWo (Thompson, Ito, Harle), A375 (Thompson, in4mer, Harle), MelJuso (Ito, in4mer), RPE-1 (Thompson, CHyMErA), MEL202 (Ito, Chan). The line-level call is:
- SL if SL in ≥1 study and no study is clear_negative for that line;
- clear_negative if all testable studies are clear_negative;
- `line_discordant` otherwise, which counts as borderline.

**(c) Counting.** n_lines counts distinct Cellosaurus IDs. n_studies counts distinct publications, with in4mer prototype and Inzolia as one study.

**(d) Audit.** Report positive rate by source study and by n_testable bin (round-1 D3 audit) before freezing.

### R2-8 (must). Per-gene expansion: exact mapping rules (implementable now).

**Yeast (UP000002311, 2026_03):**
- `gene2acc` is keyed by SGD ID, not ORF name.
- Map Costanzo systematic names via `idmapping` `Gene_OrderedLocusName`.
- Five accessions carry two ORFs: P02309, P02994, P10081, P32324, P61830.
- Use the ORF name as the gene ID. The copies are identical, so which SGD ID goes with which ORF does not matter.

**Human (UP000005640):**
- `gene2acc` has 172,773 rows (it includes non-reference accessions); 20,803 rows fall on the 20,652 FASTA accessions.
- Of these, 89 accessions encode ≥2 HGNC genes, giving 151 extra gene records. Histone H4, P62805, maps to 14 HGNC IDs and 14 ENSG.
- 110 FASTA accessions have no HGNC ID (gene column `-`).
- 2 HGNC IDs map to 2 reference accessions.

**Rules:**
- gene ID = HGNC ID, falling back to the accession when there is none.
- Restrict to accessions in the reference FASTA.
- An HGNC ID with 2 reference accessions keeps the one whose gene2acc third column equals itself; otherwise the longer sequence; log it.
- Screen symbols → HGNC via the pinned HGNC complete set, including previous and alias symbols; drop ambiguous ones. **The HGNC file is not on disk;** add it to DATA.md.
- Write `results/labels/mapping_report.tsv` with losses per source and step.

### R2-9 (should). Dean 2008 and Kuzmin 2020 need their own rule; keep Dean independent for one QC run (round-1 D16).

**Dean** gives only point fitness values: `ORF1_fit ORF2_fit ORF1ORF2_fit`, growth-rate based, with no SE or p-value. The v2 sentence "ε and DMF … meet the same rule" therefore cannot be applied, because p < 0.05 cannot be computed and the scale is not colony size.

```yaml
dean2008:
  testable: W_A * W_B > 0.5
  positive: W_AB == 0              # tetrad-lethal
  negative: W_AB - W_A*W_B >= -0.05
  ambiguous: otherwise
  parser_assert: 201 measured pairs / 49 lethal, as reported in the paper (the survey parse gave 197/48; fail loudly)
```

**Kuzmin** Table S5 (manual download pending): apply R2-1 with p-values if S5 has them; if not, use Dean-style point rules.

**Source conflicts:**
- Dean lethal + SGA NEG → positive, flagged `dean_override` (tetrad dissection is the direct test).
- Dean negative + SGA POS → ambiguous.

**QC:** build the yeast labels once without Dean and report the SGA-call vs Dean agreement table before merging.

### R2-10 (should). Make the yeast tiers nested.

- **Strict:** v2's strict tier is "DMF ≤ 0.2 **or** ε < −0.35". With E > 0.85, ε < −0.35 can coexist with DMF > 0.5, so a pair can be strict-positive and primary-ambiguous. Define strict = primary positive ∧ (DMF ≤ 0.2 or Dean W_AB = 0).
- **Dennler-comparable tier:** keep "ε < −0.35, p < 0.05" as a separate, named tier.
- **Lenient** (ε < −0.12, p < 0.05, no DMF rule): use the same testability floor and negatives as primary. With E > 0.5 and DMF ≤ 0.5, every primary positive has ε < −0.08. It is not automatically < −0.12, so state that lenient ⊉ primary, or define lenient as primary ∪ {ε < −0.12, p < 0.05}.

### R2-11 (should). Linked yeast pairs (round-1 D11, still missing).

- Same-chromosome pairs within a pre-set window (e.g. 50 kb; verify against the Baryshnikova 2010 linkage filter) are flagged `sga_linked`.
- Under SGA they are untestable unless Dean or BioGRID-LT gives a tetrad call.
- Coordinates need a pinned SGD GFF (not on disk; label construction only).
- Report how many U-pairs are lost. Tandem duplicates are the highest-identity pairs, so this interacts directly with the B1 baseline.

### R2-12 (should). Known-answer tests for the label parsers (several verified now).

| Test | Expected value | Status |
|---|---|---|
| Thompson D28 rra&t<0.1 & filter | 44/40/57 equal to the SD5 lists | [V] |
| Dede zdLFC < −3 | 17/21/19 per line, 24 pairs, 403 pairs total | [V] |
| CHyMErA S8 per-sheet rows | 692/692/691/691 rows; duplicated pairs 4/4/3/3 | [V] |
| Ito | 4,859 GEMINI pairs; 5,065 DKO and 3,284 single rows in Table 4 | [V] |
| Costanzo strain classes | sn 4,071 / dma 3,893 / tsq 1,167 / damp 819 / tsa 792 / S 117 | [V] |
| Harle | 117 SL-in-≥1-line / 355 never | — |
| Chan | 105 | — |
| Cusack | 51/110 and 190/110 | — |
| Dean | 201/49 | — |

### R2-13 (should). Acquisition status and pre-registered fallbacks (decide before label code runs).

**On disk:** Costanzo (thecellmap zip, sha256 recorded), BioGRID 5.0.262, Dede, Thompson, Ito, CHyMErA, in4mer figshare (Supp_table1–9), Chou, Harle, Chan, Flister preprint, Dennler, Cusack, Tischler (.doc), PomBase interactions, UniProt 2026_03 for six species.

**Missing:**

| Missing source | Fallback |
|---|---|
| Parrish S4 (manual) | ZdLFC fallback, R2-6 |
| Kuzmin S5 (manual) | Drop Kuzmin; SGA already covers most WGD pairs |
| Klingbeil | Exclude from primary |
| Ryan 2012 pombe E-MAP (manual) | Pombe not an AUPRC set (v2 already conditional) |
| Lloyd & Meinke 2012 single-mutant classes | Report Arabidopsis RD4 with all 110 negatives, flagged as "single-mutant lethality unchecked" |
| Ensembl Compara | — (needed for R2-4) |
| HGNC complete set | — (needed for R2-8) |
| SGD GFF | — (needed for R2-11) |

Also: Costanzo was fetched from thecellmap, not from the BioGRID mirror that round 1 suggested as the MIT source of record. Record the license basis in DATA.md, and if possible confirm that the BioGRID mirror has the same checksum.

### R2-14 (could). Extra label audits.

- Positive rate by E bin (yeast) and by (a + b)/T bin (human).
- Label rate vs number of measurements (yeast) and vs n_testable lines (human).
- These are the direct checks that R2-1, R2-3, R2-5 and R2-7 removed the confounds. Commit them with the labels.

### R2-15 (could). Add a penetrance floor to the "≥2 studies" route.

"SL in ≥2 lines and ≥2 studies" lets a pair that is SL in 1/11 Ito lines plus 1/27 Harle lines (penetrance 0.05) become a robust positive. Add penetrance ≥ 0.1 to that route.

---

## 3. Risks the design does not address

1. **Harmonised LFCs are a re-analysis.** in4mer `Supp_table1` was produced by the Hart lab from counts, so its timepoint and normalisation choices are theirs. Thompson's timepoint in that sheet is unverified.
   - Mitigation: use it only for the DKO depth and testability, keep the authors' significance calls, and record the provenance per line.
   - Check that `exp_mean` = `gene1_mean` + `gene2_mean` holds in every sheet (it does in the rows I saw).
2. **Many human screens have a non-zero single-gene median.** Ito ≈ +0.13 to +0.30; CHyMErA ≈ −0.12 to +0.23. Any rule written on raw LFC (e.g. "obs ≤ −1") instead of LFC − m will be line-biased. R2-5 handles this; keep it in the parser test.
3. **CHyMErA timepoint sensitivity is large.** RPE-1 has 111 `both_gi` calls early vs 54 late. Pre-register one rule; the "either timepoint GI, late-sheet DKO" rule in R2-6 is the least arbitrary option.
4. **Ito's DKO table has 206 pairs without GEMINI scores** (5,065 vs 4,859). Treat them as untested, not negative.
5. **Expression flags.** Harle/Chan testability uses expression flags (`target*__is_expressed`). This is fine for label construction, but those columns must never reach the feature registry. Add them to the forbidden-column list used by the quarantine/whitelist reader.
6. **Prevalence shift between species becomes larger.** With the floor, yeast negatives will be cleaner but rarer among sick-single genes. Human negatives become more numerous with the borderline allowance. Report prevalence per species next to every AP.

## 4. Fine as it is

- The single SL definition (significant negative GI ∧ strong DKO defect ∧ both singles testable) is the right target. R2-1 and R2-5 only make "testable" operational in the same way for both species.
- The null-allele-only primary yeast tier; TS/DAmP kept as an exploratory stratum.
- Flister outside the primary tier. Chou rescoring only as sensitivity. No SLKB/SynLethDB. Dennler used only for a report-only cross-check against our lenient tier.
- Gene, not accession, as the unit. It is implementable now from the on-disk UniProt files (R2-8).
- BioGRID-LT as a separate tier with and without BioGRID-only positives.
- Commit policy (derived calls only; non-open sources flagged) and quarantine of files with forbidden columns.
- The pre-registered thresholds, not calibrated on Dean W = 0. With R2-1 the DMF floor question is moot: the observed minimum DMF is ≈ −0.15, and 0.5 is far from the floor.
