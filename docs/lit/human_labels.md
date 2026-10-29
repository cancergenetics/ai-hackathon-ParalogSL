# Human paralog synthetic-lethality (SL) label sources — survey and labelling recommendation

Survey date: 2026-10-27. Author: label-survey agent (for the PI). Scope: human paralog SL labels for training and testing a
sequence-only paralog SL scorer (metric: AUPRC). Yeast labels are covered elsewhere.

**Verification legend.** [V] = verified directly in this survey: I read the paper/PMC/EPMC full text, the GitHub/Zenodo/figshare
API listing, or I inspected the file. [P] = taken from a secondary source (abstract, search snippet, another paper's
description), not verified at the primary source. [U] = unverified or unknown. Pair and hit counts are as reported by authors
unless I recomputed them from an inspected file (marked "recomputed").

---

## 1. Key take-aways

1. **There is no gold standard.** Every comparison paper says so (De Kegel 2021; Esmaeili Anvar 2024; Chou 2025; Ajmal 2025).
   Agreement across screens is modest. Pairwise screen-to-screen hit overlap is about 31% (De Kegel 2021 notebook [V]). Of 388
   pairs that were hits in one of five early studies, only 26 were hits in more than one study, and 13 form a "candidate gold
   standard" (Esmaeili Anvar 2024 SD1 [V]). The best achievable cross-experiment agreement in the same cell line is roughly
   precision 0.40 / recall 0.32 (Kebabci 2026 [V]). This sets the ceiling on label quality, and any AUPRC we report must be read
   against it.
2. **Tested negatives exist** for most combinatorial screens: per-pair scores for all screened pairs are published. This is the
   main reason to build labels from the **primary screen tables**, not from SL databases. SynLethDB holds mostly positives, and
   SLKB's negatives are defined by a relative top-10% rule.
3. **Screened pairs are a strongly biased sample of paralog space.** Median identity is 54% for screened pairs vs 32% for all
   Ensembl pairs; 89% of screened pairs are the closest paralog and 73% are WGD, vs 13% and 19% overall (De Kegel 2021 [V]).
   Our test distribution will not resemble "all paralog pairs in a proteome". The only near-unbiased human sample is the
   Flister 2025 all-paralog screen (36,648 pairs, 2 cell lines; full per-pair table not yet verified).
4. **SL is highly context-specific.** About 68–86% of hits are found in a single cell line (De Kegel 2021). Only about 4 pairs
   are SL in ≥25 of 27 lines (Harle 2025). A sequence-only, context-free tool can only target **robust/recurrent SL**. The
   primary label should therefore require recurrence, and single-line hits should be treated as ambiguous rather than positive.
5. **Ready-made CC-BY compilations exist.** Two are especially useful:
   - Ryan lab "SL-lenient" human set: 5,316 pairs, 403 SL, 5 screens, Ensembl 111 IDs, CC BY 4.0, Zenodo 14973633 (Dennler & Ryan 2025).
   - Hart lab ZdLFC rescoring of 4 in4mer screens plus a 12-screen hit matrix, CC BY 4.0 (Chou 2025).

   These are good cross-checks, but we should rebuild labels from primary tables so that every label carries its evidence counts.

---

## 2. Source table — primary combinatorial (double-knockout) screens

"Neg avail." = are per-pair results for non-hits published, so that tested-negative pairs can be derived?

| # | Source (DOI) | Platform | Paralog pairs screened | Cell lines | Scoring / SL call (authors) | Hits | IDs | Neg avail. | License | Data location (reachability) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Dede et al. 2020 Genome Biol, 10.1186/s13059-020-02173-2 | enCas12a | 403 (400 computational + lit) [V, recomputed 403] | 3: A549, HT29, OVCAR8 | zdLFC (dLFC Z-scored after trimming); SL = zdLFC < −3 | 24 SL pairs (recomputed: 24; per line 17/21/19); 79% recur in ≥2 lines | HGNC symbols (Ensembl 95 paralogs) | **Yes** — Table S2 has zdLFC for all 403 pairs × 3 lines [V, inspected] | CC BY 4.0; BMC data CC0 waiver | Addl. file 3 (13 kB txt) [V 200]; SRA PRJNA664967; figshare 12275642 (code) |
| 2 | Thompson et al. 2021 Nat Commun, 10.1038/s41467-021-21478-9 | Cas9 dual gRNA | 1,191 pairs total, of which 645 are paralog pairs [V, recomputed class counts]; the rest are mutual-exclusivity (447), SynLethDB (95) and other pairs | 3: A375, MeWo, RPE-1 | Bliss residual, LOESS + variance smoothing; t-test + RRA, FDR < 0.1; pairs with a singly essential gene removed (MAGeCK/BAGEL) | 105 SL (72 paralog); 40–57 per line | Symbols + Ensembl IDs (SD3) | **Yes** — SD1 lists FDRs for all 1,191 pairs × 3 lines × 2 timepoints plus essentiality filters [V, inspected] | CC BY 4.0 | Nature ESM MOESM4 (SD1), MOESM6 (SD3), MOESM8 (SD5) [V 200]; raw data at EGA EGAD00001006648 (controlled access) |
| 3 | Parrish et al. 2021 Cell Rep, 10.1016/j.celrep.2021.109597 (pgPEN) | Cas9 paired gRNA | 1,030 (≥50% identity, 2-member families) | 2: PC9, HeLa | GI = observed − expected DKO; SL = GI < −0.5 and FDR < 0.1 | PC9 87, HeLa 70; 122 unique (12%) | Symbols | **Yes** — GI scores for all pairs (Table S4 per text) [P]; counts in GEO | Article CC BY-NC-ND 4.0 [V EPMC] | GEO GSE178179 (pgPEN_counts.txt.gz, 1.3 MB) [V]; Zenodo 5081113 (MIT code); Cell Rep mmc tables [U] |
| 4 | Ito et al. 2021 Nat Genet, 10.1038/s41588-021-00967-z | Cas9 dual (SpCas9/SaCas9) | 5,065 designed; 4,858 unique pairs scored (Kebabci notebook) [V] | 11: MELJUSO, GI1, MEL202, PK1, MEWO, HS944T, IPC298, A549, HSC5, HS936T, PATU8988S | GEMINI; authors use FDR < 0.05 (lenient: 435–705 hits per line per SLKB; Flister maps 1,741 "lethal" pairs) | Many, but threshold-dependent. Kebabci: 958 SL of 41,244 pair×line points after GEMINI ≥90th percentile plus an LFC filter (10 lines, 4,170 pairs) | Symbols (Ensembl 91); needs HGNC remapping | **Yes** — LFC + GEMINI FDR for all pairs × lines (Supp. Tables, incl. S4/S9 per Kebabci) [P] | **Not OA**; no open license on supplementary (Springer Nature copyright) [V EPMC] | Supp. Tables 1–13 xlsx, 6.37 MB [V reachable, not downloaded, >5 MB]; SRA PRJNA745952; GitHub sellerslab/ParalogV1_DUSP46 (no license) |
| 5 | Gonatopoulos-Pournatzis et al. 2020 Nat Biotech (CHyMErA), 10.1038/s41587-020-0437-z | Cas9–Cas12a hybrid | ~672 paralog pairs (688 pairs in Table S8 incl. a few controls) [V, recomputed] | 2: HAP1, RPE-1 (+ HAP1 Torin1) | Observed vs expected LFC, per-orientation FDR; authors' "gi_type". The authors' negative-GI calls include weak effects: 40.6% hit rate raw; De Kegel re-thresholded to hit and LFC ≤ −0.9, giving 10.9% | 72 SL (De Kegel filtered) | Symbols | **Yes** — Table S8 lists effect sizes and FDRs for all pairs [V, inspected] | **Not OA**; no open license | Nature ESM MOESM10 (Table S8, 1.2 MB) [V]; GEO GSE144281; lab site crispr.ccbr.utoronto.ca/chymera returns **404** [V] |
| 6 | Esmaeili Anvar et al. 2024 Nat Commun (in4mer/Inzolia), 10.1038/s41467-024-47795-3 | Cas12a 4-guide arrays | Prototype: 2,082 pairs; Inzolia: 4,435 pairs + triples/quads; 4,540 pairs scored overall (recomputed from Chou T1) [V] | 4: K562, A549 (prototype); MelJuso, A375 (Inzolia) | dLFC < −1 (authors) → superseded by ZdLFC (row 7) | Hit counts not given in the paper; SD1 = 13-pair "gold standard" + 26 multi-study hits + 388 union of hits from 5 earlier studies [V, inspected] | Symbols | **Yes** — LFC tables (SD6–7; figshare) [V listing] | CC BY 4.0 (article + figshare) [V] | figshare 10.6084/m9.figshare.24243832 [V API]; Nature ESM MOESM4 = SD1 (18 kB) [V] |
| 7 | Chou et al. 2025 Genome Biol, 10.1186/s13059-025-03660-0 (re-analysis of 6) | — | 4,540 (K562/A549: 2,070; MelJuso/A375: 4,414) [V, recomputed] | 4 in4mer + re-scored Dede, Parrish, Thompson (12 screens total) | ZdLFC < −2 and ZLFC < −2 (alternatives: RdLFC < −0.7; dLFC < −1 and LFC < −1 gives 16–75 per line) | Supp. T2: 371 hit pairs across 12 screens; 180 hit in ≥2 screens; 80 hit in ≥2 *studies* (recomputed) [V] | Symbols | **Yes** for in4mer (dLFC/ZdLFC/RdLFC for all 4,540 pairs) [V, inspected]. T2 lists **hits only**; there 0 = "not a hit or not tested" (ambiguous) | CC BY 4.0 (article + figshare) [V] | figshare 10.6084/m9.figshare.28195865 [V API]: Supp. Table 1/2, LFC tables, 13-SL reference list |
| 8 | Klingbeil et al. 2024 Cancer Discov, 10.1158/2159-8290.CD-23-1529 | Cas9 dual gRNA (domain-focused) | 2,529 paralog DKOs in text; Fig. 1 says n = 2,726 (inconsistent) [V PMC]; signalling + epigenetic regulators | 22 cancer lines (PDAC, SCLC, AML, …; listed in Table S3) | GEMINI (R v1.4.0); no explicit hit threshold in the paper | Not given; Kebabci uses FDR / GEMINI | Symbols + domain tags (pairs repeat per domain combination) | **Yes** — GEMINI score/FDR + LFC tables for all pairs × 22 lines (S2, S4–S6, S10) [P via Kebabci notebook] | Article CC BY-NC-ND 4.0 [V]; Zenodo 10042504 code CC BY 4.0 [V] | GEO GSE242517 (raw tar, 906 MB) [V]; AACR supp. tables [U, not fetched] |
| 9 | Harle et al. 2025 Genome Biol, 10.1186/s13059-025-03737-w | Cas9 dual guide (Sanger) | 472 pairs, of which 262 are paralogs (≥45% identity, lethal fly/worm ortholog) [V] | 27: 8 melanoma, 10 NSCLC, 9 PDAC | Bliss GI, LOESS, t-test; SL = mean norm. GI < −0.5 and FDR < 0.01 and Bassik-method hit; pairs with a singly essential gene excluded | 117 pairs SL in ≥1 line (882 pair×line hits); 355 never SL (Kebabci notebook) [V]; 4 pairs SL in ≥25/27 lines | Symbols (Ensembl 92) | **Yes** — Table S4/S5: GI, FDR, single-depletion flags for 472 × 27 [V via Kebabci notebook] | CC BY 4.0; code + figshare MIT [V] | Addl. file 1 xlsx (32 MB, not downloaded) [V size]; ENA PRJEB60853; GitHub team113sanger/harle_dgCRISPR_paralogs (MIT) [V]; figshare 25954027 (mostly imaging, MIT) [V] |
| 10 | Chan et al. 2025 Nat Genet, 10.1038/s41588-025-02222-1 | Cas9 dual guide (Sanger, same platform as 9) | 514 pairs (262 paralog; 210 putative SL; 42 uveal-specific) [V PMC] | 10 uveal melanoma lines | Same as 9 (GI < −0.5, FDR < 0.01, single-essential filtered) | 105 pairs SL in ≥1 line; 23 shared by ≥6 lines [V] | Symbols | Probably — Supp. Tables 5–7 [P]; full coverage unverified | CC BY 4.0 [V] | Nature ESM Supp. Tables 1–17 [U, not fetched] |
| 11 | Flister et al. 2025 Cell Rep, 10.1016/j.celrep.2025.116512 (preprint 10.1101/2024.07.16.603642) | enAsCas12a, 4-guide arrays | **36,648** (all Ensembl paralog pairs ≥20% identity) | 2: NCI-H1299, MDA-MB-231 (genome-wide); hybrid screen of ~5,318 (or 5,116; text inconsistent) pairs in 20 lines | zL2FC; SL = median DKO < −1 zL2FC, dzL2FC > 1.64 (genome-wide) or > 2 (hybrid); expression > 2.5 TPM | Genome-wide: 144 + 45 SL (<0.4%); hybrid: 780 SL interactions in ≥1 line (14.6%) | Symbols | **Unverified.** Preprint Table S3 [V, inspected] lists only "Lethal" calls for 35,108 pairs (their 2 lines + mapped Ito/Dede/Parrish/Thompson/CHyMErA/Klingbeil). NaN conflates "not SL" with "not tested". Tables S1/S2/S5 are >5 MB and not inspected; they may hold full scores | Journal: CC BY-NC-ND [V EPMC]; preprint license [U] | bioRxiv supp. media-1…10.xlsx [V reachable]; raw-data accession [U] |
| 12 | Burgold et al. 2025 Nat Commun, 10.1038/s41467-025-67256-9 | Cas9 tRNA dual guide | 100,136 guide pairs (paralogs + other GIs; paralog count [U]) | 1: HT-29 | [U] | Validated CNOT7/8, ASF1A/B | [U] | Probably (Zenodo analysis data) [P] | CC BY 4.0; Zenodo 17191951 CC BY 4.0 [V] | Zenodo 17191951 (89 MB zip); figshare 25533091 (1.15 GB raw) [V API] |

**Not paralog combinatorial screens, but often confused with them.**
- **Köferle et al. 2022** Cell Rep 39:110636 (10.1016/j.celrep.2022.110636; CC BY per EPMC [V]). This is a DepMap/DRIVE
  correlation analysis ("PaCT"): paralog A expression vs paralog B dependency across 108,406 testable pairs, giving 2,040 hits
  (3 SD and p < 0.05). It also includes single-gene validation screens in 7 lines. It is **not** a double-knockout screen. Its
  labels are inferred from cancer-cell-line covariation.
- **Fielden et al. 2025 Nature** (DNA-damage-response GI map) and **Horlbeck 2018 / Han 2017 / Shen 2017 / Najm 2018 /
  Zhao 2018** are combinatorial screens, but not paralog-focused. They contain only a handful of paralog pairs. Not surveyed in depth.

## 3. Inferred (non-combinatorial) label sets, meta-analyses and databases

| Source | What it is | Size | Negatives? | License / location | Use for us |
|---|---|---|---|---|---|
| De Kegel et al. 2021 Cell Syst, 10.1016/j.cels.2021.08.006 | (a) DepMap-inferred "robust SL" set: A2 loss (HD, LOF or no expression) associated with A1 dependency across ~762 lines, 10% FDR. (b) "Screen consensus" set from 4 screens | (a) 3,634 pairs: 126 SL, 3,508 non-SL (p > 0.05); 240 ambiguous dropped [V preprint text]. (b) 545 pairs screened in ≥2 studies: 50 SL in ≥2 studies, 407 never SL [V] | Yes, explicit | Article CC BY (EPMC) [V]; preprint CC BY-NC-ND; code github cancergenetics/paralog_SL_prediction (no license file) [V] | (a) is an orthogonal, inferred label. Keep it as a separate auxiliary tier, never mixed with screen labels. Its selection rule requires frequent loss of one paralog. (b) is the template for our consensus rule |
| Dennler & Ryan 2025 NAR Genom Bioinform 7:lqaf051 (Zenodo 10.5281/zenodo.14973633) | Ready-made human and yeast paralog datasets (Ensembl 111) | Human "SL-lenient": 5,316 pairs, 403 SL (5 screens: Thompson, Dede, Parrish, CHyMErA, Ito; SL if SL in any screen, non-SL if 0 in all) [V notebook]. Human "SL" (DepMap): 23,734 pairs, 116 SL [V notebook] | Yes | **CC BY 4.0** [V Zenodo API]; ens111_human_SL-lenient.csv 12.7 MB (first 30 kB inspected: 134 cols incl. A1/A2 Ensembl + Entrez IDs, `SL`) | Cross-check / baseline replication. Binary only, no per-screen evidence. Includes data-derived feature columns, which we must not use as features |
| Ajmal et al. 2025 NAR Genom Bioinform, 10.1093/nargab/lqaf129 | Benchmark of 5 GI scoring methods (Gemini-Sensitive, Parrish score, Orthrus, zdLFC, …) on Dede, CHyMErA, Ito, Parrish, Thompson (21 screens) vs De Kegel (3,634 pairs, 3.5% SL) and Köferle (48,095 pairs, 0.92% SL) benchmarks | — | — | CC BY 4.0 [V]; code github cancergenetics/Benchmarking-GI-Scores (license NOASSERTION) [V] | Gemini-Sensitive and the Parrish score do best overall; method choice changes hits. Supports rescoring as a sensitivity analysis |
| Kebabci et al. 2026 bioRxiv, 10.64898/2026.01.19.700065 | Cell-line-specific SL classifier trained on Ito (GEMINI + LFC dual threshold), tested on Klingbeil, Parrish, Harle | Training: 958 SL / 40,286 non-SL pair×line points | Yes | CC BY 4.0 preprint; figshare 31058182 (predictions, CC BY 4.0); code MIT [V] | Labelling recipe (dual threshold: GEMINI ≥ 90th pct within screen AND DKO LFC below a screen-specific essential/non-essential midpoint). Sequence identity alone: AP ≈ 0.06 vs 0.04 chance on Klingbeil, so the baseline is weak on screened pairs |
| Esmaeili Anvar 2024 SD1; Chou 2025 reference list | 13-pair "candidate gold standard" (CNOT7/8, PITPNA/B, TIA1/TIAL1, SAR1A/B, PTP4A1/2, GSK3A/B, CSNK2A1/2, CSNK1D/E, MAPK1/3, ARFGEF1/2, HDAC1/2, ASF1A/B, SLC25A28/37) [V file] | 13 | No | CC BY 4.0 | Positive controls / sanity checks only |
| BaCoN (Rohde et al. 2025 Mol Syst Biol, 10.1038/s44320-025-00103-7) | Expression-vs-fitness buffering network from DepMap | 808 high-confidence buffering pairs (899 in preprint) [P] | No explicit negatives | [U] | Optional inferred tier; not a screen |
| SLKB (Gökbağ et al. 2024 NAR, 10.1093/nar/gkad806) | 11 CDKO studies, 22 lines, rescored with 5 methods | 16,059 SL / 264,424 non-SL as reported; 13,173 / 248,785 by majority vote | Yes, but "SL" = top 10% per method (relative) | Article CC BY-NC; **data GPL-3.0+** (figshare 22902839, 241 MB SQL dump) [V API]; site reachable, last updated 05/2023 | Useful index of screens. Do **not** derive committed labels from it (GPL, relative thresholds, mostly non-paralog). Use primary tables instead |
| SynLethDB 2.0 (Wang et al. 2022 Database, 10.1093/database/baac030) | Curated + predicted + text-mined SL | 35,943 human SL; few non-SL (5,798 across species) | Few | Article CC BY-NC; DB license not stated [V] | Not suitable: mixed evidence including predictions, and few tested negatives |
| BioGRID (Ito publication 240353) | Curated GIs from screens | 687 interactions for Ito [V] | No | MIT | Not needed; primary tables are better |
| "SLorth" | Not found under that name [V search]; Syn-Lethality / SynLethDB use yeast→human ortholog mapping | — | — | — | Do not use (cross-species leakage risk) |

## 4. Recommended human labelling scheme

### 4.1 Unit and identifiers
- **Unit:** unordered gene pair. Map symbols to HGNC ID → Ensembl gene, using HGNC previous/alias symbols (several screens use
  outdated symbols, e.g. PALM2AKAP2), then to one canonical protein per gene, consistent with our FASTA input.
- **Restriction:** keep only pairs that are paralogs under our own sequence-based paralog definition. Non-paralog pairs from
  mixed libraries (Thompson, Harle, Chan) are excluded. Record pairs that fail mapping.

### 4.2 Per screen × cell line call (from primary tables; authors' scores, harmonised thresholds)
Each pair×line observation gets one of three states:
- **SL:** the authors' significance criterion **and** an effect-size criterion. The double knockout must cause a real fitness
  defect (DKO LFC or z-score below a screen-specific cutoff), not just a negative interaction score. This matches Kebabci 2026
  and Chou 2025 (ZdLFC < −2 and ZLFC < −2).
- **not-SL:** tested, passes QC, interaction score clearly non-negative or non-significant.
- **untestable:** either single gene is essential in that line (pre-filtered by Thompson and Harle), guides failed, or no data.

Recommended per-screen rules:
- Dede: zdLFC < −3.
- Thompson / Harle / Chan: authors' filtered calls.
- Parrish: GI < −0.5 and FDR < 0.1.
- CHyMErA: authors' hit and LFC ≤ −0.9 (De Kegel).
- Ito / Klingbeil: GEMINI FDR < 0.05 **plus** the Kebabci dual threshold. Without it, Ito alone contributes about 1,700 "SL"
  pairs and dominates.
- in4mer: Chou ZdLFC rule.
- Flister: authors' calls, if the full table proves available.

Sensitivity analysis: re-derive calls with one uniform score (ZdLFC) where LFC tables allow it.

### 4.3 Pair-level aggregation
- n_lines_tested, n_lines_SL, n_studies_tested, n_studies_SL, penetrance = n_lines_SL / n_lines_tested.
- **Positive (primary):** SL in ≥2 cell lines, **and** either the lines come from ≥2 independent studies or penetrance ≥ 0.25.
  This targets recurrent SL, the only kind a context-free sequence model can predict.
- **Negative (primary):** tested in ≥2 cell lines (testable observations only) **and** never SL in any line **and** no line
  shows a borderline score. Borderline = within a margin of the threshold, e.g. zdLFC < −2 or GI < −0.3.
- **Ambiguous (excluded from primary AUPRC; kept as a separate stratum):**
  - SL in exactly one line;
  - discordant pairs (SL in some study but with low penetrance);
  - pairs tested in only one line;
  - pairs untestable in all lines.
- **Discordant calls between screens:** do not majority-vote. Use the rules above (recurrence + penetrance). Record every
  per-study call in the label file so stricter or looser variants can be rebuilt by script.
- **Secondary labels:**
  1. "Lenient": SL in ≥1 line vs never SL, matching Dennler & Ryan, for comparability with the literature.
  2. Continuous penetrance, for a ranking-correlation check.
  3. DepMap-inferred De Kegel set, as an independent auxiliary evaluation only.
- **Minimum evidence:** at least one study and two testable cell-line observations for any labelled pair. Positive controls
  (the 13-pair set) must come out positive; report it if any do not.

### 4.4 Expected size and positive rate (estimates; must be recomputed by script)

Union of pairs in screens I inspected:

| Screens | Unique pairs | Notes |
|---|---|---|
| Dede, in4mer, CHyMErA, Thompson-paralog | 5,193 | Recomputed from inspected files; overlaps are small (e.g. Dede∩in4mer 348, CHyMErA∩Thompson 171) |
| + Ito (4,858), Parrish (1,030), Klingbeil (~2,500, families partly overlap Ito), Harle/Chan (262) | ~8,000–10,000 | Estimate before ortholog/ID cleaning |
| + Flister genome-wide table, if obtainable | ~35,000 | Near all paralog pairs, but only 2 lines and very low hit rate |

Expected positive rates:
- **Lenient** (≥1 line): about 7–12% of screened pairs. Comparable sets: Dennler 403/5,316 = 7.6%; De Kegel per-screen hit
  rates 6–15%.
- **Primary recurrent rule:** roughly 300–600 positives against ~3,000–5,000 confident negatives, with ~3,000+ ambiguous
  pairs excluded. Positive rate about 5–10%; this is a rough estimate. For reference, De Kegel's ≥2-study consensus gave
  50/457 (10.9%) from 4 screens, and Chou found 80 pairs hit in ≥2 studies across 12 screens.
- **With Flister genome-wide pairs as negatives:** the positive rate falls below 1%. Keep those pairs as a separate
  "realistic-prevalence" evaluation set, not mixed into the main set.

## 5. Train / test split recommendation
- **Split by paralog family, not by pair.** Take connected components of the paralog graph built from our own sequence
  clustering, so that no gene appears in both train and test. Screens overlap heavily in pairs and families. Pair-level
  splits would leak, because genes like SAR1A/B, SEC23A/B and CNOT7/8 recur everywhere.
- **Final held-out test:**
  - About 20–25% of families, stratified by positive presence and by identity bin.
  - Selected once with a fixed seed and frozen before any modelling.
  - Evaluated only on primary-label pairs.
  - Report AUPRC with a bootstrap CI over families, alongside the sequence-identity baseline. The baseline is weak on screened
    pairs (AP ≈ 0.06 vs 0.04 chance in Kebabci), but likely stronger across all paralogs.
- **Secondary, study-held-out test:** leave out the newest independent platform. Harle + Chan (Sanger, 37 lines) is the
  natural choice; Flister is the alternative if its full table is usable. Evaluate on its pairs from test families only. This
  measures robustness to platform and selection bias.
- **Never use these for test labels:** SLKB / SynLethDB-derived calls, or the DepMap-inferred set. Report the DepMap-inferred
  set separately if at all.
- Keep yeast and human test families aligned if the yeast survey uses orthogroups. Otherwise human test genes can leak via
  yeast orthologs in training; the PI should check this.

## 6. License constraints for committing curated labels
- **Safe to commit derived labels** (pair IDs, per-study calls, counts, citation), with attribution, from CC BY sources:
  Dede, Thompson, in4mer, Chou, Harle, Chan, Ajmal, Kebabci, Dennler–Ryan (Zenodo CC BY 4.0), and Burgold. Code and
  repositories under MIT are also fine.
- **Copyrighted supplementary tables without an open license:** Ito (Nat Genet), CHyMErA (Nat Biotech).
  - Do not commit the raw tables.
  - Fetch them by script from the publisher URL and record checksums in DATA.md.
  - Committing derived binary calls (facts about gene pairs) is likely acceptable with citation. This is a judgement call, not
    verified legal advice; flag it in DATA.md.
- **CC BY-NC-ND articles:** Parrish, Klingbeil, Flister (journal).
  - The ND/NC terms restrict redistribution of adapted material; the project is non-commercial research.
  - Prefer deriving from the deposited data (GEO GSE178179 counts for Parrish; GEO GSE242517 for Klingbeil).
  - Do not commit copies of the supplementary tables.
  - Commit only derived calls, with attribution, and mark them in DATA.md.
- **SLKB data are GPL-3.0+:** do not commit anything derived from SLKB unless the label file is also released under GPL.
  Prefer the primary sources.
- **SynLethDB:** no stated database license; do not redistribute.
- **Raw downloads** stay in `data/raw/` (git-ignored), fetched by script with source, version, checksum and license in DATA.md.

## 7. Files inspected in this survey (`data/raw/_survey/`, all <5 MB; not for commit)
- `dede2020_tableS2_zdlfc.txt`: 403 pairs × 3 lines zdLFC (all pairs; 24 with zdLFC < −3).
- `thompson2021_MOESM4.xlsx` (SD1: all 1,191 pairs, FDRs and essentiality filters), `thompson2021_MOESM6.xlsx` (SD3:
  Ensembl IDs, pair origin), `thompson2021_SD5_SLpairs.xlsx` (hits only).
- `chymera2020_tableS8.xlsx`: 688 pairs; HAP1/RPE1 early/late effect sizes, FDRs, gi_type.
- `in4mer2024_SD1_goldstd.xlsx`: 388 / 26 / 13 hit lists.
- `chou2025_Supplementary_Table_1.xlsx` (dLFC/ZdLFC/RdLFC, 4,540 pairs × 4 lines), `chou2025_Supplementary_Table_2.xlsx`
  (371 hits × 12 screens; hits only), `chou2025_reference_paralog_13synthetic_lethals_221222.txt`.
- `flister2024_preprint_TableS3.xlsx` ("Lethal" calls only; no explicit negatives) and `flister2024_preprint_TableS4.xlsx`
  (guide library).
- `head_ens111_human_SL-lenient.csv`, `head_ens111_human_SL.csv`: first 30 kB of the Dennler–Ryan Zenodo files (columns only).
- Not downloaded (>5 MB): Ito supplementary tables (6.4 MB), Harle Additional file 1 (32 MB), Flister Tables S1/S2/S5, full
  Dennler files, SLKB dump.

**Open items for the data agent:**
- Verify Flister full per-pair tables (Table S1/S2) and the journal version's supplementary data.
- Get Klingbeil AACR supplementary tables.
- Get Parrish Table S4 (or recompute from GEO counts).
- Get Chan Supp. Tables 5–7 coverage.
- Get the Burgold paralog subset.

## 8. References
1. Dede M, McLaughlin M, Kim E, Hart T. Multiplex enCas12a screens detect functional buffering among paralogs otherwise masked in monogenic Cas9 knockout screens. Genome Biol 2020;21:262. doi:10.1186/s13059-020-02173-2
2. Thompson NA, Ranzani M, van der Weyden L, et al. Combinatorial CRISPR screen identifies fitness effects of gene paralogues. Nat Commun 2021;12:1302. doi:10.1038/s41467-021-21478-9
3. Parrish PCR, Thomas JD, Gabel AM, Kamlapurkar S, Bradley RK, Berger AH. Discovery of synthetic lethal and tumor suppressor paralog pairs in the human genome. Cell Rep 2021;36:109597. doi:10.1016/j.celrep.2021.109597
4. Ito T, Young MJ, Li R, et al. Paralog knockout profiling identifies DUSP4 and DUSP6 as a digenic dependence in MAPK pathway-driven cancers. Nat Genet 2021;53:1664–72. doi:10.1038/s41588-021-00967-z
5. Gonatopoulos-Pournatzis T, Aregger M, Brown KR, et al. Genetic interaction mapping and exon-resolution functional genomics with a hybrid Cas9–Cas12a platform. Nat Biotechnol 2020;38:638–48. doi:10.1038/s41587-020-0437-z
6. Köferle A, Schlattl A, Hörmann A, et al. Interrogation of cancer gene dependencies reveals paralog interactions of autosome and sex chromosome-encoded genes. Cell Rep 2022;39:110636. doi:10.1016/j.celrep.2022.110636
7. Esmaeili Anvar N, Lin C, Ma X, et al. Efficient gene knockout and genetic interaction screening using the in4mer CRISPR/Cas12a multiplex knockout platform. Nat Commun 2024;15:3577. doi:10.1038/s41467-024-47795-3
8. Chou J, Esmaeili Anvar N, Elghaish R, Chen J, Hart T. Z-scores outperform similar methods for analyzing CRISPR paralog synthetic lethality screens. Genome Biol 2025. doi:10.1186/s13059-025-03660-0
9. Klingbeil O, Skopelitis D, Tonelli C, et al. MARK2/MARK3 kinases are catalytic codependencies of YAP/TAZ in human cancer. Cancer Discov 2024;14:2471–88. doi:10.1158/2159-8290.CD-23-1529
10. Harle V, Offord V, Gökbağ B, et al. A compendium of synthetic lethal gene pairs defined by extensive combinatorial pan-cancer CRISPR screening. Genome Biol 2025;26:284. doi:10.1186/s13059-025-03737-w
11. Chan PY, Alexander D, Mehta I, et al. The synthetic lethal interaction between CDS1 and CDS2 is a vulnerability in uveal melanoma and across multiple tumor types. Nat Genet 2025;57:1672–83. doi:10.1038/s41588-025-02222-1
12. Flister MJ, Verduzco D, Petiwala S, et al. Genetic mapping and predictive modeling of paralog synthetic lethality. Cell Rep 2025;44:116512. doi:10.1016/j.celrep.2025.116512 (preprint doi:10.1101/2024.07.16.603642)
13. Burgold T, Karakoc E, Gonçalves E, et al. A next-generation dual guide CRISPR system for genetic interaction library screening. Nat Commun 2025. doi:10.1038/s41467-025-67256-9
14. De Kegel B, Quinn N, Thompson NA, Adams DJ, Ryan CJ. Comprehensive prediction of robust synthetic lethality between paralog pairs in cancer cell lines. Cell Syst 2021;12:1144–59. doi:10.1016/j.cels.2021.08.006
15. Dennler O, Ryan CJ. Evaluating sequence and structural similarity metrics for predicting shared paralog functions. NAR Genom Bioinform 2025;7:lqaf051. Data: doi:10.5281/zenodo.14973633
16. Ajmal H, Nandi S, Kebabci N, Ryan CJ. Benchmarking genetic interaction scoring methods for identifying synthetic lethality from combinatorial CRISPR screens. NAR Genom Bioinform 2025;7:lqaf129. doi:10.1093/nargab/lqaf129
17. Kebabci N, Ajmal H, Adams DJ, Ryan CJ. A predicted cancer dependency map for paralog pairs. bioRxiv 2026. doi:10.64898/2026.01.19.700065
18. Gökbağ B, Tang S, Fan K, et al. SLKB: synthetic lethality knowledge base. Nucleic Acids Res 2024;52:D1418+ (doi:10.1093/nar/gkad806). Data doi:10.6084/m9.figshare.22902839
19. Wang J, Wu M, Huang X, et al. SynLethDB 2.0: a web-based knowledge graph database on synthetic lethality for novel anticancer drug discovery. Database 2022:baac030. doi:10.1093/database/baac030
20. Rohde T, et al. BaCoN (Balanced Correlation Network) improves prediction of gene buffering. Mol Syst Biol 2025. doi:10.1038/s44320-025-00103-7 [P]
21. Fielden J, et al. Comprehensive interrogation of synthetic lethality in the DNA damage response. Nature 2025;640:1093–102 [P; not paralog-focused]
