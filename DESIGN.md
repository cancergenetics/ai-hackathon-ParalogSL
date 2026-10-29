# DESIGN — sequence-only prediction of synthetic-lethal paralog pairs

Version: **v3** (after review round 2) · 2026-10-28 · Author: PI agent.
Reviews: `reviews/round1/*.md`, `reviews/round2/*.md`; item-by-item PI responses in `reviews/round<N>/PI_response.md`.
**Normative implementation specs** (**precedence, D4:** a normative item wins over any sentence in this file unless this file explicitly names the item it overrides; `reviews/round2/statistics.md` and `external_pi.md` accepted items and data R2-13…R2-15 are also normative): label rules `reviews/round2/data.md` R2-1…R2-12; search mechanics, divergence calibration, duplication classification, masking and panel rule `reviews/round2/algorithms_comparative.md` §2.2–2.5 (changes 1–9 and should-items); gates `reviews/round2/agent_failure_modes.md` G1–G15 (subset adopted in §8).
Fixed by PROMPT.md (not open to revision): goal, input/output contract, sequence-only features, yeast + human labels, unseen-species generalisation, AUPRC, beating identity, hardware/disk limits, deliverables (installable CLI + tests + curated labels + human-proteome runtime/memory), reproducibility rules.

## 1. Problem statement and contract
- **Input:** one protein FASTA, one protein per gene, any eukaryote, possibly unannotated, possibly with fragments/allelic copies.
- **Output:** a TSV with every paralog pair the tool calls (`gene_a, gene_b, sl_score, sl_rank`, key features, QC flags such as `near_identical` for ≥99% identity / possible allelic or isoform copies, `fragment_suspect`), sorted by `sl_score`. Pair order inside a row is canonical (lexicographic); the score is symmetric by construction.
- **Allowed information:** the input sequences; a fixed, versioned reference panel of other species' proteomes (sequences only; headers stripped of names/annotations at setup); pretrained sequence/structure models (PLMs, ProstT5). **Forbidden as features or auxiliary targets (D1):** expression, PPI, essentiality, genetic interactions, GO/Pfam or any curated annotation, labels of orthologous pairs in other species, and any hand-entered knowledge of specific genes/pairs. Annotations (YGOB/OHNOLOGS WGD status, ribosomal-protein status) may be used **only for stratified evaluation**, never as inputs.
- **Supervision:** SL labels from S. cerevisiae and human only.

## 2. Data

### 2.1 Proteomes (gene is the unit)
- UniProt reference proteomes from one pinned release (recorded in DATA.md): human UP000005640, yeast UP000002311; unseen-species test proteomes Arabidopsis UP000006548, S. pombe UP000002485, C. elegans UP000001940 (D. melanogaster UP000000803 only if a test set materialises).
- **Per-gene FASTA:** UniProt collapses identical-sequence genes into one accession (verified: yeast HHF1/2, HHT1/2, TEF1/2, TIF1/2, EFT1/2; human histone genes). `scripts/data/build_gene_proteome.py` expands each accession into one record per gene (gene identifiers from the release's `.gene2acc` / `.idmapping`: Ensembl gene for human, ordered locus name for yeast, AGI for Arabidopsis, systematic ID for pombe, WormBase gene for worm). Identical copies are kept as separate genes. Mapping losses are logged per source.

### 2.2 Reference panel
- **Core:** QfO 2026_02 eukaryotes (51 species).
- **Gap filling by a written rule** (`scripts/data/build_panel.py`, `config/panel.yaml`), exactly as `reviews/round2/algorithms_comparative.md` change 8: UniProt reference proteomes per under-sampled clade, quality filters (BUSCO D ≤10%, S/F ladder, ≤50k proteins), one proteome per genus/order/family unit, explicit pre-/post-WGD Saccharomycetaceae genus lists, fixed quotas → +29 species, 80 species, ≈1.55 M proteins. Units already in the QfO core (incl. test-species congeners such as *Arabidopsis lyrata*/*Capsella*, *Caenorhabditis*) are not added, so the Arabidopsis and C. elegans tests resemble the "no close relative" use case. No non-*cerevisiae Saccharomyces* added (masks cover the close-relative regime). Metadata snapshots (JSON + `x-uniprot-release`) saved and checksummed; selection runs from the snapshot. Known gap: no gymnosperm/lycophyte passes quality filters. Choice uses only proteome metadata, never labels or results.
- Panel FASTA headers are reduced to `species|index` at setup; checksums are pinned. Core panel QfO 2026_02, gap-fill and training proteomes UniProt 2026_03 (recorded in DATA.md).
- **Self-exclusion:** panel species with median RBH identity ≥0.98 to the query are removed and logged — this removes self, conspecific strains **and sister species** (for human: chimpanzee 0.993, gorilla 0.990; macaque 0.973 kept). No sequence-only rule can separate a strain from a sister species; the behaviour is the same at training and inference.

### 2.3 Labels (pre-registered in `config/labels.yaml` before any label code is run; frozen before any model is fit)
Labels are built from primary tables by scripts (`scripts/labels/`), never from SLKB or SynLethDB. Unit: unordered pair of **genes** in one species. Each pair keeps per-source evidence (study, line/orientation, scores, call) so tiers can be regenerated. Thresholds are never tuned against the 13-pair human candidate list, against agreement with Dennler & Ryan, or against any model result.

**Label universe (R2-4; external-PI R2 M1):** labelled pairs = tested pairs ∩ paralog pairs in one pinned Ensembl Compara release (all within-species paralog homology types), for both human and yeast. Ensembl Compara is used **only** to define which tested pairs count as paralog pairs for labelling, never as a feature; caller recall is measured against this universe.

**One biological definition for both species:** *SL = a significant negative genetic interaction AND a strong absolute fitness defect of the double mutant, in a context where both single mutants are viable/testable.* "Synthetic sick or lethal" is the operational meaning in both species; this is stated in the report.

**Yeast (Costanzo 2016 raw SGA; null deletion alleles only in the primary tier; exact rules R2-1…R2-3, R2-9…R2-11):**
- Allele whitelist by strain-ID suffix; TS/DAmP alleles (essential genes) and `_S` suppressor-carrying strains excluded from the primary tier.
- **Testability floor:** a measurement is untestable if expected DMF = SMF_A × SMF_B ≤ 0.5 (a negative call would be near-impossible there).
- Measurement-level positive: ε < −0.08, p < 0.05 and DMF ≤ 0.5. Measurement-level negative: (ε ≥ −0.04, or ε ≥ −0.08 with p ≥ 0.05) and DMF > 0.5 (R2-1). Else borderline.
- Pair-level aggregation that does not reward measurement count, conflicts resolved exactly as R2-3 (flag column `orientation_conflict` kept; metrics also reported without flagged pairs).
- Dean 2008 / Kuzmin 2020 pairs labelled by their own rule (R2-9) from their fitness values; Dean kept independent for one QC run. Linked (adjacent) gene pairs handled per R2-11.
- Nested tiers exactly as R2-10 (strict ⊂ primary ⊂ lenient positives). BioGRID low-throughput SL only in a separate tier; AUPRC reported with/without BioGRID-only positives.
- Note: DMF ≤ 0.5 makes the yeast primary tier *stricter* than ε < −0.35 for pairs with healthy singles; the v2 claim that it raises positive counts was wrong (corrected; thresholds not re-tuned).

**Human (combinatorial screens; exact rules R2-5…R2-7, R2-15):**
- Per study × cell line call from **harmonised single- and double-knockout LFCs** (in4mer `Supp_table1.xlsx` for Dede, Thompson, Parrish, Ito, CHyMErA; Chou gene-level LFCs for in4mer; Harle/Chan/Klingbeil from their own tables): one uniform, screen-scaled rule — SL = the study's significance criterion **and** a DKO defect relative to the line's single-gene median (LFC − m) beyond a screen-scaled threshold; untestable by R2-5's screen-scaled single-gene rule. Gene IDs: Ensembl gene (override of R2-8's HGNC choice, D4). Ito's 90th-percentile cut dropped. Per-study exceptions (R2-6): Thompson authors' rule (reproduces 44/40/57); CHyMErA requires effect < 0, either-timepoint GI with late-sheet DKO, Torin1 excluded; Ito ERKi screen and in4mer selumetinib screen excluded; Ito pairs without GEMINI scores = untested. Expression-flag columns used for Harle/Chan testability never reach features (forbidden-column list).
- **Primary (robust) tier:** positive = SL in ≥2 distinct testable cell lines **and** (≥2 studies with penetrance ≥0.10 **or** penetrance ≥0.25 of testable lines); negative = testable in ≥2 cell lines, SL in none, with at most floor(0.1·n) borderline lines; same cell line across studies merged by Cellosaurus ID; else ambiguous.
- Sensitivity tiers: lenient (SL in ≥1 line vs never); ≥2-studies (De Kegel-style consensus).
- Sources: Dede 2020, Thompson 2021, Parrish 2021, Ito 2021, CHyMErA 2020, in4mer/Chou 2025, Klingbeil 2024, Harle 2025, Chan 2025. Flister 2025 is excluded from the primary tier (no verified non-hit table); De Kegel 2021 DepMap-inferred labels are not mixed in.
- Ascertainment controls: positive rate by source; a "source-library prediction" control (can the features predict which library a pair came from? if yes, report the confound).

**Unseen species (test only; never used for training, model selection or threshold tuning).** Endpoint semantics are stated per set:
1. Arabidopsis — Cusack 2021: strict positive = double mutant lethal with normal single mutants; negatives re-checked for lethal single mutants (excluded if so); broad tier as sensitivity.
2. S. pombe — E-MAP data with tested negatives (positive = strong negative score and severe double-mutant defect by the authors' thresholds) **only if** the full score matrix is obtained (currently behind a browser wall; manual download pending). Otherwise dropped and stated.
3. C. elegans — Tischler 2006 combinatorial RNAi (16/143).
4. Insects — no AUPRC test set exists; stated as a limitation.

**Label QC (scripted; `results/labels/`):** counts per tier/source/species; known-answer tests reproducing published counts (Dede 24 SL pairs; Harle 117/355; Costanzo threshold counts; Cusack 51/110); agreement with Dennler & Ryan's SL-Lenient set compared with *our lenient tier* (report only, never tuned); 13-pair candidate list reported (not tuned on); composition (identity, family size, closest-pair fraction) vs all de novo pairs.

**Commit policy:** derived labels (pair, tier, per-source evidence codes) committed with attribution. Rows whose only evidence comes from sources without an open license (Ito, CHyMErA) or CC BY-NC-ND (Parrish, Klingbeil) are committed as derived binary calls, flagged in DATA.md and in a column; raw tables never committed. Files containing forbidden columns (e.g. Dennler & Ryan feature columns) live in `data/raw/_quarantine/`; label code reads them only through an explicit column whitelist.

## 3. Paralog calling (part of the tool)
- MMseqs2 all-vs-all of the input proteome, `-s 7.5`, `--max-seqs 1000`, **`-a` (backtrace; required for exact identity)**; queries that hit the cap are re-searched with a larger cap. A candidate pair is defined by the bitscore-based E′ of algorithms_comparative change 3 (E′ ≤ 1e-3 with a fixed reference database size N_ref in `config/caller.yaml`, so calling does not depend on proteome size) in either direction (no identity or coverage floor; those are features). Scoring settings (gap, composition-bias correction) are pinned and identical in self- and panel searches, and included in the model's schema hash.
- Both orientations are realigned (MMseqs2 `align` on the pair list) so every pair has both directional statistics; self-scores are computed explicitly for every protein.
- **Label universe is independent of the caller** (tested pairs ∩ pinned Ensembl Compara paralogs, §2.3). Caller recall on it is reported per species/tier. Thresholds may be adjusted once on development labels; frozen in `config/`.
- End-to-end evaluation (§6.4) gives missed labelled pairs the minimum score.

## 4. Features (identical code for training and inference — training features are produced by running the shipped CLI's feature stage on the whole proteome)
All pair features are symmetric (order-invariant transforms: min, max, mean, |difference|). A **feature registry** (`src/slparalog/features/registry.py`) lists every allowed column; models accept only registered columns, and a test checks the registry against the forbidden-source list. No raw E-values (they scale with database size).

Tier A — self-search:
1. Similarity: identity (local; min/max over directions), Ensembl-style identity = identical positions / length of each protein (min, max), global (Needleman–Wunsch, parasail) identity, coverage of each protein, bitscore / self-bitscore (min, max), length ratio.
2. Family context: number of hits per gene at several bitscore-ratio levels (min/max over the pair), partner rank in each gene's hit list (min/max), reciprocal-best-hit flag, number of shared third paralogs, identity of each gene to its best *other* paralog and the gap to the partner, local family size defined by **mutual-coverage** links (both coverages ≥0.5) rather than single-linkage components. Proteome-percentile versions are computed over the full called-pair set (same at training and inference).
3. Per-gene sequence properties: length, compositional-bias/low-complexity fraction, `fragment_suspect` (length < 50% of the median homolog length).

Tier B — reference panel (MMseqs2 **one search per panel species**, `-a`; presence thresholds by bitscore, ≥50 bits, not E-value; ≈45–70 min for human, to be re-measured):
4. **Divergence calibration, no species tree needed:** pseudo-RBHs per species; core gene set G* (RBH in ≥80% of retained species); d_s = Kimura protein distance from the median identity over G*; d_s logged per run (algorithms_comparative change 5). Fixed bin edges B1–B6 in d (frozen before yeast d_s is computed). Validation without labels: Spearman of d_s with TimeTree divergence.
5. **Duplication classification per species** (change 6): pre / post / undetermined / absent from query-side bitscores; depth bounds L (max d over pre) and U (min d over post), depth, U − L, discordant fraction, censoring flags (`younger_than_nearest`, `older_than_panel`), relative rate d(A,B)/(2·depth). Conservation per gene/pair as fractions of retained species with a hit.
6. Divergence asymmetry: for species where A and B share an outgroup best hit O, |id(A,O) − id(B,O)| and its normalisation; identity of the more conserved copy to the outgroup.
7. "Paralogs rank" (Dennler): number of same-proteome proteins ranked above the partner (from Tier A), and "species rank": fraction of retained panel species with a protein scoring higher to A than B does.
- **No feature may be constant within a proteome** (registry test): per-bin features and bin-missingness patterns are excluded from the model (they would act as a species indicator with two training species); only pair-varying summaries are used. Species-identification control: a classifier predicting species from the model's features on the natural panel; its AUROC is reported.
- **Panel-density robustness (conditional adoption):** mask ladder by absolute ceilings c ∈ {none, 0.03, 0.12, 0.35, 0.70} (remove species with d_s < c; post-processing via a hidden CLI option `--panel-dmin`, everything downstream recomputed per mask; identical rows deduplicated with weight spread). First run the panel-thinning test without masking; **masking augmentation is adopted only if** thinning causes a development AUPRC drop that masking recovers. With masking: all copies of a pair share a fold, weights per pair, `min_samples_leaf` in pairs, no row-random early stopping, evaluation on one row per pair (natural panel).

Tier C — PLM (opt-in `--plm` tier unless the measured default stays ≤2 h; ablation decides; within-identity-strata evaluation mandatory; adoption also requires the human default run to stay ≤2 h):
8. Embedding similarity (cosine, Euclidean) of mean-pooled embeddings, starting with **one** model (ESM-2 150M); others (ESM C 300M, ProtT5-XL — Dennler found ProtT5 > ESM-2) only if it shows identity-conditional signal. Pooling variants tested only then.
(Per-gene embedding PCA components dropped — memorisation risk; listed under remaining ideas.)

Tier D — predicted structure (optional; only if Tier C signal suggests structure helps): ProstT5 encoder + CNN head → 3Di; Foldseek alignment of candidate pairs in 3Di+AA mode giving 3Di-alignment bitscore/identity (no TM-score/LDDT: no Cα coordinates exist). Evaluated first on the labelled-gene subset.

**WGD-proxy validation (prior-work M3b, comparative M6):** at M2, the duplication-depth and co-orthology features are compared against YGOB (yeast) and OHNOLOGS v2 (human) WGD status — evaluation only, no labels involved — and reported.

Not planned: rate correlation from MSAs; predicted heteromerisation (opposite-sign evidence in yeast vs human; PLM-interact trained on measured PPI → excluded).

## 5. Models
- Baselines: (B0) random; (B1) **identity alone** — the best of the pre-registered identity variants {local min, local max, Ensembl-style min/max/mean, global, identity×coverage}, chosen on development folds, scored on the same pairs as the model; (B2) identity + family size + closest-pair flag (logistic regression) — pre-registered secondary comparator; (B3) gene-prior leakage baseline.
- Main model: gradient-boosted trees (sklearn HistGradientBoosting, or LightGBM in a separate process from torch — libomp conflict), species-agnostic features only, sample weights giving each training species equal total weight **and equal weighted prevalence**; no row-random early stopping (iterations chosen by grouped CV). Logistic regression on the same features as comparator. Hyper-parameters from a small pre-declared grid; **selection metric = development species-mean ΔAP vs B1, one-standard-error rule**.
- Auto-research (all feature/model exploration) happens on development data **before** the single test access, driven by a backlog in RESEARCH_TREE.md with a stop rule: stop when two consecutive ideas fail to improve dev species-mean AP beyond one SE, or the backlog is empty.
- Every development configuration is logged in `results/experiments.tsv` (id, git hash, features, params, dev metrics per species) with its adoption rule declared before the run.
- The released model is refit on all labelled human + yeast data **after** the final test evaluation is recorded; the refit model is never scored on the test set.

## 6. Evaluation protocol

### 6.1 Splits (built once by script, frozen, checksummed, committed)
- Graph over all labelled genes of both species. **Each labelled pair is first contracted into one supernode** (a labelled pair can never span the split). Edges = MMseqs2 hits within and across species with a stricter rule (bitscore ratio and mutual coverage thresholds pre-registered in `config/splits.yaml`); default = connected components on these edges. Components above a pre-registered size cap are partitioned by Leiden on bitscore-ratio-weighted edges (fixed seed and resolution rule in config) — fallback only.
- **Test:** fraction of groups set at M1 by the power simulation (below), between 30% and 40%, stratified so each species keeps that fraction of its positives and yeast has ≥30 test positives. The rule moving very large groups into development is reported with what it removes.
- **Leakage measurement:** for every test pair, the maximum bitscore ratio of its genes to any training gene; ΔAP reported as a dose–response by leakage tertile with a trend CI (the no-hit subset is a different population, reported descriptively only).
- **Development:** 5-fold grouped CV over development groups, repeated over the 5 fixed seeds.
- **Leave-one-species-out (LOSO)** within development: yeast → human, human → yeast.
- **Power:** at M1, before any feature exists and without computing any score on test pairs (identity included), a simulation (`scripts/eval/power_sim.py`; inputs pre-registered: score correlation ρ = 0.8, identity-baseline AUROC from the literature, real group sizes and positive counts) computes power of the full decision rule at T = 0.04 AP. The test fraction is the smallest in [0.30, 0.40] giving power ≥0.8; if none does, 0.40 is used and the primary comparison is pre-declared underpowered. A non-significant result is reported as "not demonstrated", not as "no effect". **No nested-CV fallback** (it would be biased by development-time choices).

### 6.2 Metrics
- AUPRC = `sklearn.metrics.average_precision_score` (pinned version). Prevalence and lift always shown. Secondary: AUROC, precision at top-k, lenient/strict/≥2-studies tier AUPRC, WGD vs SSD and with/without ribosomal pairs (YGOB/OHNOLOGS, UniProt keywords — evaluation only).
- **Primary pair population:** labelled pairs (primary tier) that the caller calls; B1 and the model are scored on exactly the same pairs (enforced in code).
- Uncertainty: group bootstrap, 10,000 resamples, seed from `config/success.yaml`; groups resampled with replacement within strata defined by species composition (human-only, yeast-only, mixed) × (group contains ≥1 positive or not); percentile 95% CI; paired for Δ; leave-one-group-out influence check on the largest groups.
- **Reporting hierarchy:** (1) primary: robust tier, test, human + yeast, T vs B1; (2) secondary: Δ vs B2, LOSO, unseen species, lenient/strict tiers, end-to-end; (3) exploratory: WGD/SSD, ribosomal, leakage tertiles, mask ladder, robustness.
- **Primary success criterion (pre-registered):** T = mean over {human, yeast} of ΔAP_s = AP(model) − AP(B1) on the frozen test set; pass if the 95% CI of T excludes 0 **and** ΔAP_s > 0 in each species. "Clearly beats identity" additionally requires mean relative gain (AP_model/AP_B1 − 1) ≥ 0.20 over the two species, judged on the point estimate. Secondary: Δ vs B2; LOSO ΔAP; unseen-species ΔAP with CI (pre-registered as underpowered; no pass/fail).
- Expected precision at genome-wide prevalence (reweighted negatives; human ≈0.5% per Flister 2025).

### 6.3 Test-set discipline
- **Split mechanics (D4):** `config/splits.yaml` (protected) fixes edge rule, component size cap (components holding >10% of a species' labelled supernodes go to the Leiden fallback), Leiden seed/resolution, and **one draw per candidate test fraction (no redraws)**; `config/mde.yaml` (protected) fixes the power-simulation inputs. If yeast has <75 primary positives, the yeast test target becomes ⌈0.40·n_pos⌉ and the yeast comparison is pre-declared underpowered. `scripts/labels/make_splits.py`, `scripts/eval/power_sim.py` and `tests/test_splits.py` are allowlisted to read sealed counts.
- Test labels, the split-assignment file and unseen-species labels live in `data/labels/sealed/` (with `MANIFEST.sha256`); development code only sees `data/labels/dev/`. Raw unseen-species sources are read only by their label-builder scripts (allowlisted), which write only to `sealed/` and print counts. `load_labels()` raises on sealed paths (G1).
- Before the test is opened, the final model and config are committed and tagged `final-model-v1`; `final_test.py` checks its own preconditions (G15: tag, clean tree, model SHA, sealed manifest, prior accesses committed, exactly one primary access) and appends every invocation (git hash, config hash, timestamp, outcome) to `results/test_access_log.tsv`.
- One planned access. A re-run requires a reviewed bug, a development-only test that exposes it, and reporting both numbers.

### 6.4 End-to-end evaluation
- AUPRC over all labelled pairs with caller-missed pairs at the minimum score, alongside AUPRC over called pairs.

### 6.5 Controls and verification (scripts + pytest)
- Toy proteome with planted duplicates and a toy panel → known expected values for caller and every feature.
- Shuffled-label negative control (AUPRC ≈ prevalence). Gene-prior baseline B3. Source-library prediction control.
- Leakage audit script (gene overlap must be 0; homology leakage quantified as in §6.1).
- Identity-stratified analysis for Tier C/D features.
- Invariance tests: input order, ID renaming, `score(A,B)=score(B,A)`, determinism across thread counts; golden-file regression on a small proteome.
- **Robustness test (pre-registered):** human and yeast development proteomes degraded by simulated truncation (fragments), split and fused genes, 10–20% missing genes and haplotig-like near-identical copies; report AUPRC drop.
- Panel-thinning test (§4).
- Independent review: every key result is re-checked by an agent that did not produce it (including a re-run from a clean checkout) and the review is recorded in `reviews/results/` before the tree node is marked [x] or [!].

## 7. Compute plan and budget
- Per-proteome run, human (measured in Low Power Mode, to be re-measured): MMseqs2 self-search ≈10 min; per-species panel search ≈45–70 min (80 species); ESM-2 150M ≈1.7 h if Tier C is adopted (candidate-pair genes hold 88% of residues, so restriction saves little); 3Di optional (≈11 h). Target: default tier ≤2 h and ≤8 GB peak for human; if Tier C would push the default beyond 2 h it ships as an opt-in `--plm` tier. Runtime per 10k proteins documented; no silent subsampling.
- Disk: raw data a few GB, panel ≤1.5 GB, model weights ≤25 GB (unused weights pruned), interim ≤20 GB. `du -sh` at each milestone; total ≪150 GB.
- Final runtime/memory of the CLI on the full human proteome with `/usr/bin/time -l`, recording power state (`pmset -g`), thread count, and versions.

## 8. Software deliverable and reproducibility
- Package `slparalog` (src layout, `pyproject.toml`, console script `slparalog predict proteome.fasta -o pairs.tsv`; `slparalog setup` fetches the panel and model weights with pinned checksums). Model bundle stores feature-schema hash and search parameters, checked at run time.
- Pinned environment: `requirements-lock.txt` (pip freeze), MMseqs2/Foldseek versions and checksums in `config/tools.yaml`, Hugging Face model revisions pinned.
- Tests: unit tests on toy data, label-parser fixtures with known answers, CLI smoke test, invariance tests, sealed-path guard, gene-name scanner (fails on hard-coded gene symbols/accessions outside fixtures and config label-source files).
- Numbers in REPORT.md and the tree are generated from `results/*.json` and checked by `scripts/check_report_numbers.py`.
- Protected files (label/split/seed/success/caller configs, metric code, final_test, leakage audit, controls, split files, sealed manifest, guard tests, `tests/conftest.py`, `pyproject.toml`) are listed with SHA-256 in `config/protected.sha256`; a test fails if they change without a matching decision note (hash prefix) in RESEARCH_TREE.md (G3).
- **Gates adopted** from `reviews/round2/agent_failure_modes.md`: M1 — G1 data access guard, G2 data registry, G3 protected hashes, G4 label known answers (with tolerances; unreproducible counts as strict xfail with note, never deleted), G5 metric toy tests and same-pair enforcement, G6 split reproducibility and invariants, G7 `scripts/milestone_check.py` (tests pass, no unlisted skip/xfail, clean tree, disk < 150 GB, lock files present; the git-ancestry ordering checks are applied as a documented manual check, not code). M2 — G8 feature registry, G9 narrow ID scanner, G10 train/inference parity at mask none, G11–G12 controls, G14 review-record format. M4 — G15.
- Every delegation to a sub-agent includes `docs/AGENT_RULES.md` (extended in round 2: no editing failing tests/tolerances/pytest config/skip markers, no git-history rewriting, no editing RESEARCH_TREE/DESIGN, no reading raw unseen-species files, every reported number cites its result file and commit).

## 9. AI-agent failure modes and mitigations (summary; see review `reviews/round1/agent_failure_modes.md`)
Test peeking → sealed directory + tag + access log. Silent changes → protected hashes + decision notes. Fabricated numbers → results JSON + report number check + independent reviews. Hard-coded knowledge → gene-name scanner, registry of feature columns, feature cards written before results are seen. Train/inference skew → training features produced by the CLI. Iterative dev overfitting → experiment registry + one-SE rule + single final model. Context loss → research tree.

## 10. Milestones
M1 data + per-gene proteomes + panel + labels + splits + MDE (frozen) → M2 tool skeleton + Tier A/B features + baselines + controls → M3 model development (Tier C/D ablations, robustness, panel thinning) → M4 freeze, tag, final test, CLI packaging, human benchmark → M5 REPORT.md.
