# Round 1 review: AI-agent failure modes

Reviewer role: auditor of the ways an autonomous PI agent and its sub-agents could fail while running this design.
Inputs read: PROMPT.md, DESIGN.md v1, LITERATURE.md, docs/lit/human_labels.md (parts), DATA.md, LESSONS.md,
RESEARCH_TREE.md, .gitignore, `git ls-files`, config/seeds.yaml, scripts/bench/*, data/raw/{proteomes,qfo}/README.md,
data/raw/qfo/species_list.tsv, a listing of data/raw/_survey/, and the header of data/raw/_survey/head_ens111_human_SL-lenient.csv.
Nothing was run beyond `ls`, `head`, version flags and `pip freeze`.

## 1. Overall verdict

The design names the right failure modes (§9, §6.3, §6.5), but almost every safeguard is a promise the same agents are
trusted to keep: the test-access log only works if the agents choose to go through `final_test.py`, "no gene names in
code" has no check, and "an independent agent re-checks" has no record format. Some rules are already broken before M1:
DATA.md is empty although about 40 raw files are on disk, benchmark logs sit in the gitignored `logs/` directory so the
compute numbers in LITERATURE.md are not committed, and a file holding forbidden data-derived features plus SL labels
(the Dennler head CSV) lies unquarantined in `data/raw/_survey/`. The design is sound once the promises below are turned
into scripts and gates that fail loudly. Most of them cost hours, not days.

## 2. Concrete changes (priority · change · justification)

### A. Test-set discipline

1. **MUST: physically separate the test labels, and make tests enforce the separation.** `scripts/labels/make_splits.py` writes
   `data/labels/dev.tsv` and `data/labels/sealed/test.tsv` (plus `data/labels/sealed/unseen_*.tsv`). A pytest
   (`tests/test_sealed_access.py`) greps `src/`, `scripts/` and notebooks for the string `sealed/` and for the sealed file
   names, and fails on any hit outside `scripts/eval/final_test.py` and the split builder. The dev loaders
   (`slparalog.data.load_labels`) raise if they are given a sealed path.
   *Why:* §6.3's access log is opt-in. Any agent can run `pd.read_csv("data/labels/splits.tsv")` and filter to the test rows
   without leaving a trace.

2. **MUST: freeze the model before the test is opened, and prove it.** Before the first test access, commit
   `results/frozen/model.txt`, `feature_list.json` and `config_hash`, and tag the commit (e.g. `git tag test-freeze-v1`).
   `final_test.py` refuses to run unless HEAD is that tag, or a descendant whose diff touches only `results/`. It logs the tag
   SHA, the model SHA-256 and the sealed-file SHA-256 to `results/test_access_log.tsv`. A reviewer agent checks
   `git log --follow` of the log and of `sealed/` at M4.
   *Why:* this makes "no model or feature selection after first access" auditable from git, not from memory.

3. **MUST: tighten the "re-run if a bug is found" loophole.** A re-run is allowed only when (a) an agent that did not write
   the code reviews the bug and writes it up in `results/reviews/`, (b) a dev-only test exposes the bug and fails before
   the fix, and (c) REPORT.md shows both the original and the re-run numbers.
   *Why:* after seeing a disappointing test number, an agent can relabel almost any change as a "bug fix".

4. **MUST: the released model is never scored on the test set.** Every metric JSON records the model SHA-256 and the
   list of training-pair hashes. `scripts/report/make_tables.py` refuses to put a test metric into REPORT.md when that
   model's training set overlaps the test pairs. The refit model (§5) gets only a sanity check that its scores correlate
   with the frozen model's on dev pairs, and that number is labelled as such.
   *Why:* this prevents a contaminated number from being reported as a test result, whether by accident or on purpose.

5. **SHOULD: build the split once and pin it.** `make_splits.py` uses seed 11 only and is deterministic. It writes
   `splits.sha256`, and the split files get exactly one commit. Any later change needs a decision note D-n and a
   `splits_v2` file, with v1 results kept.
   *Why:* this stops "split shopping", i.e. regenerating the split until the test fold looks easy, including through the
   >5%-cluster rule and the stratification retry loop.

### B. Labels, thresholds, metrics: no silent changes

6. **MUST: put the protected files under hash control.** `config/protected.sha256` lists `config/labels.yaml`,
   `config/splits.yaml`, `src/slparalog/eval/metrics.py`, `scripts/eval/final_test.py`, `scripts/eval/leakage_audit.py`,
   the shuffled-label control, and the success-criterion constants. `scripts/check_protected.py` runs in pytest and at
   every milestone. A mismatch fails unless a RESEARCH_TREE decision note cites the new hash.
   *Why:* the commonest reward hack is editing the check (the audit, the tolerance, the metric) rather than the code it
   checks. Making the rule mechanical catches it.

7. **MUST: freeze all thresholds before the label code first runs, and never tune them on QC outputs.** Commit
   `config/labels.yaml`, with every per-study threshold and its citation (page/table), *before* any label is built.
   Agreement with Dennler & Ryan and the 13-pair reference list are reported as they come out. If they look wrong, the
   action is a bug hunt that is reviewed, documented and run with the threshold values untouched. Any threshold change
   becomes a new tier (`labels_v2`), with both tiers reported.
   *Why:* the design's QC requirement that the 13 pairs "must come out positive or ambiguous" invites an agent to loosen
   thresholds until the well-known pairs pass. That fits labels to gene-level knowledge, which PROMPT.md forbids.

8. **MUST: test the label code against known answers from the papers.** Add pytest checks that reproduce published counts
   from the primary tables before our harmonised rules are applied. Examples: Harle 2025 gives 117 pairs SL in ≥1 line and
   355 never SL; Kebabci's Ito recipe gives 958 SL of about 41k pair×line points; Costanzo's stringent ε<−0.12 count; Cusack
   strict gives 51/110. Tolerance must be stated, e.g. ±2% or an exact match. Where a count cannot be reproduced, the
   discrepancy is logged and reviewed, not hidden.
   *Why:* parser and orientation bugs (query/array swaps, duplicated domain-tag rows in Klingbeil, symbol remapping) are
   silent. A published count is the cheapest known-answer control there is.

9. **MUST: pin the metric implementation.** AUPRC = `sklearn.metrics.average_precision_score`, as stated in
   `metrics.py`, with no trapezoidal alternative. A toy test with a hand-computed AP checks it, and so does a test that
   identical scores (ties) are handled identically for the model and for B1. Prevalence is computed from the same
   pair set the AP uses.
   *Why:* switching AP and trapezoid AUPRC, or tie handling, changes Δ AUPRC by amounts comparable to the expected effect.

10. **SHOULD: freeze the identity-baseline menu now.** List the B1 variants in DESIGN §5: min/max/mean directional
    identity, local vs alignment-length-normalised, bitscore/min-self-score. B1 picks among these on dev folds only, and
    the baseline is scored on *exactly* the pair set the model is scored on, in both the candidate-only and the
    end-to-end views (missed pairs at the minimum for both).
    *Why:* the cheapest way to "beat identity" is to give it a weak variant or a different pair set. The prompt's headline
    criterion depends on B1 being fair.

11. **SHOULD: pre-commit to reporting both evaluation views and both label tiers.** REPORT.md tables are generated for
    candidate-only and end-to-end, primary and lenient, every species, whichever looks better.
    *Why:* this blocks selective reporting of the favourable view.

### C. Data and knowledge leakage, including LLM priors

12. **MUST: quarantine the Dennler/Ryan files and check feature-matrix schemas.** Move `head_ens111_human_SL-lenient.csv`
    and any future Dennler Zenodo download to `data/raw/_quarantine/` with a README. That file contains GTEx correlations,
    PPI sets, complex essentiality, `has_essential_cerevisiae_ortholog`, Dennler's own rank features and the `SL` label.
    The training-matrix builder accepts only columns registered in `src/slparalog/features/registry.py`. A pytest asserts
    that `train_matrix.columns ⊆ registry` and that no column name matches `gtex|ppi|essential|complex|ceres|go_|pfam|SL$`.
    *Why:* an agent "speeding up" feature work by merging columns from a ready-made table is the most likely way forbidden
    data enters. It is a one-line `merge`.

13. **MUST: add a scanner for hard-coded gene knowledge.** `tests/test_no_hardcoded_genes.py` scans `src/`, `scripts/` and
    `config/` for HGNC-like symbols (via a regex for gene-symbol tokens, checked against the HGNC list loaded at test
    time), yeast systematic ORFs (`Y[A-P][LR]\d{3}[WC]`), UniProt accessions, AGI codes and WormBase IDs. Anything outside
    `tests/fixtures/` and an explicit allowlist fails the test. The allowlist is reviewed at each milestone.
    *Why:* §9 says "no gene names in code" but has no check. Agents routinely add "sanity" special cases.

14. **MUST: mask gene identity for feature-development agents.** Agents who design or ablate features get dev pairs
    with opaque IDs (hash of the accession pair) and no symbols. A *feature card* (`docs/features/<name>.md`) states the
    rationale, the literature source and the decision rule *before* the feature's dev result is computed.
    *Why:* an LLM looking at "HDAC1–HDAC2, MAPK1–MAPK3 are false negatives" will, from its priors, engineer features that
    fit known families (length windows, composition, kinase-like signals). Once pairs are named, LLM gene knowledge
    becomes a feature that passes no code check.

15. **SHOULD: derive the species tree and phylostrata by script, not from model memory.** Feature 4's "fixed species
    tree for the panel, distances stored in `config/`" must be generated from NCBI Taxonomy (or TimeTree) files listed in
    DATA.md, by `scripts/panel/make_tree.py`.
    *Why:* an agent typing divergence times from memory is a silent, unverifiable assumption, and it will be partly wrong.

16. **SHOULD: verify the self-exclusion rule empirically, and check its train/inference asymmetry.** The panel contains
    chimpanzee, gorilla and macaque. Under the "≥50% of proteins at ≥98% identity" rule, the human run will probably also
    drop chimp and gorilla, possibly macaque, while yeast's nearest panel relative is Candida albicans. Add
    `scripts/panel/exclusion_report.py`, which logs the excluded species per run and prints the per-species fraction of
    near-identical hits, and pytest controls. *Positive control:* the human, yeast, Arabidopsis, S. pombe and C. elegans
    runs each drop themselves. *Negative control:* a panel species is not dropped when the query is a distant species.
    Also add one "novel-genome simulation" evaluation that masks the query's whole clade (for example all mammals for
    human). The comparative-genomics reviewer should judge the biology.
    *Why:* the rule is a silent assumption with a numeric threshold. Panel-composition differences between the training
    species and a novel plant or insect are a direct train/inference skew in every Tier B feature.

### D. Train/inference skew

17. **MUST: generate training features with the shipped CLI, on whole proteomes.** Training matrices are produced by
    `slparalog features <full proteome FASTA>` (the same entry point as `predict`), and only then joined to labels.
    Features are never computed on a labelled subset or for pairs the caller did not emit. A pytest runs `predict` on a
    500-protein subsample and asserts its features are bit-identical to the training-path features for the same pairs.
    *Why:* family-context features (percentiles, component size, ranks) depend on the whole proteome. Subset computation
    or "forced" alignment of labelled-but-uncalled pairs creates skew that cross-validation cannot see.

18. **SHOULD: decide now whether labelled pairs the caller misses can be training examples** (recommended: no, because
    they cannot occur at inference) and write it in DESIGN §3.
    *Why:* if this is left open, an agent will decide it silently.

19. **SHOULD: add a robustness check for proteome quality.** Re-run dev evaluation on a perturbed proteome: 10%
    fragmented proteins, 5% duplicated alleles at 99% identity, a few hundred extra isoforms. The perturbation is scripted
    and seeded. Report the AUPRC change.
    *Why:* the target use case (newly sequenced genomes) violates the clean "one canonical protein per gene" assumption
    of the UniProt training proteomes. Near-identical allele pairs will sit at the top of an identity-driven ranking.

### E. Overfitting the development set and gaming the comparison with identity

20. **MUST: keep an experiment registry and pre-specify decision rules.** Every dev evaluation appends to
    `results/experiments.jsonl` (git SHA, feature set, hyper-parameters, seeds, per-species dev AUPRC with CI). A feature
    tier is adopted only if a rule written *before* the ablation passes: for example, mean Δ dev AUPRC > 0 in both species,
    the paired 95% CI over clusters excludes 0 in at least one species, and LOSO does not worsen. REPORT.md states the
    total number of configurations tried.
    *Why:* with roughly 300–600 human and roughly 1k yeast positives, dozens of ablations over seven feature groups and a
    hyper-parameter grid will "find" dev gains that are selection noise. The registry makes the garden of forking paths
    visible.

21. **SHOULD: add two pipeline controls beside the shuffled-label one.** (a) *Identity-only model through the full
    LightGBM pipeline* must reproduce B1's AUPRC within CI; a large excess means leakage or a bug. (b) *Planted-signal
    control:* synthetic labels drawn from a known function of two features at the real prevalence must be recovered, with
    AP near the oracle. Give the shuffled-label control an explicit tolerance: AP within the bootstrap 95% interval of
    prevalence.
    *Why:* these turn "the model beats identity" from a single number into a comparison with known expected answers, and
    catch target leakage through features such as "number of labelled partners".

22. **SHOULD: write the stopping rule and the negative-result commitment into DESIGN §10.** Ideas are queued in
    RESEARCH_TREE with an expected-value note. Work stops after K (e.g. 3) consecutive pre-registered ideas fail the
    adoption rule, or at the time budget. If the success criterion fails, REPORT.md says so as the headline.
    *Why:* "keep going until you run out of ideas" plus a headline goal of "beat identity" is the textbook setup for
    endless dev-set mining or for quietly redefining success.

### F. Fabricated or unverified numbers, unreviewed results, context loss

23. **MUST: generate every number in REPORT.md and RESEARCH_TREE.md from results.** Numbers in REPORT.md and in RESEARCH_TREE
    status summaries come from `results/**/*.json` through `scripts/report/render.py` (template keys such as
    `{{human.test.ap}}`). `scripts/check_numbers.py` flags any decimal number in REPORT.md that does not trace to a JSON
    key or to a cited literature value. Sub-agent summaries must cite the JSON path and commit SHA for every number they
    report. The PI spot-checks with `jq`.
    *Why:* PI context is kept lean by reading summaries. A sub-agent's mistyped or invented number would then propagate
    unchecked.

24. **MUST: give independent review a record format.** Each [x]/[!] node in RESEARCH_TREE links to
    `results/reviews/<id>.md`, which states the reviewer agent, the commit reviewed, what was *re-run from a clean
    checkout*, and whether the numbers matched within a stated tolerance. `scripts/check_tree.py` fails if a [x]/[!] node
    lacks the link. The reviewer gets the code and the results, not the producer's conclusion.
    *Why:* "re-checked by a second agent" without a re-run usually degrades into a second agent agreeing with a summary.

25. **MUST: write a standard sub-agent preamble.** `docs/AGENT_RULES.md` (about one page) is pasted verbatim into every
    delegation. It covers forbidden features and sources (D1), the sealed-file rule, no gene names or knowledge, the
    protected files, the required report format (paths, SHAs, JSON keys), the disk/CPU lock, and "do not edit RESEARCH_TREE
    or DESIGN". Only the PI edits RESEARCH_TREE.md.
    *Why:* sub-agents start with no context. Rules that live only in DESIGN.md and the PI's head get violated by the third
    delegation, and after a context reset the PI needs the same document.

26. **SHOULD: verify the [U] literature facts that drive design numbers before they are used.** There are 20 [U] marks in
    docs/lit. Items that set thresholds or reweighting are the priority: Flister's genome-wide prevalence (<0.5%), which is
    used to reweight to genome-wide precision; the Flister non-hit table; the S. pombe E-MAP negatives. Also fix a provenance
    inconsistency: DESIGN §2.3 cites the 13-pair list as "Esmaeili Anvar 2024 SD1", but the file on disk is
    `chou2025_reference_paralog_13synthetic_lethals_221222.txt`; and the Flister files on disk are 2024 *preprint* tables
    while the design cites Cell Rep 2025.
    *Why:* unverified secondary numbers, once in DESIGN.md, become "facts" for later agents.

### G. Reproducibility

27. **MUST: backfill DATA.md now, by script, and gate M1 on it.** `scripts/data/register.py` hashes every file under
    `data/raw/` and `models_cache/` and writes or updates DATA.md rows (URL, version, date, SHA-256, license, committed?).
    `scripts/data/check_registry.py` fails if any file on disk lacks a row. Today DATA.md is empty while `data/raw/_survey/`
    holds about 25 files, several of them unlicensed or restricted (CHyMErA, Flister, FlyBase). The QfO and proteome
    READMEs hold checksums but are not DATA.md.
    *Why:* PROMPT.md requires it, and the survey phase already broke the rule. Survey downloads made ad hoc should either be
    re-fetched by `scripts/data/download_*.py` or registered with their URLs.

28. **MUST: pin the environment.** Commit `requirements.lock` (from `pip freeze`; currently torch 2.14.1, lightgbm 4.7.0,
    sklearn 1.9.1, …), `python --version`, and `config/binaries.yaml` with the MMseqs2 commit (8cc5ce3…), the Foldseek
    commit (941cd33…), DIAMOND 2.2.8, binary SHA-256s and download URLs. `bin/` is gitignored, so these are recorded nowhere
    today. Record the HuggingFace model revisions (commit hashes) for each PLM used, and load them with `revision=`.
    *Why:* "end-to-end reproducible" cannot be checked without them. PLM checkpoints on the Hub can change under a fixed name.

29. **MUST: commit benchmark and result outputs that the reports cite.** `logs/` is gitignored, yet LITERATURE.md §4
    quotes runtimes from `logs/bench/*.jsonl`. Copy small summaries (jsonl and time files, not the .npy embeddings) to
    `results/bench/` and commit them. Rule: any number cited in a document lives in a committed file under `results/`.
    *Why:* PROMPT.md says "commit results to git", and right now the cited numbers cannot be traced.

30. **SHOULD: make the seed policy concrete and test determinism.** Map each stochastic step to its seeds in
    `config/seeds.yaml` (CV folds: all 5; bootstrap: seed 11 only; split: 11). Set LightGBM `deterministic=True`, a fixed
    `num_threads` and `seed`. Run the toy and 500-protein pipeline twice in pytest and diff the features. MPS fp16 embeddings
    are not bitwise deterministic, so set a stated tolerance (cos ≥ 0.9999) and report whether the final model's AP moves.
    MMseqs2 `--max-seqs` truncation under ties should be checked the same way. Remove hard-coded
    `ROOT = Path(".")` from scripts in favour of a repo-relative root.
    *Why:* determinism claims fail silently. Absolute paths break the "installable tool" and any clean-checkout re-run by
    a reviewer.

31. **SHOULD: commit at every result.** RESEARCH_TREE.md is currently modified but uncommitted. Use a milestone checklist
    script (`scripts/milestone_check.sh`) that runs check_protected, check_registry, check_tree, pytest, `du -sh` and
    `git status --porcelain` (must be empty).
    *Why:* uncommitted state is lost on context reset and cannot be audited.

### H. Wasted compute and resource contention

32. **SHOULD: serialise heavy jobs.** Use a single job runner (`scripts/run_heavy.sh`, wrapping `flock` on
    `logs/heavy.lock`) for MMseqs2 panel searches, PLM embedding and ProstT5. Agents running in parallel never run two at
    once. Record `pmset -g` and Low Power Mode state with every timing.
    *Why:* with 10 cores and 24 GB, two agents embedding at the same time causes swapping, and the slow timings then get
    reported as the tool's runtime. The bench phase already ran entirely in Low Power Mode.

33. **SHOULD: cache embeddings and searches by content hash.** Key embeddings by (SHA-256 of sequence, model id, HF
    revision, pooling, chunking) in `data/interim/cache/`. Gate the expensive tiers (C beyond 35M, D) on a pre-registered
    subset experiment with a time budget written in DESIGN §7, e.g. ≤ 6 h of compute per tier before a go/no-go decision.
    Prune `models_cache` (currently 24 GB) to the models actually used.
    *Why:* re-embedding a proteome after a pooling change costs 1–3 h each time. Agents re-run more than they think.

## 3. Risks the design does not address

- **Evaluator and producer are the same agent.** Nothing says the agent that builds `final_test.py`, the metric and the
  leakage audit differs from the agent that develops the model. Recommend that a separate "evaluation owner" sub-agent
  writes and reviews these files and that modelling agents cannot edit them (covered mechanically by item 6).
- **Survey agents have already seen the unseen-species test labels.** Cusack, FlyBase and Tischler-related counts appear in
  LITERATURE.md, and `cusack2021_*` files are on disk. Seeing counts is harmless, but feature-development agents must not
  open these files (item 1 should cover `data/raw/_survey/cusack*` and `flybase*` as sealed sources too).
- **Label-source selection after seeing agreement.** Sources are listed, but nothing stops an agent from dropping a
  "noisy" screen after seeing that it hurts AUPRC. The source list should be frozen in labels.yaml (item 7), and dropping a
  source should count as a label change.
- **The Dennler cross-check becomes a target.** If agreement is low, an agent may tune harmonisation until it matches
  Dennler, which turns their labels into our thresholds by the back door. Report agreement once; do not iterate on it
  (item 7).
- **Ascertainment shortcut.** The model can learn "looks like a screened pair" (high identity, WGD-like, closest pair). On
  the screened candidate set this is indistinguishable from SL signal. The design reports composition, but no agent is
  asked to check whether gains persist within identity strata for Tier A/B features too, not just for PLM/3Di. Extend the
  identity-stratified analysis to every adopted tier.
- **Feature-code licensing/data from pretrained resources.** No issue for ESM/ProstT5 weights, but record their licenses
  in DATA.md along with the datasets.
- **Disk accounting includes `models_cache`.** It is already 24 GB of the 150 GB budget. `du -sh` at milestones should
  cover the whole working directory (`du -sh .`), not only `data/`.

## 4. Things that are fine as they are

- Homology-cluster split across species with LOSO and genuinely held-out unseen species: the right structure against the
  most damaging leakage.
- Test-access log, single planned access, refit-after-test policy: correct intent. The changes above only make them
  enforceable.
- The gene-prior leakage baseline (B3), the shuffled-label control, the random-split leakage quantification, and the
  pre-registered success criterion with paired bootstrap.
- D1, which bans measured data as auxiliary targets as well as features.
- End-to-end evaluation with missed pairs at the minimum score, which blocks the "only output easy pairs" shortcut.
- Generic, non-species-specific self-exclusion in principle (item 16 asks only that it be verified).
- Fixed seed list, bench scripts with fixed sampling seeds, and the checksummed QfO/proteome READMEs.
