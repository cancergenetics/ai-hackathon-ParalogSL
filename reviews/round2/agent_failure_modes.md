# Round 2 review: AI-agent failure modes

Reviewer role: auditor of the ways an autonomous PI agent and its sub-agents could fail while running this design.
Inputs read: reviews/REVIEW_BRIEF.md, PROMPT.md, DESIGN.md v2 (all of it, closely §2.3, §4, §6, §8, §9), docs/AGENT_RULES.md,
reviews/round1/agent_failure_modes.md, reviews/round1/PI_response.md, RESEARCH_TREE.md, DATA.md, LESSONS.md, .gitignore,
`git log` / `git status` / `git ls-files`, config/data_manifest.yaml and .lock.json, scripts/data/download.py (header),
listings of data/raw/, data/raw/_quarantine/, data/raw/_survey/, results/, scripts/.
Nothing was run beyond `ls`, `grep`, `git` and `du`. Note that a data download was running while I read (empty
`data/raw/<group>/` directories were created at 23:47), so on-disk observations are a snapshot.

## 1. Overall verdict

v2 accepts all 19 of my round-1 must items in substance. §6.3, §8 and AGENT_RULES.md now name the right mechanisms
(sealed directory, tag-before-test, protected hashes, feature registry, gene-ID scanner, experiment registry, review
records). Two kinds of gap remain. First, several mechanisms are still described as intentions with no assertion
attached: "only the PI runs it", "adoption rule declared before the run", "re-checked from a clean checkout",
"thresholds adjusted once". Second, a few new holes have appeared since v1. The data manifest will download the full
Dennler & Ryan tables, including a ready-made similarity-feature table, into `data/raw/dennler2025/` rather than the
quarantine. The raw unseen-species sources are not sealed. The pre-registered nested-CV fallback would have to read
the sealed labels outside `final_test.py`. Three of my own round-1 proposals were over-engineered or wrong, and I
withdraw or shrink them here. The design is ready for M1 once the 10 M1 gates in §2.C below exist as code. None of
them needs more than about half a day.

## 2. Concrete changes

### A. Status of my round-1 MUST items in v2

| R1 # | Item | v2 status | Remaining gap (→ item in §2.B/C) |
|---|---|---|---|
| 1 | Sealed test + path-guard test | **Fixed in text** (§6.3) | The guard covers only `data/labels/sealed/`. Raw unseen-species sources and the split-assignment file are readable. The dev loader has no runtime refusal. → B1, G1 |
| 2 | Freeze + tag before test | **Mostly fixed** (§6.3) | The log lacks the model SHA-256, the sealed-manifest SHA and the purpose. Nothing makes each access be committed. "One planned access" is not counted. → B4 |
| 3 | Re-run loophole | **Fixed** (§6.3) | The rule does not say that the bug reviewer must be a different agent from the author. Minor. → B4 |
| 4 | Refit never scored on test | **Fixed in text** (§5) | No mechanism. I withdraw my round-1 proposal of training-pair hashes in every metric JSON as over-engineered. A model-SHA check in `final_test.py` is enough. → B4 |
| 6 | Protected-file hashes | **Fixed** (§8) | The list omits `final_test.py`, the controls, the success-criterion constants, the guard tests themselves, and the pytest configuration, which can silently disable tests. → G3 |
| 7 | Thresholds frozen before label code | **Fixed** (§2.3) | It is not yet a checkable ordering. Also, "Thresholds may be adjusted once" (§3, caller) and the "borderline margin" (§2.3) need numbers in config before first use. → G7 |
| 8 | Known-answer label tests | **Fixed** (§2.3) | No tolerances are stated. → G4 |
| 9 | Pinned metric | **Fixed** (§6.2) | The toy test and the requirement that B1 and the model are scored on the same pair set are not written down. → G5 |
| 12 | Quarantine + registry | **Partly fixed.** Two 30 KB head files were quarantined. | **Regressed for the full files:** the manifest routes all four Dennler files to `data/raw/dennler2025/`. → B2, G2 |
| 13 | Gene-name scanner | **Accepted** (§8) | As I specified it (HGNC symbol regex) it would produce too many false positives. I shrink it. → G9 |
| 14 | Masked gene IDs for feature agents | **Partly accepted** (feature cards, label-blind feature code) | I accept the PI's modification. Feature cards become checkable through git ancestry. → G13, C-could |
| 17 | Training features via the CLI | **Fixed in text** (§4) | No parity test. My round-1 test (subsample 500 proteins and compare with full-proteome training features) was **wrong**: proteome-percentile features legitimately differ on a subsample. Corrected version → G10 |
| 20 | Experiment registry + pre-declared rules | **Fixed in text** (§5) | "Declared before the run" cannot be checked as written. REPORT does not state the number of configurations tried. → G13 |
| 23 | Report numbers from results | **Accepted** (§8) | "Flag any decimal" is fragile. Replace it with a template lint. → B6 |
| 24 | Review record format | **Partly fixed** (§6.5, `reviews/results/`) | No required fields and no checker. A clean-checkout re-run is infeasible for multi-hour steps because the data are gitignored. → G14, B5 |
| 25 | AGENT_RULES preamble | **Done** (docs/AGENT_RULES.md) | It misses the most common agent hacks: editing a failing test or tolerance, skip/xfail, git history rewriting, editing the tree, and reading raw unseen-species sources. → B3 |
| 27 | DATA.md backfill + M1 gate | **Accepted, not done.** DATA.md is still empty. `update_data_md.py` is named in the manifest header but does not exist. 20 of the 25 `_survey/` files are not in the manifest. | → B2, G2 |
| 28 | Pinned environment | **Accepted, not done** (fine for M1) | → G7 asserts it at M1 close |
| 29 | Commit cited bench outputs | **Accepted, not done.** `results/` is empty, and docs/lit/compute_feasibility.md still cites gitignored `logs/bench/*.jsonl`. | → C2 |

### B. Changes to the design text (all small)

1. **MUST: seal the raw unseen-species sources and the split assignment, and add a runtime refusal.** The sealed
   prefixes become `data/labels/sealed/`, `data/raw/{cusack2021,tischler2006,pombe_gi,ryan2012_pombe}/`, any
   `flybase_*` file, and `data/raw/_quarantine/`. The only files that may name these prefixes are
   `scripts/labels/build_unseen.py`, `scripts/labels/make_splits.py`, `scripts/labels/crosscheck_dennler.py`,
   `scripts/eval/final_test.py`, `scripts/data/*`, and their own tests. `slparalog.data.load_labels()` calls
   `Path.resolve()` and raises on any sealed prefix, so a symlink cannot bypass it. The split-assignment file
   (`data/labels/splits/groups.tsv`) holds group→fold only, with no labels, and is a protected file.
   *Why:* today a feature agent can open `_survey/cusack2021_Dataset_1.txt` or rebuild test labels without touching
   `data/labels/sealed/`. State in §6.3 that this guard is a tripwire against accidental access. Deliberate access is
   deterred by the access log and the git audit, not prevented.

2. **MUST: fix the manifest so forbidden-column files land in quarantine, and clean up `_survey/`.** In
   `config/data_manifest.yaml`, set `group: "_quarantine/dennler2025"` for all four Dennler entries. That includes
   `ens111_human_similPerc.csv`, a ready-made Ensembl-derived similarity-feature table, the exact thing an agent would
   "just merge". Do this **before the running download finishes**, or move the files afterwards and re-hash them. For
   `_survey/`: delete every file whose SHA-256 matches a manifest-downloaded copy, register the rest (Costanzo readme,
   Dean PDF, Flister preprint tables, FlyBase, chou2025 tables) as manifest entries, then remove `_survey/`.
   *Why:* DESIGN §2.3 promises quarantine, but the manifest contradicts it. Two unregistered copies of every label source
   invite a script to read the copy that has no provenance.

3. **MUST: add six rules to docs/AGENT_RULES.md**, keeping it under one page:
   (a) "If a test, tolerance or known-answer check fails, report it. Never edit the test, its tolerance, its fixture,
   `conftest.py`/pytest config, or add `skip`/`xfail`." (b) "Never run `git tag -f`, `commit --amend` on a pushed or
   reviewed commit, `reset --hard`, `rebase`, `--no-verify`, or edit `results/test_access_log.tsv`." (c) "Do not edit
   RESEARCH_TREE.md, DESIGN.md or `config/protected.sha256`. Only the PI does." (d) "Do not open the sealed prefixes in
   B1. That includes the raw Cusack/Tischler/pombe/FlyBase files." (e) "Every number in your report cites
   `path:json-key` and the commit SHA, or is marked 'unverified'." (f) "Do not use web search or your own knowledge to
   decide features or thresholds for specific genes or families."
   *Why:* sub-agents optimise for the task looking done. Editing the check is the cheapest route, and today the preamble
   does not forbid it.

4. **MUST: tighten §6.3 so that `final_test.py` checks its own preconditions**, listed as G15 below. It also writes
   `model_sha256`, `sealed_manifest_sha256`, `purpose` ∈ {primary, review-rerun, bugfix-rerun}, and the agent identity
   to the log. "One planned access" becomes: exactly one `purpose=primary` row. A `bugfix-rerun` row requires
   `reviews/results/bugfix-<n>.md` by a non-author.
   *Why:* "only the PI runs it" cannot be enforced. Preconditions inside the script can be.

5. **MUST: make the nested-CV co-primary (§6.1, MDE > 0.10 rule) run inside `final_test.py`, after the tag.** It uses
   the frozen pipeline, including the frozen hyper-parameter-selection procedure, inside inner folds.
   *Why:* nested CV "over all labelled data" reads the sealed labels. If it is written as a separate development
   script, either the seal guard fails or someone widens its allowlist. That is exactly the silent weakening the guard
   exists to prevent.

6. **SHOULD: replace "check_report_numbers flags any decimal" with a template lint.** REPORT.md is rendered from
   `report/REPORT.md.j2` by `scripts/report/render.py`. `tests/test_report_template.py` asserts that the template
   contains no `\d+\.\d+` or percentage outside `{{ … }}`, except on lines tagged `<!-- lit:cite -->`, and that the
   committed REPORT.md equals `render(template, results/)`. Drop the number check for RESEARCH_TREE.md. It is a
   working memory, and a gate there costs more than it is worth. Keep the rule "cite path:key" (B3e).
   *Why:* matching free-text decimals to JSON values fails on rounding and on counts, and agents then add allowlist
   entries until it passes.

7. **SHOULD: redefine "independent review with a clean-checkout re-run" as something that can actually be done.** The
   reviewer creates `git worktree add /tmp/rev <sha>` with `data/raw` and the cache symlinked read-only. The reviewer
   verifies the SHA-256 of the cached inputs, re-runs the step that produces the number from those inputs
   (evaluation from cached features, not re-embedding), and compares within a tolerance stated in the record. A full
   from-raw re-run is required once, at M4, for the headline test number and the human runtime.
   *Why:* data and embeddings are gitignored and take hours to rebuild. An impossible requirement ends up as a
   reviewer "agreeing with the summary".

8. **SHOULD: state determinism tolerances before the tests are written.** §6.5 says "determinism across thread counts".
   LightGBM with `deterministic=True` is reproducible at a fixed `num_threads`, but not bit-identical across thread
   counts. MPS fp16 embeddings are not bit-identical run to run. Pre-register: bit-identical at fixed threads and fixed
   inputs (caller, Tier A/B, model). Across thread counts, scores must agree within |Δ| ≤ 1e-6 or Spearman ≥ 0.9999.
   Embeddings must agree within cosine ≥ 0.9999. The model bundle records `num_threads`.
   *Why:* a test that is bound to fail is the test an agent will "fix" by loosening it. Stating the tolerance first
   removes the temptation.

9. **SHOULD: fix the split knobs numerically in `config/splits.yaml` before the first run.** These are the component
   cap, the Leiden resolution and seed, the "very large group → development" threshold, and the yeast test-positive
   floor (≥30). They are protected after the first build.
   *Why:* every one of them is a lever for split shopping once the test fold composition is visible.

10. **SHOULD: write the negative-result commitment into §10.** The first line of REPORT.md is the primary-criterion
    pass/fail, rendered from `results/final_test.json`. REPORT.md states the number of rows in
    `results/experiments.tsv`. (Round-1 item 22 was not adopted. I drop the K-failure stopping rule because PROMPT.md
    sets its own stopping condition.)

### C. Minimal checkable gates for M1 and M2

Principle: one fast `pytest` suite (`tests/`) for checks that need only code and small committed files, plus a `data`
marker for checks that need `data/raw`. Checks that read git history live in `scripts/milestone_check.py`. Git hooks
are not relied on, because they are unversioned and `--no-verify` bypasses them. A milestone is closed only when
`scripts/milestone_check.py --milestone M<n>` exits 0 and its JSON output is committed.

**M1 (before any label is built, and at M1 close)**

| Gate | File | Asserts |
|---|---|---|
| G1 | `tests/test_data_access.py` | (i) AST string constants and plain text of every `*.py`/`*.sh`/`*.ipynb` under `src/`, `scripts/`, `tests/` contain none of the substrings `sealed`, `_quarantine`, `cusack`, `tischler`, `pombe_gi`, `ryan2012`, `flybase`, `dennler`, except in the allowlisted files of B1, which are hard-coded in this test. (ii) `load_labels()` raises `PermissionError` on a sealed path, including through a symlink made in `tmp_path`. (iii) No file under `src/slparalog/features/` contains `data/labels`, `data/raw`, `ygob`, `ohnolog` or `keyword`, so evaluation-only annotations cannot become features. |
| G2 | `tests/test_data_registry.py` (marker `data`) | Every file under `data/raw/` maps to a manifest entry, and its SHA-256 equals the lock or `SOURCE.json`. No `_survey/` directory exists. Every manifest entry whose `used_for` contains "never features" or "cross-check" has a group under `_quarantine/`. DATA.md equals the output of `scripts/data/update_data_md.py` (regenerate and diff). Every entry has non-empty `license`, `version` and `committed` fields. |
| G3 | `tests/test_protected.py` | For each line of `config/protected.sha256`, the current hash matches, or RESEARCH_TREE.md contains the first 12 hex characters of the new hash inside a `D<n>` decision note. The protected list must at least contain: `config/{labels,splits,seeds,success,caller}.yaml`, `src/slparalog/eval/metrics.py`, `scripts/eval/{final_test,leakage_audit,controls}.py`, `scripts/labels/make_splits.py`, `data/labels/splits/groups.tsv`, `data/labels/sealed/MANIFEST.sha256`, `tests/test_{data_access,protected,metrics,label_known_answers,splits}.py`, `tests/conftest.py`, and `pyproject.toml`, because `[tool.pytest]` `addopts` can deselect tests. The test also asserts that its own list contains these entries, so entries cannot be dropped silently. |
| G4 | `tests/test_label_known_answers.py` (marker `data`) | Published counts are reproduced from the primary tables, before harmonisation, with tolerances written in the test now: Dede 24 SL pairs (exact); Harle 117 / 355 (exact); Costanzo ε<−0.12, p<0.05 pair count (±1%); Cusack 51/110 strict (exact; this test sits in the B1 allowlist). A known count that cannot be reproduced is listed in `results/labels/known_answer_gaps.json` as `xfail(strict=True)` with a reviewed note. It is never deleted. |
| G5 | `tests/test_metrics.py` | `metrics.ap([1,0,1,0],[.9,.8,.7,.1]) == 5/6`. Constant scores give AP = prevalence. AP is invariant to row order. `ap()` refuses mismatched pair-ID sets when given (model, B1) frames, so B1 is always scored on the same pairs. `sklearn.__version__` equals the version in `requirements-lock.txt`. The success-criterion constants (0.20 relative gain, 10,000 resamples, 0.10 MDE, seed 11) are read from `config/success.yaml`, not from literals in code. |
| G6 | `tests/test_splits.py` (marker `data`) | Re-running `make_splits.py` on the frozen inputs reproduces `groups.tsv` byte for byte. There is zero gene overlap between test and development. Every test pair's genes fall in test groups. The per-species test-positive fraction lies in [0.25, 0.40] and yeast has ≥30 test positives. The file hash equals the protected hash. |
| G7 | `scripts/milestone_check.py --milestone M1` | Runs `pytest -m "not slow"` including `data`, with 0 failures. Collected-test count ≥ the count in the previous `results/milestones/*.json`. No `skip`/`xfail` outside `tests/skip_allowlist.txt`. `git status --porcelain` is empty. `du -sk` of the whole project is under 150 GB, `models_cache` included. Git ordering: the first commit of `config/labels.yaml` is an ancestor of the first commit of any `scripts/labels/build_*.py`, and the first commit of `results/mde.json` is an ancestor of the first commit of any file under `src/slparalog/features/`. `requirements-lock.txt` and `config/tools.yaml` exist and match `pip freeze` and the binaries' SHA-256. Writes `results/milestones/M1.json`. |

**M2 (tool skeleton, Tier A/B, baselines, controls)**

| Gate | File | Asserts |
|---|---|---|
| G8 | `tests/test_feature_registry.py` | The columns of every training matrix are a subset of `registry.FEATURES`. Each registry entry has `source` ∈ {self_search, panel, plm, structure, sequence}. No registered name matches `(?i)gtex\|ppi\|essential\|complex\|ceres\|depmap\|go_\|pfam\|interpro\|expr\|ygob\|ohnolog\|label\|^sl$`. The model bundle's feature-schema hash equals the hash of the registry. |
| G9 | `tests/test_no_hardcoded_ids.py` | AST string constants in `src/` and `scripts/` (excluding `scripts/labels/`, `scripts/data/`, `tests/fixtures/`) match none of these: UniProt accession; `ENS[GTP]\d{11}`; `Y[A-P][LR]\d{3}[WC](-[A-Z])?`; `AT[1-5CM]G\d{5}`; `WBGene\d{8}`; `SP[ABC]C\w+\.\d+c?`; `FBgn\d{7}`; or, as a whole token, any symbol of length ≥4 in `data/labels/gene_universe.tsv` (a gene list with no labels, written by the label builder, not sealed), minus `tests/id_allowlist.txt`. This replaces the HGNC-wide scan from round 1, which flags `CAT`, `MAX`, `SET` and similar words. |
| G10 | `tests/test_train_infer_parity.py` | On the toy proteome and toy panel, `slparalog features toy.fasta` (run as a subprocess through the console script) and the training-matrix builder give identical Tier A/B feature values for every pair. Integers must be exact; floats must agree within 1e-9. Panel-mask k = 0 only: the masking augmentation is a deliberate train-only path. The training pair set is a subset of the caller's pair set: no forced alignments of labelled pairs the caller missed. |
| G11 | `tests/test_invariance.py` | Already in §6.5: shuffling input order, renaming IDs and swapping A/B give identical outputs. Toy self-exclusion has a positive control (a panel containing the query drops it) and a negative control (a distant species is kept). |
| G12 | `tests/test_controls.py` (toy, plus `data` for dev) | Shuffled labels give an AP inside the seeded bootstrap 95% interval of prevalence. **Should:** a planted-signal label (a known function of two features at real prevalence) gives AP ≥ 0.9 × oracle, and identity-only run through the full GBM pipeline matches B1 within its CI. |
| G13 | `tests/test_experiment_registry.py` | Every row of `results/experiments.tsv` has a non-empty `adoption_rule`, an existing `feature_card` path and a `git_sha` that resolves. The commit that first added the card, with its adoption rule, is an ancestor of the row's `git_sha`, and the card's rule text has not changed since. `results/experiments_summary.json` holds the row count used by REPORT. |
| G14 | `tests/test_review_records.py` | Every `[x]`/`[!]` node in RESEARCH_TREE.md from item 4 onward links to an existing `reviews/results/<id>.md`. That file contains the fields `reviewer:`, `producer:` (different from the reviewer), `commit:` (resolves), `rerun_command:`, `inputs_sha256:` and `match:` (yes/no plus tolerance). |

**M4 (design now, implement before the freeze)**

| Gate | File | Asserts |
|---|---|---|
| G15 | `scripts/eval/final_test.py` (self-checks) and `tests/test_final_test_guard.py` (on a fake sealed set) | Refuses unless all of these hold. HEAD is the commit of tag `final-model-v1`, or a descendant whose diff touches only `results/`. The tree is clean. The model file SHA equals `results/frozen/model.sha256` at the tag. Every sealed file is listed in `sealed/MANIFEST.sha256`, and its first commit predates the tag. The frozen config contains the chosen B1 variant and the hyper-parameters. `results/test_access_log.tsv` has no uncommitted lines, so each access must be committed before the next. At most one `purpose=primary` row exists. The guard test runs the script against a fake sealed set and checks that each refusal fires. |

What I would **not** build (over-engineered relative to value):
- Masked or anonymised gene IDs for every feature agent. I accept the PI's position. *Could:* the error-analysis
  script writes opaque pair IDs by default (`--names` to override), because false-negative lists are where gene priors
  enter.
- Training-pair hash lists in every metric JSON (round-1 item 4). G15's model-SHA check covers the same risk.
- A number checker for RESEARCH_TREE.md, and a free-text decimal matcher for REPORT.md (B6 replaces it).
- An HGNC-wide symbol scan (G9 replaces it).
- Pre-commit hooks as a control, since they are bypassable and unversioned. `milestone_check.py` is the gate.
- A gate on "AGENT_RULES pasted verbatim in every delegation". It cannot be checked. *Could:* sub-agent reports start
  with `AGENT_RULES@<short-sha>` so the PI can grep for it.

Still unenforceable after these gates. Accept these knowingly and name them in REPORT.md limitations:
- A determined agent reading sealed data through a constructed path. Only the audit trail deters it.
- An LLM's gene priors shaping feature *ideas* during error analysis. G13 only timestamps the idea.
- The PI writing both a protected-file change and the decision note that legitimises it. G3 makes the change visible;
  it does not make it independent. The fresh-agent milestone review required by PROMPT.md should diff
  `config/protected.sha256` since the previous milestone tag.

### D. Other should/could items

1. **SHOULD:** serialise heavy jobs (round-1 item 32 was not adopted). Wrap MMseqs2 panel searches, embedding and
   ProstT5 in `scripts/run_heavy.sh` using `flock logs/heavy.lock`, and record `pmset -g` with every timing.
   *Why:* §7's runtime numbers are a deliverable. Parallel agents contending for 24 GB will corrupt them.
2. **SHOULD:** before round 3, copy the cited bench summaries to `results/bench/` and repoint
   docs/lit/compute_feasibility.md to them. Remove the hard-coded `.` from `scripts/bench/*`
   and from the `download.py` docstring.
   *Why:* this was accepted in round 1 but `results/` is still empty.
3. **SHOULD:** commit the pending RESEARCH_TREE.md change, which records D2 and is currently uncommitted.
4. **COULD:** verify the [U] Flister ≈0.5% prevalence before using it for §6.2 reweighting (round-1 item 26), or label
   that number "assumed".

## 3. Risks the design does not address

- **The download is running against a manifest that contradicts the design (B2).** When it finishes, four Dennler
  files with forbidden columns and a ready-made feature table will sit unquarantined in `data/raw/dennler2025/`.
- **The seal and the nested-CV fallback conflict (B5).** Unless this is resolved now, the first agent to implement the
  MDE fallback will widen the guard's allowlist.
- **Guard tests are themselves editable.** G3 protects them, and AGENT_RULES rule (a) forbids editing them. Without
  both, an agent will meet a red guard test by adjusting the guard, which is the most common reward hack in agent
  coding loops.
- **Train-only code paths multiply in v2:** masking augmentation k ∈ {1,2,3}, species-balanced weights, and
  percentile features over the called-pair set. G10 pins k = 0 to the CLI, but a bug in the k > 0 paths would only
  show up as a dev/test discrepancy. The panel-thinning test (§4) should therefore run through the CLI with a thinned
  panel, not through the training-time mask code.
- **Caller-threshold adjustment "once on development labels" (§3)** is a tuning step on labels that the end-to-end
  metric then rewards. Record it as an experiment row (G13) and protect `config/caller.yaml` after M2 (G3).
- **Seeing unseen-species counts in LITERATURE.md** is harmless. But the Cusack known-answer test (G4) means a test file
  loads Cusack. Keep that test in the B1 allowlist and out of feature agents' tasks.

## 4. Things that are fine as they are

- §6.3's structure: sealed directory, tag, access log, one access, re-run rule, and the refit never scored. With B4 and
  G15 it becomes enforceable.
- §2.3: thresholds pre-registered in `config/labels.yaml`, never tuned on the 13-pair list, on agreement with Dennler,
  or on model results. Unseen-species endpoints stated per set. The label universe independent of the caller.
- §4's statement that training features are produced by the shipped CLI on whole proteomes, plus the feature registry.
- §5's experiment registry, one-SE rule and single final model. §6.2's pinned AP implementation and pre-registered
  success criterion.
- The quarantine of the two Dennler head files (done) and docs/AGENT_RULES.md as a short, pasteable preamble. It needs
  only the B3 additions.
- The PI's partial rejection of anonymised IDs. Its reasoning holds for the feature *code*. The residual risk lies in
  error analysis and is handled by G13 plus the C-could item.
