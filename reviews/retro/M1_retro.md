# M1 retrospective (fresh process consultant)

Scope: PROMPT.md, RESEARCH_TREE.md (D0-D11), LESSONS.md, docs/AGENT_RULES.md, `git log` (42 commits), reviews/ (round1-3, results/*). No data or label files read. As of M1 no model has been scored; `results/experiments.tsv` shows no rows, so M1 evidence is about process, not science.

## 1. What should we stop doing?

1. **Writing numeric pre-registered rules without a label-blind dry run on the real data structure.** Evidence: D9 (split size cap counted supernodes, so a 1,611-pair family reached test; regenerated), D4 (14 contradictions in DESIGN v3), D6 (strict xfail asserted post-mapping counts), D5 (G1 allowlist missed sanctioned readers). Four of eleven decisions repair our own specs. Before freezing a rule, have the author run it label-blind on the real structure (pair counts per group, family sizes, caller flags) and review the output.
2. **Adopting reviewer specs "by reference" with a precedence escape clause.** DESIGN v3 said "verbatim unless this file says otherwise", and stale v2 sentences overrode the adopted specs (round 3). It cost a whole extra round (~2.2k words, plus D4). Have one agent write a single normative spec and delete the superseded text, not layer it.
3. **Running production artefacts from uncommitted code or on the shared git index.** Evidence: yeast review finding 3 (labels built at HEAD ea97950, code only committed in a79394b); commits 2c7ab4e and 415356b (shared-index race). Rule: commit the code, then run it, and record that hash. Builders commit with `git commit <paths>`.
4. **Launching full-scale runs before a small-scale smoke test.** Evidence: D11, where the human Tier A run took 2.3 h and the capped-query rule flagged 67% of human queries, so the rule is being changed after the run. Run on 5% of genes, check runtime and flag rates, then scale up.
5. **Leaving reviewers and builders room to touch protected material.** D10: the reviewer computed an identity AP on unsplit labels. The LESSONS fix is in place; stop relying on the brief alone. Give reviewers a label-view without split membership, or a script that refuses scoring.
6. **Putting hash prefixes, counts and gate bookkeeping in RESEARCH_TREE decision prose.** D8-D11 are 100-200 words each with sha prefixes (the tree is ~10 KB, and PROMPT says to keep it lean). Hashes already live in `config/protected.sha256` and `results/*.json`; the tree should hold decision, why and pointer.
7. **Timing anything under Low Power Mode.** Lessons 4 and D7: all compute-feasibility numbers are provisional, and agents cannot change it. Ask the user once to turn it off, or label timings invalid.
8. **Rule 6 (log each experiment in `results/experiments.tsv`) is currently unexercised**, so no adoption rule is on record. Before M2 modelling, create the file header and have the first row precede the first run. Otherwise stop claiming the registry exists.

## 2. What worked well (keep)

- **Reviewers who measure on real files.** Round-1/2 found the 52% single-linkage component, 49% panel saturation, missing `-a` bias, yeast testability floor (Lessons 1).
- **Independent re-implementation reviews of labels.** Yeast: 0 differences over 9.1 M pairs and an explained README gap; human: found the NBSP Cellosaurus bug and the Harle/Chan 22% fallback; both ACCEPT WITH FIXES with no label changed afterwards.
- **Pre-declared sensitivity analyses before any result** (D8 B4 ascertainment baseline AP 0.127 vs prevalence 0.027; D10 ribosomal / tie-resolved strata). These make later AUPRC interpretable.
- **Transparent incident logging** (D9, D10): the decision was made before any score existed, and the cause was stated as a spec bug.
- **Structural integrity controls:** sealed directory + G1 guard, protected hashes, quarantine of Dennler tables, label universe independent of caller (Compara), known-answer tests (12/12 human exact).
- **Label-blind feature agents** (D7) and disjoint file ownership kept parallel work flowing without leakage.
- **Decision notes with an explicit "Why"** let the tree work as memory.

## 3. Proposed LESSONS.md (replace; 9 rules, one line each)

1. Ask reviewers to measure claims on the real files, not just read the design; evidence: rounds 1/2 found the 52% giant component, 49% panel saturation, `-a` identity bias and the yeast testability floor by quick measurements.
2. Write one normative spec (no "unless this file says otherwise"), then do one consistency pass before implementation; evidence: round 3 found 14 contradictions, D4.
3. Before freezing any threshold, cap or split rule, dry-run it label-blind on the real structure and cap the unit you actually fear (pairs, not supernodes); evidence: D9 (1,611-pair family took 18% of human pairs into test).
4. Every guard or allowlist must name all sanctioned readers at creation; evidence: G1 failed on the unseen-species known-answer tests (D5).
5. Forbid scoring labelled pairs in every builder/reviewer brief AND give them no split membership; evidence: the yeast reviewer computed an identity AP on unsplit labels (D10).
6. Quarantine forbidden-column files at download time through the manifest; evidence: Dennler feature tables landed on a normal path and were caught only by round 2.
7. Commit code before running a production step, record that hash, and commit with `git commit <paths>`; evidence: 2c7ab4e/415356b index race; yeast labels built from uncommitted code (review finding 3).
8. Smoke-test new pipelines on a ~5% subset for runtime and flag rates before the full run; evidence: D11 (2.3 h human run, 67% of queries flagged, rule revised after).
9. Record machine power state with each benchmark; evidence: all feasibility timings taken in Low Power Mode and unusable.

Pruning notes: old 5 and the commit part of 7 merge into new 7; old 6+7 become new 4 and 5; old 8 generalises into new 3; old 3, 4 kept as-is (6 and 9); new 8 is added. Suggested M2 prompt for the next retro: check which of 3, 7 and 8 were violated again.
