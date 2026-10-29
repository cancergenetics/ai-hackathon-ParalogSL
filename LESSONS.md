# LESSONS — process rules for running this project efficiently

Each rule: one line, with evidence. Pruned at every milestone (last: M3/M4, from reviews/retro/M3_retro.md; previous: M1_retro.md).

1. Write one normative spec (no "unless this file says otherwise"), dry-run every threshold, cap or split rule label-blind on the real structure, and cap the unit you actually fear; evidence: round 3 found 14 contradictions (D4); a 1,611-pair family took 18% of human pairs into test (D9).
2. Quote the pre-registered section in every brief and D-note; a recalled definition is not a spec; evidence: the final_test review found the brief's primary-population premise wrong (DESIGN §6.2) and D19 had to correct D17.
3. Forbid scoring labelled pairs in every builder/reviewer brief, give them no split membership, and list every sanctioned sealed reader in the guard at creation; evidence: D10 (identity AP on unsplit labels), D5 (G1 failed on the known-answer tests).
4. Commit code, then run production steps from a clean tree and record that hash (scripts should refuse a dirty tree); evidence: D13 and D18 (runs from uncommitted code), the Tier C meta `src_dirty=True` (re-run 00991c6), the yeast labels review finding 3.
5. One owner per file and one agent per working tree at a time; commit with `git commit <paths>`; evidence: 415356b shared-index race, D18's imported uncommitted edit, final_test review R7.
6. Stamp each artifact with a hash of its own column or parameter definitions, not of the whole registry; evidence: 239c3fd had to restamp all 12 tables, and final_test review R7 found the global hash still checked.
7. Declare rules and candidates in a commit before the run: adoption rule, secondaries, ablation sets, and per-feature missingness checks; evidence: 2ab485d, 6e3c00a, 048c11a, D17(c); the unchecked `tb_asym_*` forced a full re-run (6bb9438).
8. Run an adversarial, measure-on-real-data review before adopting any gain and before any irreversible step (one-shot test access, freeze); evidence: D17 (PARTLY ARTEFACT, later confirmed by the -0.020 within-cell test result), final_test review R1-R3 (fixed before the access), rounds 1-2 measurements.
9. Report every headline gain next to a confound-controlled version (within ascertainment cells, strata, source library, yeast-only transfer, unseen species); evidence: D12(b), D17, D20 (T = 0.297 headline vs -0.020 within cells).
10. Smoke-test new pipelines on ~5% before a full run, gate expensive tiers on a dev increment check, and record power state and load average with every timing; evidence: D11 (2.3 h run, 67% flagged), Tier C not adopted (D17b), D13 (load 95-154).
