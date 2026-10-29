# Round 3 review: consistency (statistics and data curation)

Reviewer role: fresh consistency reviewer. I read REVIEW_BRIEF.md, PROMPT.md, DESIGN.md v3, reviews/round1/PI_response.md, all of reviews/round2/*.md, and the round-1 must-lists. To check what is already being built, I also read `config/{success,gene_ids,labels_unseen,protected}.yaml|sha256` and `tests/test_data_access.py`. I did not reopen any settled choice; every item below is a contradiction, an omission, or a blocker for M1–M2. Line numbers refer to DESIGN.md v3.

## 1. Overall verdict

v3 is close to implementable, and its substance matches round 2. Most problems share one root cause. The header (l.5) adopts the round-2 specs "verbatim **unless this file says otherwise**", but several sentences carried over from v2 still say otherwise. As written, those stale sentences override the specs that the PI response claims were adopted: the label universe, the caller threshold, the yeast negative rule, and human testability. A second gap is that `reviews/round2/statistics.md` and `external_pi.md` are not in the normative list, so any accepted item that v3 did not transcribe is lost; R2-13…R2-15 are in the same position. Fifteen must-items follow. All are text-level fixes, and nearly all have to be settled before labels and splits freeze at M1.

## 2. Concrete changes

### Must

0. **Fix the precedence rule (l.5).** Replace "unless this file says otherwise" with "where this file's text and a listed normative item differ, the normative item wins unless this file names the item it overrides". Then make items 1–8 explicit.
   *Why:* without this, an implementer who reads only DESIGN.md builds the v2 rules.

1. **§3 l.62 still says the label universe is "all tested pairs from the sources".** That is the exact text that external-PI M1 and data R2-4 objected to, and it contradicts §2.3 l.29. Replace it with "U = tested pairs ∩ pinned Ensembl Compara paralogs (§2.3); caller recall is measured on U".

2. **§3 l.60 calls candidate pairs at "E ≤ 1e-3".** algorithms_comparative change 3 (must, claimed adopted) calls pairs by E′ = qlen·N_ref·2^(−bits) ≤ 1e-3, with N_ref fixed in config. Raw E depends on proteome size: yeast has about 6k proteins and a crop can exceed 40k, so the caller and the primary pair population would differ by species. Write E′ into §3, and add N_ref to `config/caller.yaml` and to the schema hash.

3. **Yeast measurement-level negative (l.36) omits "DMF > 0.5", which R2-1 requires.** At E just above 0.5, a measurement with ε ≥ −0.04 and DMF ≤ 0.5 is negative under v3 but MID under R2-1. Use the R2-1 rule verbatim.

4. **The yeast tiers are not nested, although l.39 says strict ⊂ primary ⊂ lenient.** Lenient is defined as "ε < −0.12, p < 0.05". A primary positive with sick singles can have −0.12 ≤ ε < −0.08 (R2-10; statistics 10), so lenient does not contain primary.
   - Define lenient = primary ∪ {testable, ε < −0.12, p < 0.05}, with the same negatives as primary.
   - Define strict = primary ∧ (DMF ≤ 0.2 ∨ Dean W_AB = 0).
   - Either list the Dennler-comparable tier (ε < −0.35, p < 0.05) that R2-10 asks for, or state that it is dropped.

5. **Yeast orientation conflicts (l.37).** v3 calls them "ambiguous". R2-3 calls a pair positive when n_POS ≥ n_NEG, so one positive and one negative orientation gives a positive. Choose one rule and write it out; R2-3 is the count-neutral one.

6. **Human testability (l.43).** v3 says "untestable = either single gene essential in that line". R2-5 (must) says testable iff (a + b) > −T, which is screen-scaled and explicitly needs no essential-gene list. Replace the v3 phrase with R2-5. Also state R2-5's consequence: a study whose testability cannot be computed enters the lenient tier only. This applies to Harle/Chan if their single-depletion columns are not confirmed, and to Klingbeil if it is not obtained.

7. **Human "≥2 studies" route (l.44).** The text cites R2-15, but R2-15 lies outside the normative range R2-1…12, and the rule as written has no penetrance floor. Write the floor explicitly: "≥2 studies **and** penetrance ≥ 0.10". Otherwise SL in 1/11 Ito lines plus 1/27 Harle lines makes a robust positive.

8. **Human gene ID.** §2.1 l.18 and `config/gene_ids.yaml` use the Ensembl gene ID; R2-8 (normative) uses the HGNC ID. Keep Ensembl, because it is native to Compara (R2-4), and state the override. Keep R2-8's other rules: restrict to reference-FASTA accessions, the multi-gene accession expansion, and screen symbol → HGNC (including previous and alias symbols) → Ensembl through the pinned HGNC complete set, dropping ambiguous symbols. `mapping_report.tsv` must log the losses.

9. **M1 inputs and fallbacks are missing from §2.**
   - **Inputs to add:** Ensembl Compara (pinned release; R2-4), the HGNC complete set (R2-8), the SGD GFF plus a numeric linkage window (R2-11), and YGOB/OHNOLOGS (evaluation only). Add each to §2 and DATA.md.
   - **Fallbacks:** adopt the R2-13 fallback table (Parrish S4 → ZdLFC; Kuzmin S5 → drop; Klingbeil → excluded from primary; pombe → dropped), plus external-PI S8's time box: any source not parsed and QC'd by the M1 freeze is dropped and documented.
   - **Record two decisions already made in config but missing from §2.3:**
     - Arabidopsis negatives are re-checked with Cusack's own lethality proxy, because Lloyd & Meinke is not on disk.
     - The C. elegans set is blocked until Tischler's Table 4 positives exist (`labels_unseen.yaml` status: partial).

10. **Split and power knobs must be numeric before the first split run (statistics 2–3, agent B9).** v3 says "pre-registered" but gives no values, and the two round-2 reviews disagree on the Leiden trigger (statistics: cap = max(5% of positives, 5% of pairs) per group; external PI: Leiden only if one component holds > 20% of positives). Fix in `config/splits.yaml`:
    - the edge thresholds;
    - the cap;
    - the γ grid, the seed (first seed in `config/seeds.yaml`) and the edge weight bitscore/min(self-bitscore);
    - the rule that a supernode still above the cap goes to development.

    Fix the fraction search too:
    - a fixed grid of fractions in [0.30, 0.40];
    - **one** deterministic draw per fraction with the first seed, and no redraws (otherwise this is split shopping);
    - the power simulation runs on exactly that draw.

    Commit the simulation inputs in a protected `config/mde.yaml`: binormal scores, B1 AUROC 0.65 for human and 0.70 for yeast, ρ = 0.8 as the deciding value (0.7 and 0.9 also reported), a cluster random effect, B ≥ 1,000 and ≥ 500 simulated datasets.

11. **No fallback if yeast cannot reach 30 test positives.** At a test share of 0.40, at least 75 yeast primary positives are needed. The DMF ≤ 0.5 rule, the testability floor, null alleles only, the Compara universe and the loss of linked pairs all shrink the yeast set (statistics 9; external PI S7a). Pre-register what happens if fewer than 75 exist: use 0.40 and declare the yeast comparison underpowered. Make G6 assert this rule instead of the bare count, or G6 can never pass.

12. **The bootstrap and the success rule are incomplete (statistics must 6 and should 12; external PI S10).** l.109 stratifies only by species composition. Add the "group contains ≥ 1 positive" axis, and log and drop resamples that have no positives in a species. Fix the CI type (percentile, with BCa as a check). State that the 0.20 relative-gain clause is judged on the **point estimate**. `config/success.yaml` is already protected and carries none of this, so amend it now with a decision note, before any split exists.

13. **Sealed-path layout contradicts G1/G3.**
    - **Where the file lives:** §6.3 l.115 puts "the split-assignment file" in `sealed/`. B1/G3 place `data/labels/splits/groups.tsv` (group → fold only, no labels, protected) outside `sealed/`. Development CV needs the fold assignment, so follow B1 and change l.115.
    - **Allowlist gaps:** the implemented allowlist in `tests/test_data_access.py` lacks `scripts/labels/make_splits.py`, `scripts/eval/power_sim.py` (or a statement that `make_splits` calls it in-process) and `tests/test_splits.py`. All three must read sealed test-label counts: G6 asserts yeast ≥ 30 test positives, and the power simulation needs real test group sizes and positive counts. Today the M1 split step cannot be written without widening the guard ad hoc.

14. **Round-1 must-item lost: WGD-proxy validation.** This is prior-work M3(b) and comparative M6, accepted in round 1 and re-raised as external-PI S6. v3 has the WGD/SSD *strata* but not the evaluation-only check of how well the Tier-B depth and 1:2-pattern features separate YGOB/OHNOLOGS ohnologs from SSD pairs (development data only). Add it to §6.5 as an M2 check; it is the cheapest test that Tier B works at all.

### Should

15. **G5 (adopted verbatim) asserts a "0.10 MDE" constant that v3 abolished.** State that G5 reads the v3 constants: `t_target_ap`, `power_target`, `test_fraction_range` (already in `success.yaml`).
16. **G8's forbidden-name regex matches the planned feature names.** `complex` matches the Tier-A "low-complexity fraction" (l.71), and `ppi` matches words like "overlapping". Anchor the regex tokens, or name features to avoid them, before the M2 registry exists.
17. **§4.6 (l.76) is v2 text.** Replace it with change 10: the outgroup is the RBH in the nearest *post* species, distances are computed on shared aligned columns, Kimura distances are used, and the median over the 3 nearest is taken. Also replace the `fragment_suspect` definition (l.71) with change 11: the median length of A's best hits in the retained panel. That makes it a panel feature, not a Tier-A one.
18. **Quantify the masking adoption rule (l.79).** "A drop that masking recovers" needs numbers. Use statistics 16 / external-PI S1b: the thinning drop exceeds 1 SE, masking recovers more than 1 SE, and natural-panel T is not worse by more than 1 SE.
19. **Tier C cannot be the default under §7's own numbers.** Self-search (10 min) + panel search (45–70 min) + ESM-2 (~100 min) is more than 2 h. Pre-declare that the single primary test access scores the default-tier model and that any opt-in `--plm` model is reported as secondary in the same access. Otherwise there are two "final models".
20. **Missing specifications for secondary analyses.**
    - End-to-end B1 on caller-missed pairs: direct pairwise alignment, with no alignment scoring 0 (statistics 5).
    - LOSO with joint groups: train on species X in folds ≠ k, test on species Y in fold k (statistics 14).
    - The SE used by the one-SE rule: cluster-bootstrap SD of pooled out-of-fold T (statistics 15).
    - The stop rule (l.92) uses AP but selection uses ΔAP vs B1. Use ΔAP for both.
21. **The power simulation at M1 uses labelled-universe counts, but the primary population is caller-called pairs,** and caller thresholds may change once at M2. Freeze the E′ threshold at M1 and run the simulation on called pairs, or state the recall assumption.
22. **Determinism across thread counts (l.127) needs the tolerances from agent B8.** As written, the M2 invariance test is bound to fail. **Memory (l.135):** `/usr/bin/time -l` reports the peak RSS of a single process, not the MMseqs2 child tree or MPS allocations (algorithms round-1 22). Use a process-tree sampler, because the human memory figure is a PROMPT deliverable.
23. **Restore the omitted should-items:**
    - `config/report_plan.yaml` (statistics 11): a CI only for strata with ≥ 10 positives, interaction CIs, post-hoc labels;
    - n_lines strata if the ascertainment control is positive (statistics 13);
    - the label-noise ceiling (statistics 14, external-PI S9, data D16);
    - SMF and E-bin strata for yeast (statistics 9).
24. **Reviews (l.130) say "re-run from a clean checkout" for every key result.** That conflicts with B7/G14: a worktree on cached, hashed inputs, with a full from-raw re-run once at M4. Align the text.
25. **Extend the protected list (l.142)** with `config/{mde,panel,splits,caller,gene_ids,labels_unseen}.yaml`.

### Could

26. **`near_identical` (l.10):** add "over ≥ 95% of the shorter protein" (comparative M5).
27. **Panel wording (l.22).** "Units already in the QfO core (incl. test-species congeners …)" reads as if *A. lyrata*, *Capsella* and *Caenorhabditis* are in the core; it is their units (Brassicales, Rhabditidae) that are. Also record the Riccia vs Marchantia choice.
28. **G6's [0.25, 0.40] vs §6.1's [0.30, 0.40]:** state that the first is a tolerance. **Dean 197/48 vs 201/49:** state that G4's strict xfail governs, not R2-9's "fail loudly".
29. **Label Flister's ≈ 0.5% prevalence (l.112) "assumed"** until it is verified (agent D4).

## 3. Risks not addressed

- **Yeast primary power.** Item 11 is likely to bind. The report should be prepared for "not demonstrated in yeast".
- **Precedence ambiguity (item 0) is the main implementation risk.** Sub-agents receive DESIGN.md, not the review files.
- **The human harmonised LFCs** (in4mer `Supp_table1`) are a third-party re-analysis with unverified timepoints (data risk 1). Testability, and so negatives, depend on them.
- **The training-time panel contains the unseen test species' proteomes** (A. thaliana, S. pombe and C. elegans are in the QfO core). This is sequences only and allowed, but §2.3's "never used for training" should say "labels never used".

## 4. Fine as it is

- The primary criterion (species-mean ΔAP vs B1, CI excludes 0, ΔAP > 0 in each species), with no nested-CV fallback.
- Same-pair scoring of B1, enforced in code.
- Contraction of labelled pairs into supernodes, with connected components as the default.
- The leakage dose–response.
- The masking statistical conditions (shared folds, pair weights, leaf size in pairs, no row-random early stopping, one row per pair at evaluation).
- The ban on within-proteome-constant features, plus the species-identification control.
- Per-species panel search with backtrace, core-set D_s, loss-aware duplication classification, and the written panel rule.
- Self-exclusion wording, the sealed test, the tag, G15, and the refit that is never scored.
- The commit and quarantine policy.
