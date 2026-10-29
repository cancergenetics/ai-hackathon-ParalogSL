# Round 1 review: external PI

Reviewer role: external PI (senior computational biologist). I looked at priorities, the order of work, scope fit, what is over-built or missing, and whether the success criterion and deliverables match PROMPT.md.
Read: PROMPT.md, DESIGN.md v1, LITERATURE.md, docs/lit/{methods,human_labels,yeast_and_other_species_labels,compute_feasibility}.md, data/raw/qfo/species_list.tsv.

## 1. Overall verdict

The design is sound and unusually disciplined. Its labels come from primary tables with tested negatives, it splits by homology across species, it runs leave-one-species-out (LOSO), it controls test access, and it uses the right baselines. I would fund it.

Its weak spots are not in machine learning but in the *deployment setting* the prompt cares about: a newly sequenced plant or insect proteome. The design does not test robustness to messy input proteomes (allelic or isoform duplicates, fragments). Its panel features depend on how densely the panel samples the query's relatives, and that density varies a great deal (human sees macaque; Arabidopsis sees only rice, maize and moss). Unseen-species performance carries almost no weight in the success criterion.

Some label rules are inconsistent between species, in ways that will set the model's ceiling. They must be fixed before the labels are frozen, because they cannot be fixed afterwards.

The expensive tiers are 3Di and per-gene embedding PCs. Their expected payoff is low, so they should be gated or cut.

## 2. Concrete changes

### Labels: these must be settled before the M1 freeze

1. **MUST: make the yeast negatives consistent with the yeast positive definition.**
   - The primary positive includes strong synthetic *sick* pairs (ε < −0.35).
   - The design adds Kuzmin 2020 and Dean 2008 *viable* double mutants as tested negatives. A viable double mutant can still have ε ≪ 0.
   - Fix: compute a multiplicative ε from their single-mutant and double-mutant fitness (Kuzmin Table S5; Dean W values), then apply the same SGA thresholds: positive if ε < −0.35; negative only if the ε that would make it a negative is not met. If ε cannot be computed, use those pairs only in a strict "lethal vs viable" tier.
   - *Why:* otherwise the negatives contain true positives by construction.

2. **MUST: treat essential-gene pairs the same way in both species.**
   - Human labels mark a pair "untestable" when either gene is singly essential.
   - Yeast labels include pairs tested with TS/DAmP alleles of essential genes. These are hypomorph interactions, not SL.
   - Fix: exclude from primary yeast labels any pair whose member is essential as a single deletion. Report those pairs as a separate "essential-member" stratum.
   - *Why:* SL is undefined when a single mutant is already lethal. Mixing the two definitions teaches the model different targets in each species.

3. **MUST: add a penetrance clause to the human positive rule** (the label survey §4.3 recommended it, but it was dropped).
   - Change "SL in ≥2 testable lines (any studies)" to "SL in ≥2 testable lines **and** (≥2 independent studies **or** penetrance ≥ 0.25)".
   - Keep "number of lines tested" as a QC column.
   - *Why:* pairs tested in Klingbeil's 22 or Harle's 27 lines get many chances to pass "≥2 lines" by noise alone, while Dede/Thompson pairs get 2–3 chances. Under the current rule, positive status becomes a function of which library a pair was in, and library choice correlates with identity and family size.

4. **SHOULD: Arabidopsis labels.**
   - Require viable single mutants in *negatives*, not only in positives. Cusack's "non-redundant" class includes pairs where a single mutant already has the phenotype.
   - Add GABI-DUPLO (Bolle 2013). It is a systematic set of about 200 double knockouts with about 29 lethal, and its design is much less ascertainment-biased than the literature-curated Lloyd & Meinke positives. Deduplicate it against Cusack.
   - *Why:* Arabidopsis is the only genuinely unseen-lineage test, and right now it is small (51/110) and positive-biased.

5. **SHOULD: timebox label acquisition and close the source list at the freeze.**
   - Any source not parsed and QC'd by the M1 freeze is dropped and documented. It is not added after modelling starts.
   - Start the parsers with the sources already in `data/raw/_survey/` (Dede, Thompson, CHyMErA, in4mer/Chou) plus Ito. Klingbeil, Harle/Chan and Flister follow.
   - *Why:* labels must be frozen before fitting, and ten heterogeneous sources is where the schedule will slip.

6. **SHOULD: make obtaining Flister 2025's full per-pair table a priority.**
   - If it is obtainable, use its held-out-cluster pairs as a co-reported "near-genome-wide human" test. Positive = SL in both lines; negative = tested in both lines and SL in neither.
   - *Why:* it is the only human sample close to the tool's real use case (all paralog pairs, prevalence < 0.5%). Every other human label lives in the screened-pair space, which is 89% closest pairs and 73% WGD.

### Evaluation and success criterion

7. **MUST: put unseen-species performance into the pre-registered success criterion.** PROMPT.md asks for good performance "on both species and on unseen species". Currently only "Δ > 0 point estimate on Arabidopsis" counts. Add:
   - (a) LOSO, yeast→human and human→yeast: the model beats B1 with the paired-bootstrap CI excluding 0 in at least one direction, and Δ > 0 in both.
   - (b) Arabidopsis and C. elegans: report lift over prevalence and Δ vs B1 with CIs. State in advance that these sets are underpowered and that a CI including 0 is "not demonstrated", not "failed".
   - *Why:* the criterion should measure the stated goal, not only within-species held-out clusters.

8. **MUST: run a power simulation before freezing the split, and choose the holdout design from it.**
   - Using the label counts once they are computed, simulate the cluster-bootstrap CI half-width of ΔAUPRC for a 20% holdout.
   - If the expected human test set is about 60–120 positives and the CI half-width exceeds the plausible effect (about 0.03–0.06 AP), switch to one of:
     - a larger holdout (30%), or
     - cluster-grouped nested CV, with the outer folds evaluated once through `final_test.py` under the same access log.
   - *Why:* a criterion that is unachievable for statistical reasons, not methodological ones, wastes the project's one test access.

9. **MUST: audit the single-linkage cluster split before freezing it.**
   - Single linkage at E ≤ 1e-3 and coverage ≥ 0.3 over human ∪ yeast will very likely produce a giant component (kinases, small GTPases, WD40/ankyrin-repeat proteins).
   - Forcing clusters larger than 5% into development would then push most of the large families out of the test set, so the test would cover only small families.
   - Fix:
     - Report the cluster-size distribution and the fraction of labelled positives in the top clusters.
     - If one component holds more than about 20% of positives, re-cluster at a stricter linkage (e.g. E ≤ 1e-10 and bidirectional coverage ≥ 0.5) and run a residual-similarity audit: report test results stratified by maximum identity of each test gene to any training gene, following Hamp & Rost.
   - *Why:* the test set must represent the family-size spectrum, because family-context features are the main hope for beating identity.

10. **SHOULD: select models on a scale-free cross-species score.** Use the mean over species of ΔAUPRC vs B1, or of log-lift, not mean AUPRC. *Why:* yeast prevalence (about 20–35%) is several times human prevalence (about 5–10%), so mean AUPRC is dominated by yeast.

11. **SHOULD: add a label-noise ceiling.**
    - On pairs tested by two independent studies, compute the AUPRC of one study's calls (or scores) at predicting the other's.
    - Do the same for SGA orientation A×B vs B×A.
    - *Why:* the replicate ceiling (precision ≈ 0.40 / recall ≈ 0.32) is needed to interpret any AUPRC the model reaches, and it is cheap.

12. **SHOULD: report an accuracy-vs-runtime table for each feature tier.** Report AUPRC per species with human wall time and peak RSS for Tier A, A+B, A+B+C and A+B+C+D. Benchmark with Low Power Mode off. *Why:* PROMPT.md says "report runtime alongside accuracy", and the tier choice should come from that trade-off.

### Deployment realism: the main gap

13. **MUST: guard against allelic and isoform duplicates, and test robustness to input quality.** New plant and insect assemblies often contain haplotigs (allelic copies), unmerged isoforms, split genes and fragments. These appear as ~99–100%-identity "paralogs", and an identity-driven model will put them at the top of its ranking.
    - (a) Add a documented, generic pre-filter or flag for near-identical pairs, for example identity ≥ 98% over ≥ 90% of the shorter sequence, and fragment-containment pairs (one sequence nearly contained in the other). Output them with a flag rather than silently dropping them.
    - (b) Add a robustness experiment on development data. Degrade the human and yeast proteomes: fragment 10% of proteins, inject 5% isoform or allelic copies, and swap in longest-isoform-per-gene from Ensembl `pep.all`. Report the change in AUPRC and in the top-k composition.
    - *Why:* this is the failure most likely to embarrass the tool on its stated application, and the current evaluation, which uses curated UniProt reference proteomes, cannot see it.

14. **MUST: make panel features robust to how densely the panel samples the query's relatives, and test that robustness.**
    - QfO 2026_02 is vertebrate-heavy. It has 15+ vertebrates including macaque, and gorilla and chimp will probably be removed by the 98% self-exclusion rule. It has 4 arthropods and 5 plants.
    - It has no pre-WGD Saccharomycetaceae. Yeast's nearest panel relatives are Candida, Yarrowia and Aspergillus, so the yeast WGD-vs-SSD age proxy will be weak.
    - Count-type features ("number of species with a hit", "fraction of species where best hit A = best hit B") therefore mean different things for human, yeast, Arabidopsis and a new insect.
    - (a) Prefer divergence-relative definitions over species counts. Examples: is id(A,B) greater than id(A,O) and id(B,O) for each outgroup O; the age bin of the duplication, scaled by outgroup distance; asymmetry relative to the nearest shared outgroup.
    - (b) Add a **panel-thinning simulation** on development data. Recompute human features with all panel vertebrates except one distant one removed, and yeast features with all non-Saccharomycotina fungi removed. Report AUPRC; a large drop means the features will not transfer to sparsely sampled lineages.
    - (c) Optionally balance the panel by a stated taxonomic rule, for example adding UniProt reference proteomes so that each major clade (Viridiplantae, Arthropoda, other Metazoa, Fungi, protists) has at least 5–8 species. Fix the panel before any test access. Cost is about 11 s per species.
    - *Why:* the deployment targets, plants and insects, are exactly the sparsely sampled lineages.

15. **SHOULD: validate the WGD and duplication-age proxies against held-out annotations, for evaluation only.**
    - Report the AUROC of the Tier B age proxy against YGOB yeast ohnologs and OHNOLOGS v2 human ohnologs. Never use those annotations as features.
    - *Why:* WGD status is one of the best-documented determinants of SL (Guan 2007; De Kegel 2021). It is the cheapest check that Tier B captures it, and if the check fails, the feature code is debugged before modelling.

### Scope and cost: cut or gate

16. **SHOULD: drop per-gene embedding PCA features (Tier C item 8).** *Why:* they encode lineage and family identity. That is the memorisation shortcut found by Bernett 2024 and TripleThreat. It will look fine in human↔yeast LOSO but shift for plants and insects, and it cannot be checked with the available unseen-species data.

17. **SHOULD: demote Tier D (ProstT5→3Di) to "remaining ideas" unless a cheap pilot justifies it.**
    - The pilot is yeast candidate-pair proteins only, at most 1,000 aa. It must show an identity-conditional gain whose development-fold CI excludes 0.
    - *Why:* the cost is about 19 h for human, against a 2 h default target. Dennler found Foldseek fident and LDDT highly correlated with identity (ρ 0.86 / 0.78), and human SL is the setting where identity was already "exceptionally" strong.

18. **SHOULD: use ESM-2 (MIT) as the default PLM; keep ESM C 300M as an ablation only.** *Why:* the ESM C license must be re-checked before release (compute survey), and a shipped tool needs clean redistribution terms.

19. **SHOULD: add explicit decision gates to M2/M3, and put all auto-research before the single test access.**
    - G1: after Tier A+B, if the model does not beat B2 in LOSO, investigate before adding tiers.
    - G2: a tier is adopted only if the development and LOSO gain, conditional on identity, has a CI excluding 0 **and** the runtime fits the budget.
    - Keep a ranked idea backlog in RESEARCH_TREE.md with a stopping rule. The final test runs only after the backlog is exhausted or time is nearly up.
    - *Why:* PROMPT.md asks for continued auto-research, but the design puts the test at M4 with nothing after it. Any idea that comes up after the test access cannot be evaluated honestly.

20. **COULD: drop the Drosophila enrichment check.** *Why:* the data are positive-only and curation-biased, so no defensible test is possible. State the missing insect test as a limitation instead.

21. **COULD: define what "adjust thresholds once to maximise recall" in paralog calling is actually maximising.** Recall alone drives toward the loosest threshold. Use, for example, maximum recall subject to a cap on pairs per gene, or recall at a fixed output size. Report output size for human: expect several hundred thousand pairs, with zinc-finger and olfactory-receptor families dominating.

22. **COULD: report published external predictions for comparison.** Where licenses allow, show Dennler & Ryan 2025 predictions on overlapping test pairs as a reference point (not as features), with the caveat that their training set may contain our test pairs.

## 3. Risks the design does not address

- **Input-quality shift:** haplotigs, isoforms and fragments in new assemblies (change 13). This is the highest-impact unaddressed risk.
- **Panel-density shift** between the training species (dense vertebrate neighbourhood for human) and the deployment species (sparse for plants and insects) (change 14).
- **Polyploid plants and lineage-specific WGDs.** Crops carry recent homeologs at very high identity and in large numbers. Training sees only the yeast WGD (~100 My) and the vertebrate 2R (~500 My). Scores for recent polyploid homeologs are extrapolations; state this as a limitation and look at the score distribution on rice and maize, without labels.
- **Tandem-duplicate blind spot.** SGA under-tests physically linked pairs, so yeast labels under-represent tandem SSDs, which are common in plants. This adds to the ascertainment bias, and screened human pairs share it.
- **Gene-prior baseline B3 is trivial under the cluster split** (test genes are never seen). It is only informative on the random split, so report it there and say so.
- **Committing labels derived from non-open sources** (Ito, CHyMErA, Parrish, Klingbeil). The fallback "regenerate by script" is the right default; if it is used, make sure the committed label file still satisfies "deliver the curated labels used for training", for example by committing the CC BY subset plus a deterministic rebuild script for the rest.
- **The release model** is refit on all data after test, so its performance is not measured directly. Say so explicitly in REPORT.md.

## 4. Things that are fine as they are

- The sequence-only constraint as stated (D1), including the ban on auxiliary targets.
- Labels built from primary tables with tested negatives only, ambiguous pairs excluded, and a lenient tier kept as a sensitivity analysis.
- Cross-species homology-cluster splits, LOSO, a random split reported only to measure leakage inflation, and paired cluster bootstrap for Δ.
- Baselines B1 (best identity variant chosen on development data) and B2 (identity + family size + closest pair). B2 is the honest bar for showing that the extra tiers matter.
- End-to-end evaluation with caller misses scored at the minimum.
- Test-access logging through a single script; toy proteome with planted duplicates; shuffled-label control.
- Identity-stratified evaluation of the PLM and 3Di features.
- The generic self-exclusion rule for the panel.
- MMseqs2 for both searches, LightGBM plus a logistic-regression comparator, species-balanced weights, no species indicator.
- CLI layout (`slparalog setup` / `predict`), checksummed weights and panel, `/usr/bin/time -l` benchmark with power state recorded.
