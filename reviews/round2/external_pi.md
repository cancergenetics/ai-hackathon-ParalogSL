# Round 2 review: external PI

Read: PROMPT.md, DESIGN.md v2, reviews/round1/*.md (all seven) and PI_response.md, RESEARCH_TREE.md.

## 1. Overall verdict

v2 is a real improvement. The label inconsistencies are fixed: one SL definition for both species, essential-member pairs treated the same way, Kuzmin/Dean pairs labelled by ε, and a penetrance clause in the human rule. Deployment realism is now taken seriously, and test discipline is enforced in code. I would fund it.

What remains is mostly over-specification. v2 adds three mechanisms (Leiden partitioning, panel-masking augmentation and a nested-CV co-primary). Each sounds rigorous, but each has an unstated failure mode, and a simpler rule does the job. One round-1 must item (D6, the paralog label universe) is marked accepted but is not actually implemented.

### Round-1 must items: coverage check

| Item | v2 status |
|---|---|
| Labels: unified SL definition, essential-member pairs, Kuzmin/Dean by ε, human penetrance clause, BioGRID tier, Flister out, gene as unit (ext PI 1–3; data D1–D5, D7–D9; prior work M1–M2) | Addressed |
| **Label universe independent of the caller (data D6)** | **Not addressed.** "All tested pairs from the sources" is the whole SGA matrix in yeast, not a paralog universe. See M1. |
| Power / MDE / test size (statistics 1–4; ext PI 8) | Addressed, but the nested-CV fallback is unsound. See M2. |
| Cluster split (statistics 5–6; algorithms 7; ext PI 9) | Addressed via Leiden. Labelled pairs can be cut. See M3. |
| Panel truncation, gap filling, divergence calibration, thinning (comparative M1–M4; algorithms 8–9; ext PI 14) | Addressed. Masking adds a leakage path. See M4. |
| Robustness to input quality (ext PI 13; comparative M5) | Addressed (§6.5 and QC flags). |
| Prior-work baselines and stratified reporting (prior work M3, M5) | Mostly addressed. The age-proxy-vs-annotation check (M3b) is not explicit. See S6. |
| Algorithms musts (search caps, symmetric features, no raw E-values, self-scores, no TM-score/LDDT, libomp, pinning, percentile definition) | Addressed |
| Agent-failure musts (1–4, 6–9, 12–13, 17, 20, 23–25, 27–29) | Addressed. Item 14 was partially accepted with a reason; I agree with that. |
| Unseen species in the success criterion (ext PI 7) | Downgraded to secondary with CIs, which is acceptable given power. Pre-register the wording; see S3. |

## 2. Concrete changes

### Must

**M1. Pre-register the paralog label universe for both species, independent of the shipped caller (data D6 is still open).**
- §3 says "label universe = all tested pairs from the sources". In yeast, SGA tested about 23M gene pairs, nearly all of them non-paralogous. As written, the yeast labelled set and its prevalence are undefined.
- Fix: in `config/labels.yaml`, define a fixed homology rule used only for building labels, for example a BLASTP/MMseqs2 search at a pinned version with E ≤ 1e-5 and coverage ≥ 0.3 of the shorter protein. Ensembl Compara paralogs at a pinned release are an alternative; using them to define the universe is not a feature. Apply the same rule to the human screen pairs, which drops Thompson's mutual-exclusivity pairs and the non-paralog controls.
- The caller's recall is then measured against this universe.
- *Why:* prevalence, B1's AP and caller recall all depend on this choice. It must be fixed before labels are built, not decided implicitly by the parser.

**M2. Drop the nested-CV co-primary fallback (§6.1, last bullet).**
- Model and tier choices (Tier C/D adoption, the hyper-parameter grid, caller thresholds) are made by agents over weeks on the development data. They cannot be re-run inside outer folds. Outer-fold estimates that include development data would therefore be optimistically biased, and a second primary analysis adds a forking path.
- The round-1 simulations put yeast MDE near or above 0.10 AP, so the switch would very likely fire. This would make a biased estimate the de facto headline.
- Replace with: compute the MDE at M1. The only allowed lever is the test fraction (and the yeast test share), set at M1 before any feature exists. The held-out test stays the single confirmatory analysis. A ΔAP whose CI includes 0 is reported as "not demonstrated". Repeated grouped development CV is reported as a secondary result and labelled optimistic.

**M3. Labelled pairs must never span the split, and connected components should be the default instead of Leiden.**
- As written, Leiden partitions a graph whose edges include labelled pairs. It can cut a labelled pair, leaving a pair with one gene in test and one in development, and the design does not say what happens to such pairs.
- *Must:* contract every labelled pair (and every called pair between labelled genes) into a single node before any partitioning.
- *Should, as a simplification:* the default should be connected components over stricter edges (E ≤ 1e-10 and both coverages ≥ 0.5) plus labelled-pair edges. Round 1 measured a largest human component of 941 genes under that rule, which is small enough to need no community detection. Use Leiden only as a documented fallback, if the largest component still holds more than 20% of positives.
- In both cases, keep the residual-homology reporting already in §6.1. If Leiden is used, pin `leidenalg`/`igraph`, the resolution and the seed in `config/`.
- *Why:* components have no resolution or cap parameter to fork on. Their residual homology is the standard, interpretable quantity (Hamp & Rost).

**M4. Ban proteome-constant features. Masking augmentation creates a species-indicator leak.**
- D_s, the divergence bins and "empty bin" missingness are properties of the whole proteome. With 2 training species × 4 masks, any feature that is constant within a proteome takes 8 values. Per-bin fractions that are NaN for every pair when a bin is empty are an example. Such a feature is a condition ID, which violates the "no species indicator" rule. LightGBM will use it to fit species-specific prevalence.
- Fix: the feature registry test fails any column that has zero within-proteome variance, or a missingness pattern that is constant within a proteome, on the toy and training proteomes. Encode panel information only in pair-relative form. Examples: the fraction of available panel species whose RBH identity to the query is below id(A,B); the D-scaled depth of the duplication. Do not use bin-indexed slots.
- Extend the source-library control to a "species/mask prediction" control on pair features.

### Should

**S1. Simplify the panel-masking augmentation and make it conditional.**
- (a) Define masks by absolute divergence: remove every panel species with D_s above t, for t ∈ {none, 0.75, 0.55}, fixed now. Do not use "the k closest bins", which is rank-based: for yeast, "3 closest bins" removes nearly everything, while for human it removes only mammals.
- (b) Run the panel-thinning test first, with no augmentation. Adopt the augmentation only if thinning lowers development ΔAP by more than its bootstrap SE and augmentation recovers it. Masks reuse the same search output, so this is cheap, but it should be justified by data, not assumed.

**S2. State the model-selection metric.** §5 names the one-SE rule but not what it is applied to. Use the development out-of-fold species-mean ΔAP vs B1, the same statistic as the primary criterion (my round-1 #10). Mean AP would be dominated by yeast prevalence.

**S3. Put runtime into tier adoption, and shrink Tier C.**
- Adoption rule (round-1 #19, G2): a tier is adopted only if its identity-conditional development gain has a CI excluding 0 **and** it fits ≤2 h / ≤8 GB on human. ProtT5-XL (≈11 h) and 3Di (≈11 h) can at most be opt-in flags, never the default.
- Tier C: start with one PLM (ESM-2 150M, mean-pooled cosine; MIT licence). Escalate to more models or pooling variants only on signal. The current model × pooling × metric matrix is many forks for a small effect.
- Drop the per-gene embedding PCs (item 9). LOSO human↔yeast cannot detect the lineage memorisation they encode, so the gate does not protect anything.

**S4. Make auto-research precede the single test access, with a stop rule.**
- M3 should be an explicit loop: a ranked idea backlog in RESEARCH_TREE.md, gate G1 (if A+B does not beat B2 in LOSO, investigate before adding tiers), and a stopping rule (backlog exhausted, or N consecutive non-adopted ideas, or time).
- M4 (test) comes only after that.
- *Why:* PROMPT.md asks to keep going until ideas run out, and nothing done after the test can be evaluated honestly.

**S5. Pre-register the unseen-species wording.** LOSO (human→yeast, yeast→human) is the best-powered proxy for unseen species. Pre-register that a generalisation claim requires its paired CI to exclude 0 in at least one direction, with Δ > 0 in both. Arabidopsis, C. elegans and S. pombe remain descriptive. Report each test species' nearest retained panel D_s, because gap filling may put a close relative of Arabidopsis in the panel.

**S6. Validate the age and WGD proxy before modelling.** Report the AUROC of the Tier B duplication-depth feature against YGOB and OHNOLOGS ohnolog status, for evaluation only (prior work M3b; my round-1 #15). This is the cheapest debug check that Tier B works at all.

**S7. Yeast labels: correct one rationale and add one missing rule.**
- (a) PI_response says the DMF ≤ 0.5 rule "raises yeast positive counts". With healthy single mutants (SMF ≈ 1), DMF ≤ 0.5 implies ε ≤ about −0.5, which is *stricter* than the old ε < −0.35. Do not re-tune; the rule is justified biologically. Correct the claim, and let the M1 MDE step show the power consequence.
- (b) Data D11 is missing from v2: mark physically linked yeast pairs (within the SGA linkage window) untestable unless Dean or BioGRID-LT gives a tetrad call. *Why:* tandem duplicates are high-identity, so the linkage artefact correlates with the main feature.

**S8. Close the label source list.** Add a time box and a rule that any source not parsed and QC'd by the M1 freeze is dropped and documented (my round-1 #5). Record the decision on GABI-DUPLO (round-1 #4), which the PI response did not mention.

**S9. Add the label-noise ceiling.** Cross-study and SGA cross-orientation AP (round-1 #11) is still missing. It is cheap, and it is needed to interpret every AUPRC in the report.

**S10. Fix the bootstrap stratification wording.** "Stratified by species and label" is not possible for a cluster bootstrap, since clusters contain both labels. Stratify by species and by "cluster contains a positive" (statistics 9).

### Could

- **C1.** Robustness test: keep fragments, near-identical copies and missing genes. Split/fused genes can be deferred; they cost more to simulate correctly than they tell us.
- **C2.** Commit `labels_open.tsv` (CC BY/CC0/MIT sources only) and report the human headline result on it as a sensitivity analysis (data D14).
- **C3.** Define the caller-threshold objective, for example maximum recall of the M1 universe subject to a cap on pairs per gene (round-1 #21).
- **C4.** Add a per-tier accuracy vs runtime/RSS table, benchmarked with Low Power Mode off.

## 3. Risks the design does not address

- **Proteome-level confounding** (M4) is the main new risk introduced by v2.
- **Yeast power** may fall rather than rise under the DMF rule (S7a). The project should be ready to report "not demonstrated in yeast".
- **Leiden residual homology:** if the fallback is used, cut edges guarantee non-zero test↔development homology. The strict no-hit subset may then be tiny, so report its n.
- **Schedule:** 9 human sources, 4 yeast sources, 3 unseen-species sets, an 80–90-species panel and four feature tiers. Labels (S8) and Tier C/D (S3) are where time will go.
- Unchanged from round 1: recent polyploid homeologs are an extrapolation; the yeast labels have a blind spot for tandem duplicates; B3 is only informative on the random split.

## 4. Things that are fine as they are

- The unified SL definition, the strict and lenient tiers, the human penetrance rule with a borderline margin, and Flister kept out of the primary tier.
- Per-gene proteome expansion; header stripping; self-exclusion at median RBH ≥ 0.98, provided M4 and S1 are applied.
- Divergence calibration through D_s, which needs no species tree; per-species hit guarantees; gap filling by a written, label-blind rule.
- The primary criterion (species-mean ΔAP vs B1, CI excludes 0, Δ > 0 in each species, relative gain ≥ 0.20) and the pre-registered B1 variant set.
- The sealed test, tag, access log, protected hashes, feature registry, gene-name scanner, quarantine and report-number checker.
- The robustness test and the `near_identical` / `fragment_suspect` flags.
- The refit model is never scored on the test set.
