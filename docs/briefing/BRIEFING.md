# Predicting synthetic-lethal paralog pairs from protein sequence alone

*A research briefing on a project run by an AI principal-investigator agent and a team of Claude agents*

**Status of this document.** Everything described here was planned, run, reviewed and written up by AI agents. **No human scientist has yet verified the results, re-run the analyses or checked the labels.** Every number below is copied from a file committed in the project's git repository; the file is cited next to the number (paths are relative to the repository root, with the commit where it helps). This briefing was itself written by an agent from that record, without re-running anything. Commit hashes (for example a88cff1) refer to the project's development repository. Its history is not published, because it contains data files that cannot be redistributed; the public repository `ai-hackathon-ParalogSL` holds the same files as of project close, with no data, as a fresh history.

## Lay summary

Many genes come in pairs of related copies, called paralogs, that arose when a gene was duplicated during evolution. Often a cell can lose either copy and survive, because the other copy covers for it, but it dies if it loses both. This "synthetic lethality" is useful in cancer research: if a tumour has already lost one copy, a drug against the other could kill tumour cells and spare healthy ones. Knowing which paralog pairs behave like this is also basic knowledge about how genomes stay robust. Today such pairs are found by expensive experiments in a few species. This project asked whether a computer program can predict them for any plant, animal or fungus, using only the protein sequences of its genes. An AI "principal investigator" agent, directing 49 other AI agents, built the tool, trained it on known pairs from baker's yeast and human cells, and tested it once on pairs it had never seen. The tool ranks true synthetic-lethal pairs far better than the usual rule of thumb ("the more similar the two copies, the more they back each other up"). But much of the human improvement comes from learning which pairs scientists chose to test, not only from biology. The tool also works in a worm, but not in a plant. The whole tool runs on a laptop in under two hours for the human genome.

> **Key results at a glance**
>
> - **Headline.** On a held-out test set opened once, the tool beats sequence identity in both training species: human average precision (AP) 0.349 against 0.106 for identity, yeast 0.541 against 0.190; the pre-registered statistic is T = 0.297 [95% CI 0.131, 0.395] (`results/final/final_test.md`, a88cff1).
> - **Main caveat (ascertainment).** Among human pairs tested by exactly the same set of studies, the tool does *not* beat identity: ΔAP −0.020 [−0.098, 0.060] (`results/final/final_test.md`). Much of the human headline reflects which pairs were chosen for screening. An adversarial review predicted this before the test (`reviews/results/tierB_gain_review.md`, 457ed8c).
> - **Evidence for real biology.** A model trained on yeast alone lifts human AP from 0.106 to 0.259 (`results/final/final_test.md`).
> - **Unseen species.** Gain in the nematode *C. elegans* (ΔAP +0.178 [0.064, 0.349], 16 positive pairs); no gain in the plant *Arabidopsis* (ΔAP −0.005 [−0.111, 0.097]) (`results/final/final_test.md`). Plants are not supported.
> - **Statistics.** The test was declared underpowered before any score existed (power 0.588 at the design effect size; `results/mde.json`). The yeast confidence interval is wide.
> - **Runtime.** Full human proteome (20,770 proteins, 686,225 pairs): 5,990 s wall time, 2.52 GB peak memory on a 10-core laptop (`results/bench/e2e_human.json`).
> - **Process.** 49 delegated agents (`results/process/agents.tsv`), 3 design-review rounds, independent review of every key result, one logged test access (`results/final/access_log.tsv`). Several process slips are documented (Section 2.5).

### Glossary

| term | plain meaning |
|----|-------------|
| **paralog** | One of two or more genes in the same genome that descend from a single ancestral gene through duplication. |
| **ortholog** | The "same" gene in two different species, separated by speciation rather than duplication. |
| **whole-genome duplication (WGD)** | An event in which an organism's entire genome was copied. Its surviving gene pairs are called *ohnologs*. The ancestor of baker's yeast underwent one, and early vertebrates underwent the duplications whose surviving human pairs are called "2R" ohnologs. Duplications of single genes are *small-scale duplications* (SSD). |
| **synthetic lethality (SL)** | Losing either gene alone is tolerated, but losing both kills the cell or makes it grow very poorly ("synthetic sick or lethal"). |
| **genetic interaction, ε** | How much worse (or better) a double mutant grows than expected from the two single mutants. A strongly negative ε means the genes back each other up. |
| **combinatorial CRISPR screen** | An experiment that knocks out two genes at once in cultured human cells, for many pairs and often several cell lines. |
| **proteome** | The full set of protein sequences encoded by a genome; here, one protein per gene. |
| **sequence identity** | The fraction of aligned positions at which two proteins carry the same amino acid. The standard simple predictor of paralog backup. |
| **reference panel** | A fixed set of 80 other species' proteomes that each query is compared against. |
| **ascertainment** | The process by which pairs came to be measured at all. If scientists chose which pairs to test, the tested set is not a random sample, and a model can learn the choice instead of the biology. |
| **prevalence** | The fraction of pairs in a set that are positive (here, SL). |
| **AUPRC / average precision (AP)** | A score for a ranked list: it is high when the true positives sit near the top. A random ranking scores about the prevalence. Section 4 explains it. |
| **ΔAP** | AP of the model minus AP of a baseline on the same pairs. |
| **bootstrap confidence interval (CI)** | A range of plausible values obtained by recomputing the statistic on many resampled versions of the data (here 10,000), resampling whole groups of related genes. |
| **homology-grouped split** | Dividing data into training and test sets so that genes with detectable sequence similarity always land on the same side, which prevents a model from memorising gene families. |
| **development / test set** | Data used freely to build and choose the model (with *cross-validation*: training on part, evaluating on the rest, repeatedly), versus data sealed until one final evaluation. |
| **pre-registration** | Writing down, and committing, the analysis rules before seeing the data they will be applied to. |
| **gradient-boosted trees (GBDT)** | A machine-learning model made of many small decision trees, each correcting the errors of the previous ones. |
| **protein language model (PLM)** | A neural network trained on millions of protein sequences that turns each protein into a vector of numbers. |
| **agent / PI agent** | An instance of a large language model (here Claude) working with tools. The PI agent plans and decides; delegated agents carry out tasks and report back. |

---

## 1. The question and why paralog synthetic lethality matters

**What was asked.** The task was given as a one-page prompt (`PROMPT.md`). The tool must take the proteome of any eukaryote as one protein FASTA file, find every paralog pair, and give each an SL score. Every feature must be computed from the sequences alone, because the intended use is on newly sequenced, unannotated genomes such as plants and insects. Measured data such as expression, protein interactions or essentiality, and annotations such as Gene Ontology terms, were forbidden as inputs. Training and test labels had to come from known SL pairs in budding yeast (*Saccharomyces cerevisiae*) and human cells. The metric was AUPRC, and the tool had to "clearly beat sequence identity alone". It also had to run on one Apple M4 laptop with 10 cores and 24 GB memory, and its runtime and memory on a full human proteome had to be reported (`PROMPT.md`).

![Figure 1. How the tool works. A protein FASTA file (one protein per gene) is searched against itself with MMseqs2; the hits are the candidate paralog pairs. Tier A features come from this self-search (identity, alignment coverage, family context). Tier B features come from comparing each protein with other species' proteomes in a reference panel (conservation, single-copy orthologs, relative evolutionary rate, timing of the duplication). A gradient-boosted tree model combines them into an SL score for every candidate pair. Tier C (protein language model, dashed) was tested and not adopted. Schematic with no numbers; the frozen model uses 50 Tier A and 23 Tier B features and an 80-proteome panel (`REPORT.md` §4). Source: `docs/briefing/figures/SOURCES.md`.](figures/fig1_pipeline.png)

**Why it matters.** Paralogs are less often essential than single-copy genes, because a surviving copy can buffer the loss of the other. In DepMap, a large catalogue of gene-knockout effects in cancer cell lines, paralog buffering explained about 13–17% of variable essentiality (De Kegel & Ryan 2019, as summarised in `docs/lit/methods.md` §1.1). Paralog SL pairs are therefore candidate drug targets in tumours that have lost one copy, and they are a window on how genomes evolve redundancy. Experiments to find them are costly and have been done systematically only in yeast and in a few human cell lines; the literature survey found no labelled paralog-SL data set for insects or crops (`LITERATURE.md` §1, "Gaps").

**Why identity is the bar to beat.** The oldest and simplest rule is that more similar copies back each other up more often. Gu et al. 2003 found that compensation in yeast correlates strongly with sequence similarity (`docs/lit/methods.md` §2.1). In human data, identity alone reached AUROC about 0.73 and AP about 0.10 in De Kegel et al. 2021, whose full model with interaction and expression features reached AP about 0.23 (`LITERATURE.md` §1). Dennler & Ryan 2025, the closest prior work, reported that identity "performed exceptionally well" on human SL and that single-modality sequence or structure models only matched it (`LITERATURE.md` §1). Beating identity with sequence-only features was therefore known to be hard, but not impossible.

**What was open.** The survey listed three gaps: no labelled paralog-SL data for insects or crops, no study of yeast-to-human transfer of paralog-SL predictors, and no sequence-only tool that calls paralogs itself and is evaluated end to end (`LITERATURE.md` §1).

---

## 2. How the agentic project was organised

![Figure 2. The 49 delegated agents over time. One dot per agent, placed at its launch time (hours after the first agent, launched 27 October 2026, 20:44 UTC) and grouped by task category: literature (4), design review (13), data and labels (8), features and tool (6), evaluation and experiments (4), independent review and QC (10), release and write-up (4). Filled pink: Opus (26); hollow pink: Opus inherited from the lead session (4); green: Sonnet (19). The file records only launch times, not durations. Categories were assigned by a keyword map on the task description (listed in `docs/briefing/figures/SOURCES.md`). Source: `results/process/agents.tsv` (016001a).](figures/fig2_process.png)

### 2.1 The PI and its team

The prompt made an AI agent the principal investigator (PI). It told the PI to delegate, review, integrate and decide, and to choose between two models, Claude Opus 5.5 and Claude Sonnet 5.5, for its agents (`PROMPT.md`). The PI kept three documents as its memory: a research tree with dated decision notes D0–D22 to project close, plus D23 on the public export the next day (`RESEARCH_TREE.md`), a list of process lessons (`LESSONS.md`), and the design (`DESIGN.md`). The prompt asked the PI to keep its own context small by reading agents' summaries rather than raw logs, and to be able to resume from the tree after a context reset (`PROMPT.md`).

The roster of delegated agents was reconstructed from the PI transcript after the project ended (`results/process/agents.tsv`, 016001a). It lists 49 agents: 26 marked `opus`, 4 `opus (inherited)` (running on the PI's own model) and 19 `sonnet` (counted from the rows of `results/process/agents.tsv`; Figure 2). Opus did the literature surveys on methods and labels, all design reviews, the human and yeast label builders, the splits, the caller and Tier A and Tier B features, the evaluation harness, the label and harness reviews, the adversarial Tier B review, the review of the test script, and the report draft. Sonnet did the compute benchmark, downloads, the sealed unseen-species labels, packaging, PLM embeddings, the Tier B and Tier C code reviews, the thinning and robustness experiments, the test script itself, the re-review of its fixes that returned READY, the release refit and benchmark, and both retrospectives (`results/process/agents.tsv`).

The research tree records that by 10:00 on 28 October the weekly plan usage stood at 47%, after which "routine agents" were moved to Sonnet and reviews kept short (`RESEARCH_TREE.md` D14). Model choice was therefore partly a budget decision.

**Timeline.** The first commit is dated 2026-10-27 20:44 UTC and the project-close commit 2026-10-28 18:12 UTC (`git log`, da27709 and e901438). The single test access ran from 11:59:51 to 12:00:18 UTC on 28 October (`results/final/access_log.tsv`). The development repository's git history holds 102 commits up to the roster commit (`git log`).

**The human role.** The human gave the initial prompt and afterwards sent only short status questions (progress, time and usage), one "continue", and the request for this write-up. They reviewed or steered no scientific decision. The only action recorded in the research tree is that they switched off the laptop's Low Power Mode at 09:42 on 28 October; agents could not change it themselves (`RESEARCH_TREE.md`, M1 status summary and D7). Some data sources needed a manual browser download and were never obtained (`REPORT.md` §9, "missing source").

### 2.2 Literature first, then design by critique

Four agents surveyed prior methods, human label sources, yeast and other-species labels, and the compute feasibility of candidate tools on this laptop (`results/process/agents.tsv` rows 1–4; `docs/lit/*.md`). The PI condensed them into `LITERATURE.md` and drafted `DESIGN.md`.

The prompt required the design to be "reviewed and revised in at most three rounds" by agents representing an external PI and experts in prior work, statistics, data, algorithms and compute, AI-agent failure modes, and other needed fields (`PROMPT.md`). Round 1 had seven reviewers, round 2 five, and round 3 one consistency reviewer (`RESEARCH_TREE.md` tree items 2.1–2.3). Reviewers were asked to measure claims on the real files, not only to read the text, and this paid off:

- An algorithms reviewer found that single-linkage clustering at the planned threshold put **10,650 human genes (52% of the proteome)** into one component, which would have made homology-grouped splitting meaningless (`reviews/round1/algorithms_compute.md`).
- The same review found that **10,132 of 20,617 human queries (49%)** hit the 100-hit cap in the panel search, which censors exactly the multi-copy families of interest (`reviews/round1/algorithms_compute.md`).
- Round 3 found **14 contradictions** in DESIGN v3, which would have produced inconsistent labels (`LESSONS.md` rule 1; `RESEARCH_TREE.md` D4).

The round-1 reviewers converged on label inconsistency, degenerate clusters, panel-density artefacts and unenforced test discipline; their fixes became the backbone of the design (`RESEARCH_TREE.md` D2).

### 2.3 Pre-registration and the single test access

The design fixed, before any score existed: the label rules, the baselines, the primary statistic, the pass criterion, the bootstrap, the split rule and the test fraction (`DESIGN.md` §2.3, §5, §6). Label rules were committed before any label code ran (`REPORT.md` §3). Every development experiment was logged in `results/experiments.tsv` with its adoption rule declared beforehand (`docs/AGENT_RULES.md` rule 6), often in a separate "declare-only" commit such as 2ab485d, 6e3c00a or 048c11a (`LESSONS.md` rule 7).

Test labels and the unseen-species labels sat in a sealed directory. Only one script, `scripts/eval/final_test.py`, could read them, and only the PI could run it (`docs/AGENT_RULES.md` rule 3). The script checks its own preconditions before opening anything: the model is tagged `final-model-v1`, the working tree is clean, the model hash and frozen specification match, the sealed manifest matches, and no earlier access exists. It then writes an `access_started` record to disk before reading sealed labels (`REPORT.md` §4 "Guards"; `reviews/results/final_test_review.md` Q4).

### 2.4 Reviews and guards

The prompt required that "an agent that didn't produce a result reviews it before you rely on it" (`PROMPT.md`). Every key result has a review file in `reviews/results/`. The most consequential catches:

| review | verdict | what it caught | source |
|---|---|---|---|
| human labels | ACCEPT WITH FIXES | One cell line not merged across studies (a non-breaking-space bug in a cell-line identifier); a fallback rule behind 22% of positives that needed a flag and a sensitivity analysis; strong ascertainment | `reviews/results/human_labels_review.md` |
| yeast labels | ACCEPT WITH FIXES | Independent re-implementation agreed on all 9.1 M pairs; flagged 7 tie-broken positives and ribosomal-protein pairs (23% of positives) for sensitivity analyses | `reviews/results/yeast_labels_review.md` |
| evaluation harness | ACCEPT WITH FIXES | The human gain of Tier A over identity "disappears once ascertainment strata are respected"; the pre-registered pooled statistic was about half the per-fold one | `reviews/results/eval_harness_review.md` |
| Tier B features | ACCEPT WITH FIXES | Two asymmetry features were missing for 94–98% of pairs; the yeast WGD could not be dated by distance; 25.7% of worm pairs had no panel evidence | `reviews/results/tierB_review.md` |
| Tier C features | ACCEPT WITH FIXES | Features were produced from uncommitted code; an identifier-mapping step had no test | `reviews/results/tierC_review.md` |
| Tier B gain (adversarial) | PARTLY ARTEFACT | No leakage, but most of the human gain ranks pairs between screening-ascertainment cells | `reviews/results/tierB_gain_review.md` |
| test script | READY WITH FIXES, then READY | R1: nothing checked that the frozen model was the selected one; R2: eight pre-registered secondary analyses were missing; R3: a failing secondary would have discarded the primary result while still using up the single access | `reviews/results/final_test_review.md` |

The last row deserves emphasis. All three test-script defects could only have been fixed before the access. They were fixed in 563b15d at 11:55 UTC, the re-review returned READY at 11:57 (29e8fc8), and the access ran at 11:59 (`git log`; `results/final/access_log.tsv`). The access completed with no secondary errors (`results/final/final_test.json`, field `secondary_errors`).

Guards were code, not only instructions: `load_labels()` refuses sealed paths, tests check the hashes of protected files, a registry test checks that no feature comes from a forbidden source, and feature builders never see labels (`REPORT.md` §4 "Guards"). The M3 retrospective found no test-discipline incident after D10 (`reviews/retro/M3_retro.md` §2 item 7). The guards also produced false alarms: the sealed-path guard initially failed on legitimate known-answer tests and had to be extended (`RESEARCH_TREE.md` D5).

### 2.5 What went wrong

The research tree logs its own incidents (`REPORT.md` §9; `RESEARCH_TREE.md`):

- **Test-discipline incident (D10).** A label reviewer computed one identity-only AP on unsplit labels, including future test pairs, before splits were frozen. The value was not recorded and informed no decision; briefs were changed to forbid scoring labelled pairs.
- **Specification bug (D9).** The split size cap counted the wrong unit, so one 238-gene human family (1,611 labelled pairs, 18% of human labelled pairs) would have landed in the test set. It was fixed before any score existed.
- **Runs from uncommitted code (D13, D18).** In D18 the thinning run imported another agent's uncommitted edit; it was rescued only because the reference condition reproduced the stored value 0.1432 exactly. Tier C was re-run from a clean tree (00991c6).
- **Shared working tree.** Agents editing one working tree caused an index race (415356b) and a review of still-changing code (`reviews/retro/M3_retro.md` §1 item 2).
- **Misremembered specification (D17, D19).** A decision note misstated the primary test population; the test-script review caught it by quoting the design (`LESSONS.md` rule 2).
- **Scaling before smoke-testing (D11).** A full human run took 2.3 h and flagged 67% of queries under a rule that then had to be changed.
- **Skipped test suite (D22)** for two release scripts, and **timings under load** (load averages 95–154 in D13).

Two retrospectives by fresh agents asked "what should we stop doing?" and "what worked well?" (`reviews/retro/M1_retro.md`, `M3_retro.md`); the distilled rules are in `LESSONS.md`.

### 2.6 A rule refined after seeing development results

One decision needs particular transparency. After the first Tier A development report showed that the human "gains" were learned ascertainment, the PI added a human stratum guard to the adoption rule, changed the selection statistic to a per-fold paired one, and added a deconfounding option (`RESEARCH_TREE.md` D12); the condition was later restated (95ac397). Only development data had been seen, and both changes made adoption stricter. The report states that every adoption decision is the same under both wordings (`REPORT.md` §9). It is still data-informed rule-making, and the retrospective lists it as a pattern to avoid: "Freeze the adoption rule and the candidate list before the first batch of a new family" (`reviews/retro/M3_retro.md` §1 item 7). By then 144 configurations had been tried (`reviews/retro/M3_retro.md`).

---

## 3. Building the benchmark

### 3.1 One definition of SL for two very different experiments

Yeast and human data come from different technologies. Yeast labels come from systematic crosses of deletion mutants (Costanzo et al. 2016, the "SGA" global genetic-interaction map); human labels come from combinatorial CRISPR screens in cancer cell lines. Round-1 reviewers asked for one biological definition (`RESEARCH_TREE.md` D2): **SL = a significant negative genetic interaction plus a strong absolute double-mutant defect, in a context where both single mutants are testable** (`REPORT.md` §3). The second half matters: a negative interaction between two genes whose single mutants are already sick says little about backup.

- **Yeast primary positive:** ε < −0.08, p < 0.05 and double-mutant fitness ≤ 0.5; negatives need double-mutant fitness > 0.5 with a near-zero interaction (`REPORT.md` §3). Dean et al. 2008 data enter under their own rule, and BioGRID low-throughput SL records form a separate tier only.
- **Human primary positive:** SL in at least 2 cell lines, plus either at least 2 studies with penetrance (fraction of testable lines in which the pair is SL) ≥ 0.10, or penetrance ≥ 0.25. Negative: testable in at least 2 cell lines and SL in none (`REPORT.md` §3). Calls come from nine screens: Dede, Thompson, Parrish, Ito, CHyMErA, in4mer/Chou, Klingbeil, Harle and Chan (`REPORT.md` §3).
- **Only tested pairs are labelled.** An untested pair is never a negative. Labels were built from primary tables, never from SL databases (SLKB, SynLethDB), because those mix scoring methods and mostly lack tested negatives (`LITERATURE.md` §3).

**Label universe.** A labelled pair must be both tested and a paralog pair in Ensembl Compara release 116. Compara only decides which tested pairs count as paralogs; it is never a feature, and it is independent of the tool's own paralog caller, so caller recall can be measured honestly (`REPORT.md` §3; `DESIGN.md` §2.3).

**Checks.** The human builder reproduced 12 of 12 published counts; an independent re-implementation of the yeast rules agreed on all 9.1 M measured pairs (`REPORT.md` §3; `RESEARCH_TREE.md` D8, D10).

![Figure 3. Sizes of the curated label sets (log scale; dark = SL positives, grey = negatives). (a) All labelled pairs before splitting, by tier: human counts are over all 11,125 tested pairs (primary 222 / 8,040; lenient 907 / 9,401; two-study 152 / 2,196), yeast counts within the Compara paralog universe (primary 104 / 2,049; lenient 183 / 2,047; strict 68 / 2,049; BioGRID 227 / 1,998). (b) Development split (human 123 / 4,297; yeast 62 / 1,143), the called test pairs actually scored (human 82 / 2,723; yeast 40 / 703), and the sealed unseen species (*Arabidopsis* strict 51 / 107 and broad 190 / 107; worm 16 / 120). Sources: `results/labels/human_label_qc.json`, `results/labels/yeast_label_qc.json`, `results/labels/unseen_label_counts.json`, `results/splits/split_summary.json`, `results/final/final_test.json`; see `docs/briefing/figures/SOURCES.md`.](figures/fig3_labels.png)

| set | positives | negatives | source |
|---|---|---|---|
| human, primary, in universe | 205 | 7,026 | `results/splits/split_summary.json` |
| yeast, primary, in universe | 104 | 2,049 | `results/splits/split_summary.json` |
| development split, human / yeast | 123 / 62 | 4,297 / 1,143 | `RESEARCH_TREE.md` D9 |
| test split, human / yeast | 82 / 42 | 2,729 / 906 | `RESEARCH_TREE.md` D9 |
| *Arabidopsis* (Cusack 2021), strict / broad | 51 / 190 | 107 / 107 | `results/labels/unseen_label_counts.json` (via `REPORT.md` §3) |
| *C. elegans* (Tischler 2006) | 16 | 120 | same |

*S. pombe* was planned but its interaction matrix needed a manual browser download and was never obtained; no insect data set with tested negatives exists (`REPORT.md` §3).

### 3.2 The ascertainment problem appears early

The human label reviewer found that the number of cell lines in which a pair was testable, a fact about the experiment and not the biology, by itself gives AUPRC 0.127 against a prevalence of 0.027 (`RESEARCH_TREE.md` D8). This is ascertainment: screens chose which pairs to test, and they chose pairs they expected to be SL. The literature had quantified this: screened human pairs are 73% WGD pairs against 19% of all paralog pairs, 89% "closest pairs" against 13%, with median family size 3 against 11 (De Kegel et al. 2021, as summarised in `LITERATURE.md` §2 item 4). Any feature that predicts "this pair would be chosen for a screen" can therefore look predictive of SL.

Because of this, the PI pre-declared, before any model result: an ascertainment baseline B4 (number of testable lines, a diagnostic only, never a tool input); stratified reporting by the number of testing studies (1, 2, ≥3); and a sensitivity analysis without labels that depend on the Harle/Chan fallback (`RESEARCH_TREE.md` D8).

### 3.3 Why homology-grouped splits

If a test pair shares a gene, or a close relative of a gene, with a training pair, a model can score well by recognising the family rather than learning general rules. In the protein-interaction field, random pair splits let deep models learn similarity shortcuts, and performance became near-random once leakage was removed (Bernett et al. 2024, `docs/lit/methods.md` §3; Park & Marcotte 2012, §4.1). For paralogs the risk is acute: pairs A–B and A–C from the same family share a gene, and yeast and human families are often homologous to each other (`docs/lit/methods.md` §4.1).

The split was therefore built from a homology graph over all labelled genes of *both* species, with each labelled pair first merged into one node, and whole connected groups assigned to development or test (`REPORT.md` §4 "Splits"). Groups holding more than 10% of a species' labelled pairs were forced into development; this moved one 238-gene human family with 1,611 pairs and no primary positives (`REPORT.md` §4). The test fraction was 0.40, drawn once (`results/mde.json`).

**Power, declared in advance.** A simulation run before any feature existed estimated the probability that the full decision rule would detect a modest true gain (T = 0.04 AP). At the largest allowed test fraction it was 0.588, below the 0.8 target, so the primary comparison was **pre-declared underpowered** (`results/mde.json`, fields `power_at_chosen` and `underpowered`). A non-significant result was to be reported as "not demonstrated", not "no effect" (`DESIGN.md` §6.1).

### 3.4 Unseen species

Two species with tested negatives were sealed and never used for training, selection or tuning: *Arabidopsis thaliana* double mutants (Cusack et al. 2021) and *C. elegans* combinatorial RNAi (Tischler et al. 2006) (`DESIGN.md` §2.3). Close relatives of these species already in the panel were not supplemented, so the tests resemble the "no close relative" use case (`DESIGN.md` §2.2). They were pre-declared as descriptive only, with no pass/fail (`REPORT.md` §6). The *Arabidopsis* labels measure genetic redundancy in a plant, a related but different endpoint, at high prevalence (0.323 for the strict tier; `results/final/final_test.md`).

---

## 4. Baselines and the evaluation metric

### 4.1 What AUPRC measures, intuitively

Imagine the tool ranks all candidate pairs from most to least likely SL, and you walk down the list. Each time you reach a true SL pair, you note the *precision* so far: the fraction of pairs above and including it that are truly SL. **Average precision (AP)**, the project's AUPRC (`sklearn.metrics.average_precision_score`; `docs/AGENT_RULES.md` rule 5), is the mean of those precisions. It is 1 if every SL pair sits above every non-SL pair.

**Why prevalence matters.** A random ranking scatters the positives, so its precision hovers around the overall fraction of positives, the prevalence. Random AP is therefore about the prevalence, not 0.5. In the human test set the prevalence is 0.029 and in yeast 0.054 (`results/final/final_test.md`); a human AP of 0.349 is therefore far above chance, whereas a similar AP in a set where a third of pairs are positive would be barely above chance. AP cannot be compared across sets with different prevalence (Saito & Rehmsmeier 2015 and McDermott et al. 2024, `docs/lit/methods.md` §4.4). For that reason every AP in the project is reported with its prevalence and counts, and methods are compared only on identical pairs (`docs/AGENT_RULES.md` rule 5; `DESIGN.md` §6.2). The *Arabidopsis* strict set, with prevalence 0.323, shows the point: an AP of 0.593 there is a much weaker result than it looks (`results/final/final_test.md`).

There is a second, practical consequence. Across the whole human genome, fewer than 0.5% of paralog pairs are SL (Flister et al. 2025, `LITERATURE.md` §1). Precision on a screened set, where prevalence is several times higher, overstates what a user would see when ranking a whole proteome. The design asked for expected precision at genome-wide prevalence, but it was not computed in the single access (Section 6.4).

### 4.2 The baselines

| baseline | definition | role |
|---|---|---|
| B0 | random | sanity check |
| **B1** | **sequence identity**: the best of seven pre-registered identity variants on development folds, which was global (Needleman–Wunsch) identity | primary comparator |
| B2 | logistic regression on identity + family size + reciprocal-best-hit flag | stronger "trivial" comparator |
| B3 | per-gene prior | leakage check; uninformative by construction under homology-grouped splits |
| B4 | number of testable cell lines | human ascertainment diagnostic, never a tool input |

(`REPORT.md` §4 "Baselines"; selected B1 in `results/dev/dev_primary_AB_ABC_noasym_6bb9438.md`.)

Choosing the best identity variant on development data, and adding B2, followed the methods survey's advice (`docs/lit/methods.md` §4.5).

### 4.3 The primary statistic

The primary statistic **T** is the mean, over human and yeast, of ΔAP (model minus B1) on the test set. The tool passes if the 95% CI of T excludes 0 and ΔAP > 0 in each species; "clearly beats" additionally requires a mean relative gain of at least 0.20 (`REPORT.md` §4 "Statistics"). CIs come from 10,000 bootstrap resamples of homology groups, stratified by species composition and by whether the group holds a positive (`REPORT.md` §4). During development, selection used **T_pf**, the mean of per-fold paired ΔAPs over five-fold grouped cross-validation repeated with five fixed seeds (`RESEARCH_TREE.md` D12a).

---

## 5. Model development, step by step

**Figure 4** summarises the development record.

![Figure 4. Development results for the selected configuration of each model family. (a) T_pf, the mean per-fold paired gain in AP over identity (B1), with 95% group-bootstrap CI. (b) Human AP averaged over the three ascertainment strata (1, 2, ≥3 testing studies); the dashed line is identity's value, the adoption guard. Open circles fail the guard: every Tier A variant, including deconfounded training and rank-averaging with identity, has a positive T_pf but falls below the guard. Filled circles pass. The star marks the adopted GBDT on Tier A + B. B2 is shown for reference. Sources: `results/dev/dev_primary_AB_ABC_noasym_6bb9438.json` (e6849e8), `dev_primary_ablation_anchor_95ac397.json` (6574597), `dev_primary_dc_tierA_68ca19a.json` (157aa80); see `docs/briefing/figures/SOURCES.md`.](figures/fig4_dev.png)

All development runs used the same population: 4,415 called human pairs (123 positive, prevalence 0.0279) and 996 yeast pairs (60 positive, prevalence 0.0602) (`results/dev/dev_primary_AB_ABC_noasym_6bb9438.md`). By the M3 retrospective, `results/experiments.tsv` held 663 rows (`reviews/retro/M3_retro.md` §1 item 7).

### 5.1 The paralog caller

Because the tool must work on unannotated proteomes, it calls paralogs itself. It runs MMseqs2 all-against-all and calls a pair when a bitscore-based E-value, computed with a fixed reference database size so that calling does not depend on proteome size, is ≤ 1e-3 in either direction; there is no identity or coverage floor (`REPORT.md` §4). It calls 686,225 human and 18,014 yeast pairs, recovering 89% (human) and 82% (yeast) of Compara paralog pairs, and 98–100% above 30% identity (`REPORT.md` §4; `RESEARCH_TREE.md` D11, D13).

### 5.2 Tier A and the ascertainment trap

**Tier A** contains 50 model features from the self-search alone (`REPORT.md` §4): several identity measures, coverage, bitscore ratios, length ratio, and *family context*, meaning how many other relatives each gene has at several similarity levels, the partner's rank among them, whether the two are reciprocal best hits, how many third paralogs they share, and per-gene length and low-complexity fraction.

On development data, logistic regression (LR) on Tier A reached human AP 0.179 against 0.118 for identity, and T_pf 0.058 [0.026, 0.099] (`results/dev/dev_primary_AB_ABC_noasym_6bb9438.md`). That looked like a gain. The harness review showed otherwise: within testing-study strata, the gain vanished (`reviews/results/eval_harness_review.md`, finding M2). The human stratum-averaged AP of LR Tier A was 0.180, *below* identity's 0.201 (`results/dev/dev_primary_AB_ABC_noasym_6bb9438.md`). Tier A had learned which pairs get screened many times, and heavily screened pairs are more often positive.

The logical chain was: because the gain disappeared within strata, the PI added a guard requiring that a candidate not lower the human stratum-averaged AP below identity (`RESEARCH_TREE.md` D12b; Section 2.6). An ablation then located the signal: similarity features alone gave a negative T_pf, and the family-context group carried the apparent gain while still failing the guard (`REPORT.md` §5; `results/dev/dev_primary_ablation_anchor_95ac397.json`). Family context is exactly what distinguishes screened from unscreened pairs (Section 3.2).

Two repair attempts failed:

- **Deconfounded training.** The number of testing studies was given to the model during training only, then fixed to a constant at prediction, so that the model could attribute that variation to the nuisance variable (`RESEARCH_TREE.md` D12c). LR and GBDT reached stratum-averaged AP 0.179 and 0.168, still below 0.201 (`results/dev/dev_primary_dc_tierA_68ca19a.json`, via `REPORT.md` §5).
- **Identity-anchored models.** A GBDT forced to be monotone in identity and rank-ensembles of the model with B1 also failed the guard. The closest was the LR rank-ensemble, 0.200 against 0.201 (`REPORT.md` §5; `results/dev/dev_primary_ablation_anchor_95ac397.json`).

No Tier A configuration was adopted (`results/dev/dev_primary_AB_ABC_noasym_6bb9438.md`, "Selection").

### 5.3 Tier B: comparative features from a reference panel

**Tier B** compares each protein with 80 reference proteomes (about 1.55 M proteins), chosen by a written rule that filled gaps in the Quest-for-Orthologs core set with more fungi, plants, nematodes and insects (`REPORT.md` §4; `DESIGN.md` §2.2). It contributes 23 model features (`REPORT.md` §4). The key ideas, in plain words:

- **Divergence calibration without a species tree.** For each panel species, the tool estimates how far it is from the query (d_s) from the median identity of a core set of reciprocal best hits. Validated without labels, d_s correlates with published divergence times at Spearman 0.78 (human), 0.87 (yeast), 0.77 (worm) and 0.59 (*Arabidopsis*) (`results/features/tierB_validation.json`, 3d5209d).
- **Did the duplication happen before or after each species split off?** If, in a given panel species, both copies best match the *same* single gene, the duplication probably happened after that species diverged; if they match different genes, before. From this come bounds on the duplication's age.
- **`tb_n_post_single`:** the number of panel species with a single post-duplication ortholog, a "one gene there, two here" pattern.
- **Relative rate (`tb_rel_rate`):** how fast the pair has diverged relative to the age of the duplication.
- **Conservation:** in how many panel species each gene has a detectable relative.

These features recover the human vertebrate 2R ohnologs without any labels (AUROC 0.92). They do *not* recover the yeast WGD by depth: the pre-WGD species *Torulaspora* sits closer to *S. cerevisiae* than several post-WGD species, a long-branch effect that a distance axis cannot fix; only `tb_n_post_single` separates yeast ohnologs (AUROC 0.91) (`REPORT.md` §5; `reviews/results/tierB_review.md` §4).

On development data, adding Tier B passed the adoption rule for the first time. GBDT A+B reached human AP 0.361 and yeast AP 0.591, T_pf 0.143 [0.095, 0.187], and stratum-averaged human AP 0.260 against the 0.201 guard (`results/dev/dev_primary_AB_ABC_noasym_6bb9438.md`). LR A+B also passed, at T_pf 0.099 (same file). Two caveats were already visible: the per-fold yeast ΔAP of GBDT A+B was only 0.040, and a classifier could identify the species from the A+B features with AUROC 0.998 (`results/dev/dev_primary_AB_ABC_noasym_6bb9438.md`; `REPORT.md` §5 "Controls").

### 5.4 The adversarial review of the Tier B gain

Because the gain was large and Tier A had already fooled the harness once, the PI did not adopt A+B outright. It marked A+B as a provisional candidate and commissioned an adversarial review, "leakage, ascertainment, shuffles, seeds, family concentration", before adoption (`RESEARCH_TREE.md` D15). An Opus agent that had not built the features or the harness ran it (`results/process/agents.tsv` row 39; `reviews/results/tierB_gain_review.md`). Its verdict was **PARTLY ARTEFACT**:

**No leakage, and the signal is real.**
- Shuffling Tier B rows within species returned T_pf to 0.040, the Tier A level, so the gain needs the right values on the right pairs.
- Regrouping folds by any homology hit gave T_pf 0.150; dropping the five largest positive-holding families kept human ΔAP at 0.241–0.262.
- A+B beat A in 9 of 9 human source libraries, with the CI excluding 0 in 8.
- A model trained on yeast alone lifted human per-fold ΔAP over B1 from 0.083 (A) to 0.144 (A+B).

**But ascertainment carries much of it.**
- Within *exact study-set cells*, groups of pairs tested by exactly the same set of studies and therefore with identical ascertainment, A+B beat identity by only +0.010 [−0.057, 0.063].
- Given the ascertainment covariates plus identity, the A+B score added only +0.015.
- Tier B alone predicted *whether a pair was screened at all* with AUROC 0.84.
- The reviewer concluded: "Treat about 0.01–0.05 as the realistic transferable gain, not 0.14."

(All from `reviews/results/tierB_gain_review.md`, 457ed8c.)

**Which features drive it, and do they make biological sense?** These features rank among the most important in both species (Figure 7, Section 6.5), and their standardised positive-versus-negative differences agree in direction across species: more one-to-two single-ortholog patterns (`tb_n_post_single` +1.08 human, +1.60 yeast), slower relative rate (`tb_rel_rate` −0.70, −0.63), and higher conservation (`tb_cons_min_pct` +0.38, +0.28). The idea that older duplicates are more often SL was *not* supported: the depth feature reversed direction between species and between human strata (`reviews/results/tierB_gain_review.md` §5).

**Decision.** GBDT A+B (configuration `GBDT_AB_17`) passed the pre-declared rule and was adopted *by rule*; the review's caveat was carried into the report, and the within-cell analysis was made a pre-declared test secondary (`RESEARCH_TREE.md` D17a, D17c). This is a judgement call worth noting: the rule said adopt, the review said the gain was mostly artefact, and the PI resolved the tension by adopting while pre-registering the test that could expose the artefact.

### 5.5 Tier C not adopted

**Tier C** added similarity of ESM-2 (150M-parameter) protein-language-model embeddings. Over Tier A it added ΔT_pf 0.006 [−0.000, 0.014] for LR and 0.013 [−0.010, 0.035] for GBDT; over A+B, about +0.014 in human and nothing in yeast; and it cost 5,982 s of embedding time for the human proteome in Low Power Mode (`REPORT.md` §5; `results/bench/plm_esm2_150m_human.json`). GBDT A+B+C actually scored slightly higher on development (T_pf 0.154) but was not adopted, because its increment's CI included 0 and the runtime cost was large (`RESEARCH_TREE.md` D17b). It remains an optional install. The retrospective noted that this tier was fully built and reviewed before a cheap increment check, which wasted effort (`reviews/retro/M3_retro.md` §1 item 5). The literature had warned that PLM similarity is often identity in disguise (TripleThreat, ICLR 2026, `docs/lit/methods.md` §3).

### 5.6 Panel thinning: what if no close relative is in the panel?

A new plant or insect may have no close relative in the panel. The panel-thinning test removed panel species closer than a ceiling c at prediction time only. At c = 0.35, which cut the human panel from 77 to 65 usable species, T_pf fell from 0.143 to 0.127, a drop of 0.016 [0.003, 0.030] (`results/dev/dev_thinning_GBDT_AB_17.md`, eea23bc). Training on equally thinned panels restored 0.145, so the loss comes from a mismatch between training and prediction panels, not from missing information (same file). Training with random masking recovered only 34% of the drop, short of the pre-declared 50%, so masking was **not adopted** (same file; `RESEARCH_TREE.md` D18). Figure 8a (Section 8) plots these conditions.

### 5.7 The frozen model

The frozen model `AB-gbdt-7ff069ac` is a scikit-learn histogram GBDT with learning rate 0.03, 100 iterations, 7 leaves, minimum leaf size 50 and L2 penalty 1.0, on 73 features (50 Tier A, 23 Tier B), seed 11, trained on all 5,411 primary development pairs at clean commit 563b15d (`REPORT.md` §4; `RESEARCH_TREE.md` D19). Each species gets equal total weight and equal weighted prevalence, so yeast is not swamped by human pairs (`REPORT.md` §4). Its score ranks pairs; it is not a calibrated probability.

### 5.8 Where the working ideas came from

The final method is a gradient-boosted tree on self-search similarity and family-context features plus reference-panel comparative features. Almost every ingredient was proposed or documented in prior work, as catalogued in `docs/lit/methods.md` and `LITERATURE.md`:

| ingredient in the final method | prior work credited in the project's literature files |
|---|---|
| Sequence identity as the main predictor of buffering | Gu et al. 2003; De Kegel & Ryan 2019; De Kegel et al. 2021 (`docs/lit/methods.md` §1.1–1.2, §2.1) |
| Family size, closest-paralog / reciprocal-best-hit status | De Kegel & Ryan 2019; De Kegel et al. 2021; Cusack et al. 2021 (reciprocal best match among top features) (`docs/lit/methods.md` §1.2, §2.2, §6.1) |
| Partner rank among a gene's relatives ("paralogs" rank) and species-rank features | Dennler & Ryan 2025, who found the "paralogs" rank nearly orthogonal to identity (`LITERATURE.md` §1; `docs/lit/methods.md` §1.5) |
| Conservation and gene age against other species | De Kegel et al. 2021 (conservation and age "very predictive"); Li, Yuan & Zhang 2010 (long-retained redundant pairs evolve slowly) (`docs/lit/methods.md` §1.2, §2.1) |
| Single-ortholog ("one-to-two") patterns as a WGD proxy; WGD pairs more often SL | Proposed as a proteome-only WGD proxy in `docs/lit/methods.md` §5; WGD effect from Guan et al. 2007 and De Kegel & Ryan 2019; low family size of ohnologs from Makino & McLysaght 2010 |
| Relative or asymmetric divergence rate lowering redundancy | VanderSluis et al. 2010; Byrne & Wolfe 2007; Kellis et al. 2004; Kuzmin et al. 2020 (`docs/lit/methods.md` §2.1) |
| Gradient-boosted trees, species-agnostic features, proteome percentiles | Recommended in `docs/lit/methods.md` §6.1; Dennler & Ryan 2025 used XGBoost, De Kegel et al. 2021 a random forest |
| Family-grouped evaluation, identity on the same pairs, AUPRC with prevalence | Park & Marcotte 2012; Hamp & Rost 2015; Bernett et al. 2024; Feng et al. 2024; Saito & Rehmsmeier 2015 (`docs/lit/methods.md` §4) |
| Awareness of screen ascertainment | De Kegel et al. 2021 (`LITERATURE.md` §2 item 4); ELISL 2024 explicitly targets selection bias (`docs/lit/methods.md` §1.7) |
| Search engine and PLM | MMseqs2 (doi:10.1038/nbt.3988) and ESM-2 (Lin et al. 2023) (`docs/lit/methods.md` §3) |

**What was new here.** As far as the literature survey found, the new parts are the combination and the evaluation, not the features: (i) a strictly sequence-only tool that calls paralogs itself and is scored end to end; (ii) comparative features computed against a fixed panel with a label-free divergence calibration, so the same code runs on any species; (iii) joint training on yeast and human with homology groups built across both species; (iv) an explicit test of yeast-to-human transfer of paralog-SL predictors, which the survey listed as a gap (`LITERATURE.md` §1, "Gaps"); and (v) a pre-registered, within-ascertainment-cell analysis of the human gain. The direction of the three leading Tier B effects matches the literature priors (`reviews/results/tierB_gain_review.md` §5); in that sense the project confirmed known biology more than it discovered new biology. Results are not directly comparable with De Kegel et al. 2021 or Dennler & Ryan 2025, which used different labels (for example DepMap-inferred SL), different splits and features that are forbidden here.

---

## 6. The frozen test result and every pre-declared secondary

The single access ran at 11:59:51–12:00:18 UTC on 28 October 2026 from tag `final-model-v1` (29e8fc8), and was logged with the git head, model hash, configuration hash and sealed-manifest hash (`results/final/access_log.tsv`). No model or analysis changed afterwards (`RESEARCH_TREE.md` D20).

![Figure 5. Primary test result. Average precision of the frozen model, identity (B1) and identity + family size + closest pair (B2) on the held-out test pairs, with prevalence (the chance level) indicated, for human (n = 2,805, 82 positive) and yeast (n = 743, 40 positive); and the gain over B1 with 95% group-bootstrap CIs. T = 0.297 [0.131, 0.395]. Source: `results/final/final_test.json` (a88cff1); see `docs/briefing/figures/SOURCES.md`.](figures/fig5_test.png)

### 6.1 Primary result

| species | n | positives | prevalence | AP model | AP B1 | AP B2 | ΔAP vs B1 [95% CI] | ΔAP vs B2 [95% CI] |
|---|---|---|---|---|---|---|---|---|
| human | 2,805 | 82 | 0.029 | 0.349 | 0.106 | 0.148 | 0.244 [0.156, 0.335] | 0.201 [0.122, 0.285] |
| yeast | 743 | 40 | 0.054 | 0.541 | 0.190 | 0.434 | 0.351 [0.050, 0.523] | 0.108 [0.010, 0.187] |

(`results/final/final_test.md`, a88cff1.)

T = **0.297 [0.131, 0.395]**, relative gain 2.077 against the 0.20 threshold; T against B2 is 0.154 [0.087, 0.210]. Verdict: **CLEARLY BEATS identity** (`results/final/final_test.md`). The observed effect is far larger than the design target of 0.04 for which the test was underpowered (`results/mde.json`). Labelled pairs the caller missed are excluded from the primary population, as pre-registered: 6 human pairs (no positives) and 205 yeast pairs (2 positives) (`REPORT.md` §6; `DESIGN.md` §6.2).

**A note on yeast.** The model's yeast test AP (0.541) is close to its development value (0.591), but identity's AP fell from 0.567 on development to 0.190 on test (`results/dev/dev_primary_AB_ABC_noasym_6bb9438.md`; `results/final/final_test.md`). On development, the yeast per-fold gain over identity was only 0.040 [−0.024, 0.097] (`results/dev/dev_thinning_GBDT_AB_17.md`, row N→N). Much of the large yeast test margin therefore reflects identity performing unusually poorly on this particular test draw. The cause is unexplained and was not investigated after the access (`REPORT.md` §6, point 3). The yeast CI, with a lower bound of 0.050, is wide.

### 6.2 Pre-declared secondaries

![Figure 6. Pre-declared secondary analyses and unseen species: gain in AP over identity (B1) on the test set, with 95% bootstrap CIs and the number of positives per row. Groups: how heavily each human pair was screened (studies, cell lines, within identical study sets, without Harle/Chan-dependent calls); cross-species transfer (single-species-trained versus final model); label tiers; similarity of test genes to development genes (human tertiles and pairs with no similar development gene); the unseen species; and other checks (identity tertiles, without ribosomal pairs, without orientation-conflict positives). Rows with fewer than 10 positives have no CI in the file and are omitted. The near-zero gain within identical study sets is the key ascertainment result. Source: `results/final/final_test.json` (a88cff1); omissions listed in `docs/briefing/figures/SOURCES.md`.](figures/fig6_secondaries.png)

**Ascertainment (human).**

| analysis | n | positives | prevalence | AP model | AP B1 | ΔAP vs B1 [95% CI] |
|---|---|---|---|---|---|---|
| 1 testing study | 1,885 | 10 | 0.005 | 0.010 | 0.005 | 0.005 [0.001, 0.022] |
| 2 testing studies | 659 | 23 | 0.035 | 0.333 | 0.280 | 0.054 [−0.146, 0.225] |
| ≥3 testing studies | 261 | 49 | 0.188 | 0.569 | 0.492 | 0.077 [−0.046, 0.192] |
| **within exact study-set cells (33 cells)** | 2,054 | 67 | 0.033 | 0.566 | 0.585 | **−0.020 [−0.098, 0.060]** |

(`results/final/final_test.md`.)

This is the most important secondary result. The adversarial development review predicted it (+0.010 [−0.057, 0.063] on development; `reviews/results/tierB_gain_review.md`), and the test confirmed it. Within each stratum the CIs for the 2 and ≥3 strata include 0. The stratum-averaged human AP is 0.304 for the model against 0.259 for identity and 0.102 for the ascertainment baseline B4 (`results/final/final_test.json`, `secondary.human_strata_D8`). The honest human claim is therefore: a large improvement in ranking pairs *as they are screened today*, but **no demonstrated gain among pairs screened under identical designs** (`REPORT.md` §6).

**Cross-species transfer.**

| analysis | AP model | AP B1 | ΔAP vs B1 [95% CI] |
|---|---|---|---|
| human, yeast-only-trained model | 0.259 | 0.106 | 0.153 [0.085, 0.244] |
| yeast, human-only-trained model | 0.497 | 0.190 | 0.307 [−0.031, 0.498] |

(`results/final/final_test.md`.)

A model that never saw a human label, and so never saw human screening ascertainment, lifts human AP from 0.106 to 0.259. This is the strongest single piece of evidence that part of the signal is transferable biology. It is not proof: yeast pairs may carry their own selection biases, and the features that predict screening in human (Section 5.4) could partly overlap with features that predict SL in yeast.

**Sensitivity analyses and tiers.**

| analysis | positives | AP model | AP B1 | ΔAP vs B1 [95% CI] |
|---|---|---|---|---|
| human, without Harle/Chan-dependent labels | 70 | 0.287 | 0.092 | 0.195 [0.109, 0.289] |
| yeast, without orientation-conflict positives | 39 | 0.529 | 0.181 | 0.348 [0.054, 0.517] |
| human, lenient tier | 337 | 0.223 | 0.134 | 0.088 [0.049, 0.136] |
| human, ≥2-studies tier | 61 | 0.406 | 0.210 | 0.196 [0.094, 0.286] |
| yeast, lenient tier | 58 | 0.507 | 0.195 | 0.312 [0.044, 0.461] |
| yeast, strict tier | 26 | 0.574 | 0.210 | 0.363 [0.030, 0.561] |
| yeast, BioGRID tier | 86 | 0.643 | 0.261 | 0.383 [0.026, 0.531] |

(`results/final/final_test.md`.) Without ribosomal-protein pairs, yeast AP is 0.460 against 0.106 for identity (31 positives; ΔAP 0.355 [0.057, 0.513]) (`results/final/final_test.json`, `secondary.yeast_and_ribosomal_D10`). The gain survives every label variant. Note that the human lenient tier, which includes pairs SL in only one cell line, shows a much smaller gain (0.088).

**End-to-end and leakage.** Scoring caller-missed labelled pairs at the minimum score gives yeast AP 0.518 instead of 0.541, and leaves human unchanged at 0.349 (`results/final/final_test.md`). Splitting test pairs into tertiles by their maximum similarity to any development gene gives ΔAP 0.218, 0.194 and 0.134, with a trend slope of −0.042 [−0.130, 0.066] (same file): no evidence that the gain depends on homology to training genes.

### 6.3 Unseen species

| species / tier | n | positives | prevalence | AP model | AP B1 | AP B2 | ΔAP vs B1 [95% CI] |
|---|---|---|---|---|---|---|---|
| *Arabidopsis*, strict (headline) | 158 | 51 | 0.323 | 0.593 | 0.598 | 0.692 | −0.005 [−0.111, 0.097] |
| *Arabidopsis*, broad | 297 | 190 | 0.640 | 0.778 | 0.757 | 0.810 | 0.021 [−0.023, 0.062] |
| *C. elegans* | 135 | 16 | 0.119 | 0.358 | 0.180 | 0.187 | 0.178 [0.064, 0.349] |

(`REPORT.md` §6, from `results/final/final_test.json`; ΔAP and AP values also in `results/final/final_test.md`.)

In *Arabidopsis*, the tool shows no gain over identity and scores below B2 on both tiers. Plausible contributors, none of them tested: the plant part of the panel is sparse, d_s agrees less well with divergence times in *Arabidopsis* (0.59), and the labels measure a different endpoint at high prevalence (`REPORT.md` §8 item 2). The non-gain was deliberately not used to retune anything (`RESEARCH_TREE.md` D20). The worm gain is clear but rests on 16 positive pairs. It is notable because 25.7% of worm pairs had no panel evidence at all (`reviews/results/tierB_review.md` §3).

### 6.4 Pre-registered items that were not done

The design asked for test AP stratified by WGD versus small-scale duplication, and for expected precision at genome-wide prevalence (`DESIGN.md` §6.2). Neither was computed in the single access, and neither can now be added without a second access (`REPORT.md` §8 item 11). Tier D (predicted structure) was not built. The access log lives at a different path from the one named in the design (`REPORT.md` §9). The release model, refit on all 8,959 called primary pairs from development plus test, has **no independent evaluation** (`REPORT.md` §8 item 7; `results/release/release_model.json`).

### 6.5 What the model relies on

Figure 7 shows held-out permutation importance on development data. Both the Tier A block (family context such as the gap to each gene's best other paralog and partner rank, plus identity) and the Tier B block (relative rate, conservation, single-ortholog counts) matter in both species: shuffling the Tier B block costs 0.239 AP in human and 0.271 in yeast, the Tier A block 0.286 and 0.463 (`reviews/results/tierB_gain_review.md` §5). Importance gives size, not direction; the directions are in Section 5.4.

![Figure 7. Which features the model relies on, human versus yeast. Held-out permutation importance (fall in development AP when one feature is shuffled; whiskers = SD over 25 seed × fold splits) for the ten single features with the highest mean importance across the two species, plus the whole Tier A and Tier B blocks. Leading features include Tier A family-context measures (gap to the best other paralog, partner rank) and the Tier B features `tb_rel_rate`, `tb_cons_min_pct` and `tb_n_post_single`. Importance shows size only, not direction; directions are given in the text. Source: `results/review/tierB_gain/importance.json` (457ed8c); see `docs/briefing/figures/SOURCES.md`.](figures/fig7_features.png)

---

## 7. Runtime and memory

The prompt required runtime and memory on a full human proteome (`PROMPT.md`). The end-to-end benchmark ran `slparalog predict` with 10 threads from an empty work directory and an empty panel-hit cache (`results/bench/e2e_human.json`, commit 211c264 recorded in the file):

| item | value |
|---|---|
| input | 20,770 proteins; 686,225 pairs scored |
| wall time | 5,990 s, of which the Tier B panel search took 5,277 s |
| CPU time | 47,240 s user + 535 s system |
| peak memory | 2.524 GB (process tree), 2.387 GB largest single process |
| peak scratch disk | 1.516 GB |
| machine | Apple M4, 10 cores (4 performance + 6 efficiency), Low Power Mode off |

(`results/bench/e2e_human.json`.)

The design target of at most 2 hours and 8 GB was met (`REPORT.md` §7). One-off panel setup took 16 s and 0.16 GB, excluding download time (`results/bench/e2e_setup.json`, via `REPORT.md` §7). Yeast takes 1,115 s and 1.46 GB (`results/bench/e2e_yeast.json`, via `REPORT.md` §7). The shipped pipeline reproduces the frozen model's scores on the evaluated feature tables exactly, with a maximum difference of 0.0 (`results/bench/e2e_consistency_*.json`, via `REPORT.md` §7).

The benchmark was not on an idle machine: the one-minute load average was 4.6 at the start and 48.6 at the end, mostly the run's own MMseqs2 threads but also system daemons and one unrelated user process (`results/bench/e2e_human.json`, fields `loadavg_start`, `loadavg_end`; `REPORT.md` §7). The panel search dominates, so the main lever for speed is restricting or parallelising it (`REPORT.md` §10 item 7).

---

## 8. Robustness to degraded proteomes

A newly sequenced genome rarely comes with clean gene models. The design pre-registered a robustness test (`DESIGN.md` §6.5); it was declared in b58dbc2 and run after the test access on development labels only, so it is exploratory and touches no test result (`RESEARCH_TREE.md` D21). Each proteome was degraded with one defect (20% of proteins truncated to 50–80% of their length; 5% of genes split in two; 2% fused with another gene; 15% of genes removed; near-identical "haplotig" copies of 5% of genes) or with all combined, and the shipped `predict` recomputed all features (`REPORT.md` §6a). Pairs that lost a gene or were no longer called were scored at the minimum for both the model and identity.

| species | defect | model AP, natural → degraded | identity AP, natural → degraded |
|---|---|---|---|
| yeast | none (control) | 0.626 → 0.626 | 0.586 → 0.586 |
| yeast | fragments | 0.626 → 0.621 | 0.586 → 0.538 |
| yeast | split genes | 0.626 → 0.609 | 0.586 → 0.546 |
| yeast | fusions | 0.626 → 0.628 | 0.586 → 0.591 |
| yeast | 15% missing | 0.626 → 0.442 | 0.586 → 0.422 |
| yeast | haplotigs | 0.626 → 0.604 | 0.586 → 0.586 |
| yeast | combined | 0.626 → 0.393 | 0.586 → 0.329 |
| human | combined | 0.387 → 0.228 | 0.140 → 0.111 |

(`results/dev/dev_robustness_GBDT_AB_17.md`, a84a4f1. AP here is the per-fold mean, so it differs from pooled development values.)

Figure 8 shows both robustness tests: panel thinning (a; Section 5.6) and degraded proteomes (b).

![Figure 8. Robustness. (a) Thinning the reference panel: T_pf (gain over identity, 95% CI) when panel species closer than c are removed at prediction only, at training and prediction, or with masking augmentation; the shaded band is the natural-panel CI. (b) Degrading the proteome: development AP (per-fold mean) of the model (blue) and identity (grey) on natural (open) and degraded (filled) proteomes, for each simulated defect in yeast and all defects together in yeast and human. Sources: `results/dev/dev_thinning_GBDT_AB_17.json` (eea23bc), `results/dev/dev_robustness_GBDT_AB_17.json` (a84a4f1); see `docs/briefing/figures/SOURCES.md`.](figures/fig8_robustness.png)

**Reading.**
- Fragments and split genes hurt identity more than the model; fusions do nothing.
- Haplotig copies cost the model 0.022 and identity nothing, plausibly because extra near-identical partners disturb family-context features.
- Missing genes dominate the loss for both methods; a pair with an absent gene cannot be scored.
- In human, with all defects combined, the margin over identity falls from 0.247 to 0.118 (drop 0.129 [0.070, 0.184]), roughly halving, although the model still ranks above identity (`results/dev/dev_robustness_GBDT_AB_17.md`).
- Over both species, T_pf falls from 0.143 to 0.090, a drop of 0.053 [0.008, 0.096] (same file).

The defect rates were chosen in the design, not estimated from real assemblies, so the effect on real draft proteomes is unknown (`REPORT.md` §6a). Human was run only for the combined case, for compute reasons.

---

## 9. Take-away messages

### 9.1 For biology and bioinformatics

1. **Comparative features carry real, transferable signal about paralog backup.** A pair is more likely SL when many other species have a single gene corresponding to both copies, when the pair evolves slowly for its age, and when both genes are widely conserved. These directions hold in yeast and human (`reviews/results/tierB_gain_review.md` §5), and a yeast-trained model improves human ranking (`results/final/final_test.md`). This confirms earlier observations (Section 5.8) in a stricter setting; it is not a new mechanism.
2. **Human combinatorial-screen labels are heavily ascertained, and sequence features can learn the ascertainment.** Family-context features predict which pairs were screened (AUROC 0.84 for Tier B alone; `reviews/results/tierB_gain_review.md`), and within identical screening designs the tool's gain over identity disappears (`results/final/final_test.md`). Anyone training or evaluating paralog-SL predictors on these screens should report a within-design comparison next to the headline.
3. **"Older duplicates are more often SL" was not supported** as a simple rule: the depth feature's direction differed between species and strata (`reviews/results/tierB_gain_review.md` §5), consistent with the conflicting literature on whether redundancy decays with divergence (`LITERATURE.md` §1).
4. **Transfer to plants failed in the one test available.** Before using such a tool on crop genomes, the panel needs more plant species and the *Arabidopsis* failure needs a label-free diagnosis (`REPORT.md` §10 item 2).
5. **Protein language-model similarity added little** beyond what self-search and panel features already capture (`REPORT.md` §5).
6. **The yeast WGD cannot be dated by distance alone** because of unequal evolutionary rates among yeast lineages; gene-tree reconciliation is the proposed fix (`reviews/results/tierB_review.md` §4; `REPORT.md` §10 item 3).
7. **Scores are rankings, not probabilities,** and precision on screened sets overstates precision genome-wide (Section 4.1).

### 9.2 For agentic science

1. **An agent team can carry a project from a one-page prompt to a pre-registered, single-access test in under a day.** The first commit is dated 20:44 UTC on 27 October and the test access began at 11:59 UTC on 28 October (`git log`, da27709; `results/final/access_log.tsv`). The deliverables (installable tool, tests, curated labels, benchmark and report) exist in the repository.
2. **Structure beat intelligence.** The most valuable safeguards were procedural: sealed data behind code guards, declare-before-run commits, a script that refuses to run unless its preconditions hold, and a written adoption rule that bound the PI in both directions. The rules rejected every Tier A configuration and rejected masking augmentation (`results/dev/dev_primary_AB_ABC_noasym_6bb9438.md`; `RESEARCH_TREE.md` D18).
3. **Adversarial reviews by agents that did not build the thing were the single most useful practice.** The Tier B gain review predicted the key test caveat (+0.010 on development within cells, −0.020 on test), and the test-script review caught three defects that could not have been fixed after the access (`reviews/results/tierB_gain_review.md`; `reviews/results/final_test_review.md`; `LESSONS.md` rule 8). Asking reviewers to *measure* on real files, not just read, found the problems that mattered in the design rounds (`reviews/retro/M1_retro.md` §2).
4. **Written lessons did not prevent repeat mistakes; code did.** "Commit before running" was a lesson after M1 and was broken again in M2 (`reviews/retro/M3_retro.md` §1 item 1). Agents sharing one working tree repeatedly interfered with each other. The retrospective's recommendation is to enforce rules in scripts (refuse a dirty tree) and isolate agents in separate worktrees.
5. **Agents misremember specifications.** A decision note and a brief both misstated the primary test population; a reviewer caught it by quoting the design (`LESSONS.md` rule 2). Quoting the pre-registered text in every brief is cheap insurance.
6. **Data-informed rule changes are a real risk even with good intentions.** The D12 guard was right and made adoption stricter, but it was added after seeing development results, with 144 configurations already tried (`reviews/retro/M3_retro.md` §1 item 7).
7. **Speed creates pressure on deliberation.** The test-script fixes were committed at 11:55, the re-review returned READY at 11:57, and the access began at 11:59 (`git log`; `results/final/access_log.tsv`). Nothing went wrong, but a human PI would rarely move from fixes to an irreversible step within minutes.

### 9.3 The limits of judging an AI-run project

- **No human has verified the results.** The labels, code, statistics and conclusions have been checked only by other agents.
- **"Independent" reviewers were not fully independent.** They were separate agent instances that had not written the code, but they were Claude models given briefs written by the same PI, and they may share blind spots. After a budget note, the PI moved routine agents to Sonnet and kept reviews short (`RESEARCH_TREE.md` D14); the Tier B and Tier C code reviews and the READY re-review of the test script were Sonnet agents (`results/process/agents.tsv` rows 36, 37, 45).
- **The record was written by the system being judged.** The tree, lessons, incident list and report were written by the PI agent; incidents no agent noticed would not appear. This briefing is also agent-written from that record and inherits its blind spots.
- **Pre-registration constrains but does not certify.** It shows that the test was not tuned; it does not show that the labels are right, that the code does what it claims, or that the interpretation is sound.

---

## 10. Caveats, in one place

| caveat | consequence | source |
|---|---|---|
| Human gain is mostly between screening-ascertainment cells | Within identical designs, ΔAP −0.020 [−0.098, 0.060]; human AP is not precision on unscreened pairs | `results/final/final_test.md` |
| *Arabidopsis* shows no gain | Plants not supported; B2 beats the model there | `results/final/final_test.md` |
| Primary comparison pre-declared underpowered | Power 0.588 at T = 0.04 | `results/mde.json` |
| Yeast test CI wide; identity unusually weak on the yeast test draw | Yeast ΔAP 0.351 [0.050, 0.523]; dev yeast gain was 0.040 | `results/final/final_test.md`; `results/dev/dev_thinning_GBDT_AB_17.md` |
| Small unseen sets | Worm result rests on 16 positives | `results/final/final_test.md` |
| Release model not evaluated | Users of the release model rely on the frozen model's test as a proxy | `REPORT.md` §8 item 7 |
| Adoption rule refined after Tier A development results | Stricter, decisions unchanged, but data-informed | `RESEARCH_TREE.md` D12; `REPORT.md` §9 |
| Pre-registered items missing | WGD/SSD stratified test AP and expected precision at genome-wide prevalence not computed | `REPORT.md` §8 item 11 |
| Process slips | Runs from uncommitted code, a test-discipline incident, shared-tree races, skipped test suite | `REPORT.md` §9; `RESEARCH_TREE.md` D10, D13, D18, D22 |
| Simulated degradation, not real assemblies | Effect on real draft genomes unknown; human margin halves under combined defects | `results/dev/dev_robustness_GBDT_AB_17.md` |
| Licensing | 44 of 222 human primary positives rest only on non-open evidence; licence check needed before public release | `REPORT.md` §3 |
| **No human verification** | All results are agent-checked only | this briefing |

## 11. What would come next

The project's follow-ups (`REPORT.md` §10) are: ascertainment-robust human training and evaluation; more plant proteomes in the panel; tree-based duplication dating; calibrated probabilities; more unseen species; and a faster panel search. The first step should be an independent human re-check of the labels, the test run and the claims made here.

---

*Repository: `ai-hackathon-ParalogSL`. Main record: `REPORT.md`, `RESEARCH_TREE.md`, `LESSONS.md`, `reviews/`, `results/`. Figures and their data sources: `docs/briefing/figures/` and `docs/briefing/figures/SOURCES.md`.*
