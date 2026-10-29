# Predicting synthetic-lethal paralog pairs from protein sequences alone

> This is an exercise in assessing the bioinformatics capabilities of agentic AI. This particular exercise was run by the Ryan lab at University College Dublin as part of a larger hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects. All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication and bioinformatics tool.
>
> Models and machine: the project was run in Claude Code. The lead (PI) agent ran on Claude Opus 5.5 (`claude-opus-5-5`); it planned the work, delegated it, reviewed the results and made every decision. During the project it delegated 49 tasks to sub-agents (`results/process/agents.tsv`). 30 of them used Claude Opus 5.5: the literature review, the design reviewers, the label builders, the feature tiers, the evaluation harness, the adversarial and pre-test reviews, and the report draft. The other 19 used Claude Sonnet 5.5 (`claude-sonnet-5-5`): data acquisition, the compute-feasibility tests, the package skeleton, the embeddings, code reviews and their fixes, the CLI plumbing, the panel-thinning and robustness experiments, the release refit and benchmark, and the retrospectives. For this write-up, an Opus agent drafted the briefing and Sonnet agents produced the figures and the repository export (`results/process/agents_including_writeup.tsv`). Everything ran on one Apple MacBook Air laptop (M4, 10 CPU cores, 24 GB unified memory, no discrete GPU).
>
> Currently none of the results have been verified by a human. We plan to verify the agent's methodology and benchmarking, and we will update this section accordingly.
>
> The rest of this repository is entirely AI generated. It includes code for implementation of the predictive model and a AI-written report in the style of a journal article, detailing the final method, the path followed to reach it, and benchmarking against other methods.

## The task

Two paralogs (gene copies created by a duplication) are synthetic lethal (SL) when a cell survives losing either one but not both. The agent was asked to build an installable command-line tool that takes only a eukaryote proteome FASTA (one protein per gene) and returns every paralog pair it finds, each with an SL score.
- **Features:** computed from sequences alone. Expression, interactions, essentiality, annotations, and the agent's own knowledge of specific genes were all forbidden.
- **Training data:** known SL pairs from budding yeast and human.
- **Evaluation:** AUPRC on a held-out test set, which had to clearly beat sequence identity alone. Unseen species (plants, insects) were the eventual target.
- **Constraint:** everything had to run on one laptop.

The full prompt is in [PROMPT.md](PROMPT.md).

## Where to start

| document | what it is |
|---|---|
| [docs/briefing/BRIEFING.pdf](docs/briefing/BRIEFING.pdf) ([Markdown](docs/briefing/BRIEFING.md)) | Research briefing for a broad scientific audience: what was done, why, what came out, and the lessons for agentic science. Start here. |
| [REPORT.md](REPORT.md) | The agent's final technical report: approach, development and test results, runtime, limitations, follow-ups, reproduction |
| [PROMPT.md](PROMPT.md) | The task as given to the agent, verbatim |
| [LITERATURE.md](LITERATURE.md) | Literature review written before the design (detailed notes in `docs/lit/`) |
| [DESIGN.md](DESIGN.md) | Pre-registered study design (v3, after three review rounds): labels, features, splits, metric, success criterion |
| [RESEARCH_TREE.md](RESEARCH_TREE.md) | The research record: milestone summaries, every decision note (D0–D23), and the tree of ideas tried |
| [LESSONS.md](LESSONS.md) | Process rules the agent learned, each with its evidence |
| [DATA.md](DATA.md) | Every dataset used: source, version, checksum, licence |
| [reviews/](reviews/) | Design reviews (three rounds), independent result reviews, and retrospectives |
| [docs/AGENT_RULES.md](docs/AGENT_RULES.md) | The rules every delegated agent had to follow |
| [docs/USAGE.md](docs/USAGE.md) | Tool documentation: install, setup, predict, output columns |
| [mvp/README.md](mvp/README.md) | Minimal runnable version of the final method, with an example |

## Repository layout

- **Root:** the research record (prompt, literature, design, research tree, lessons, data sources, report).
- `src/slparalog/`, `scripts/`, `tests/`, `config/`, `pyproject.toml`: the tool and every script used to build labels, features, models and evaluations, with their tests and pinned configuration.
- `results/`: small result files only (JSON/TSV summaries of every development run, the single test access, benchmarks, and the frozen and release models). Every number in the report and briefing cites one of these files.
- `reviews/`: design reviews, result reviews and retrospectives.
- `docs/`: the briefing with its figures (PNG and SVG, plus the script that draws them from `results/`), the literature notes, and the agent rules.
- `mvp/`: the release model with a minimal run script.

Commit hashes cited in the documents (for example `a88cff1`) refer to the project's development repository. Its history is not published, because it contains data files that cannot be redistributed. This repository holds the same files as of project close, with no data, as a fresh history.

## Key findings

| finding | evidence | strength |
|---|---|---|
| On the held-out test set the tool clearly beats sequence identity: T = 0.297 [0.131, 0.395]; human AP 0.349 vs 0.106, yeast 0.541 vs 0.190 | `results/final/final_test.json` | **Strong.** Pre-registered criterion, model frozen beforehand, single test access. The comparison was pre-declared underpowered, and the yeast CI is wide |
| Much of the human gain reflects which pairs the screens chose to test: within identical study-set cells the gain is −0.020 [−0.098, 0.060] | `results/final/final_test.json`; `reviews/results/tierB_gain_review.md` | **Moderate.** Predicted by the development review and confirmed on test, but with only 67 positives in comparable cells |
| Part of the signal is transferable biology: a model trained on yeast alone lifts human AP from 0.106 to 0.259 (ΔAP 0.153 [0.085, 0.244]) | `results/final/final_test.json` | **Moderate.** One transfer direction; human→yeast is not significant |
| Transfers to an unseen animal: *C. elegans* ΔAP +0.178 [0.064, 0.349] | `results/final/final_test.json` | **Weak to moderate.** Only 16 positives |
| Does **not** work on an unseen plant: *Arabidopsis* ΔAP −0.005 [−0.111, 0.097] (strict labels); identity + family size does better | `results/final/final_test.json` | **Negative result** (51 positives). Plants are not supported |
| What helps is comparison with 80 other species' proteomes: single-copy orthologs elsewhere, similar evolutionary rates, conservation (dev T_pf 0.046 → 0.143) | `results/dev/dev_primary_AB_ABC_noasym_6bb9438.json`; `results/review/tierB_gain/importance.json` | **Moderate.** Development folds; these features follow ideas from prior work (see the briefing) |
| Self-search features alone (identity, family context) mainly learn screening ascertainment: they fail the within-stratum guard | `results/dev/dev_primary_ablation_anchor_95ac397.json` | **Negative result**, development data |
| Protein language-model embeddings (ESM-2) add little (≈ +0.01) at ~1.7 h extra compute | `results/dev/dev_primary_AC_tierA_tierC_a1917e0.json` | **Negative result.** Not adopted |
| Identity-anchored models, deconfounded training and panel-masking augmentation did not help | `results/dev/*.json` (ablation, dc, thinning) | **Negative results**, development data |
| Full human proteome (20,770 proteins): 99.8 min, 2.52 GB peak memory, on a laptop | `results/bench/e2e_human.json` | **Measured once.** The machine was not otherwise idle |
| Tolerant of fragmented or split gene models. With all simulated defects combined, the human margin over identity roughly halves | `results/dev/dev_robustness_GBDT_AB_17.json` | **Exploratory.** Simulated defects, development labels |
| Agentic process: independent reviews caught real problems before they mattered (ascertainment, a wrong premise and three fixes before the test access), while guards and pre-registration kept the test single-use. Slips still happened: runs on uncommitted code, and a rule refined after development results | `reviews/`, `RESEARCH_TREE.md`, `LESSONS.md` | **Qualitative.** One project, not yet human-verified |

## Data, licences and reproduction

**No data files are included in this repository.** Proteomes, labels, raw screen tables, panel databases, model caches and feature tables are all excluded by `.gitignore`. The public export also drops result files that contain per-pair label rows, and it redacts the per-pair rows of one review table. Only small summary results and the two trained models (decision trees, no data rows) are included.

[DATA.md](DATA.md) lists every source with its licence. These sources do **not** allow redistribution of their tables, so the derived label calls built from them are not shared either:
- Parrish et al. 2021 and Klingbeil et al. 2024: CC BY-NC-ND 4.0.
- Ito et al. 2021 and CHyMErA (Gonatopoulos-Pournatzis et al. 2020): publisher copyright, with no open licence on the supplementary tables.
- Flister et al. 2025: preprint supplementary tables.
- Costanzo et al. 2016 (TheCellMap): no licence stated; used for research only.
- Ryan et al. 2012, *S. pombe* E-MAP: publisher copyright.

All other sources are open (CC BY 4.0, CC0, MIT); see DATA.md.

**Use the final method** (details in [mvp/README.md](mvp/README.md)). Requirements: Python 3.12 and MMseqs2 on `PATH` or at `bin/mmseqs`.

```bash
pip install -e .
slparalog setup --panel              # downloads and checksums the 80 reference proteomes, builds search databases
slparalog setup --model mvp/model    # installs the release model
slparalog predict proteome.fasta -o pairs.tsv -t 8
```

**Reproduce the study.** The full pipeline is in [REPORT.md §11](REPORT.md):
1. Run `scripts/data/download.py`. Some sources need a manual browser download; DATA.md says which.
2. Build the labels and splits.
3. Compute the features.
4. Run the development reports.
5. Train the frozen model.

Run the test-suite with `pytest -m "not data"`. Tests marked `data` need the downloaded data. The single test-set access (`scripts/eval/final_test.py`) has already been made, and by design it refuses to run again.
