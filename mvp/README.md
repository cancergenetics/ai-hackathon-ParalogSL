# slparalog MVP: synthetic-lethal paralog pair ranking from a proteome FASTA

`slparalog` takes one proteome FASTA (one protein per gene) and returns every paralog pair it calls in that
proteome, each with an **SL score**: a ranking of how likely the pair is to be synthetic lethal. It uses only
sequence-derived features (self-search of the proteome with MMseqs2, global alignment identity, and
duplication-age features from a panel of 80 reference proteomes). It uses no expression, interaction,
essentiality or annotation data. The model shipped here is the winning method of the project: a gradient-boosted
tree model on Tier A + Tier B features (`model_type: gbdt`, `feature_set: A+B`, 73 features). It is the **release
model**: the frozen, tested specification refit on all human and yeast labels (development + test). The refit itself
has no held-out evaluation. Its expected performance is that of the identical frozen model (`results/frozen/model`),
which was trained on development labels only and tested once. The full method, validation and limits are in
`REPORT.md`.

## Requirements

* Python 3.12.
* The MMseqs2 binary, either at **`bin/mmseqs` under the repository root** (`mmseqs.binary` in
  `config/caller.yaml`) or on your `PATH` (used when `bin/mmseqs` is absent). The build pinned for the reported results is in
  `config/tools.yaml`. A different build still runs, but the output header then carries the warning
  "MMseqs2 binary differs from the pinned version/checksum".
* About 2.2 GB of disk for the built reference panel (`~/.slparalog` or `$SLPARALOG_DATA`; setup refuses to start with less than 8 GB free), plus about 1.5 GB of
  scratch space per human-size run. Peak memory in the benchmark runs was about 2.5 GB.

## Install

From the repository root:

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -e .
slparalog --help
```

## One-off setup

```bash
# 1. Reference panel: downloads the 80 reference proteomes (UniProt / Quest for Orthologs; about 0.3 GB compressed,
#    the sizes in config/panel_sources.json sum to 0.31 GB),
#    checks every file against the SHA-256 recorded in config/panel_sources.json, builds the MMseqs2
#    databases (about 2.2 GB on disk). -t sets the threads for the database build.
slparalog setup --panel -t 8

# 2. Install the trained model that ships with this repository.
slparalog setup --model mvp/model
```

`slparalog setup` with no flag would also fetch ESM-2 weights (`--plm`). The winning model does not use them
(Tier C is not part of it), so skip that. Use `--data-dir DIR` on both `setup` and `predict` to keep the data
somewhere other than `~/.slparalog`.

## Predict

```bash
slparalog predict your.fasta -o pairs.tsv -t 8
```

`your.fasta` has one protein per gene (use one isoform per gene), and the first word of each header is the gene
ID that appears in the output. The model installed last by `setup --model` is used; pass `--model mvp/model` to
choose it explicitly. Other options (`--data-dir`, `--workdir`, `--keep-work`, `--features-only`,
`--allow-sklearn-mismatch`, `-v`) are listed by `slparalog predict --help`.

`mvp/run_example.sh` downloads the yeast reference proteome (checksummed) and runs setup and predict end to end.

## Output

`pairs.tsv` is a tab-separated table. Lines starting with `#` are a header block (tool version, model ID, input
checksum, run time, QC warnings), so read it with `pandas.read_csv(path, sep="\t", comment="#")`. One row per
called pair, sorted by rank:

| column | meaning |
|---|---|
| `gene_a`, `gene_b` | the two paralogs (record IDs from your FASTA; `gene_a` sorts before `gene_b`) |
| `sl_score` | ranking score in [0, 1]; higher means more likely synthetic lethal |
| `sl_rank` | 1 = highest score; ties share the best rank |
| `ident_global` | global Needleman-Wunsch identity of the pair (BLOSUM62, gap open/extend 11/1, end gaps penalised): identical positions / alignment columns |
| `fragment_suspect` | 1 if either protein is shorter than half the median length of its called partners (likely a fragment) |
| `near_identical` | 1 if the local identity is at least 0.99 (possible allelic or isoform copy) |
| `n_dir_called` | number of search directions (1 or 2) in which the pair was called |
| `tb_n_outgroups` | number of outgroup triplets used for the asymmetry features (0 to 3) |
| `tb_fragment_suspect` | 1 if either protein is shorter than half the median length of its best hits in the panel species |

**`sl_score` ranks pairs; it is not a probability.** It is not calibrated, and the same value means a different
absolute likelihood in different species, because the fraction of paralog pairs that are synthetic lethal
(prevalence) differs between species and label sources. Use it to prioritise pairs within one proteome. Do not
read 0.8 as "80% chance".

## Runtime

The cost is dominated by the Tier B panel search. On a 10-core Apple M4 laptop (`-t 10`), the full human proteome
(20,770 proteins, 686,225 called pairs) took 5,990 s, about 100 minutes
(`results/bench/e2e_human.json`; the machine was also loaded by other jobs during this run). The yeast proteome
(6,071 proteins, 18,014 pairs) took about 1,115 s, about 19 minutes (`results/bench/e2e_yeast.json`). Time scales
with the number of called pairs, so expect a small proteome to take minutes.

## Scope and honest limits

* **Validated in human and yeast** (frozen model, held-out test labels, single access). The numbers, intervals and
  caveats are in `REPORT.md` and `results/final/final_test.md`. A large part of the human gain reflects which pairs
  screens chose to test (see the briefing). The release model shipped here is the same specification refit with the
  test labels added, and has not been evaluated separately.
* **Worm (C. elegans): acceptable.** The gain over a simple identity baseline transfers to worm, but the worm test
  set is small (16 positives).
* **Plants are not supported.** In Arabidopsis the model shows no gain over identity, and a simple baseline scores
  higher on some tiers. Do not use the scores for plant proteomes.
* Species not listed here have not been evaluated. The score depends on a reference panel of 80 proteomes, so very
  young duplications in species far from the panel are unresolved (the output header warns about this).
* The labels behind the model come from published double-knockout screens (human) and double-mutant data (yeast).
  The score predicts the label as defined there, not a guaranteed phenotype in your system.
