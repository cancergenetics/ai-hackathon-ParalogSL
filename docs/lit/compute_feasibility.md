# Compute feasibility on Apple M4 (10 cores, 24 GB) — measured 2026-10-27

All numbers are **measured on this machine** unless marked "estimate". Scripts: `scripts/bench/` (raw logs: `logs/bench/`).

**Important caveat — macOS Low Power Mode was ON for every benchmark** (`pmset -g` -> `lowpowermode 1`; the first runs, i.e. all
MMseqs2/DIAMOND runs, were on battery; later ones on AC, Low Power Mode still on). Effective MPS GEMM throughput was only
~0.75 TFLOPs (4096^2 matmul, fp32 and fp16 alike) and `torch.get_num_threads()` = 4 (performance cores). Two system daemons (`duetexpertd`,
`siriknowledged`) also used ~1 core each. I did not change system settings. Treat all times as **pessimistic**; re-run
`scripts/bench/run_all_embed.sh` with Low Power Mode off (plugged in) before locking the compute budget. Rule of thumb measured here for
transformer encoders on MPS: time ~= 2 x params x residues / 0.75e12 FLOP/s (e.g. 650M model x 11.5M residues ~ 5.4 h).

## 1. Summary table (full human proteome = 20,652 proteins, 11.46 M residues)

| Tool / model | Version | Install route | License | Time, human proteome (measured -> extrapolated) | Peak memory |
|---|---|---|---|---|---|
| MMseqs2 all-vs-all human, `-s 5.7` (default), `--max-seqs 1000`, 10 thr | 18-8cc5c | GitHub release `mmseqs-osx-universal.tar.gz` -> `bin/mmseqs` (arm64 slice) | MIT | **278 s** (1.09 M hits) | 1.15 GB RSS |
| MMseqs2 human vs human `-s 7.5`, `--max-seqs 1000` | " | " | " | **615 s** (1.29 M hits) | 1.48 GB |
| MMseqs2 human vs human `-s 7.5`, `--max-seqs 200` | " | " | " | **251 s** (0.75 M hits) | 1.31 GB |
| DIAMOND human vs human `--sensitive` / `--very-sensitive` / `--ultra-sensitive`, `-k 200` | 2.2.8 | GitHub release `diamond-macos.tar.gz` (universal, has arm64) -> `bin/diamond` | GPL-3.0 | **37 s / 60 s / 133 s** | 0.61 / 0.81 / 1.9 GB |
| MMseqs2 human vs QfO eukaryote panel (51 spp, 976,393 prot.), `-s 5.7 --max-seqs 100` | " | " | MIT | **560 s** (1.5 M hits) | 4.2 GB |
| DIAMOND human vs same panel, `--very-sensitive -k 200` | " | " | GPL-3.0 | **1734 s (29 min)** | 8.5 GB RSS, **19.1 GB footprint** (near the limit) |
| Foldseek | 10-941cd33 | GitHub release `foldseek-osx-universal.tar.gz` -> `bin/foldseek` | GPL-3.0 | not benchmarked for search (only ProstT5 3Di, below) | — |
| ESM-2 t6 8M | HF transformers 5.21, torch 2.14.1 (pip, in venv) | `facebook/esm2_t6_8M_UR50D` | MIT | MPS fp16 28 min; CPU 24 min | <1.5 GB RSS |
| ESM-2 t12 35M | " | `..._t12_35M_UR50D` | MIT | MPS fp16 66 min; CPU 65 min | <2 GB |
| ESM-2 t30 150M | " | `..._t30_150M_UR50D` | MIT | MPS fp16 **104 min**; fp32 113 min; CPU 230 min | MPS driver 1.5 GB; RSS 1.3 GB |
| ESM-2 t33 650M | " | `..._t33_650M_UR50D` (2.6 GB) | MIT | MPS fp16 **326 min (5.4 h)**; fp32 388 min; CPU (n=50) 561 min | MPS driver 2.2 GB (fp16), RSS 5 GB at load |
| ESM C 300M | " (native `esmc` in transformers; no remote code) | `biohub/ESMC-300M` (1.3 GB; ungated) | HF card: MIT + third-party notice (earlier EvolutionaryScale releases: Cambrian Open (300M) / Non-Commercial (600M) — **re-check before release**) | MPS fp16 **158 min**; fp32 216 min; CPU (n=50) 390 min | RSS 2.8 GB; MPS driver 1.9 GB |
| ESM C 600M | " | `biohub/ESMC-600M` (2.3 GB; ungated) | as above | MPS fp16 249 min (4.1 h) | RSS 4.3 GB |
| ProtT5-XL-half (encoder) | " | `Rostlab/prot_t5_xl_half_uniref50-enc` (2.4 GB) | AFL-3.0 per ProtTrans repo (HF card has no tag) | MPS fp16 **653 min (10.9 h)**; CPU fp32 (n=30) 1250 min | RSS 3.5 GB; MPS driver 6.5 GB |
| ProstT5 encoder (embeddings only) | " | `Rostlab/ProstT5_fp16` (5.6 GB; the fp32 repo is 11.3 GB, not downloaded) | MIT | MPS fp16 660 min (11 h) | RSS 3.4 GB; MPS driver 7.6 GB |
| ProstT5 AA->3Di via HF `generate` (sampling, bs<=16) | " | same | MIT | MPS fp16 **0.05 seq/s -> ~125 h** (n=100) — not viable | MPS driver 13.5 GB |
| ProstT5 AA->3Di via **Foldseek `createdb --prostt5-model`** (GGUF f16, Metal) | Foldseek 10 | `foldseek databases ProstT5 ...` (2.1 GB `prostt5-f16.gguf` -> `models_cache/foldseek/`) | GPL-3.0 (Foldseek); MIT (weights) | **667 s / 200 proteins = 0.30 seq/s -> ~19 h** extrapolated (n=200) | 2.6 GB RSS |
| ESMFold v1 / ESMFold2(-Fast) | — | `facebook/esmfold_v1` 8.4 GB; `biohub/ESMFold2-Fast` 26 GB (**not** downloaded, would not fit 24 GB in fp32) | MIT | **not run; estimate: infeasible** (see section 4) | — |

Embedding benchmark details: 200 proteins drawn at random (seed 0) from the human proteome (real length distribution; 107k residues), mean-pooled
last hidden state excluding special tokens, length-sorted token-budget batches, one process per config. Full matrix and per-run JSON in
`logs/bench/embed_results_main.jsonl`. seq/s: t6 13.4 (MPS fp32), t12 5.4, t30 3.4, t650M 1.1, ESMC-300M 2.3, ESMC-600M 1.4, ProtT5 0.55, ProstT5-enc 0.54.
Notes: (a) CPU rows for the big models use n=30–50 so extrapolation is noisy; (b) varying `tok_budget` 1k-4k changed nothing for t6/t12 (compute-bound);
(c) **fp16 vs fp32 on MPS gives identical embeddings** (mean per-protein cosine >= 0.99999, pairwise-similarity Spearman 1.0000; checked for t6, t12, t30, t650M,
ESMC-300M; `scripts/bench/compare_emb.py`) and is 8–35% faster, so use fp16; (d) small ESM-2 models are *not* proportionally faster (t6 only 4x faster than t30; CPU = MPS speed
for t6/t12) because they are overhead/attention dominated.

**Length policy used** (11.2% of human proteins > 1022 aa, 2.25% > 2000 aa, max 34,350): split any protein longer than the model window
(1022 aa for ESM-2/ProtT5/ProstT5, 2046 for ESM C) into ceil(L/window) *equal-length* chunks (23,730 chunks for 20,652 proteins at 1022; 21,202 at 2046), embed independently, length-weighted
average of chunk means. No residues are dropped, so cost scales with total residues (no truncation savings). Alternative truncation to the first 1022 aa would save only ~12% of tokens
and loses the C-terminus; not recommended. ProtT5/ProstT5 (relative positions) could take longer windows but attention memory grows quadratically; chunking is safer on 24 GB.

## 2. Homology search notes

* Human proteome (UP000005640, release 2026_03, 20,652 proteins; yeast UP000002311, 6,066 proteins) downloaded to `data/raw/proteomes/` (README with URLs, date, sha256).
* MMseqs2 all-vs-all human is cheap (4–10 min, ~1.5 GB); cost is dominated by the alignment stage and grows with `--max-seqs`/sensitivity because large families
  (zinc fingers, olfactory receptors, kinases) generate huge hit lists. For paralog detection `-s 7.5 --max-seqs 200` (4 min) is a sensible default; verify recall against `--max-seqs 1000` later (not yet done).
* DIAMOND is 5–10x faster within-species (`--very-sensitive` 1 min), but against the 976k-protein panel it was 3x slower than MMseqS2 and its peak footprint was 19 GB -> use MMseqs2 for the panel. DIAMOND remains an option for the self-search (note: `-k 200` caps hits).
* Foldseek has an official macOS universal binary (works, Metal-enabled ProstT5 inference). No bioconda needed for any of the three.
* **QfO reference proteomes** (release 2026_02, https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/): 81 species total (51 eukaryotes, 30 prokaryotes);
  full tarball `QfO_release_2026_02.tar.gz` 2.2 GB (+ 20 GB `_additional` isoforms; not needed). I downloaded only the 51 eukaryote canonical `*.fasta.gz` = **309 MB gz / 588 MB, 976,393 proteins** (matches README counts) into `data/raw/qfo/`
  (README, SHA256SUMS, species_list.tsv). Species mix: plants/algae (Arabidopsis, rice, maize, moss, Chlamydomonas), fungi (S. cerevisiae, S. pombe, Aspergillus, Neurospora, Candida...), few arthropods (Drosophila, Anopheles, Tribolium, Ixodes), C. elegans, many vertebrates, several protists. Panel is vertebrate-heavy; for plants/insects add more
  invertebrate/plant reference proteomes from UniProt if needed (cost scales ~linearly: ~11 s average per species at s5.7 -> 100 species ~20 min). Panel includes human itself (exclude the query species in real use).
* Panel search cost for one new proteome: **~9.5 min, 4.2 GB** with MMseqs2 `-s 5.7 --max-seqs 100`. For the training species (yeast, human) the panel search is a one-off (~10 min + ~3 min).

## 3. Protein LM / 3Di notes

* ESM C weights moved to the `biohub` HF org; `AutoModel.from_pretrained("biohub/ESMC-300M")` works without login or `trust_remote_code` on transformers >= 5 (remote-code shim only wraps the `esm` pip package; I read it before deciding not to use it).
* ProtT5-half and ProstT5 are 4–8x slower than ESM-2 650M per residue (1.2B-param T5 encoder; MPS fp16 T5 produced no NaNs). They are only worth it if they add accuracy over ESM C/ESM-2.
* Full-proteome 3Di prediction is the single most expensive feature: Foldseek-ProstT5 ~19 h (human), ~5.6 h (yeast, 6,066 proteins) at Low-Power-Mode speeds; HF `generate` is 6.5x slower still. Mitigations: predict 3Di only for proteins that have a paralog hit
  (typically a fraction of the proteome, to be measured), cap length (e.g. <= 1000 aa; autoregressive cost is superlinear), run with Low Power Mode off, and cache by sequence hash. Not tested: a non-autoregressive 3Di head on ProstT5/ESM embeddings.

## 4. ESMFold feasibility (not run; estimate)

ESMFold v1 = ESM-2 3B stem (6 GFLOP/token; 11.5 M residues x 6 GFLOP = 69 PFLOP -> ~25 h at the measured 0.75 TFLOPs for the language-model stem *alone*) plus a 48-block folding trunk with L x L x 128 pair
tensors (memory ~L^2; long proteins > ~1,000 aa will exceed 24 GB). Published GPU timings (~0.1–15 s per protein on A100/V100, i.e. 50–100x our effective compute) imply roughly 1–5 min per median (411 aa) protein here ->
**hundreds of hours for a proteome; infeasible**, and only borderline for a few hundred proteins. Weights are 8.4 GB (fp32 .bin); ESMFold2-Fast is 26 GB (contains the 6B ESM C) — would not fit. Recommendation: do not use predicted structures; use ProstT5 3Di (sequence-only structure tokens) instead.

## 5. Recommendations

1. **Homology**: MMseqs2 only (MIT, 18-8cc5c). Self-search `-s 7.5 --max-seqs 200` (~4 min human); panel search `-s 5.7 --max-seqs 100` against a 50–100-species QfO-style panel (~10–20 min, ~4 GB). Keep DIAMOND as optional fast pre-check.
2. **Embedding workhorse**: ESM C 300M or ESM-2 150M (t30) in fp16 on MPS: 1.7–2.6 h for a full human proteome; embed only proteins that appear in >= 1 candidate pair (likely ~half) and cache by sequence hash -> ~1 h. Use ESM-2 t6/t12 (0.5–1.1 h) for development loops. ESM-2 650M (5.4 h) / ESM C 600M (4.1 h) only as an ablation on the candidate subset; ProtT5 (11 h) probably not worth it.
3. **Structure tokens**: Foldseek `createdb --prostt5-model` (GGUF) is the only feasible route (0.3 seq/s); budget ~19 h/human proteome unless restricted to candidate-pair proteins/length cap; treat 3Di as an optional "heavy" feature tier and report accuracy with and without it. Verify with Low Power Mode off — expected speed-up unknown (possibly 1.5–2x).
4. **Total budget** for a new proteome, fast tier (MMseqs self + panel, ESM-2 150M on candidates): roughly 1–1.5 h and < 5 GB; full tier (+ ProstT5 3Di on candidates, ESM C 600M): ~10–25 h.
5. Peak memory is never the binding constraint (all < 8.5 GB RSS except DIAMOND vs panel, 19 GB footprint); **time is**. Keep batch sizes by token budget (4k–16k tokens) and MPS fp16.
6. TODO for follow-up: re-benchmark with Low Power Mode off; test MMseqs2 recall of `--max-seqs 200` vs 1000 on known paralog pairs; measure the fraction of proteins having a paralog hit (to size the restricted-embedding approach); benchmark Foldseek search with ProstT5-derived 3Di.

## 6. Reproduce

```
scripts/bench/bench_mmseqs.sh  Q.fasta T.fasta SENS TAG        # MAXSEQS=200 env to change --max-seqs
scripts/bench/bench_diamond.sh Q.fasta T.fasta --very-sensitive TAG
scripts/bench/run_all_embed.sh                                 # full embedding matrix (~3 h) -> logs/bench/embed_results.jsonl
scripts/bench/bench_embed.py --model esm2_t30_150M --device mps --dtype fp16 --n 200
scripts/bench/bench_foldseek_prostt5.sh data/interim/bench_foldseek/sample200.fasta s200
scripts/bench/compare_emb.py ; scripts/bench/make_sample_fasta.py 200 OUT.fasta
```
Disk: `models_cache/` 24 GB (HF 22 GB incl. unused ESM-2 .bin duplicates, ProstT5 5.6 GB, ProtT5 2.4 GB, ESM-2 650M, ESMC; Foldseek GGUF 2.3 GB), `data/raw` 348 MB, `data/interim` ~1.2 GB; free volume 522 GB.
Binaries in `bin/` (89 MB): mmseqs, foldseek, diamond (macOS universal; sha256 of the tarballs: mmseqs 158102d4…, foldseek 6d8b07e1…, diamond 7d17c9c7…).
