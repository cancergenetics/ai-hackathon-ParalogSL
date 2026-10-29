# Tier B independent review

Reviewer: independent agent (did not write the code). No source files modified. No SL labels read; no metric computed against SL labels.
Commits: 7a61fda, 48cd3fe, 3d5209d. Scope: src/slparalog/panel.py, features/tier_b.py, registry `_TIER_B`, `--tier-b` hook, tests/test_tier_b_toy.py, config/panel.yaml, tierB outputs for human, yeast, arabidopsis, worm. Panel searches were not re-run; cached hit tables were used. Scratch scripts are in the session scratchpad (`align_check.py`, `recompute.py`, `ds.py`, `yeast_wgd.py`).

## VERDICT: ACCEPT WITH FIXES

The implementation matches the normative spec (algorithms_comparative.md changes 5, 6, 9, 10, 11), and every number I recomputed agreed. The required fixes are test and robustness gaps, not feature-value errors. The yeast long-branch problem is a genuine limitation of a distance axis, not a bug.

## 1. Row alignment and ID mapping (measured)

- **Row alignment:** tierB.parquet and tierB.pairbits.parquet have the same `gene_a`/`gene_b` as tierA.parquet, in the same row order, in all 4 species (human 686,225; yeast 18,014; arabidopsis 914,960; worm 260,434).
  - There are no duplicate pairs, and `gene_a < gene_b` for every row.
  - The `aligned_to` sha256 in each Tier B meta matches the current tierA.parquet.
- **Query-ID mapping:** I rebuilt the canonical order (sort by sequence, then id) and the cache key from the proteome FASTA.
  - For every panel hit of the sampled queries, `qlen` equals the length of the mapped proteome sequence: 0 mismatches in 1,666,784 human and 1,148,467 yeast hits.
  - The top hit of each query in the self species (HUMAN, YEAST) is a byte-identical sequence for 99.4% (human, first 500 queries) and 99.7% (yeast, first 300). The residual is presumably identical-sequence ties, which I did not investigate. This validates the query-index to gene and target-index to sequence mapping end to end.
- **Would the tests catch a scrambled ID?**
  - Query-side scramble: yes. I mutated a copy of `load_panel_data` in the scratchpad (`np.roll` of the per-species vectors by one query index). The 19 tests that pass on the unmodified code gave 8 failures: planted_histories P/Z/X, pre_post_bounds, outgroup_asymmetry, masking, bin_nan_vs_zero and no_proteome_constant.
  - Panel-target scramble: invisible by design. Targets are only compared for equality (a == b), which is permutation-invariant. The outgroup path reads sequences by index and is covered by the asymmetry known-answer test.
  - Gap: `scripts/features/run_tierB.py` (the script that produced the training features) has no test. Its `gene_a/gene_b` to `pos` mapping and the Tier A row alignment are tested only by my measurement above.

## 2. Independent recomputation

Method: pandas-only re-implementation from the cached hit parquet files (no import of `panel.classify` or `reduce_species`), 200 random pairs each for human and yeast. It uses `d_s` from the meta and the query-side pair bits from the pairbits sidecar. The sidecar comes from the pipeline's own realignment, so it is not independently derived.

| Feature | Human | Yeast |
|---|---|---|
| `tb_depth` non-conflict | 162/162 | 187/187 |
| `tb_depth` conflict (L ≥ U, bin-weighted discordance, brute force) | 38/38 | 13/13 |
| `tb_frac_pre`, `tb_frac_post` (pre/post classification) | 200/200 | 200/200 |
| `tb_n_post_single` | 200/200 | 200/200 |
| `tb_cons_both` | 200/200 | 200/200 |
| `tb_species_rank_min`, `tb_species_rank_max` | 200/200 | 200/200 |

`d_s` itself was recomputed for all species from the RBH, G* and Kimura definitions. Human: core 1065, 77 usable species. Yeast: core 934, 79 usable species. Maximum |Δd_s| is 5e-7, and the self-excluded sets match (human: GORGO, HUMAN, PANTR; yeast: YEAST).

Not independently checked: the outgroup asymmetry features (CIGAR decoding), `tb_dupability_*`, `tb_fragment_suspect`. The asymmetry decoding is covered by the toy known-answer test and by the logged `n_cigar_nident_mismatch` counter.

## 3. Rule compliance

- **Inputs only from the proteome plus panel proteomes:** yes. There is no use of labels, annotation, ygob, ohnolog or timetree in `src/` (grep clean). `config/panel.yaml` genus lists (Sacch pre/post-WGD) drive panel construction only (`scripts/data/build_panel.py`), not feature computation.
- **Annotations only in evaluation:** YGOB, OHNOLOGS and TimeTree appear only in `scripts/eval/tierB_validation.py` (`evaluation_only: true` in its JSON output).
- **Constant features:** none constant or all-NaN in any of the 4 proteomes (27 columns each).
  - The registry test covers the toy and yeast cases.
  - Not measured: the DESIGN species-identification control (a classifier predicting species from the features). I did not find it among the Tier B outputs.
- **Self-exclusion:** one code path (`calibrate`, median RBH identity ≥ 0.98), run in `run_stage` for both training and predict. Human removes GORGO and PANTR as well as HUMAN. Nearest remaining human species: MACMU 0.011.
- **Species-specific hard-coding:** none in panel.py or tier_b.py (no species names; grep clean). The only human-derived constant is the bin-edge table. It is used only for the 1/n_bin discordance weights, and the effect is small.
  - Bin edges: "from the reviewer's human core-set values", frozen before yeast `d_s`.
  - On an unseen species the edges are applied to the query's own `d_s` scale, which is acceptable.
- **Unseen plant or insect proteome:** the code path is species-agnostic, and arabidopsis and worm ran end to end.
  - Risk 1: the core-set fallback (down to a 0.1 species fraction) can leave `|G*| < 200` with only a log entry. I saw no hard failure or warning flag.
  - Risk 2: worm has `tb_no_evidence` = 25.7% (human 0.0%, yeast 1.5%, arabidopsis 0.7%).
    - Cause: the panel's nearest nematodes are far (nearest worm species CYLNA, d = 0.38). Most pairs have no pre or post species, so depth is NaN.
    - This is a train/test distribution shift. The model must tolerate NaN-rich depth, and it should be reported.
  - Risk 3: asymmetry features (`tb_asym_*`, `tb_outgroup_ident_max`) are NaN for 94 to 98% of pairs in every species, so they are effectively unusable.

## 4. Yeast long-branch issue: genuine limitation, not a bug

Measured from the cached hits. Nearest species by `d_s`:

| Species | `d_s` | WGD status |
|---|---|---|
| TORDC | 0.319 | pre-WGD |
| CANGA | 0.335 | post-WGD |
| KAZAF | 0.340 | post-WGD |
| HUIN7 | 0.368 | post-WGD |
| TX1245769 | 0.395 | pre-WGD |
| EREGS | 0.410 | pre-WGD |

Pre- and post-WGD species interleave inside 0.32 to 0.41 (the B3/B4 edge is at 0.35), so a distance-ordered axis cannot encode the WGD split.

On the 522 called YGOB ohnolog pairs (evaluation only):

- Class in TORDC: pre 0.8%, post 21.6%, undetermined 76.4%.
- Class in CANGA, KAZAF, HUIN7: pre only 6 to 10%, post 24 to 33%.
- In those three post-WGD species 76 to 80% of ohnolog pairs have both members hitting the same gene (a == b). The median bits of A to the species' ortholog is 405 to 433, against 370 for the A-B bits, so they are classified undetermined or fail the pre test.
- Only 17% of ohnologs have any pre species. 83% are `younger_than_nearest`, so depth = U/2 and the median ohnolog depth is 0.198.
- The depth estimate is a distance-order artefact: it is not a pre-/post-WGD boundary estimate.

Two distinct causes:

1. Species-level rate heterogeneity: fast-evolving post-WGD lineages (Kazachstania, Huiozyma, C. glabrata) sit farther from S. cerevisiae than Torulaspora, even though Torulaspora branched earlier. A median-identity distance cannot fix this.
2. Gene-level: ohnologs evolve faster or more asymmetrically than G*, and gene loss in post-WGD species is common. The classification is correct by the spec and works as designed.

Tier B still carries the WGD signal through counts rather than depth. Ohnolog versus other AUROC: `tb_n_post_single` 0.91, `tb_frac_post` 0.89, `tb_frac_pre` 0.16 (inverted, i.e. 0.84), against 0.20 for `tb_depth` (inverted, 0.80) (results/features/tierB_validation.json). The identity-stratified AUROC for `tb_n_post_single` is 0.90 (yeast YGOB).

Within-scope proposals (no code written):
- **(a) Rate-corrected split depth from the panel alone.** Compute panel-panel RBH distances once at setup, on a fixed ortholog set. Replace d(q,s) with the Gromov product (d(q,s) + d(q,O) − d(s,O)) / 2, with O drawn from a fixed set of distant outgroups (median over O). This removes species-private branch length, the three-point/relative-rate correction, and needs no annotation.
- **(b) Topology-derived order.** Build an NJ or UPGMA tree on the panel species (panel-only, fixed). Order species by branching order from the query's attachment point, and define L and U on that rank axis instead of d.
- **(c) Cluster-level calls.** Collapse panel species into fixed distance clusters (from panel-panel distances) and classify per cluster by majority. This reduces noise and interleaving. Add rank-based features such as "pre/post state in the nearest 1, 2 and 3 clusters".
- **(d) Pair-rate-aware `post`.** Use the existing `tb_rel_rate` to condition `post` on a gene-level rate estimate, to separate fast ohnologs from loss.
- **Not recommended:** a Saccharomycetaceae WGD genus list inside the features. It is species-specific and would not transfer to the unseen test species.

## 5. Training versus fresh `predict` differences

Direct test: I ran `pipeline.run(yeast.fasta, tier_b=True, features_only=True, threads=4)`, the same stage as `slparalog predict --tier-b`, against the stored tierA.parquet and tierB.parquet.
- It matched on all 18,014 pairs, on all 27 Tier B columns and all numeric Tier A columns, with tolerance 1e-9 and no differing cells. Only the pairs and values are compared; the row order of the two files was not.
- It used the cached panel hits (80/80 species cache hits), so it did not test search determinism or thread-count dependence.

Remaining risks:
1. **Cache dependence.** The training features were produced by `scripts/features/run_tierB.py` with `species_jobs=3` and 8 threads. The cache key covers the sequences, settings and tool sha, but not the thread count. MMseqs2 tie-breaking under `--max-seqs 50` with different thread counts was not tested by this review.
2. **Capped queries.** 2,110 yeast queries remain capped after one re-search round (`max_research_rounds: 1`). This affects second-hit and `dup_count` for large families, identically at training and inference.
3. **Proteome-wide dependence by design.** `d_s`, G* and the percentile columns (`*_pct`) depend on the whole input proteome. A subset or differently isoformed FASTA changes every feature value. This is by design, but undocumented in the CLI.
4. **Fresh-machine requirements.** The panel DB and hit cache sit under `PROJECT_ROOT/data/interim`. A fresh install needs the panel FASTA and DBs. The first run costs about 55 min for yeast (panel_search_s 3127 s, results/bench/tierB_yeast.json). `pipeline.run` accepts `panel_cache`, but the CLI does not pass it.

## Required fixes

1. Add a test that runs `scripts/features/run_tierB.py`'s alignment logic on a toy proteome whose FASTA order differs from the canonical order, and checks the Tier B rows against `pipeline.run`. This covers the untested row alignment and `pos` mapping.
2. Add an end-to-end regression test or script that compares `pipeline.run --tier-b` against stored tierB.parquet on a fixed small proteome. Log the thread and `species_jobs` values used for each run in the Tier B meta.
3. Raise a warning flag when the core set `|G*| < 200` after the fallback, and when the nearest retained species has `d_s > 0.3` (worm-like proteomes). Write it into the output meta and the QC columns.
4. Report the Tier B NaN and missingness rates per species (`tb_no_evidence`: worm 25.7%) and the asymmetry feature coverage (2 to 6%) in the dev report. Either drop `tb_asym_*` from the model features or document it as near-empty.
5. Document the yeast long-branch limitation in DESIGN or REPORT: `tb_depth` does not encode the WGD split, and the WGD signal is carried by `tb_n_post_single` and `tb_frac_post`. Evaluate proposal (a) or (c) as an optional follow-up. It would need a pre-declared adoption rule and a schema change.
6. Implement and report the DESIGN species-identification control (species AUROC on the Tier B features).
7. Note in the code or DESIGN that the bin edges are human-derived and used only for the discordance weights.

## Unverified

- Search determinism across thread counts.
- The outgroup asymmetry CIGAR decoding, `tb_dupability_*` and `tb_fragment_suspect`, independently of the code under review.
- Species-identification AUROC.
- Whether the `d_s`-vs-TimeTree Spearman is adequate for the unseen species: it is 0.78 for human, 0.87 for yeast, 0.77 for worm and only 0.59 for arabidopsis (results/features/tierB_validation.json).
