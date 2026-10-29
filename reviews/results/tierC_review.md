# Tier C review (independent reviewer)

Verdict: **ACCEPT WITH FIXES** (no correctness defect found in the features; fixes are test coverage / documentation).

No labels read, no metrics computed. Script: scratchpad/rev_c.py (numpy only, does not import tier_c.py).

## 1. Independent recomputation (full-proteome float64 rows)
Yeast, 800 random pairs; human, 300 random pairs (seed 12345). Species worm/arabidopsis: IDs, cache, row alignment checked (100 pairs, same script).
| feature | yeast max abs diff | human max abs diff |
|---|---|---|
| tc_cos / tc_euc / tc_cos_c | 7e-16 / 0 / 1e-15 | 8e-16 / 0 / 9e-16 |
| tc_rank_min, tc_rank_max, tc_mnn | 0 | rank_max: 1 pair off by 1 (see below); others 0 |
| tc_above_frac_min/max | 0 | 0 / 4.8e-5 (same pair) |
| tc_gap_min/max | 9.6e-7 | 9.8e-7 (float32 noise) |
| *_pct (all 4, full set) | 0 | 0 |
The single human mismatch: nearest competitor lies 1.9e-8 from the (cos + 1e-6) tie threshold; the implementation compares in float32, the reference in float64. Benign, but the "tie within 1e-6" definition is only reproducible to ~1e-7 (float32 noise); features are not bit-reproducible across BLAS/thread layouts at such boundary cases.

## 2. ID alignment
- For all four species: ids.txt == fasta record ids (order), unique; every .npy row equals (atol 1e-3, float16 cache) the cache shard row keyed by sha256 of that row's own sequence: 0 mismatches, 0 missing.
- tierC.parquet gene_a/gene_b equal tierA.parquet gene_a/gene_b row-for-row (all species); gene_a<gene_b, no duplicate pairs.
- Hand-recomputed features from the gene IDs via my own id->row dict match the stored values (item 1), so no off-by-one / scrambling.
- meta: registry_hash matches current registry for all species (output printed above in session); cache_check "identical"; meta git src_dirty=True and HEAD differs yeast vs human (provenance weakness only).

## 3. Rule compliance
Inputs are the proteome embeddings only; centring mean and neighbourhoods use the same species' proteome; no labels, annotations, other-species data, per-gene coordinates in output. No feature constant (nunique >> 1; rank/pct features 201 levels, mnn 2). Percentiles computed over all called pairs (rank average, pct) - matches registry/Tier A convention (verified exactly). torch not imported (test passes).

## 4. Rank/denominator/cap/symmetry
Ties share best rank (eps 1e-6), denominator n-2, cap 201 verified numerically (min 1, max 201). rank_min<=rank_max, mnn==(rank_max==1), above_frac==0 iff rank_max==1: all hold. Min/max construction is symmetric; swap-orientation test passes.
Observation: cap saturates heavily in human (74.7% of rank_max, 59% of rank_min == 201; tc_rank_max_pct has 62.7% of mass at one value; arabidopsis/worm similar). Not a bug (above_frac is uncapped) but rank features carry little information there.

## 5. Tests (9 pass, 84 s)
Good: brute-force reference, toy known answers, ties, cap, shift, block/thread, row-permutation, orientation, no-torch, pipeline hook vs cache.
**Would NOT catch a scrambled-ID bug in scripts/features/run_tierC.py**: the id->row mapping (pos dict, ta gene_a/gene_b -> i/j) is untested; the permutation test relabels rows and indices consistently, so it is blind to a wrong mapping. The embedded ids/hash checks in run_tierC guard the npy side at run time but are not tested. The pipeline hook test is partially protective (compares against a cache lookup by sequence hash).

## Required fixes
1. Add a test (toy fasta + synthetic embeddings, or a run_tierC function refactor) asserting that features for a named pair equal compute() on the rows of those two ids, and that a row-shifted/shuffled .npy or ids.txt is rejected.
2. Add a run_tierC output assertion that tierC gene_a/gene_b equal Tier A row-for-row (currently true by construction, not asserted/tested).
3. Document the float32 tie-boundary non-reproducibility (or evaluate the eps comparison in float64) and heavy rank saturation in human.
4. Re-run Tier C from a clean tree (meta shows src_dirty=True, differing HEADs) or record the source hash.
Cosmetic: registry description of tc_above_frac_min is phrased as a count then "/ (n-2)".
