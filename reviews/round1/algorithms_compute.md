# Round 1 review — Algorithms & compute

Reviewer role: algorithms and compute (homology search, MMseqs2/Foldseek/ProstT5, PLMs on Apple silicon, gradient boosting, CLI packaging).
Inputs read: PROMPT.md, DESIGN.md v1, LITERATURE.md, docs/lit/compute_feasibility.md, docs/lit/methods.md (Dennler/De Kegel sections), the existing benchmark outputs in `data/interim/bench_mmseqs/`, and the venv.

**New measurements.** I ran a few light checks on outputs that already exist. There was no new search and no download. The scripts are in my scratchpad and can be re-run with `.venv/bin/python -I`:
- `<scratchpad>` (panel saturation)
- `<scratchpad>` (components)
- `<scratchpad>` (coverage of candidates)
- `<scratchpad>` (self-exclusion rule)

The PI should have a project script reproduce any number below before relying on it.

| Check (human, MMseqs2 18-8cc5c, E ≤ 1e-3) | Result |
|---|---|
| Unordered non-self pairs, `-s 7.5 --max-seqs 1000` | **691,190** pairs. Only 578,796 (84%) are found in both directions; 112,394 are one-directional. |
| Same with `--max-seqs 200` | 454,609 pairs: **236,581 pairs (34%) lost**. 1,515 genes saturate at 199 hits. `--max-seqs 200` is unusable. |
| Largest per-gene hit list at `--max-seqs 1000` | 728 hits. The *aligned* cap is not reached, but prefilter saturation is unknown. |
| Pairs involving a gene with > 100 paralogs | **599,401 / 691,190 (87%)**: the candidate set is dominated by giant families. |
| Pairs with paralog rank ≤ 5 in either direction | 46,733 |
| Genes in ≥ 1 candidate pair | 16,933 / 20,652 genes, which hold **88% of residues** |
| Single-linkage components, all E ≤ 1e-3 edges | Largest component = **10,650 genes (52% of the proteome)** |
| Same, both coverages ≥ 0.3 (the §6.1 rule) | Largest component = **5,812 genes (28%)**; next is 433 |
| Same, both coverages ≥ 0.5 | Largest component = 941 genes (4.6%) |
| Proteins with no self-hit in the self-search | 28 |
| Human vs QfO panel, `--max-seqs 100`: queries at the 100-hit cap | **10,132 / 20,617 (49%)** |
| Self-exclusion rule (≥ 50% of proteins with a hit at ≥ 98% id and ≥ 90% cov) for a human query | HUMAN 99.8%, **PANTR 66.9%, GORGO 57.0% (both excluded)**, MACMU 31.7% (kept; median best-hit id 96.8%), MOUSE 7.5% |
| Same E-value hit, self-search vs panel search | 1.5e-73 vs 6.2e-72: E-value scales with database size |
| `import torch` then `LGBMClassifier.fit()` (multi-threaded) in this venv | **Segfault (exit 139)**. It works if lightgbm is imported first, with `n_jobs=1`, or with sklearn `HistGradientBoostingClassifier`. Two libomp copies are loaded (torch/lib and sklearn/.dylibs). |

## 1. Overall verdict

The architecture is sound and fits the machine: MMseqs2 self-search plus panel search, cheap family/rank/conservation features, optional PLM and 3Di tiers, a GBDT, and one feature code path for training and inference. The run sizes that matter (~700k candidate pairs for human, ~1 M panel hits) are small for this hardware.

However, several specified parameters and rules break on real data:
- The panel search cap truncates conservation features for half the proteome.
- Single-linkage "families" collapse half of human into one component, which affects both the component-size feature and the §6.1 splits.
- The self-exclusion rule leaves the human query with a 97%-identical neighbour, while yeast's nearest neighbour is about 40% identical.
- The planned 3Di TM-score/LDDT features cannot exist without Cα coordinates.
- Torch and LightGBM crash when used in the same process.

All of these can be fixed cheaply and within scope.

## 2. Concrete changes

### Paralog calling

1. **(must) Keep `--max-seqs ≥ 1000` for the self-search. Detect prefilter saturation and re-search only the saturated queries with a 5× higher cap.**
   - Why: `--max-seqs 200` (the compute doc's suggested default) loses 34% of pairs. Plant and insect proteomes have larger families than human: Arabidopsis F-box ~700, maize 39k proteins, Trichomonas 50k.
   - How: use `mmseqs prefilter` result counts, or compare hit counts against the cap. Re-search the saturated queries with `createsubdb`. Log the number of saturated genes in the run log and freeze the policy in `config/`.

2. **(must) Make every pair feature invariant to pair order, and complete both directions explicitly.**
   - Why: 16% of human pairs are hit in one direction only, so "both directions" columns would be partly missing. A missing column would then mean "weak or asymmetric hit" in a way that depends on sensitivity and cap.
   - How: build the union of directed pairs, add the swapped orientation (`mmseqs swapresults`), and re-score all of them with `mmseqs align` using fixed parameters. This costs minutes.
   - Then feed the model only symmetric transforms: min/max, sum, |diff|, and "the longer vs the shorter protein". Never feed columns that depend on "A = lexicographically smaller accession".
   - Add a unit test that `score(A,B) == score(B,A)`.

3. **(must) Do not use the raw E-value as a feature, or as a cross-species-comparable threshold beyond calling.**
   - Why: E-values scale with database size. Measured: a 30× change between the 20k-protein self-search and the 1 M-protein panel search for the same alignment. Yeast (6k proteins) and a 100k-protein plant would sit on different scales.
   - How: use bitscore, bitscore / min(self-scores), and bitscore / max(self-scores). The E ≤ 1e-3 calling threshold can stay; optionally recompute E from bitscore at a fixed effective database size so that calling is species-comparable.

4. **(must) Guarantee self-scores.**
   - Why: 28 human proteins have no self-hit, so "bitscore normalised by self-score" is undefined for them.
   - How: use `--add-self-matches 1` or align each sequence to itself explicitly. Add a test.

5. **(should) Add a global-identity variant (identity over the full length) to the features and to B1.**
   - Why: MMseqs `pident` is identity over the local alignment columns. 263,625 directed human hits have both coverages < 0.3; these are domain-only matches with inflated local identity. The literature's "identity" (Ensembl/De Kegel) is closer to a global measure.
   - How, cheapest: `nident / max(qlen, tlen)` and `nident / min(...)` from `--format-output ...,nident,qlen,tlen`.
   - How, better: a SIMD Needleman–Wunsch on the candidate pairs. About 2·10¹¹ DP cells for human, which is minutes on 10 cores; benchmark it before adopting. B1 must be allowed to pick this variant on development folds. Otherwise "beats identity" is measured against a weak baseline.

6. **(should) Add an input pre-check for near-identical sequences.**
   - Why: an unannotated plant or insect proteome often contains un-collapsed alleles or haplotigs, isoforms or fragments. These appear as ~99–100%-identity "paralogs" and would sit at the top of the SL ranking.
   - How: detect pairs with ≥ 98% identity where the coverage of the shorter protein is ≥ 95% (`mmseqs linclust`/`clusthash`, or reuse the self-search). Flag them in the output (`possible_allele_or_isoform`) and report their count. Do not silently drop them; that would violate the "every pair" contract.

### Features

7. **(must) Replace connected-component size, and the "giant cluster forced into development" logic, because single-linkage components are degenerate.**
   - Evidence: at E ≤ 1e-3, 52% of human genes form one component; at both coverages ≥ 0.3 it is 28%. This is before adding yeast to the union graph.
   - Effect on the feature: it then mostly encodes "contains a promiscuous domain".
   - Effect on the splits: the > 5% rule would move most labelled kinase/GTPase/transporter families into development, so the test set would hold only small, isolated families.
   - Feature replacement: local family size, i.e. the number of paralogs above bitscore-ratio thresholds of 0.2/0.3/0.5. Optionally MCL/Leiden cluster size on the bitscore-ratio graph.
   - Splits (decision belongs to the statistics reviewer): measure the fraction of labelled pairs in the largest component before freezing. If it is > 20%, use a stricter linkage (both coverages ≥ 0.5 gives a largest component of 941 genes), and report the *residual* cross-split homology in the leakage audit. Residual homology = the fraction of test genes with an E ≤ 1e-3 hit to any training gene. Do not claim zero.

8. **(must) Fix the panel search so per-species quantities are not truncated.**
   - Problem: `--max-seqs 100` is global across 51 species, and 49% of human queries hit the cap. Saturated queries then show hits in a median of only 23 species, and the query's own species uses about 5 of the slots. So "number of species with a hit", "best hit of A in species S" and the Dennler "species" rank are censored, and the censoring is strongest exactly for multi-copy families.
   - Fix: search each panel species as its own database with a small cap (e.g. `--max-seqs 50` per species). Alternatively use about 8 clade databases with `--max-seqs 300`. Keep only the top-k hits per species.
   - Cost: about the same total time (the doc measured ~11 s per species). In `slparalog setup`, precompute each species' database and index once with `createindex`.
   - With per-species best hits, compute "species rank" per species ("does B outscore A's best hit in S?"). That is cap-robust.

9. **(must) Define conservation, age and duplication-age features relative to the query's measured divergence, not a fixed species tree, and make the self-exclusion rule divergence-based.**
   - The tool receives only a FASTA, so it cannot know where the query sits in the panel tree. Any phylostratum must come from sequence.
   - The current rule gives very different panel neighbourhoods for the two training species:
     - human: chimp and gorilla excluded, macaque kept at a median best-hit identity of 96.8%, ~10 mammals at 85–90%;
     - yeast: nearest neighbour at about 40%;
     - a new beetle or crop: somewhere in between.
     The same panel features would then mean different things in training and at inference.
   - Proposed rule, pre-registered:
     (a) run a pre-search of 1,000 sampled query proteins;
     (b) compute each species' median best-hit identity d_S to the query;
     (c) exclude every species with d_S ≥ τ. Pick a single τ, e.g. 90%, chosen *before* modelling and justified by typical within-genus divergence. It should be one value, not tuned on labels.
   - Then express conservation as features of the d_S distribution. Examples: "lowest d_S among species with a hit for A" (age proxy), and counts of species with a hit in d_S bins [<40, 40–60, 60–90). This is continuous, tree-free and comparable across queries. Write d_S per species to the run log.

10. **(must) Drop the TM-score/LDDT features from Tier D.** Foldseek databases built with `createdb --prostt5-model` contain only 3Di and amino-acid sequences, without Cα. Foldseek's README states that TM-score, LDDT and `--alignment-type 1` are unavailable for them. Tier D can only provide 3Di+AA alignment bits, fident and E-value, and Dennler found these strongly correlated with identity (ρ ≈ 0.86).

11. **(should) Use the ProstT5 encoder (HF, MPS fp16) plus the published CNN 3Di head instead of Foldseek's GGUF route.**
   - The ProstT5 paper's fast mode predicts 3Di with a small CNN on encoder embeddings.
   - Measured here: the HF encoder took 11 h for human, against 19 h for Foldseek-GGUF.
   - One pass then yields both the ProstT5 embedding (a Tier-C candidate) and the 3Di string, which is written to a Foldseek `_ss` database for 3Di+AA search.
   - Pin the CNN checkpoint by sha256; it is hosted at rostlab.org, not HF.
   - Before adopting it, verify that 3Di strings agree with the Foldseek route on the 200-protein benchmark sample (per-residue token agreement and search hits).

12. **(must) Correct the compute plan: restricting per-protein PLM/3Di work to "candidate proteins" saves only ~12%, not ~50%.** 88% of residues belong to genes with at least one candidate pair, and even the genes in rank-≤5 pairs number 16,886.
   - Revised Low-Power-Mode estimates for human:
     - ESM-2 150M: ~92 min;
     - ESM C 300M: ~140 min;
     - 3Di: ~17 h (Foldseek route) or ~10 h (route in change 11).
   - Restricting to a subset only works for *pairwise* expensive features, e.g. 3Di alignment on the 46,733 rank-≤5 pairs. Apply that rule identically at training and inference.
   - For development ablations, embed only the genes in labelled pairs (a few thousand) and compare 35M / 150M / 650M / ESM C there. PLM features are per-protein and do not depend on the rest of the proteome, so this is exact, provided no proteome-percentile transform is applied to them.

13. **(must) Define proteome-percentile transforms over the CLI's full called-pair set, and compute them the same way at training time.** Do not compute them over the labelled set. Store the search parameters used (sensitivity, caps, MMseqs version) in the model bundle. The CLI must refuse to apply a model whose feature schema or search parameters differ: assert the feature names, their order and a schema hash.

14. **(should) Drop per-gene embedding PCA components from the main plan; keep them only as an explicitly flagged ablation.** With a few hundred robust positives they will mostly memorise families. Grouped CV with the residual homology in change 7 cannot fully exclude this.

15. **(could) Add zero-cost PLM features from the same forward pass.** Examples: mean per-residue entropy, or the native-residue log-probability from the LM head (`AutoModelForMaskedLM`), as sequence-only constraint proxies for each gene. Also the per-pair difference of these values, as a divergence-asymmetry signal. Evaluate them identity-stratified like the other PLM features.

16. **(could) Compute the triangle feature in chunks.**
   - A full sparse A·A on the human graph touches ~5.7·10⁸ wedges, which risks several GB of memory.
   - Instead, multiply row blocks masked by A, or count triangles only on the stricter bitscore-ratio ≥ 0.3 graph, which is also biologically more meaningful.
   - Store the pair table as float32 Parquet in chunks, so that 10⁷-pair proteomes (polyploid crops) do not exhaust 24 GB.

### Models

17. **(must) Fix the torch/LightGBM libomp conflict.** Measured: a segfault in this venv when torch is imported before a multi-threaded LightGBM fit. Choose one of:
   (a) run PLM and 3Di inference in a separate subprocess, which also isolates MPS memory;
   (b) use sklearn `HistGradientBoostingClassifier`. It supports missing values natively, as well as monotone constraints and sample weights, needs no extra OpenMP, and I measured it as working after importing torch.

   I recommend (b) as the shipped model, with LightGBM kept for development only if (a) is implemented. Either way, add a CLI test that runs the whole pipeline in one process.

18. **(should) Constrain the GBDT capacity to the label size.**
   - Settings: about 15 leaves or fewer, `min_samples_leaf` ≥ 30–50, about 20–40 features, early stopping inside grouped CV.
   - Always report the logistic-regression comparator next to the GBDT.
   - Class weighting does not help AUPRC (it ranks the same); keep only the species-balancing weights.
   - The CLI's `sl_score` is not calibrated for a new species' prevalence. Say so in the output header and rely on `sl_rank`.

19. **(should) Handle the family-size mismatch between training and inference.**
   - 87% of output pairs involve genes with more than 100 paralogs, whereas screened pairs have a median family size of 3.
   - Add an output column `in_training_domain`, e.g. local family size and identity within the training range.
   - Report score distributions by family-size bin on the full human output, as a check that giant families are not ranked high by tree extrapolation.

### Packaging, tests, runtime

20. **(must) Pin and record every external artefact.**
   - The MMseqs2 version: the CLI checks `mmseqs version` and warns or refuses on a mismatch, since hit lists and features change between versions. The same applies to Foldseek if it is used.
   - HF model revisions by commit hash, plus sha256.
   - The panel release, plus sha256.
   - The CNN checkpoint sha256.
   - Write all of them into the output header.
   - `slparalog setup` should also rewrite panel FASTA headers to `accession|taxid`. The QfO headers carry GN= and protein names, which are annotation that must never reach the feature code.

21. **(must) Add these tests on top of the planned ones.**
   - Input-order and ID invariance: shuffle the FASTA order and rename every ID → identical scores. This protects against hidden ID or ordering leakage.
   - Thread determinism: MMseqs with `--threads 1` vs `10` → identical pair set and features.
   - A golden-file regression on a ~300-protein yeast subset, with tolerances.
   - A CPU-vs-MPS embedding check (cosine ≥ 0.999) using ESM-2 8M, which keeps CI fast.
   - A schema/feature-hash test of the model bundle.
   - Pure-Python feature unit tests on hand-built hit tables with known answers, separate from a skip-if-no-binary MMseqs integration test.
   - The single-process torch+GBDT test from change 17.
   - A CPU-only run, with no MPS available.

22. **(should) Measure memory correctly.**
   - `/usr/bin/time -l` reports the peak of one process, not the concurrent sum across MMseqs children. It does not include MPS driver allocations, which reached 1.5–7.6 GB in the benchmarks.
   - Use a 1-s sampler of the whole process tree (`psutil`, RSS plus `phys_footprint`) and log `torch.mps.driver_allocated_memory()`.
   - Also record `pmset -g`, the power source, and thermal state.
   - Re-benchmark with Low Power Mode off before freezing the "≤ 2 h, ≤ 8 GB" target. All current numbers were taken with Low Power Mode on and 4 torch threads.

23. **(should) Ship explicit tiers (`--tier fast|standard|full`), with one model trained per tier on exactly that tier's features.**
   - Do not ship a single model that receives NaN for absent tiers.
   - Suggested defaults:
     - `fast`: Tiers A+B, about 25 min for human;
     - `standard`: + ESM-2 150M or 35M, chosen by ablation;
     - `full`: + 3Di, used only if it adds AUPRC conditional on identity.
   - Report runtime and AUPRC for every tier.

24. **(could) Overlap the panel MMseqs search (CPU) with PLM embedding (MPS).** This saves about 10–20 min of wall time. Cap MMseqs at about 6 threads while the PLM runs, and check the combined peak memory.

## 3. Risks the design does not address

- **Panel-density confound.** Conservation features learned on human, with a primate- and mammal-dense panel, may transfer poorly to insects and plants, whose panel coverage is 3–5 arthropods and 5 plants or algae. Changes 8–9 reduce this, but they do not remove it.
  - Add a robustness experiment: re-run the human and yeast features with a "thinned" panel (e.g. all vertebrates removed for human) and report the AUPRC change.
  - Consider adding a few more arthropod and plant reference proteomes. The cost is about linear, ~11 s per species.
- **Sensitivity/cap drift.** Family-size, rank and percentile features all depend on `-s`, `--max-seqs` and the MMseqs version. Any later change to the search silently changes the meaning of features for a frozen model (change 13 guards this).
- **Fragmented or partial gene models in new proteomes.** These shrink coverage and global identity and inflate one-directional hits. The planned evaluation proteomes are all high-quality UniProt reference sets, so there is no test of degradation on, for example, BUSCO-incomplete annotations.
  - Could: re-evaluate the yeast and human models after randomly truncating 20% of proteins.
- **The 3Di tier may cost about 10–17 h for little gain.** Dennler found structural similarity highly correlated with identity. Run the identity-stratified ablation on labelled-pair genes first (about 7k proteins, 4–6 h, one-off) before investing in full-proteome 3Di.
- **Output size.** About 0.7 M rows for human and possibly more than 10 M for polyploid crops. Write gzip output and a top-N companion file. Keep "every pair" as the default, per PROMPT.
- **ESM C licence.** The EvolutionaryScale terms differ between the 300M and 600M models and between releases. Re-check before releasing any tier that depends on ESM C.
- **GPL-3 binaries.** Foldseek and DIAMOND are GPL-3; MMseqs2 is MIT. Calling them as subprocesses and fetching them in `setup` is fine. Do not vendor them inside the Python wheel.

## 4. Fine as it is

- MMseqs2 for both searches. DIAMOND is not needed: it reached a 19 GB footprint against the panel.
- `-s 7.5` with E ≤ 1e-3 for the self-search, and no identity or coverage floor in calling. Keeping identity and coverage as features rather than filters is right.
- Mean-pooled, chunked PLM embeddings in fp16 on MPS. fp16 was verified equivalent to fp32.
- No ESMFold or AF2, which is correctly judged infeasible.
- Identity-stratified evaluation of the PLM and 3Di features.
- Training features produced by the CLI's own feature code.
- End-to-end AUPRC in which missed pairs receive the minimum score.
- The overall shape of the CLI: `slparalog setup` / `slparalog predict`, the src layout, and the console script.
- Memory is not the binding constraint for human-size proteomes; time is.
- The species-agnostic feature set with no species indicator, and refitting the released model only after the test evaluation is recorded.

Sources: [Foldseek README (ProstT5 databases have no Cα, so no TM-score/LDDT)](https://github.com/steineggerlab/foldseek); [ProstT5 model card](https://huggingface.co/Rostlab/ProstT5); [D-SCRIPT predict_3di.py (ProstT5 CNN head checkpoint URL)](https://huggingface.co/spaces/samsl/D-SCRIPT/blob/main/predict_3di.py).
