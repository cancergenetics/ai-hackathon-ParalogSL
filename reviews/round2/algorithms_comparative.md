# Round 2 review: algorithms/compute and comparative genomics (combined)

Reviewer role: combined algorithms/compute and comparative-genomics reviewer.
Inputs read: PROMPT.md, DESIGN.md v2 (all sections; §2.2, §3, §4, §7 and §8 in detail), `reviews/round1/{algorithms_compute, comparative_genomics, PI_response}.md`, `data/raw/qfo/{README.md, species_list.tsv}`, and the existing benchmark outputs in `data/interim/bench_mmseqs/`.

## New measurements in this round

These were light checks. Nothing was downloaded into the project, and no project file other than this review was touched. Scripts and outputs are in `<scratchpad>`. A project script must reproduce any number before it is relied on.

The machine was in Low Power Mode with a load average of 5–22, because other agents were running, so every timing below is pessimistic.

| # | Check | Result |
|---|---|---|
| a | D_s for a human query, computed from the existing human-vs-QfO search (`ds.py`, `ds2.py`). Median RBH identity, both coverages ≥ 0.5. | PANTR 0.993, GORGO 0.990, MACMU 0.973, mouse 0.872, chicken 0.743, zebrafish 0.631, amphioxus 0.502, Drosophila 0.449. **All 27 non-bilaterian and non-metazoan species fall in 0.35–0.43**, in a non-phylogenetic order (Monosiga 0.425 > Arabidopsis 0.399 > yeast 0.387; C. elegans 0.414 ranks with fungi and plants). The number of RBH pairs falls from 17k (mammals) to 0.8–3k (protists), so each species' median is taken over a different, increasingly conserved gene subset. |
| b | Same, restricted to a common core gene set: human genes with a coverage ≥ 0.5 RBH in ≥ 80% of species (1,082 genes). | Range widens to 0.38–0.996, and the order becomes mostly sensible: mammals 0.94–0.96, other vertebrates 0.77–0.86, invertebrates 0.56–0.65, Monosiga/C. elegans 0.50, plants 0.47, Dikarya 0.42–0.44. Long-branch effects remain: C. elegans sits next to Monosiga, and fungi sit beyond plants. |
| c | Panel search on 400 random human queries (`bench/`, `rec.py`, `rec2.py`). Reference = one search per species (`--max-seqs 50`, `-s 5.7`). Measured: share of (query, species) cells with a hit at ≥ 50 bits that the concatenated search recovers. | Global `--max-seqs 100`: **77.8%** recovered (190/400 queries saturated). Global `--max-seqs 1000`: **94.9%** recovered (45/400 saturated). Even queries **below** the cap lose 1.1–1.3% of cells, and 179–205 cells change best hit. Losses concentrate in **distant** species (Plasmodium, Dictyostelium, Paramecium, S. pombe, Pezizomycotina): exactly where age and conservation features live. |
| d | The same per-species reference vs global search at E ≤ 1e-3. | Global search drops weak cells that the per-species search keeps (median 46–48 bits), because E-values are ~50× smaller in a 20k-protein species database. So "species has a hit" depends on how the database is partitioned unless presence is defined by bitscore. |
| e | Cost: 400 vs 1,600 human queries, 8 threads. | Per-species (51 species, `--max-seqs 50`): 150 s vs 300 s, i.e. **0.125 s per query** plus ~100 s fixed, giving **≈45 min for human** under these conditions. Global `--max-seqs 1000`: 0.05 s per query, ≈18 min. The earlier full-human `--max-seqs 100` run took 560 s. In the per-species run, 9.5% of (query, species) cells reached the 50-hit cap. |
| f | MMseqs2 `pident`/`nident` without vs with backtrace (`-a`): 400 human queries vs the human proteome, 26,656 pairs. | **Without `-a`, `nident` = 0 and `pident` is a score-based estimate.** It differs from the true identity by > 5 points for **29%** of pairs and by > 10 points for 4% (median +1.5). With `-a`, `pident` = `nident`/`alnlen` exactly. |
| g | UniProt REST metadata, release 2026_03 (`ref_*.json`, `an.py`, `genus.py`, `select.py`). | Works: `https://rest.uniprot.org/proteomes/stream?query=(proteome_type:REFERENCE) AND (taxonomy_id:T)&format=json`. Each record carries BUSCO complete/single/duplicated/fragmented/missing, `lineageDb`, the CPD status, `proteinCount` and the full `taxonLineage` (genus/family/order). Note: `proteome_type:1` returns 0 hits; the value must be the string `REFERENCE`. Canonical FASTA files are at `.../current_release/knowledgebase/reference_proteomes/Eukaryota/UPxxxx/UPxxxx_taxid.fasta.gz`. |
| h | What "highest BUSCO per clade" actually picks. | **Saccharomycetaceae**: the top 9 by BUSCO are all pre-WGD, and 8 of them are *Lachancea*. **Hymenoptera**: *Nasonia vitripennis* reports C = 96% but **D = 53%**. **Lycophytes**: *Selaginella* D = 65%. **Gymnosperms**: none usable (best S = 47%). **Monocots**: under C ≥ 90 / D ≤ 10 / F ≤ 5 the only pass is *Oryza sativa indica*, i.e. rice again. **Nematoda**: the top picks are *Caenorhabditis* congeners of the C. elegans test species. BUSCO lineage databases also differ within a clade (eudicots_odb10, fabales_odb10, brassicales_odb10, …), so scores are only roughly comparable. |
| i | My proposed rule (change 8) run on that metadata. | Adds **29 species / 576,744 proteins**, giving **80 species and 1.55 M proteins** (inside v2's targets of 80–90 species and ≤ 2.5 M proteins). The species list is in change 8. |

## 1. Overall verdict

v2 fixes most round-1 must-items of both reviews in substance: symmetric features, explicit self-scores, no raw E-values, mutual-coverage families, Leiden splits, no TM-score from 3Di, torch/GBDT isolation, robustness test, WGD/SSD strata, and the move to divergence-calibrated, mask-augmented panel features. That last change is the right idea.

However, the new elements are specified too loosely to implement or pre-register, and three of them break on real data:
- D_s, as a median over each species' own RBH set, collapses all deep species into a 0.35–0.43 band.
- Neither panel-search option in v2 meets a per-species guarantee as worded. A capped concatenated search loses 5% of (query, species) cells even at `--max-seqs 1000`.
- "Highest BUSCO per clade" picks redundant or wrong species: 8 *Lachancea*, a 53%-duplicated wasp, rice again for monocots.

Separately, every identity value from MMseqs2 must come from backtraced alignments. Below are concrete specifications that make each element implementable. With them the design is ready for M1/M2 from my side.

## 2. Concrete changes

### 2.1 Status of round-1 must-items

| Round-1 item | v2 status |
|---|---|
| Alg 1: self-search cap + re-search on saturation | **Partial.** Saturation must be detected from **prefilter** counts, not aligned-hit counts (change 4). |
| Alg 2: symmetric features, both orientations realigned | Fixed, provided the realignment uses `-a` (change 1). |
| Alg 3 / Comp S5: no raw E-values | Fixed for features. **Not fixed** for panel presence and for the self-search calling threshold (change 3). |
| Alg 4: explicit self-scores | Fixed. |
| Alg 7 / Comp M7: single-linkage families and splits | Fixed. |
| Alg 8 / Comp M1: panel truncation | **Not yet.** The "< 1% of queries saturate" criterion is the wrong criterion (measurement c). See change 2. |
| Alg 9 / Comp M2, M3: divergence-relative features, query placement, thinning | **Partial.** Right idea, but D_s, bins, depth and masks are undefined or ill-defined (changes 5–7). The self-exclusion text is wrong (change 9). |
| Alg 10: no TM-score/LDDT from 3Di | Fixed. |
| Alg 12: compute plan for PLM restriction | Fixed. |
| Alg 13: percentiles over the called-pair set, schema hash | Fixed. Percentiles must also be recomputed per mask (change 7). |
| Alg 17: libomp conflict | Fixed. Do not use HGB's internal early stopping with mask replicates (change 7d). |
| Alg 20: pin artefacts | Mostly fixed. The gap-filled panel must be snapshotted (change 8f). |
| Alg 21: extra tests | Mostly fixed (change 12). |
| Comp M4: gap filling | **Rule not implementable as written** (measurement h). See change 8. |
| Comp M5: robustness test and near-identical flag | Fixed. The flag needs exact identity (change 1). |
| Comp M6: WGD/SSD stratification | Fixed. |

### 2.2 Paralog calling and search mechanics

**1. (must) Require backtraced alignments (`-a`; `--alignment-mode 3` in `align`) wherever identity, `nident` or alignment coordinates are used. Store `nident`, `alnlen`, `qstart`, `qend`, `tstart`, `tend`, `qlen`, `tlen` and `bits`.**
- Applies to the self-search realignment (§3), the panel search (§4, items 4–6), the B1 identity variants, and the `near_identical` flag.
- Why: without backtrace, MMseqs2 reports an estimated `pident` (> 5 points wrong for 29% of pairs) and `nident` = 0 (measurement f). Ensembl-style identity, global identity, the ≥ 99% flag, D_s and the asymmetry features would all be silently wrong.
- Add a unit test that `pident` equals `nident`/`alnlen` on the toy proteome.

**2. (must) Specify the panel search as one search per panel species. Drop the "single search with a large cap" option.**
- Exact spec: per-species MMseqs2 databases are built once in `slparalog setup`. Each run uses `search -s 5.7 -a --max-seqs 50 -e 10` (permissive; presence is decided by bitscore, see change 3) and keeps all hits with bits ≥ 40.
- Why: measurement c. Even at `--max-seqs 1000`, the concatenated search loses 5.1% of (query, species) cells at ≥ 50 bits, loses cells even for unsaturated queries, and changes the best hit in ~1.6% of cells. The losses are biased towards distant species. "< 1% of queries saturate" cannot detect this.
- Acceptance benchmark (M2, pre-registered):
  - 500 human and 500 yeast queries; reference = per-species `-s 7.5 --max-seqs 300`.
  - Candidate settings: `--max-seqs` ∈ {20, 50} × `-s` ∈ {5.7, 7.5}.
  - Adopt the cheapest setting with ≥ 99% concordance of per-(query, species) best hit and ≥ 99% presence recall at ≥ 50 bits.
  - Commit the result table.

**3. (must) Define every presence, absence and calling threshold in bitscore, not E-value.**
- Panel: "A has a hit in s" ⇔ bits ≥ 50. This is about E = 1e-3 for a 500-residue query against the concatenated 1.5 M-protein panel, so it matches the current intent without depending on how the database is partitioned (measurement d).
- Self-search calling: keep the search at `-e 1e-3`, or better at `-e 10` for headroom, but call a pair by E′ ≤ 1e-3 with E′ = qlen · N_ref · 2^(−bits). N_ref = 1.1·10⁷ residues (≈ human proteome) is fixed in `config/`. This makes the call proteome-size invariant: yeast is 6k proteins, a crop can exceed 40k (Comp S5, not done in v2).
- Report caller recall of labelled pairs under both definitions once.

**4. (should) Detect saturation from prefilter counts.**
- Self-search: run with `--remove-tmp-files 0` and read per-query entry counts from the prefilter database (`pref_*` index). Re-search queries with ≥ 0.95 · `max-seqs` prefilter candidates using a 5× cap (`createsubdb`).
- Log the counts. Apply the same check per species in the panel search.
- Aligned-hit counts below the cap do not prove that the prefilter was not truncated (round-1 Alg 1).

### 2.3 Divergence calibration, bins, duplication depth, masking

**5. (must) Define D_s on a common core gene set, with an additive distance. Fix the bin edges in advance.**

Exact spec:
1. From the per-species forward table, take **pseudo-RBHs**: q → x is an RBH if x is q's top-bitscore hit in s, q is the top-bitscore query among all queries hitting x, and both coverages are ≥ 0.5. This is valid once change 2 holds; no reverse search is needed.
2. **Core set** G* = query proteins with such an RBH in ≥ 80% of the retained panel species (1,082 genes for human). Fallback if |G*| < 200: lower the threshold in steps of 10% and log it.
3. id_s = median identity (`nident`/`alnlen`) over G* ∩ RBH(s). d_s = Kimura protein distance −ln(1 − p − 0.2p²) with p = 1 − id_s. Cap at d = 3 if the argument is ≤ 0.
4. Write d_s, |G* ∩ RBH(s)| and the RBH fraction per species to the run log.

Why: the per-species-own-set median compresses all deep species into a 0.35–0.43 band with no usable order (measurements a–b). The core set gives the same genes for every species and a monotone, roughly additive scale.

Validation (M2, no labels):
- Spearman correlation of d_s with TimeTree divergence for the human and yeast queries. Report it.
- Expect long-branch outliers (C. elegans, Saccharomycotina). The bins below are coarse at the deep end for this reason.

**Bin edges** in d (from human core-set values; to be frozen before yeast D_s is computed):

| Bin | d range | Human members |
|---|---|---|
| B1 | [0, 0.03) | macaque |
| B2 | [0.03, 0.12) | eutherian mammals, opossum |
| B3 | [0.12, 0.35) | birds, amphibians, fish |
| B4 | [0.35, 0.70) | non-vertebrate bilaterians, cnidarians |
| B5 | [0.70, 0.95) | choanoflagellates, nematodes, plants, Dictyostelium, Chlamydomonas, oomycetes |
| B6 | [0.95, ∞) | Dikarya, diatom, kinetoplastids, apicomplexans, diplomonads, parabasalids, ciliates |

- A bin with no retained species is **NaN** ("not observable"). This is distinct from 0 ("species present, no hit"). Feature code must keep the two apart; add a test.

**6. (must) Replace the item 5 definition of duplication depth, which is garbled in v2 ("distinct best hits closer to each than to each other"), with a loss-aware, InParanoid-style classification on the d axis.**

For pair (A,B) and retained species s, let a = A's best hit in s and b = B's best hit in s. Compare bitscores from the query's side (bits(A,·) are comparable because A is the query in both searches; the same gap and composition-bias settings must be used in both searches, plus a test).

Classify s as:
- **pre** (the duplication predates the split from s; both copies were inherited): a ≠ b, bits(A,a) > bits(A,B) and bits(B,b) > bits(B,A).
- **post** (the duplication is younger than the split; A and B are in-paralogs relative to s): bits(A,B) > bits(A,a) and bits(B,A) > bits(B,b). This holds whether or not a = b.
- **undetermined**: everything else, including gene loss in s and strong rate asymmetry.
- **absent**: neither A nor B has a hit at ≥ 50 bits.

Depth:
- lower bound L = max d_s over **pre** species;
- upper bound U = min d_s over **post** species;
- depth = (L + U)/2 when L < U.
- When L ≥ U (conflict), set depth = argmin_t of the bin-weighted count of discordant species, i.e. pre species with d_s > t plus post species with d_s < t, with each bin's species weighted 1/n_bin. Report the discordant fraction as its own feature.

Censoring:
- no pre species → depth ∈ [0, U]; set depth = U/2 and the flag `younger_than_nearest` = 1;
- no post species → depth = L and the flag `older_than_panel` = 1.

Derived features:
- depth, U − L, the discordant fraction, the two flags, and per-bin fractions of pre/post/undetermined species (NaN for empty bins).
- **(should) Relative rate** = d(A,B) / (2·depth), using the Kimura distance on the pair's backtraced identity. This is the paralog divergence relative to the genome-average clock, an identity-orthogonal proxy for asymmetric or accelerated evolution.
- Test on a toy panel: a planted pre-split duplication, a post-split duplication, and a duplication with one copy lost in s.

**7. (must) Define the training masks as absolute d ceilings, not "k closest bins", and fix how they interact with weights, percentiles and early stopping.**

a. **Mask ladder** c ∈ {none, 0.03, 0.12, 0.35, 0.70}. Mask c removes every panel species with d_s < c. Masks are post-processing on one hit table (no re-search), so they cost seconds.
   - "k closest bins" means different things for human (bins B1–B3 occupied) and yeast (B1–B2 empty), and is ambiguous about empty bins.
   - With absolute ceilings, a mask that removes nothing for a species yields identical rows. Deduplicate them and spread that species' weight over its distinct configurations, so species balance is preserved.
   - The ladder emulates, in order: a human with no mammals, then an "insect-like" isolation (d_nearest ≈ 0.35–0.5), then a "plant-from-a-metazoan" isolation.

b. **Per-mask recomputation.** D_s is computed once on the full panel. Everything downstream is recomputed per mask: bins, depth, conservation, species rank, and proteome percentiles of Tier-B features. Do this through a hidden CLI option (`--panel-dmin c`), so that training still goes through the shipped code path.

c. **Primary evaluation conditions** (pre-register):
   - test on the natural panel (mask none, after self-exclusion);
   - report the full ladder as secondary;
   - report the unseen species (Arabidopsis, C. elegans, S. pombe) on the natural panel and on the ladder, with the realised d_nearest of each.

d. **GBDT with replicated rows.** sklearn HistGradientBoosting's built-in early stopping splits **rows** at random, so mask copies of one pair fall in both the train and validation parts. Set `early_stopping=False` and choose `max_iter` by grouped CV. Also express `min_samples_leaf` in pairs, i.e. multiply by the number of distinct mask configurations, or use weights. Otherwise a 30-row leaf can be six pairs.

e. **Species-identity leakage.** Empty-bin patterns, flags and percentiles let the model recognise human vs yeast. The two species also differ in label prevalence. Two fixes:
   - reweight so that each species has the same weighted prevalence, as well as the same total weight;
   - add a control: a classifier predicting species from the Tier-B features on the natural panel, with its AUROC reported, and Tier-B importance reported conditional on that.

**8. (must) Replace "highest BUSCO per clade" with an explicit, phylogenetically spread rule. The tested version (`select.py`) is below; freeze it in `config/panel.yaml`.**

a. **Candidates.** UniProt `proteome_type:REFERENCE` under each clade taxid. Exclude organisms, and the genus/family units, already in the core QfO panel. Require `proteinCount` ≤ 50,000 (no polyploid or isoform-inflated sets).

b. **Quality.** BUSCO D ≤ 10% (this removes *Nasonia* at D = 53% and *Selaginella* at D = 65%). Ladder (S ≥ 85, F ≤ 8) → (S ≥ 80, F ≤ 8) → (S ≥ 75, F ≤ 10), applied per clade only if its quota is unmet, and logged.

c. **Spread.** At most one proteome per unit:
   - genus for Saccharomycetaceae;
   - order for eudicots and basal plants;
   - family for monocots, insects and nematodes.
   Rank by S, then F (lower first), then UPID. Pre- and post-WGD Saccharomycetaceae need **explicit genus lists** (BUSCO ranking alone returns only pre-WGD genera):
   - post-WGD: Saccharomyces, Nakaseomyces, Kazachstania, Naumovozyma, Tetrapisispora, Vanderwaltozyma, Huiozyma, Maudiozyma, Monosporozyma, Arxiozyma, Henningerozyma;
   - pre-WGD: Kluyveromyces, Lachancea, Eremothecium, Zygosaccharomyces, Torulaspora, Zygotorulaspora.
   This is panel composition, not a feature.

d. **Quotas, and what they select today** (release 2026_03):

| Clade (quota) | Selected |
|---|---|
| Post-WGD Saccharomycetaceae (3) | *N. glabratus*, *Huiozyma naganishii*, *Kazachstania africana* |
| Pre-WGD Saccharomycetaceae (3) | *Torulaspora delbrueckii*, *Lachancea lanzarotensis*, *Eremothecium gossypii* |
| Eudicots, orders other than Brassicales (6) | *Dipteronia*, *Nicotiana attenuata*, *Prunus*, *Cephalotus*, *Daucus*, *Lithocarpus* |
| Monocots, not Poaceae (3) | *Dendrobium*, *Zostera*, *Spirodela* (via the S ≥ 75 step) |
| Basal angiosperm + liverwort (2) | *Amborella*, *Riccia* |
| Hymenoptera (2) | *Atta cephalotes*, *Vespula maculifrons* |
| Coleoptera, not Tenebrionidae (2) | *Phyllotreta*, *Aquatica* |
| Lepidoptera (2) | *Mythimna loreyi*, *Dendrolimus* |
| Hemiptera (2) | *Aphis craccivora*, *Laodelphax* |
| Polyneoptera/Palaeoptera (1) | *Zootermopsis* |
| Nematoda, not Rhabditidae (3) | *Cylicocyclus*, *Brugia pahangi*, *Thelazia* |

   Total: +29 species, +577k proteins, giving **80 species and 1.55 M proteins**.
   - No gymnosperm or lycophyte passes the quality filter. State this as a known gap.
   - Optional: *Marchantia* (S = 78.9) may replace *Riccia* if the PI prefers the model liverwort. Decide that before any result exists.

e. **Test-species congeners.** Excluding units already in the core panel keeps *Arabidopsis lyrata*/*Capsella* and *Caenorhabditis* congeners out. The Arabidopsis and C. elegans tests then resemble the "no close relative" use case instead of being made easy by a sister species. State this choice explicitly, and note that change 7c reports the ladder anyway.
   - Alternative, if the PI wants yeast to also see a close-relative regime like human/macaque: allow one non-*cerevisiae Saccharomyces* (*S. mikatae*, S = 97.6). I recommend against it, because masks already cover that regime through human.

f. **Reproducibility.** The REST API serves only the current release.
   - Save the JSON metadata responses together with the `x-uniprot-release` header (2026_03 now). Take checksums. Have `build_panel.py` select from the saved snapshot.
   - Download the canonical `UPxxxx_taxid.fasta.gz` (not `_additional`) from the same release.
   - Record that the core panel is QfO 2026_02 while gap-fill and training proteomes are 2026_03. This is harmless for self-exclusion, which works on identity, but DATA.md must say it.

**9. (must) Correct the self-exclusion statement.**
- At median RBH identity ≥ 0.98, the rule removes **chimpanzee (0.993) and gorilla (0.990)**, not only "the query species itself or a strain of it". Macaque (0.973) is kept.
- No sequence-only rule can separate a conspecific strain from a sister species about 6 My apart. Keep the 0.98 threshold, but state it as "self, conspecific strains and sister species", log which species were removed, and add the human example to DATA.md/REPORT.md.

### 2.4 Features: remaining comparative items

**10. (should) Divergence asymmetry (§4 item 6): use the nearest informative outgroup, the shared region, and distances.**
- Outgroup O = the RBH partner of A or B in the species with the smallest d_s classified **post** in change 6. That is the closest species that predates the duplication.
- Compute id(A,O) and id(B,O) only over the columns of O covered by both alignments (from the backtraced coordinates of change 1), convert to Kimura distances, and use |d_AO − d_BO| / (d_AO + d_BO).
- Do it for the 3 nearest post species and take the median.
- Cost: ≤ 3 triplets per pair from rows already in the per-species table, so no new search.
- Raw identity saturates below about 30%, where human 2R pairs and yeast pairs sit against their outgroups (Comp S4, not in v2).

**11. (should) Add the two cheapest WGD-like retention proxies and explicit missingness indicators.**
- `n_post_single`: the number of **post** species where a = b, and a has no second hit within 0.9 × bits(a's best query). This is the 1:2 orthology pattern.
- `family_duplicability`: the median over species of the number of hits ≥ 0.5 × self-bits of A.
- Indicator columns for "no Tier-B evidence", and AUPRC on pairs with complete Tier-B features (Comp S6, S8).
- `fragment_suspect` in v2 says "median homolog length". Define it as the median length of A's best hits in retained panel species. Self-paralogs share the same assembly artefacts, so they are not a fair reference.

### 2.5 Compute, tests

**12. (must) Update the §7 runtime budget, and add the missing tests.**

Runtime:
- Measured per-species panel search (change 2 setting, 51 species) is ≈ 0.125 s per query plus ~100 s, giving **≈ 45 min for human**. With the 80-species, 1.55 M-protein panel this becomes **≈ 70 min**.
- These figures were taken under Low Power Mode and contention, so expect about half that on an idle machine. v2's "≈ 15–30 min" is not supported yet.
- The one-off feature build for human, yeast, Arabidopsis, S. pombe and C. elegans (≈ 79k queries) is ≈ 3–5 h. That is acceptable.
- Re-measure with Low Power Mode off before freezing the ≤ 2 h default-tier target.
- If the budget is tight, the change-2 benchmark may adopt `--max-seqs 20` or `-s 5.7` for species with d < 0.35, where sensitivity is not the issue. Only do this if it passes the 99% concordance criterion.

Tests to add:
- the `-a`/`nident` test (change 1);
- NaN-vs-0 bin semantics (change 5);
- the toy duplication-depth triplets (change 6);
- mask-copy grouping in CV and no row-random early stopping (change 7d);
- order invariance of core-set selection: G* depends only on hits, never on FASTA order;
- a CPU-only single-process run including torch + GBDT (round-1 Alg 21, not listed in v2 §6.5/§8).

**13. (could) Run label-free distribution-shift diagnostics** on 2–3 unlabelled gap-fill proteomes run as queries (e.g. *Atta*, *Nicotiana*, *Aphis*). Report a per-feature PSI against human and yeast under the natural panel. This uses no labels and is the only direct evidence about insects and non-Arabidopsis plants (round-1 Comp S9).

## 3. Risks the design does not address

1. **The deep end of D is rate-confounded.** C. elegans and Saccharomycotina are long-branch lineages. Any query or panel species with an atypical rate is misplaced on d, and depth estimates for old duplications inherit that. Coarse deep bins and the discordance feature limit the damage, but do not remove it.
2. **Two training species give two natural-panel configurations.** Masks create intermediate densities only for human (B1–B3) and, weakly, yeast. The "plant query" geometry, with dense eudicot neighbours at d ≈ 0.1–0.3 and no metazoan bins, is never seen at training. Only the Arabidopsis test and the change-13 diagnostics speak to it.
3. **BUSCO is not comparable across lineage databases**, and UniProt's BUSCO is not a contamination check. A gap-filled proteome with bacterial or TE contamination would pass. Check the selected 29 once for a long tail of proteins that have no hits outside their own species (orphans > 30% would be a red flag).
4. **Per-species search costs scale with panel size × query size.** A 100k-protein polyploid query would take ~3–5 h for Tier B alone. Document the runtime per 10k proteins and offer `--threads`. Do not silently subsample.
5. **The bitscore comparability in change 6** assumes identical scoring settings in the self-search and the panel search. A future MMseqs2 upgrade, or a different `--comp-bias-corr` default, would shift pre/post classifications. Pin the settings and include them in the schema hash.

## 4. Things that are fine as they are

- The decision to place the query by sequence divergence (D_s) instead of a species tree, and to emulate sparse panels by masking during training. These are the right answers to round-1 Alg 9 and Comp M2/M3.
- Self-search at `-s 7.5`, `--max-seqs ≥ 1000`, both orientations realigned, explicit self-scores, symmetric transforms, no raw E-values as features.
- Mutual-coverage local family size, Leiden-partitioned splits, and per-test-pair leakage reporting.
- Tier C/D gating, 3Di-alignment-only Tier D, torch/GBDT process isolation, and the feature registry and schema hash.
- WGD/SSD and ribosomal stratification as evaluation-only, the robustness test with degraded proteomes, and the `near_identical`/`fragment_suspect` flags, given change 1.
- A target panel size of 80–90 species and ≤ 2.5 M proteins: the explicit rule lands at 80 / 1.55 M. Disk use is < 2 GB.
