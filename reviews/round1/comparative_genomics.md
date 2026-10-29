# Round 1 review: comparative and evolutionary genomics

Reviewer role: comparative and evolutionary genomicist (orthology inference, duplication dating, WGDs, phylostratigraphy, proteome quality).
Inputs read: PROMPT.md, DESIGN.md v1, LITERATURE.md, docs/lit/methods.md (sections 1, 2, 5, 6), docs/lit/compute_feasibility.md, data/raw/qfo/{README.md, species_list.tsv}.
I also ran a few quick checks on benchmark outputs that already existed (scratch scripts only; no project files changed):
- `data/interim/bench_mmseqs/hq_s5.7_ms100/out.m8`: human vs the QfO panel.
- `.../hh_s7.5/out.m8`: human vs human.
- One small re-search: 150 human proteins vs the panel at `--max-seqs 3000`, which took 79 s and 4.2 GB.

The numbers quoted below come from those checks.

## 1. Overall verdict

The feature families are the right ones: similarity, family context, conservation and age against a panel, and divergence asymmetry. Not trying to call WGDs from synteny is also correct. But the Tier B implementation details would make the duplication-age and conservation features measure *panel density around the query* rather than biology. Three things cause this:
- The panel search truncates hit lists.
- Phylostratum assignment needs the query's taxonomic position, which the tool never has.
- The nearest panel relative ranges from about 29 My (human, after self-exclusion) to more than 300 My (a new hymenopteran or lepidopteran).

A model trained on human and yeast would then learn features whose meaning shifts in exactly the target cases: plants and insects. Each of these problems has a cheap fix inside the current scope. A robustness test for fragmented or incomplete gene models is also missing and should be pre-registered.

## 2. Concrete changes

### Must

**M1. Fix hit-list truncation in the panel search; exclude the query species before the search, not after.**
The benchmarked panel search (`-s 5.7 --max-seqs 100`, human vs all 51 species) saturates for 49% of human queries: 10,132 of 20,652 have 100 hits, and self and in-paralog hits use up slots.

I re-searched 150 randomly chosen saturated queries at `--max-seqs 3000`:

| Measure | `--max-seqs 100` | `--max-seqs 3000` |
|---|---|---|
| Median panel species with a hit | 22 | 34 |
| Queries with a yeast hit | 20/150 | 65/150 |
| Queries with an Arabidopsis hit | 35/150 | 69/150 |
| Queries with a Drosophila hit | 78/150 | 119/150 |
| Queries losing ≥1 species | — | 142/150 |
| Queries losing ≥5 species | — | 106/150 |

The genes that have paralogs are exactly the ones that saturate. As configured, conservation breadth, phylostratum and "best hit shared" would systematically make multi-gene families look younger and less conserved. This bias would also be stronger for human, with its 15 vertebrate relatives, than for yeast, so it differs between the two training species.

Fixes:
- (a) Build the target DB without the self or near-identical species. Detect them in a quick pre-pass on, say, 500 random query proteins.
- (b) Guarantee per-species coverage: search each species or clade DB separately with a modest `--max-seqs`, or use `--max-seqs` ≥1000.
- (c) Benchmark per-species recall against an uncapped reference on a sample, as I did above, and pick the cheapest setting with ≥99% species-level recall.
- (d) Write the saturated-query fraction to the run log.

**M2. Specify how the query is placed on the panel tree.**
Section 4.4 ("deepest panel clade with a hit, using a fixed species tree") and any age expressed in My need to know where the query sits on the tree. The input is an unannotated FASTA with no taxon.

Pre-register a sequence-only placement method: for example, the nearest panel species by median identity of reciprocal best hits (or by RBH count), with the query attached as its sister. Validate it by leave-one-panel-species-out on about 10 panel species across clades, and log the placement in every run. Without this, phylostratum features cannot be computed for a new genome at all.

**M3. Make age and conservation features calibrated by distance, and test them under panel thinning.**
The distance from the query to its nearest retained panel relative varies a lot. The self-exclusion rule, as measured, removes chimpanzee (66.9% of human proteins pass ≥98% identity with ≥90% coverage) and gorilla (57.0%) but keeps macaque (31.7%).

| Query | Nearest retained panel relative | Approx. divergence |
|---|---|---|
| Human | Macaque | ~29 My |
| Yeast | *C. albicans* / *Yarrowia* | ≥200 My |
| Arabidopsis | Rice / maize | ~150 My |
| *C. elegans* | No other nematode | >550 My |
| New hymenopteran or lepidopteran | *Tribolium* / *Drosophila* | ~300–350 My |

Raw counts ("n species where best hits coincide") therefore change meaning between queries. Two changes:
- (a) Express duplication age on a query-independent scale. Use the divergence time (TimeTree) of the deepest panel node that the duplication predates, or a relative "speciation ladder": d(A,B) compared with d(A, ortholog in S) for panel species at known distances.
- (b) Pre-register a panel-thinning protocol on development data. Recompute Tier B with every panel species within T My of the query removed, for T ∈ {0 (self only), ~100, ~300, ~450} (human: drop primates, then mammals, then vertebrates; yeast: drop Saccharomycotina, then Ascomycota). Report AUPRC per condition, and include the thinned conditions in model selection (mean development AUPRC over conditions). Training-time augmentation with thinned panels is the natural follow-up if AUPRC falls as T grows. The plant and insect use cases are the large-T cases, and the current design never evaluates them.

**M4. Fill pre-registered taxonomic gaps in the panel before any modelling; pin it in DATA.md.**
The 51 QfO eukaryotes are vertebrate-heavy: 15 vertebrates, 2 other chordates, 12 fungi and 8 protists.
- Plants: only 4 land plants (Arabidopsis, rice, maize, *Physcomitrium*) plus *Chlamydomonas*. That is one eudicot, and no asterid, non-Brassicaceae rosid, basal angiosperm, gymnosperm, lycophyte or liverwort.
- Insects: only 3, two of them Diptera (*Drosophila*, *Anopheles*) plus *Tribolium*. There is no Hymenoptera, Lepidoptera, Hemiptera or hemimetabolous insect (*Daphnia* and *Ixodes* are non-insect arthropods).
- Yeast: no other Saccharomycetaceae. No pre-WGD outgroup (*K. lactis*, *Eremothecium*, *Lachancea*) and no post-WGD relative (*Nakaseomyces glabratus*, *Naumovozyma*). Yeast ohnologs therefore cannot be told apart from younger yeast-specific SSDs by any panel feature, and their divergence asymmetry can only be measured against distant outgroups (Byrne & Wolfe 2007 used pre-WGD species).

Add about 15–20 UniProt reference proteomes (sequences only, about 300 MB, +35–40% search time). Choose them by a written rule, for example "≥1 proteome per major lineage with a gap of >100 My on the panel tree for eudicots, monocots, land-plant grades, insect orders and Saccharomycetaceae". Do not choose them by test-species performance. Candidates:
- Plants: *Vitis*, *Solanum lycopersicum*, *Populus* or *Medicago*, *Amborella*, *Selaginella*, *Marchantia*, *Klebsormidium*.
- Insects: *Apis*, *Bombyx*, *Acyrthosiphon*, *Pediculus* or *Zootermopsis*.
- Yeasts: *K. lactis*, *E. gossypii*, *L. thermotolerans*, *N. glabratus*, *N. castellii*.

The data agent should check BUSCO completeness on each and reject fragmentary assemblies.

**M5. Pre-register a gene-model robustness test, and flag near-identical pairs in the output.**
Proteomes from newly sequenced plant and insect genomes contain several defects:
- truncated or fragmented models;
- split genes and fused (chimeric) genes;
- missing genes (BUSCO completeness 80–95%);
- uncollapsed haplotigs or retained isoforms (BUSCO "duplicated" often 5–20% in heterozygous insect assemblies);
- TE-derived ORFs.

Training proteomes (UniProt human and yeast) have none of these.

Test design: on development data (yeast and human, the clean-trained model, never the frozen test), apply seeded perturbations at 2–3 levels:

| Perturbation | Levels |
|---|---|
| (i) Truncate a fraction of genes to 30–80% of length (N-terminal, C-terminal or internal deletion) | 10%, 25% |
| (ii) Split genes into two non-overlapping models | 5% |
| (iii) Fuse random gene pairs | 2% |
| (iv) Delete genes | 5%, 10%, 20% |
| (v) Add near-identical copies (99–100% identity, small indels: haplotig or isoform mimic) | 2%, 5% |

Report:
- caller recall of labelled pairs;
- Spearman correlation of each feature (clean vs perturbed) on surviving pairs;
- AUPRC on surviving pairs, plus end-to-end AUPRC;
- overlap of the top-k predictions;
- the same metrics for the identity baseline B1.

Separately, flag in the CLI output every pair with ≥99% identity over ≥95% of the shorter protein ("possible allele, isoform or assembly duplicate"). An identity-driven model will otherwise rank these artefacts at the top in a heterozygous assembly. Human has only 509 such pairs among about 636k candidate pairs, so the flag costs nothing on clean input.

**M6. Stratify evaluation by WGD vs SSD using evaluation-only annotations, and report SSD-only AUPRC.**
Labelled sets are dominated by WGD pairs:
- human screened pairs are 73% WGD;
- the yeast sources Dean, DeLuna and Kuzmin are mostly ohnologs;
- Cusack found α-WGD origin to be a top predictor in Arabidopsis.

Insects have few WGDs, and plant WGDs are much younger than yeast or 2R WGDs. Use YGOB (yeast) and OHNOLOGS v2 (human) only to report AUPRC within the WGD and SSD strata. Use them also to validate the WGD-like proxy features (age bin, 1:2 orthology pattern, small family), with the proxy checked against these labels on development data only. Treat the SSD-only number as the best available proxy for insect-like genomes. This shows whether the model has learned "WGD-ness", which will not transfer. Arabidopsis α/β/γ annotation may be used for stratification only at the final test access.

**M7. Replace single-linkage family definitions; this affects both features and splits.**
In the human self-search at E ≤ 1e-3, single linkage puts 10,650 of 20,652 proteins into one component. With ≥0.3 coverage on both sequences it is still 5,812. With ≥0.5 coverage on both, the largest component is 941.

Domain sharing (kinase, ZNF, WD40 and similar) chains unrelated families together. As specified, "connected-component size" is therefore nearly a binary "is in the giant component" feature. The split clustering in §6.1 (E ≤ 1e-3, coverage ≥0.3, single linkage, over human + yeast) will create a giant cluster that the >5% rule forces into development, which leaves a small and biased test set.

For family-size features, use coverage-aware families (≥0.5 coverage of both sequences, or of the shorter one) or a non-chaining clustering such as MMseqs2 `cluster` or MCL-style. Keep a looser "distant relatives" count as a separate feature. Check the size of the largest cluster before freezing splits.

### Should

**S1. Weight conservation counts by clade.**
Count presence or absence over about 12 major clades (vertebrates, other deuterostomes, arthropods, nematodes plus lophotrochozoans, early metazoans plus choanoflagellates, Saccharomycotina, other fungi, land plants, green algae, SAR, Amoebozoa, excavates), or weight species by phylogeny. Raw species counts in this panel mostly measure "present in vertebrates" for human queries and "present in fungi" for yeast queries.

**S2. Compute the Dennler-style "species rank" at species level.**
Use "fraction of panel species with any protein scoring higher to A than B does" instead of a protein count. Panel species with their own recent WGDs or expansions inflate protein counts: *X. laevis* (allotetraploid), *Paramecium* (several WGDs), *Trichomonas* (50k proteins), maize.

**S3. Make duplication dating aware of gene loss.**
"A and B share the same best hit in S" is also what you see when S lost one copy of an older duplication, and differential loss is the rule after WGD. Place the duplication node at the deepest split for which some outgroup species has positive evidence: two distinct co-orthologs, each reciprocal-best with one of A and B. Treat a shared best hit only as weak evidence. RBH status can be read from the forward table, since the query with the top bitscore to O is O's best query hit, once M1 removes truncation, so no reverse search is needed.

**S4. Correct divergence asymmetry for distance and alignment region.**
- Use the nearest outgroup (the shallowest species where A and B share a reciprocal-best ortholog), not any outgroup.
- Convert identities to distances (Scoredist, or a Poisson/Kimura correction). Raw identity saturates below about 30%, which is where human 2R pairs sit against invertebrates and yeast pairs sit against distant fungi.
- Compute A–O and B–O identity over the region that all three sequences share. Otherwise one fragmented or multidomain copy produces spurious asymmetry.

**S5. Use features that do not depend on proteome size.**
MMseqs2 E-values scale with database size; proteomes here range from 5k (*S. pombe*) to more than 40k (rice, maize) and over 100k (polyploid crops). Define the calling threshold and family counts by bitscore (or by E-value at a fixed effective DB size), and drop raw E-value as a model feature. Keep percentile features only if LOSO in both directions (human and yeast differ 3.4× in size) shows they transfer.

**S6. Add explicit outgroup copy-number and pattern features.**
- The 1:2 orthology pattern: the number of outgroup species where a single-copy ortholog is reciprocal-best to A or B and has no second close homolog.
- Family "duplicability": mean per-species copy number of the family across the panel.

Both are sequence-only and are the cleanest available proxies for ohnolog-like, dosage-sensitive retention.

**S7. Add fragment and fusion indices, plus a panel-derived completeness QC.**
- Per gene: query length divided by the median length of its panel orthologs. This can be a feature and a QC flag.
- Per run: a BUSCO-like report from panel-derived near-universal single-copy families (present and single-copy in ≥90% of panel species), giving the fraction present, fragmented and duplicated.

Print both in the run log so users can tell whether predictions on their proteome are trustworthy.

**S8. Handle missing values explicitly.**
Orphans, lineage-specific families and genes with no outgroup ortholog are common in insects and plants. Their Tier B features are undefined. Missingness differs by species: yeast has more missing values than human because it has fewer close relatives. Use explicit indicator features and test whether the model relies on missingness, for example by AUPRC on pairs with complete Tier B features.

**S9. Run label-free distribution-shift diagnostics on unlabelled plant and insect proteomes.**
Compare feature distributions between the training species and unlabelled target-like proteomes, using adversarial validation or PSI per feature. Use rice, maize, *Tribolium*, *Anopheles*, and the added *Apis* and *Bombyx*. This uses no labels, so it is a legitimate input to feature selection, and it directly addresses the stated use case.

**S10. Document what the self-exclusion rule actually does.**
At 98% identity and 50% of proteins, the rule removes the great apes, not only "the query species or a strain". That is acceptable, but it should be stated, and the threshold should be justified or tied to the M3 thinning conditions. Otherwise the default human evaluation has a much denser panel (macaque at ~29 My) than any plant or insect user will have.

### Could

**C1. Higher-order redundancy flag.**
In allopolyploid crops (wheat, *Brassica napus*, cotton), homeolog triplets make pairwise SL largely meaningless. Output a flag when a third paralog is about as close as the pair partner. The design already computes the "max identity to other paralog" and "shared third paralogs" features this needs.

**C2. Two-tier paralog output.**
Report full-length homologs (coverage ≥0.5 of the shorter) separately from pairs that share only a domain. The 636k human candidate pairs are mostly domain sharing. Report labelled-pair recall per tier, and give the genome-wide prevalence estimate for the tier that matches Ensembl-like paralogs.

**C3. Calling-threshold rule.**
Do not choose the calling threshold by "maximise recall of labelled pairs" alone, because that always picks the loosest setting. Use recall at a fixed output size, or require ≥99% recall of labelled pairs at the strictest threshold that achieves it.

**C4. Phylostratigraphy caveat.**
Short, fast-evolving proteins look young through homology-detection failure (Moyers & Zhang 2015/2016). Report gene-age features with length and rate controls, and do not interpret their importance causally.

**C5. Avoid WGD peak-calling features.**
Do not build features from WGD peaks fitted to the query's own paralog-age distribution, such as mixture models on protein distance. Plant peaks overlap (Arabidopsis α/β/γ), insects have none, and the calibration per proteome is unstable. If tried at all, use it only as an evaluation-only diagnostic.

## 3. Risks the design does not address

1. **No insect labels, and plant labels are Arabidopsis only.** Arabidopsis redundancy labels are enriched for young α-WGD pairs with high identity, so the one unseen-plant test mostly rewards identity. The M3 thinning conditions and M6 SSD-only numbers are the only evidence about insects. REPORT.md should say this plainly.
2. **TE-derived and contaminant ORFs** in new plant and insect annotations can form large, high-identity pseudo-families that inflate family-size features and fill the top of the ranking. M5 (v) and S7 partly cover this. A low-complexity or repeat-family flag would cover the rest.
3. **Tandem arrays** are common in plants (about 15% of Arabidopsis genes) but are invisible without coordinates. The model cannot tell them apart from dispersed SSDs. This is a limitation to state, not something to fix.
4. **Release mismatch.** The panel is QfO 2026_02 (UniProt 2026_02) while the training proteomes are from 2026_03. Self-exclusion works by identity, so it is harmless, but the pinned versions must agree with DATA.md.
5. **Search cost after M1.** Per-species coverage raises alignment work. Saturated queries average about 634 hits at `--max-seqs 3000` versus 100. Budget this before freezing the default runtime target of ≤2 h.

## 4. Things that are fine as they are

- Not inferring WGD from synteny, and using YGOB and OHNOLOGS only for evaluation (keep it that way).
- Using only sequences of other species (no label transfer by orthology), and the "best hit shared" plus asymmetry ideas, with the refinements above.
- MMseqs2 as the homology engine for both self-search and panel search, and DIAMOND kept off the panel because of its 19 GB footprint.
- A generic, logged self-exclusion rule rather than species-specific code.
- LOSO in both directions, and Arabidopsis as the first unseen species with the frozen-test discipline.
- Ruling out synteny, Ks and CDS-based dating: these are correctly out of scope for a protein-only input.
