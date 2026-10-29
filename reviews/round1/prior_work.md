# Round 1 review: prior work (paralog SL and duplicate-gene evolution)

Reviewer role: prior-work expert. Read: PROMPT.md, DESIGN.md v1, LITERATURE.md, docs/lit/methods.md, docs/lit/human_labels.md, docs/lit/yeast_and_other_species_labels.md, data/raw/qfo/species_list.tsv, data/raw/_survey/Costanzo2016-Readme.txt. Key claims spot-checked on the web (Europe PMC records and PMC full text for Dennler & Ryan 2025 PMC12034104; Dean 2008 PLoS Genet full text; VanderSluis 2010 abstract PMC3010121; Flister 2025 Europe PMC record; Kuzmin 2020 abstract/Dryad).

## 1. Overall verdict

The design is well grounded in prior work. It splits by homology family rather than by random pair (fixing De Kegel 2021's main weakness), uses only tested negatives, avoids SLKB and SynLethDB, includes Dennler-style rank features, and conditions PLM similarity on identity. Four problems need fixing before labels are frozen:
- Two yeast label rules contradict the labels' own definition of a positive (viable Kuzmin/Dean pairs used as blanket negatives; essential-gene TS alleles treated as testable).
- Duplication mode (WGD vs SSD), the best-established determinant and the main ascertainment axis, has no evaluation-only check or stratified reporting.
- The reference panel's resolution for dating duplications differs greatly between human, yeast and the unseen species.
- The identity baseline B1 may not include the identity variant prior work actually used, and leaving it out would flatter the model.

## 2. Concrete changes

### MUST

**M1. Do not use Kuzmin 2020 or Dean 2008 *viable* pairs as blanket negatives. Label them with the same ε/fitness rule as SGA.**
- DESIGN §2.3 says "Kuzmin 2020 (viable WGD double mutants) and Dean 2008 viable pairs serve as additional tested negatives".
- But the primary yeast positive is ε < −0.35, which is *synthetic sick*, not lethal. Kuzmin's 240 WGD pairs are viable by construction. The paper says pairs that were "essential or refractory to double-mutant query strain construction" were left out. About 30% of duplicate pairs show negative digenic interactions (Kuzmin), and Dean found 34% of measured pairs had significantly reduced double-mutant fitness, beyond the 24% that were lethal.
- Calling all viable pairs negative would plant false negatives among high-identity WGD pairs, which is exactly the region where the model must beat identity.
- Fix: compute ε (or the Dean W_AB deficit) from Kuzmin Table S5 SMF/DMF and Dean W values. Apply the same positive, negative and ambiguous thresholds as for SGA. Viable-but-sick pairs become positive or ambiguous.

**M2. Treat yeast pairs with an essential member (TS/DAmP alleles only) as untestable in the primary tier, matching the human rule.**
- Human pairs where either single gene is essential are "untestable" (DESIGN §2.3). Yeast pairs instead include "TS alleles of essential genes … as tested".
- A negative interaction between a hypomorph of an essential gene and a paralog deletion is not SL of two nulls. The yeast survey says so itself (docs/lit/yeast_and_other_species_labels.md, caveat ii and line 121: report "essential-member pairs" separately). DESIGN dropped that recommendation without saying why.
- Mixing the two label semantics across training species corrupts the jointly trained model. Keep essential-member pairs as a separate, reported sensitivity set only.

**M3. Stratify all results by duplication mode, and validate the WGD-proxy features, using evaluation-only annotations.**
- Prior work (Guan 2007; Hakes 2007; De Kegel 2019/2021; Cusack 2021, where α-WGD origin was the top feature) makes WGD vs SSD the most consistent determinant. It is also the largest ascertainment axis: 73% of screened human pairs are WGD vs 19% of all pairs.
- The envisioned insect applications are mostly SSD or tandem (Drosophila has had no recent WGD). A model that is mainly a WGD-age detector could look good on yeast and human and fail there.
- Required:
  - (a) AUPRC and Δ vs B1 per stratum: WGD vs SSD, using YGOB for yeast and OHNOLOGS v2 / Makino & McLysaght for human. These are evaluation-only labels, never features. docs/lit/methods.md §5 already says this, but DESIGN omits it.
  - (b) An evaluation-only check of how well the duplication-age features (Tier B, items 4–5) separate WGD from SSD pairs.
  - (c) Yeast results with and without ribosomal-protein (RP) pairs (from SGD annotation, evaluation only). Dean 2008: RP pairs were redundant in 30/36 cases vs 39/165 non-RP pairs. RP pairs are a large share of yeast WGD positives and form many separate small families, so the homology split does not control them. They have a distinctive sequence signature (short, basic, ultra-conserved, near-identical), and human RP genes are mostly singletons.

**M4. Make the duplication-age and conservation features robust to panel density (and see S1 on densifying the panel).**
- The QfO panel (data/raw/qfo/species_list.tsv) gives human about 20 vertebrates/chordates that bracket the 2R events (Branchiostoma and Ciona before; Lepisosteus, Danio and Oryzias after). It also includes chimp, gorilla, macaque, mouse and others.
- By contrast:
  - **Yeast:** the only Saccharomycotina left after self-exclusion are *C. albicans* and *Yarrowia*. There is no close pre-WGD outgroup (Kluyveromyces, Zygosaccharomyces, Lachancea, Eremothecium), so the yeast WGD cannot be told apart from SSDs spanning >200 My.
  - **Arabidopsis:** no other eudicot remains (maize, rice, Physcomitrella, Chlamydomonas only).
  - **C. elegans:** no other nematode.
  - **Insects:** 3 insects plus Daphnia and Ixodes.
- Features counted in "number of panel species" or "deepest clade" therefore mean different things per query species. Their signal will come mostly from human, and they shift on unseen species.
- Fix:
  - Define age features in species-count-free units. Example: the fraction of the query's own RBH-ortholog identities to panel species that are lower than id(A,B), i.e., "the duplication is younger than the split from s" expressed on an identity scale. Dennler's "species" feature is a count; also normalise it by the number of panel species available.
  - Add a panel-thinning sensitivity analysis: recompute features with the nearest panel relatives removed down to a fixed divergence. This simulates a phylogenetically isolated new genome and shows how much performance depends on close relatives.

**M5. B1 must include the identity variants used by the prior work it is compared against.**
- De Kegel and Dennler use Ensembl-style identity: identical positions divided by the length of each protein, i.e., length-normalised, not local fident. De Kegel uses the max of the two directions; Dennler uses the min.
- MMseqs2 local `fident` over a short shared domain can be high for domain-only pairs, which weakens B1 and inflates "clearly beats identity".
- Add to the B1 candidate set:
  - identical positions / len(A), identical positions / len(B), and their min and max;
  - global (end-to-end) alignment identity for each candidate pair.
- Choose among them on development folds only, as planned.

### SHOULD

**S1. Densify the reference panel evenly (sequences only, UniProt reference proteomes, pinned in DATA.md).**
- Add roughly 2–4 species per under-sampled clade:
  - Saccharomycetaceae pre- and post-WGD (*K. lactis*, *Z. rouxii*, *L. thermotolerans*, *E. gossypii*, *N. glabratus*);
  - eudicots (Vitis, Populus or Solanum, Medicago, Brassica rapa);
  - insects (Apis, Bombyx, Aedes, Acyrthosiphon, Nasonia);
  - nematodes (C. briggsae, Pristionchus or Brugia).
- This is the realistic setting for a new plant or insect genome, which will have sequenced relatives. Rate-asymmetry analyses of yeast ohnologs (Byrne & Wolfe 2007) depend on a close non-WGD outgroup.
- Keep it even across clades so the training species are not specially favoured. Expected panel cost is about +30–40% sequences, which is small next to the current 9-minute panel search.

**S2. Make the 1:2 co-orthology (WGD-like) pattern and pair-level retention explicit features.**
- For each panel species s, compute:
  - whether A and B share a single best hit O_s that is itself a singleton (no close paralog of O_s in s), i.e., an ancestral singleton duplicated in the query lineage;
  - whether s has two distinct co-orthologs, one reciprocal-best to A and one to B (the pair was retained across that lineage).
- Summarise as counts and fractions over panel clades.
- Use RBH-based ortholog presence, not "any hit at E ≤ 1e-3", for the conservation and age features (item 4). Any-hit counts saturate for multi-domain proteins and large families, which is a known phylostratigraphy artefact.
- This is the sequence-computable version of De Kegel's "conservation score" and "has S. cerevisiae / S. pombe ortholog", and of the WGD proxy listed in docs/lit/methods.md feature #10. It is not in DESIGN §4.

**S3. Correct LITERATURE.md where it misstates prior work.**
- **Dennler & Ryan 2025 main human labels.** Their main human "SL" label is the De Kegel DepMap-inferred pipeline over all Ensembl 111 pairs (about 0.5% SL; 116/23,734 in the notebook), not screens. The screen-based set is their "SL-Lenient" (403/5,316). Their "identity performed exceptionally well" refers to the DepMap-based task, which LITERATURE §1 labels "Human (screens)".
- **Dennler rank features.** Only the "paralogs" rank is nearly orthogonal to identity (ρ ≈ 0.13–0.18). The "species" feature correlates at ρ ≈ 0.58–0.68. LITERATURE says "rank features weakly correlated (ρ≈0.15)".
- **Dennler embeddings.** ProtT5 generally beat ESM-2, and the pooling method mattered more than the distance metric. DESIGN Tier C does not reflect either point.
- **"Redundancy decays with divergence" is contested.**
  - Dean 2008 reports redundancy that is "pervasive and persistent" (about 100 My), with no measured divergence trend.
  - Li 2010 reports a decay proportional to divergence time.
  - Gu 2003 concerns single-gene dispensability, not SL.
- **"WGD pairs buffer more than SSD" is weaker than stated.** Among non-RP pairs, Dean found 17/120 WGD vs 4/45 SSD SL, which is not significant.
- **The heteromer evidence points in opposite directions in the two training species.** In yeast, Diss 2017 and Dandage & Landry 2019 report heteromeric paralogs as *dependent*, i.e., less buffering. In human, De Kegel 2019/2021 found physically interacting / same-complex paralogs *more* often SL. This matters for any future predicted-heteromer feature.
- **Kebabci 2026 chance level.** LITERATURE says AP 0.06 vs 0.04 chance; docs/lit/methods.md says AP 4–6% with chance about 2% on Klingbeil. Reconcile.
- **Flister 2025 model type.** The published version reports an XGBoost model; the survey describes the preprint's weighted RF. Neither reports a held-out AUPRC or an identity baseline, so this is not a comparator.

**S4. Pre-register Δ vs B2 (identity + family size + closest pair) as a secondary comparison.**
- The sequence-only subset of De Kegel's model has never been reported. B2 is the strongest cheap prior-work baseline (docs/lit/methods.md says the same).
- If the final model beats B1 but not B2, the expensive tiers add nothing, and the report must say so. This does not change the PROMPT's success criterion. It reports the comparison prior work implies.

**S5. Add a "≥2 studies" human sensitivity tier, and cross-check Dennler against the matching tier.**
- De Kegel's screen consensus required SL in ≥2 *studies*. DESIGN's primary tier allows ≥2 cell lines from one study, where library- or pipeline-specific artefacts replicate.
- Report a ≥2-studies tier for comparability with De Kegel's external test (50/407).
- Specify that the Dennler agreement table compares their SL-Lenient set (any screen) with our *lenient* tier, and their DepMap "SL" set with nothing (it is a different label type).

**S6. Make the yeast measurement rule symmetric across measurements.**
- Positive is "ε < −0.35 and p < 0.05 in any orientation" over all alleles and temperatures. Negative is "never ε < −0.08". Pairs with more measurements (both orientations, several alleles or temperatures) get more chances to become positive and fewer to stay negative. This is a multiple-testing ascertainment that correlates with gene properties.
- Restrict the primary tier to deletion alleles of nonessential genes (see M2). Use a per-pair consensus, e.g., the deletion-allele measurement with the smallest p, or require agreement when both orientations exist. Report the "any" rule as a sensitivity analysis.

**S7. Report BioGRID low-throughput "Synthetic Lethality" positives that lie outside the SGA-tested universe separately.**
- These are literature-curated and biased toward well-studied genes, and they have no matched negatives. Primary metrics with and without them.

**S8. Treat Arabidopsis "non-redundant" negatives the way human single-essential pairs are treated.**
- Cusack negatives (110, shared by RD4 and RD9) are pairs where a single mutant is at least as severe as the double mutant. Some have a single-mutant lethal phenotype, which the human rule would call untestable rather than negative.
- Report Arabidopsis AUPRC with and without negatives where a single mutant is already lethal, if the Cusack or Lloyd & Meinke class data allow it. Note that RD4 ⊂ RD9 and that the positives come from literature, so they are biased toward well-studied genes.

**S9. Keep feature models trained on measured interactions or annotations out of every tier, including "remaining ideas".**
- PLM-interact is trained on measured PPI data, so it imports forbidden data, possibly including the labelled pairs themselves. Only AF-Multimer-type models (trained on structures) fit the PROMPT.
- For "fraction disordered", prefer a composition- or structure-statistics method (an IUPred-style energy estimate, or charge/hydropathy). Avoid predictors trained on curated annotation sets.

### COULD

- **C1.** Run a ProtT5 ablation restricted to proteins in labelled pairs (Dennler found ProtT5 better than ESM-2). Also compare mean-of-residue pooling with BOS/EOS pooling. If it helps, decide whether its cost can be justified for the default run.
- **C2.** Check the SL rate in the top identity bin per human source, and for the C. elegans RNAi set. At very high identity, guides or dsRNAs can hit both paralogs, so a single knockout becomes a double. The pair then shows up as "single essential / untestable" or non-SL, producing an artificial non-monotone identity–SL curve. Guide sequences are in the supplementary tables.
- **C3.** Count tandem (adjacent) yeast paralog pairs among labelled vs all pairs, using genome coordinates for evaluation only. SGA tests linked pairs poorly. The Costanzo raw readme describes only the reciprocal-sign filter, not linkage handling, so check the SOM. Tandem SSDs are probably under-represented in training but common in plants and insects.
- **C4.** Keep Kuzmin's position-specific rate-correlation ("evolutionarily constrained" pairs) as a remaining idea. It is the one sequence-only determinant from that paper not covered. A subset version (small families only) may be affordable.

## 3. Risks the design does not address

1. **Distribution shift to the target applications.** Training pairs are dominated by WGD/2R, closest pairs and small families, while insect paralogs are mostly SSD/tandem and plant paralogs come from several WGDs plus tandem arrays. Without M3 stratification, and some SSD-only result, the "unseen species" claim rests on 51 Arabidopsis and 16 C. elegans positives.
2. **Label-noise ceiling.** Same-line replicate agreement (precision ≈ 0.40, recall ≈ 0.32) and about 31% between-screen overlap cap achievable AUPRC. Report an empirical ceiling (e.g., AUPRC of one screen's calls predicting another's on shared pairs) so a modest absolute AUPRC is not misread.
3. **Hub families in human screens.** Some families appear in many screens and many positive pairs. The B3 gene-prior baseline is good; also report a per-family-capped AUPRC (each homology cluster weighted equally) so a few large positive families do not dominate the test result.
4. **Giant component.** Single-linkage at E ≤ 1e-3 with coverage ≥ 0.3 over the human and yeast proteomes may chain kinases, WD40, ZNF and similar families into one component. Forcing >5% clusters into development would leave a test set of small, simple families, where SL is easier to predict and more common. Report test vs development composition (family size, identity, WGD fraction).
5. **Family-size truncation.** `--max-seqs 1000` truncates hit lists in large families, so family-size and rank features saturate. This is probably harmless, since SL is rare there, but it should be measured on the toy and human data.
6. **Self-exclusion side effects.** For human, the ≥98%-identity / ≥50% rule will probably also remove chimp and gorilla. That is acceptable, but the set of removed species per query should be reported, since it changes panel density (M4).

## 4. Things that are fine as they are

- The homology-cluster split across species, the random-split leakage estimate, and the gene-prior baseline B3 correctly answer Park & Marcotte and De Kegel's random CV.
- Tested negatives only; ambiguous pairs excluded; lenient tier kept as a sensitivity analysis. This matches De Kegel and Dennler practice.
- SLKB and SynLethDB are excluded as label sources. The De Kegel DepMap set is kept separate. The yeast primary threshold ε < −0.35, p < 0.05 matches Dennler & Ryan (4.8% prevalence), so their yeast set is a usable cross-check.
- The feature list covers the established sequence-computable determinants:
  - identity variants;
  - family size and closest/reciprocal-best pair;
  - paralog rank and species rank (Dennler);
  - conservation and age;
  - divergence asymmetry against an outgroup (VanderSluis 2010; Byrne & Wolfe 2007; Kuzmin 2020);
  - length ratio;
  - triangle redundancy in families of 3 or more.
- PLM and 3Di similarity are evaluated conditional on identity. This is the right response to Dennler's and TripleThreat's finding that these measures largely proxy identity.
- The end-to-end evaluation, with missed pairs scored at the minimum, and the reweighting to genome-wide prevalence (Flister: SL is "infrequent", <0.5%) are both appropriate.
- Proteome-percentile features and no species indicator are suitable for cross-species transfer.
