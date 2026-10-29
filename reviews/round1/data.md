# Round 1 review — Data curation (genetic-interaction labels)

Reviewer role: data curator for combinatorial CRISPR screens, yeast SGA, and Arabidopsis/worm/pombe genetic data.
Read: PROMPT.md, DESIGN.md v1, LITERATURE.md, docs/lit/human_labels.md, docs/lit/yeast_and_other_species_labels.md, docs/lit/methods.md §1.5/§4.3.
Checked by hand: small files in data/raw/_survey/ (Costanzo readme, yeast UniProt gene2acc/idmapping, Dede S2, Chou T1, Thompson SD1, CHyMErA S8 headers, Cusack datasets, Chou 13-pair list).

## 1. Overall verdict

The label strategy is sound in outline. It uses primary tables, tested negatives only, an excluded ambiguous middle, per-source evidence and no SL databases. But the human and yeast definitions do not yet describe the same biology, and there are four concrete defects:
- Yeast includes TS/DAmP alleles of essential genes, while human treats a singly essential gene as untestable.
- The yeast positive rule is a pure ε cut, and ε cannot reach −0.35 when the single mutants are sick. The human rule has an absolute DKO criterion, but only for Ito.
- The human "≥2 lines, any study" rule is confounded with how many lines each library was screened in.
- The label unit is the UniProt accession, which silently merges paralogs that encode identical proteins.

All four are fixable before M1. None needs new data sources.

## 2. Concrete changes

**D1 (must) — Essential-gene pairs are untestable in yeast too; null alleles only in the primary tier.**
Replace "TS alleles of essential genes included as tested" with the following rules:
- The primary yeast tier uses only deletion (null) alleles of non-essential genes. Strain suffixes: `_sn`, plus DMA arrays.
- Any pair with an essential member, i.e. one measured only through `_tsq`/`_tsa`/`_damp` alleles, is labelled `untestable_essential`.
- TS/DAmP (ExN/ExE) calls are kept as a separately reported stratum.

*Why:* the human rule (single essential → untestable) and the survey (A.1 caveat ii, A.8 pre-filter) both say SL is undefined when one single mutant is lethal. A TS-hypomorph negative GI is not null-allele SL, and including it means the two species' labels differ in meaning.

**D2 (must) — One SL definition for both species: significant negative interaction AND a strong absolute double-mutant defect.**

*Yeast.* Use ε < −0.08, p < 0.05 **and** DMF ≤ τ as the primary positive.
- Calibrate τ once, before any modelling, from the SGA DMF distribution of Dean 2008 W_AB = 0 pairs. The colony-size DMF has a non-zero floor.
- Freeze τ in `config/labels.yaml`.
- Keep "ε < −0.35, p < 0.05" as a secondary tier so results stay comparable with Dennler & Ryan (methods.md §1.5: prevalence ≈ 4.8%).
- Under the multiplicative model ε ≥ −SMF_A·SMF_B. With SMF = 0.6 each, a completely lethal double mutant has ε = −0.36, and slightly sicker singles can never pass −0.35. The pure ε rule therefore systematically excludes SL pairs whose singles are sick.
- Pairs with min(SMF) below a pre-set floor (e.g. 0.5) have too little dynamic range to be confident negatives. Mark them ambiguous.

*Human.* Apply the absolute DKO criterion (DKO LFC or Z below the screen's essential/non-essential midpoint, as in Kebabci and Chou's ZLFC < −2) to **every** study where single/DKO LFCs exist, not only Ito.
- Checked: Chou Supplementary Table 1 holds only dLFC/ZdLFC/RdLFC. ZLFC needs the figshare LFC tables.
- Checked: Dede Table S2 holds only zdLFC, with no p-value, no DKO LFC and no single-gene essentiality flag. So "authors' significance criterion" and "study's own essentiality flag" cannot be applied to Dede as written.

**D3 (must) — A human aggregation rule that does not reward being screened in many lines.**
"SL in ≥2 testable lines (any studies)" lets a pair with 2/27 SL lines (Harle) or 2/22 (Klingbeil) become a robust positive. A pair tested in Parrish's or CHyMErA's 2 lines can only be positive by being SL in both.
- Positives are then enriched for many-line libraries.
- Negatives are enriched for few-line libraries.
- Those libraries were chosen by different rules: Parrish ≥50% identity, 2-member families; Harle ≥45% identity plus a lethal fly/worm ortholog; Klingbeil signalling/epigenetic domains.
- So the label partly encodes library membership, which sequence features can recover.

Adopt the survey's rule (human_labels.md §4.3):
- Positive = SL in ≥2 *distinct* cell lines **and** (≥2 independent studies **or** penetrance ≥0.25).
- Negative = ≥2 testable distinct lines, SL in none, **and** no borderline score. Borderline means within a pre-set margin, e.g. zdLFC < −2 or GI < −0.3.
- The same line in two studies (A549: Dede/Ito/in4mer; MeWo: Thompson/Ito; A375: Thompson/in4mer; RPE-1: Thompson/CHyMErA) counts as one line with replicated evidence.

Add a scripted audit to `results/labels/`: positive rate by source study and by number of testable lines.

**D4 (must) — Keep Flister 2025 genome-wide pairs out of the primary human tier.**
DESIGN lists Flister among the primary sources. Its 36,648 pairs × 2 lines would contribute about 35k "negatives", swamp the other screens and move prevalence below 1%. As the survey recommends, keep it as a separate realistic-prevalence set, and use it only if the full per-pair non-hit table is verified. Preprint Table S3 has "Lethal" calls only, and NaN there conflates not-SL with not-tested.

**D5 (must) — The label unit is the gene, not the UniProt accession; build training proteomes per gene.**
Checked in `UP000002311_559292.gene2acc/idmapping`: five yeast ohnolog pairs encode identical proteins and share one accession. They are HHF1/HHF2 (P02309), HHT1/HHT2 (P61830), TEF1/TEF2 (P02994), TIF1/TIF2 (P10081) and EFT1/EFT2 (P32324). These are classic redundant, likely-SL pairs at 100% identity.

In human, Dede screened HIST1H4H_HIST4H4, and both genes encode histone H4 (one accession). If the unit is "pair of UniProt accessions", these pairs become self-pairs and vanish, removing exactly the top-identity positives.

The PROMPT input is "one protein per gene". So:
- Build the human, yeast and test proteomes per gene from gene2acc, duplicating identical sequences under gene IDs.
- Key labels by (species, gene_id, gene_id), with the accession as an attribute.
- Log every collapse.
- Also keep a toy test with identical paralogs (the caller and identity features must handle 100% / E = 0).

**D6 (must) — Define the yeast labelled universe independently of the caller, before caller tuning.**
In human, screens define which pairs were tested as paralogs. In yeast, SGA tests essentially all pairs, so the "paralog" universe has to come from somewhere.
- If it comes from our caller, then §3's "adjust thresholds once to maximise recall of labelled pairs" is circular: recall is 100% by construction.
- Freeze a reference paralog universe for label construction only. Examples: Ensembl Compara paralogs of a pinned release (as Dennler & Ryan used, 5,673 yeast pairs), or a fixed BLASTP rule documented in `config/labels.yaml`. Annotation used to define the universe is not a feature.
- Use the same rule to restrict human screen pairs to "paralogs". This drops Thompson's 447 mutual-exclusivity pairs, Harle/Chan non-paralog pairs and Klingbeil duplicates.
- Then measure caller recall against this universe.

**D7 (must) — "Viable" is not "negative" when positives include synthetic-sick pairs.**
DESIGN uses Kuzmin 2020 viable WGD double mutants and Dean 2008 viable pairs as extra negatives wherever SGA gives no call. Yet the primary positive (ε < −0.35) includes strong synthetic-sick pairs, and Dean reports 34% of pairs with significantly reduced fitness.
- Apply the same ε/DMF rule (D2) to Kuzmin Table S5 fitness and to Dean W values. Viable-but-strongly-sick pairs are ambiguous, not negative.
- Note also that Kuzmin's 240 pairs are the *constructible* 44% of WGD pairs. The others were left out partly because the double mutant was essential, so they are a selected sample.

**D8 (must) — BioGRID low-throughput positives: restrict, de-duplicate and stratify.**
- Use only `Throughput = Low` and `Experimental System = Synthetic Lethality`.
- Exclude records from Costanzo 2016/2010, Dean, Kuzmin and other HT publications (no double counting).
- Exclude records where either gene is essential (D1), and records whose Qualifications indicate a drug, temperature or other condition, or a non-null allele, where parseable.
- BioGRID has no negatives. Pairs added *only* by BioGRID-LT raise prevalence in exactly the strata SGA cannot test (linked, well-studied, essential-like). Tag them `positive_lt_only` and report yeast AUPRC with and without them.
- Pin BioGRID 5.0.262.

**D9 (must) — A complete per-source rule table in `config/labels.yaml` before any label is computed.**
The primary tables have free choices that are currently unspecified. Each must be fixed in advance, not chosen after seeing counts:
- *Thompson:* D14 vs D28, and RRA vs t-test FDR. SD1 has both timepoints and both tests, plus "Passes filter?" columns.
- *CHyMErA:* early vs late timepoint. Orientation (gi1 Cas9/gi2 Cas12a, "both_gi"). **Exclude the HAP1 + Torin1 drug condition** (it is in Table S8).
- *Klingbeil:* collapse repeated domain-combination rows to one pair call. Pre-register any/most-significant/median; "any" inflates SL.
- *in4mer:* prototype and Inzolia libraries are separate screens.
- *Chou 2025:* the rescored Dede/Parrish/Thompson are a sensitivity re-scoring only. Each screen enters the primary tier once, with one scoring.
- *Ito:* state that the Kebabci "GEMINI ≥90th percentile within screen" is a *relative* threshold. DESIGN rejects SLKB for using relative thresholds, so justify the difference (it is paired with an absolute LFC cut) or replace it with an absolute GEMINI score cut.
- *Dede:* the untestable rule, since Table S2 has no essentiality flag (D2).

**D10 (should) — SGA measurement multiplicity and orientation conflicts.**
"Positive if ε < cut in *any* orientation" combined with "negative only if *never* ε < −0.08" pushes heavily measured pairs out of negatives and into positives. Measurements multiply across alleles, query/array orientation, NxN at different temperatures, and Costanzo 2010 vs 2016 strains.
- Pre-register an aggregation rule: e.g. the most significant null-allele measurement per orientation, then require agreement or no strong contradiction across orientations.
- Report label rate vs number of measurements.
- Write the conflict rule exactly: e.g. one orientation passes the positive rule, the other has ε > −0.04 with p > 0.05 → `discordant`.
- Note that the raw file has already dropped reciprocal pairs with opposite-sign significant scores (Costanzo readme). That loss is invisible, so say so in DATA.md.

**D11 (should) — Mark linked yeast pairs untestable.**
In SGA, double mutants of physically linked genes are rarely recovered by meiotic recombination. This produces artifactual "negative interactions". The scoring pipeline is reported to mask a window around the query; I did not verify the window size.
- Tandem duplicates are high-identity paralogs, so this artifact would correlate directly with the main feature.
- Using chromosome coordinates from a pinned SGD/Ensembl GFF (label construction only), mark same-chromosome pairs within the documented linkage window as untestable, unless Dean or BioGRID-LT gives a tetrad-based call.
- State the consequence: yeast training labels contain almost no tandem duplicates.

**D12 (must) — Unseen-species endpoints with stated semantics; Arabidopsis negatives re-checked for lethal singles.**
- *Arabidopsis.* The Cusack RD4 negatives are "non-redundant": a single mutant is at least as severe as the double. Where that single mutant is itself lethal, the pair is *untestable* under our human/yeast definition, not a negative. Arabidopsis enters the pre-registered success criterion, so:
  - Obtain the Lloyd & Meinke 2012 single-mutant phenotype classes (Suppl. table; OUP blocks scripts, so a documented manual download).
  - Drop negatives with a lethal single mutant, or report AUPRC both ways, with the choice pre-registered.
  - Freeze RD4 (lethal) as primary and RD9 (any enhanced phenotype) as secondary.
- *C. elegans (Tischler).* Positives are "synthetic phenotype", including brood size, not lethality. Label it as such. 16 positives give a very wide CI: report it as descriptive only.
- *S. pombe.* Tested negatives need the raw E-MAP S-score matrix, which is unverified. Pre-register the threshold (e.g. S < −2.5) and a drop rule if no tested-negative matrix is obtainable.
- Map worm 2006 sequence names to current WBGene with logged merges/splits.

**D13 (should) — Identifier mapping hub and loss accounting.**
- *Human:* symbol → HGNC ID. Use previous and alias symbols from a pinned HGNC complete-set date. Dede uses pre-2020 names such as HIST1H4H; the survey notes PALM2AKAP2. Drop, never guess, aliases that resolve to more than one current gene.
- HGNC/Ensembl gene → gene2acc of the pinned UniProt release. Assert that every mapped accession is in that release's reference-proteome FASTA.
- *Yeast:* strip the allele suffix from strain IDs. Drop merged/dubious ORFs (BioGRID removed 21,039 such interactions from the HighestConfidence set). Map via `Gene_OrderedLocusName`.
- Report pairs lost per source and per step in `results/labels/mapping_report.tsv`.

**D14 (should) — Commit policy: binary calls for all sources, an open-only variant, and named sources of record.**
- Commit per-pair, per-source *binary/ternary call codes and counts* for every source, which are facts derived by our thresholds. Never commit scores, LFCs or table extracts from non-open sources (Ito, CHyMErA) or NC-ND sources (Parrish, Klingbeil, Flister journal). There is precedent: Dennler & Ryan redistribute Ito- and CHyMErA-derived SL calls under CC BY 4.0.
- Also build and commit `labels_open.tsv` from CC BY/CC0/MIT sources only. Report the headline human result on it as a sensitivity analysis, so a user who refuses the derived non-open calls can still reproduce a result.
- Note in DATA.md that the project is EU-based (Irish database right on substantial extraction). This is a judgement call, not legal advice.
- *Costanzo:* use the BioGRID mirror (MIT) as the source of record. TheCellMap's 521,571,629 B is 497.4 MiB, so it is probably the same file as BioGRID's "497.41 MB". Confirm by checksum. The Boone-lab SOM zip (431,580,254 B) is a different "representative subset" and must not be mixed in.
- *Kuzmin:* the Dryad CC0 files cannot be fetched by script (401/403). Document the manual browser step and verify the checksum in the build script.

**D15 (should) — Report how much the excluded middle inflates AUPRC.**
Excluding everything between "strong SL" and "clearly neutral" makes the task easier than whole-proteome use, where the tool scores those pairs too. In yeast, roughly 25% of duplicate pairs fall in −0.35 < ε < −0.08.
- Report a sensitivity analysis with ambiguous pairs as negatives, clearly labelled.
- Report a family-weighted AUPRC so large families do not dominate.
- This does not change the primary metric. It is a pre-registered secondary.

**D16 (should) — Label-noise ceilings in both species.**
Add yeast cross-source agreement to the label QC: SGA call vs Dean 2008 W, vs Kuzmin S5, vs BioGRID-LT, on overlapping pairs. This mirrors the human between-screen agreement. Better still, keep Dean out of the yeast positive definition for one QC run, so it acts as an independent validation of the SGA-derived calls.

**D17 (should) — Fix the source list at M1, with fallbacks.**
Several sources are unverified: Flister full table, Klingbeil AACR tables, Parrish Table S4, Chan Supp. Tables 5–7, pombe S-scores, Harle 32 MB file.
- Freeze which ones are in at M1. Any later addition needs a decision note *before* modelling.
- Pre-register the fallback for each, e.g. Parrish from GEO GSE178179 counts with our own GI scoring.

**D18 (should) — A cross-species label-semantics table in DESIGN §2.3.**
Columns: assay, readout, allele type, context, positive definition, negative definition, untestable definition, expected prevalence.
- Human: proliferation of cancer lines over 2–4 weeks; CRISPR nulls (incomplete); context-specific.
- Yeast: haploid colony size, rich medium.
- Dean: growth rate / inviability.
- Arabidopsis: organismal embryo/gametophyte/seedling lethality.
- Worm: any synthetic RNAi phenotype.

This makes explicit what "SL" means in each test set, and why prevalence and AUPRC are not comparable across species.

**D19 (could) — Flag cross-reactivity at very high identity.**
- High nucleotide identity lets one sgRNA or dsRNA hit both paralogs. The single "KO" then looks lethal, so the pair is untestable or shows no GI.
- Libraries may also have dropped near-identical pairs (ascertainment at the top of the identity range).
- Flag pairs above about 90% CDS identity where the study gives guide specificity. Check Tischler's cross-RNAi filter. Report label rates in the top identity bin.

**D20 (could) — Flag closely linked human pairs.**
Dual cutting of tandem neighbours can excise the intervening genes, which can produce false SL. Flag pairs within about 1 Mb, using coordinates for label construction only.

**D21 (could) — Smaller additions.**
- Cross-check the yeast calls with the Dennler & Ryan yeast set, which is on the same Zenodo record.
- Add GABI-DUPLO lethal / no-phenotype lines if the Wiley SOM is obtainable, de-duplicating against the 29 Bolle pairs already in Cusack.
- Cusack has strong gene reuse: 41 genes occur in more than one RD9 pair. Bootstrap the Arabidopsis CIs by family, not by pair.

## 3. Risks the design does not address

1. **Ascertainment differs fundamentally between the training species.**
   - Human labels come from curated libraries: 89% closest pairs, 73% WGD, median family size 3, and per-library identity floors (Parrish ≥50%, Harle ≥45%, Flister ≥20%).
   - Yeast labels come from a near-exhaustive SGA over a paralog universe of our choosing, dominated by large families and remote pairs.
   - Species-balanced weights equalise mass, not distribution. The model may learn "closest pair / small family → SL" from human and "low identity → negative" from yeast.
   - Report per-species metrics within shared strata: RBH/closest-pair subset; family size bins.
2. **Source is predictable.** Per-library selection rules mean "which study screened this pair" may be recoverable from sequence features. Since label rates differ by study, this is a leakage path. Audit it: train a classifier to predict source study from the features, and report AUPRC within single-source subsets.
3. **No tandem duplicates in the yeast or GABI-DUPLO labels** (SGA linkage; DUPLO is unlinked pairs only), while plant genomes are rich in tandem arrays. The envisioned plant application is out of distribution on exactly this axis. State it as a limitation.
4. **Context.** Yeast labels are from a single rich-medium condition (Musso 2008 shows condition-dependent epistasis). Human labels are cancer-line proliferation. "SL" in the deliverable means "robust in standard growth conditions". Say so in the CLI documentation.
5. **The 13-pair "gold standard" is not independent.** It was assembled from hits in the very screens we label from, so passing that check validates the code, not the biology. DESIGN already frames it correctly; keep it that way in REPORT.md.
6. **Unseen-species leakage through homology.** "Unseen species are held out entirely" does not stop a pombe or Arabidopsis pair from being orthologous to a training family. Report unseen-species AUPRC stratified by whether either gene has a significant hit to a gene in a labelled training pair. This matters most for pombe.
7. **PDF parsing of Dean 2008.** The survey parse gives 197 pairs / 48 lethal vs the paper's 201 / 49. Make the parser assert against the paper totals and fail loudly. A silent parse loss of a positive in a 48-positive set is material.

## 4. Fine as it is

- Labels come from primary tables with tested negatives only. SLKB (relative top-10%, GPL) and SynLethDB (mixed predicted evidence) are rejected. Dennler & Ryan serve only as a cross-check. The De Kegel DepMap-inferred set is kept separate.
- Per-pair, per-source evidence is retained so tiers can be regenerated. The ambiguous tier is excluded from primary metrics, and a lenient tier serves as sensitivity.
- UniProt reference proteomes and idmapping come from a single pinned release (2026_03), with mapping losses logged (D5 and D13 refine this).
- Unseen-species sets are test-only. Drosophila is an enrichment check only. The absence of insect/crop AUPRC data is stated as a limitation.
- Raw tables are never committed; raw data are fetched by script with checksums in DATA.md.
- D1 (no measured data as features or auxiliary targets) is correctly scoped. Using essentiality, coordinates or phenotypes to *construct labels* (D1, D11, D12 above) does not violate it, because none of these enter the feature vector.
- The yeast negative definition already has a margin (never ε < −0.08) between the positive and negative thresholds. Human negatives need the same margin (D3).
