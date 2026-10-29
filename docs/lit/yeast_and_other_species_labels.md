# Yeast and other-species SL paralog label sources

Survey date: 2026-10-27. Scope: (A) labels for *S. cerevisiae* paralog SL; (B) candidate held-out
"unseen species" test sets; (C) reference proteomes and identifier mapping.

Status markers:
- **[V]** I checked this myself: the URL was reachable, or I inspected the file or the paper text.
- **[L]** Reported in the literature or on a web page, but I did not inspect the data.
- **[U]** Unverified or my own estimate. Check before relying on it.

Nothing was bulk-downloaded. The small files inspected are listed in §5.

---

## 0. Decision-relevant summary

1. **Yeast primary label source: the Costanzo 2016 SGA raw pairwise file [V].**
   - It covers about 23 M double mutants, and every tested pair carries ε, p-value, both single-mutant fitnesses (SMF), the double-mutant fitness (DMF) and its SD.
   - Because every tested pair is included, **tested negatives are available**, which is essential for AUPRC.
   - Size: 521 MB zip from TheCellMap, or 497 MB from the BioGRID mirror. The license is unclear (see A.1).
   - The BioGRID "HighestConfidence" and tab2 derivatives contain **only significant interactions**, so they hold no negatives. Do not use them as the label universe.
2. **Independent yeast low-throughput sets with explicit lethal / viable double mutants.** Use them as a second yeast test, or as checks on the SGA labels.
   - **Dean et al. 2008 [V]:** 283 pairs listed, 197 with complete fitness data in the PDF parse; 48 of them have DMF = 0.
   - **Musso 2008 [L]:** 399 WGD pairs, 140 aggravating.
   - **Kuzmin 2020 [L]:** 240 WGD pairs whose double deletion is *viable* by construction.
3. **Yeast positive rate is about 25–35% for duplicates.**
   - Musso: 140/399 = 35% aggravating.
   - Dean: 49/201 = 24% strictly lethal.
   - De Kegel 2021 cites 25–35%.
   - Expected labelled yeast paralog pairs: about 1–2.5k [U]. That number depends on the paralog definition and has to be computed.
4. **Best unseen-species test sets, ranked:**
   1. *Arabidopsis*: GABI-DUPLO, plus Cusack 2021 / Lloyd & Meinke 2012. About 300–450 pairs with tested double mutants; 51 (RD4) to 190 (RD9) positives depending on definition [V for Cusack].
   2. *S. pombe*: E-MAP / BioGRID. Quantitative, but a fungus, so not very "unseen".
   3. *C. elegans*: Tischler 2006 [L]. 143 pairs tested by combinatorial RNAi, 16 synthetic.
   4. *Drosophila*: FlyBase / BioGRID curated [V]. Positives only, noisy, and no systematic paralog screen exists.
   5. *C. albicans* and mouse: small or close to human.
   6. Rice and maize: essentially no data. BioGRID has 2 and 0 genetic interactions.
   - **No systematic insect paralog double-mutant dataset was found.** This is a real gap and should be stated as a limitation.
5. **Proteomes:** UniProt reference proteomes release 2026_03 [V] give one canonical protein per gene directly (`*.fasta.gz`), with `*.gene2acc.gz` and `*.idmapping.gz` for mapping IDs. That makes them the simplest uniform choice for all species.

---

## A. *S. cerevisiae* label sources

### A.1 Costanzo et al. 2016, global SGA network

| Item | Detail |
|---|---|
| Ref | Costanzo M. et al. *Science* 353:aaf1420 (2016). DOI 10.1126/science.aaf1420. PMID 27708008 |
| Landing page | https://thecellmap.org/yeast/costanzo2016/ [V, 200]. The old `https://thecellmap.org/costanzo2016/` is now 404 [V], but the data-file paths below it still work. |
| Main file | `https://thecellmap.org/costanzo2016/data_files/Raw%20genetic%20interaction%20datasets:%20Pair-wise%20interaction%20format.zip`. HTTP 200, Content-Length 521,571,629 B, Last-Modified 2020-01-15 [V]. |
| Mirrors | Boone-lab SOM: `https://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/data_files/Data File S1_Raw genetic interaction datasets: Pair-wise interaction format.zip` (431,580,254 B [V]). This is the "representative subset" version used in the paper, not the complete 23 M set. BioGRID: https://downloads.thebiogrid.org/BioGRID/Published-Datasets/Costanzo2016/ (497.41 MB raw; 5.60 MB HighestConfidence; 17.13 MB BioGRID tab2) [V]. |
| Contents | Separate files for NxN (nonessential × nonessential), ExN/NxE, ExE (TS alleles of essential genes) and DAmP alleles, plus a strain list with SMF [V, TheCellMap page]. |
| Columns (11) [V, readme] | Query Strain ID; Query allele name; Array Strain ID; Array allele name; Arraytype/Temp (DMA or TSA, 26 °C or 30 °C); ε (genetic interaction score); P-value; Query SMF; Array SMF; Double mutant fitness; DMF SD |
| Identifiers | Strain IDs are systematic ORF names with a suffix, e.g. `YAL001C_tsq508`, `_damp`, `_sn…` [L]. Strip the suffix to get the systematic ORF name, which maps to UniProt `Gene_OrderedLocusName` / SGD. |
| Thresholds (authors' own) [V, Boone SOM page] | lenient: P < 0.05; intermediate: P < 0.05 and \|ε\| > 0.08; stringent: P < 0.05 and (ε < −0.12 or ε > 0.16). In the raw file, reciprocal pairs (A–B and B–A) that were both significant at the intermediate threshold with opposite signs were already removed [V, BioGRID readme]. |
| SL definition | The authors give no explicit "synthetic lethal" cutoff. A negative interaction is ε < 0, with SL as the extreme case [L]. A DMF-based SL tier can be built from the DMF column (our choice, see A.7). |
| Tested negatives | **Yes.** All tested pairs are present, including non-significant ones. |
| Scale | ~23 M double mutants; ~550 k negative and ~350 k positive interactions; ~90% of genes covered as query and/or array. Screens: 3,589 NxN queries × 3,892 arrays; 1,162 TS queries (804 essential genes) × DMA; ExE: 2,241 deletion + 1,108 TS queries × 792 TS arrays (561 essential genes) [V, TheCellMap page]. |
| License | TheCellMap gives none [U]. The BioGRID mirror states **MIT** for "all BioGRID project files" [V]. Cite Costanzo 2016 and Stark 2006. Treat it as redistributable via the BioGRID mirror, but store only derived labels plus the checksum of the raw file. |
| Caveats | (i) Pairs of genes physically linked on a chromosome are poorly tested by SGA because of linkage-based depletion. Tandem SSD paralogs may be systematically missing or filtered [U; check the SOM]. (ii) TS / DAmP alleles are hypomorphs, so an ExE "negative GI" is not SL of null alleles. (iii) Colony-size DMF has a non-zero floor, so DMF ≈ 0 is not reached exactly. |

### A.2 Kuzmin et al. 2020, WGD paralog digenic and trigenic analysis

| Item | Detail |
|---|---|
| Ref | Kuzmin E. et al. *Science* 368:eaaz5667 (2020). DOI 10.1126/science.aaz5667. PMID 32586993 |
| Data | Dryad doi:10.5061/dryad.g79cnp5m9, **CC0-1.0** [V via Dryad API]. Version 4 (2020-08-05) has TableS1.xlsx (61.5 MB, digenic/trigenic raw), TableS2 (0.63 MB), TableS3 (10.9 MB), and TableS4–S13 (11–567 KB). **Dryad file downloads now return 401/403 to scripted requests** [V]. Download them in a browser, or use the Science SOM: http://science.sciencemag.org/content/suppl/2020/06/24/368.6498.eaaz5667.DC1 [L]. |
| Design [V, paper text] | "240 double-mutant query yeast strains, each deleted for a pair of nonessential WGD paralog genes … dispensable paralog pairs and represent 44% (240/551) of unique WGD paralog pairs, a number of which were not included because the pair was either essential or refractory to double-mutant query strain construction." Fitness of the 240 double mutants and the 480 single mutants is in **Table S5** (colony size). Thresholds: digenic and trigenic negative/positive at \|ε\| or \|τ\| > 0.08, P < 0.05. Overall: 537,911 double mutants, 256,861 triple mutants, 4,650 negative digenic interactions, 2,466 negative trigenic interactions. |
| SL relevance | **The 240 pairs are viable double mutants by construction, so they are tested non-SL for a strict "lethal" label.** They can still be synthetic sick. The DMF in Table S5 gives the degree. **Table S13** relates to "paralogs that show a synthetic lethal genetic interaction" (the essential subset) [L; content not inspected]. The 311 WGD pairs left out contain the SL pairs but are not individually labelled. |
| IDs | ORF / standard names [U] |

### A.3 DeLuna et al. 2008

| Item | Detail |
|---|---|
| Ref | DeLuna A. et al. *Nat Genet* 40:676–681 (2008). DOI 10.1038/ng.123. PMID 18408719 |
| Data | Competition-based fitness of single and double deletions of WGD paralogs, focused on metabolism. Supplementary material is behind the Nature paywall [U]. BioGRID curates **73 interactions** from this paper [L]. |
| Size / def. | Exact number of pairs **not verified** [U]. DeLuna 2010 (PLoS Biol, PMC2846854) re-used these data and defined SSL as **ε < −0.2**, giving 27 neutral and 16 SSL pairs in minimal medium [L]. Dean 2008 rich-medium data gave 37 neutral and 18 SSL pairs. |
| Negatives | Yes, all measured pairs. Small. |

### A.4 Dean et al. 2008

| Item | Detail |
|---|---|
| Ref | Dean E.J., Davis J.C., Davis R.W., Petrov D.A. *PLoS Genet* 4:e1000113 (2008). DOI 10.1371/journal.pgen.1000113. PMC2440806. **CC BY** (PLoS) |
| Design [L, PMC text] | 289 duplicate pairs from two-member families with no third homolog (FASTA E < 0.01): 204 WGD and 85 SSD (reciprocal best hits, E < 1e-10). 32 pairs were not built because one gene is essential, and 56 were unfinished. **201 pairs measured; 49 (24%) synthetic lethal (W_AB = 0)**; 69 (34%) significantly reduced fitness. 90 singleton control pairs, none SL. Growth rate in YPD. |
| Data [V] | Dataset S1 PDF: `https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=10.1371/journal.pgen.1000113.s004` (176,698 B). Columns `ORF1 ORF2 ORF1_fit ORF2_fit ORF1ORF2_fit`, using systematic ORF IDs. My pdftotext parse: 283 pair rows, **197 with all three fitness values, of which 48 have ORF1ORF2_fit = 0**, and 86 rows with IDs only. This differs slightly from the paper's 201 / 49, probably from page-break parsing. Re-check when building labels. |
| Use | High-quality, independent of SGA, includes **SSD pairs**, and gives **explicit lethal versus viable** outcomes. Best yeast check for the SGA-derived labels. |

### A.5 Musso et al. 2008 (added, relevant)

Musso G. et al. *Genome Res* 18:1092–1099 (2008), DOI 10.1101/gr.076174.108, PMC2493398. Of **399** examinable WGD pairs, **140** showed aggravating (negative) epistasis in standard conditions, and more appeared under stress [L]. The per-pair table should be in the Genome Research supplement [U].

### A.6 BioGRID genetic interactions (low and high throughput)

| Item | Detail |
|---|---|
| Release | **BioGRID 5.0.262** (Oct 2026) [V]. `BIOGRID-ORGANISM-5.0.262.tab3.zip` (180.12 MB, all organisms split by species). Archive: https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-5.0.262/ [V]. **MIT license** [V]. |
| Yeast size | *S. cerevisiae* S288c: 600,374 raw / 457,599 non-redundant genetic interactions [V, BioGRID stats page] |
| Relevant evidence codes | `Synthetic Lethality`, `Synthetic Growth Defect`, `Negative Genetic` (quantitative screens, mostly SGA / E-MAP), `Synthetic Rescue`, `Dosage Lethality`, `Phenotypic Enhancement`. The tab3 file also has a `Throughput` column (Low / High) and `Score` / `Qualifications` fields. |
| IDs | Systematic names (ORF), Entrez Gene ID, and SWISS-PROT/TREMBL accessions in tab3 [L] |
| Negatives | **None.** BioGRID only records reported interactions, so it is positives only. |
| Use | Low-throughput `Synthetic Lethality` records are high-specificity positives. They can rescue SL pairs that SGA cannot test (linked pairs, pairs with an essential-like double mutant), but they are biased toward well-studied genes. Use them as positives only, with negatives coming from SGA or Dean tested-neutral pairs. |

### A.7 YGOB ohnolog list (for describing pairs only, not a feature)

Byrne & Wolfe *Genome Res* 15:1456 (2005), DOI 10.1101/gr.3672305: **551** S. cerevisiae ohnolog pairs [L]. Host `ygob.ucd.ie` **refused connections on port 80/443 today** [V: connection refused], so the site may be temporarily down. Fallbacks:
- the Genome Research supplement;
- the 551-pair list implied by Kuzmin 2020 (who use the same 551).

Use YGOB only to annotate pairs as WGD or SSD for stratified reporting, never as a model input.

### A.8 Recommended yeast labelling scheme (our proposal; thresholds are our choice, flagged)

**Unit.** An unordered paralog pair {A, B}, using systematic ORF names mapped to UniProt canonical accessions. The paralog universe comes from our own sequence-only homology pipeline, so the same definition is used at train and inference time. Labels are attached afterwards.

**Measurements per pair.** Collect every Costanzo raw record for A×B and B×A over all alleles and screens, and keep the deletion (`_dma`/`_sn`) alleles of nonessential genes. Add the Dean 2008 W values, Kuzmin Table S5 DMF, and BioGRID low-throughput records.

**Pre-filter (exclude from the labelled set).**
- Either gene essential as a single deletion (SMF unavailable, or only TS/DAmP alleles). SL is undefined when a single mutant is already lethal. Report these separately as "essential-member pairs".
- Dubious or merged ORFs.

**Positive (SL / SSL), primary label.** At least one measurement with ε < −0.12 and P < 0.05 (the authors' *stringent* negative threshold), **or** a Dean 2008 W_AB = 0, **or** a BioGRID low-throughput `Synthetic Lethality` record.

**Strict-SL tier (sensitivity analysis).**
- Positive additionally requires DMF ≤ 0.3, or Dean W_AB = 0, or low-throughput SL.
- Negatives are as below.
- DMF ≤ 0.3 is [U], our choice, and should be checked against the distribution of DMF for known lethal pairs (Dean W = 0) inside the SGA data.

**Negative.**
- The pair was tested in at least one orientation with a deletion allele, and **no** measurement meets even the lenient-negative criterion (ε < −0.08 and P < 0.05).
- Also negative: a Dean 2008 pair with no significant fitness reduction.
- Kuzmin pairs whose DMF is near the multiplicative expectation can serve as additional negatives, but only for the strict tier.

**Grey zone (excluded from the primary evaluation).** −0.12 ≤ ε < −0.08 with P < 0.05, or a pair that was never tested.

**Conflicts.**
- One orientation stringent-negative and the other clearly neutral (ε > −0.04 with ≥ 2 replicates): mark `discordant`. Keep it as positive in the primary set (SGA false negatives are common) but flag it, and report AUPRC with and without discordant pairs.
- Low-throughput SL versus SGA neutral: low-throughput wins, flagged.
- Dean versus SGA: Dean's lethal call wins for the strict tier.
- Never resolve conflicts by looking at model scores.

**Positive rate and counts (estimates, [U]).**
- WGD pairs: 551. About 35% aggravating (Musso), about 24% strictly lethal (Dean).
- Under a BLAST-style paralog definition (≥ 20–25% identity, SSD + WGD), yeast probably has about 1–2.5 k paralog pairs. Large families produce many pairs, so cap family size or weight by family.
- Expect about 60–80% to be covered by SGA (missing: essential members, linked pairs, untested strains).
- **Expected labelled pairs ≈ 0.8–1.8 k; positive rate ≈ 20–35% (primary), ≈ 8–15% (strict).**
- These are **not computed**. Compute them once the paralog list exists.

---

## B. Candidate UNSEEN-species test sets

| Rank | Species | Source(s) | Pairs with tested outcome | Positives | Negatives available? | SL definition | Data URL / status | License |
|---|---|---|---|---|---|---|---|---|
| 1 | *Arabidopsis thaliana* (plant) | **Cusack et al. 2021** benchmark, built on Lloyd & Meinke 2012 plus Bolle 2013 | 300 single + double mutant trios (RD9) [V] | RD9 "redundant" (DM more severe than both singles): 190 [V]. RD4 "extreme" (singles show no phenotype, DM lethal): 51 [V] | **Yes.** 110 "non-redundant" (a single mutant ≥ DM) [V] | DM lethal = Lloyd & Meinke class 2 "essential" (embryo, gametophyte or seedling lethal) | Zenodo 10.5281/zenodo.3987384, CC-BY-4.0 [V]. Dataset_1.txt = RD4 (161 pairs: 51 pos / 110 neg); Dataset_2.txt = RD9 (300 pairs: 190 pos / 110 neg). Column 1 `pair_ID` (`AT1G..._AT2G...`), column 2 `Class`; the remaining ~4,100 columns are measured-data features, **not usable**. RD4 ⊂ RD9 [V]. 552/553 genes map to UniProt UP000006548 [V]. | CC-BY-4.0 (Zenodo); code GPL-3.0 |
| 1 | *A. thaliana* | **GABI-DUPLO** (Bolle et al. *Plant J* 75:157, 2013, DOI 10.1111/tpj.12197) | 200 double knockouts of unlinked two-member paralog pairs (from 2,108 pairs) [L] | 23 lethal before germination + 6 seedling lethal = **29 lethal** [L]; 7 adult phenotype; 14 conditional or molecular | **Yes.** The remaining ~144 show no visible phenotype [L/U] | DM lethality, unbiased systematic design | Per-line table in the Wiley supplement [U, paywalled]. DUPLOdb at gabi-kat.de is **offline** (server compromised) [V]. NASC collection page lists 367 lines by DUPLO pair ID only [V]. Obtaining the pair → phenotype table needs the Wiley SOM or the authors. | Wiley copyright [U] |
| 1 | *A. thaliana* | Lloyd & Meinke 2012 *Plant Physiol* 158:1115, Suppl. Table S6 | 347 groupings / 591 genes; 188 doubles in simple clusters [L] | Positive-biased (only redundant genes with a multiple-mutant phenotype) | Partly (asymmetric / symmetric classes) | Phenotype classes: 0 / conditional / cellular / morphological / essential | `plphys_v158_3_1115_s1.zip` (1.3 MB) [L]; OUP returned 403 to scripts | OUP copyright [U] |
| 2 | *S. pombe* (fungus, no WGD) | Ryan et al. 2012 **Mol Cell** 46:691 (not *Cell*), DOI 10.1016/j.molcel.2012.05.028: genome-scale E-MAP, ~60% of the nonessential genome. Also Roguev 2008 *Science*, Dixon 2008 *PNAS* (SpSGA), BioGRID | Unknown; S. pombe has few paralogs. Probably about 50–200 tested paralog pairs [U] | BioGRID pombe: 52,872 non-redundant GIs in total [V]; 40,831 curated from Ryan 2012 [L] | Only via the **raw E-MAP S-score matrix** (Mol Cell SOM) [U]; BioGRID has positives only | E-MAP S-score (e.g. S < −2.5 negative) [U] | BioGRID tab3 (MIT); Mol Cell SOM [U] | MIT (BioGRID) |
| 3 | *C. elegans* (nematode) | Tischler, Lehner, Chen & Fraser 2006 *Genome Biol* 7:R69, DOI 10.1186/gb-2006-7-8-r69, PMC1779603 | **143** duplicate pairs (of 239 InParanoid pairs) [L] | **16** (11%) synthetic phenotype, i.e. combined phenotype > product of singles (embryonic lethality or brood size) [L] | **Yes** (127) | Synthetic RNAi phenotype, not strictly lethal | Additional file 2 (`gb-2006-7-8-r69-S2.doc`) lists the screened pairs; Table 4 lists the 16 [L]. Lehner 2006 *Nat Genet* (~65 k pairs, signalling genes) is not paralog-focused. BioGRID *C. elegans*: 2,302 GIs [V]. | CC BY (BMC) |
| 4 | *Drosophila melanogaster* (insect) | FlyBase 2026_03 curated genetic interactions and genotype→phenotype; BioGRID; combinatorial RNAi screens (Horn 2011, Fischer 2015 *eLife* 4:e05464, Heigwer 2018 *eLife* 7:e40174, all signalling-focused) | No systematic paralog screen found [V by search]. Must be mined from multi-allele genotypes | FlyBase `genotype_phenotype_data_fb_2026_03.tsv.gz`: 400,940 rows, 51,437 "lethal*" rows, 10,437 lethal rows with ≥ 2 alleles [V]. Most of these are GAL4 driver + RNAi, not double mutants, so they need heavy filtering. Gene-level GI file: 21,550 rows (suppressible / enhanceable) [V]. BioGRID Dmel: 10,764 non-redundant GIs [V] | Weak. "viable" annotations exist, but tested-viable double mutants are rarely curated | Literature curation; hypomorphs and RNAi | https://s3ftp.flybase.org/releases/current/precomputed_files/ [V] | FlyBase CC-BY-4.0 [U, check] |
| 5 | *Candida albicans* | Shapiro et al. 2018 *Nat Microbiol* 3:73 (CRISPR gene-drive double deletions: efflux pumps, adhesins); BioGRID 459 GIs [V] | Few tens [U] | Few | Few | Growth | Nature SOM | — |
| 6 | *Mus musculus* | MGI `MGI_PhenoGenoMP.rpt` (54 MB [V]): multi-gene genotypes with MP terms (lethality terms) | Could be hundreds [U] | Lethal double knockouts | "no abnormal phenotype" genotypes are sparse | MP lethality terms | https://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt [V] | MGI terms (free, cite) |
| — | Rice, maize, other crops, other insects | None found. BioGRID rice 2 GIs, maize 0 [V]. Multi-Knock (Hu et al. 2023 *Nat Plants*) is a multi-targeted CRISPR *Arabidopsis* library with no pairwise SL labels | — | — | — | — | — | — |

**Ranking rationale.**
- *Arabidopsis* is the only plant with tested double-mutant outcomes. It has explicit negatives and a lethal-phenotype readout, and is the most distant lineage from yeast and human. That makes it the most valuable held-out set.
- Recommended *Arabidopsis* label:
  - Positive = DM lethal (Cusack RD4 positives ∪ GABI-DUPLO lethal lines).
  - Negative = DM tested and not lethal (Cusack non-redundant ∪ GABI-DUPLO no-phenotype lines ∪ RD9-positive pairs whose DM is not lethal).
  - Report RD9 (any enhanced DM phenotype) as a secondary "redundancy" endpoint.
  - Expected size about 300–450 unique pairs with about 15–25% lethal positives [U]. Deduplicate overlaps between Cusack (which already includes 29 Bolle pairs) and GABI-DUPLO.
- *S. pombe* gives quantitative, SGA-like labels. It is fungal and shares many orthologs with *S. cerevisiae*, so it tests species transfer, not lineage transfer.
- *C. elegans* (Tischler) is tiny but has clean negatives and is the only systematic non-human animal paralog screen found.
- *Drosophila* is the only insect option. It is positives-only and curation-biased, so treat it as a "recall at top-k / enrichment" sanity check, not an AUPRC test.
- **Leakage control for all of them:** report performance on test pairs whose families have no close homolog among labelled train pairs. Use a cross-species family clustering, e.g. both members' best hits to labelled yeast or human pairs below an identity threshold.

---

## C. Reference proteomes (UniProt release 2026_03, released 2026-09-02 [V])

Base URL [V; all files below returned HTTP 200]:
`https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/<UPID>/<UPID>_<taxid>{.fasta.gz,.gene2acc.gz,.idmapping.gz}`

- `<UPID>_<taxid>.fasta.gz` = **canonical sequences only, one per gene.** Isoforms are in the separate `_additional.fasta.gz`.
- `.gene2acc` maps gene → accession, with MOD IDs preferred: SGD:S…, PomBase:SP…, TAIR:AT…, CG… for fly.
- `.idmapping` carries the cross-references: `Gene_OrderedLocusName`, SGD, Ensembl, FlyBase, WormBase, and so on.

| Species | UPID | Taxid | Genes (= canonical seqs) | Total entries | Assembly | fasta.gz size | Label ID → proteome mapping |
|---|---|---|---|---|---|---|---|
| *S. cerevisiae* S288c | UP000002311 | 559292 | 6,066 (6,066 FASTA records [V]) | 6,067 | GCA_000146045.2 | 2.0 MB | ORF (YAL001C) via idmapping `Gene_OrderedLocusName` (6,072 rows) or SGD ID [V] |
| *H. sapiens* | UP000005640 | 9606 | 20,652 | 147,520 | GCA_000001405.29 | 7.7 MB | Ensembl gene / Entrez / HGNC via idmapping (27.7 MB) |
| *S. pombe* 972h- | UP000002485 | 284812 | 5,129 | 5,130 | GCA_000002945.2 | 1.7 MB | PomBase systematic ID (SPAC…) in gene2acc [V] |
| *A. thaliana* | UP000006548 | 3702 | 27,496 | 39,272 | GCA_000001735.1 (TAIR10) | 7.9 MB | AGI (TAIR:AT1G01010) in gene2acc. 27,496 gene groups [V]; Cusack IDs 552/553 mapped [V] |
| *D. melanogaster* | UP000000803 | 7227 | 13,814 | 21,947 | GCA_000001215.4 | 5.0 MB | gene2acc uses CG numbers [V]; FBgn via idmapping `FlyBase` |
| *C. elegans* | UP000001940 | 6239 | 19,792 | 26,629 | GCA_000002985.3 | 5.7 MB | WormBase WBGene / sequence name via idmapping |
| *C. albicans* SC5314 | UP000000559 | 237561 | 6,040 | 6,041 | GCA_000182965.3 | 2.0 MB | CGD orf19 IDs via idmapping (check) |
| *M. musculus* | UP000000589 | 10090 | 21,860 | 54,855 | GCA_000001635.9 | 8.1 MB | MGI IDs via idmapping |
| *O. sativa* japonica | UP000059680 | 39947 | 43,678 | 48,896 | GCA_001433935.1 | 9.3 MB | RAP-DB / MSU IDs via idmapping |
| *Z. mays* | UP000007305 | 4577 | 39,346 | 63,237 | GCA_902167145.1 | 9.9 MB | Zm00001eb… |
| *Tribolium castaneum* (extra insect) | UP000007266 | 7070 | 16,568 | 18,505 | GCA_000002335.3 | 5.1 MB | TC IDs |

**Alternatives.**
- Ensembl release 116 (vertebrates, yeast): `https://ftp.ensembl.org/pub/current_fasta/<species>/pep/<Species>.<asm>.pep.all.fa.gz` [V].
- Ensembl Genomes release 63 (plants, metazoa, fungi): `https://ftp.ensemblgenomes.ebi.ac.uk/pub/{plants,metazoa,fungi}/current/fasta/<species>/pep/` [V]. Examples: `Arabidopsis_thaliana.TAIR10.pep.all.fa.gz`, `Drosophila_melanogaster.BDGP6.54.pep.all.fa.gz`, `Schizosaccharomyces_pombe.ASM294v2.pep.all.fa.gz`.
- `pep.all` contains **all isoforms**. Choosing one per gene needs the `Ensembl_canonical` tag from the GTF, or longest-isoform selection.
- The Ensembl BioMart service returned "Service unavailable" today [V].
- TAIR (arabidopsis.org) returned 403 to scripted access and is partly subscription-based [V/U].

**Recommendation.** Use UniProt reference proteome canonical FASTA for every species: one protein per gene, uniform headers, consistent IDs. Record release 2026_03 and the SHA-256 of each file in DATA.md. For a "newly sequenced proteome" simulation, also test with longest-isoform-per-gene derived from Ensembl `pep.all`, to check robustness to isoform choice. Filter sequences < 30 aa and fragments (`(Fragment)` in the header) consistently across species.

---

## 5. Files inspected in `data/raw/_survey/` (all < 5 MB)

| File | Source | What I checked |
|---|---|---|
| `Costanzo2016-Readme.txt` | BioGRID Costanzo2016 readme | Column definitions, thresholds, filtering |
| `cusack2021_Dataset_1.txt` (3.1 MB) | Zenodo 3987384 | Header and `Class` counts: 51 positive / 110 negative (RD4) |
| `cusack2021_Dataset_2_id_class.txt` | Zenodo 3987384 Dataset_2 (streamed, kept columns 1–2 only) | 190 positive / 110 negative (RD9) |
| `cusack2021_RD9_test_set.txt` | ShiuLab GitHub | 30 held-out pair IDs |
| `dean2008_DatasetS1_duplicate_pairs.pdf` | PLoS Genet SOM s004 | 283 pairs; 197 with fitness; 48 DMF = 0 |
| `flybase_gene_genetic_interactions_fb_2026_03.tsv.gz` | FlyBase | 21,550 rows; types suppressible / enhanceable |
| `flybase_allele_genetic_interactions_fb_2026_03.tsv.gz` | FlyBase | Free-text interaction statements |
| `flybase_genotype_phenotype_data_fb_2026_03.tsv.gz` | FlyBase | 400,940 rows; 10,437 lethal multi-allele rows |
| `UP000002311_559292.{fasta,gene2acc,idmapping}.gz` | UniProt 2026_03 | 6,066 records; ID types |
| `UP000006548_3702.gene2acc.gz`, `UP000002485_284812.gene2acc.gz`, `UP000000803_7227.gene2acc.gz` | UniProt 2026_03 | Gene ID formats |

The other files in `_survey/` (chymera, dede, thompson, in4mer, chou2025, flister, ens111 human) were put there by other agents and were not inspected by me.

---

## References

1. Costanzo M. et al. (2016) A global genetic interaction network maps a wiring diagram of cellular function. *Science* 353:aaf1420. doi:10.1126/science.aaf1420. Data: https://thecellmap.org/yeast/costanzo2016/ ; https://downloads.thebiogrid.org/BioGRID/Published-Datasets/Costanzo2016/
2. Costanzo M. et al. (2010) The genetic landscape of a cell. *Science* 327:425. doi:10.1126/science.1180823. Baryshnikova A. et al. (2010) *Nat Methods* 7:1017 (SGA score).
3. Kuzmin E. et al. (2020) Exploring whole-genome duplicate gene retention with complex genetic interaction analysis. *Science* 368:eaaz5667. doi:10.1126/science.aaz5667. Data: doi:10.5061/dryad.g79cnp5m9 (CC0)
4. DeLuna A. et al. (2008) Exposing the fitness contribution of duplicated genes. *Nat Genet* 40:676. doi:10.1038/ng.123
5. DeLuna A. et al. (2010) Need-based up-regulation of protein levels in response to deletion of their duplicate genes. *PLoS Biol* 8:e1000347 (PMC2846854)
6. Dean E.J., Davis J.C., Davis R.W., Petrov D.A. (2008) Pervasive and persistent redundancy among duplicated genes in yeast. *PLoS Genet* 4:e1000113. doi:10.1371/journal.pgen.1000113
7. Musso G. et al. (2008) The extensive and condition-dependent nature of epistasis among whole-genome duplicates in yeast. *Genome Res* 18:1092. doi:10.1101/gr.076174.108
8. Byrne K.P., Wolfe K.H. (2005) The Yeast Gene Order Browser. *Genome Res* 15:1456. doi:10.1101/gr.3672305
9. Stark C. et al. (2006) BioGRID. *Nucleic Acids Res* 34:D535. BioGRID 5.0.262: https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-5.0.262/
10. De Kegel B., Quinn N., Thompson N.A., Adams D.J., Ryan C.J. (2021) Comprehensive prediction of robust synthetic lethality between paralog pairs in cancer cell lines. *Cell Syst* 12:1144. (Cites 25–35% SL among yeast paralogs)
11. Cusack S.A. et al. (2021) Predictive models of genetic redundancy in *Arabidopsis thaliana*. *Mol Biol Evol* 38:3397. doi:10.1093/molbev/msab111. Data: doi:10.5281/zenodo.3987384. Code: https://github.com/ShiuLab/Manuscript_Code/tree/master/2021_Arabidopsis_redundancy_modeling
12. Lloyd J., Meinke D. (2012) A comprehensive dataset of genes with a loss-of-function mutant phenotype in Arabidopsis. *Plant Physiol* 158:1115. (PMC3291275)
13. Bolle C. et al. (2013) GABI-DUPLO: a collection of double mutants to overcome genetic redundancy in *Arabidopsis thaliana*. *Plant J* 75:157. doi:10.1111/tpj.12197
14. Hu Y. et al. (2023) Multi-Knock: a multi-targeted genome-scale CRISPR toolbox to overcome functional redundancy in plants. *Nat Plants* 9:572.
15. Ryan C.J. et al. (2012) Hierarchical modularity and the evolution of genetic interactomes across species. *Mol Cell* 46:691. doi:10.1016/j.molcel.2012.05.028
16. Roguev A. et al. (2008) Conservation and rewiring of functional modules revealed by an epistasis map in fission yeast. *Science* 322:405. Dixon S.J. et al. (2008) *PNAS* 105:16653.
17. Tischler J., Lehner B., Chen N., Fraser A.G. (2006) Combinatorial RNA interference in *C. elegans* reveals that redundancy between gene duplicates can be maintained for more than 80 million years. *Genome Biol* 7:R69. doi:10.1186/gb-2006-7-8-r69
18. Lehner B. et al. (2006) Systematic mapping of genetic interactions in *C. elegans* identifies common modifiers of diverse signaling pathways. *Nat Genet* 38:896.
19. Fischer B. et al. (2015) A map of directional genetic interactions in a metazoan cell. *eLife* 4:e05464. Heigwer F. et al. (2018) Time-resolved mapping of genetic interactions to model rewiring of signaling pathways. *eLife* 7:e40174.
20. Shapiro R.S. et al. (2018) A CRISPR–Cas9-based gene drive platform for genetic interaction analysis in *Candida albicans*. *Nat Microbiol* 3:73.
21. FlyBase release FB2026_03 precomputed files: https://s3ftp.flybase.org/releases/current/precomputed_files/
22. MGI reports: https://www.informatics.jax.org/downloads/reports/
23. UniProt Consortium. Reference proteomes release 2026_03: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/
