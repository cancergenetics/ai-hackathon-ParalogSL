# Literature review: methods and features for paralog SL / redundancy prediction, and evaluation pitfalls

Scope: what is known about predicting synthetic lethality (SL) and genetic redundancy between paralogs, which features are computable from protein sequence alone, and how to evaluate fairly. Searched up to Oct 2026. DOIs were checked against Crossref or Europe PMC unless marked **[unverified]**. "Seq-only" means the feature can be computed from a proteome FASTA, plus homology search against other proteomes, PLMs, or structure prediction. It excludes expression, PPI, essentiality, GI, GO, Pfam and other curated data.

Abbreviations: AP = average precision (≈AUPRC); prev = positive prevalence, which is the chance level of AP; WGD = whole-genome duplication; SSD = small-scale duplication; RF = random forest.

---

## 1. Paralog SL prediction in human cancer cell lines

### 1.1 De Kegel & Ryan 2019, PLoS Genet — doi:10.1371/journal.pgen.1008466
- Data: DepMap essentiality in 558 cell lines. This is a descriptive study, not a classifier.
- Paralogs are less often essential than singletons. The effect is stronger for genes with **several paralogs** or with **highly similar paralogs**.
- **WGD paralogs** have more variable essentiality than SSD paralogs. WGD pairs and pairs in protein complexes are more likely to be SL.
- Paralog buffering explains about 13–17% of variable essentiality.

### 1.2 De Kegel et al. 2021, Cell Systems — doi:10.1016/j.cels.2021.08.006 (bioRxiv 10.1101/2020.12.16.423022; numbers below are from the preprint v2 full text)
- **Paralog universe:** Ensembl 93 protein-coding paralog pairs with HGNC "gene with protein product", ≥20% identity in at least one direction. This gives ~36.6k pairs. Training-set filters removed families with >20 genes, tissue-specific pairs, and pairs not in DepMap. The authors note that **different Ensembl releases give different pair sets**: Ensembl 93 overlapped the screened pairs better than Ensembl 102.
- **Labels:** In DepMap (762 lines), a pair is SL if loss of A2 (deletion, LoF mutation or silencing) is associated with higher dependency on A1 in a regression with lineage covariate, at 10% FDR. This gave 126 SL pairs and 3,508 non-SL pairs (p>0.05); ambiguous pairs were dropped. Total 3,634 pairs, **prev ≈ 3.5%**. In independent Behan et al. screens, 47% of the SL pairs replicated and none of the non-SL pairs did.
- **Features (22), with our seq-only tag:**
  - Sequence identity — yes.
  - WGD vs SSD from Makino & McLysaght 2010 / OHNOLOGS v2. These use synteny, so not directly; only via proxies (see §5).
  - Family size — yes.
  - Closest paralog (each other's most similar paralog) — yes.
  - Conservation score: number of 9 species with an ortholog. Yes, via homology search.
  - Gene age from ProteinHistorian / phylostratigraphy — yes, via a homology-search panel.
  - Protein length ratio — yes.
  - Has an S. cerevisiae / S. pombe ortholog — yes.
  - Essential ortholog — no.
  - Shared Pfam domains — no (annotation).
  - Shared PPI, number of PPIs, PPI mean essentiality, physical interaction between A1 and A2, complex membership, complex essentiality — no.
  - Co-localisation — no.
  - GTEx mean expression and co-expression — no.
- **Single-feature results (training data, 4-fold CV):** Identity AUROC ~0.73, AP ~0.10. Shared PPI is about the same. Conservation score and gene age are "very predictive". SL pairs come from smaller families (median 5 vs 8) and are more often closest pairs. Having any S. pombe ortholog gave AP 6.2%; having an *essential* one gave AP 5.3%.
- **Model:** scikit-learn RF with 600 trees, max_depth 3, max_features 0.5. It was regularised heavily on purpose.
- **Splits:** stratified random 4-fold CV on pairs, repeated. There is **no gene or family grouping**.
- **Performance:**
  - CV: RF AUROC ~0.80, AP ~0.23, against identity at 0.73 and 0.10.
  - External test on the consensus of four combinatorial screens: 50 positives (SL in ≥2 screens), 407 negatives, prev ~11%. RF AP ~0.49 vs identity AP ~0.20.
  - Removing pairs seen in training gave AP ~0.46.
  - Hits in ≥1 screen: AP 0.38 vs identity 0.12.
- **Relevance to us:** Most of the lift over identity comes from PPI, complex and expression features, which we cannot use. The seq-only features that rank highly are identity, conservation/age, family size and closest-paralog status. A seq-only subset of this model has not been reported anywhere I found.

### 1.3 Dandage & Landry 2019, Mol Syst Biol — doi:10.15252/msb.20198871
- Data: CRISPR robustness landscape over >450 cell lines.
- **Heteromeric paralogs** (paralogs that physically interact) buffer each other less. They are more deleterious on loss because they have higher expression, more PPI partners, and stricter dosage balance.
- This links to Diss et al. 2017 (§2): physically dependent paralogs. For our purposes it means "high identity" does not imply SL when the paralogs form an obligate heteromer.
- Physical interaction between paralogs is measured data, but it could in principle be *predicted* from sequence or structure (AF-Multimer, PLM-based PPI). That is expensive and unvalidated for this purpose.

### 1.4 Flister et al. 2025, Cell Reports — doi:10.1016/j.celrep.2025.116512 (bioRxiv 10.1101/2024.07.16.603642)
- **Screen:** enAsCas12a digenic screen of **all 36,648 human paralog pairs** in 2 lines (NCI-H1299, MDA-MB-231), followed by a hybrid screen of ~5.1–5.3k pairs in 20 lines.
- **Genome-wide SL rate is very low.** The abstract reports <0.5% of all pairs as SL (~0.4% in the Discussion). The hybrid subset, which was enriched for likely pairs, had a 14.6% hit rate. **This is the best available estimate of prevalence for whole-proteome scoring.** AUPRC on all de novo pairs will have a chance level well below 1%.
- **Model:** weighted RF on 223 positives against random negatives, averaged over 10,000 models. Features are De Kegel's 23 plus DepMap correlation, TCGA co-expression and AlphaFold RMSD. Top features were PPI extent and essentiality of the shared PPI network, followed by expression. **No held-out AUROC/AUPRC and no identity baseline were reported** in the preprint I read. The published version may differ **[not checked in full]**.

### 1.5 Dennler & Ryan 2025, NAR Genom Bioinform — doi:10.1093/nargab/lqaf051 — **the closest prior work to our setting**
- **Question:** do structure and PLM similarity predict shared paralog function beyond identity? Labels were shared PPI and SL, in human and yeast.
- **Paralogs:** Ensembl 111, giving 107,103 human pairs after filtering (124,813 raw) and 5,673 yeast pairs.
- **SL labels:**
  - Human: De Kegel pipeline, plus a lenient variant. Prev ≈ 0.5% of all human pairs.
  - Yeast: SGA ε ≤ −0.35 and P ≤ 0.05 gives prev ≈ 4.8%. A lenient "neg. GI" set used ε ≤ −0.08.
  - **Only screened pairs are used**; unscreened pairs are dropped rather than treated as negatives.
- **Features (36, all seq-only):**
  - 9 Foldseek features on AlphaFold structures: fident, 3 TM-scores, LDDT, E-value, bits, prob, alignment length.
  - 20 PLM distances: ESM2-15B and ProtT5, using cosine, Euclidean, Manhattan and TS-SS distances.
  - 6 "similarity-search" features from MMseqs2 vs Swiss-Prot and Foldseek vs AFDB. **"Paralogs"** is the number of same-species proteins ranked above the partner. **"Species"** is the number of species with hits ranked above the partner. Both take the minimum over the two directions.
  - Spearman correlations with identity in human: fident 0.86, LDDT 0.78, TM 0.74, mean PLM distance 0.62, "species" 0.58/0.68, "paralogs" 0.18/0.13. **The "paralogs" and "species" rank features are almost orthogonal to identity.**
- **Model:** XGBoost on all 36 features, compared with identity. Splits were 4-fold CV and a **family-split 4-fold CV**, where each family sits wholly in train or test. The text does not say how families were defined.
- **Results (numbers are in figures only):**
  - The combined model beats identity on all datasets.
  - **Exception: on human SL, identity "performed exceptionally well"**, and the structure-only and PLM-only classifiers only *matched* it.
  - PLM features helped most in yeast.
  - Under the family split, performance drops but feature integration still helps. The similarity-search classifier and the all-feature classifier were the most robust.
- **Implications for us:**
  1. Beating identity on human SL with seq-only features is **hard but has been done** with the combined model.
  2. The rank and context features ("paralogs", "species") carry signal not in identity and are cheap.
  3. ESM2-15B is out of reach on our hardware; we must test whether ESM2-650M or 3B, or ProtT5, retain the signal.
  4. Features and predictions are on Zenodo (records 14973633 and 14975581). These are useful for sanity-checking our pipeline, **not as inputs**.

### 1.6 Kebabci, Ajmal, Adams & Ryan 2026, bioRxiv — doi:10.64898/2026.01.19.700065
- **Task:** cell-line-specific paralog SL, trained on 4,170 pairs × 10 lines from Ito, Parrish and Klingbeil.
- **Model:** RF using cell-line features (expression, essentiality, copy number, mutations) plus STRING/BioGRID and GO features. Not transferable to us.
- **Evaluation design is the useful part.** Four CV schemes: random; grouped by cell line; grouped by pair; and unseen pair × unseen line. AUROC falls from 0.92–0.93 (random) to 0.87 (strategy IV). On external Klingbeil data, AP was 0.32–0.33 for seen pairs and 0.22–0.24 for unseen pairs, with chance at ~2%.
- **Identity baseline in this context-specific task:** AUROC 0.54, AP ~4–6%. The context-agnostic classifier reached AUROC 0.68–0.73, AP 13–15%.
- **Lesson:** grouping by pair ID still lets a family leak across folds. The authors did not do family-level splits.

### 1.7 Other 2022–2026 human work (brief)
- **Combinatorial screens used as labels:**
  - Dede et al. 2020 — doi:10.1186/s13059-020-02173-2. ~400 pairs, 3 lines, 24 SL; 79% replicate in ≥2 lines.
  - Thompson et al. 2021 — doi:10.1038/s41467-021-21478-9. 1,191 pairs, 105 hits.
  - Parrish et al. 2021 — doi:10.1016/j.celrep.2021.109597. pgPEN, >2,000 paralogs; 12% SL in ≥1 context.
  - Gonatopoulos-Pournatzis et al. 2020 — doi:10.1038/s41587-020-0437-z. CHyMErA.
  - Ito et al. 2021 — doi:10.1038/s41588-021-00967-z. 5,065 pairs, 815 families, 11 lines.
  - Köferle et al. 2022 — doi:10.1016/j.celrep.2022.110636.
  - Esmaeili Anvar et al. 2024 — doi:10.1038/s41467-024-47795-3. in4mer/Inzolia, ~4,000 paralog pairs. The authors also identify "background-independent" paralog SLs across prior screens.
- **Aggregator:** SLKB, Gökbağ et al. 2024 — doi:10.1093/nar/gkad806. 11 CDKO experiments, 22 lines, 16,059 SL and 264,424 non-SL pairs. **The five scoring methods share only 1.21% of their top-10% SL pairs.**
- **Scoring benchmarks:** Ajmal et al. 2025 bioRxiv — doi:10.1101/2025.03.31.645224 — found "no single method performs best across all screens". Chou et al. 2024 bioRxiv — doi:10.1101/2024.08.19.608665 — found ZdLFC and RdLFC more consistent across lines than raw dLFC.
- **Review:** Ryan et al. 2023, Trends Cancer — doi:10.1016/j.trecan.2023.02.002.
- **Structure / PLM-based general SL models:**
  - Struct2SL 2025, CSBJ — doi:10.1016/j.csbj.2025.04.012. AlphaFold2 + sequence + PPI graph.
  - ELISL 2024, Bioinformatics — doi:10.1093/bioinformatics/btad764. Protein-sequence embeddings + omics. It explicitly targets **selection bias** in SL data and says sequence similarity is "underexplored".
  - MiT4SL 2025 bioRxiv — doi:10.1101/2025.04.20.649694.
  - Cilantro-SL 2026 bioRxiv — doi:10.64898/2026.02.25.708096. Single-cell foundation model on expression.
  - All of these rely on PPI or omics, and none is paralog-specific. "ESM4SL" is cited by Cilantro-SL, but **I could not find or verify the paper**.
- **No paper I found uses PLM or structure features for paralog SL with a sequence-only constraint and a cross-species evaluation.** That niche looks open.

---

## 2. Yeast (and other species): determinants of duplicate redundancy and SL

### 2.1 Core yeast findings
- **Gu et al. 2003, Nature — doi:10.1038/nature01198.** Duplicates compensate more often than singletons. Compensation frequency **correlates strongly with sequence similarity**. Deleting the more highly expressed copy is more severe.
- **Guan, Dunham & Troyanskaya 2007, Genetics — doi:10.1534/genetics.106.064329.** WGD paralogs share more interactors and functions than SSD paralogs, independent of identity. WGD paralogs have **a higher SL rate** and fewer essential genes.
- **Hakes et al. 2007, Genome Biol — doi:10.1186/gb-2007-8-10-r209.** WGD duplicates are more dispensable than SSD duplicates, but only for complex members. This is consistent with dosage balance.
- **Ihmels et al. 2007, MSB — doi:10.1038/msb4100127.** Duplicates contribute ~25% of robustness to null mutations. Even buffering pairs have largely non-overlapping deletion phenotypes.
- **DeLuna et al. 2008, Nat Genet — doi:10.1038/ng.123.** Competition fitness of single and double deletions of WGD metabolic paralogs. Functional overlap explains single-gene dispensability.
- **Dean et al. 2008, PLoS Genet — doi:10.1371/journal.pgen.1000113.** Redundancy is pervasive and persists ~100 My, for both WGD and SSD pairs.
- **Musso et al. 2008, Genome Res — doi:10.1101/gr.076174.108.** 140 of 399 WGD pairs are epistatic in standard conditions, and more are in stress media. **Labels are condition-dependent.**
- **VanderSluis et al. 2010, MSB — doi:10.1038/msb.2010.82.** Duplicates have fewer GIs. GI degree is **asymmetric** between copies (one copy diverges faster), and the asymmetry predicts other divergence properties.
- **Li, Yuan & Zhang 2010, PLoS Genet — doi:10.1371/journal.pgen.1001187.** ~500 WGD and SSD pairs. Backup capacity **decays roughly in proportion to divergence time**. Selected, long-retained redundant pairs share complexes and partners and evolve slowly.
- **Byrne & Wolfe 2007, Genetics — doi:10.1534/genetics.106.066951.** 56% of yeast ohnolog pairs show significant **rate asymmetry** against a non-WGD outgroup.
- **Kellis et al. 2004, Nature — doi:10.1038/nature02424.** 95% of cases of accelerated evolution involve only one copy.
- **Diss et al. 2017, Science — doi:10.1126/science.aai7685.** In 22 of 56 pairs one paralog compensates for lost PPIs. In about as many pairs the paralogs are *dependent*: they interact physically, regulate each other's abundance, and descend from **self-interacting ancestors**. Such pairs give fragility rather than buffering.
- **Marchant et al. 2019, eLife — doi:10.7554/elife.46754.** Heteromerisation is frequent among duplicated homomers and correlates with functional similarity. Structural pleiotropy slows divergence.
- **Soria, McGary & Rokas 2014, MBE — doi:10.1093/molbev/msu050.** Combines 8 functional-divergence measures over 499 ohnolog pairs. The **sequence-level measures** were DIVERGE type I/II rate-shift coefficients and SUPERFAMILY domain architecture. GI-profile divergence was found in 96% of pairs; redundancy was rare in the strict sense.
- **Kuzmin et al. 2020, Science — doi:10.1126/science.aaz5667.**
  - ~30% of duplicate pairs show a negative digenic interaction, against ~3.6% of random pairs. 240 dispensable WGD pairs were screened (44% of the 551 unique WGD pairs).
  - The trigenic fraction, a measure of retained redundancy, is negatively correlated with digenic-degree asymmetry and with sequence divergence rate (r ≈ −0.22).
  - "Evolutionarily constrained" pairs have correlated position-specific rate patterns. This is computable from an alignment with a pre-WGD outgroup and is **seq-only**.
  - These features belong to a model of "structural and functional entanglement".
- **Costanzo et al. 2016, Science — doi:10.1126/science.aaf1420.** SGA global network with ~550k negative interactions. Essential genes are dense GI hubs. Thresholds are commonly ε<−0.08 (lenient) or ε<−0.12 (stringent), at p<0.05.

**Seq-computable determinants that recur across studies:**
1. Identity or divergence (strong, monotone).
2. Duplication age; WGD > SSD in yeast.
3. Rate asymmetry against an outgroup, which signals divergence.
4. Family size and closest-paralog status.
5. Conservation breadth.
6. Ancestral self-interaction or heteromer formation (structural; reduces buffering).

Expression and PPI are also strong but are disallowed as features.

### 2.2 Predictive models of redundancy outside human cancer
- **Chen, Bandyopadhyay, Shasha & Birnbaum 2010, BMC Evol Biol — doi:10.1186/1471-2148-10-357.** Arabidopsis; SVM on multiple attributes is ~2× more precise than BLAST E-value or co-expression alone. Many predicted redundant pairs are old (Ks>1).
- **Cusack et al. 2021, MBE (Shiu lab) — doi:10.1093/molbev/msab111.**
  - Arabidopsis; 300 benchmark pairs labelled from double-mutant phenotypes, with 9 redundancy definitions tested.
  - ~4,100 features: GO/AraCyc/Pfam, expression, epigenetic, network, **evolutionary** (Ka, Ks, Ka/Ks, identity, duplication mode, reciprocal best match, retention across 15 species, family size) and **protein properties** (length, pI).
  - Best results: RD4 CV AUROC 0.84 / AUPRC 0.82; holdout 0.73 / 0.62, with roughly balanced classes. Top features: **α-WGD origin**, **reciprocal best match**, TF annotation, stress expression.
  - The definition of redundancy, the transformations and the algorithm all strongly affected performance.
  - **No cross-species transfer was tested.** Test pairs may share genes with training pairs; the paper does not describe handling this.
- **Tischler et al. 2006, Genome Biol — doi:10.1186/gb-2006-7-8-r69.** C. elegans combinatorial RNAi: 16 of ~143 duplicate pairs redundant, 14 of them older than 80 My.
- **Insects and plants beyond Arabidopsis:** I found **no systematic, labelled paralog-SL dataset or predictive model** for Drosophila, other insects, or crop plants. Unseen-species generalisation can therefore only be tested by leave-one-species-out between yeast and human, plus small literature or curation sets.

### 2.3 Cross-species conservation of GIs / paralog SL
- Dixon et al. 2008, PNAS — doi:10.1073/pnas.0806261105: ~29% of GIs are shared between S. cerevisiae and S. pombe in a directly compared gene set.
- Tischler, Lehner & Fraser 2008, Nat Genet — doi:10.1038/ng.114: yeast GIs are **poorly conserved** in worms.
- So **direct label transfer by orthology is weak**. What should transfer are *determinants* (divergence, age, family context, asymmetry), not specific pairs. I found no study that maps yeast paralog-SL pairs to human paralog pairs and measures concordance **[gap; not found]**.

---

## 3. General GI/SL prediction from sequence (transferable parts only)
- Most SL predictors (SLMGAE, KG4SL, GCATSL, PTGNN and others; see Feng et al. 2024) are **graph or knowledge-graph methods built on known SL edges, PPI and GO**. They are not usable in our setting, and their gains come mostly from memorising each gene's SL-degree.
- Pair-input lessons from sequence-based PPI prediction transfer directly:
  - **Bernett, Blumenthal & List 2024, Brief Bioinform — doi:10.1093/bib/bbae076.** With random splits, deep sequence-PPI models learn **sequence similarity and node degree shortcuts**. Once leakage is removed by similarity-aware splits, performance becomes **near-random**. Simple similarity and degree baselines are competitive.
  - PLM-interact (Nat Commun 2025, doi:10.1038/s41467-025-64512-w **[DOI from search result, not Crossref-checked]**) fine-tunes ESM-2 on protein pairs for physical PPI. It could serve as a seq-only proxy for "paralogs physically interact" (Diss/Dandage heteromer effect). That is speculative and costly.
- **PLM embeddings for homology and function similarity (seq-only tools):**
  - ESM-2 / ESMFold, Lin et al. 2023 — doi:10.1126/science.ade2574.
  - ProtT5, Elnaggar et al. 2022 — doi:10.1109/TPAMI.2021.3095381.
  - ProstT5 (sequence→3Di, no folding needed), Heinzinger et al. 2024 — doi:10.1093/nargab/lqae150.
  - Foldseek, van Kempen et al. 2024 — doi:10.1038/s41587-023-01773-0.
  - pLM-BLAST, Kaminski et al. 2023 — doi:10.1093/bioinformatics/btad579.
  - TM-Vec, Hamamsy et al. 2024 — doi:10.1038/s41587-023-01917-2.
  - Search and clustering: MMseqs2 — doi:10.1038/nbt.3988; DIAMOND — doi:10.1038/s41592-021-01101-x.
  - Orthology and gene trees: OrthoFinder — doi:10.1186/s13059-019-1832-y.
- **TripleThreat (ICLR 2026; Misra, Nayar, Lundberg, Altman) [abstract only; no DOI]:** pLM embeddings separate functionally divergent paralogs from conserved orthologs *worse as identity, length and species confounders are controlled*. Alignment-based pLMs encode species more than function. This is a warning that PLM-similarity features may be largely identity in disguise, so **always test them conditionally on identity**.

---

## 4. Evaluation pitfalls

### 4.1 Pair-input leakage (gene, family, cross-species)
- **Park & Marcotte 2012, Nat Methods — doi:10.1038/nmeth.2259.** Pair-input CV is inflated because test pairs share components with training pairs. Report C1 (both genes seen), C2 (one seen) and C3 (neither seen) separately. C3 is closest to real-world use.
- **Hamp & Rost 2015, Bioinformatics — doi:10.1093/bioinformatics/btu857.** In C3, test proteins must also be *non-similar* to training proteins, both positives and negatives. The apparent best method depends on test–train similarity and on the assumed negative:positive ratio.
- **Feng et al. 2024, Nat Commun — doi:10.1038/s41467-024-52900-7** (author correction doi:10.1038/s41467-025-56368-x).
  - Benchmarks 12 SL methods under CV1 (random pairs), CV2 (one gene unseen) and CV3 (both unseen).
  - F1 of the best method falls from >0.877 (CV1) to 0.779 (CV2) to ~0.73 (CV3); others fall below 0.70. NDCG@10 collapses from ~0.27 to <0.01.
  - Performance also degrades as the negative ratio grows (1:1 → 1:50).
- **Paralog-specific leakage.** Pairs from the same family share genes and near-identical sequences; for example A–B in train and A–C or B–C in test, which forms a triangle. *Random pair splits leak family membership.* De Kegel 2021 used random stratified CV. Kebabci 2026 grouped by pair only. Dennler 2025 showed that a family split lowers scores.
- **Cross-species leakage (specific to us).** Human and yeast families are often homologous. A model trained on a yeast family and tested on its human orthologous family is not an "unseen-species/unseen-family" test. Families must be grouped **across species** using shared homology clusters or orthogroups.
- General taxonomy: Kapoor & Narayanan 2023, Patterns — doi:10.1016/j.patter.2023.100804.

### 4.2 Degree / hub bias
- In SL data, a few genes take part in many positive pairs. Examples: big complexes, recurrent drug targets, and genes from families of 3–10 where all pairs are SL.
- Graph SL models exploit per-gene positive rates. Bernett et al. 2024 showed explicitly that sequence-PPI models exploit node degree.
- Feng et al. did not analyse degree bias; I found no SL-specific quantification **[gap]**.
- **Mitigations:**
  1. Report a "gene-prior" baseline: the mean training label of each gene, or the gene's count of positives. If it is competitive, the split leaks.
  2. Group splits by family.
  3. Optionally cap or weight pairs per gene in training and in metrics.
  4. Report performance on C3-type pairs.
  5. Use a per-family bootstrap for CIs.

### 4.3 Label noise, ascertainment, context
- **Screen agreement is low** (De Kegel 2021): hit rates are ~6–15% per screen. 68–86% of hits in Thompson, Parrish and CHyMErA occur in only one cell line. Pairwise between-screen hit overlap averages ~31%. **In RPE-1, which two studies screened, only 2 of 26 hits agreed (~8%).**
- SLKB: the 5 scoring methods share 1.21% of their top-10% pairs. Ajmal 2025: no scoring method is best everywhere.
- **Ascertainment bias in screened pairs** (De Kegel 2021 Fig. 1): compared with all paralog pairs, screened pairs have higher identity and smaller families (median 3 vs 11). They are 73% WGD vs 19%, and 89% closest pairs vs 13%. Consequences:
  1. Identity's range is restricted within screened pairs, so it is *less* discriminative there.
  2. Models trained on screens learn little about non-closest pairs in large families, which make up most of a whole proteome.
  3. Prevalence in screens (~11%) is far above the genome-wide rate (<0.5%; Flister 2025).
- **Yeast labels:**
  - Prevalence depends strongly on the threshold: ε≤−0.35 gives ~4.8% (Dennler), while lenient negative GI gives ~30% of duplicate pairs (Kuzmin).
  - The deletion SGA cannot test pairs that include an essential gene without TS alleles.
  - GIs are condition-dependent (Musso 2008).
  - Untested pairs must not be called negatives. Dennler dropped unscreened pairs.
- **Recommendation:** define positives as "robust" (≥2 screens or cell lines in human; stringent ε in yeast). Draw negatives only from tested non-hits; consider dropping ambiguous pairs, as De Kegel and Dennler did. Keep a lenient label set for sensitivity analysis only, and never tune labels on test results.

### 4.4 Class imbalance and AUPRC
- AP's chance level equals prevalence, so AP is not comparable across datasets with different prevalence (yeast vs human, screened set vs all pairs). Saito & Rehmsmeier 2015, PLoS One — doi:10.1371/journal.pone.0118432 — argue for PR curves under imbalance. McDermott et al. 2024 (arXiv:2401.06091) caution that AUPRC weights high-score regions and can favour subpopulations with higher prevalence. This is relevant when pooling species.
- **Report per species and per split:**
  - AP, prevalence and AP/prevalence (lift).
  - AP of the identity baseline on *the same pairs*.
  - ΔAP with a paired bootstrap CI, resampling families rather than pairs.
  - Optionally precision at fixed recall.
- Avoid pooling yeast and human into one AUPRC unless each species is reweighted to a common prevalence.

### 4.5 What is a fair "sequence identity only" baseline
- Identity has several definitions:
  - Ensembl gives a value in each direction. De Kegel used "≥20% in at least one direction"; Dennler used the minimum identity of the two directions.
  - Alternatives are global (Needleman–Wunsch) identity, local identity × coverage, and bitscore normalised by self-score.
- **A fair baseline:** compute identity with the *same* aligner and pair set as the tool. Choose the best variant (min, max or mean direction; local or global) **on training data only**, then fix it. Optionally fit a monotone calibration (isotonic or logistic) per species on training folds; this does not change AUPRC within a species.
- **A stronger "trivial" baseline to also report:** logistic regression on identity + family size + closest-paralog flag. De Kegel and Cusack show these are cheap and strong. Beating identity alone while losing to this baseline would mean our expensive features add nothing.
- The identity baseline must be evaluated on the **same candidate set**, including pairs our caller finds that the label sets lack, and with the same handling of labelled pairs our caller misses.

---

## 5. How paralog pairs are defined in label datasets, and implications for de novo calling

| Source | Paralog definition |
|---|---|
| De Kegel 2019/2021, Kebabci 2026 | Ensembl 93 Compara paralogs (gene-tree based), ≥20% identity in ≥1 direction; HGNC protein-coding. WGD labels from Makino & McLysaght 2010 (doi:10.1073/pnas.0914697107) / OHNOLOGS v2 (Singh & Isambert 2020, doi:10.1093/nar/gkz909), both synteny-based |
| Dennler & Ryan 2025 | Ensembl 111 (human 124,813 raw → 107,103; yeast 5,673 pairs) |
| Flister 2025 | "36,648 known paralog pairs" (≈ the Ensembl-derived De Kegel universe; exact definition not stated in the preprint) |
| Combinatorial screens | Hand-selected candidates (high identity, small families, WGD, closest pairs), typically drawn from Ensembl |
| Yeast WGD studies (Kuzmin, DeLuna, Musso, Byrne & Wolfe) | YGOB synteny-curated ohnologs (Byrne & Wolfe 2006, doi:10.1093/nar/gkj041); ~551 unique WGD pairs |
| Yeast SGA (Costanzo 2016) | No paralog definition; GI scores for all tested pairs; paralogs are added afterwards (Ensembl or YGOB) |

- **Ensembl Compara** (current documentation): HMM classification into families (the documentation names the TreeFam HMM library) with BLAST+/hcluster_sg for unclassified genes. Then MSA and a TreeBeST gene tree, reconciled to the species tree. **Within-species paralogs are pairs whose last common ancestor is a duplication node.** Clusters larger than 1,500 genes are split. Gene-tree methods: Vilella et al. 2009 — doi:10.1101/gr.073585.107. Ensembl gives identity relative to each protein, so values are asymmetric **[field semantics from Ensembl BioMart conventions; not re-verified this session]**. Paralogy therefore includes very old, low-identity pairs (many <20%), and pair sets change between releases.
- **Implications for our de novo caller (proteome-only):**
  1. **Recall against label universes:** report the fraction of labelled positive and negative pairs that our caller recovers at each threshold. A pair we never call scores the minimum, which hurts recall. Evaluate on (a) labelled pairs our caller finds and (b) all labelled pairs with misses scored at the minimum. Report both.
  2. Fix thresholds a priori and do not tune them on test labels. Ensembl-like sensitivity likely needs E-value ≤1e-5 or so plus a coverage requirement on the shorter protein. A ≥20% identity floor matches the De Kegel universe. The family graph should use a *looser* threshold than the pair list, so that family size and redundancy-context features see distant relatives.
  3. **WGD status cannot be called from a proteome**, since there are no coordinates or synteny. Proxies:
     - duplication age from species overlap in gene trees (OrthoFinder with a reference panel);
     - "1:2 orthology to an outgroup" patterns;
     - low family size (ohnologs resist SSD; Makino & McLysaght 2010).
     The proxy must be validated against YGOB/OHNOLOGS labels as an *evaluation-only* check; those labels must not become features.
  4. **No CDS means no Ks/Ka.** Use protein distances: identity, ML protein distance, or distance to an outgroup ortholog.
  5. Alternative isoforms are excluded by the one-protein-per-gene input; still guard against near-identical fragments and gene splits.
  6. **Large families:** Ensembl splits clusters; De Kegel dropped families >20. Our tool must score all pairs, so it needs explicit family-size features and must make tractable decisions about all-vs-all pair enumeration (n² within families).

---

## 6. Implications for our design

### 6.1 Prioritised feature ideas

Cost is for a full human proteome (~20k proteins) on an M4 with 24 GB, rough estimate.

| # | Feature group | Seq-only? | Est. cost | Evidence / rationale |
|---|---|---|---|---|
| 1 | **Pairwise similarity suite**: identity (min/max/mean direction; local and global), query/target coverage, bitscore / self-bitscore, E-value, length ratio, aligned-length fraction | yes | low (MMseqs2/DIAMOND all-vs-all: minutes) | Identity is the strongest single feature (AUROC ~0.73 in De Kegel; "exceptionally well" for human SL in Dennler); Gu 2003 |
| 2 | **Family context**: family size at 2–3 thresholds, closest-paralog / reciprocal-best-hit flag, rank of partner among each gene's paralogs ("paralogs" feature, min over directions), number of third paralogs similar to both (triangle redundancy), max identity to any *other* paralog | yes | low | De Kegel (family size, closest pair); Dennler (rank features ⟂ identity, ρ≈0.15); Cusack (reciprocal best match top feature) |
| 3 | **Conservation / age vs a fixed reference proteome panel** (~20–50 eukaryotes spanning fungi, metazoa incl. insects, plants, protists): per-gene phyletic breadth, per-pair conservation score, phylostratum, duplication-node depth / "species" rank feature | yes | moderate (MMseqs2 20k × ~0.5–1.5M seqs; tens of minutes to ~2 h; panel ~1–3 GB) | De Kegel (conservation score and gene age "very predictive"); Dennler ("species" feature); Li 2010 (backup decays with divergence time) |
| 4 | **Divergence asymmetry vs outgroup ortholog**: |d(A,O) − d(B,O)| / (d(A,O)+d(B,O)), using best panel hit O outside the duplication | yes | moderate (needs ortholog assignment; reuses #3 search) | VanderSluis 2010; Byrne & Wolfe 2007 (56% asymmetric); Kellis 2004; Kuzmin 2020 (asymmetry anti-correlates with redundancy) |
| 5 | **Per-gene "importance" proxies from sequence**: identity to most distant panel ortholog (slow evolution ≈ essential), phyletic breadth, protein length, low-complexity/disorder fraction (sequence-only predictor) | yes | low–moderate | De Kegel (any ortholog in S. pombe more informative than essential ortholog); conservation ≈ essentiality |
| 6 | **PLM embedding similarity** (ESM2-650M first; 3B if time allows; ProtT5-XL-half as alternative): cosine/Euclidean of mean-pooled embeddings, plus per-gene PCA-reduced embeddings | yes | moderate–high (ESM2-650M on 20k proteins on MPS/CPU: ~1–3 h est.; 3B ~4–5× more; truncate >1–2k aa) | Dennler (PLM adds beyond identity, mostly yeast; used 15B); TripleThreat warns of confounding. Per-gene embeddings risk memorising genes, so test only with family-grouped CV |
| 7 | **Structural similarity without folding**: ProstT5 → 3Di, then Foldseek 3Di+AA alignment (TM-score proxy, LDDT-like, fident) | yes | moderate (ProstT5 ~ same order as #6) | Dennler (Foldseek features; LDDT strong for yeast PPI). ESMFold for whole proteomes is likely too slow on CPU; AFDB is unavailable for new species and bulk download is disallowed |
| 8 | **Position-specific rate correlation** between paralogs (from MSA with outgroup orthologs) | yes | high (per-family MSAs + rate inference) | Kuzmin 2020 ("evolutionarily constrained" pairs); exploratory |
| 9 | **Predicted physical interaction / heteromer between paralogs** (PLM-interact-style or AF-Multimer on top candidates only) | yes (predicted) | high/very high | Diss 2017, Dandage 2019, Marchant 2019: dependent heteromeric paralogs buffer less. Low priority; only on a shortlist if cheap |
| 10 | **WGD-like proxy**: duplication age bin + 1:2 orthology pattern to outgroup + low family size | yes | low given #3 | WGD > SSD for SL (Guan 2007; De Kegel 2019/2021); validate only against YGOB/OHNOLOGS |
| — | Expression, PPI, complex membership, essentiality, GO, Pfam, ortholog essentiality, label transfer from orthologous pairs | **no** | — | Excluded by the project rules. Orthologous-pair label transfer is also weak (Tischler 2008; Dixon 2008), and the fact that a pair is SL in another species is measured data |

Gray area to flag to the PI: using measured PPI or shared-PPI as an **auxiliary training target**, never as an input feature. An example is multi-task learning that predicts "shared PPI" from sequence features. Inference stays seq-only, but this needs an explicit decision against PROMPT.md.

**Suggested order:** 1 → 2 → 3 → 10 → 4 → 5 → 7 → 6 → 8/9. Groups 1–3 are cheap and are the most likely route to beating identity under family-grouped CV. Groups 6–7 are where the literature reports gains beyond identity (Dennler), but they need ablation conditional on identity.

**Model:** a gradient-boosted tree or regularised RF, with monotone constraint on identity as an option, trained jointly on yeast and human with a species-agnostic feature set. No species ID feature: the tool must work on unseen species. Normalise family-size and rank features within each proteome (as percentiles), because human families are larger than yeast families (median 4 vs smaller in yeast; De Kegel).

### 6.2 Evaluation recommendations
1. **Primary split: family-grouped and cross-species-grouped.** Build homology clusters over the union of the yeast and human proteomes (MMseqs2 at a loose threshold, or OrthoFinder orthogroups). Assign each cluster wholly to one fold so that no test pair shares a gene *or* a homologous family with training. Use repeated grouped K-fold with fixed seeds.
2. **Secondary splits:**
   - (a) Leave-one-species-out: train on yeast and test on human, and the reverse. This is the proxy for unseen plants and insects.
   - (b) Park–Marcotte C1/C2/C3 breakdown.
   - (c) Random pair split, reported only to show how much leakage inflates scores.
3. **Report per species:** AP, prevalence, lift, and AUROC as secondary. Give identity-baseline AP and the trivial-baseline AP (identity + family size + closest flag) on identical pairs, ΔAP with family-level bootstrap CIs, and a **gene-prior leakage baseline**.
4. **Label sets:**
   - Human: robust positives (≥2 screens/lines, or the De Kegel DepMap set) and negatives from tested non-hits. Keep screen-consensus pairs as a held-out external test, as in De Kegel 2021.
   - Yeast: SGA at a pre-registered stringent threshold, tested pairs only, with a lenient threshold as sensitivity analysis.
   - Document how each candidate set is biased (identity, family size, WGD and closest fractions against the full de novo pair set).
5. **Whole-proteome realism:** also evaluate with all de novo pairs as candidates, treating unlabelled pairs as unknown and excluding them. Report expected precision at genome-wide prevalence (<0.5% human; Flister 2025) by reweighting negatives. Report de novo caller recall of labelled pairs.
6. **Controls:**
   - Shuffled-label negative control, which should give AP ≈ prevalence.
   - A positive-control toy proteome with planted duplicates.
   - An identity-matched check that new features add signal *within identity strata*. This guards against PLM features that only proxy identity (TripleThreat).
7. **Freeze before testing:** thresholds for paralog calling, label definitions, and the baseline identity variant are chosen on training folds only and recorded before any test evaluation.

---

## References (DOIs verified via Crossref / Europe PMC unless marked)

- Ajmal H, Nandi S, Kebabci N, Ryan CJ (2025) Benchmarking genetic interaction scoring methods for identifying synthetic lethality from combinatorial CRISPR screens. bioRxiv. doi:10.1101/2025.03.31.645224
- Bernett J, Blumenthal DB, List M (2024) Cracking the black box of deep sequence-based protein–protein interaction prediction. Brief Bioinform. doi:10.1093/bib/bbae076
- Buchfink B, Reuter K, Drost HG (2021) Sensitive protein alignments at tree-of-life scale using DIAMOND. Nat Methods. doi:10.1038/s41592-021-01101-x
- Byrne KP, Wolfe KH (2006) Visualizing syntenic relationships among the hemiascomycetes with the Yeast Gene Order Browser. Nucleic Acids Res. doi:10.1093/nar/gkj041
- Byrne KP, Wolfe KH (2007) Consistent patterns of rate asymmetry and gene loss indicate widespread neofunctionalization of yeast genes after WGD. Genetics. doi:10.1534/genetics.106.066951
- Chen HW, Bandyopadhyay S, Shasha DE, Birnbaum KD (2010) Predicting genome-wide redundancy using machine learning. BMC Evol Biol. doi:10.1186/1471-2148-10-357
- Chou J, Esmaeili Anvar N, Elghaish R, Chen J, Hart T (2024) Optimal methods for analyzing targeted pairwise knockout screens. bioRxiv. doi:10.1101/2024.08.19.608665
- Costanzo M et al. (2016) A global genetic interaction network maps a wiring diagram of cellular function. Science. doi:10.1126/science.aaf1420
- Cusack SA et al. (2021) Predictive models of genetic redundancy in Arabidopsis thaliana. Mol Biol Evol. doi:10.1093/molbev/msab111
- Dandage R, Landry CR (2019) Paralog dependency indirectly affects the robustness of human cells. Mol Syst Biol. doi:10.15252/msb.20198871
- De Kegel B, Ryan CJ (2019) Paralog buffering contributes to the variable essentiality of genes in cancer cell lines. PLoS Genet. doi:10.1371/journal.pgen.1008466
- De Kegel B, Quinn N, Thompson NA, Adams DJ, Ryan CJ (2021) Comprehensive prediction of robust synthetic lethality between paralog pairs in cancer cell lines. Cell Syst. doi:10.1016/j.cels.2021.08.006 (preprint doi:10.1101/2020.12.16.423022)
- Dean EJ, Davis JC, Davis RW, Petrov DA (2008) Pervasive and persistent redundancy among duplicated genes in yeast. PLoS Genet. doi:10.1371/journal.pgen.1000113
- Dede M, McLaughlin M, Kim E, Hart T (2020) Multiplex enCas12a screens detect functional buffering among paralogs otherwise masked in monogenic Cas9 knockout screens. Genome Biol. doi:10.1186/s13059-020-02173-2
- DeLuna A et al. (2008) Exposing the fitness contribution of duplicated genes. Nat Genet. doi:10.1038/ng.123
- Dennler O, Ryan CJ (2025) Evaluating sequence and structural similarity metrics for predicting shared paralog functions. NAR Genom Bioinform 7(2):lqaf051. doi:10.1093/nargab/lqaf051
- Diss G et al. (2017) Gene duplication can impart fragility, not robustness, in the yeast protein interaction network. Science. doi:10.1126/science.aai7685
- Dixon SJ et al. (2008) Significant conservation of synthetic lethal genetic interaction networks between distantly related eukaryotes. PNAS. doi:10.1073/pnas.0806261105
- Elnaggar A et al. (2022) ProtTrans. IEEE TPAMI. doi:10.1109/TPAMI.2021.3095381
- Emms DM, Kelly S (2019) OrthoFinder. Genome Biol. doi:10.1186/s13059-019-1832-y
- Esmaeili Anvar N et al. (2024) Efficient gene knockout and genetic interaction screening using the in4mer CRISPR/Cas12a multiplex knockout platform. Nat Commun. doi:10.1038/s41467-024-47795-3
- Feng Y et al. (2024) Benchmarking machine learning methods for synthetic lethality prediction in cancer. Nat Commun. doi:10.1038/s41467-024-52900-7 (correction doi:10.1038/s41467-025-56368-x)
- Flister MJ et al. (2025) Genetic mapping and predictive modeling of paralog synthetic lethality. Cell Rep. doi:10.1016/j.celrep.2025.116512 (preprint doi:10.1101/2024.07.16.603642)
- Gökbağ B et al. (2024) SLKB: synthetic lethality knowledge base. Nucleic Acids Res. doi:10.1093/nar/gkad806
- Gonatopoulos-Pournatzis T et al. (2020) Genetic interaction mapping and exon-resolution functional genomics with a hybrid Cas9–Cas12a platform. Nat Biotechnol. doi:10.1038/s41587-020-0437-z
- Gu Z et al. (2003) Role of duplicate genes in genetic robustness against null mutations. Nature. doi:10.1038/nature01198
- Guan Y, Dunham MJ, Troyanskaya OG (2007) Functional analysis of gene duplications in Saccharomyces cerevisiae. Genetics. doi:10.1534/genetics.106.064329
- Hakes L et al. (2007) All duplicates are not equal: the difference between small-scale and genome duplication. Genome Biol. doi:10.1186/gb-2007-8-10-r209
- Hamamsy T et al. (2024) Protein remote homology detection and structural alignment using deep learning (TM-Vec). Nat Biotechnol. doi:10.1038/s41587-023-01917-2
- Hamp T, Rost B (2015) More challenges for machine-learning protein interactions. Bioinformatics. doi:10.1093/bioinformatics/btu857
- Heinzinger M et al. (2024) Bilingual language model for protein sequence and structure (ProstT5). NAR Genom Bioinform. doi:10.1093/nargab/lqae150
- Hua K, Haber E, Ma J (2026) Uncertainty-aware synthetic lethality prediction with pretrained foundation models. bioRxiv. doi:10.64898/2026.02.25.708096
- Huang Y et al. (2025) Struct2SL. Comput Struct Biotechnol J. doi:10.1016/j.csbj.2025.04.012
- Ihmels J et al. (2007) Backup without redundancy: genetic interactions reveal the cost of duplicate gene loss. Mol Syst Biol. doi:10.1038/msb4100127
- Ito T et al. (2021) Paralog knockout profiling identifies DUSP4 and DUSP6 as a digenic dependence in MAPK pathway-driven cancers. Nat Genet. doi:10.1038/s41588-021-00967-z
- Kaminski K et al. (2023) pLM-BLAST. Bioinformatics. doi:10.1093/bioinformatics/btad579
- Kapoor S, Narayanan A (2023) Leakage and the reproducibility crisis in machine-learning-based science. Patterns. doi:10.1016/j.patter.2023.100804
- Kebabci N, Ajmal H, Adams DJ, Ryan CJ (2026) A predicted cancer dependency map for paralog pairs. bioRxiv. doi:10.64898/2026.01.19.700065
- Kellis M, Birren BW, Lander ES (2004) Proof and evolutionary analysis of ancient genome duplication in the yeast S. cerevisiae. Nature. doi:10.1038/nature02424
- Köferle A et al. (2022) Interrogation of cancer gene dependencies reveals paralog interactions of autosome and sex chromosome-encoded genes. Cell Rep. doi:10.1016/j.celrep.2022.110636
- Kuzmin E et al. (2020) Exploring whole-genome duplicate gene retention with complex genetic interaction analysis. Science. doi:10.1126/science.aaz5667
- Li J, Yuan Z, Zhang Z (2010) The cellular robustness by genetic redundancy in budding yeast. PLoS Genet. doi:10.1371/journal.pgen.1001187
- Lin Z et al. (2023) Evolutionary-scale prediction of atomic-level protein structure with a language model (ESM-2/ESMFold). Science. doi:10.1126/science.ade2574
- Makino T, McLysaght A (2010) Ohnologs in the human genome are dosage balanced and frequently associated with disease. PNAS. doi:10.1073/pnas.0914697107
- Marchant A et al. (2019) The role of structural pleiotropy and regulatory evolution in the retention of heteromers of paralogs. eLife. doi:10.7554/eLife.46754
- McDermott M et al. (2024) A closer look at AUROC and AUPRC under class imbalance. arXiv:2401.06091 (NeurIPS 2024)
- Misra M, Nayar G, Lundberg E, Altman R (2026) TripleThreat: benchmarking functional sensitivity in protein representations with paralog–ortholog triplets. ICLR 2026. **[abstract only seen; no DOI]**
- Musso G et al. (2008) The extensive and condition-dependent nature of epistasis among whole-genome duplicates in yeast. Genome Res. doi:10.1101/gr.076174.108
- Park Y, Marcotte EM (2012) Flaws in evaluation schemes for pair-input computational predictions. Nat Methods. doi:10.1038/nmeth.2259
- Parrish PCR et al. (2021) Discovery of synthetic lethal and tumor suppressor paralog pairs in the human genome. Cell Rep. doi:10.1016/j.celrep.2021.109597
- PLM-interact (2025) Extending protein language models to predict protein–protein interactions. Nat Commun. doi:10.1038/s41467-025-64512-w **[DOI from URL; not Crossref-checked]**
- Ryan CJ, Mehta I, Kebabci N, Adams DJ (2023) Targeting synthetic lethal paralogs in cancer. Trends Cancer. doi:10.1016/j.trecan.2023.02.002
- Saito T, Rehmsmeier M (2015) The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. PLoS One. doi:10.1371/journal.pone.0118432
- Singh PP, Isambert H (2020) OHNOLOGS v2. Nucleic Acids Res. doi:10.1093/nar/gkz909
- Soria PS, McGary KL, Rokas A (2014) Functional divergence for every paralog. Mol Biol Evol. doi:10.1093/molbev/msu050
- Steinegger M, Söding J (2017) MMseqs2. Nat Biotechnol. doi:10.1038/nbt.3988
- Tao et al. (2025) MiT4SL. bioRxiv. doi:10.1101/2025.04.20.649694
- Tepeli YI, Seale C, Gonçalves JP (2024) ELISL: early–late integrated synthetic lethality prediction in cancer. Bioinformatics. doi:10.1093/bioinformatics/btad764
- Thompson NA et al. (2021) Combinatorial CRISPR screen identifies fitness effects of gene paralogues. Nat Commun. doi:10.1038/s41467-021-21478-9
- Tischler J, Lehner B, Chen N, Fraser AG (2006) Combinatorial RNA interference in C. elegans reveals that redundancy between gene duplicates can be maintained for more than 80 million years. Genome Biol. doi:10.1186/gb-2006-7-8-r69
- Tischler J, Lehner B, Fraser AG (2008) Evolutionary plasticity of genetic interaction networks. Nat Genet. doi:10.1038/ng.114
- van Kempen M et al. (2024) Fast and accurate protein structure search with Foldseek. Nat Biotechnol. doi:10.1038/s41587-023-01773-0
- VanderSluis B et al. (2010) Genetic interactions reveal the evolutionary trajectories of duplicate genes. Mol Syst Biol. doi:10.1038/msb.2010.82
- Vilella AJ et al. (2009) EnsemblCompara GeneTrees. Genome Res. doi:10.1101/gr.073585.107

**Not verified / not found:**
- "ESM4SL" (cited by Cilantro-SL).
- Any study of yeast→human paralog-SL concordance by orthology.
- Any labelled paralog-SL dataset or model for insects or crop plants.
- Exact AUROC/AUPRC values in Dennler & Ryan 2025; they are in figures only and were not read.
- Held-out metrics in Flister et al. 2025 (published version not read in full).
