# LITERATURE — synthesis

Detailed reviews (with full reference lists and verification marks) are in `docs/lit/`:
- `methods.md` — prediction methods, features, evaluation pitfalls (≈60 refs, DOIs checked against Crossref/Europe PMC)
- `human_labels.md` — human combinatorial-screen label sources, licenses, labelling scheme
- `yeast_and_other_species_labels.md` — yeast labels, unseen-species test sets, reference proteomes
- `compute_feasibility.md` — measured runtimes/memory of candidate tools on this M4 (scripts in `scripts/bench/`)

This file condenses what matters for the design. Items marked [U] were not verified at the primary source by the survey agents.

## 1. Prior work on predicting paralog SL

| Work | Setting | Features | Sequence-only? | Result | Split |
|---|---|---|---|---|---|
| De Kegel & Ryan 2019 PLoS Genet (10.1371/journal.pgen.1008466) | Human, DepMap | paralog buffering of essentiality | — | Paralogs are less often essential; family size, identity matter | — |
| De Kegel et al. 2021 Cell Syst (10.1016/j.cels.2021.08.006) | Human, DepMap-inferred robust SL (126/3,634, 3.5%) | 22 features: identity, family size, closest pair, conservation, gene age, PPI, complexes, expression… | partly (identity, family size, closest-pair, conservation, age) | RF AUROC ≈0.80, AP ≈0.23; identity alone AUROC ≈0.73, AP ≈0.10 | random pair CV (optimistic) |
| Dandage & Landry 2019 MSB (10.15252/msb.20198871) | Human | heteromer status, divergence | partly | heteromeric paralogs buffer less | — |
| **Dennler & Ryan 2025 NAR Genom Bioinform (10.1093/nargab/lqaf051)** — closest prior work | Human: main label DepMap-inferred SL over all Ensembl 111 pairs (116/23,734 ≈0.5%); screen-based "SL-Lenient" set (403/5,316) secondary. Yeast set too | Foldseek on AFDB, ESM2-15B / ProtT5 embedding distances (ProtT5 generally beat ESM-2; pooling mattered more than distance metric), MMseqs2/Foldseek rank features ("paralogs", "species") | yes (but used AFDB structures and 15B PLM) | Identity "performed exceptionally well" on the **DepMap-inferred** human task; single-modality models only matched it; combined XGBoost beat it, also under family-split CV. Only the "paralogs" rank feature is nearly orthogonal to identity (ρ≈0.13–0.18); "species" correlates ρ≈0.58–0.68 | family split |
| Flister et al. 2025 Cell Rep (10.1016/j.celrep.2025.116512) | Human, all 36,648 Ensembl paralog pairs, 2 lines | screen + XGBoost (published version; preprint used weighted RF) | no | genome-wide SL prevalence <0.5%; no held-out AUPRC or identity baseline reported → not a comparator | — |
| Kebabci et al. 2026 bioRxiv (10.64898/2026.01.19.700065) | Human, cell-line-specific, trained on Ito | omics | no | identity alone barely above chance on the Klingbeil screen (survey docs disagree on exact numbers: AP≈0.06 vs 0.04, or 0.04–0.06 vs ≈0.02 — to be recomputed, not relied on) | external screen |
| Cusack et al. 2021 MBE (10.1093/molbev/msab111) | Arabidopsis genetic redundancy | expression, conservation, sequence, duplication | partly | reciprocal-best-match among top features | — |
| TripleThreat (ICLR 2026) [abstract only] | PLM representations | PLM similarity | yes | PLM similarity largely confounded with identity | — |

Yeast determinants of duplicate redundancy (Gu 2003; Dean 2008; DeLuna 2008; Guan 2007; Hakes 2007; VanderSluis 2010; Li 2010; Kuzmin 2020): asymmetric divergence (one copy faster) anti-correlates with redundancy; small families and closest pairs are enriched for SL. Two claims are weaker than often stated (round-1 prior-work review): (i) "WGD buffers more than SSD" — among non-ribosomal pairs Dean 2008 found 17/120 WGD vs 4/45 SSD SL, not significant; ribosomal-protein pairs dominate (30/36 redundant vs 39/165 others); (ii) "redundancy decays with divergence" — Li 2010 reports decay, Dean 2008 reports pervasive, persistent redundancy (~100 My) without a divergence trend. Heteromer evidence points in opposite directions: in yeast heteromeric paralogs are more often *dependent* (Diss 2017; Dandage & Landry 2019), in human physically interacting/same-complex paralogs are *more* often SL (De Kegel 2019/2021). Cross-species conservation of *specific* genetic interactions is modest (Tischler 2008; Dixon 2008) — and in any case using another species' labels as features is forbidden here.

**Gaps:** no labelled paralog-SL dataset for insects or crops; no study of yeast→human transfer of paralog-SL *predictors*; no sequence-only tool that calls paralogs de novo and is evaluated end-to-end.

## 2. Evaluation pitfalls that shape our protocol
1. **Pair-input leakage** (Park & Marcotte 2012; Bernett 2024): random pair splits share genes/families between train and test and inflate performance, sometimes to near-random once removed. → split by homology cluster built across all species jointly.
2. **Hub/degree bias**: a few genes appear in many positive pairs; a per-gene prior can look predictive. → report a gene-prior baseline and per-gene-capped variants.
3. **Label noise**: screens share ~31% of hits; same-line replicate agreement ≈ precision 0.40 / recall 0.32; scoring method changes hits (Ajmal 2025; SLKB). → robust consensus positives, tested negatives only, ambiguous excluded; sensitivity analyses with lenient labels.
4. **Ascertainment bias**: screened human pairs are 73% WGD (vs 19% of all), 89% closest pairs (vs 13%), median family size 3 (vs 11), median identity 54% (vs 32%). AUPRC on screened pairs ≠ AUPRC genome-wide. → report candidate-set composition; report expected precision at genome-wide prevalence.
5. **AUPRC baseline = prevalence**; prevalence differs by species and label tier. → always report prevalence and lift; compare methods on identical pairs.
6. **Fair identity baseline**: choose the best identity variant (directional min/max/mean, local/global) on development data only.

## 3. Label sources (decision-relevant)

**Human** (all primary tables publish non-hits, i.e. tested negatives): Dede 2020 (403 pairs × 3 lines, CC BY), Thompson 2021 (645 paralog pairs × 3, CC BY), Parrish 2021 (1,030 × 2, CC BY-NC-ND), Ito 2021 (4,858 × 11, no open license), CHyMErA 2020 (~672 × 2, no open license), in4mer/Esmaeili Anvar 2024 + Chou 2025 rescoring (4,540 pairs × 4, CC BY), Klingbeil 2024 (~2,500 × 22, CC BY-NC-ND), Harle 2025 (262 paralog × 27, CC BY), Chan 2025 (262 paralog × 10, CC BY), Flister 2025 (genome-wide, 2 lines; full non-hit table [U]). Ready-made cross-check: Dennler & Ryan 2025 Zenodo 14973633 (human 403 SL / 5,316, yeast set, CC BY 4.0; contains data-derived columns we must not use as features). Databases (SLKB GPL-3 with relative top-10% calls; SynLethDB mostly positives) are unsuitable as label sources. Expected union ≈8–10k screened pairs; robust tier ≈300–600 positives at 5–10% prevalence [estimate].

**Yeast**: Costanzo 2016 SGA raw pairwise file (TheCellMap, 521 MB; ε, p, single and double mutant fitness for every tested pair → tested negatives available; authors' thresholds: lenient p<0.05, intermediate |ε|>0.08, stringent ε<−0.12 / ε>0.16). Independent: Dean 2008 (197 pairs with fitness, 48 lethal; CC BY). Kuzmin 2020 (240 WGD pairs, viable double mutants; Dryad CC0, browser-only download). BioGRID low-throughput "Synthetic Lethality" (MIT) for positives only. Expected ≈0.8–1.8k labelled paralog pairs [estimate].

**Unseen species (held out):** Arabidopsis — Cusack 2021 Zenodo (CC BY): strict set 51 lethal / 110 non-lethal double mutants; broad 190/110; 552/553 genes map to UniProt. S. pombe — E-MAP genetic interaction data (Ryan 2012 and others; negatives need raw score matrix). C. elegans — Tischler 2006 combinatorial RNAi, 143 pairs, 16 positive. Drosophila — no systematic paralog screen; FlyBase lethal double-allele phenotypes without negatives (enrichment check at most). No insect or crop AUPRC test set exists — a stated limitation.

**Proteomes:** UniProt reference proteomes (one canonical sequence per gene; `.gene2acc` and `.idmapping` files for mapping). Reference panel: Quest-for-Orthologs 2026_02, 51 eukaryotes, 976k proteins, 309 MB (already downloaded, `data/raw/qfo/`).

## 4. Compute facts measured on this machine (Low Power Mode was ON — pessimistic)
- MMseqs2 human self-search: 4–10 min, ≤1.5 GB; human vs QfO panel: ~9 min, 4.2 GB. DIAMOND faster for self-search but 19 GB footprint on the panel.
- Embedding the full human proteome (MPS fp16): ESM-2 35M 66 min; 150M 104 min; 650M 5.4 h; ESM C 300M 2.6 h; ProtT5-XL 10.9 h. fp16 ≡ fp32 (cos ≥0.99999).
- ProstT5→3Di via Foldseek: ~19 h for human — only viable restricted to candidate-pair proteins or as an optional tier.
- ESMFold / AF2 for whole proteomes: infeasible here.

## 5. What this implies (carried into DESIGN.md)
- Beating identity is hard on the DepMap-inferred human task (Dennler) and identity is only weakly above chance on some screens (Kebabci); the most promising sequence-only signals beyond identity are family-context/rank features, conservation/age against an outgroup panel, divergence asymmetry, and (to be tested conditional on identity) PLM and 3Di similarity.
- Evaluation must be homology-cluster-grouped across species, with leave-one-species-out and genuinely unseen species (Arabidopsis first), with prevalence-aware reporting and leakage baselines.
