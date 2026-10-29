# PI response to review round 1 → DESIGN v2

Seven reviews: external PI, prior work, statistics, data, algorithms/compute, AI-agent failure modes, comparative genomics. All "must" items were accepted (in some cases modified as noted); "should/could" items were accepted where cheap or decision-relevant. Rejections and modifications are explained.

## Labels
- Yeast and human SL now share one definition (significant negative GI + absolute double-mutant defect, both singles testable) — data D2, external PI 1–2, prior work M1–M2. **Modified:** yeast primary = ε<−0.08, p<0.05 and DMF ≤ 0.5 (pre-registered; not calibrated on Dean W=0 pairs, to avoid tuning on any outcome); ε<−0.35 / DMF ≤ 0.2 is the strict tier, ε<−0.12 the lenient tier. This also raises yeast positive counts (statistics item 2).
- Essential-member (TS/DAmP-only) yeast pairs untestable in the primary tier — accepted.
- Kuzmin/Dean viable pairs labelled by the same ε/DMF rule, not as blanket negatives — accepted.
- Human robust tier: ≥2 lines AND (≥2 studies OR penetrance ≥0.25), negatives with borderline margin; ≥2-studies tier added; Flister out of primary — accepted.
- Gene, not accession, is the unit; per-gene proteome expansion — accepted (data D5).
- Label universe independent of the caller — accepted (data D6).
- Per-source rules pre-registered in `config/labels.yaml` before label code runs; known-answer count tests; no tuning on the 13-pair list or Dennler agreement — accepted (data D9, agent 7–8).
- BioGRID-LT only in a separate tier — accepted (data D8). Arabidopsis negatives re-checked for lethal singles — accepted (data D12).

## Splits and statistics
- Single-linkage clusters replaced: homology graph over labelled genes, large components partitioned by bitscore-weighted Leiden; leakage quantified per test pair and a strict no-hit subset reported — accepted (statistics 5–6, comparative M7, algorithms 4, external PI 6).
- Test fraction raised to ~30% with species-stratified positives; MDE simulated at M1 from the real composition, with a pre-registered switch to co-primary nested grouped CV if MDE > 0.10 AP — accepted (statistics 1–2, external PI 5).
- Primary criterion: species-mean ΔAP vs B1 with stratified cluster bootstrap (10k) CI excluding 0 and ΔAP>0 in each species; "clearly beats" = mean relative gain ≥ 0.20 — accepted (statistics 1, 3, 9). Arabidopsis pass/fail clause replaced by reported ΔAP + CI — accepted (statistics 4). LOSO and B2 as pre-registered secondary comparisons — accepted (external PI 4, prior work S4).
- Identity baseline: pre-registered variant set incl. Ensembl-style and global identity — accepted (prior work M5, statistics 8, algorithms B1).
- Experiment registry, one-SE rule, single final model — accepted (statistics 10, agent 20).
- Source-library prediction control — accepted (statistics 7).
- WGD vs SSD and ribosomal-pair stratified reporting, annotations used for evaluation only — accepted (prior work M3, comparative M6).

## Features, panel, compute
- Panel gap filling by a written rule (Saccharomycetaceae, eudicots, monocots, non-dipteran insects, nematodes) — accepted (comparative M4, external PI 8, prior work M4).
- Panel-density independence: divergence calibration D_s per species (median RBH identity), divergence-binned features, training-time masking of the closest bins, pre-registered panel-thinning test — accepted, and this also answers comparative M2 (query placement without a species tree) and agent 16 (training/inference asymmetry from primates in the panel).
- Self-exclusion rule changed to median RBH identity ≥0.98 (removes self/strain only); close relatives handled by masking augmentation rather than exclusion — **modified** from algorithms 6 (they proposed a divergence rule; we use masking so the model sees both settings).
- Panel search with per-species hit guarantees — accepted (comparative M1, algorithms 5). Self-search `--max-seqs ≥1000` with re-search on saturation; both orientations realigned; explicit self-scores; no raw E-values; symmetric features; percentile features over the full called-pair set — accepted (algorithms 1–3, 9).
- Tier D: no TM-score/LDDT from 3Di-only databases; ProstT5 encoder+CNN — accepted (algorithms 7). Tier D gated behind Tier C signal — accepted (external PI).
- torch/LightGBM process isolation — accepted (algorithms 10). Pinning of tools, models, environment — accepted (algorithms 11, agent 28).
- Robustness test on degraded proteomes and `near_identical`/`fragment_suspect` flags — accepted (comparative M5, external PI 7).
- ProtT5 > ESM-2 and pooling matters (Dennler) — Tier C updated to include ProtT5 conditional on small-model signal (prior work S3); full ProtT5 costs ≈11 h on human, so it is only run if cheaper PLMs show identity-conditional signal.

## Process (AI-agent failure modes)
- Sealed test directory + path-guard test, tag-before-test, access log, re-run rule, refit-never-scored — accepted (agent 1–4).
- Protected file hashes, feature registry, gene-name scanner, quarantine of files with forbidden columns (done 2026-10-27: `data/raw/_quarantine/`), report-number checker, review records with clean-checkout re-runs, `docs/AGENT_RULES.md` preamble — accepted (agent 6, 9, 12–13, 23–25).
- Anonymised pair IDs for feature developers (agent 14) — **partially accepted:** features are computed proteome-wide by code that never sees labels or gene names; feature cards (hypothesis + adoption rule) are written before results are seen. Anonymising IDs adds friction without protection beyond that.
- DATA.md backfill by script; benchmark outputs committed (summary copied to `results/bench/`) — accepted (agent 27, 29).

## Literature corrections
- LITERATURE.md fixed for Dennler's main label type and rank-feature correlations, contested "WGD > SSD" and "redundancy decays" claims, heteromer sign difference, Kebabci chance level, Flister model type (prior work S3).
