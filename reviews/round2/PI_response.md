# PI response to review round 2 → DESIGN v3

Five reviews: external PI, statistics, data, algorithms+comparative genomics, AI-agent failure modes. All must-items accepted; the reviewers' exact specifications are adopted as normative (referenced from DESIGN v3 header) rather than re-transcribed, to avoid copy errors.

## Accepted
- **Label universe** = tested pairs ∩ pinned Ensembl Compara paralogs, both species; Compara only defines the universe (external PI M1, data R2-4).
- **Yeast:** testability floor on expected DMF, negative margin, strain-suffix whitelist excluding `_S` suppressor strains, count-neutral aggregation, Dean/Kuzmin own rule, nested tiers, linked pairs (data R2-1/2/3/9/10/11). The v2 claim that DMF ≤ 0.5 raises positive counts is withdrawn (external PI, statistics); thresholds are not re-tuned.
- **Human:** uniform harmonised-LFC rule, per-study exceptions, borderline allowance floor(0.1n), Cellosaurus merging, Ito 90th-percentile cut dropped (data R2-5/6/7).
- **Splits:** labelled pairs contracted to supernodes; connected components on stricter edges by default, Leiden only as capped fallback with pre-registered parameters; leakage as dose–response by tertile (external PI M3, statistics 3–4).
- **Power:** pre-registered simulation inputs, no score on test pairs, test fraction chosen in [0.30, 0.40] for power ≥0.8 at T = 0.04 AP (statistics 1–2).
- **Nested-CV fallback dropped** (external PI M2). Statistics item 8 (how to run it) is therefore moot.
- Primary pair population and same-pair scoring of B1; bootstrap strata by species composition; reporting hierarchy (statistics 5–6, should-items).
- **Masking:** absolute d ceilings, adopted only if the thinning test shows a recoverable drop, with all statistical conditions (statistics 7, algorithms change 7, external PI should-item).
- **Ban on within-proteome-constant features** + species-identification control (external PI M4, algorithms 7e).
- Search mechanics: `-a` backtrace, per-species panel search, bitscore thresholds, pinned scoring settings (algorithms 1–3).
- Divergence calibration on a core gene set with Kimura distance, fixed bins, loss-aware duplication classification (algorithms 4–6).
- Panel gap-filling rule as tested by the reviewer, 80 species (algorithms 8); congeners of test species not added; no extra Saccharomyces.
- Self-exclusion text corrected — 0.98 also removes chimp and gorilla (algorithms 9).
- Tier C starts with one PLM; embedding PCs dropped; runtime budget part of adoption (external PI should-items).
- Auto-research before the single test access, with backlog and stop rule (external PI).
- Data hygiene: full Dennler files moved to `data/raw/_quarantine/dennler2025/` with manifest updated (done 2026-10-28); sealed directory covers splits and unseen labels; loader refuses sealed paths; AGENT_RULES extended (agent B1–B3).
- `final_test.py` preconditions (agent B4/G15).

## Modified
- Gates: G1–G6, G8–G12, G14, G15 adopted as specified; G7's git-ancestry ordering checks are kept as a documented manual check at milestone review (cost/benefit), the rest of G7 implemented; G13 (git-ancestry check of experiment registry) replaced by committing the registry row before each run.
- `_survey/` files: kept but registered in the manifest as survey copies (not used by code), rather than deleted, to preserve the audit trail of the literature agents.

## Not adopted
- None of the must-items rejected.

## Decision on round 3
A single consistency reviewer (fresh agent) checks DESIGN v3 against both rounds for contradictions or omissions, while M1 data work that does not depend on open points starts. Substantive redesign is not expected.
