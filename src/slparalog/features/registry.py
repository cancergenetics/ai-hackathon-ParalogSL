"""Feature registry (DESIGN v3 §4; gate G8 in reviews/round2/agent_failure_modes.md).

Every column a model may use is listed here with its name, source and description. Models accept
only registered columns. ``source`` is one of ALLOWED_SOURCES; no registered name may match
guards.FORBIDDEN_NAME_RE (checked by tests/test_feature_registry.py). ``kind`` is "feature" (model input)
or "qc" (output flag; also sequence-derived, may be used as a feature after ablation).
``pct`` marks features that also get a proteome-percentile version ``<name>_pct`` (rank percentile
over the full called-pair set of the run, identical at training and inference).
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass

from ..guards import FORBIDDEN_NAME_RE

ALLOWED_SOURCES = ("self_search", "panel", "plm", "structure", "sequence")


@dataclass(frozen=True)
class Feature:
    name: str
    source: str
    description: str
    dtype: str = "float64"
    kind: str = "feature"
    tier: str = "A"
    pct: bool = False


_TIER_A: list[Feature] = [
    # --- 4.1 similarity (self-search realignment, both orientations, backtraced) ---
    Feature("ident_local_min", "self_search", "local identity nident/alnlen, min over the two orientations", pct=True),
    Feature("ident_local_max", "self_search", "local identity nident/alnlen, max over the two orientations", pct=True),
    Feature("ident_ens_min", "self_search", "Ensembl-style identity: identical positions / length of each protein, min (= / longer); identical positions = max nident over orientations", pct=True),
    Feature("ident_ens_max", "self_search", "Ensembl-style identity, max (= identical positions / shorter protein)", pct=True),
    Feature("ident_global", "self_search", "global Needleman-Wunsch identity (BLOSUM62, 11/1, end gaps penalised): identical / alignment columns", pct=True),
    Feature("cov_min", "self_search", "aligned-span coverage of each protein (max over orientations), min over the pair"),
    Feature("cov_max", "self_search", "aligned-span coverage of each protein (max over orientations), max over the pair"),
    Feature("bits_pair", "self_search", "pair bitscore = max over orientations (not an E-value)"),
    Feature("bsr_min", "self_search", "bits_pair / self-bitscore of each protein, min", pct=True),
    Feature("bsr_max", "self_search", "bits_pair / self-bitscore of each protein, max", pct=True),
    Feature("len_ratio", "sequence", "shorter / longer protein length"),
    # --- 4.2 family context (called-pair graph; per-gene lists use bits(g->p) with g as query) ---
    Feature("hits_all_min", "self_search", "number of called partners per gene, min over the pair", dtype="int32"),
    Feature("hits_all_max", "self_search", "number of called partners per gene, max over the pair", dtype="int32", pct=True),
    Feature("hits_r10_min", "self_search", "partners with bits(g->p) >= 0.10 x selfbits(g), min", dtype="int32"),
    Feature("hits_r10_max", "self_search", "partners with bits(g->p) >= 0.10 x selfbits(g), max", dtype="int32"),
    Feature("hits_r25_min", "self_search", "partners with bits(g->p) >= 0.25 x selfbits(g), min", dtype="int32"),
    Feature("hits_r25_max", "self_search", "partners with bits(g->p) >= 0.25 x selfbits(g), max", dtype="int32"),
    Feature("hits_r50_min", "self_search", "partners with bits(g->p) >= 0.50 x selfbits(g), min", dtype="int32"),
    Feature("hits_r50_max", "self_search", "partners with bits(g->p) >= 0.50 x selfbits(g), max", dtype="int32", pct=True),
    Feature("hits_r75_min", "self_search", "partners with bits(g->p) >= 0.75 x selfbits(g), min", dtype="int32"),
    Feature("hits_r75_max", "self_search", "partners with bits(g->p) >= 0.75 x selfbits(g), max", dtype="int32"),
    Feature("rank_min", "self_search", "rank of the partner in each gene's called list by bits(g->p) (1 = best; ties share the best rank), min", dtype="int32"),
    Feature("rank_max", "self_search", "partner rank, max over the pair (proteome 'paralogs rank' = rank - 1)", dtype="int32", pct=True),
    Feature("rbh", "self_search", "reciprocal best hit: both partner ranks == 1", dtype="int8"),
    Feature("shared_third", "self_search", "number of genes called with both a and b", dtype="int32", pct=True),
    Feature("shared_third_jacc", "self_search", "shared_third / |N(a) u N(b) minus {a,b}| (0 if empty)"),
    Feature("other_ident_min", "self_search", "local identity of each gene to its best other partner (by bits; 0 if none), min"),
    Feature("other_ident_max", "self_search", "local identity of each gene to its best other partner, max"),
    Feature("other_gap_min", "self_search", "ident(g->partner) - ident(g->best other partner), min over the pair"),
    Feature("other_gap_max", "self_search", "ident(g->partner) - ident(g->best other partner), max over the pair"),
    Feature("mc_deg_min", "self_search", "mutual-coverage degree (links with both coverages >= 0.5), min", dtype="int32"),
    Feature("mc_deg_max", "self_search", "mutual-coverage degree, max", dtype="int32"),
    Feature("fam_mc_size", "self_search", "local family size |{a,b} u Nmc(a) u Nmc(b)| on mutual-coverage links", dtype="int32", pct=True),
    # --- 4.3 per-gene sequence properties ---
    Feature("len_min", "sequence", "protein length, min over the pair", dtype="int32", pct=True),
    Feature("len_max", "sequence", "protein length, max over the pair", dtype="int32"),
    Feature("lowcomp_frac_min", "sequence", "fraction of residues masked by MMseqs2 tantan (compositional bias), min"),
    Feature("lowcomp_frac_max", "sequence", "fraction of residues masked by MMseqs2 tantan, max"),
    # --- QC flags ---
    Feature("fragment_suspect", "self_search", "1 if either protein is shorter than 0.5 x the median length of its called partners (provisional self-search definition; panel-based definition per algorithms_comparative change 11 pending Tier B)", dtype="int8", kind="qc"),
    Feature("near_identical", "self_search", "1 if ident_local_max >= 0.99 (possible allelic/isoform copy)", dtype="int8", kind="qc"),
    Feature("n_dir_called", "self_search", "number of search directions with E' <= threshold (1 or 2)", dtype="int8", kind="qc"),
]


# Tier B (DESIGN v3 §4 items 4-7; algorithms_comparative changes 5, 6, 10, 11). Computed on the retained
# panel species after self-exclusion and the optional --panel-dmin mask (src/slparalog/features/tier_b.py).
# Only pair-varying summaries are registered: no per-bin feature and no bin-missingness pattern
# (DESIGN §4: no feature may be constant within a proteome).
_TIER_B: list[Feature] = [
    # --- duplication classification / depth on the d axis (change 6) ---
    Feature("tb_depth", "panel", "duplication depth on the d (Kimura) axis: (L+U)/2; U/2 if no pre species; L if no post species; bin-weighted discordance minimiser if L >= U; NaN if no pre and no post", tier="B", pct=True),
    Feature("tb_depth_range", "panel", "U - L (U = min d over post species, L = max d over pre; L := 0 if no pre; NaN if no post); negative under conflict", tier="B"),
    Feature("tb_discordant_frac", "panel", "bin-weighted (1/n_bin) fraction of pre+post species discordant with the depth (0 unless L >= U; NaN if no pre and no post)", tier="B"),
    Feature("tb_younger_than_nearest", "panel", "censoring flag: no pre species (duplication younger than the nearest classified species)", dtype="int8", tier="B"),
    Feature("tb_older_than_panel", "panel", "censoring flag: no post species (duplication older than the panel)", dtype="int8", tier="B"),
    Feature("tb_no_evidence", "panel", "indicator: no pre and no post species (depth undefined)", dtype="int8", tier="B"),
    Feature("tb_rel_rate", "panel", "relative rate d(A,B) / (2 depth); d(A,B) = Kimura distance of the mean backtraced pair identity over both orientations", tier="B"),
    Feature("tb_frac_pre", "panel", "fraction of retained species classified pre", tier="B"),
    Feature("tb_frac_post", "panel", "fraction of retained species classified post", tier="B"),
    Feature("tb_frac_undet", "panel", "fraction of retained species classified undetermined (loss, rate asymmetry)", tier="B"),
    Feature("tb_frac_absent", "panel", "fraction of retained species with no hit >= 50 bits for either member", tier="B"),
    Feature("tb_n_post_single", "panel", "number of post species with a = b whose higher-scoring member has no second hit in s within 0.9 x its best bits (1:2 orthology pattern)", dtype="int32", tier="B"),
    # --- conservation (change 6 derived) ---
    Feature("tb_cons_min", "panel", "per-gene fraction of retained species with a hit >= 50 bits, min over the pair", tier="B", pct=True),
    Feature("tb_cons_max", "panel", "per-gene conservation fraction, max over the pair", tier="B"),
    Feature("tb_cons_both", "panel", "fraction of retained species where both members have a hit >= 50 bits", tier="B"),
    # --- species rank (item 7) ---
    Feature("tb_species_rank_min", "panel", "per gene g with partner p: fraction of retained species whose best hit scores higher to g than p does (query-side bits), min", tier="B"),
    Feature("tb_species_rank_max", "panel", "species rank, max over the pair", tier="B", pct=True),
    # --- family duplicability (change 11) ---
    Feature("tb_dupability_min", "panel", "per gene: median over retained species with a hit of the number of hits >= 0.5 x self-bits, min", tier="B"),
    Feature("tb_dupability_max", "panel", "family duplicability, max over the pair", tier="B"),
    # --- divergence asymmetry (change 10) ---
    Feature("tb_asym_norm", "panel", "median over the <=3 nearest post species of |d_AO - d_BO| / (d_AO + d_BO), O = RBH partner of A or B, identities over O columns covered by both alignments", tier="B"),
    Feature("tb_asym_abs", "panel", "median |d_AO - d_BO| over the same outgroups", tier="B"),
    Feature("tb_outgroup_ident_max", "panel", "median identity of the more conserved copy to the outgroup (shared O columns)", tier="B"),
    # --- QC ---
    Feature("tb_n_outgroups", "panel", "number of outgroup triplets used for the asymmetry features (0-3)", dtype="int8", kind="qc", tier="B"),
    Feature("tb_fragment_suspect", "panel", "1 if either protein is shorter than 0.5 x the median length of its best hits in retained panel species (change 11 definition)", dtype="int8", kind="qc", tier="B"),
]


# Tier C (DESIGN v3 §4 item 8): pairwise similarity of mean-pooled PLM embeddings (ESM-2 150M,
# config/plm.yaml) and embedding-neighbourhood context over the whole input proteome
# (src/slparalog/features/tier_c.py). Pairwise only: per-gene embedding coordinates/PCs are excluded
# (design review: memorisation of screened gene types). Neighbourhoods use proteome-mean-centred cosine.
_TIER_C: list[Feature] = [
    Feature("tc_cos", "plm", "cosine similarity of the two mean-pooled embeddings", tier="C", pct=True),
    Feature("tc_euc", "plm", "Euclidean distance between the two mean-pooled embeddings (identical for proteome-mean-centred embeddings)", tier="C", pct=True),
    Feature("tc_cos_c", "plm", "cosine after centring each embedding by the input-proteome mean embedding (mean over all proteins)", tier="C", pct=True),
    Feature("tc_rank_min", "plm", "rank of the partner among all other proteome proteins by centred cosine to g (1 = nearest; ties within 1e-6, compared in float32, share the best rank; not bit-reproducible for competitors within ~1e-7 of that boundary), capped at k+1 = 201, min over the pair (saturates at 201 for 59% of human pairs, reviews/results/tierC_review.md)", dtype="int32", tier="C"),
    Feature("tc_rank_max", "plm", "embedding partner rank (capped at 201; same float32 tie caveat), max over the pair (saturates at 201 for 74.7% of human pairs; tc_above_frac_* carries the uncapped information)", dtype="int32", tier="C", pct=True),
    Feature("tc_mnn", "plm", "mutual nearest neighbours in centred-embedding space (both partner ranks == 1; same float32 tie caveat)", dtype="int8", tier="C"),
    Feature("tc_above_frac_min", "plm", "fraction of proteome proteins other than the pair whose centred cosine to g exceeds the pair's by > 1e-6 (count / (n_proteome - 2); exact, uncapped; float32 tie caveat as tc_rank_min), min over the pair", tier="C"),
    Feature("tc_above_frac_max", "plm", "fraction of proteome proteins closer to g than its partner (centred cosine), max over the pair", tier="C"),
    Feature("tc_gap_min", "plm", "centred cos(g, partner) - max centred cos(g, x) over all other proteome proteins x, min over the pair", tier="C"),
    Feature("tc_gap_max", "plm", "centred-cosine gap to the best other proteome protein, max over the pair", tier="C"),
]


def _with_pct(feats: list[Feature]) -> list[Feature]:
    out: list[Feature] = []
    for f in feats:
        out.append(f)
    for f in feats:
        if f.pct:
            out.append(Feature(f"{f.name}_pct", f.source,
                               f"proteome percentile (average rank / n over called pairs) of {f.name}",
                               dtype="float64", kind=f.kind, tier=f.tier, pct=False))
    return out


FEATURES: list[Feature] = _with_pct(_TIER_A) + _with_pct(_TIER_B) + _with_pct(_TIER_C)
BY_NAME: dict[str, Feature] = {f.name: f for f in FEATURES}
TIER_A_COLUMNS: list[str] = [f.name for f in FEATURES if f.tier == "A"]
TIER_B_COLUMNS: list[str] = [f.name for f in FEATURES if f.tier == "B"]
TIER_C_COLUMNS: list[str] = [f.name for f in FEATURES if f.tier == "C"]
MODEL_FEATURES: list[str] = [f.name for f in FEATURES if f.kind == "feature"]


def registry_hash(features: list[Feature] | None = None) -> str:
    feats = FEATURES if features is None else features
    payload = json.dumps([asdict(f) for f in feats], sort_keys=True).encode()
    return hashlib.sha256(payload).hexdigest()


COLUMN_DEFINITION_FIELDS = ("name", "source", "dtype", "kind", "tier", "pct")


def column_definition(f: Feature) -> list:
    """The part of a feature's definition that determines the computed column (descriptions excluded)."""
    return [f.name, f.source, f.dtype, f.kind, f.tier, bool(f.pct)]


def columns_hash(columns, features: list[Feature] | None = None) -> str:
    """Hash of the column definitions (name, source, dtype, kind, tier, pct flag) of ``columns``.

    Order-independent (sorted by name) and independent of descriptions and of every other registered column, so
    registering a new tier or rewording a description leaves it unchanged. Unregistered names raise ValueError.
    ``features`` defaults to the registry; an alternative list (e.g. an older registry version) may be passed.
    """
    by_name = BY_NAME if features is None else {f.name: f for f in features}
    names = sorted(set(columns))
    missing = [c for c in names if c not in by_name]
    if missing:
        raise ValueError(f"columns not in the registry: {missing}")
    payload = json.dumps([column_definition(by_name[c]) for c in names]).encode()
    return hashlib.sha256(payload).hexdigest()


def table_columns_hash(columns, features: list[Feature] | None = None) -> str:
    """``columns_hash`` of the registered columns among ``columns`` (a table's gene_a/gene_b/etc. are skipped)."""
    by_name = BY_NAME if features is None else {f.name: f for f in features}
    return columns_hash([c for c in columns if c in by_name], features)


def validate(features: list[Feature] | None = None) -> list[str]:
    """Return a list of registry violations (empty = OK)."""
    feats = FEATURES if features is None else features
    errs = []
    names = [f.name for f in feats]
    if len(set(names)) != len(names):
        errs.append("duplicate feature names")
    for f in feats:
        if f.source not in ALLOWED_SOURCES:
            errs.append(f"{f.name}: source {f.source!r} not allowed")
        if FORBIDDEN_NAME_RE.search(f.name):
            errs.append(f"{f.name}: matches forbidden-name regex")
        if f.kind not in ("feature", "qc"):
            errs.append(f"{f.name}: bad kind {f.kind!r}")
    return errs
