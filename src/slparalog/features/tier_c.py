"""Tier C features (DESIGN v3 §4 item 8): pairwise similarity of mean-pooled PLM embeddings.

Input: an (n x d) embedding matrix covering the WHOLE input proteome (one row per protein; the
cached ESM-2 150M mean-pooled embeddings, config/plm.yaml) and the called pairs as row indices.
Only pairwise quantities are produced; no per-gene embedding coordinate or PC is ever a feature
(design review: such features let a model memorise which gene types were screened).

Definitions (registry descriptions in registry._TIER_C):
* ``tc_cos`` / ``tc_euc``: cosine and Euclidean distance of the raw mean embeddings (float64).
* ``tc_cos_c``: cosine after centring every embedding by the input-proteome mean embedding (the mean
  over all n rows, recomputed at inference on the input proteome). The Euclidean distance of centred
  embeddings equals ``tc_euc`` (translation invariance) and is therefore not a separate feature.
* Neighbourhood context, in centred-cosine space over all n proteome proteins. For gene g with
  partner p, ``above(g, p)`` = number of proteins x not in {g, p} with cos_c(g, x) > cos_c(g, p) + eps
  (exact, computed against the full proteome row by blocked float32 matrix products; eps absorbs
  float32 noise so identical sequences tie). Rank = above + 1 (ties share the best rank, as in Tier A),
  capped at ``k + 1``. ``tc_mnn`` = both ranks == 1 (mutual nearest neighbours; analogue of Tier A
  ``rbh``). ``tc_above_frac`` = above / (n - 2) (uncapped, proteome-size normalised).
  ``tc_gap`` = cos_c(g, p) - max over x not in {g, p} of cos_c(g, x) (analogue of Tier A other_gap).
  Pair-level min/max over the two genes, so every feature is symmetric in (a, b).

Ranks are computed exactly (full-row counts), so no top-k approximation is involved; the cap only
bounds the rank scale.

Caveats (independent review, reviews/results/tierC_review.md):
* Reproducibility at the tie boundary: the eps comparison is done on float32 cosines from BLAS
  matrix products, so a competitor lying within ~1e-7 of (cos(g, p) + eps) can fall on either side
  depending on precision / BLAS / block layout. ``tc_rank_*``, ``tc_mnn`` and ``tc_above_frac_*`` are
  therefore not bit-reproducible across precisions or BLAS builds in such cases (review: 1 of 300
  human pairs differed by one rank against a float64 reference; the competitor was 1.9e-8 from the
  threshold). ``tc_gap_*`` carries float32 noise (~1e-6).
* Saturation: the cap at k + 1 = 201 binds for most pairs of large families; in human 74.7% of
  ``tc_rank_max`` and 59% of ``tc_rank_min`` values equal 201 (arabidopsis/worm similar; see also
  ``frac_rank_max_capped`` in results/features/tierC_sanity.json), so the rank features carry little
  information there. ``tc_above_frac_*`` is the uncapped counterpart. Everything is numpy/numba on CPU: this module never imports torch, so it
is safe in a process that fits LightGBM (the embedding step runs in its own process, ``plm.py``).
"""

from __future__ import annotations

import hashlib
import json

import numpy as np
import pandas as pd
from numba import njit, prange

from .registry import BY_NAME, FEATURES, TIER_C_COLUMNS

DEFAULTS = {
    "rank_cap_k": 200,          # ranks are capped at k + 1
    "tie_eps": 1e-6,            # cosine differences <= eps count as ties (shared best rank)
    "neighbour_space": "proteome_mean_centred_cosine",
    "block_rows": 1024,         # query rows per float32 matrix product block
}


def settings_hash(settings: dict | None = None) -> str:
    s = {**DEFAULTS, **(settings or {})}
    s.pop("block_rows", None)  # does not change results (tested)
    return hashlib.sha256(json.dumps(s, sort_keys=True).encode()).hexdigest()


@njit(parallel=True, cache=True)
def _row_top2(S, top1_idx, top1_val, top2_val):
    """Per row: index/value of the largest entry and the second-largest value (self already -inf)."""
    for r in prange(S.shape[0]):
        b1 = -np.inf
        b2 = -np.inf
        i1 = -1
        for x in range(S.shape[1]):
            v = S[r, x]
            if v > b1:
                b2 = b1
                b1 = v
                i1 = x
            elif v > b2:
                b2 = v
        top1_idx[r] = i1
        top1_val[r] = b1
        top2_val[r] = b2


@njit(parallel=True, cache=True)
def _count_above(S, erow, edst, eps, above, val):
    """For each directed edge (row r of S, partner column p): number of entries > S[r, p] + eps."""
    for e in prange(erow.shape[0]):
        r = erow[e]
        v = S[r, edst[e]]
        thr = v + eps
        c = 0
        for x in range(S.shape[1]):
            if S[r, x] > thr:
                c += 1
        above[e] = c
        val[e] = v


def _unit_rows(X: np.ndarray) -> np.ndarray:
    nrm = np.linalg.norm(X, axis=1, keepdims=True)
    return X / np.maximum(nrm, 1e-12)


def _pair_rowdot(A: np.ndarray, lo: np.ndarray, hi: np.ndarray, chunk: int = 65536) -> np.ndarray:
    out = np.empty(len(lo), np.float64)
    for s in range(0, len(lo), chunk):
        a, b = A[lo[s:s + chunk]], A[hi[s:s + chunk]]
        out[s:s + chunk] = np.einsum("ij,ij->i", a, b)
    return out


def _pair_dist(X: np.ndarray, lo: np.ndarray, hi: np.ndarray, chunk: int = 65536) -> np.ndarray:
    out = np.empty(len(lo), np.float64)
    for s in range(0, len(lo), chunk):
        d = X[lo[s:s + chunk]] - X[hi[s:s + chunk]]
        out[s:s + chunk] = np.sqrt(np.einsum("ij,ij->i", d, d))
    return out


def neighbourhood(Xn32: np.ndarray, src: np.ndarray, dst: np.ndarray, eps: float,
                  block_rows: int = 1024) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Directed-edge neighbourhood statistics in the space of unit rows ``Xn32`` (float32).

    Returns (above, cos, best_other) per edge g=src -> p=dst: count of proteins x not in {g, p}
    with cos(g, x) > cos(g, p) + eps, the float32 cos(g, p), and max cos(g, x) over x not in {g, p}
    (-inf if the proteome has only the two proteins).
    """
    n = Xn32.shape[0]
    m = len(src)
    above = np.zeros(m, np.int64)
    val = np.zeros(m, np.float32)
    best_other = np.full(m, -np.inf, np.float64)
    if m == 0:
        return above, val.astype(np.float64), best_other
    order = np.argsort(src, kind="mergesort")
    s_src, s_dst = src[order], dst[order]
    genes = np.unique(s_src)
    XT = np.ascontiguousarray(Xn32.T)
    for b0 in range(0, len(genes), block_rows):
        blk = genes[b0:b0 + block_rows]
        S = Xn32[blk] @ XT                      # (B x n) float32 cosines
        S[np.arange(len(blk)), blk] = -np.inf   # exclude self
        e0 = np.searchsorted(s_src, blk[0], side="left")
        e1 = np.searchsorted(s_src, blk[-1], side="right")
        erow = np.searchsorted(blk, s_src[e0:e1]).astype(np.int64)
        edst = s_dst[e0:e1].astype(np.int64)
        a = np.empty(e1 - e0, np.int64)
        v = np.empty(e1 - e0, np.float32)
        _count_above(S, erow, edst, np.float32(eps), a, v)
        t1i = np.empty(len(blk), np.int64)
        t1v = np.empty(len(blk), np.float32)
        t2v = np.empty(len(blk), np.float32)
        _row_top2(S, t1i, t1v, t2v)
        bo = np.where(t1i[erow] == edst, t2v[erow], t1v[erow]).astype(np.float64)
        idx = order[e0:e1]
        above[idx], val[idx], best_other[idx] = a, v, bo
    return above, val.astype(np.float64), best_other


def compute(pairs: pd.DataFrame, emb: np.ndarray, settings: dict | None = None,
            threads: int | None = None) -> tuple[pd.DataFrame, dict]:
    """Return (feature frame aligned with ``pairs`` rows, summary).

    ``pairs`` has integer columns i, j (row indices into ``emb``, any orientation, no self pairs).
    ``emb`` is the (n x d) embedding matrix of the whole input proteome (centring and neighbourhoods
    use all n rows, not only the genes in pairs).
    """
    st = {**DEFAULTS, **(settings or {})}
    if threads is not None:
        import numba
        numba.set_num_threads(max(1, min(int(threads), numba.config.NUMBA_NUM_THREADS)))
    X = np.asarray(emb, dtype=np.float64)
    n = X.shape[0]
    if not np.all(np.isfinite(X)):
        raise ValueError("embedding matrix has non-finite values")
    pi = pairs["i"].to_numpy(np.int64)
    pj = pairs["j"].to_numpy(np.int64)
    if np.any(pi == pj):
        raise ValueError("self pairs are not allowed")
    if len(pi) and (min(pi.min(), pj.min()) < 0 or max(pi.max(), pj.max()) >= n):
        raise ValueError("pair index outside the embedding matrix")
    if n < 3:
        raise ValueError("need at least 3 proteins for neighbourhood features")
    lo, hi = np.minimum(pi, pj), np.maximum(pi, pj)
    npairs = len(lo)

    F: dict[str, np.ndarray] = {}
    Xu = _unit_rows(X)
    F["tc_cos"] = _pair_rowdot(Xu, lo, hi)
    F["tc_euc"] = _pair_dist(X, lo, hi)
    mu = X.mean(axis=0)
    Xc = _unit_rows(X - mu)
    F["tc_cos_c"] = _pair_rowdot(Xc, lo, hi)
    del Xu

    src = np.concatenate([lo, hi])
    dst = np.concatenate([hi, lo])
    above, _, best = neighbourhood(Xc.astype(np.float32), src, dst, float(st["tie_eps"]),
                                   int(st["block_rows"]))
    k = int(st["rank_cap_k"])
    rank = np.minimum(above + 1, k + 1)
    frac = above / float(n - 2)
    cosc = np.concatenate([F["tc_cos_c"], F["tc_cos_c"]])
    gap = cosc - best

    def mm(v):
        a, b = v[:npairs], v[npairs:]
        return np.minimum(a, b), np.maximum(a, b)

    F["tc_rank_min"], F["tc_rank_max"] = mm(rank)
    F["tc_mnn"] = (F["tc_rank_max"] == 1).astype(np.int8)
    F["tc_above_frac_min"], F["tc_above_frac_max"] = mm(frac)
    F["tc_gap_min"], F["tc_gap_max"] = mm(gap)

    out = pd.DataFrame(F, index=pairs.index)
    for f in FEATURES:
        if f.tier == "C" and f.name.endswith("_pct") and f.name[:-4] in BY_NAME:
            out[f.name] = out[f.name[:-4]].rank(method="average", pct=True)
    missing = [c for c in TIER_C_COLUMNS if c not in out.columns]
    extra = [c for c in out.columns if c not in BY_NAME]
    if missing or extra:
        raise RuntimeError(f"registry mismatch: missing={missing} extra={extra}")
    out = out[TIER_C_COLUMNS].astype({c: BY_NAME[c].dtype for c in TIER_C_COLUMNS})
    involved = np.unique(src)
    summary = {
        "n_proteome": int(n),
        "embedding_dim": int(X.shape[1]),
        "n_pairs": int(npairs),
        "n_genes_in_pairs": int(len(involved)),
        "n_rank_capped_edges": int((above + 1 > k + 1).sum()),
        "frac_mnn": float(F["tc_mnn"].mean()) if npairs else None,
        "settings": {k_: v for k_, v in st.items()},
        "settings_hash": settings_hash(st),
    }
    return out, summary


# --------------------------------------------------------------------------- embedding step
def embed_subprocess(seqs: list[str], workdir, cache_dir, model_key: str | None = None) -> np.ndarray:
    """Embed ``seqs`` (row order kept) in a separate process (``python -m slparalog.plm``), so that
    torch is never imported here. Returns the (n x d) float32 matrix."""
    from pathlib import Path

    from .. import fasta, plm  # plm's module level is torch-free

    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    fa = wd / "plm_input.fasta"
    fasta.write_fasta(fa, ((str(k), s) for k, s in enumerate(seqs)))
    out = wd / "plm_input.npy"
    mat = plm.embed_fasta_subprocess(fa, cache_dir, out, model_key)
    if mat.shape[0] != len(seqs):
        raise RuntimeError(f"embedding rows {mat.shape[0]} != sequences {len(seqs)}")
    return mat
