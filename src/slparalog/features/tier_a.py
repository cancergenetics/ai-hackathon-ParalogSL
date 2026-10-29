"""Tier A features (DESIGN v3 §4.1-4.3): self-search similarity, family context, per-gene properties.

Inputs are indexed by the pipeline's internal protein index (canonical, sequence-sorted order).
Every pair feature is a symmetric transform (min/max/union/intersection) of per-orientation or
per-gene quantities, so f(a, b) == f(b, a) by construction; ``compute`` accepts pairs in either
orientation (tests/test_tier_a.py swaps them).

Conventions (see registry descriptions):
* Directional statistics come from the realignment of both orientations (``search.realign``).
  If one orientation is missing from the realignment, the other one is used for both (counted in
  the run summary as ``n_pairs_missing_orientation``).
* Coverage of a protein = aligned span / length, max over the two orientations.
* Per-gene hit lists use bits(g->p) with g as the query. Rank ties share the best rank.
* "Best other partner" = the partner with the highest bits(g->p), ties broken by higher identity;
  if g has no other partner, other identity = 0 and the gap = ident(g->p).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from numba import njit, prange

from .. import nw
from .registry import BY_NAME, FEATURES, TIER_A_COLUMNS


@njit(parallel=True, cache=True)
def _intersect_counts(indptr, indices, pi, pj, out):
    for k in prange(pi.shape[0]):
        a0, a1 = indptr[pi[k]], indptr[pi[k] + 1]
        b0, b1 = indptr[pj[k]], indptr[pj[k] + 1]
        c = 0
        while a0 < a1 and b0 < b1:
            x = indices[a0]
            y = indices[b0]
            if x == y:
                c += 1
                a0 += 1
                b0 += 1
            elif x < y:
                a0 += 1
            else:
                b0 += 1
        out[k] = c


def _csr(n: int, i: np.ndarray, j: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Undirected graph (edges i-j) as CSR with sorted neighbour lists."""
    src = np.concatenate([i, j]).astype(np.int64)
    dst = np.concatenate([j, i]).astype(np.int64)
    order = np.lexsort((dst, src))
    src, dst = src[order], dst[order]
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(np.bincount(src, minlength=n), out=indptr[1:])
    return indptr, dst


def _intersections(n, ei, ej, pi, pj):
    indptr, indices = _csr(n, ei, ej)
    out = np.zeros(len(pi), dtype=np.int64)
    if len(pi):
        _intersect_counts(indptr, indices, np.asarray(pi, np.int64), np.asarray(pj, np.int64), out)
    deg = np.diff(indptr)
    return out, deg


def _orient(aln: pd.DataFrame, q: np.ndarray, t: np.ndarray) -> pd.DataFrame:
    key = pd.DataFrame({"query": q, "target": t})
    cols = ["query", "target", "bits", "nident", "alnlen", "qstart", "qend", "tstart", "tend"]
    return key.merge(aln[cols], on=["query", "target"], how="left", sort=False)


def compute(pairs: pd.DataFrame, aln: pd.DataFrame, seqs: list[str], lowcomp: np.ndarray,
            cfg: dict, threads: int = 1) -> tuple[pd.DataFrame, dict]:
    """Return (feature frame aligned with ``pairs`` rows, summary). ``pairs`` has columns i, j
    (internal indices, any orientation, no duplicates as unordered pairs) and n_dir_called."""
    fcfg = cfg["features"]
    n = len(seqs)
    lens = np.array([len(s) for s in seqs], dtype=np.int64)
    pi = pairs["i"].to_numpy(np.int64)
    pj = pairs["j"].to_numpy(np.int64)
    if np.any(pi == pj):
        raise ValueError("self pairs are not allowed")
    lo, hi = np.minimum(pi, pj), np.maximum(pi, pj)  # canonical orientation for all computations
    npairs = len(pi)

    # --- self scores ---
    selfrows = aln[aln["query"] == aln["target"]]
    selfbits = np.full(n, np.nan)
    selfbits[selfrows["query"].to_numpy()] = selfrows["bits"].to_numpy(np.float64)
    involved = np.unique(np.concatenate([lo, hi]))
    n_missing_self = int(np.isnan(selfbits[involved]).sum())

    # --- both orientations ---
    ab = _orient(aln, lo, hi)   # lo as query
    ba = _orient(aln, hi, lo)   # hi as query
    miss_ab = ab["bits"].isna().to_numpy()
    miss_ba = ba["bits"].isna().to_numpy()
    n_missing_orient = int((miss_ab ^ miss_ba).sum())
    n_missing_both = int((miss_ab & miss_ba).sum())

    def arr(df, c):
        return df[c].to_numpy(np.float64)

    bits_ab, bits_ba = arr(ab, "bits"), arr(ba, "bits")
    id_ab = arr(ab, "nident") / arr(ab, "alnlen")
    id_ba = arr(ba, "nident") / arr(ba, "alnlen")
    # coverage of lo / hi in each orientation (1-based inclusive coordinates)
    cov_lo_ab = (arr(ab, "qend") - arr(ab, "qstart") + 1) / lens[lo]
    cov_hi_ab = (arr(ab, "tend") - arr(ab, "tstart") + 1) / lens[hi]
    cov_hi_ba = (arr(ba, "qend") - arr(ba, "qstart") + 1) / lens[hi]
    cov_lo_ba = (arr(ba, "tend") - arr(ba, "tstart") + 1) / lens[lo]
    # fill a missing orientation from the other one
    bits_ab = np.where(miss_ab, bits_ba, bits_ab)
    bits_ba = np.where(miss_ba, bits_ab, bits_ba)
    id_ab = np.where(miss_ab, id_ba, id_ab)
    id_ba = np.where(miss_ba, id_ab, id_ba)
    nid_ab, nid_ba = arr(ab, "nident"), arr(ba, "nident")

    with np.errstate(invalid="ignore", divide="ignore"):
        cov_lo = np.fmax(cov_lo_ab, cov_lo_ba)
        cov_hi = np.fmax(cov_hi_ab, cov_hi_ba)
        nid = np.fmax(nid_ab, nid_ba)
        bits_pair = np.fmax(bits_ab, bits_ba)
        lmin, lmax = np.minimum(lens[lo], lens[hi]), np.maximum(lens[lo], lens[hi])
        bsr_lo = bits_pair / selfbits[lo]
        bsr_hi = bits_pair / selfbits[hi]

    F: dict[str, np.ndarray] = {}
    F["ident_local_min"] = np.fmin(id_ab, id_ba)
    F["ident_local_max"] = np.fmax(id_ab, id_ba)
    F["ident_ens_min"] = nid / lmax
    F["ident_ens_max"] = nid / lmin
    g = nw.global_identity(seqs, lo, hi, cfg["global_alignment"]["gap_open"],
                           cfg["global_alignment"]["gap_extend"], threads=threads)
    F["ident_global"] = g["identity"]
    F["cov_min"] = np.fmin(cov_lo, cov_hi)
    F["cov_max"] = np.fmax(cov_lo, cov_hi)
    F["bits_pair"] = bits_pair
    F["bsr_min"] = np.fmin(bsr_lo, bsr_hi)
    F["bsr_max"] = np.fmax(bsr_lo, bsr_hi)
    F["len_ratio"] = lmin / lmax

    # --- per-gene hit lists (directed edges, g as query) ---
    E = pd.DataFrame({
        "src": np.concatenate([lo, hi]), "dst": np.concatenate([hi, lo]),
        "bits": np.concatenate([bits_ab, bits_ba]), "ident": np.concatenate([id_ab, id_ba]),
        "pair": np.concatenate([np.arange(npairs), np.arange(npairs)]),
    })
    E["bsr"] = E["bits"].to_numpy() / selfbits[E["src"].to_numpy()]
    E["rank"] = E.groupby("src")["bits"].rank(method="min", ascending=False)
    hits = {"all": np.bincount(E["src"], minlength=n)}
    for lev in fcfg["bsr_levels"]:
        sel = E["bsr"].to_numpy() >= float(lev)
        hits[f"r{int(round(lev * 100)):02d}"] = np.bincount(E["src"][sel], minlength=n)
    # best and second-best partner per gene by (bits desc, ident desc, dst asc)
    Es = E.sort_values(["src", "bits", "ident", "dst"], ascending=[True, False, False, True],
                       kind="mergesort")
    pos = Es.groupby("src").cumcount().to_numpy()
    top1 = Es[pos == 0].set_index("src")
    top2 = Es[pos == 1].set_index("src")
    src = E["src"].to_numpy()
    t1_dst = top1["dst"].reindex(src).to_numpy()
    t1_id = top1["ident"].reindex(src).to_numpy()
    t2_id = top2["ident"].reindex(src).fillna(0.0).to_numpy()
    other = np.where(t1_dst == E["dst"].to_numpy(), t2_id, t1_id)
    E["other_ident"] = other
    E["other_gap"] = E["ident"].to_numpy() - other
    first, second = E.iloc[:npairs], E.iloc[npairs:]

    def mm(col):
        a, b = first[col].to_numpy(np.float64), second[col].to_numpy(np.float64)
        return np.fmin(a, b), np.fmax(a, b)

    def gene_mm(v):
        a, b = v[lo], v[hi]
        return np.minimum(a, b), np.maximum(a, b)

    for key, v in hits.items():
        F[f"hits_{key}_min"], F[f"hits_{key}_max"] = gene_mm(v)
    F["rank_min"], F["rank_max"] = mm("rank")
    F["rbh"] = (F["rank_max"] == 1).astype(np.int8)
    F["other_ident_min"], F["other_ident_max"] = mm("other_ident")
    F["other_gap_min"], F["other_gap_max"] = mm("other_gap")

    # --- shared third paralogs and mutual-coverage local family ---
    inter, deg = _intersections(n, lo, hi, lo, hi)
    F["shared_third"] = inter
    union = deg[lo] + deg[hi] - inter - 2
    F["shared_third_jacc"] = np.where(union > 0, inter / np.maximum(union, 1), 0.0)
    thr = float(fcfg["mutual_coverage_min"])
    link = (cov_lo >= thr) & (cov_hi >= thr)
    inter_mc, deg_mc = _intersections(n, lo[link], hi[link], lo, hi)
    F["mc_deg_min"], F["mc_deg_max"] = gene_mm(deg_mc)
    F["fam_mc_size"] = deg_mc[lo] + deg_mc[hi] - inter_mc + np.where(link, 0, 2)

    # --- per-gene sequence properties ---
    F["len_min"], F["len_max"] = lmin, lmax
    F["lowcomp_frac_min"], F["lowcomp_frac_max"] = gene_mm(np.asarray(lowcomp, np.float64))

    # --- QC ---
    med = E.assign(plen=lens[E["dst"].to_numpy()]).groupby("src")["plen"].median()
    med_len = np.full(n, np.nan)
    med_len[med.index.to_numpy()] = med.to_numpy()
    frag_gene = lens < float(fcfg["fragment_ratio"]) * med_len
    F["fragment_suspect"] = (frag_gene[lo] | frag_gene[hi]).astype(np.int8)
    F["near_identical"] = (F["ident_local_max"] >= float(fcfg["near_identical_min"])).astype(np.int8)
    F["n_dir_called"] = pairs["n_dir_called"].to_numpy(np.int8)

    out = pd.DataFrame(F, index=pairs.index)
    for f in FEATURES:
        if f.tier == "A" and f.name.endswith("_pct") and f.name[:-4] in BY_NAME:
            out[f.name] = out[f.name[:-4]].rank(method="average", pct=True)
    missing = [c for c in TIER_A_COLUMNS if c not in out.columns]
    extra = [c for c in out.columns if c not in BY_NAME]
    if missing or extra:
        raise RuntimeError(f"registry mismatch: missing={missing} extra={extra}")
    out = out[TIER_A_COLUMNS]
    out = out.astype({c: BY_NAME[c].dtype for c in TIER_A_COLUMNS})
    summary = {
        "n_pairs_missing_orientation": n_missing_orient,
        "n_pairs_missing_both_orientations": n_missing_both,
        "n_genes_missing_selfscore": n_missing_self,
        "n_genes_in_pairs": int(len(involved)),
        "nw_cells": int((lens[lo] * lens[hi]).sum()),
    }
    return out, summary
