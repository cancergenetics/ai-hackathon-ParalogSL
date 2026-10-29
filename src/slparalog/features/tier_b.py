"""Tier B features: reference-panel duplication depth, conservation, species rank, asymmetry
(DESIGN v3 §4 items 4-7; reviews/round2/algorithms_comparative.md changes 5, 6, 10, 11, normative).

Inputs are indexed by the pipeline's internal protein index. Every pair feature is computed on the
canonical orientation (lo < hi) from symmetric rules (both members enter every condition in the
same way, ties are resolved without reference to which member is "A"), so f(a, b) == f(b, a).

Duplication classification of retained species s for pair (A, B), with a = A's best hit in s and
b = B's best hit in s (presence: bits >= presence_min_bits; bitscores are query-side and come from
searches with identical scoring settings):

* absent: neither A nor B has a hit;
* pre:  a != b, bits(A,a) > bits(A,B) and bits(B,b) > bits(B,A);
* post: bits(A,B) > bits(A,a) and bits(B,A) > bits(B,b) (a missing hit counts as 0 bits);
* undetermined: everything else (gene loss in s, rate asymmetry).

Depth on the d axis (d_s from panel.calibrate): L = max d over pre, U = min d over post;
L < U -> depth = (L+U)/2; no pre -> depth = U/2 and ``younger_than_nearest``; no post -> depth = L
and ``older_than_panel``; L >= U (conflict) -> depth = argmin_t of the bin-weighted (1/n_bin)
count of discordant species (pre with d > t plus post with d < t), t over the classified species'
d values (the cost is piecewise constant and minimal at one of them; ties -> lower median of the
tied values); the discordant fraction is that weighted count / weighted number of pre+post species.
Neither pre nor post -> depth NaN and ``tb_no_evidence`` = 1.

The bin edges (config/panel.yaml ``tier_b.bins``) are human-derived (frozen from human core-set d
values) and are used ONLY for the 1/n_bin discordance weights above.

Per-bin fractions are computed internally only (``bin_fractions``); they are not registered model
features (a bin-missingness pattern is constant within a proteome).
"""

from __future__ import annotations

import logging
import time

import numpy as np
import pandas as pd
from numba import njit, prange

from .. import panel as P
from .registry import BY_NAME, FEATURES, TIER_B_COLUMNS

CHUNK = 20000
log = logging.getLogger(__name__)


# ----------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------

def directional(aln: pd.DataFrame, lo: np.ndarray, hi: np.ndarray, n: int):
    """bits and identity for lo->hi and hi->lo from the realignment table; a missing orientation
    is filled from the other one (as in Tier A). Also self-bits per protein."""
    cols = ["query", "target", "bits", "nident", "alnlen"]
    a = (aln[cols].sort_values(["query", "target", "bits"], ascending=[True, True, False], kind="mergesort")
         .drop_duplicates(["query", "target"]))

    def orient(q, t):
        k = pd.DataFrame({"query": q, "target": t})
        m = k.merge(a, on=["query", "target"], how="left", sort=False)
        return m["bits"].to_numpy(np.float64), (m["nident"] / m["alnlen"]).to_numpy(np.float64)

    b_ab, i_ab = orient(lo, hi)
    b_ba, i_ba = orient(hi, lo)
    b_ab2 = np.where(np.isnan(b_ab), b_ba, b_ab)
    b_ba2 = np.where(np.isnan(b_ba), b_ab, b_ba)
    i_ab2 = np.where(np.isnan(i_ab), i_ba, i_ab)
    i_ba2 = np.where(np.isnan(i_ba), i_ab, i_ba)
    selfrows = aln[aln["query"] == aln["target"]]
    selfbits = np.full(n, np.nan)
    selfbits[selfrows["query"].to_numpy()] = selfrows["bits"].to_numpy(np.float64)
    return b_ab2, b_ba2, i_ab2, i_ba2, selfbits


def classify(bAa, bBb, a, b, bAB, bBA, pA, pB):
    """Per-species classification (module docstring). bAa/bBb: best-hit bits of A/B per species,
    a/b: best-hit targets, bAB/bBA: query-side pair bits (broadcastable), pA/pB: presence.
    Returns boolean arrays (pre, post, undetermined, absent)."""
    absent = ~pA & ~pB
    pre = pA & pB & (a != b) & (bAa > bAB) & (bBb > bBA)
    post = ~absent & (bAB > np.where(pA, bAa, 0)) & (bBA > np.where(pB, bBb, 0))
    undet = ~absent & ~pre & ~post
    return pre, post, undet, absent


@njit(cache=True)
def _conflict_depth(d, w, pre, post):
    """Depth and discordant fraction for one pair with L >= U (see module docstring)."""
    n = d.shape[0]
    cand = np.empty(n)
    m = 0
    tot = 0.0
    for s in range(n):
        if pre[s] or post[s]:
            cand[m] = d[s]
            m += 1
            tot += w[s]
    cand = np.sort(cand[:m])
    best = np.inf
    costs = np.empty(m)
    for k in range(m):
        t = cand[k]
        c = 0.0
        for s in range(n):
            if pre[s] and d[s] > t:
                c += w[s]
            elif post[s] and d[s] < t:
                c += w[s]
        costs[k] = c
        if c < best - 1e-12:
            best = c
    # tied minima (distinct values), lower median
    ties = np.empty(m)
    nt = 0
    for k in range(m):
        if costs[k] <= best + 1e-12 and (nt == 0 or cand[k] != ties[nt - 1]):
            ties[nt] = cand[k]
            nt += 1
    depth = ties[(nt - 1) // 2]
    return depth, best / tot if tot > 0 else np.nan


@njit(parallel=True, cache=True)
def _conflicts(rows, d, w, pre, post, out_depth, out_disc):
    for k in prange(rows.shape[0]):
        r = rows[k]
        dep, disc = _conflict_depth(d, w, pre[r], post[r])
        out_depth[r] = dep
        out_disc[r] = disc


@njit(parallel=True, cache=True)
def _shared_identity(qseq, qoff, oseq, ooff, cig, cigoff, qstart, tstart, aln_gene,
                     ta, tb, to, out_na, out_nb, out_nboth, out_full):
    """For triplet k (alignment ta[k] of gene A to O, alignment tb[k] of gene B to O, O index
    to[k]): identical columns of each alignment restricted to O positions covered by both
    alignments. MMseqs2 CIGAR: M = aligned column, I = query residue against a gap (query
    advances), D = target residue against a gap (target advances). out_full[k, 0/1] = identical
    columns over the whole alignment (checked against nident by the caller)."""
    for k in prange(ta.shape[0]):
        o = to[k]
        L = ooff[o + 1] - ooff[o]
        marks = np.zeros((2, L), np.int8)
        for side in range(2):
            a = ta[k] if side == 0 else tb[k]
            g = aln_gene[a]
            qp = qstart[a] - 1
            tp = tstart[a] - 1
            qL = qoff[g + 1] - qoff[g]
            num = 0
            full = 0
            for c in range(cigoff[a], cigoff[a + 1]):
                ch = cig[c]
                if ch >= 48 and ch <= 57:
                    num = num * 10 + (ch - 48)
                    continue
                if ch == 77:  # 'M'
                    for _ in range(num):
                        if tp < L and qp < qL:
                            same = qseq[qoff[g] + qp] == oseq[ooff[o] + tp]
                            marks[side, tp] = 2 if same else 1
                            if same:
                                full += 1
                        qp += 1
                        tp += 1
                elif ch == 73:  # 'I'
                    qp += num
                else:  # 'D'
                    tp += num
                num = 0
            out_full[k, side] = full
        na = 0
        nb = 0
        nboth = 0
        for p in range(L):
            if marks[0, p] > 0 and marks[1, p] > 0:
                nboth += 1
                if marks[0, p] == 2:
                    na += 1
                if marks[1, p] == 2:
                    nb += 1
        out_na[k] = na
        out_nb[k] = nb
        out_nboth[k] = nboth


def _encode(seqs: list[str]) -> tuple[np.ndarray, np.ndarray]:
    off = np.zeros(len(seqs) + 1, np.int64)
    np.cumsum([len(s) for s in seqs], out=off[1:])
    flat = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8).copy()
    return flat, off


# ----------------------------------------------------------------------------------------------
# main computation
# ----------------------------------------------------------------------------------------------

def compute(pairs: pd.DataFrame, aln: pd.DataFrame, seqs: list[str], pdata: P.PanelData,
            cal: P.Calibration, mask: P.Mask, pcfg: dict, threads: int = 1,
            return_internal: bool = False):
    """Tier B feature frame aligned with ``pairs`` rows (columns i, j: internal indices, either
    orientation). Returns (frame, summary[, internal dict with per-bin fractions])."""
    t0 = time.time()
    t = pcfg["tier_b"]
    n = len(seqs)
    lens = np.array([len(s) for s in seqs], np.int64)
    pi = pairs["i"].to_numpy(np.int64)
    pj = pairs["j"].to_numpy(np.int64)
    if np.any(pi == pj):
        raise ValueError("self pairs are not allowed")
    lo, hi = np.minimum(pi, pj), np.maximum(pi, pj)
    npairs = len(lo)
    b_ab, b_ba, i_ab, i_ba, selfbits = directional(aln, lo, hi, n)

    cols = mask.cols
    S = len(cols)
    d = mask.d.astype(np.float64)
    w = mask.weight.astype(np.float64)
    nbins = len(t["bins"]["edges"])
    presence = float(t["presence_min_bits"])
    cap = float(t["kimura"]["cap"])

    BB = pdata.best_bits[:, cols].astype(np.float32)          # n x S (0 = absent)
    BT = pdata.best_target[:, cols]
    SB = pdata.second_bits[:, cols]
    RB = pdata.rbh[:, cols]
    present = BB >= presence

    # --- per-gene quantities ---
    with np.errstate(invalid="ignore", divide="ignore"):
        cons_gene = present.mean(1) if S else np.full(n, np.nan)
        dc = np.where(present, pdata.dup_count[:, cols].astype(np.float64), np.nan)
        dup_gene = np.full(n, np.nan)
        anyp = present.any(1)
        if S:
            dup_gene[anyp] = np.nanmedian(dc[anyp], axis=1)
        tl = np.where(present, pdata.best_tlen[:, cols].astype(np.float64), np.nan)
        med_len = np.full(n, np.nan)
        if S:
            med_len[anyp] = np.nanmedian(tl[anyp], axis=1)
    frag_gene = lens < float(t["fragment_ratio"]) * med_len   # NaN median -> False

    F = {k: np.full(npairs, np.nan) for k in (
        "tb_depth", "tb_depth_range", "tb_discordant_frac", "tb_rel_rate", "tb_frac_pre",
        "tb_frac_post", "tb_frac_undet", "tb_frac_absent", "tb_cons_both",
        "tb_species_rank_min", "tb_species_rank_max", "tb_asym_norm", "tb_asym_abs",
        "tb_outgroup_ident_max")}
    for k in ("tb_younger_than_nearest", "tb_older_than_panel", "tb_no_evidence",
              "tb_n_post_single", "tb_n_outgroups"):
        F[k] = np.zeros(npairs, np.int64)
    binfrac = {c: np.full((npairs, nbins), np.nan) for c in ("pre", "post", "undet", "absent")} \
        if return_internal else None
    trip = []   # (pair, col_position, O, use_A_alignment_gene, ...) gathered per chunk
    n_conflict = 0
    d_pair = P.kimura(0.5 * (i_ab + i_ba), cap)

    for c0 in range(0, npairs, CHUNK):
        c1 = min(c0 + CHUNK, npairs)
        A, B = lo[c0:c1], hi[c0:c1]
        bAa, bBb = BB[A], BB[B]
        pA, pB = present[A], present[B]
        a, b = BT[A], BT[B]
        bAB = b_ab[c0:c1, None]
        bBA = b_ba[c0:c1, None]
        pre, post, undet, absent = classify(bAa, bBb, a, b, bAB, bBA, pA, pB)
        m = c1 - c0
        if S:
            F["tb_frac_pre"][c0:c1] = pre.mean(1)
            F["tb_frac_post"][c0:c1] = post.mean(1)
            F["tb_frac_undet"][c0:c1] = undet.mean(1)
            F["tb_frac_absent"][c0:c1] = absent.mean(1)
            F["tb_cons_both"][c0:c1] = (pA & pB).mean(1)
            srA = (pA & (bAa > bAB)).mean(1)
            srB = (pB & (bBb > bBA)).mean(1)
            F["tb_species_rank_min"][c0:c1] = np.minimum(srA, srB)
            F["tb_species_rank_max"][c0:c1] = np.maximum(srA, srB)
        has_pre, has_post = pre.any(1), post.any(1)
        L = np.where(has_pre, np.max(np.where(pre, d, -np.inf), axis=1, initial=-np.inf), np.nan)
        U = np.where(has_post, np.min(np.where(post, d, np.inf), axis=1, initial=np.inf), np.nan)
        depth = np.full(m, np.nan)
        rng = np.full(m, np.nan)
        disc = np.full(m, np.nan)
        both = has_pre & has_post
        ok = both & (L < U)
        depth[ok] = 0.5 * (L[ok] + U[ok])
        rng[both] = U[both] - L[both]
        disc[ok] = 0.0
        conf = np.flatnonzero(both & ~(L < U))
        n_conflict += len(conf)
        if len(conf):
            _conflicts(conf.astype(np.int64), d, w, pre, post, depth, disc)
        yn = ~has_pre & has_post
        depth[yn] = 0.5 * U[yn]
        rng[yn] = U[yn]
        disc[yn] = 0.0
        op = has_pre & ~has_post
        depth[op] = L[op]
        disc[op] = 0.0
        F["tb_depth"][c0:c1] = depth
        F["tb_depth_range"][c0:c1] = rng
        F["tb_discordant_frac"][c0:c1] = disc
        F["tb_younger_than_nearest"][c0:c1] = yn
        F["tb_older_than_panel"][c0:c1] = op
        F["tb_no_evidence"][c0:c1] = ~has_pre & ~has_post
        with np.errstate(invalid="ignore", divide="ignore"):
            F["tb_rel_rate"][c0:c1] = np.where(depth > 0, d_pair[c0:c1] / (2.0 * depth), np.nan)
        # n_post_single: post, a == b, and the member scoring higher to a (both on a tie) has
        # no second hit in s within ratio x its best bits
        r = float(t["n_post_single_ratio"])
        sA, sB = SB[A], SB[B]
        condA = sA < r * bAa
        condB = sB < r * bBb
        cond = np.where(bAa > bBb, condA, np.where(bBb > bAa, condB, condA & condB))
        F["tb_n_post_single"][c0:c1] = (post & pA & pB & (a == b) & cond).sum(1)
        if return_internal:
            bins = mask.bins
            for name, M in (("pre", pre), ("post", post), ("undet", undet), ("absent", absent)):
                for k in range(nbins):
                    sel = bins == k
                    if sel.any():
                        binfrac[name][c0:c1, k] = M[:, sel].mean(1)
        # outgroup triplets: the n nearest post species (columns are sorted by d)
        if S:
            kmax = int(t["outgroup"]["n_nearest_post"])
            rank = np.cumsum(post, axis=1)
            rr, cc = np.nonzero(post & (rank <= kmax))
            Ar, Br = A[rr], B[rr]
            rA, rB = RB[Ar, cc], RB[Br, cc]
            tA, tB = BT[Ar, cc], BT[Br, cc]
            bA, bB = BB[Ar, cc], BB[Br, cc]
            # O = RBH partner of A or B; if both have one, the one with higher bits (tie: lower O)
            useA = rA & (~rB | (bA > bB) | ((bA == bB) & (tA <= tB)))
            useB = rB & ~useA
            sel = useA | useB
            O = np.where(useA, tA, tB)
            trip.append(np.stack([rr[sel] + c0, cc[sel], O[sel]], axis=1))
    # --- divergence asymmetry over shared outgroup columns ---
    asym_stats = _asymmetry(trip, lo, hi, seqs, pdata, cols, pcfg, F, npairs, cap)

    F["tb_cons_min"] = np.minimum(cons_gene[lo], cons_gene[hi])
    F["tb_cons_max"] = np.maximum(cons_gene[lo], cons_gene[hi])
    F["tb_dupability_min"] = np.fmin(dup_gene[lo], dup_gene[hi])
    F["tb_dupability_max"] = np.fmax(dup_gene[lo], dup_gene[hi])
    F["tb_fragment_suspect"] = (frag_gene[lo] | frag_gene[hi]).astype(np.int8)

    out = pd.DataFrame(F, index=pairs.index)
    for f in FEATURES:
        if f.tier == "B" and f.name.endswith("_pct") and f.name[:-4] in BY_NAME:
            out[f.name] = out[f.name[:-4]].rank(method="average", pct=True)
    missing = [c for c in TIER_B_COLUMNS if c not in out.columns]
    extra = [c for c in out.columns if c not in BY_NAME]
    if missing or extra:
        raise RuntimeError(f"registry mismatch: missing={missing} extra={extra}")
    out = out[TIER_B_COLUMNS].astype({c: BY_NAME[c].dtype for c in TIER_B_COLUMNS})
    summary = {"n_pairs": int(npairs), "n_species_used": int(S), "n_conflict_pairs": int(n_conflict),
               "frac_no_evidence": float(out["tb_no_evidence"].mean()) if npairs else None,
               "frac_younger_than_nearest": float(out["tb_younger_than_nearest"].mean()) if npairs else None,
               "frac_older_than_panel": float(out["tb_older_than_panel"].mean()) if npairs else None,
               **asym_stats, "features_s": round(time.time() - t0, 2)}
    if return_internal:
        return out, summary, {"bin_fractions": binfrac, "d_pair": d_pair, "bits_lo_hi": b_ab,
                              "bits_hi_lo": b_ba, "lo": lo, "hi": hi}
    return out, summary


def _asymmetry(trip, lo, hi, seqs, pdata, cols, pcfg, F, npairs, cap):
    t = pcfg["tier_b"]
    min_cols = int(t["outgroup"]["min_shared_columns"])
    T = np.concatenate(trip, axis=0) if trip else np.zeros((0, 3), np.int64)
    if len(T) == 0:
        return {"n_triplets": 0, "n_triplets_valid": 0, "n_cigar_nident_mismatch": 0,
                "n_cigar_nident_off_by_one": 0}
    qflat, qoff = _encode(seqs)
    vals = np.full((len(T), 4), np.nan)   # asym_norm, asym_abs, ident_max, valid
    n_mismatch = 0
    n_off1 = 0
    n_missing_aln = 0
    for ci in np.unique(T[:, 1]):
        sel = np.flatnonzero(T[:, 1] == ci)
        code = pdata.codes[cols[ci]]
        pr, O = T[sel, 0], T[sel, 2]
        A, B = lo[pr], hi[pr]
        need = pd.DataFrame({"query": np.concatenate([A, B]), "target": np.concatenate([O, O])})
        need = need.drop_duplicates()
        h = P.load_hits(pdata.hit_paths[code],
                        ["query", "target", "bits", "nident", "qstart", "tstart", "cigar"])
        h = h.merge(need, on=["query", "target"], how="inner")
        h = (h.sort_values(["query", "target", "bits"], ascending=[True, True, False], kind="mergesort")
              .drop_duplicates(["query", "target"]).reset_index(drop=True))
        key = pd.Series(np.arange(len(h)), index=pd.MultiIndex.from_arrays([h["query"], h["target"]]))
        ia = key.reindex(pd.MultiIndex.from_arrays([A, O])).to_numpy()
        ib = key.reindex(pd.MultiIndex.from_arrays([B, O])).to_numpy()
        okk = ~np.isnan(ia) & ~np.isnan(ib)
        n_missing_aln += int((~okk).sum())
        if not okk.any():
            continue
        ia, ib = ia[okk].astype(np.int64), ib[okk].astype(np.int64)
        so = sel[okk]
        Ou = np.unique(O[okk])
        oseqs = P.read_species_seqs(pcfg, code, Ou)
        oloc = {int(x): k for k, x in enumerate(Ou)}
        oflat, ooff = _encode([oseqs[int(x)] for x in Ou])
        cig = h["cigar"].tolist()
        cflat, coff = _encode(cig)
        to = np.array([oloc[int(x)] for x in O[okk]], np.int64)
        m = len(ia)
        na, nb, nboth = np.zeros(m, np.int64), np.zeros(m, np.int64), np.zeros(m, np.int64)
        full = np.zeros((m, 2), np.int64)
        _shared_identity(qflat, qoff, oflat, ooff, cflat, coff, h["qstart"].to_numpy(np.int64),
                         h["tstart"].to_numpy(np.int64), h["query"].to_numpy(np.int64),
                         ia, ib, to, na, nb, nboth, full)
        nid = h["nident"].to_numpy(np.int64)
        # MMseqs2's nident is occasionally one lower than the identical columns of its own CIGAR
        # (coordinates agree); count those separately from real decoding mismatches
        dif = np.abs(np.concatenate([full[:, 0] - nid[ia], full[:, 1] - nid[ib]]))
        n_mismatch += int((dif > 1).sum())
        n_off1 += int((dif == 1).sum())
        good = nboth >= min_cols
        with np.errstate(invalid="ignore", divide="ignore"):
            idA = na / nboth
            idB = nb / nboth
            dA, dB = P.kimura(idA, cap), P.kimura(idB, cap)
            s = dA + dB
            norm = np.where(s > 0, np.abs(dA - dB) / s, 0.0)
        vals[so[good], 0] = norm[good]
        vals[so[good], 1] = np.abs(dA - dB)[good]
        vals[so[good], 2] = np.maximum(idA, idB)[good]
        vals[so[good], 3] = 1.0
    valid = vals[:, 3] == 1.0
    pr = T[valid, 0]
    V = vals[valid]
    if len(pr):
        df = pd.DataFrame({"p": pr, "n": V[:, 0], "a": V[:, 1], "i": V[:, 2]})
        g = df.groupby("p")
        med = g[["n", "a", "i"]].median()
        idx = med.index.to_numpy()
        F["tb_asym_norm"][idx] = med["n"].to_numpy()
        F["tb_asym_abs"][idx] = med["a"].to_numpy()
        F["tb_outgroup_ident_max"][idx] = med["i"].to_numpy()
        F["tb_n_outgroups"][idx] = g.size().to_numpy()
    return {"n_triplets": int(len(T)), "n_triplets_valid": int(valid.sum()),
            "n_triplets_missing_alignment": n_missing_aln,
            "n_cigar_nident_mismatch": n_mismatch, "n_cigar_nident_off_by_one": n_off1}


# ----------------------------------------------------------------------------------------------
# end-to-end stage (used by the pipeline hook and by scripts/features/run_tierB.py)
# ----------------------------------------------------------------------------------------------

def run_stage(mm, ccfg: dict, pcfg: dict, qdb: str, seqs: list[str], pairs: pd.DataFrame,
              aln: pd.DataFrame, threads: int = 1, panel_dmin: float | None = None,
              cache_root=None, species_jobs: int = 1) -> tuple[pd.DataFrame, dict]:
    """Panel search (cached) -> reduction -> calibration -> mask -> features. ``aln`` must hold
    both orientations of every pair and the self-alignments (``search.realign``)."""
    t0 = time.time()
    n = len(seqs)
    paths, smeta = P.search_panel(mm, ccfg, pcfg, qdb, seqs, cache_root=cache_root,
                                  species_jobs=species_jobs)
    t1 = time.time()
    lo = np.minimum(pairs["i"].to_numpy(np.int64), pairs["j"].to_numpy(np.int64))
    hi = np.maximum(pairs["i"].to_numpy(np.int64), pairs["j"].to_numpy(np.int64))
    _, _, _, _, selfbits = directional(aln, lo[:0], hi[:0], n)
    pdata = P.load_panel_data(paths, n, selfbits, pcfg)
    cal = P.calibrate(pdata, pcfg)
    mask = P.apply_mask(cal, pcfg, panel_dmin)
    qc = P.qc_flags(cal, pcfg)   # output metadata only, never a feature column
    for w in qc["warnings"]:
        log.warning("Tier B QC: %s", w)
    t2 = time.time()
    feats, fsum = compute(pairs, aln, seqs, pdata, cal, mask, pcfg, threads=threads)
    meta = {"settings_hash": P.settings_hash(pcfg, ccfg, smeta["db_stamp"]["panel_fasta_sha256"]),
            "panel_dmin": panel_dmin,
            "paths": {"db_dir": str(pcfg["tier_b"]["db_dir"]), "cache_dir": str(pcfg["tier_b"]["cache_dir"]),
                      "panel_fasta": str(pcfg["panel"]["fasta"])},  # provenance only, outside the hash
            "threads": int(threads), "species_jobs": int(species_jobs), "qc": qc,
            "search": smeta, "calibration": cal.log, "mask": mask.log, "features": fsum,
            "timings": {"panel_search_s": round(t1 - t0, 1), "reduce_calibrate_s": round(t2 - t1, 1),
                        "features_s": round(time.time() - t2, 1)}}
    return feats, meta
