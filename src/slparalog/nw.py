"""Global (Needleman-Wunsch, Gotoh affine-gap) identity for protein pairs.

parasail (the DESIGN's choice) has no macOS arm64 wheel on PyPI (checked 2026-10-28: sdist only,
which needs a C toolchain build of libparasail), so this is the documented fallback: a numba JIT
kernel, parallel over pairs (``prange``), linear memory per pair.

Method: the DP carries, for every cell and state (M/E/F), the number of identical columns and the
number of aligned (residue-residue) columns of the path that produced the cell's score. Ties are
broken deterministically (diagonal > horizontal gap > vertical gap; gap open >= gap extend), so the
result is the identity of one well-defined optimal global alignment, obtained without a traceback
matrix (O(min) memory: ~titin x titin needs < 1 MB). No banding is needed.

Scoring: BLOSUM62 (20 standard residues + X), gap of length k costs ``gap_open + k * gap_extend``
(BLAST/NCBI convention, 11/1), end gaps penalised like internal gaps (classic NW, as parasail
``nw``). Identity = identical columns / alignment columns, where alignment columns =
len(a) + len(b) - aligned columns. X never counts as identical.

Symmetry: callers pass pairs in a canonical sequence order (the pipeline's internal index order is
sorted by sequence), so the value does not depend on which gene is "a".
"""

from __future__ import annotations

import numpy as np
import numba
from numba import njit, prange

# Use numba's own thread pool, not OpenMP: keeps libomp out of this package so it can share a
# process with LightGBM/torch (AGENT_RULES 8: libomp conflicts). Results do not depend on it.
numba.config.THREADING_LAYER = "workqueue"

ALPHABET = "ARNDCQEGHILKMFPSTWYVX"
X_CODE = ALPHABET.index("X")

# NCBI BLOSUM62, rows/cols in ALPHABET order (20 standard residues + X). Verified against
# Biopython's BLOSUM62 in tests/test_nw.py.
_BLOSUM62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3 -1
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3 -1
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -2
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2 -1
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2 -1
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3 -1
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -1
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2 -1
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -2
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -1
 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1
"""
BLOSUM62 = np.array([[int(v) for v in r.split()] for r in _BLOSUM62_ROWS.strip().splitlines()],
                    dtype=np.int32)

_LUT = np.full(256, X_CODE, dtype=np.uint8)
for _k, _c in enumerate(ALPHABET):
    _LUT[ord(_c)] = _k


def encode(seq: str) -> np.ndarray:
    return _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


@njit(cache=True, nogil=True)
def _nw_one(a, b, S, go, ge, xcode):
    # payload of a path = matches * 2**32 + aligned columns (both < 2**31 for any protein)
    n = a.shape[0]
    m = b.shape[0]
    NEG = np.int32(-(1 << 29))
    ONE = np.int64(1)
    MATCH = np.int64((1 << 32) + 1)
    goe = np.int32(go + ge)
    H = np.empty(m + 1, np.int32)
    P = np.zeros(m + 1, np.int64)
    F = np.empty(m + 1, np.int32)
    FP = np.zeros(m + 1, np.int64)
    H[0] = 0
    for j in range(1, m + 1):
        H[j] = -(go + j * ge)
        F[j] = NEG
    for i in range(1, n + 1):
        ai = a[i - 1]
        row = S[ai]
        dH = H[0]
        dP = P[0]
        H[0] = -(go + i * ge)
        P[0] = 0
        E = NEG
        EP = np.int64(0)
        hl = H[0]
        pl = P[0]
        for j in range(1, m + 1):
            hj = H[j]
            pj = P[j]
            fo = hj - goe
            fe = F[j] - ge
            if fo >= fe:
                fj = fo
                fpj = pj
            else:
                fj = fe
                fpj = FP[j]
            F[j] = fj
            FP[j] = fpj
            eo = hl - goe
            ee = E - ge
            if eo >= ee:
                E = eo
                EP = pl
            else:
                E = ee
            bj = b[j - 1]
            best = dH + row[bj]
            if ai == bj and ai != xcode:
                bp = dP + MATCH
            else:
                bp = dP + ONE
            dH = hj
            dP = pj
            if E > best:
                best = E
                bp = EP
            if fj > best:
                best = fj
                bp = fpj
            H[j] = best
            P[j] = bp
            hl = best
            pl = bp
    return H[m], P[m] >> 32, P[m] & 0xFFFFFFFF


@njit(parallel=True, cache=True)
def _nw_batch(buf, offs, lens, pi, pj, S, go, ge, xcode, out_s, out_m, out_d):
    for k in prange(pi.shape[0]):
        x = pi[k]
        y = pj[k]
        a = buf[offs[x]: offs[x] + lens[x]]
        b = buf[offs[y]: offs[y] + lens[y]]
        s, mm, dd = _nw_one(a, b, S, go, ge, xcode)
        out_s[k] = s
        out_m[k] = mm
        out_d[k] = dd


def global_identity(seqs: list[str], pi: np.ndarray, pj: np.ndarray, gap_open: int = 11,
                    gap_extend: int = 1, threads: int = 1) -> dict[str, np.ndarray]:
    """NW identity for pairs (seqs[pi[k]], seqs[pj[k]]). Returns score, matches, columns, identity."""
    pi = np.ascontiguousarray(pi, dtype=np.int64)
    pj = np.ascontiguousarray(pj, dtype=np.int64)
    lens = np.array([len(s) for s in seqs], dtype=np.int64)
    offs = np.zeros(len(seqs), dtype=np.int64)
    if len(seqs) > 1:
        offs[1:] = np.cumsum(lens)[:-1]
    buf = encode("".join(seqs)) if seqs else np.zeros(0, np.uint8)
    npairs = pi.shape[0]
    out_s = np.zeros(npairs, np.int32)
    out_m = np.zeros(npairs, np.int32)
    out_d = np.zeros(npairs, np.int32)
    if npairs:
        # interleave work so static prange chunks get similar total DP cells (timing only;
        # results are per-pair and independent of the order and of the thread count)
        cells = lens[pi] * lens[pj]
        order = np.argsort(-cells, kind="stable")
        nthreads = max(1, min(int(threads), numba.config.NUMBA_NUM_THREADS))
        perm = np.concatenate([order[t::nthreads] for t in range(nthreads)])
        prev = numba.get_num_threads()
        numba.set_num_threads(nthreads)
        try:
            s, m_, d = (np.zeros(npairs, np.int32) for _ in range(3))
            _nw_batch(buf, offs, lens, pi[perm], pj[perm], BLOSUM62, np.int32(gap_open),
                      np.int32(gap_extend), np.uint8(X_CODE), s, m_, d)
        finally:
            numba.set_num_threads(prev)
        out_s[perm], out_m[perm], out_d[perm] = s, m_, d
    cols = lens[pi] + lens[pj] - out_d
    ident = np.where(cols > 0, out_m / np.maximum(cols, 1), np.nan)
    return {"score": out_s, "matches": out_m, "columns": cols, "identity": ident}
