"""Global (NW) identity kernel: matrix, known answers, scores vs Biopython, symmetry, threads."""

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

from slparalog.nw import ALPHABET, BLOSUM62, global_identity

AA = "ACDEFGHIKLMNPQRSTVWY"


def test_blosum62_matches_biopython():
    ref = substitution_matrices.load("BLOSUM62")
    for i, x in enumerate(ALPHABET):
        for j, y in enumerate(ALPHABET):
            assert BLOSUM62[i, j] == ref[x][y], (x, y)


def test_known_answers():
    s = "MKVLAAGIVGLLLAAPQSHA"
    r = global_identity([s, s], np.array([0]), np.array([1]))
    assert r["identity"][0] == 1.0 and r["columns"][0] == 20
    # one substitution, no gaps -> 19/20
    t = s[:5] + ("W" if s[5] != "W" else "Y") + s[6:]
    r = global_identity([s, t], np.array([0]), np.array([1]))
    assert r["identity"][0] == 19 / 20
    # prefix fragment: end gaps penalised; 10 identical columns over 20 columns
    r = global_identity([s, s[:10]], np.array([0]), np.array([1]))
    assert r["columns"][0] == 20 and r["matches"][0] == 10
    assert r["score"][0] == sum(BLOSUM62[ALPHABET.index(c), ALPHABET.index(c)] for c in s[:10]) - (11 + 10)


def test_x_never_identical():
    r = global_identity(["AXXA", "AXXA"], np.array([0]), np.array([1]))
    assert r["matches"][0] == 2 and r["columns"][0] == 4


def _random_pairs(seed=0, n=40):
    rng = np.random.default_rng(seed)
    seqs = ["".join(rng.choice(list(AA), int(rng.integers(5, 200)))) for _ in range(n)]
    for k in range(10):  # related pairs with substitutions and an indel
        s = list(seqs[k])
        for p in rng.choice(len(s), len(s) // 4, replace=False):
            s[p] = AA[rng.integers(20)]
        del s[2:4]
        seqs.append("".join(s))
    pi = np.concatenate([np.arange(n - 1), np.arange(10)])
    pj = np.concatenate([np.arange(1, n), np.arange(n, n + 10)])
    return seqs, pi, pj


def test_scores_equal_biopython_global_affine():
    seqs, pi, pj = _random_pairs()
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.mode = "global"
    al.open_gap_score = -12  # Biopython: first gap residue; = 11 + 1 x 1 in our convention
    al.extend_gap_score = -1
    r = global_identity(seqs, pi, pj)
    for k, (x, y) in enumerate(zip(pi, pj)):
        assert r["score"][k] == al.score(seqs[x], seqs[y])


def test_score_symmetric_and_threads_deterministic():
    seqs, pi, pj = _random_pairs(seed=1)
    a = global_identity(seqs, pi, pj, threads=1)
    b = global_identity(seqs, pj, pi, threads=1)
    c = global_identity(seqs, pi, pj, threads=4)
    assert (a["score"] == b["score"]).all()
    for key in ("score", "matches", "columns"):
        assert (a[key] == c[key]).all()
    assert np.array_equal(a["identity"], c["identity"])
