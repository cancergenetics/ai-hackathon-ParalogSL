"""Toy query proteome + toy reference panel with planted duplication histories (Tier B tests).

Deterministic (numpy seed 20261029). Substitutions only (no indels), never in the first/last 10
residues, each to a different residue, and all mutation sets of one family DISJOINT, so every
identity below is exact over the full length (local alignments span the whole protein).

Panel species: S0 = identical copy of every query protein (must be self-excluded); S1, S2, S3 with
core-gene identities 0.90, 0.75, 0.60 -> d = Kimura(0.90) < Kimura(0.75) < Kimura(0.60).

Query proteins and their panel homologs (identities to the query copies):
* core K01..K30 (300 aa): one ortholog per species at 0.90 / 0.75 / 0.60.
* P1/P2 (pre-split duplication, 400 aa): P1 = R + 15%, P2 = R + 15% (id 0.70). S1: P1+4%, P2+4%;
  S2: P1+8%, P2+8% (pre in S1, S2); S3: one copy R+25% (id 0.60 < 0.70: post).
  -> L = d_S2, U = d_S3, depth = (d_S2 + d_S3)/2.
* Y1/Y2 (post-split, young): Y1 = R+1%, Y2 = R+4% (id 0.95); one copy per species R+10/20/35%
  -> post everywhere, younger_than_nearest, depth = d_S1/2; asymmetry from 3 outgroups.
* E1/E2 (old): E1 = R+25%, E2 = R+25% (id 0.50); S1/S2/S3 have E1+5/10/15% and E2+5/10/15%
  -> pre everywhere, older_than_panel, depth = d_S3.
* Q1/Q2 (one copy lost in S2): as P, but S2 keeps only Q1+8% -> S1 pre, S2 undetermined, S3 post
  -> depth = (d_S1 + d_S3)/2.
* Z1/Z2 (absent outside S1): as P in S1 (pre); no homolog in S2, S3 -> absent, older_than_panel.
* X1/X2 (conflict): id 0.70; S1 one copy R+25% (post); S3 X1+5%, X2+5% (pre); nothing in S2
  -> L = d_S3 >= U = d_S1; bin-weighted discordance ties at both -> depth = d_S1, discordant 1/2.
* F1/F2 (fragment): F1 = R (400 aa), F2 = residues 100-249 of R + 10%; every species has R+10%
  -> F2 is shorter than 0.5 x its panel best-hit median length (fragment_suspect).
* U01..U05: unrelated singletons.
"""

from __future__ import annotations

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
FREQ = np.array([8.25, 1.38, 5.46, 6.72, 3.86, 7.07, 2.27, 5.91, 5.80, 9.65,
                 2.41, 4.06, 4.74, 3.93, 5.53, 6.65, 5.36, 6.86, 1.10, 2.92])
FREQ = FREQ / FREQ.sum()
SEED = 20261029
SPECIES = ["S0", "S1", "S2", "S3"]
CORE_ID = {"S1": 0.90, "S2": 0.75, "S3": 0.60}


class _Fam:
    """A root sequence and a pool of unused mutable positions (keeps mutation sets disjoint)."""

    def __init__(self, rng, length):
        self.rng = rng
        self.root = "".join(rng.choice(list(AA), size=length, p=FREQ))
        self.pool = list(rng.permutation(np.arange(10, length - 10)))

    def branch(self):
        """Independent lineage: same root, a copy of the remaining pool (sister lineages may reuse
        positions; mutations along one lineage stay disjoint from everything it is compared to)."""
        b = _Fam.__new__(_Fam)
        b.rng, b.root, b.pool = self.rng, self.root, list(self.pool)
        return b

    def mut(self, seq, frac):
        k = int(round(frac * len(self.root)))
        pos, self.pool = self.pool[:k], self.pool[k:]
        assert len(pos) == k, "position pool exhausted"
        s = list(seq)
        for p in pos:
            choices = [c for c in AA if c != self.root[p]]
            s[p] = choices[self.rng.integers(len(choices))]
        return "".join(s)


def build():
    """Return (query_records, panel {species: [seq, ...]}, expected dict)."""
    rng = np.random.default_rng(SEED)
    query: list[tuple[str, str]] = []
    panel: dict[str, list[str]] = {s: [] for s in SPECIES}

    for k in range(1, 31):
        f = _Fam(rng, 300)
        query.append((f"K{k:02d}", f.root))
        for s, ident in CORE_ID.items():
            panel[s].append(f.branch().mut(f.root, 1 - ident))

    def pair_like_p(name, s1=True, s2="both", s3=True):
        f = _Fam(rng, 400)
        a, b = f.mut(f.root, 0.15), f.mut(f.root, 0.15)
        query.extend([(f"{name}1", a), (f"{name}2", b)])
        if s1:
            g = f.branch()
            panel["S1"] += [g.mut(a, 0.04), g.mut(b, 0.04)]
        if s2 == "both":
            g = f.branch()
            panel["S2"] += [g.mut(a, 0.08), g.mut(b, 0.08)]
        elif s2 == "one":
            panel["S2"] += [f.branch().mut(a, 0.08)]
        if s3:
            panel["S3"].append(f.branch().mut(f.root, 0.25))

    pair_like_p("P")
    pair_like_p("Q", s2="one")
    pair_like_p("Z", s2="none", s3=False)

    f = _Fam(rng, 400)
    y1, y2 = f.mut(f.root, 0.01), f.mut(f.root, 0.04)
    query.extend([("Y1", y1), ("Y2", y2)])
    for s, m in (("S1", 0.10), ("S2", 0.20), ("S3", 0.35)):
        panel[s].append(f.branch().mut(f.root, m))

    f = _Fam(rng, 400)
    e1, e2 = f.mut(f.root, 0.25), f.mut(f.root, 0.25)
    query.extend([("E1", e1), ("E2", e2)])
    for s, m in (("S1", 0.05), ("S2", 0.10), ("S3", 0.15)):
        g = f.branch()
        panel[s] += [g.mut(e1, m), g.mut(e2, m)]

    f = _Fam(rng, 400)
    x1, x2 = f.mut(f.root, 0.15), f.mut(f.root, 0.15)
    query.extend([("X1", x1), ("X2", x2)])
    panel["S1"].append(f.branch().mut(f.root, 0.25))
    g = f.branch()
    panel["S3"] += [g.mut(x1, 0.05), g.mut(x2, 0.05)]

    f = _Fam(rng, 400)
    frag = list(f.root[100:250])
    pos = rng.choice(np.arange(10, 140), size=15, replace=False)
    for p in pos:
        frag[p] = [c for c in AA if c != frag[p]][rng.integers(19)]
    query.extend([("F1", f.root), ("F2", "".join(frag))])
    for s in ("S1", "S2", "S3"):
        panel[s].append(f.branch().mut(f.root, 0.10))

    for k in range(1, 6):
        query.append((f"U{k:02d}", "".join(rng.choice(list(AA), size=250, p=FREQ))))

    panel["S0"] = [s for _, s in query]   # identical to the query: must be self-excluded

    def kim(i):
        p = 1 - i
        return float(-np.log(1 - p - 0.2 * p * p))

    d = {s: kim(i) for s, i in CORE_ID.items()}
    expected = {
        "d": d,
        "self_excluded": ["S0"],
        "P": {"depth": (d["S2"] + d["S3"]) / 2, "younger": 0, "older": 0,
              "pre": 2, "post": 1, "undet": 0, "absent": 0},
        "Y": {"depth": d["S1"] / 2, "younger": 1, "older": 0, "pre": 0, "post": 3},
        "E": {"depth": d["S3"], "younger": 0, "older": 1, "pre": 3, "post": 0},
        "Q": {"depth": (d["S1"] + d["S3"]) / 2, "younger": 0, "older": 0,
              "pre": 1, "post": 1, "undet": 1},
        "Z": {"depth": d["S1"], "older": 1, "absent": 2},
        "X": {"depth": d["S1"], "disc": 0.5},
        # Y outgroup identities in S1/S2/S3: Y1 = 1 - 0.01 - m, Y2 = 1 - 0.04 - m
        "Y_outgroup": [(1 - 0.01 - m, 1 - 0.04 - m) for m in (0.10, 0.20, 0.35)],
        "P_rel_rate": kim(0.70) / (2 * (d["S2"] + d["S3"]) / 2),
    }
    return query, panel, expected


def write(tmp, query, panel):
    """Write query FASTA, panel FASTA (headers '<code>|<index>') and species table under tmp."""
    from pathlib import Path

    tmp = Path(tmp)
    qf = tmp / "query.fasta"
    qf.write_text("".join(f">{n}\n{s}\n" for n, s in query))
    pf = tmp / "panel.fasta"
    with open(pf, "w") as fh:
        for code in SPECIES:
            for i, s in enumerate(panel[code]):
                fh.write(f">{code}|{i}\n{s}\n")
    st = tmp / "species.tsv"
    st.write_text("code\ttaxid\tname\n" + "".join(f"{c}\t{k + 1}\ttoy {c}\n" for k, c in enumerate(SPECIES)))
    return qf, pf, st
