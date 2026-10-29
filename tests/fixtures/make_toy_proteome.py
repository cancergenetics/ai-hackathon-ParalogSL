"""Generate the toy proteome with planted duplicates (tests/fixtures/toy_proteome.fasta).

Deterministic (numpy seed 20261028). Planted structure, substitutions only (no indels), never in
the first/last 10 residues, each substitution to a different residue:
* toyA1/toyA2: 300 aa, 60 substitutions -> identity exactly 0.80 over the full length.
* toyF1/toyF2/toyF3: 3-member family, 400 aa. F2 and F3 each carry 80 substitutions relative to F1
  at DISJOINT positions -> F1-F2 = F1-F3 = 0.80, F2-F3 = 0.60.
* toyG1/toyG2: G2 = residues 101-250 of G1 (150 aa) with 15 substitutions -> fragment, local
  identity 0.90 over 150 aa; Ensembl-style 135/500 and 135/150.
* toyH1/toyH2: identical 200 aa copies.
* toyS01..toyS12: unrelated random singletons.
Run: python tests/fixtures/make_toy_proteome.py  (rewrites the FASTA; a test checks it is current).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
# approximate UniProt background frequencies (A C D E F G H I K L M N P Q R S T V W Y)
FREQ = np.array([8.25, 1.38, 5.46, 6.72, 3.86, 7.07, 2.27, 5.91, 5.80, 9.65,
                 2.41, 4.06, 4.74, 3.93, 5.53, 6.65, 5.36, 6.86, 1.10, 2.92])
FREQ = FREQ / FREQ.sum()
SEED = 20261028
OUT = Path(__file__).with_name("toy_proteome.fasta")


def _rand(rng, n):
    return "".join(rng.choice(list(AA), size=n, p=FREQ))


def _mutate(rng, seq, positions):
    s = list(seq)
    for p in positions:
        choices = [c for c in AA if c != s[p]]
        s[p] = choices[rng.integers(len(choices))]
    return "".join(s)


def build() -> list[tuple[str, str]]:
    rng = np.random.default_rng(SEED)
    recs = []
    a1 = _rand(rng, 300)
    pos = rng.choice(np.arange(10, 290), size=60, replace=False)
    recs += [("toyA1", a1), ("toyA2", _mutate(rng, a1, pos))]

    f1 = _rand(rng, 400)
    pos = rng.choice(np.arange(10, 390), size=160, replace=False)
    recs += [("toyF1", f1), ("toyF2", _mutate(rng, f1, pos[:80])),
             ("toyF3", _mutate(rng, f1, pos[80:]))]

    g1 = _rand(rng, 500)
    frag = g1[100:250]
    pos = rng.choice(np.arange(10, 140), size=15, replace=False)
    recs += [("toyG1", g1), ("toyG2", _mutate(rng, frag, pos))]

    h1 = _rand(rng, 200)
    recs += [("toyH1", h1), ("toyH2", h1)]

    for k in range(1, 13):
        recs.append((f"toyS{k:02d}", _rand(rng, int(rng.integers(150, 451)))))
    return recs


def render(recs) -> str:
    out = []
    for rid, seq in recs:
        out.append(f">{rid} planted toy record")
        out += [seq[i:i + 60] for i in range(0, len(seq), 60)]
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    OUT.write_text(render(build()))
    print(f"wrote {OUT}")
