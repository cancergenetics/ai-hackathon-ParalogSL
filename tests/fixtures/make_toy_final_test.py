"""Writes tests/fixtures/final_test_toy/: synthetic test-format labels and feature tables with hand-computable
answers for `scripts/eval/final_test.py --dry-run-toy` (tests/test_final_test_guard.py states the answers).

Pairs are listed with the rank order each method gives them; scores are descending ranks. Columns used as
stand-ins (all registered feature names): ident_local_max = B1, bits_pair = model, cov_min = B2,
len_min = yeast-only model (fam_mc_size and rbh are filler so that a real LR model can be fitted on the toy tables). No biology; run: .venv/bin/python tests/fixtures/make_toy_final_test.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).parent / "final_test_toy"


def scores(order: list[str], names: list[str]) -> dict[str, float]:
    n = len(order)
    s = {k: (n - i) / (n + 1) for i, k in enumerate(order)}
    return {k: s[k] for k in names}


def features(sp: str, names: list[str], genes: dict[str, tuple[str, str]], orders: dict[str, list[str]]) -> pd.DataFrame:
    """Stand-in columns plus filler for every other registered Tier A/B model feature, so that the real backend
    (a GBDT on A+B fitted by slparalog.train on these tables) can run on the toy repository."""
    from slparalog.features import registry

    rng = np.random.default_rng(len(names) + len(sp))
    rows = []
    for i, n in enumerate(names):
        row = {"gene_a": genes[n][0], "gene_b": genes[n][1], **{col: scores(orders[col], names)[n] for col in orders},
               "fam_mc_size": np.int32(2 + i % 3), "rbh": np.int32(i % 2)}
        for f in registry.MODEL_FEATURES:
            if registry.BY_NAME[f].tier in ("A", "B") and f not in row:
                row[f] = float(rng.random())
        row["ident_global"] = row["ident_local_max"]  # the real backend's B1 variant (ident_global) = the toy B1
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    (OUT / "features").mkdir(parents=True, exist_ok=True)
    cols = ["ident_local_max", "bits_pair", "cov_min", "len_min"]  # B1, model, B2, yeast-only stand-in

    # ---- human test pairs P0..P7 (called), P8 (labelled, not called), P9 (not in the primary tier)
    P = [f"P{i}" for i in range(10)]
    g = {p: (f"th_{p}a", f"th_{p}b") for p in P}
    y = {"P0": 1, "P1": 0, "P2": 1, "P3": 0, "P4": 0, "P5": 1, "P6": 0, "P7": 0, "P8": 1}  # P8: a caller-missed positive
    called = P[:8]
    orders = {
        "bits_pair": ["P0", "P2", "P5", "P1", "P3", "P4", "P6", "P7"],
        "ident_local_max": ["P0", "P1", "P2", "P3", "P4", "P5", "P6", "P7"],
        "cov_min": ["P1", "P0", "P3", "P2", "P5", "P4", "P6", "P7"],
        "len_min": ["P2", "P0", "P1", "P3", "P5", "P4", "P6", "P7"],
    }
    (OUT / "features" / "human").mkdir(exist_ok=True)
    features("human", called, g, orders).to_parquet(OUT / "features" / "human" / "tierA.parquet")
    nst = {"P0": 3, "P1": 3, "P2": 3, "P3": 3, "P4": 3, "P5": 2, "P6": 2, "P7": 2, "P8": 3, "P9": 3}
    studies = {"P0": "s1,s2", "P1": "s3,s1,s2", "P2": "s1,s2,s3", "P3": "s1,s2,s3", "P4": "s1,s2", "P5": "s4,s5",
               "P6": "s1,s2", "P7": "s2,s1", "P8": "s1,s2", "P9": "s1,s2"}
    h = pd.DataFrame({
        "gene_a": [g[p][0] for p in P], "gene_b": [g[p][1] for p in P],
        "label_primary": [y.get(p, float("nan")) for p in P],
        "n_studies_testable": [nst[p] for p in P], "studies": [studies[p] for p in P],
        "group": ["gh01" if p in ("P0", "P1") else f"gh{p}" for p in P],
        "n_testable": [5, 2, 6, 2, 1, 4, 2, 1, 3, 3],
        "harle_chan_dependent": ["True" if p == "P0" else "False" for p in P],
        "label_lenient": [1 if p in ("P0", "P2", "P4", "P5", "P8") else (float("nan") if p == "P9" else 0) for p in P],
        "label_two_studies": [1 if p in ("P0", "P2") else (float("nan") if p == "P9" else 0) for p in P],
        "leakage_bsr": [0.0, 0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.1, 0.9],
    })
    h.to_csv(OUT / "human_test.tsv", sep="\t", index=False)

    # ---- yeast test pairs Q0..Q7 (called), Q8 (labelled, not called)
    Q = [f"Q{i}" for i in range(9)]
    gq = {q: (f"ty_{q}a", f"ty_{q}b") for q in Q}
    yq = {q: int(q in ("Q1", "Q3")) for q in Q}
    qc = Q[:8]
    orders_y = {
        "bits_pair": ["Q1", "Q0", "Q3", "Q2", "Q4", "Q5", "Q6", "Q7"],
        "ident_local_max": ["Q0", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7"],
        "cov_min": ["Q1", "Q3", "Q0", "Q2", "Q4", "Q5", "Q6", "Q7"],
        "len_min": ["Q7", "Q6", "Q5", "Q4", "Q3", "Q2", "Q1", "Q0"],
    }
    (OUT / "features" / "yeast").mkdir(exist_ok=True)
    features("yeast", qc, gq, orders_y).to_parquet(OUT / "features" / "yeast" / "tierA.parquet")
    pd.DataFrame({"gene_a": [gq[q][0] for q in Q], "gene_b": [gq[q][1] for q in Q],
                  "primary": ["pos" if yq[q] else "neg" for q in Q], "group": [f"gy{q}" for q in Q],
                  "lenient": ["pos" if q in ("Q1", "Q3", "Q5") else "neg" for q in Q],
                  "strict": ["pos" if q == "Q1" else "neg" for q in Q],
                  "orientation_conflict": [1 if q == "Q1" else 0 for q in Q],
                  "leakage_bsr": [0.0, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.0, 0.05]}
                 ).to_csv(OUT / "yeast_test.tsv", sep="\t", index=False)

    # ---- unseen: Arabidopsis-like A0..A6 (A0 and A1 share a gene; A5 strict NA; A6 not called)
    A = [f"A{i}" for i in range(7)]
    ga = {"A0": ("ta1", "ta2"), "A1": ("ta2", "ta3"), "A2": ("ta4", "ta5"), "A3": ("ta6", "ta7"),
          "A4": ("ta8", "ta9"), "A5": ("ta10", "ta11"), "A6": ("ta12", "ta13")}
    ac = A[:6]
    orders_a = {
        "bits_pair": ["A3", "A0", "A1", "A2", "A4", "A5"],
        "ident_local_max": ["A1", "A0", "A2", "A3", "A4", "A5"],
        "cov_min": ["A0", "A1", "A2", "A3", "A4", "A5"],
        "len_min": ["A0", "A1", "A2", "A3", "A4", "A5"],
    }
    (OUT / "features" / "arabidopsis").mkdir(exist_ok=True)
    features("arabidopsis", ac, ga, orders_a).to_parquet(OUT / "features" / "arabidopsis" / "tierA.parquet")
    strict = {"A0": 1, "A1": 0, "A2": 0, "A3": 1, "A4": 0, "A5": "NA", "A6": 0}
    broad = {"A0": 1, "A1": 0, "A2": 0, "A3": 1, "A4": 1, "A5": 0, "A6": 0}
    pd.DataFrame({"gene_a": [ga[a][0] for a in A], "gene_b": [ga[a][1] for a in A],
                  "label_strict": [strict[a] for a in A], "label_broad": [broad[a] for a in A],
                  "single_lethal_n": [1 if a == "A0" else 0 for a in A], "note": ""}).to_csv(OUT / "unseen_arabidopsis.tsv", sep="\t", index=False)

    # ---- unseen: worm-like W0..W5 (W5 not called); B1 constant (all tied)
    W = [f"W{i}" for i in range(6)]
    gw = {w: (f"tw_{w}a", f"tw_{w}b") for w in W}
    yw = {"W0": 1, "W1": 0, "W2": 1, "W3": 0, "W4": 0, "W5": 0}
    wc = W[:5]
    orders_w = {
        "bits_pair": ["W0", "W2", "W1", "W3", "W4"],
        "ident_local_max": ["W0", "W1", "W2", "W3", "W4"],
        "cov_min": ["W1", "W0", "W3", "W2", "W4"],
        "len_min": ["W0", "W1", "W2", "W3", "W4"],
    }
    fw = features("worm", wc, gw, orders_w)
    fw["ident_local_max"] = 0.5
    (OUT / "features" / "worm").mkdir(exist_ok=True)
    fw.to_parquet(OUT / "features" / "worm" / "tierA.parquet")
    pd.DataFrame({"gene_a": [gw[w][0] for w in W], "gene_b": [gw[w][1] for w in W],
                  "label": [yw[w] for w in W]}).to_csv(OUT / "unseen_worm.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
