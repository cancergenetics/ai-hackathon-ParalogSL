"""Candidate paralog pairs from the self-search (DESIGN v3 §3, algorithms_comparative change 3).

A pair {a, b} (a != b) is called when E' <= evalue_prime_max in EITHER search direction, with
E' = qlen * N_ref * 2**(-bits) and N_ref fixed in config/caller.yaml, so the call does not depend
on proteome size. Pairs are unordered; internally they are stored as (i, j) with i < j on the
internal protein index, and in outputs as (gene_a, gene_b) with gene_a < gene_b lexicographically.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# The formula implemented by evalue_prime(); tests assert it equals config/caller.yaml.
EVALUE_PRIME_FORMULA = "qlen * N_ref * 2**(-bits)"


def evalue_prime(qlen, bits, n_ref: float):
    qlen = np.asarray(qlen, dtype=np.float64)
    bits = np.asarray(bits, dtype=np.float64)
    return qlen * n_ref * np.exp2(-bits)


def call_pairs(hits: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, dict]:
    """Return unordered called pairs (columns i < j, n_dir_found, n_dir_called) and stats."""
    c = cfg["caller"]
    if c["evalue_prime_formula"] != EVALUE_PRIME_FORMULA:
        raise ValueError("config E' formula differs from the implemented one")
    if c.get("directions", "either") != "either":
        raise ValueError("only directions: either is implemented")
    h = hits[hits["query"] != hits["target"]]
    ep = evalue_prime(h["qlen"].to_numpy(), h["bits"].to_numpy(), float(c["N_ref"]))
    q = h["query"].to_numpy(np.int64)
    t = h["target"].to_numpy(np.int64)
    called = ep <= float(c["evalue_prime_max"])
    i, j = np.minimum(q, t), np.maximum(q, t)
    d = pd.DataFrame({"i": i, "j": j, "q": q, "called": called})
    # one row per (direction, pair); a re-searched query cannot duplicate a direction, but be safe
    d = d.groupby(["i", "j", "q"], sort=False)["called"].max().reset_index()
    g = d.groupby(["i", "j"], sort=True).agg(n_dir_found=("q", "size"),
                                              n_dir_called=("called", "sum")).reset_index()
    pairs = g[g["n_dir_called"] > 0].reset_index(drop=True)
    pairs = pairs.astype({"i": np.int64, "j": np.int64, "n_dir_found": np.int8,
                          "n_dir_called": np.int8})
    n = len(pairs)
    stats = {
        "n_pairs": int(n),
        "n_pairs_called_one_direction": int((pairs["n_dir_called"] == 1).sum()),
        "n_pairs_found_one_direction": int((pairs["n_dir_found"] == 1).sum()),
        "frac_called_one_direction": float((pairs["n_dir_called"] == 1).mean()) if n else 0.0,
        "frac_found_one_direction": float((pairs["n_dir_found"] == 1).mean()) if n else 0.0,
    }
    return pairs, stats


def canonical_gene_order(ids_i: np.ndarray, ids_j: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(gene_a, gene_b) with gene_a < gene_b lexicographically."""
    a = np.where(ids_i < ids_j, ids_i, ids_j)
    b = np.where(ids_i < ids_j, ids_j, ids_i)
    return a, b
