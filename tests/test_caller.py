"""Paralog caller: E' definition (fixed N_ref, either direction), canonical pair order."""

import numpy as np
import pandas as pd

from slparalog import caller
from slparalog.config import CONFIG_DIR, load_yaml

CFG = load_yaml(CONFIG_DIR / "caller.yaml")


def test_config_formula_is_the_implemented_one():
    assert CFG["caller"]["evalue_prime_formula"] == caller.EVALUE_PRIME_FORMULA
    assert float(CFG["caller"]["N_ref"]) == 1.1e7
    assert float(CFG["caller"]["evalue_prime_max"]) == 1e-3
    assert CFG["mmseqs"]["sensitivity"] == 7.5 and CFG["mmseqs"]["max_seqs"] == 1000
    assert CFG["mmseqs"]["backtrace"] is True


def test_evalue_prime_values():
    # qlen * N_ref * 2^-bits
    assert np.isclose(caller.evalue_prime(100, 10, 1e6), 100 * 1e6 / 1024)
    # threshold bits for qlen 300 at N_ref 1.1e7: log2(300 * 1.1e7 / 1e-3)
    b = np.log2(300 * 1.1e7 / 1e-3)
    assert np.isclose(caller.evalue_prime(300, b, 1.1e7), 1e-3)


def _hits(rows):
    return pd.DataFrame(rows, columns=["query", "target", "bits", "qlen"])


def test_either_direction_and_self_excluded():
    b_hi, b_lo = 60.0, 30.0  # 60 bits passes for qlen 100 (needs ~40), 30 does not
    hits = _hits([
        (0, 0, 500.0, 100),          # self hit: ignored
        (0, 1, b_hi, 100), (1, 0, b_lo, 100),   # called in one direction only
        (2, 3, b_hi, 100), (3, 2, b_hi, 100),   # both directions
        (4, 5, b_lo, 100), (5, 4, b_lo, 100),   # found, never called
        (6, 7, b_hi, 100),                      # found and called in one direction only
    ])
    pairs, st = caller.call_pairs(hits, CFG)
    got = {(r.i, r.j): (r.n_dir_found, r.n_dir_called) for r in pairs.itertuples()}
    assert got == {(0, 1): (2, 1), (2, 3): (2, 2), (6, 7): (1, 1)}
    assert st["n_pairs"] == 3
    assert st["n_pairs_called_one_direction"] == 2
    assert st["n_pairs_found_one_direction"] == 1


def test_call_is_independent_of_database_size():
    # no database-size term: the same hit is called identically whatever else is in the table
    base = _hits([(0, 1, 45.0, 200)])
    big = pd.concat([base, _hits([(k, k + 1, 20.0, 200) for k in range(2, 2000, 2)])])
    p1, _ = caller.call_pairs(base, CFG)
    p2, _ = caller.call_pairs(big, CFG)
    assert p1[["i", "j"]].values.tolist() == p2[["i", "j"]].values.tolist()


def test_canonical_gene_order():
    a, b = caller.canonical_gene_order(np.array(["g2", "g1"], dtype=object),
                                       np.array(["g1", "g3"], dtype=object))
    assert list(a) == ["g1", "g1"] and list(b) == ["g2", "g3"]
