"""Toy known-answer tests for the robustness-test proteome degrader (scripts/robustness/degrade_proteome.py)."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("degrade_proteome", ROOT / "scripts" / "robustness" / "degrade_proteome.py")
D = importlib.util.module_from_spec(_spec)
sys.modules["degrade_proteome"] = D
_spec.loader.exec_module(D)

N = 400


def toy(n=N, seed=0):
    rng = np.random.default_rng(seed)
    return {f"g{i:04d}": "".join(rng.choice(list(D.AA), size=int(rng.integers(60, 300)))) for i in range(n)}


PROT = {f"g{i:04d}" for i in range(0, N, 5)}  # 20% "labelled"
SEED = 11


def run(kind, recs=None, seed=SEED):
    recs = recs or toy()
    return recs, *D.degrade(recs, kind, seed, PROT, species="toy")


def by_ops(rows, prefix):
    return [r for r in rows if r["role"] == "gene" and any(o.startswith(prefix) for o in r["ops"].split(";"))]


def test_fragments_fraction_length_and_end():
    recs, out, rows, c = run("fragments")
    assert c["n_fragments"] == round(0.2 * N) == 80
    assert len(out) == N and set(out) == set(recs)
    frag = by_ops(rows, "fragment")
    assert len(frag) == 80
    ends = set()
    for r in frag:
        o, n = recs[r["orig_id"]], out[r["orig_id"]]
        assert 0.5 - 1e-9 <= len(n) / len(o) <= 0.8 + 1 / len(o)  # rounding of the target length
        assert o.startswith(n) or o.endswith(n)
        ends.add("C" if o.startswith(n) else "N")
        assert r["degraded_len"] == len(n)
    assert ends == {"N", "C"}  # both ends occur
    for g in set(recs) - {r["orig_id"] for r in frag}:
        assert out[g] == recs[g]


def test_split_fraction_pieces_and_ids():
    recs, out, rows, c = run("split")
    assert c["n_split"] == 20 and len(out) == N + 20
    for r in by_ops(rows, "split"):
        g = r["orig_id"]
        nid = g + D.SPLIT_SUFFIX
        keep, other = out[g], out[nid]
        o = recs[g]
        assert o in (keep + other, other + keep)         # a clean split into two pieces
        assert len(keep) >= len(other) >= D.MIN_PIECE    # longer piece keeps the gene ID, internal point
    derived = [r for r in rows if r["role"] == "derived_split_piece"]
    assert len(derived) == 20 and all(r["orig_id"] + D.SPLIT_SUFFIX == r["degraded_id"] for r in derived)


def test_split_tie_keeps_n_terminal_piece():
    recs = {"a": "A" * 10 + "C" * 10, "b": "D" * 50}
    # a has length 20 = 2 * MIN_PIECE: the only split point is 10, a tie
    st = D._State(recs)
    st.n0 = 40  # force round(0.05 * 40) = 2 genes: both
    D.op_split(st, SEED)
    assert st.seqs["a"] == "A" * 10 and st.seqs["a" + D.SPLIT_SUFFIX] == "C" * 10


def test_fused_donor_removed_unlabelled_and_appended():
    recs, out, rows, c = run("fused")
    assert c["n_fused_targets"] == c["n_fused_donors_removed"] == round(0.02 * N) == 8
    assert len(out) == N - 8
    for r in by_ops(rows, "fused:"):
        donor = r["ops"].split("donor=")[1]
        assert donor not in PROT and donor not in out
        assert out[r["orig_id"]] == recs[r["orig_id"]] + recs[donor]
    removed = [r for r in rows if r["role"] == "gene" and r["degraded_id"] == ""]
    assert len(removed) == 8 and all(r["ops"].startswith("fused_donor_removed") for r in removed)


def test_fused_needs_enough_unlabelled_donors():
    recs = toy(50)
    with pytest.raises(D.DegradeError):
        D.degrade(recs, "fused", SEED, set(recs), species="toy")  # everything labelled


def test_missing_fraction():
    recs, out, rows, c = run("missing")
    assert c["n_missing"] == 60 and len(out) == N - 60
    assert all(out[g] == recs[g] for g in out)
    gone = [r for r in rows if r["degraded_id"] == ""]
    assert len(gone) == 60 and set(r["orig_id"] for r in gone) == set(recs) - set(out)


def test_haplotig_copy_near_identical():
    recs, out, rows, c = run("haplotig")
    assert c["n_haplotig_copies"] == 20 and len(out) == N + 20
    subs = tot = 0
    for r in rows:
        if r["role"] != "derived_haplotig":
            continue
        o, cp = recs[r["orig_id"]], out[r["degraded_id"]]
        assert r["degraded_id"] == r["orig_id"] + D.HAP_SUFFIX and len(o) == len(cp) and out[r["orig_id"]] == o
        subs += sum(a != b for a, b in zip(o, cp))
        tot += len(o)
    assert 0.003 < subs / tot < 0.03   # ~1% substitutions in expectation (20 genes, ~3.5k residues)


def test_mutate_rate_and_always_different_residue():
    rng = np.random.default_rng(1)
    s = "ACDEFGHIKL" * 2000
    m, n = D.mutate(s, rng, 0.01)
    assert sum(a != b for a, b in zip(s, m)) == n          # every hit changes the residue
    assert 120 < n < 280                                    # 20000 * 0.01 = 200
    assert D.mutate("AAAA", rng, 0.0) == ("AAAA", 0)


def test_combined_disjoint_and_counts():
    recs, out, rows, c = run("combined")
    assert (c["n_missing"], c["n_fused_targets"], c["n_fused_donors_removed"], c["n_split"], c["n_fragments"],
            c["n_haplotig_copies"]) == (60, 8, 8, 20, 80, 20)
    gene_rows = [r for r in rows if r["role"] == "gene"]
    # missing / fused / split / fragments act on disjoint genes (a gene's ops contain at most one of them)
    for r in gene_rows:
        kinds = {o.split(":")[0] for o in r["ops"].split(";") if o and not o.startswith("haplotig")}
        assert len(kinds) <= 1
    assert len(out) == N - 60 - 8 + 20 + 20
    # labelled genes are never fusion donors
    assert not any(r["ops"].startswith("fused_donor") and r["orig_id"] in PROT for r in gene_rows)


def test_deterministic_seeded_and_order_independent():
    recs = toy()
    a = D.degrade(recs, "combined", SEED, PROT)
    b = D.degrade(dict(reversed(list(recs.items()))), "combined", SEED, PROT)
    key = lambda r: (r["orig_id"], r["role"])  # noqa: E731
    assert a[0] == b[0] and sorted(a[1], key=key) == sorted(b[1], key=key)  # same content; row order follows the file
    c = D.degrade(recs, "combined", 23, PROT)
    assert a[0] != c[0]
    assert D.degrade(recs, "combined", SEED, PROT)[0] == a[0]


def test_identity_kind_unchanged_and_mapping_identity():
    recs, out, rows, c = run("identity")
    assert out == recs and D.mapping_dict(rows) == {g: g for g in recs}


def test_mapping_dict_and_new_ids_unique_no_collision():
    recs, out, rows, _ = run("combined")
    m = D.mapping_dict(rows)
    assert set(m) == {r["orig_id"] for r in rows if r["role"] == "gene" and r["degraded_id"]}
    assert all(v in out for v in m.values())
    ids = [r["degraded_id"] for r in rows if r["degraded_id"]]
    assert len(ids) == len(set(ids)) == len(out)
    with pytest.raises(D.DegradeError):
        recs2 = {f"h{i}": "A" * 40 for i in range(10)}
        recs2.update({f"h{i}" + D.SPLIT_SUFFIX: "C" * 15 for i in range(10)})  # too short to split: every eligible gene collides
        D.degrade(recs2, "split", SEED, set())


def test_write_kind_roundtrip(tmp_path, monkeypatch):
    from slparalog import fasta

    recs = toy(100)
    src = tmp_path / "p.fasta"
    fasta.write_fasta(src, recs.items())
    meta = D.write_kind("toy", "split", src, tmp_path / "out", SEED, set())
    back, _ = fasta.read_proteome(tmp_path / "out" / "toy" / "split" / "proteome.fasta")
    assert len(back) == 100 + 5 and meta["counts"]["n_split"] == 5
    import pandas as pd
    mp = pd.read_csv(tmp_path / "out" / "toy" / "split" / "mapping.tsv", sep="\t", keep_default_na=False)
    assert list(mp.columns) == D.MAPPING_COLUMNS and len(mp) == 105
