"""Caller + Tier A features on the planted toy proteome (tests/fixtures/make_toy_proteome.py).

Known answers come from the planted construction (substitution-only copies, disjoint mutation
sets), not from a previous run. Invariance: input order, ID renaming, A/B swap, thread count.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from slparalog import fasta, pipeline, search
from slparalog.config import CONFIG_DIR, load_yaml
from slparalog.features import registry, tier_a

FIX = Path(__file__).parent / "fixtures"
TOY = FIX / "toy_proteome.fasta"
CFG = load_yaml(CONFIG_DIR / "caller.yaml")

EXPECTED_PAIRS = {
    ("toyA1", "toyA2"), ("toyF1", "toyF2"), ("toyF1", "toyF3"), ("toyF2", "toyF3"),
    ("toyG1", "toyG2"), ("toyH1", "toyH2"),
}


def _features(df):
    return df.set_index(["gene_a", "gene_b"]).sort_index()


@pytest.fixture(scope="module")
def base(tmp_path_factory):
    out = tmp_path_factory.mktemp("base") / "toy.parquet"
    df, meta = pipeline.run(TOY, out, threads=1, features_only=True)
    return df, meta, out


def test_fixture_is_current():
    sys.path.insert(0, str(FIX))
    try:
        import make_toy_proteome as mk
    finally:
        sys.path.pop(0)
    assert TOY.read_text() == mk.render(mk.build())


def test_capped_rule_counts_called_hits_only():
    from slparalog.search import capped_queries

    # qlen 100: 60 bits passes E' <= 1e-3 (needs ~40), 30 bits does not
    rows = [(0, 0, 500.0), (0, 1, 60.0), (0, 2, 60.0), (0, 3, 30.0),
            (1, 1, 500.0), (1, 0, 60.0), (1, 2, 30.0), (1, 3, 30.0),
            (2, 2, 500.0)]
    hits = pd.DataFrame(rows, columns=["query", "target", "bits"]).assign(qlen=100)
    assert capped_queries(hits, CFG, 3) == {0}        # 3 called hits >= 0.95 x 3
    assert capped_queries(hits, CFG, 4) == set()      # 3 < 3.8; query 1 has 4 hits but 2 called
    assert capped_queries(hits, CFG, 2) == {0, 1}


def test_called_pairs(base):
    df, meta, _ = base
    assert set(zip(df.gene_a, df.gene_b)) == EXPECTED_PAIRS
    assert (df.gene_a < df.gene_b).all()
    assert meta["caller"]["n_pairs"] == 6
    assert meta["n_proteins"] == 21
    assert meta["saturation"]["n_capped_first_round"] == 0


def test_planted_identities(base):
    f = _features(base[0])
    exact = {
        ("toyA1", "toyA2"): 0.8, ("toyF1", "toyF2"): 0.8, ("toyF1", "toyF3"): 0.8,
        ("toyF2", "toyF3"): 0.6, ("toyH1", "toyH2"): 1.0,
    }
    for pair, ident in exact.items():
        for col in ("ident_local_min", "ident_local_max", "ident_ens_min", "ident_ens_max",
                    "ident_global"):
            assert f.loc[pair, col] == pytest.approx(ident, abs=1e-12), (pair, col)
        assert f.loc[pair, "cov_min"] == 1.0
    g = f.loc[("toyG1", "toyG2")]
    assert g.ident_local_max == pytest.approx(0.9, abs=1e-12)        # 135 / 150
    assert g.ident_ens_min == pytest.approx(135 / 500, abs=1e-12)
    assert g.ident_ens_max == pytest.approx(135 / 150, abs=1e-12)
    assert g.ident_global == pytest.approx(135 / 500, abs=1e-12)     # end gaps -> 500 columns
    assert g.cov_min == pytest.approx(0.3) and g.cov_max == 1.0
    assert g.len_ratio == pytest.approx(0.3)


def test_family_context_and_flags(base):
    f = _features(base[0])
    fam = [("toyF1", "toyF2"), ("toyF1", "toyF3"), ("toyF2", "toyF3")]
    for p in fam:
        assert f.loc[p, "shared_third"] == 1
        assert f.loc[p, "fam_mc_size"] == 3
        assert f.loc[p, "hits_all_min"] == f.loc[p, "hits_all_max"] == 2
    # F2-F3 (0.6) is never a best hit: F1 (0.8) is better for both
    assert f.loc[("toyF2", "toyF3"), "rbh"] == 0
    assert f.loc[("toyF2", "toyF3"), "rank_min"] == 2
    assert f.loc[("toyF2", "toyF3"), "other_ident_min"] == pytest.approx(0.8)
    assert f.loc[("toyF2", "toyF3"), "other_gap_max"] == pytest.approx(-0.2)
    assert f.loc[("toyF1", "toyF2"), "other_gap_max"] == pytest.approx(0.2)  # F2: 0.8 vs F3 0.6
    for p in [("toyA1", "toyA2"), ("toyG1", "toyG2"), ("toyH1", "toyH2")]:
        assert f.loc[p, "rbh"] == 1
        assert f.loc[p, "shared_third"] == 0
        assert f.loc[p, "fam_mc_size"] == 2
        assert f.loc[p, "other_ident_max"] == 0.0
        assert f.loc[p, "hits_all_max"] == 1
    assert f.loc[("toyG1", "toyG2"), "mc_deg_max"] == 0  # G1 coverage 0.3 < 0.5: no mc link
    assert f.loc[("toyA1", "toyA2"), "mc_deg_max"] == 1
    # QC flags
    assert f["near_identical"].to_dict() == {p: int(p == ("toyH1", "toyH2")) for p in f.index}
    assert f["fragment_suspect"].to_dict() == {p: int(p == ("toyG1", "toyG2")) for p in f.index}
    assert (f["n_dir_called"] == 2).all()
    assert f.loc[("toyH1", "toyH2"), "bsr_min"] == 1.0


def test_percentiles_over_called_pairs(base):
    f = _features(base[0])
    assert f["ident_global_pct"].between(0, 1).all()
    assert f.loc[("toyH1", "toyH2"), "ident_global_pct"] == 1.0
    assert f.loc[("toyG1", "toyG2"), "ident_global_pct"] == pytest.approx(1 / 6)


def test_output_metadata_and_schema(base):
    df, meta, out = base
    import pyarrow.parquet as pq

    md = json.loads(pq.read_schema(out).metadata[b"slparalog"])
    assert md["schema_hash"] == pipeline.schema_hash(CFG) == meta["schema_hash"]
    assert md["tool"]["mmseqs_version"]
    assert "sl_score" not in df.columns
    assert list(df.columns[2:]) == registry.TIER_A_COLUMNS
    assert json.loads(Path(str(out) + ".meta.json").read_text())["schema_hash"] == md["schema_hash"]


def _assert_same(a, b):
    pd.testing.assert_frame_equal(_features(a), _features(b), check_exact=True)


def test_input_order_invariance(base, tmp_path):
    recs = fasta.parse_fasta(TOY.read_text().splitlines(True))
    rng = np.random.default_rng(7)
    shuffled = [recs[k] for k in rng.permutation(len(recs))]
    p = tmp_path / "shuffled.fasta"
    fasta.write_fasta(p, shuffled)
    df, _ = pipeline.run(p, None, threads=1, features_only=True)
    _assert_same(base[0], df)


def test_id_renaming_invariance(base, tmp_path):
    recs = fasta.parse_fasta(TOY.read_text().splitlines(True))
    # reverse-sorting names: renaming also flips which gene is gene_a
    names = [r[0] for r in recs]
    new = {old: f"q{999 - k:03d}" for k, old in enumerate(sorted(names))}
    p = tmp_path / "renamed.fasta"
    fasta.write_fasta(p, [(new[r], s) for r, s in recs])
    df, _ = pipeline.run(p, None, threads=1, features_only=True)
    back = {v: k for k, v in new.items()}
    a = df.gene_a.map(back)
    b = df.gene_b.map(back)
    df = df.assign(gene_a=np.where(a < b, a, b), gene_b=np.where(a < b, b, a))
    _assert_same(base[0], df)


def test_thread_count_determinism(base):
    df, _ = pipeline.run(TOY, None, threads=4, features_only=True)
    _assert_same(base[0], df)


def test_saturation_research_path(base, tmp_path):
    # cap 2: queries with >= 0.95 x 2 hits passing E' (self + >= 1 called partner) are capped (D11);
    # singletons (self only passes) are not, whatever their prefilter count
    cfg = load_yaml(CONFIG_DIR / "caller.yaml")
    cfg["mmseqs"]["max_seqs"] = 2
    cfg["mmseqs"]["research_cap_factor"] = 50
    p = tmp_path / "caller_small_cap.yaml"
    p.write_text(yaml.safe_dump(cfg))
    df, meta = pipeline.run(TOY, None, threads=1, features_only=True, config_path=p)
    sat = meta["saturation"]
    assert sat["rule"] == "called_hits"
    assert sat["n_capped_first_round"] == 9  # A1 A2 F1 F2 F3 G1 G2 H1 H2
    assert sat["n_capped_final"] == 0
    # the old prefilter-count rule would also have flagged singletons (17 of 21 fill a cap of 2)
    assert sat["rounds"][0]["n_prefilter_saturated"] > 9
    assert len(sat["rounds"]) >= 2
    _assert_same(base[0], df)


@pytest.fixture(scope="module")
def toy_alignments(tmp_path_factory):
    wd = tmp_path_factory.mktemp("aln")
    seqs_by_id, _ = fasta.read_proteome(TOY)
    ids = fasta.canonical_order(seqs_by_id)
    seqs = [seqs_by_id[k] for k in ids]
    mm = search.MMseqs(binary=search.resolve_binary(CFG), threads=1, workdir=wd)
    fasta.write_fasta(wd / "input.fasta", ((str(k), s) for k, s in enumerate(seqs)))
    search.createdb(mm, wd / "input.fasta", "db")
    hits, _ = search.all_vs_all(mm, CFG, "db", len(seqs))
    n = len(seqs)
    q = np.concatenate([np.repeat(np.arange(n), n)])
    t = np.concatenate([np.tile(np.arange(n), n)])
    aln = search.realign(mm, CFG, "db", q, t)
    return ids, seqs, hits, aln


def test_backtrace_identity_is_exact(toy_alignments):
    # change 1: with -a, fident == nident / alnlen. MMseqs2 prints fident truncated (not rounded)
    # to 3 decimals, so 0 <= nident/alnlen - fident < 1e-3 (checked: 15/43 = 0.3488 -> 0.348).
    _, _, hits, aln = toy_alignments
    for df in (hits, aln):
        assert (df["nident"] > 0).all()
        diff = df["nident"] / df["alnlen"] - df["fident"]
        assert (diff > -1e-9).all() and (diff < 1e-3).all()


def test_realign_returns_every_requested_pair(toy_alignments):
    ids, seqs, _, aln = toy_alignments
    n = len(seqs)
    assert len(aln) == n * n  # every ordered pair incl. unrelated ones and self pairs (-e inf)
    selfrows = aln[aln["query"] == aln["target"]]
    assert (selfrows["nident"] == selfrows["qlen"]).all()


def test_feature_symmetry_under_swap(toy_alignments):
    ids, seqs, hits, aln = toy_alignments
    from slparalog import caller

    pairs, _ = caller.call_pairs(hits, CFG)
    lowcomp = np.linspace(0, 0.3, len(seqs))
    f1, _ = tier_a.compute(pairs, aln, seqs, lowcomp, CFG)
    swapped = pairs.rename(columns={"i": "j", "j": "i"})[["i", "j", "n_dir_found", "n_dir_called"]]
    f2, _ = tier_a.compute(swapped, aln, seqs, lowcomp, CFG)
    pd.testing.assert_frame_equal(f1, f2, check_exact=True)


def test_cli_smoke(tmp_path):
    out = tmp_path / "pairs.tsv"
    r = subprocess.run([sys.executable, "-m", "slparalog.cli", "predict", str(TOY), "-o", str(out),
                        "-t", "2"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 6 and df["sl_score"].isna().all()
    assert {"gene_a", "gene_b", "near_identical", "fragment_suspect"} <= set(df.columns)


def test_cli_relative_workdir(tmp_path):
    # mmseqs runs with the workdir as cwd; a relative --workdir must still work
    out = tmp_path / "pairs.tsv"
    r = subprocess.run([sys.executable, "-m", "slparalog.cli", "predict", str(TOY), "-o", str(out),
                        "-t", "2", "--workdir", "rel_work"], capture_output=True, text=True, cwd=tmp_path)
    assert r.returncode == 0, r.stderr
    assert len(pd.read_csv(out, sep="\t")) == 6 and (tmp_path / "rel_work").is_dir()


def test_resolve_binary_falls_back_to_path(tmp_path, monkeypatch):
    from slparalog import search

    fake = tmp_path / "mmseqs"
    fake.write_text("#!/bin/sh\n")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    assert search.resolve_binary({"mmseqs": {"binary": "bin/does_not_exist"}}) == fake
