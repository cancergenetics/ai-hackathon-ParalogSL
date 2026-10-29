"""Tier C (PLM embedding-pair) features on hand-built embeddings with known answers.

Known answers come from the geometric construction (unit vectors at chosen angles, proteome mean 0),
plus a brute-force float64 reference on random embeddings. Invariances: pair orientation, proteome
row permutation, block size, thread count, constant shift of all embeddings (centred features).
"""

import math

import numpy as np
import pandas as pd
import pytest

from slparalog.features import registry, tier_c
from slparalog.guards import FORBIDDEN_NAME_RE

TH = math.radians(20)
# 6 proteins in 3-d, mean exactly 0 => centred == raw
TOY = np.array([
    [1.0, 0.0, 0.0],                      # 0
    [math.cos(TH), math.sin(TH), 0.0],    # 1  (20 deg from 0)
    [-1.0, 0.0, 0.0],                     # 2  (= -0)
    [-math.cos(TH), -math.sin(TH), 0.0],  # 3  (= -1)
    [0.0, 0.0, 1.0],                      # 4
    [0.0, 0.0, -1.0],                     # 5
])
TOY_PAIRS = pd.DataFrame({"i": [0, 4, 2, 0], "j": [1, 0, 3, 2]})


def _run(emb, pairs, **kw):
    return tier_c.compute(pairs, emb, **kw)


def test_registry_tier_c_block():
    cols = registry.TIER_C_COLUMNS
    assert cols and all(registry.BY_NAME[c].tier == "C" for c in cols)
    assert all(registry.BY_NAME[c].source == "plm" for c in cols)
    assert not any(FORBIDDEN_NAME_RE.search(c) for c in cols)
    # pairwise only: no per-gene embedding coordinates / PCs
    assert not any(c.startswith(("tc_pc", "tc_dim", "tc_emb")) for c in cols)
    for c in ("tc_cos_pct", "tc_cos_c_pct", "tc_euc_pct", "tc_rank_max_pct"):
        assert c in cols


def test_known_answers_toy():
    df, s = _run(TOY, TOY_PAIRS)
    assert list(df.columns) == registry.TIER_C_COLUMNS
    c20 = math.cos(TH)
    # pair (0,1): nearest neighbours of each other
    r = df.iloc[0]
    assert r.tc_cos == pytest.approx(c20) and r.tc_cos_c == pytest.approx(c20)
    assert r.tc_euc == pytest.approx(math.sqrt(2 - 2 * c20))
    assert (r.tc_rank_min, r.tc_rank_max, r.tc_mnn) == (1, 1, 1)
    assert r.tc_above_frac_min == 0 and r.tc_above_frac_max == 0
    # best other for 0 (excl. 1) is 4/5 at cos 0; same for 1
    assert r.tc_gap_min == pytest.approx(c20) and r.tc_gap_max == pytest.approx(c20)
    # pair (4,0): for 4, all of 0..3 tie at 0 -> rank 1; for 0, only 1 is above 0 (5 ties) -> rank 2
    r = df.iloc[1]
    assert r.tc_cos == pytest.approx(0.0, abs=1e-12)
    assert (r.tc_rank_min, r.tc_rank_max, r.tc_mnn) == (1, 2, 0)
    assert r.tc_above_frac_min == 0 and r.tc_above_frac_max == pytest.approx(1 / 4)
    # gap for 0: 0 - cos20 ; for 4: 0 - 0 (ties with 1,2,3)
    assert r.tc_gap_min == pytest.approx(-c20) and r.tc_gap_max == pytest.approx(0.0, abs=1e-6)
    # pair (2,3): mirror of (0,1)
    assert (df.iloc[2].tc_rank_max, df.iloc[2].tc_mnn) == (1, 1)
    # pair (0,2): antipodal -> last of the 4 others for both
    r = df.iloc[3]
    assert r.tc_cos == pytest.approx(-1.0)
    assert (r.tc_rank_min, r.tc_rank_max) == (5, 5)
    assert r.tc_above_frac_min == 1.0 and r.tc_above_frac_max == 1.0
    assert r.tc_euc == pytest.approx(2.0)
    assert s["n_proteome"] == 6 and s["frac_mnn"] == pytest.approx(0.5)
    # percentiles: rank pct over the called pairs
    assert df["tc_rank_max_pct"].tolist() == pytest.approx([0.375, 0.75, 0.375, 1.0])


def test_rank_cap():
    df, s = _run(TOY, TOY_PAIRS, settings={"rank_cap_k": 2})
    assert df.iloc[3].tc_rank_min == 3 and df.iloc[3].tc_rank_max == 3   # capped at k+1
    assert df.iloc[3].tc_above_frac_min == 1.0                          # fraction is uncapped
    assert s["n_rank_capped_edges"] == 2


def test_centring_shift():
    # adding a constant vector changes raw cosine but not centred features / Euclidean distance
    shift = np.array([5.0, -3.0, 2.0])
    base, _ = _run(TOY, TOY_PAIRS)
    sh, _ = _run(TOY + shift, TOY_PAIRS)
    centred = [c for c in registry.TIER_C_COLUMNS if not c.startswith("tc_cos") or c.startswith("tc_cos_c")]
    pd.testing.assert_frame_equal(base[centred], sh[centred], atol=1e-6, check_exact=False)
    a, b = TOY[0] + shift, TOY[1] + shift
    assert sh.iloc[0].tc_cos == pytest.approx(a @ b / np.linalg.norm(a) / np.linalg.norm(b))
    assert sh.iloc[0].tc_cos != pytest.approx(base.iloc[0].tc_cos)


def test_identical_embeddings_tie():
    X = np.vstack([TOY, TOY[[0]], TOY[[2]]])   # 6 == 0 and 7 == 2 (mean stays 0)
    df, _ = _run(X, pd.DataFrame({"i": [0, 0], "j": [6, 1]}))
    c20 = math.cos(TH)
    r = df.iloc[0]   # identical copies tie at cos 1 -> mutual nearest neighbours
    assert (r.tc_rank_min, r.tc_rank_max, r.tc_mnn) == (1, 1, 1)
    assert r.tc_gap_min == pytest.approx(1 - c20, abs=1e-6)
    r = df.iloc[1]   # for 0, copy 6 is above 1 (rank 2); for 1, 0 and 6 tie (rank 1)
    assert (r.tc_rank_min, r.tc_rank_max, r.tc_mnn) == (1, 2, 0)


def _brute(X, lo, hi, k):
    Xc = X - X.mean(0)
    Xc = Xc / np.linalg.norm(Xc, axis=1, keepdims=True)
    C = Xc @ Xc.T
    n = len(X)
    out = {"rank": [], "frac": [], "gap": []}
    for a, b in zip(lo, hi):
        rs, fs, gs = [], [], []
        for g, p in ((a, b), (b, a)):
            others = [x for x in range(n) if x not in (g, p)]
            ab = sum(C[g, x] > C[g, p] + 1e-6 for x in others)
            rs.append(min(ab + 1, k + 1))
            fs.append(ab / (n - 2))
            gs.append(C[g, p] - max(C[g, x] for x in others))
        out["rank"].append(sorted(rs))
        out["frac"].append(sorted(fs))
        out["gap"].append(sorted(gs))
    return out


def test_random_vs_bruteforce_and_invariances():
    rng = np.random.default_rng(0)
    n, d = 60, 16
    X = rng.normal(size=(n, d)) + 2.0
    lo = rng.integers(0, n, 300)
    hi = rng.integers(0, n, 300)
    keep = lo != hi
    pairs = pd.DataFrame({"i": lo[keep], "j": hi[keep]}).drop_duplicates(ignore_index=True)
    a = np.minimum(pairs.i, pairs.j)
    b = np.maximum(pairs.i, pairs.j)
    pairs = pd.DataFrame({"i": a, "j": b}).drop_duplicates(ignore_index=True)
    k = 10
    df, _ = _run(X, pairs, settings={"rank_cap_k": k, "block_rows": 7})
    ref = _brute(X, pairs.i.to_numpy(), pairs.j.to_numpy(), k)
    np.testing.assert_array_equal(df[["tc_rank_min", "tc_rank_max"]].to_numpy(), np.array(ref["rank"]))
    np.testing.assert_allclose(df[["tc_above_frac_min", "tc_above_frac_max"]].to_numpy(), np.array(ref["frac"]))
    np.testing.assert_allclose(df[["tc_gap_min", "tc_gap_max"]].to_numpy(), np.array(ref["gap"]), atol=1e-5)
    Xu = X / np.linalg.norm(X, axis=1, keepdims=True)
    np.testing.assert_allclose(df.tc_cos, np.einsum("ij,ij->i", Xu[pairs.i], Xu[pairs.j]))
    np.testing.assert_allclose(df.tc_euc, np.linalg.norm(X[pairs.i] - X[pairs.j], axis=1))
    # every feature varies within the proteome
    for c in registry.TIER_C_COLUMNS:
        assert df[c].nunique() > 1, c

    # orientation swap
    sw, _ = _run(X, pairs.rename(columns={"i": "j", "j": "i"}), settings={"rank_cap_k": k})
    pd.testing.assert_frame_equal(df, sw, atol=1e-6, check_exact=False)
    # block size and thread count
    b1, _ = _run(X, pairs, settings={"rank_cap_k": k, "block_rows": 1}, threads=1)
    pd.testing.assert_frame_equal(df, b1, atol=1e-6, check_exact=False)
    # proteome row permutation (relabel indices)
    perm = rng.permutation(n)
    inv = np.argsort(perm)
    Xp = X[perm]
    pp = pd.DataFrame({"i": inv[pairs.i], "j": inv[pairs.j]})
    dp, _ = _run(Xp, pp, settings={"rank_cap_k": k})
    pd.testing.assert_frame_equal(df, dp, atol=1e-6, check_exact=False)


def test_rejects_bad_input():
    with pytest.raises(ValueError):
        _run(TOY, pd.DataFrame({"i": [0], "j": [0]}))
    with pytest.raises(ValueError):
        _run(TOY, pd.DataFrame({"i": [0], "j": [6]}))
    bad = TOY.copy()
    bad[0, 0] = np.nan
    with pytest.raises(ValueError):
        _run(bad, TOY_PAIRS)


def test_no_torch_import():
    import subprocess
    import sys

    code = ("import sys, numpy as np, pandas as pd; from slparalog.features import tier_c; "
            "tier_c.compute(pd.DataFrame({'i':[0],'j':[1]}), np.eye(4)); "
            "assert 'torch' not in sys.modules")
    subprocess.run([sys.executable, "-c", code], check=True)


def test_pipeline_hook_toy(tmp_path):
    """`predict --tier-c` path: embedding in a child process (dev model esm2_8m), then Tier C columns."""
    from pathlib import Path

    from slparalog import fasta, pipeline, plm

    toy = Path(__file__).parent / "fixtures" / "toy_proteome.fasta"
    cache = tmp_path / "cache"
    df, meta = pipeline.run(toy, tmp_path / "out.parquet", threads=1, features_only=True,
                            workdir=tmp_path / "work", tier_c=True, plm_cache=cache, plm_model="esm2_8m")
    assert set(registry.TIER_C_COLUMNS) <= set(df.columns)
    assert not df[registry.TIER_C_COLUMNS].isna().any().any()
    assert meta["tier_c"]["n_proteome"] == meta["n_proteins"]
    # same features from the cache-only (torch-free) matrix in canonical order
    seqs_by_id, _ = fasta.read_proteome(toy)
    ids = fasta.canonical_order(seqs_by_id)
    emb = plm.load_cached_matrix([seqs_by_id[k] for k in ids], cache, "esm2_8m")
    pos = {g: k for k, g in enumerate(ids)}
    pairs = pd.DataFrame({"i": df.gene_a.map(pos), "j": df.gene_b.map(pos)})
    ref, _ = tier_c.compute(pairs, emb)
    pd.testing.assert_frame_equal(df[registry.TIER_C_COLUMNS].reset_index(drop=True),
                                  ref.reset_index(drop=True), atol=1e-6, check_exact=False)


# --------------------------------------------------------------------------- run_tierC id -> row mapping
def _load_run_tierc():
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "scripts" / "features" / "run_tierC.py"
    spec = importlib.util.spec_from_file_location("run_tierC_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture()
def mapping_setup(tmp_path):
    """Toy proteome + synthetic embeddings in a sequence-hash cache (dev model key esm2_8m, d=320)."""
    import hashlib
    import json

    from slparalog import plm

    rng = np.random.default_rng(7)
    aa = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    ids = [f"g{k}" for k in range(9)]
    seqs = ["".join(rng.choice(aa, size=int(rng.integers(30, 80)))) for _ in ids]
    prot = tmp_path / "toy.fasta"
    prot.write_text("".join(f">{g} desc\n{s}\n" for g, s in zip(ids, seqs)))
    cfg = plm.load_plm_config()
    spec = plm.model_spec("esm2_8m", cfg)
    cache_dir = tmp_path / "cache"
    cache = plm.EmbeddingCache(cache_dir, plm._cache_meta(spec, cfg["pooling"]))
    by_gene = {g: rng.normal(size=spec["dim"]).astype(np.float16).astype(np.float32) for g in ids}
    # cache shard written in a different (shuffled) order than the FASTA
    order = rng.permutation(len(ids))
    cache.add_shard([plm.seq_sha256(seqs[k]) for k in order], np.stack([by_gene[ids[k]] for k in order]))
    npy = tmp_path / "toy.npy"
    np.save(npy, np.stack([by_gene[g] for g in ids]))
    idsf = tmp_path / "toy.ids.txt"
    idsf.write_text("\n".join(ids) + "\n")
    ta = pd.DataFrame({"gene_a": ["g0", "g1", "g2", "g4"], "gene_b": ["g3", "g7", "g5", "g8"]})
    ta_path = tmp_path / "tierA.parquet"
    ta.to_parquet(ta_path, index=False)
    sha = hashlib.sha256(prot.read_bytes()).hexdigest()
    (tmp_path / "tierA.parquet.meta.json").write_text(json.dumps({"input": {"sha256": sha}}))
    return {"prot": prot, "npy": npy, "idsf": idsf, "ta": ta_path, "cache": cache_dir,
            "by_gene": by_gene, "ids": ids, "pairs": ta}


def test_run_tierc_mapping_end_to_end(mapping_setup):
    m = _load_run_tierc()
    s = mapping_setup
    out, meta = m.build("toy", s["prot"], s["npy"], s["idsf"], s["ta"], s["cache"], "esm2_8m", threads=1)
    assert meta["cache_check"] == "identical"
    pd.testing.assert_frame_equal(out[["gene_a", "gene_b"]], s["pairs"])
    # per named pair, raw cosine / distance equal the embeddings looked up by gene name
    for _, r in out.iterrows():
        a, b = s["by_gene"][r.gene_a].astype(np.float64), s["by_gene"][r.gene_b].astype(np.float64)
        assert r.tc_cos == pytest.approx(a @ b / np.linalg.norm(a) / np.linalg.norm(b), abs=1e-12)
        assert r.tc_euc == pytest.approx(np.linalg.norm(a - b), abs=1e-9)
    # full frame equals compute() on an independently built name -> row mapping
    X = np.stack([s["by_gene"][g] for g in s["ids"]])
    pos = {g: k for k, g in enumerate(s["ids"])}
    ref, _ = tier_c.compute(pd.DataFrame({"i": s["pairs"].gene_a.map(pos), "j": s["pairs"].gene_b.map(pos)}), X)
    pd.testing.assert_frame_equal(out[registry.TIER_C_COLUMNS], ref, atol=1e-9, check_exact=False)


@pytest.mark.parametrize("how", ["roll", "swap2"])
def test_run_tierc_rejects_scrambled_npy(mapping_setup, how):
    m = _load_run_tierc()
    s = mapping_setup
    X = np.load(s["npy"])
    X = np.roll(X, 1, axis=0) if how == "roll" else X[[1, 0] + list(range(2, len(X)))]
    np.save(s["npy"], X)
    with pytest.raises(SystemExit, match="cache entries"):
        m.build("toy", s["prot"], s["npy"], s["idsf"], s["ta"], s["cache"], "esm2_8m", threads=1)


def test_run_tierc_rejects_scrambled_ids(mapping_setup):
    m = _load_run_tierc()
    s = mapping_setup
    ids = list(s["ids"])
    ids[3], ids[5] = ids[5], ids[3]
    s["idsf"].write_text("\n".join(ids) + "\n")
    with pytest.raises(SystemExit, match="ids differ"):
        m.build("toy", s["prot"], s["npy"], s["idsf"], s["ta"], s["cache"], "esm2_8m", threads=1)


def test_run_tierc_rejects_other_proteome(mapping_setup):
    import json

    m = _load_run_tierc()
    s = mapping_setup
    (s["ta"].parent / "tierA.parquet.meta.json").write_text(json.dumps({"input": {"sha256": "0" * 64}}))
    with pytest.raises(SystemExit, match="different proteome"):
        m.build("toy", s["prot"], s["npy"], s["idsf"], s["ta"], s["cache"], "esm2_8m", threads=1)
