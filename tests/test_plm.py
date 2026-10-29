"""PLM embedding: shape, determinism, chunking consistency, cache hits (ESM-2 8M, CPU fp32)."""

import hashlib
import json
import sys

import numpy as np
import pytest

from slparalog import plm

S1 = "MKTAYIAKQRQISFVKSHFSRQ"  # 22 aa
S2 = "GSHMLEDPVDAFKQLAEGHRA"  # 21 aa
S3 = "MSTNPKPQRKTKRNTNRRPQDVKFPGG"  # 27 aa
MODEL = "esm2_8m"


def _embed(seqs, cache, **kw):
    return plm.embed_sequences(seqs, cache, MODEL, device="cpu", log=lambda m: None, **kw)


def test_hash_is_sha256_of_canonical_sequence():
    assert plm.seq_sha256("mktay\n*") == hashlib.sha256(b"MKTAY").hexdigest()


def test_chunk_sequence_equal_chunks():
    assert [len(c) for c in plm.chunk_sequence("A" * 2044, 1022)] == [1022, 1022]
    assert [len(c) for c in plm.chunk_sequence("A" * 2045, 1022)] == [682, 682, 681]
    assert plm.chunk_sequence("ACD", 1022) == ["ACD"]


def test_shape_determinism_cache_and_chunking(tmp_path):
    spec = plm.model_spec(MODEL)
    assert len(spec["revision"]) == 40  # pinned commit hash
    seqs = [S1, S2, S3]
    c1 = tmp_path / "c1"
    m1, st1 = _embed(seqs, c1)
    assert m1.shape == (3, spec["dim"]) and m1.dtype == np.float32
    assert np.isfinite(m1).all()
    assert st1["n_computed"] == 3 and st1["n_cache_hits"] == 0
    assert not np.allclose(m1[0], m1[1])

    # determinism: fresh cache, recomputed -> identical
    m1b, st1b = _embed(seqs, tmp_path / "c1b")
    assert st1b["n_computed"] == 3
    np.testing.assert_array_equal(m1, m1b)

    # row order aligned to input order (reverse input -> reversed rows)
    mr, _ = _embed(seqs[::-1], tmp_path / "c1r")
    np.testing.assert_allclose(mr[::-1], m1, atol=2e-3)

    # cache hits: rerun computes nothing and is bitwise identical; duplicates computed once
    m2, st2 = _embed(seqs + [S2], c1)
    assert st2["n_computed"] == 0 and st2["n_cache_hits"] == 3 and st2["n_unique"] == 3
    assert "device" not in st2  # model never loaded
    np.testing.assert_array_equal(m2[:3], m1)
    np.testing.assert_array_equal(m2[3], m1[1])
    # cache stores float16 shards
    shards = sorted(c1.glob("shard_*.npy"))
    assert shards and all(np.load(s).dtype == np.float16 for s in shards)
    # torch-free cache reader returns the same matrix
    np.testing.assert_array_equal(plm.load_cached_matrix(seqs, c1, MODEL), m1)

    # partial cache: one new sequence -> only it is computed
    m3, st3 = _embed(seqs + [S3 + "AG"], c1)
    assert st3["n_computed"] == 1 and st3["n_cache_hits"] == 3
    np.testing.assert_array_equal(m3[:3], m1)

    # chunking consistency: S1+S1 with window = len(S1) splits at the boundary into two identical chunks,
    # whose length-weighted mean must equal the single-chunk embedding of S1.
    w = len(S1)
    a, sta = _embed([S1], tmp_path / "w_a", window=w)
    b, stb = _embed([S1 + S1], tmp_path / "w_b", window=w)
    assert sta["n_chunks"] == 1 and stb["n_chunks"] == 2
    np.testing.assert_allclose(a, b, atol=2e-3)
    # and chunking really changes things relative to embedding the doubled sequence whole
    whole, _ = _embed([S1 + S1], tmp_path / "w_c")
    assert not np.allclose(whole, b, atol=1e-4)


def test_cache_meta_mismatch_refused(tmp_path):
    _embed([S1], tmp_path / "c", window=40)
    with pytest.raises(ValueError, match="different settings"):
        _embed([S1], tmp_path / "c", window=41)


def test_empty_sequence_rejected(tmp_path):
    with pytest.raises(ValueError):
        _embed(["", S1], tmp_path / "c")


def test_embed_fasta_alignment_and_cli_subprocess(tmp_path):
    import subprocess

    fa = tmp_path / "p.fasta"
    fa.write_text(f">g3 desc\n{S3}\n>g1\n{S1[:11]}\n{S1[11:]}\n>g2\n{S2}\n")
    out = tmp_path / "emb.npy"
    r = subprocess.run(
        [sys.executable, "-m", "slparalog.cli", "embed", str(fa), "--cache", str(tmp_path / "cache"),
         "--out", str(out), "--model", MODEL, "--device", "cpu", "--lock", str(tmp_path / "lock"),
         "--bench-json", str(tmp_path / "b.json")],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    mat = np.load(out)
    assert out.with_suffix(".ids.txt").read_text().split() == ["g3", "g1", "g2"]
    ref, _ = _embed([S3, S1, S2], tmp_path / "ref")
    np.testing.assert_array_equal(mat, ref)  # multi-line FASTA parsed, record order kept
    b = json.loads((tmp_path / "b.json").read_text())
    assert b["n_records"] == 3 and b["peak_rss_gb"] > 0 and "lowpowermode" in b
