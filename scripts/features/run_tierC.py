"""Tier C (PLM embedding-pair) features for an existing Tier A run, aligned row-by-row to its pairs.

    .venv/bin/python scripts/features/run_tierC.py <species> [--threads 4] [--model esm2_150m]

Reads the cached proteome embeddings data/interim/plm/<model>/<species>.npy + .ids.txt (written by
``python -m slparalog.plm``, record order of data/interim/proteomes/<species>.fasta) and the called
pairs from data/interim/features/<species>/tierA.parquet. Checks that (i) the Tier A run was made on
the same proteome file, (ii) the embedding ids equal the proteome record ids, and (iii) the .npy rows
equal the sequence-hash cache entries of the current proteome sequences (unless --skip-cache-check).
No torch import: features are numpy/numba only (src/slparalog/features/tier_c.py).

Output: data/interim/features/<species>/tierC.parquet (gene_a, gene_b + registered Tier C columns in
Tier A row order) + .meta.json (input checksums, git head, settings, runtime, peak RSS), and
results/bench/tierC_<species>.json (runtime / memory).
"""

from __future__ import annotations

import argparse
import os
import sys

if __name__ == "__main__":  # thread caps must be set before numpy/numba load (not when imported by tests)
    _THREADS = "4"
    for _a, _v in zip(sys.argv, sys.argv[1:]):
        if _a == "--threads":
            _THREADS = _v
    for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS",
               "NUMBA_NUM_THREADS"):
        os.environ[_k] = _THREADS

import hashlib  # noqa: E402
import json  # noqa: E402
import platform  # noqa: E402
import resource  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from slparalog import __version__, pipeline, plm  # noqa: E402
from slparalog.features import registry, tier_c  # noqa: E402


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def peak_rss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round((r if sys.platform == "darwin" else r * 1024) / 1e9, 3)


def git_head() -> dict:
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "src", "scripts/features"], cwd=ROOT,
                               capture_output=True, text=True).stdout.strip()
        return {"head": head or None, "src_dirty": bool(dirty)}
    except OSError:
        return {"head": None, "src_dirty": None}


def _rel(p: Path) -> str:
    try:
        return str(Path(p).resolve().relative_to(ROOT))
    except ValueError:
        return str(p)


def build(species: str, prot: Path, npy: Path, idsf: Path, ta_path: Path, cache_dir: Path | None,
          model_key: str | None = None, threads: int = 4) -> tuple[pd.DataFrame, dict]:
    """Checked id -> row mapping + Tier C features in Tier A row order. ``cache_dir=None`` skips the
    sequence-hash check of the .npy rows against the embedding cache (then only ids are checked)."""
    t0 = time.time()
    spec = plm.model_spec(model_key)
    timings: dict[str, float] = {}
    ta_meta = json.loads(Path(str(ta_path) + ".meta.json").read_text())
    prot_sha = sha256(prot)
    if ta_meta["input"]["sha256"] != prot_sha:
        raise SystemExit("Tier A was computed on a different proteome file")
    ids = Path(idsf).read_text().split("\n")[:-1]
    fids, fseqs = plm.read_fasta(prot)
    if ids != fids:
        raise SystemExit("embedding ids differ from the proteome record ids")
    if len(set(ids)) != len(ids):
        raise SystemExit("duplicate ids in the embedding id list")
    t = time.time()
    emb = np.load(npy)
    if emb.shape != (len(ids), int(spec["dim"])):
        raise SystemExit(f"embedding shape {emb.shape} != ({len(ids)}, {spec['dim']})")
    cache_check = "skipped"
    if cache_dir is not None:  # each row must be the cache entry keyed by sha256 of that record's sequence
        ref = plm.load_cached_matrix(fseqs, cache_dir, spec["key"])
        if not np.array_equal(ref, emb):
            raise SystemExit("embedding .npy rows differ from the cache entries of the current sequences")
        cache_check = "identical"
        del ref
    timings["load_check_s"] = round(time.time() - t, 2)

    pos = {g: k for k, g in enumerate(ids)}
    ta = pd.read_parquet(ta_path, columns=["gene_a", "gene_b"])
    i = ta["gene_a"].map(pos)
    j = ta["gene_b"].map(pos)
    if i.isna().any() or j.isna().any():
        raise SystemExit("pair genes missing from the embedding ids")
    pairs = pd.DataFrame({"i": i.to_numpy(np.int64), "j": j.to_numpy(np.int64)})
    t = time.time()
    feats, summary = tier_c.compute(pairs, emb, threads=threads)
    timings["features_s"] = round(time.time() - t, 2)
    out = pd.concat([ta[["gene_a", "gene_b"]].reset_index(drop=True), feats.reset_index(drop=True)], axis=1)
    # output rows must be the Tier A rows, in order
    if not (len(out) == len(ta)
            and (out["gene_a"].to_numpy() == ta["gene_a"].to_numpy()).all()
            and (out["gene_b"].to_numpy() == ta["gene_b"].to_numpy()).all()):
        raise SystemExit("Tier C rows are not aligned with Tier A rows")
    timings["total_s"] = round(time.time() - t0, 2)
    src = {_rel(f): sha256(f) for f in (ROOT / "src/slparalog/features/tier_c.py",
                                        ROOT / "src/slparalog/features/registry.py", Path(__file__))}
    meta = {
        "slparalog_version": __version__,
        "species": species,
        "registry_hash": registry.registry_hash(),
        "tier_c_columns": registry.TIER_C_COLUMNS,
        "git": git_head(),
        "source_sha256": src,
        "model": {k: spec[k] for k in ("key", "hf_id", "revision", "dim", "window")},
        "inputs": {
            "embeddings": {"path": _rel(npy), "sha256": sha256(npy)},
            "embedding_ids": {"path": _rel(idsf), "sha256": sha256(idsf)},
            "proteome": {"path": _rel(prot), "sha256": prot_sha},
            "tierA": {"path": _rel(ta_path), "sha256": sha256(ta_path),
                      "n_rows": int(len(ta)), "tierA_schema_hash": ta_meta.get("schema_hash")},
        },
        "cache_check": cache_check,
        "summary": summary,
        "threads": threads,
        "timings": timings,
    }
    return out, meta


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("species")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--model", default=None, help="key in config/plm.yaml (default: default_model)")
    ap.add_argument("--skip-cache-check", action="store_true")
    a = ap.parse_args()
    spec = plm.model_spec(a.model)
    pdir = ROOT / "data/interim/plm" / spec["key"]
    fdir = ROOT / "data/interim/features" / a.species
    out, meta = build(a.species, ROOT / "data/interim/proteomes" / f"{a.species}.fasta",
                      pdir / f"{a.species}.npy", pdir / f"{a.species}.ids.txt", fdir / "tierA.parquet",
                      None if a.skip_cache_check else pdir / "cache", spec["key"], a.threads)
    timings, summary = meta["timings"], meta["summary"]
    rss = peak_rss_gb()
    meta["peak_rss_gb"] = rss
    dst = fdir / "tierC.parquet"
    pipeline.write(out, dst, meta)
    bench = {
        "species": a.species,
        "step": "tier_c_features (from cached embeddings; excludes the embedding step itself)",
        "n_proteome": summary["n_proteome"],
        "n_pairs": summary["n_pairs"],
        "threads": a.threads,
        "wall_s": timings["total_s"],
        "features_s": timings["features_s"],
        "load_check_s": timings["load_check_s"],
        "peak_rss_gb": rss,
        "lowpowermode": plm.lowpower_state(),
        "loadavg": list(os.getloadavg()),
        "machine": {"platform": platform.platform(), "cpu_count": os.cpu_count()},
        "git_head_at_run": meta["git"]["head"],
        "output": _rel(dst),
        "embedding_bench": f"results/bench/plm_{spec['key']}_{a.species}.json",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    bp = ROOT / "results/bench" / f"tierC_{a.species}.json"
    bp.parent.mkdir(parents=True, exist_ok=True)
    bp.write_text(json.dumps(bench, indent=1) + "\n")
    print(json.dumps({"species": a.species, "n_pairs": len(out), "frac_mnn": summary["frac_mnn"],
                      "timings": timings, "peak_rss_gb": rss}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
