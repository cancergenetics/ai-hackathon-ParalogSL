"""Tier B features for an existing Tier A run, aligned row-by-row to its pairs.

    .venv/bin/python scripts/features/run_tierB.py <species> [--threads 10] [--panel-dmin c]

Reads data/interim/proteomes/<species>.fasta and data/interim/features/<species>/tierA.parquet
(checks that the Tier A run was made on the same proteome), rebuilds the query DB in canonical
order, realigns the called pairs in both orientations plus self-alignments (``search.realign``,
identical to the pipeline), and runs the shipped Tier B stage (``tier_b.run_stage``: cached
per-species panel search, self-exclusion, d_s calibration, optional mask, features).

Output: data/interim/features/<species>/tierB[_dmin<c>].parquet with gene_a, gene_b and the
registered Tier B columns in the Tier A row order, plus <output>.meta.json (calibration log:
per-species d_s, RBH counts, core-set size, self-excluded species; search/saturation log; timings;
schema hash; sha256 of the Tier A file it is aligned to; threads, species_jobs; proteome-level QC
flags/warnings under ``qc`` -- output metadata only, never a feature column). The row alignment
(gene_a/gene_b -> canonical positions) lives in ``pair_positions``/``compute_tierB`` (tested).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from slparalog import __version__, fasta, panel, pipeline, search  # noqa: E402
from slparalog.config import CONFIG_DIR, load_yaml  # noqa: E402
from slparalog.features import registry, tier_b  # noqa: E402


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pair_positions(ids: list[str], gene_a, gene_b) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """Map Tier A ``gene_a``/``gene_b`` (rows in Tier A order) to positions in the canonical order
    ``ids``. Returns (i, j, pairs) with i/j the positions of gene_a/gene_b per row and ``pairs``
    the (lo, hi) = (min, max) position frame in the same row order. Raises on unknown genes."""
    pos = {g: k for k, g in enumerate(ids)}
    i = pd.Series(np.asarray(gene_a, dtype=object)).map(pos)
    j = pd.Series(np.asarray(gene_b, dtype=object)).map(pos)
    if i.isna().any() or j.isna().any():
        raise ValueError("Tier A pairs contain genes that are not in the proteome")
    i, j = i.to_numpy(np.int64), j.to_numpy(np.int64)
    return i, j, pd.DataFrame({"i": np.minimum(i, j), "j": np.maximum(i, j)})


def compute_tierB(seqs_by_id: dict[str, str], ta: pd.DataFrame, ccfg: dict, pcfg: dict, wd: Path,
                  threads: int = 1, panel_dmin: float | None = None, species_jobs: int = 1,
                  cache_root=None):
    """Tier B rows for the Tier A pairs ``ta`` (columns gene_a, gene_b), in ``ta`` row order.

    Rebuilds the canonical-order query DB in ``wd`` (must exist), realigns, runs
    ``tier_b.run_stage``. Returns (out, pairbits, meta, tool, timings): ``out`` = gene_a, gene_b +
    registered Tier B columns; ``pairbits`` = query-side pair bitscores sidecar. Factored out of
    ``main`` so the row alignment / pos mapping can be tested on a proteome whose FASTA order
    differs from the canonical order (tests/test_tier_b_toy.py)."""
    ids = fasta.canonical_order(seqs_by_id)
    seqs = [seqs_by_id[k] for k in ids]
    n = len(ids)
    i, j, pairs = pair_positions(ids, ta["gene_a"], ta["gene_b"])
    timings: dict = {}
    mm = search.MMseqs(binary=search.resolve_binary(ccfg), threads=threads, workdir=wd)
    fasta.write_fasta(wd / "input.fasta", ((str(k), s) for k, s in enumerate(seqs)))
    search.createdb(mm, wd / "input.fasta", "db")
    t = time.time()
    q = np.concatenate([pairs["i"], pairs["j"], np.arange(n)])
    tt = np.concatenate([pairs["j"], pairs["i"], np.arange(n)])
    aln = search.realign(mm, ccfg, "db", q, tt)
    timings["realign_s"] = round(time.time() - t, 1)
    feats, meta = tier_b.run_stage(mm, ccfg, pcfg, "db", seqs, pairs, aln, threads=threads,
                                   panel_dmin=panel_dmin, species_jobs=species_jobs,
                                   cache_root=cache_root)
    tool = search.tool_info(search.resolve_binary(ccfg))
    # evaluation sidecar (not a feature file): query-side pair bitscores, so label-free
    # validation can re-derive per-species classes (scripts/eval/tierB_validation.py)
    lo, hi = pairs["i"].to_numpy(), pairs["j"].to_numpy()
    b_lh, b_hl, _, _, _ = tier_b.directional(aln, lo, hi, n)
    a_is_lo = i == lo
    pairbits = pd.DataFrame({"gene_a": ta["gene_a"].to_numpy(), "gene_b": ta["gene_b"].to_numpy(),
                             "bits_a_to_b": np.where(a_is_lo, b_lh, b_hl),
                             "bits_b_to_a": np.where(a_is_lo, b_hl, b_lh)})
    out = pd.concat([ta[["gene_a", "gene_b"]].reset_index(drop=True), feats.reset_index(drop=True)], axis=1)
    return out, pairbits, meta, tool, timings


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("species")
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--panel-dmin", type=float, default=None)
    ap.add_argument("--keep-work", action="store_true")
    ap.add_argument("--species-jobs", type=int, default=1,
                    help="panel species searched concurrently (each with --threads); no effect on results")
    ap.add_argument("--search-only", action="store_true",
                    help="only fill the per-species panel-hit cache (no Tier A needed)")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
    t0 = time.time()
    prot = ROOT / "data/interim/proteomes" / f"{a.species}.fasta"
    fdir = ROOT / "data/interim/features" / a.species
    ta_path = fdir / "tierA.parquet"
    ccfg = load_yaml(CONFIG_DIR / "caller.yaml")
    pcfg = panel.load_panel_config()

    seqs_by_id, _ = fasta.read_proteome(prot)
    ids = fasta.canonical_order(seqs_by_id)
    seqs = [seqs_by_id[k] for k in ids]
    n = len(ids)
    if a.search_only:
        wd = fdir / "work_tierB_search"
        shutil.rmtree(wd, ignore_errors=True)
        wd.mkdir(parents=True)
        mm = search.MMseqs(binary=search.resolve_binary(ccfg), threads=a.threads, workdir=wd)
        fasta.write_fasta(wd / "input.fasta", ((str(k), s) for k, s in enumerate(seqs)))
        search.createdb(mm, wd / "input.fasta", "db")
        _, smeta = panel.search_panel(mm, ccfg, pcfg, "db", seqs, species_jobs=a.species_jobs)
        shutil.rmtree(wd, ignore_errors=True)
        print(json.dumps({"species": a.species, "cache": smeta["cache_dir"],
                          "search_wall_s": smeta["search_wall_s"]}))
        return 0
    ta_meta = json.loads(Path(str(ta_path) + ".meta.json").read_text())
    if ta_meta["input"]["sha256"] != sha256(prot):
        raise SystemExit("Tier A was computed on a different proteome file")
    ta = pd.read_parquet(ta_path, columns=["gene_a", "gene_b"])
    tag = "" if a.panel_dmin is None else f"_dmin{a.panel_dmin:g}"
    wd = fdir / f"work_tierB{tag}"
    shutil.rmtree(wd, ignore_errors=True)
    wd.mkdir(parents=True)
    try:
        out, pairbits, meta, tool, timings = compute_tierB(
            seqs_by_id, ta, ccfg, pcfg, wd, threads=a.threads, panel_dmin=a.panel_dmin,
            species_jobs=a.species_jobs)
    finally:
        if not a.keep_work:
            shutil.rmtree(wd, ignore_errors=True)
    timings.update(meta.pop("timings"))
    timings["total_s"] = round(time.time() - t0, 1)
    full = {
        "slparalog_version": __version__,
        "species": a.species,
        "schema_hash": pipeline.schema_hash(ccfg, meta["settings_hash"]),
        "registry_hash": registry.registry_hash(),
        "tier_b_columns": registry.TIER_B_COLUMNS,
        "aligned_to": {"path": str(ta_path.relative_to(ROOT)), "sha256": sha256(ta_path),
                       "n_rows": int(len(ta)), "tierA_schema_hash": ta_meta.get("schema_hash")},
        "input": {"path": str(prot.relative_to(ROOT)), "sha256": ta_meta["input"]["sha256"]},
        "panel_config": pcfg["tier_b"],
        "tool": tool,
        **meta,   # includes threads, species_jobs and the QC flags/warnings (output metadata only)
        "timings": timings,
    }
    dst = fdir / f"tierB{tag}.parquet"
    pipeline.write(out, dst, full)
    pairbits.to_parquet(fdir / f"tierB{tag}.pairbits.parquet", index=False)
    cal = meta["calibration"]
    print(json.dumps({"species": a.species, "n_pairs": len(out), "core_size": cal["core_size"],
                      "self_excluded": cal["self_excluded"], "n_usable": cal["n_usable"],
                      "qc_warnings": meta["qc"]["warnings"],
                      "timings": timings}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
