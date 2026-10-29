#!/usr/bin/env python
"""Tier A+B features of the degraded proteomes (robustness test, DESIGN §6.5).

For every requested species/kind with a degraded proteome under data/interim/robustness/<species>/<kind>/, runs

    .venv/bin/slparalog predict <abs proteome> -o <abs features.parquet> --features-only --tier-b -t <threads> --workdir <abs work>

(existing panel databases; the per-query panel-hit cache is written under data/interim/panel/hits) and writes
runtime.json next to the outputs: wall time of the command, the tool's own timing split, peak RSS of the process
tree (``/usr/bin/time -l``) and the panel-hit cache directories created by the run. ``--cleanup`` removes the
work directories and those hit-cache directories (keeps proteome.fasta, mapping.tsv, features.parquet, runtime.json).

  .venv/bin/python scripts/robustness/run_degraded_features.py --species yeast --kinds all --threads 8
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RDIR = ROOT / "data" / "interim" / "robustness"
HITS = ROOT / "data" / "interim" / "panel" / "hits"
SLPARALOG = ROOT / ".venv" / "bin" / "slparalog"
MAX_THREADS = 8
KINDS = {"yeast": ["identity", "fragments", "split", "fused", "missing", "haplotig", "combined"], "human": ["combined"]}


def run_one(species: str, kind: str, threads: int, force: bool) -> dict:
    kd = RDIR / species / kind
    out = kd / "features.parquet"
    if out.exists() and not force:
        print("exists, skipped:", out)
        return json.loads((kd / "runtime.json").read_text()) if (kd / "runtime.json").exists() else {}
    work = kd / "work"
    before = {p.name for p in HITS.iterdir()} if HITS.exists() else set()
    cmd = ["/usr/bin/time", "-l", str(SLPARALOG), "predict", str(kd / "proteome.fasta"), "-o", str(out),
           "--features-only", "--tier-b", "-t", str(threads), "--workdir", str(work)]
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    wall = time.time() - t0
    (kd / "predict.stdout").write_text(proc.stdout)
    (kd / "predict.stderr").write_text(proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"{species}/{kind}: slparalog predict failed ({proc.returncode}); see {kd / 'predict.stderr'}")
    rss = re.search(r"(\d+)\s+maximum resident set size", proc.stderr)
    meta = json.loads(Path(str(out) + ".meta.json").read_text())
    rt = {"species": species, "kind": kind, "command": " ".join(cmd[2:]), "threads": threads, "wall_s": wall,
          "tool_timings": meta["timings"], "n_proteins": meta["n_proteins"], "n_pairs": meta["caller"]["n_pairs"],
          "peak_rss_bytes_process_tree_max_child": int(rss.group(1)) if rss else None,
          "new_panel_hit_dirs": sorted({p.name for p in HITS.iterdir()} - before) if HITS.exists() else []}
    (kd / "runtime.json").write_text(json.dumps(rt, indent=1))
    print(species, kind, f"wall={wall:.0f}s", proc.stdout.strip())
    return rt


def cleanup(species: str, kind: str) -> None:
    kd = RDIR / species / kind
    shutil.rmtree(kd / "work", ignore_errors=True)
    rt = json.loads((kd / "runtime.json").read_text()) if (kd / "runtime.json").exists() else {}
    for name in rt.get("new_panel_hit_dirs", []):
        shutil.rmtree(HITS / name, ignore_errors=True)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", required=True, choices=["human", "yeast"])
    p.add_argument("--kinds", default="all")
    p.add_argument("--threads", type=int, default=MAX_THREADS)
    p.add_argument("--force", action="store_true")
    p.add_argument("--cleanup", action="store_true", help="remove work dirs and the panel-hit caches created by the runs")
    a = p.parse_args(argv)
    if not 1 <= a.threads <= MAX_THREADS:
        raise SystemExit(f"--threads must be 1..{MAX_THREADS}")
    kinds = KINDS[a.species] if a.kinds == "all" else a.kinds.split(",")
    for k in kinds:
        if a.cleanup:
            cleanup(a.species, k)
        else:
            run_one(a.species, k, a.threads, a.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
