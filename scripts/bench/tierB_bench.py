"""Run scripts/features/run_tierB.py for a species and record wall time, peak RSS of the whole
process tree (python + mmseqs children, psutil sampling as in tierA_bench.py) and Low Power Mode.

Usage: .venv/bin/python scripts/bench/tierB_bench.py <species> [--threads 10] [--panel-dmin c]
Writes results/bench/tierB_<species>[_dmin<c>].json. A run whose panel hits were already cached
reports n_species_cached > 0: its search time is then not a cold-run measurement.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import psutil

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).parent))
from tierA_bench import git_head, lowpowermode, tree_rss  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("species")
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--interval", type=float, default=0.5)
    ap.add_argument("--panel-dmin", type=float, default=None)
    ap.add_argument("--species-jobs", type=int, default=1)
    ap.add_argument("--search-only", action="store_true",
                    help="measure the panel-search phase alone (fills the hit cache)")
    a = ap.parse_args()
    tag = ("" if a.panel_dmin is None else f"_dmin{a.panel_dmin:g}") + ("_search" if a.search_only else "")
    fdir = ROOT / "data/interim/features" / a.species
    cmd = [sys.executable, str(ROOT / "scripts/features/run_tierB.py"), a.species,
           "--threads", str(a.threads), "--species-jobs", str(a.species_jobs)]
    if a.panel_dmin is not None:
        cmd += ["--panel-dmin", str(a.panel_dmin)]
    if a.search_only:
        cmd += ["--search-only"]
    log = open(fdir / f"tierB{tag}.log", "w")
    lpm_start, load_start = lowpowermode(), os.getloadavg()
    t0 = time.time()
    proc = subprocess.Popen(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    ps = psutil.Process(proc.pid)
    peak_total = peak_single = max_procs = n_samples = 0
    while proc.poll() is None:
        total, single, nproc = tree_rss(ps)
        peak_total, peak_single = max(peak_total, total), max(peak_single, single)
        max_procs = max(max_procs, nproc)
        n_samples += 1
        time.sleep(a.interval)
    wall = time.time() - t0
    log.close()
    meta = {}
    mp = fdir / f"tierB{tag}.parquet.meta.json"
    if proc.returncode == 0 and mp.exists():
        meta = json.loads(mp.read_text())
    s = meta.get("search", {})
    per = s.get("per_species", {})
    res = {
        "species": a.species,
        "command": " ".join(cmd[1:]),
        "returncode": proc.returncode,
        "threads": a.threads,
        "species_jobs": a.species_jobs,
        "panel_dmin": a.panel_dmin,
        "wall_s": round(wall, 2),
        "peak_rss_tree_gb": round(peak_total / 1e9, 3),
        "peak_rss_single_process_gb": round(peak_single / 1e9, 3),
        "max_processes_in_tree": max_procs,
        "rss_samples": n_samples,
        "sample_interval_s": a.interval,
        "lowpowermode_start": lpm_start,
        "lowpowermode_end": lowpowermode(),
        "loadavg_start": load_start,
        "loadavg_end": os.getloadavg(),
        "machine": {"platform": platform.platform(), "cpu_count": os.cpu_count(),
                    "mem_total_gb": round(psutil.virtual_memory().total / 1e9, 1)},
        "git_head_at_run": git_head(),
        "timings": meta.get("timings"),
        "note_realign": "realign_s re-aligns Tier A pairs (needed only because Tier B runs after a separate Tier A run; the integrated pipeline reuses its realignment)",
        "panel_search_wall_s": s.get("search_wall_s"),
        "n_species_cached": s.get("n_species_cached"),
        "panel_search_s_sum_per_species_logs": round(sum(v.get("search_s", 0) for v in per.values()), 1) if per else None,
        "n_capped_first_round_total": s.get("n_capped_first_round_total"),
        "n_capped_final_total": s.get("n_capped_final_total"),
        "n_pairs": meta.get("aligned_to", {}).get("n_rows"),
        "calibration_summary": {k: meta.get("calibration", {}).get(k) for k in
                                ("self_excluded", "n_retained", "n_usable", "core_size", "core_frac_used")},
        "features_summary": meta.get("features"),
        "schema_hash": meta.get("schema_hash"),
        "output": str((fdir / f"tierB{tag}.parquet").relative_to(ROOT)) if mp.exists() else None,
    }
    dst = ROOT / "results/bench" / f"tierB_{a.species}{tag}.json"
    dst.write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ("species", "returncode", "wall_s", "peak_rss_tree_gb",
                                           "lowpowermode_start")}))
    return proc.returncode


if __name__ == "__main__":
    sys.exit(main())
