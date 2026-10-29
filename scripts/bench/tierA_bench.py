"""Run `slparalog predict --features-only` (caller + Tier A) on a per-gene proteome and record
wall time and peak memory of the whole process tree (python + mmseqs children, sampled with psutil;
`/usr/bin/time -l` alone misses the children), plus the Low Power Mode state.

Usage: .venv/bin/python scripts/bench/tierA_bench.py human [--threads 10]
Writes data/interim/features/<species>/tierA.parquet (+ .meta.json) and results/bench/tierA_<species>.json.
Peak RSS is the maximum over samples of the summed RSS of all processes in the tree (shared pages
are counted once per process, so this is an upper bound).
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

import psutil

ROOT = Path(__file__).resolve().parents[2]


def lowpowermode() -> str | None:
    try:
        out = subprocess.run(["pmset", "-g"], capture_output=True, text=True).stdout
    except OSError:
        return None
    for line in out.splitlines():
        if "lowpowermode" in line:
            return line.split()[-1]
    return None


def tree_rss(proc: psutil.Process) -> tuple[int, int, int]:
    procs = [proc]
    try:
        procs += proc.children(recursive=True)
    except psutil.Error:
        pass
    total, biggest = 0, 0
    for p in procs:
        try:
            r = p.memory_info().rss
        except psutil.Error:
            continue
        total += r
        biggest = max(biggest, r)
    return total, biggest, len(procs)


def git_head() -> str | None:
    r = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True)
    return r.stdout.strip() or None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("species")
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--interval", type=float, default=0.5)
    ap.add_argument("--keep-work", action="store_true")
    a = ap.parse_args()

    fasta = ROOT / "data/interim/proteomes" / f"{a.species}.fasta"
    outdir = ROOT / "data/interim/features" / a.species
    out = outdir / "tierA.parquet"
    work = outdir / "work"
    if work.exists():
        shutil.rmtree(work)
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, "-m", "slparalog.cli", "predict", str(fasta), "-o", str(out),
           "--features-only", "-t", str(a.threads), "--workdir", str(work), "-v"]
    log = open(outdir / "tierA.log", "w")
    lpm_start, load_start = lowpowermode(), os.getloadavg()
    t0 = time.time()
    proc = subprocess.Popen(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    ps = psutil.Process(proc.pid)
    peak_total = peak_single = max_procs = 0
    n_samples = 0
    while proc.poll() is None:
        total, single, nproc = tree_rss(ps)
        peak_total, peak_single = max(peak_total, total), max(peak_single, single)
        max_procs = max(max_procs, nproc)
        n_samples += 1
        time.sleep(a.interval)
    wall = time.time() - t0
    log.close()
    rc = proc.returncode
    meta = {}
    meta_path = Path(str(out) + ".meta.json")
    if rc == 0 and meta_path.exists():
        meta = json.loads(meta_path.read_text())
    if rc == 0 and not a.keep_work:
        shutil.rmtree(work, ignore_errors=True)
    res = {
        "species": a.species,
        "command": " ".join(cmd[1:]),
        "returncode": rc,
        "threads": a.threads,
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
        "output": str(out.relative_to(ROOT)),
        "n_proteins": meta.get("n_proteins"),
        "caller": meta.get("caller"),
        "saturation": meta.get("saturation"),
        "n_search_hits_nonself": meta.get("n_search_hits_nonself"),
        "features_summary": meta.get("features"),
        "pipeline_timings": meta.get("timings"),
        "input": {k: meta.get("input", {}).get(k) for k in
                  ("sha256", "n_records", "n_kept", "n_dropped_empty", "n_terminal_stops_stripped",
                   "n_internal_stops_to_x", "n_nonstandard_to_x")},
        "schema_hash": meta.get("schema_hash"),
        "tool": meta.get("tool"),
        "tool_matches_config": meta.get("tool_matches_config"),
    }
    dst = ROOT / "results/bench" / f"tierA_{a.species}.json"
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():  # keep the previous run's record (one level) for comparison
        prev = json.loads(dst.read_text())
        prev.pop("previous", None)
        res["previous"] = prev
    dst.write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ("species", "returncode", "wall_s", "peak_rss_tree_gb",
                                           "lowpowermode_start")}))
    return rc


if __name__ == "__main__":
    sys.exit(main())
