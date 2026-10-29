"""End-to-end runtime and memory benchmark of the shipped CLI (generalises tierA_bench.py).

Measures the whole process tree (python + mmseqs children; psutil sampling, summed RSS = upper bound as shared
pages are counted per process) and CPU time (RUSAGE_CHILDREN) of
  setup    `slparalog setup --panel --data-dir D -t T`   (one-off; the panel source files are first copied from the
           already-downloaded local copies into D/panel/downloads - `setup` has no local-source option but skips
           any download whose sha256 already matches - so network time is NOT included; the copy is untimed)
  predict  `slparalog predict <proteome> -o <out> --model M --data-dir D -t T --workdir W`  from an empty workdir
           and an empty panel hit cache (Tier A self-search + Tier B panel search + scoring)

Usage:
  .venv/bin/python scripts/bench/e2e_bench.py setup
  .venv/bin/python scripts/bench/e2e_bench.py predict yeast
  .venv/bin/python scripts/bench/e2e_bench.py predict human
Writes results/bench/e2e_setup.json / results/bench/e2e_<species>.json. The workdir is deleted afterwards
(its peak size is recorded). Idle-machine check: refuses to start if another process uses > --max-other-cpu % CPU
unless --force; the top processes at start are recorded either way.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import resource
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

import psutil

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))


def sh(*cmd: str) -> str:
    try:
        return subprocess.run(list(cmd), capture_output=True, text=True, cwd=ROOT).stdout.strip()
    except OSError:
        return ""


def machine() -> dict:
    return {
        "platform": platform.platform(), "cpu_model": sh("sysctl", "-n", "machdep.cpu.brand_string"),
        "logical_cores": os.cpu_count(), "physical_cores": psutil.cpu_count(logical=False),
        "performance_cores": sh("sysctl", "-n", "hw.perflevel0.logicalcpu"),
        "efficiency_cores": sh("sysctl", "-n", "hw.perflevel1.logicalcpu"),
        "mem_total_gb": round(psutil.virtual_memory().total / 1e9, 1),
        "pmset_g": sh("pmset", "-g"),
        "lowpowermode": next((ln.split()[-1] for ln in sh("pmset", "-g").splitlines() if "lowpowermode" in ln), None),
    }


def top_procs(n: int = 6) -> list[dict]:
    rows = []
    for p in psutil.process_iter(["pid", "name"]):
        try:
            rows.append((p.cpu_percent(interval=None), p.info["pid"], p.info["name"]))
        except psutil.Error:
            pass
    time.sleep(1.0)
    rows = []
    for p in psutil.process_iter(["pid", "name"]):
        try:
            rows.append({"pid": p.info["pid"], "name": p.info["name"], "cpu_pct": p.cpu_percent(interval=None)})
        except psutil.Error:
            pass
    return sorted(rows, key=lambda r: -r["cpu_pct"])[:n]


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


def dir_size(path: Path) -> int:
    tot = 0
    for dp, _, fns in os.walk(path):
        for f in fns:
            try:
                tot += os.lstat(os.path.join(dp, f)).st_size
            except OSError:
                pass
    return tot


class DirPeak(threading.Thread):
    """Samples the size of ``paths`` every ``every`` seconds (peak of the sum)."""

    def __init__(self, paths: list[Path], every: float = 5.0):
        super().__init__(daemon=True)
        self.paths, self.every, self.peak, self._stop_ev = paths, every, 0, threading.Event()

    def run(self):
        while not self._stop_ev.is_set():
            self.peak = max(self.peak, sum(dir_size(p) for p in self.paths if p.exists()))
            self._stop_ev.wait(self.every)

    def stop(self):
        self._stop_ev.set()
        self.join()
        self.peak = max(self.peak, sum(dir_size(p) for p in self.paths if p.exists()))


def measured(cmd: list[str], log_path: Path, interval: float, size_paths: list[Path]) -> dict:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log = open(log_path, "w")
    load0, ru0 = os.getloadavg(), resource.getrusage(resource.RUSAGE_CHILDREN)
    dp = DirPeak(size_paths)
    dp.start()
    t0 = time.time()
    proc = subprocess.Popen(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    ps = psutil.Process(proc.pid)
    peak_total = peak_single = max_procs = n = 0
    while proc.poll() is None:
        total, single, nproc = tree_rss(ps)
        peak_total, peak_single, max_procs = max(peak_total, total), max(peak_single, single), max(max_procs, nproc)
        n += 1
        time.sleep(interval)
    wall = time.time() - t0
    log.close()
    dp.stop()
    ru1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "command": " ".join(cmd), "returncode": proc.returncode, "wall_s": round(wall, 2),
        "cpu_user_s": round(ru1.ru_utime - ru0.ru_utime, 1), "cpu_sys_s": round(ru1.ru_stime - ru0.ru_stime, 1),
        "peak_rss_tree_gb": round(peak_total / 1e9, 3), "peak_rss_single_process_gb": round(peak_single / 1e9, 3),
        "max_processes_in_tree": max_procs, "rss_samples": n, "sample_interval_s": interval,
        "peak_work_size_gb": round(dp.peak / 1e9, 3), "loadavg_start": load0, "loadavg_end": os.getloadavg(),
        "log": str(log_path.relative_to(ROOT)),
    }


def stage_downloads(data_dir: Path) -> dict:
    """Copy the already-downloaded panel source files to <data>/panel/downloads (setup then finds them verified)."""
    from slparalog import setup_data as S

    src, _ = S.panel_sources()
    dl = data_dir / "panel" / "downloads"
    dl.mkdir(parents=True, exist_ok=True)
    n = 0
    for s in src:
        dest = dl / s["filename"]
        if dest.exists():
            continue
        found = subprocess.run(["find", str(ROOT / "data/raw"), "-name", s["filename"]], capture_output=True,
                               text=True).stdout.split()
        good = [f for f in found if S._sha256_file(Path(f)) == s["sha256"]]  # same name may exist in several dirs
        if not good:
            raise SystemExit(f"no local copy of {s['filename']} under data/raw with the pinned sha256")
        shutil.copy2(good[0], dest)
        n += 1
    return {"n_panel_files": len(src), "copied": n}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["setup", "predict"])
    ap.add_argument("species", nargs="?")
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--interval", type=float, default=0.5)
    ap.add_argument("--data-dir", default="data/release/slparalog_data")
    ap.add_argument("--model", default="results/frozen/model")
    ap.add_argument("--out-dir", default="results/release")
    ap.add_argument("--max-other-cpu", type=float, default=150.0, help="refuse to start if the busiest other process exceeds this %% CPU")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    dd = (ROOT / a.data_dir).resolve()
    py = sys.executable
    head = sh("git", "rev-parse", "HEAD")
    dirty = bool(sh("git", "status", "--porcelain", "--untracked-files=no"))
    top = top_procs()
    if top and top[0]["cpu_pct"] > a.max_other_cpu and not a.force:
        raise SystemExit(f"machine not idle: {top[:3]}")
    common = {"threads": a.threads, "machine": machine(), "git_head_at_run": head, "git_tracked_changes_at_run": dirty,
              "top_processes_at_start": top, "loadavg_at_start": os.getloadavg(),
              "disk_free_gb_at_start": round(shutil.disk_usage(ROOT).free / 1e9, 1),
              "python": sys.version.split()[0]}
    dst_dir = ROOT / "results/bench"
    dst_dir.mkdir(parents=True, exist_ok=True)

    if a.step == "setup":
        staged = stage_downloads(dd)
        cmd = [py, "-m", "slparalog.cli", "setup", "--panel", "--data-dir", str(dd), "-t", str(a.threads), "-v"]
        res = measured(cmd, dd.parent / "setup.log", a.interval, [dd])
        res.update({"step": "setup", "staged_local_files": staged,
                    "note": "one-off; source proteomes were copied from the local downloads beforehand (network time "
                            "not included); setup would otherwise download the 80 checksummed proteome files "
                            "(config/panel_species.tsv, config/panel_sources.json)",
                    "data_dir_size_gb_after": round(dir_size(dd) / 1e9, 3), **common})
        (dst_dir / "e2e_setup.json").write_text(json.dumps(res, indent=1))
        print(json.dumps({k: res[k] for k in ("returncode", "wall_s", "peak_rss_tree_gb", "data_dir_size_gb_after")}))
        return res["returncode"]

    sp = a.species
    fasta = ROOT / "data/interim/proteomes" / f"{sp}.fasta"
    out_dir = ROOT / a.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{sp}_predictions.tsv"
    work = dd.parent / f"bench_work_{sp}"
    hits = dd / "panel" / "hits"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)  # fresh and empty
    if hits.exists():
        shutil.rmtree(hits)  # no cached panel hits
    hits.mkdir()
    cmd = [py, "-m", "slparalog.cli", "predict", str(fasta), "-o", str(out), "--model", str(ROOT / a.model),
           "--data-dir", str(dd), "-t", str(a.threads), "--workdir", str(work), "-v"]
    res = measured(cmd, out_dir / f"{sp}_predict.log", a.interval, [work, hits])
    meta = {}
    mp = Path(str(out) + ".meta.json")
    if res["returncode"] == 0 and mp.exists():
        meta = json.loads(mp.read_text())
    res.update({"step": "predict", "species": sp, "model_dir": a.model, "data_dir": a.data_dir,
                "output": str(out.relative_to(ROOT)), "n_pairs_scored": (sum(1 for ln in open(out) if not ln.startswith("#")) - 1) if out.exists() else None,
                "n_proteins": meta.get("n_proteins"), "pipeline_timings": meta.get("timings"),
                "model_id": (meta.get("model") or {}).get("model_id"), "workdir_empty_at_start": True,
                "panel_hit_cache_empty_at_start": True, "setup_record": "results/bench/e2e_setup.json", **common})
    shutil.rmtree(work, ignore_errors=True)
    res["workdir_deleted"] = not work.exists()
    (dst_dir / f"e2e_{sp}.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ("species", "returncode", "wall_s", "peak_rss_tree_gb",
                                           "peak_rss_single_process_gb", "peak_work_size_gb")}))
    return res["returncode"]


if __name__ == "__main__":
    sys.exit(main())
