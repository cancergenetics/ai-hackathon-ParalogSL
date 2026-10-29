"""Append-only development experiment registry ``results/experiments.tsv`` (DESIGN v3 §5)."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pandas as pd

from slparalog.config import PROJECT_ROOT

EXPERIMENTS_TSV = PROJECT_ROOT / "results" / "experiments.tsv"
COLUMNS = [
    "exp_id", "run_id", "timestamp", "git_hash", "git_dirty_paths", "tier", "family", "name",
    "features", "n_features", "params", "adoption_rule",
    "ap_human", "ap_human_se", "auroc_human", "auroc_human_se", "prev_human", "npos_human", "nneg_human",
    "ap_yeast", "ap_yeast_se", "auroc_yeast", "auroc_yeast_se", "prev_yeast", "npos_yeast", "nneg_yeast",
    "T_vs_B1", "T_vs_B1_se", "selected", "notes",
    # D12 additions
    "T_perfold_vs_B1", "T_perfold_vs_B1_se", "T_perfold_lo", "T_perfold_hi", "human_strat_avg_ap",
    "incumbent", "incumbent_human_strat_avg_ap", "guard_pass", "adopted", "deconfound", "label_files_sha256",
    "feature_tables",
]


def git_state(paths: list[str] | None = None) -> tuple[str, str]:
    """(HEAD hash, comma-joined list of modified tracked/untracked files under ``paths``)."""
    try:
        h = subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, capture_output=True, text=True,
                           check=True).stdout.strip()
        cmd = ["git", "status", "--porcelain"] + (["--"] + paths if paths else [])
        st = subprocess.run(cmd, cwd=PROJECT_ROOT, capture_output=True, text=True, check=True).stdout
        dirty = ",".join(sorted(ln[3:] for ln in st.splitlines() if ln.strip()))
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown", "unknown"
    return h, dirty


def append(rows: list[dict], path: Path = EXPERIMENTS_TSV) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    for c in COLUMNS:
        if c not in df.columns:
            df[c] = ""
    for c in ("params",):
        df[c] = df[c].map(lambda v: v if isinstance(v, str) else json.dumps(v, sort_keys=True))
    df = df[COLUMNS]
    if path.exists():
        old = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if list(old.columns) != COLUMNS:  # schema grew: rewrite with the new header, old rows kept verbatim
            for c in COLUMNS:
                if c not in old.columns:
                    old[c] = ""
            unknown = [c for c in old.columns if c not in COLUMNS]
            if unknown:
                raise ValueError(f"experiments.tsv has unknown columns {unknown}")
            old[COLUMNS].to_csv(path, sep="\t", index=False)
    header = not path.exists()
    df.to_csv(path, sep="\t", index=False, header=header, mode="a", float_format="%.6g")
