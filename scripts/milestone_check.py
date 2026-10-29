#!/usr/bin/env python
"""Milestone gate (G7). Usage: .venv/bin/python scripts/milestone_check.py --milestone M1

Checks (git-ancestry ordering checks are a documented manual check, not code):
  1. pytest -m "not slow" has 0 failures; collected count >= previous milestone's.
  2. No skip/xfail outside tests/skip_allowlist.txt (strict xfail allowed only for tests listed in
     results/labels/known_answer_gaps.json); pytest config has no deselecting options.
  3. `git status --porcelain` is empty.
  4. `du -sk` of the whole project < 150 GB.
  5. requirements-lock.txt == `pip freeze --exclude-editable`.
  6. config/tools.yaml exists, records mmseqs and foldseek versions, and sha256 of bin/ binaries match.
Writes results/milestones/<M>.json (always; `passed` says whether the gate holds). Exit 0 iff all pass.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
MAX_KB = 150 * 1024 * 1024  # 150 GB in KiB
SKIP_ALLOWLIST = ROOT / "tests" / "skip_allowlist.txt"
GAPS_JSON = ROOT / "results" / "labels" / "known_answer_gaps.json"
LOCK = ROOT / "requirements-lock.txt"
TOOLS = ROOT / "config" / "tools.yaml"
BAD_ADDOPTS = re.compile(r"(^|\s)(-k\b|--deselect|--ignore|--ignore-glob|-m\b|--co\b|--collect-only|-p\s+no:)")


def run(cmd, **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, **kw)


def check_venv() -> None:
    if Path(sys.prefix).resolve() != (ROOT / ".venv").resolve():
        sys.exit(f"must run with {ROOT}/.venv/bin/python (running {sys.executable})")


def read_allowlist(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip() and not ln.strip().startswith("#")]


def flatten_strings(obj) -> set[str]:
    out: set[str] = set()
    if isinstance(obj, str):
        out.add(obj)
    elif isinstance(obj, dict):
        for k, v in obj.items():
            out.add(str(k))
            out |= flatten_strings(v)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            out |= flatten_strings(v)
    return out


def nodeid(case: ET.Element) -> str:
    cls = case.get("classname", "")
    path = cls.replace(".", "/")
    # classname is module path (maybe with a class suffix); the file is the longest prefix that exists
    parts = cls.split(".")
    for i in range(len(parts), 0, -1):
        cand = "/".join(parts[:i]) + ".py"
        if (ROOT / cand).exists():
            path = cand
            rest = parts[i:]
            return "::".join([path, *rest, case.get("name", "")])
    return f"{path}::{case.get('name', '')}"


def check_pytest(milestone: str) -> dict:
    res: dict = {"name": "pytest"}
    problems: list[str] = []
    with tempfile.TemporaryDirectory() as td:
        xml = Path(td) / "junit.xml"
        proc = run([sys.executable, "-m", "pytest", "-m", "not slow", "-q", f"--junitxml={xml}", "-p", "no:cacheprovider"])
        res["returncode"] = proc.returncode
        res["tail"] = "\n".join(proc.stdout.strip().splitlines()[-15:])
        if not xml.exists():
            return {**res, "passed": False, "problems": ["pytest produced no junit xml"]}
        suite = ET.parse(xml).getroot()
        if suite.tag == "testsuites":
            suite = suite[0]
        total = int(suite.get("tests", 0))
        failures = int(suite.get("failures", 0)) + int(suite.get("errors", 0))
        allow = read_allowlist(SKIP_ALLOWLIST)
        gaps = flatten_strings(json.loads(GAPS_JSON.read_text())) if GAPS_JSON.exists() else set()
        skipped = []
        for case in suite.iter("testcase"):
            sk = case.find("skipped")
            if sk is None:
                continue
            nid = nodeid(case)
            kind = "xfail" if "xfail" in (sk.get("type", "") + sk.get("message", "")).lower() else "skip"
            ok = any(nid == a or nid.startswith(a + "::") or nid.split("::")[0] == a for a in allow)
            if kind == "xfail" and not ok:
                ok = any(nid in g or g in nid or case.get("name", "") == g for g in gaps if g)
            skipped.append({"test": nid, "kind": kind, "allowed": ok})
            if not ok:
                problems.append(f"unlisted {kind}: {nid}")
        res.update(collected=total, failures=failures, skipped=skipped)
        if proc.returncode != 0 or failures:
            problems.append(f"pytest failed: returncode={proc.returncode}, failures+errors={failures}")
        # previous milestone count
        prev = []
        for f in sorted((ROOT / "results" / "milestones").glob("M*.json")):
            if f.stem != milestone:
                try:
                    prev.append(json.loads(f.read_text())["checks"]["pytest"]["collected"])
                except Exception:
                    pass
        if prev and total < max(prev):
            problems.append(f"collected {total} < previous milestone count {max(prev)}")
    # pytest config must not deselect tests or hide markers
    cfg = ROOT / "pyproject.toml"
    if cfg.exists():
        import tomllib

        ini = tomllib.loads(cfg.read_text()).get("tool", {}).get("pytest", {}).get("ini_options", {})
        addopts = ini.get("addopts", "")
        addopts = " ".join(addopts) if isinstance(addopts, list) else addopts
        if BAD_ADDOPTS.search(addopts):
            problems.append(f"pytest addopts deselects/ignores tests: {addopts!r}")
        if ini.get("xfail_strict") is not True:
            problems.append("pytest xfail_strict is not true")
    for extra in ("pytest.ini", "tox.ini", "setup.cfg"):
        if (ROOT / extra).exists():
            problems.append(f"unexpected pytest config file {extra}")
    return {**res, "passed": not problems, "problems": problems}


def check_git_clean() -> dict:
    out = run(["git", "status", "--porcelain"]).stdout
    head = run(["git", "rev-parse", "HEAD"]).stdout.strip()
    lines = out.splitlines()
    return {"name": "git_clean", "passed": not lines, "head": head, "dirty": lines[:50]}


def check_disk() -> dict:
    out = run(["du", "-sk", str(ROOT)]).stdout.split()
    kb = int(out[0])
    return {"name": "disk", "passed": kb < MAX_KB, "du_sk": kb, "limit_kb": MAX_KB, "gib": round(kb / 1024 / 1024, 2)}


def check_lock() -> dict:
    freeze = run([sys.executable, "-m", "pip", "freeze", "--exclude-editable"]).stdout
    if not LOCK.exists():
        return {"name": "requirements_lock", "passed": False, "problems": ["requirements-lock.txt missing"]}
    a = sorted(x for x in LOCK.read_text().splitlines() if x.strip())
    b = sorted(x for x in freeze.splitlines() if x.strip())
    diff = sorted(set(a) ^ set(b))
    return {"name": "requirements_lock", "passed": not diff, "diff": diff}


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_tools() -> dict:
    problems = []
    if not TOOLS.exists():
        return {"name": "tools", "passed": False, "problems": ["config/tools.yaml missing"]}
    cfg = yaml.safe_load(TOOLS.read_text()) or {}
    for tool in ("mmseqs", "foldseek"):
        entry = cfg.get(tool)
        if not entry or not entry.get("version") or not entry.get("sha256"):
            problems.append(f"{tool}: version/sha256 missing in tools.yaml")
    for tool, entry in cfg.items():
        binary = ROOT / (entry or {}).get("path", f"bin/{tool}")
        if not binary.exists():
            problems.append(f"{tool}: {binary} missing")
            continue
        if sha256_file(binary) != entry.get("sha256"):
            problems.append(f"{tool}: sha256 of {binary.relative_to(ROOT)} differs from tools.yaml")
        ver = run([str(binary), "version"]).stdout.strip().splitlines()
        if not ver or ver[0].strip() != str(entry.get("version")):
            problems.append(f"{tool}: reported version {ver[:1]} != tools.yaml {entry.get('version')}")
    return {"name": "tools", "passed": not problems, "problems": problems}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--milestone", required=True, choices=["M1"])
    args = ap.parse_args()
    check_venv()
    # git/disk are evaluated before this script writes its own output file
    checks = {}
    for fn in (check_git_clean, check_disk):
        r = fn()
        checks[r["name"]] = r
    for fn in (lambda: check_pytest(args.milestone), check_lock, check_tools):
        r = fn()
        checks[r["name"]] = r
    passed = all(c["passed"] for c in checks.values())
    out = {
        "milestone": args.milestone,
        "passed": passed,
        "timestamp_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "git_head": checks["git_clean"]["head"],
        "checks": checks,
        "manual_checks_not_in_code": [
            "git ancestry: config/labels.yaml first commit precedes any scripts/labels/build_*.py",
            "git ancestry: results/mde.json first commit precedes any src/slparalog/features/ file",
        ],
    }
    dest = ROOT / "results" / "milestones" / f"{args.milestone}.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2) + "\n")
    for name, c in checks.items():
        print(f"{'PASS' if c['passed'] else 'FAIL'}  {name}  {'; '.join(c.get('problems', []) + c.get('diff', []) + c.get('dirty', [])[:5])}")
    print(f"{args.milestone}: {'PASSED' if passed else 'FAILED'} -> {dest.relative_to(ROOT)}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
