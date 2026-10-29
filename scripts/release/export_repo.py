#!/usr/bin/env python
"""Export a clean, history-free, public-ready copy of the project.

Run (from the source repo, with the project venv):
    .venv/bin/python scripts/release/export_repo.py [--commit]

Copies an explicit ALLOWLIST of tracked files (``git ls-files``) into ``export/ai-hackathon-ParalogSL/`` (deleted and
recreated on every run), sanitises text files (absolute paths, user name, host name, e-mail addresses), writes
``EXPORT_MANIFEST.tsv`` (path, sha256, sanitised yes/no) and prints a summary with every exclusion.
With ``--commit`` it also runs ``git init -b main`` and makes one commit in the export (no remote, no push; the git
identity is whatever the user's global git config says).

Included: root docs, pyproject/requirements-lock, config/, src/, scripts/, tests/ (incl. verified-synthetic fixtures),
reviews/, docs/ (AGENT_RULES.md, lit/, briefing/ if present), small results (*.json, *.md, *.tsv <= 1.5 MB, minus any
file with more than 20 gene identifiers), the two trained models. Everything else (all of data/, models_cache/, bin/,
logs/, FASTA/parquet/npy/zip/gz/xlsx, predictions) is not copied. Labels are never opened.

Sanitising replaces values only (never removes JSON keys), so nothing a report cites as a key disappears; JSON files
are re-parsed after sanitising. The identifiers in the sanitising patterns are built from fragments so that this
file survives its own sanitising when it is copied into the export.
"""

from __future__ import annotations

import argparse
import getpass
import hashlib
import json
import re
import shutil
import socket
import subprocess
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[2]
NAME = "ai-hackathon-ParalogSL"
OUT = SRC / "export" / NAME
TEMPLATES = SRC / "scripts" / "release" / "templates"

ROOT_DOCS = ["PROMPT.md", "LITERATURE.md", "DESIGN.md", "RESEARCH_TREE.md", "LESSONS.md", "DATA.md", "REPORT.md"]
ROOT_FILES = ["pyproject.toml", "requirements-lock.txt"]
DIR_PREFIXES = ["config/", "src/", "scripts/", "tests/", "reviews/", "docs/lit/"]
DOC_FILES = ["docs/AGENT_RULES.md"]
BRIEFING = "docs/briefing"  # copied from the file system (it may not be tracked yet)
MODEL_DIRS = ["results/frozen/model", "results/release/model"]
SKIP_PREFIXES = ["scripts/release/templates/"]  # templates are copied to their export locations by this script
RESULT_EXTS = (".json", ".md", ".tsv")
RESULT_MAX_BYTES = int(1.5 * 1024 * 1024)
ID_LIMIT = 20  # more than this many gene identifiers in a result file => per-pair/per-gene rows => excluded
BLOCK_EXTS = (".fasta", ".fa", ".faa", ".parquet", ".npy", ".npz", ".zip", ".gz", ".xlsx", ".joblib", ".log")
EXPLICIT_RESULT_EXCLUDES = ("results/review/tierB_gain/",)  # *.tsv.gz / *.npz there (plus anything else not scanned)
BLOCK_PREFIXES = ["data/", "models_cache/", "bin/", "logs/", ".venv/", "export/"]
# synthetic test fixtures allowed past the extension block: file -> generator script that makes it (all toy)
FIXTURE_GENERATORS = {
    "tests/fixtures/toy_proteome.fasta": "tests/fixtures/make_toy_proteome.py",
    "tests/fixtures/expected_tierB_toy.parquet": "tests/fixtures/make_expected_tierB_toy.py",
    "tests/fixtures/final_test_toy/": "tests/fixtures/make_toy_final_test.py",
}

ID_RE = re.compile(r"ENSG\d+|\bY[A-P][LR]\d{3}[WC]\b|AT\dG\d+|WBGene\d+")
_U = getpass.getuser()  # local identity is read at run time, never written into the code
_H = str(Path.home())
_S = str(SRC)
_HOST = socket.gethostname()
EMAIL_RE = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)*\.[A-Za-z]{2,}")
KEEP_EMAIL = "noreply@anthropic.com"
SANITISE = [  # (compiled regex, replacement), applied in order
    (re.compile(re.escape(_S) + "/"), ""),
    (re.compile(re.escape(_S) + r"(?![A-Za-z0-9_])"), "."),
    (re.compile(r"/(?:private/)?" + "tmp/claude-501/" + r"[^\s\"'`)\]}>,;:]*"), "<scratchpad>"),
    (re.compile(re.escape(_H)), "<home>"),
    (re.compile(re.escape(_U)), "<user>"),
    (re.compile(re.escape(_HOST)), "<host>"),
    (re.compile(re.escape(_HOST.removesuffix(".local"))), "<host>"),
]
TEXT_EXTS = {".py", ".md", ".txt", ".json", ".tsv", ".yaml", ".yml", ".toml", ".sh", ".sha256", ".cfg", ".csv", ".html",
             ".svg", ".css", ".js", ".ini", ".tex", ".bib", ".log", ""}


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sanitise_text(s: str) -> str:
    for rx, rep in SANITISE:
        s = rx.sub(rep, s)
    return EMAIL_RE.sub(lambda m: m.group(0) if m.group(0) == KEEP_EMAIL else "<email>", s)


# Review tables that quote per-pair label calls and interaction scores from a source without a
# redistribution licence (Costanzo 2016 / TheCellMap; DATA.md). Rows starting with a yeast ORF are
# replaced by one note; the review's verdicts and text are kept (PI decision D23).
REDACT_PAIR_ROWS = {"reviews/results/yeast_labels_review.md": re.compile(r"^\| *Y[A-P][LR]\d{3}[WC].*$", re.M)}
REDACT_NOTE = ("| *rows redacted in the public export: per-pair calls and interaction scores from a source "
               "without a redistribution licence (see DATA.md)* | | | |")


# JSON keys whose values are full local system dumps (e.g. `pmset -g`, which lists running apps); the
# power state used in the analysis is kept in its own field (lowpowermode_*).
REDACT_JSON_KEYS = {"pmset_g"}
REDACT_JSON_NOTE = "<redacted in public export: full local power-settings dump; Low Power Mode is in its own field>"


def redact_json_keys(s: str) -> str:
    d = json.loads(s)
    hit = False

    def walk(x):
        nonlocal hit
        if isinstance(x, dict):
            for k in list(x):
                if k in REDACT_JSON_KEYS:
                    x[k], hit = REDACT_JSON_NOTE, True
                else:
                    walk(x[k])
        elif isinstance(x, list):
            for v in x:
                walk(v)

    walk(d)
    return json.dumps(d, indent=1) + "\n" if hit else s


def redact_pair_rows(rel: str, s: str) -> str:
    rx = REDACT_PAIR_ROWS.get(rel)
    if rx is None:
        return s
    out, noted = [], False
    for line in s.splitlines(keepends=True):
        if rx.match(line.rstrip("\n")):
            if not noted:
                out.append(REDACT_NOTE + "\n")
                noted = True
            continue
        out.append(line)
    return "".join(out)


def n_ids(data: bytes) -> int:
    return len(ID_RE.findall(data.decode("utf-8", errors="ignore")))


def git_files() -> list[str]:
    out = subprocess.run(["git", "ls-files", "-z"], cwd=SRC, capture_output=True, check=True).stdout
    return [f for f in out.decode().split("\0") if f and (SRC / f).is_file()]


def fixture_generator(rel: str) -> str | None:
    for k, v in FIXTURE_GENERATORS.items():
        if rel == k or (k.endswith("/") and rel.startswith(k)):
            return v
    return None


def check_fixture(rel: str, excluded: list[tuple[str, str]]) -> bool:
    """A binary-extension test fixture is copied only if its generator exists and its rows carry no real identifier."""
    gen = fixture_generator(rel)
    if gen is None or not (SRC / gen).is_file():
        excluded.append((rel, "fixture without a toy generator script: not verified synthetic"))
        return False
    p = SRC / rel
    try:
        if p.suffix == ".parquet":
            import pandas as pd

            df = pd.read_parquet(p)
            text = df.to_csv(sep="\t") + "\n".join(map(str, df.columns))
        else:
            text = p.read_text(errors="ignore")
    except Exception as e:  # pragma: no cover
        excluded.append((rel, f"fixture could not be inspected ({type(e).__name__}): not verified synthetic"))
        return False
    k = len(ID_RE.findall(text))
    if k:
        excluded.append((rel, f"fixture contains {k} real-looking gene identifiers"))
        return False
    return True


def check_model_json(p: Path) -> list[str]:
    """model.json must hold metadata only: no gene identifiers, no long data-like lists."""
    problems = []
    raw = p.read_bytes()
    if n_ids(raw):
        problems.append(f"{n_ids(raw)} gene identifiers")
    d = json.loads(raw)

    def walk(x, path=""):
        if isinstance(x, dict):
            for k, v in x.items():
                walk(v, f"{path}/{k}")
        elif isinstance(x, list):
            if len(x) > 200 or (len(x) > 20 and any(isinstance(e, (dict, list)) for e in x)):
                problems.append(f"long or nested list at {path} ({len(x)} items)")
            for i, v in enumerate(x[:5]):
                walk(v, f"{path}[{i}]")

    walk(d)
    return problems


def check_joblib(p: Path) -> list[str]:
    """model.joblib: trained trees only. Byte scan for identifiers/paths, then a load to inspect the top-level types."""
    problems = []
    raw = p.read_bytes()
    if n_ids(raw):
        problems.append(f"{n_ids(raw)} gene identifiers in bytes")
    if p.stat().st_size > 5 * 1024 * 1024:
        problems.append("larger than 5 MB: not a small tree model")
    try:
        import joblib

        obj = joblib.load(p)
        desc = {}

        def walk(o, name, depth=0):
            tn = type(o).__module__ + "." + type(o).__name__
            if tn.startswith("pandas.") and "DataFrame" in tn:
                problems.append(f"DataFrame at {name}")
            if depth < 3 and isinstance(o, dict):
                for k, v in o.items():
                    walk(v, f"{name}/{k}", depth + 1)
            elif depth < 3 and hasattr(o, "__dict__"):
                for k, v in vars(o).items():
                    walk(v, f"{name}.{k}", depth + 1)
            desc[name] = tn

        walk(obj, "obj")
    except Exception as e:
        problems.append(f"could not load for inspection ({type(e).__name__}: {e})")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--commit", action="store_true", help="git init -b main + one commit in the export (no remote)")
    args = ap.parse_args()

    tracked = git_files()
    selected: dict[str, Path] = {}  # export path -> source path
    excluded: list[tuple[str, str]] = []  # (path, reason) for everything inside the allowlist that was dropped
    notes: list[str] = []
    n_tracked_unlisted = 0

    for rel in tracked:
        low = rel.lower()
        in_models = any(rel.startswith(m + "/") for m in MODEL_DIRS)
        if in_models:
            if Path(rel).name not in ("model.json", "model.joblib"):
                excluded.append((rel, "model dir file other than model.json/model.joblib"))
                continue
            selected[rel] = SRC / rel
            continue
        if rel.startswith("results/release/") and "predictions" in rel:
            excluded.append((rel, "release predictions"))
            continue
        if rel.startswith(tuple(EXPLICIT_RESULT_EXCLUDES)):
            if rel.endswith((".tsv.gz", ".npz")) or not rel.endswith(RESULT_EXTS):
                excluded.append((rel, "explicit exclusion (results/review/tierB_gain per-pair data)"))
                continue
        if rel.startswith(tuple(BLOCK_PREFIXES)) or rel.startswith(tuple(SKIP_PREFIXES)):
            if rel.startswith(tuple(BLOCK_PREFIXES)):
                n_tracked_unlisted += 1
            continue
        if rel in ROOT_DOCS or rel in ROOT_FILES or rel in DOC_FILES:
            selected[rel] = SRC / rel
        elif rel.startswith("results/"):
            if not rel.endswith(RESULT_EXTS):
                excluded.append((rel, "results file type not in {json, md, tsv}"))
                continue
            size = (SRC / rel).stat().st_size
            if size > RESULT_MAX_BYTES:
                excluded.append((rel, f"larger than 1.5 MB ({size / 1e6:.2f} MB)"))
                continue
            k = n_ids((SRC / rel).read_bytes())
            if k > ID_LIMIT:
                excluded.append((rel, f"{k} gene identifiers (> {ID_LIMIT}): per-pair/per-gene rows"))
                continue
            selected[rel] = SRC / rel
        elif rel.startswith(tuple(DIR_PREFIXES)):
            if low.endswith(BLOCK_EXTS):
                if rel.startswith("tests/") and fixture_generator(rel):
                    if check_fixture(rel, excluded):
                        selected[rel] = SRC / rel
                else:
                    excluded.append((rel, "blocked file type"))
                continue
            selected[rel] = SRC / rel
        else:
            n_tracked_unlisted += 1

    # docs/briefing: whatever is on disk now (the briefing is written separately)
    bdir = SRC / BRIEFING
    if bdir.is_dir():
        for p in sorted(bdir.rglob("*")):
            if not p.is_file() or p.name == ".DS_Store" or "__pycache__" in p.parts:
                continue
            rel = p.relative_to(SRC).as_posix()
            if p.suffix.lower() in BLOCK_EXTS or p.stat().st_size > 16 * 1024 * 1024:
                excluded.append((rel, "briefing file of blocked type or > 16 MB"))
                continue
            if p.suffix.lower() in TEXT_EXTS and n_ids(p.read_bytes()) > ID_LIMIT:
                notes.append(f"{rel}: more than {ID_LIMIT} gene identifiers in a briefing file (copied; PI to review)")
            selected[rel] = p
    else:
        notes.append("docs/briefing/ not found: nothing copied for it")

    # model verification (before copying anything)
    for m in MODEL_DIRS:
        for f, chk in (("model.json", check_model_json), ("model.joblib", check_joblib)):
            rel = f"{m}/{f}"
            if rel in selected:
                bad = chk(SRC / rel)
                if bad:
                    del selected[rel]
                    excluded.append((rel, "model verification failed: " + "; ".join(bad)))
            else:
                notes.append(f"{rel} missing")
    if f"results/release/model/model.json" in selected:
        selected["mvp/model/model.json"] = SRC / "results/release/model/model.json"
    if f"results/release/model/model.joblib" in selected:
        selected["mvp/model/model.joblib"] = SRC / "results/release/model/model.joblib"

    # build the export tree
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    manifest: list[tuple[str, str, str]] = []
    n_sanitised = 0
    json_failures = []
    for rel in sorted(selected):
        src, dst = selected[rel], OUT / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        changed = False
        raw = src.read_bytes()
        if src.suffix.lower() in TEXT_EXTS and b"\0" not in raw:
            try:
                txt = raw.decode("utf-8")
                new = redact_pair_rows(rel, sanitise_text(txt))
                if src.suffix.lower() == ".json":
                    new = redact_json_keys(new)
                changed = new != txt
                raw = new.encode("utf-8")
            except UnicodeDecodeError:
                pass
        dst.write_bytes(raw)
        shutil.copymode(src, dst)
        if changed:
            n_sanitised += 1
        if dst.suffix == ".json":
            try:
                json.loads(raw)
            except Exception as e:
                json_failures.append((rel, str(e)))
        manifest.append((rel, sha256_file(dst), "yes" if changed else "no"))

    # generated files (templates filled from the project's own config; nothing invented)
    import yaml

    man = yaml.safe_load((SRC / "config/data_manifest.yaml").read_text())
    entries = man.get("entries") or man.get("files") or next(v for v in man.values() if isinstance(v, list))
    ent = next(e for e in entries if e.get("id") == "uniprot_2026_03_yeast_fasta")
    lock = json.loads((SRC / "config/data_manifest.lock.json").read_text())["uniprot_2026_03_yeast_fasta"]
    generated = {
        ".gitignore": (TEMPLATES / "gitignore.txt").read_text(),
        "README.md": (TEMPLATES / "README_public.md").read_text(),
        "docs/USAGE.md": sanitise_text((SRC / "README.md").read_text()),  # the tool README becomes usage docs
        "mvp/README.md": (TEMPLATES / "mvp_README.md").read_text(),
        "mvp/run_example.sh": (TEMPLATES / "mvp_run_example.sh").read_text()
        .replace("@@YEAST_URL@@", ent["url"]).replace("@@YEAST_SHA256@@", lock["sha256"])
        .replace("@@YEAST_FILENAME@@", ent["filename"]),
    }
    for rel, txt in generated.items():
        dst = OUT / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(txt)
        if rel.endswith(".sh"):
            dst.chmod(0o755)
        manifest.append((rel, sha256_file(dst), "no"))
    manifest.sort()
    (OUT / "EXPORT_MANIFEST.tsv").write_text(
        "path\tsha256\tsanitised\n" + "".join(f"{a}\t{b}\t{c}\n" for a, b, c in manifest))

    # summary
    by_top: dict[str, int] = {}
    for rel, _, _ in manifest:
        top = rel.split("/")[0] if "/" in rel else "(root)"
        by_top[top] = by_top.get(top, 0) + 1
    size = sum(p.stat().st_size for p in OUT.rglob("*") if p.is_file())
    print(f"export: {OUT}")
    print(f"files copied/generated: {len(manifest)} (+ EXPORT_MANIFEST.tsv); sanitised: {n_sanitised}; "
          f"size {size / 1e6:.2f} MB")
    print("by top-level dir:", ", ".join(f"{k}={v}" for k, v in sorted(by_top.items())))
    print(f"tracked files outside the allowlist, not copied (data/, bin/, results/release/, ...): {n_tracked_unlisted}")
    print(f"exclusions inside the allowlist: {len(excluded)}")
    for rel, why in sorted(excluded):
        print(f"  EXCLUDED {rel}: {why}")
    for n in notes:
        print("  NOTE", n)
    print(f"JSON parse failures after sanitising: {len(json_failures)}")
    for rel, e in json_failures:
        print("  FAIL", rel, e)

    if args.commit:
        env_ok = subprocess.run(["git", "config", "--global", "user.email"], capture_output=True, text=True).stdout.strip()
        if not env_ok:
            print("no global git identity: not committing")
            return 1
        subprocess.run(["git", "init", "-q", "-b", "main"], cwd=OUT, check=True)
        subprocess.run(["git", "add", "-A"], cwd=OUT, check=True)
        subprocess.run(["git", "commit", "-q", "-m",
                        "Initial public export of ai-hackathon-ParalogSL\n\n"
                        "Sanitised export built by scripts/release/export_repo.py from the development repository "
                        "(no history, no data, no labels).\n\n"
                        "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"], cwd=OUT, check=True)
        tracked_out = set(subprocess.run(["git", "ls-files"], cwd=OUT, capture_output=True, text=True,
                                         check=True).stdout.split("\n")) - {""}
        expect = {r for r, _, _ in manifest} | {"EXPORT_MANIFEST.tsv"}
        print(f"git: committed {len(tracked_out)} files; in export but ignored by .gitignore: "
              f"{sorted(expect - tracked_out)}; tracked but not in manifest: {sorted(tracked_out - expect)}")
    return 1 if json_failures else 0


if __name__ == "__main__":
    sys.exit(main())
