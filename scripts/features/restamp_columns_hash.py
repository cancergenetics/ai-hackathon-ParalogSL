"""Add ``columns_hash`` to the sidecar meta of existing feature tables without touching the parquet files.

Feature tables built before the column-definition hash existed carry only the whole-registry ``registry_hash``,
which changes whenever any description or any other tier is edited (so Tier A tables are refused after Tier B/C
were registered although no Tier A column definition changed). For each table
``data/interim/features/<species>/tier{A,B,C}.parquet`` (species human, yeast, arabidopsis, worm) this script:

1. Finds the commit the table was built at: ``git.head`` (or ``git.commit``) of its meta if present and clean;
   otherwise the commit that introduced the meta file (``git log --diff-filter=A``); the meta files are not
   tracked here, so the last resort is the newest commit at or before the meta file's modification time.
   Every fallback is reported.
2. Reads ``src/slparalog/features/registry.py`` at that commit with ``git show`` (our own repository history),
   imports it from a temp file (after a scan that it imports nothing but the standard library and ``..guards``,
   which is stubbed), and computes the column-definition hash of the table's registered columns.
3. Cross-checks that the old registry's own ``registry_hash`` equals the ``registry_hash`` in the table meta.
   If it does not, the commit identification is unreliable (e.g. the table was built from an uncommitted working
   tree) and the registry version whose hash does equal the meta's is used instead (reported as such); if there
   is none, the table is not restamped.
4. Restamps only if that hash equals the current registry's hash for the same columns AND the table's columns
   equal the registered columns of its tier (plus gene_a/gene_b). Then ``columns_hash`` and a ``restamped``
   provenance record are added to ``<table>.parquet.meta.json``. The parquet is never opened for writing;
   ``slparalog.model.table_meta`` reads ``columns_hash`` from the sidecar.

Differences are reported and the table is left alone. Summary: results/features/restamp_report.json.
Idempotent: an already stamped table is verified, not rewritten.
"""

from __future__ import annotations

import argparse
import datetime as dt
import importlib
import json
import re
import subprocess
import sys
import tempfile
from dataclasses import fields
from pathlib import Path

import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from slparalog.features import registry  # noqa: E402

SPECIES = ("human", "yeast", "arabidopsis", "worm")
TIERS = ("A", "B", "C")
FEATURE_DIR = ROOT / "data" / "interim" / "features"
REPORT = ROOT / "results" / "features" / "restamp_report.json"
REG_PATH = "src/slparalog/features/registry.py"
ALLOWED_IMPORT = re.compile(r"^\s*(from\s+(__future__|dataclasses|\.\.guards)\s+import|import\s+(hashlib|json)\b|"
                            r"from\s+(typing)\s+import)")
ID_COLS = {"gene_a", "gene_b"}


def git(*args: str, check: bool = True) -> str:
    r = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)
    if check and r.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)} failed: {r.stderr.strip()}")
    return r.stdout


def registry_source(commit: str) -> str | None:
    r = subprocess.run(["git", "show", f"{commit}:{REG_PATH}"], cwd=ROOT, capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None


def load_old_registry(commit: str, source: str):
    """Import a historic registry.py from a temp package (``..guards`` stubbed; imports scanned first)."""
    bad = [ln for ln in source.splitlines() if re.match(r"^\s*(import|from)\s", ln) and not ALLOWED_IMPORT.match(ln)]
    if bad:
        raise RuntimeError(f"registry.py at {commit[:8]} has unexpected imports {bad}; refusing to import it")
    pkg = f"slold_{commit[:12]}"
    with tempfile.TemporaryDirectory() as td:
        base = Path(td) / pkg
        (base / "features").mkdir(parents=True)
        (base / "__init__.py").write_text("")
        (base / "features" / "__init__.py").write_text("")
        (base / "guards.py").write_text("import re\nFORBIDDEN_NAME_RE = re.compile('$^')\n")
        (base / "features" / "registry.py").write_text(source)
        sys.path.insert(0, td)
        try:
            return importlib.import_module(f"{pkg}.features.registry")
        finally:
            sys.path.remove(td)


def as_current_features(old) -> tuple[list[registry.Feature], list[str]]:
    """Old Feature objects as current ``Feature``s (missing fields take the current defaults, and are reported)."""
    cur_fields = [f.name for f in fields(registry.Feature)]
    missing = sorted({f for o in old.FEATURES for f in cur_fields if not hasattr(o, f)})
    feats = []
    for o in old.FEATURES:
        kw = {f.name: getattr(o, f.name) for f in fields(registry.Feature) if hasattr(o, f.name)}
        kw.setdefault("description", "")
        feats.append(registry.Feature(**kw))
    return feats, missing


def meta_commit(meta: dict) -> tuple[str | None, str]:
    """Commit recorded in the table meta, if any and clean; second item says why not / where it came from."""
    g = meta.get("git")
    if not isinstance(g, dict):
        return None, "meta has no git record"
    commit = g.get("head") or g.get("commit")
    dirty = g.get("src_dirty") if "src_dirty" in g else g.get("dirty")
    if not commit or commit == "unknown":
        return None, "meta git commit missing"
    if dirty is not False:
        return None, f"meta git commit {commit[:8]} recorded with a dirty/unknown tree ({dirty!r})"
    return commit, "meta git.head" if "head" in g else "meta git.commit"


def fallback_commit(meta_path: Path) -> tuple[str | None, str]:
    rel = str(meta_path.relative_to(ROOT))
    out = git("log", "--diff-filter=A", "--format=%H", "--", rel, check=False).split()
    if out:
        return out[-1], "commit that introduced the meta file (git log --diff-filter=A)"
    ts = dt.datetime.fromtimestamp(meta_path.stat().st_mtime, dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    c = git("rev-list", "-1", f"--before={ts}", "HEAD", check=False).strip()
    if c:
        return c, f"meta file untracked (data/interim is ignored): newest commit at or before its mtime {ts}"
    return None, "no commit found"


def registry_versions() -> list[str]:
    """Commits that changed registry.py, newest first."""
    return git("log", "--format=%H", "--", REG_PATH).split()


def analyse(species: str, tier: str) -> dict:
    pq_path = FEATURE_DIR / species / f"tier{tier}.parquet"
    meta_path = Path(str(pq_path) + ".meta.json")
    rec: dict = {"table": str(pq_path.relative_to(ROOT)), "species": species, "tier": tier}
    if not pq_path.exists() or not meta_path.exists():
        return {**rec, "action": "not_restamped", "reasons": ["table or meta.json missing"]}
    meta = json.loads(meta_path.read_text())
    cols = pq.read_schema(pq_path).names  # schema only; the parquet is never modified
    reg_cols = [c for c in cols if c in registry.BY_NAME]
    extras = sorted(c for c in cols if c not in registry.BY_NAME and c not in ID_COLS)
    tier_cols = {"A": registry.TIER_A_COLUMNS, "B": registry.TIER_B_COLUMNS, "C": registry.TIER_C_COLUMNS}[tier]
    cur_hash = registry.columns_hash(reg_cols)
    rec.update(n_columns=len(cols), n_registered_columns=len(reg_cols), meta_registry_hash=meta.get("registry_hash"),
               current_registry_hash=registry.registry_hash(), current_columns_hash=cur_hash)
    reasons: list[str] = []

    prior = meta.get("columns_hash")
    if prior is not None:  # built with a columns_hash, or restamped by an earlier run: verify only
        ok = prior == cur_hash
        rec.update(action="already_stamped" if ok else "not_restamped", existing_columns_hash=prior,
                   existing_restamped_record=meta.get("restamped"),
                   reasons=[] if ok else ["existing columns_hash differs from the current definitions of the table's columns"])
        return rec

    columns_equal_tier = set(reg_cols) == set(tier_cols) and not extras and ID_COLS <= set(cols)
    rec["columns_equal_registered_tier"] = columns_equal_tier
    if not columns_equal_tier:
        reasons.append(f"table columns differ from registered Tier {tier} columns "
                       f"(missing {sorted(set(tier_cols) - set(reg_cols))}, extra {sorted(set(reg_cols) - set(tier_cols))}, "
                       f"unregistered {extras})")

    # --- which commit's registry does the table correspond to
    commit, how = meta_commit(meta)
    rec["meta_commit_note"] = how
    if commit is None:
        commit, how = fallback_commit(meta_path)
        rec["fallback"] = how
        if "git" in meta:
            rec["fallback_reason"] = rec["meta_commit_note"]
    rec["basis_commit"], rec["commit_source"] = commit, how
    old = old_feats = None
    src = registry_source(commit) if commit else None
    if src is None:
        reasons.append(f"registry.py not readable at commit {commit}")
    else:
        old = load_old_registry(commit, src)
        rec["commit_registry_hash"] = old.registry_hash()
        if old.registry_hash() != meta.get("registry_hash"):
            rec["commit_registry_hash_matches_meta"] = False
            alt = None
            for c in registry_versions():  # the table was built from a tree whose registry we can identify by hash
                s2 = registry_source(c)
                m2 = load_old_registry(c, s2) if s2 else None
                if m2 is not None and m2.registry_hash() == meta.get("registry_hash"):
                    alt, old, src = c, m2, s2
                    break
            if alt is None:
                old = None
                reasons.append(f"registry at {commit[:8]} hashes to {rec['commit_registry_hash'][:12]}, meta says "
                               f"{str(meta.get('registry_hash'))[:12]}, and no registry.py version in history matches")
            else:
                rec.update(basis_commit=alt, commit_source=f"registry version whose registry_hash equals the meta's "
                                                           f"(commit {commit[:8]} from '{how}' did not)")
        else:
            rec["commit_registry_hash_matches_meta"] = True
    if old is not None:
        old_feats, missing_fields = as_current_features(old)
        rec["old_feature_fields_missing"] = missing_fields
        absent = [c for c in reg_cols if c not in {f.name for f in old_feats}]
        if absent:
            reasons.append(f"columns absent from the registry at {rec['basis_commit'][:8]}: {absent}")
        else:
            rec["old_columns_hash"] = registry.columns_hash(reg_cols, old_feats)
            rec["definitions_equal"] = rec["old_columns_hash"] == cur_hash
            if not rec["definitions_equal"]:
                diff = [c for c in reg_cols
                        if registry.column_definition(next(f for f in old_feats if f.name == c))
                        != registry.column_definition(registry.BY_NAME[c])]
                rec["changed_columns"] = diff
                reasons.append(f"column definitions changed since {rec['basis_commit'][:8]}: {diff}")

    if reasons:
        rec.update(action="not_restamped", reasons=reasons)
        return rec
    rec.update(action="restamp", reasons=[])
    rec["_meta"], rec["_meta_path"] = meta, meta_path
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dry-run", action="store_true", help="report only; do not write any meta.json")
    a = ap.parse_args()
    prev = {}
    if REPORT.exists():  # a verify-only re-run must not erase the record of what an earlier run restamped
        prev = {t["table"]: t for t in json.loads(REPORT.read_text()).get("tables", [])}
    results = []
    for sp in SPECIES:
        for t in TIERS:
            r = analyse(sp, t)
            if r["action"] == "restamp":
                meta, mp = r.pop("_meta"), r.pop("_meta_path")
                if not a.dry_run:
                    meta["columns_hash"] = r["current_columns_hash"]
                    meta["restamped"] = {
                        "utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                        "script": "scripts/features/restamp_columns_hash.py",
                        "basis_commit": r["basis_commit"], "commit_source": r["commit_source"],
                        "table_registry_hash": r["meta_registry_hash"],
                        "current_registry_hash": r["current_registry_hash"],
                        "note": "columns_hash added after verifying that the definitions (name, source, dtype, kind, "
                                "tier, pct) of this table's columns at the basis commit equal the current ones; "
                                "parquet unchanged",
                    }
                    mp.write_text(json.dumps(meta, indent=1, default=str))
                r["action"] = "restamped" if not a.dry_run else "would_restamp"
            if r["action"] == "already_stamped" and prev.get(r["table"], {}).get("action") == "restamped":
                r = {**prev[r["table"]], "verified_again_columns_hash_equals_current": True}
            results.append(r)
    summary = {k: sorted(r["table"] for r in results if r["action"] == k)
               for k in ("restamped", "would_restamp", "already_stamped", "not_restamped")}
    report = {"script": "scripts/features/restamp_columns_hash.py", "head": git("rev-parse", "HEAD").strip(),
              "dry_run": a.dry_run, "current_registry_hash": registry.registry_hash(),
              "summary": {k: v for k, v in summary.items() if v or k == "not_restamped"}, "tables": results}
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=1))
    for r in results:
        print(f"{r['table']}: {r['action']}" + (f" ({'; '.join(r['reasons'])})" if r.get("reasons") else "")
              + (f" [{r.get('commit_source')}]" if r["action"] in ("restamped", "would_restamp") else ""))
    return 1 if summary["not_restamped"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
