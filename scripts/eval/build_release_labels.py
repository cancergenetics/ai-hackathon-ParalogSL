#!/usr/bin/env python
"""Build the release label tables data/labels/release/{human,yeast}_release.tsv (dev + test pairs).

The release model is refit on every labelled pair of the label universe (development AND held-out test
pairs of human and yeast). The tables are derived from data/labels/{human,yeast}_labels_all.tsv with the
same universe restriction and "labelled in any tier" rule as the dev tables (scripts/labels/make_splits.py
``load_species_labels``, config/splits.yaml); they are written in the dev-table format (original columns +
``group`` + ``fold_s<seed>``). This script never reads data/labels/sealed/.

Columns ``group`` / ``fold_s*``: rows that are in the dev table keep their dev values unchanged; the other
rows (the held-out test pairs) get group ``HELDOUT_TEST`` and fold -1 because their real group IDs live only
in sealed files. The release model ignores both columns (no cross-validation is run on release data).

Checks (the script exits non-zero if one fails, writing nothing):
  * release rows are a superset of the dev rows and every dev row is identical in all dev columns;
  * the number of release rows not in dev equals n_all_universe_labelled - n_dev (per species).
Writes data/labels/release/<sp>_release.tsv and results/release/release_labels_summary.json.

    .venv/bin/python scripts/eval/build_release_labels.py
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SPECIES = ("human", "yeast")
OUT_DIR = ROOT / "data" / "labels" / "release"
SUMMARY = ROOT / "results" / "release" / "release_labels_summary.json"


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_make_splits():
    spec = importlib.util.spec_from_file_location("make_splits", ROOT / "scripts" / "labels" / "make_splits.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # module level only defines functions (no sealed access)
    return mod


def main() -> int:
    ms = _load_make_splits()
    cfg = ms.load_cfg()
    for sp in SPECIES:
        assert "sealed" not in cfg["labels"][sp]["file"]
    out_tables, summary = {}, {}
    for sp in SPECIES:
        allq = ms.load_species_labels(cfg, sp)  # universe + labelled-in-any-tier, strings
        dev_path = ROOT / "data" / "labels" / "dev" / f"{sp}_dev.tsv"
        dev = pd.read_csv(dev_path, sep="\t", dtype=str, keep_default_na=False)
        dev_cols = list(dev.columns)
        fold_cols = [c for c in dev_cols if c.startswith("fold_s")]
        orig = [c for c in dev_cols if c not in fold_cols and c != "group"]
        if list(allq.columns[1:1 + len(orig)]) != orig:  # allq: species first, then original columns
            raise SystemExit(f"{sp}: original column order differs between all-labels and dev tables")
        rel = allq[orig].copy()
        key = ["gene_a", "gene_b"]
        if rel.duplicated(key).any() or dev.duplicated(key).any():
            raise SystemExit(f"{sp}: duplicate ordered pairs")
        m = rel.merge(dev[key + ["group"] + fold_cols], on=key, how="left", indicator=True)
        in_dev = (m["_merge"] == "both").to_numpy()
        m["group"] = m["group"].where(in_dev, "HELDOUT_TEST")
        for c in fold_cols:
            m[c] = m[c].where(in_dev, "-1")
        rel = m[dev_cols].reset_index(drop=True)
        # --- checks
        if in_dev.sum() != len(dev):
            raise SystemExit(f"{sp}: only {in_dev.sum()} of {len(dev)} dev pairs found in the release universe")
        sub = rel[in_dev].reset_index(drop=True)
        d2 = dev.sort_values(key).reset_index(drop=True)
        s2 = sub.sort_values(key).reset_index(drop=True)
        if not s2.equals(d2):
            raise SystemExit(f"{sp}: dev rows differ between release and dev tables")
        out_tables[sp] = rel
        summary[sp] = {"n_all_labelled_universe": int(len(rel)), "n_dev_rows": int(len(dev)),
                       "n_added_heldout_rows": int((~in_dev).sum()),
                       "dev_rows_identical": True, "release_superset_of_dev": True,
                       "all_labels_file": cfg["labels"][sp]["file"],
                       "all_labels_sha256": sha256_file(ROOT / cfg["labels"][sp]["file"]),
                       "dev_file_sha256": sha256_file(dev_path)}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    for sp, rel in out_tables.items():
        p = OUT_DIR / f"{sp}_release.tsv"
        rel.to_csv(p, sep="\t", index=False, lineterminator="\n")
        summary[sp]["release_file"] = str(p.relative_to(ROOT))
        summary[sp]["release_sha256"] = sha256_file(p)
    SUMMARY.write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
