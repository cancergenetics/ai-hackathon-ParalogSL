#!/usr/bin/env python
"""Restrict yeast pair labels to the paralog label universe (R2-4; config/labels_yeast.yaml `universe`).

Re-runnable separate step. Reads data/interim/labels/yeast_pair_labels_unrestricted.tsv (from
build_yeast_labels.py) and the Ensembl Compara paralog universe data/interim/universe/yeast_paralog_pairs.tsv
(built by another step). If the universe file is absent, nothing is written to data/labels/ and the QC JSON
records universe.status = 'pending'; re-run this script once the file exists:
    .venv/bin/python scripts/labels/apply_universe.py
Writes data/labels/yeast_labels_all.tsv and updates the `universe` section of results/labels/yeast_label_qc.json
(counts before/after, per-tier counts within the universe, universe pairs never tested).
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
ROOT = HERE.parents[1]
CFG = yaml.safe_load((ROOT / "config" / "labels_yeast.yaml").read_text())


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_universe(path: Path, cols: list[str]) -> pd.DataFrame:
    hdr = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    if not set(cols) <= set(hdr):
        raise ValueError(f"universe file {path} lacks columns {cols}; has {hdr}")
    u = pd.read_csv(path, sep="\t", usecols=cols, dtype=str)
    a, b = u[cols[0]].str.strip(), u[cols[1]].str.strip()
    u = pd.DataFrame({"gene_a": a.where(a < b, b), "gene_b": b.where(a < b, a)})
    u = u[u["gene_a"] != u["gene_b"]].drop_duplicates()
    return u


def main() -> int:
    out = CFG["qc"]["outputs"]
    qc_path = ROOT / out["qc_json"]
    qc = json.loads(qc_path.read_text()) if qc_path.exists() else {}
    upath = ROOT / CFG["universe"]["file"]
    if not upath.exists():
        qc["universe"] = {"status": "pending", "file": CFG["universe"]["file"],
                          "note": "universe file absent; data/labels/yeast_labels_all.tsv not written; re-run apply_universe.py"}
        qc_path.write_text(json.dumps(qc, indent=1, default=str))
        print(f"universe file {upath} absent: nothing written (status pending)")
        return 0
    import build_yeast_labels as B  # noqa: E402  (shared QC code)

    u = read_universe(upath, list(CFG["universe"]["pair_columns"]))
    parts, n_before, genes_tested = [], 0, set()
    for ch in pd.read_csv(ROOT / out["pairs_unrestricted"], sep="\t", chunksize=1_000_000, keep_default_na=False,
                          na_values={"min_smf": [""], "linked": [""], "same_chrom_gap_bp": [""]},
                          dtype={"dean_call": str, "source": str, "evidence": str, "license_note": str}):
        n_before += len(ch)
        m = ch["n_meas"] > 0
        genes_tested |= set(ch.loc[m, "gene_a"]) | set(ch.loc[m, "gene_b"])
        parts.append(ch.merge(u, on=["gene_a", "gene_b"], how="inner"))
    r = pd.concat(parts, ignore_index=True)
    untested_u = u.merge(r[["gene_a", "gene_b"]], on=["gene_a", "gene_b"], how="left", indicator=True)
    untested_u = untested_u[untested_u["_merge"] == "left_only"]
    dest = ROOT / out["labels"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    r.to_csv(dest, sep="\t", index=False)
    qc["universe"] = {
        "status": "applied", "file": CFG["universe"]["file"], "file_sha256": sha256_file(upath),
        "n_universe_pairs": int(len(u)), "n_pairs_before": int(n_before), "n_pairs_after": int(len(r)),
        "universe_pairs_not_tested_unlabelled": int(len(untested_u)),
        "universe_pairs_not_tested_both_genes_in_sga_measured_pairs": int((untested_u["gene_a"].isin(genes_tested)
                                                                     & untested_u["gene_b"].isin(genes_tested)).sum()),
        "restricted": B.restricted_report(CFG, r),
        "labels_file_sha256": sha256_file(dest),
    }
    qc_path.write_text(json.dumps(qc, indent=1, default=str))
    print(f"universe applied: {n_before} -> {len(r)} pairs; wrote {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
