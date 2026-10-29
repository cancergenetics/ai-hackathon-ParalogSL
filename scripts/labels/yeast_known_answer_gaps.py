#!/usr/bin/env python
"""Write the yeast entries of results/labels/known_answer_gaps.json (merged with other species' entries).

Recomputes the observed values of the yeast known-answer checks that did not reproduce (xfail(strict=True) in
tests/test_label_known_answers_yeast.py) plus per-file diagnostics. Expected values and tolerances are those
pre-registered in config/labels_yeast.yaml `known_answers`; nothing here changes a threshold.
Run: .venv/bin/python scripts/labels/yeast_known_answer_gaps.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import yeast_parsers as P  # noqa: E402

ROOT = HERE.parents[1]
CFG = yaml.safe_load((ROOT / "config" / "labels_yeast.yaml").read_text())
OUT = ROOT / CFG["qc"]["outputs"]["known_answer_gaps"]
T = "tests/test_label_known_answers_yeast.py::"


def main() -> None:
    c = CFG["sources"]["costanzo2016"]
    files = P.ensure_unpacked(ROOT / c["zip"], ROOT / c["unpack_to"], c["files_labels"] + c["files_known_answer_only"])
    per = {}
    for name in c["files_labels"] + c["files_known_answer_only"]:
        r = {"stringent": 0, "int_neg": 0, "int_pos": 0, "dmf_min_all_rows": np.inf, "dmf_min_sn_dma_rows": np.inf}
        for ch in P.iter_sga(files[name], 2_000_000):
            e, p, d = ch["eps"].to_numpy(float), ch["p"].to_numpy(float), ch["dmf"].to_numpy(float)
            with np.errstate(invalid="ignore"):
                r["stringent"] += int(((p < 0.05) & ((e < -0.12) | (e > 0.16))).sum())
                r["int_neg"] += int(((p < 0.05) & (e < -0.08)).sum())
                r["int_pos"] += int(((p < 0.05) & (e > 0.08)).sum())
            m = ((P.suffix_class(ch["query_strain"]) == "sn") & (P.suffix_class(ch["array_strain"]) == "dma")).fillna(False).to_numpy()
            r["dmf_min_all_rows"] = min(r["dmf_min_all_rows"], float(np.nanmin(d)))
            if m.any():
                r["dmf_min_sn_dma_rows"] = min(r["dmf_min_sn_dma_rows"], float(np.nanmin(d[m])))
        per[name] = {k: (None if v == np.inf else v) for k, v in r.items()}
    tot = {k: sum(per[f][k] for f in per) for k in ["stringent", "int_neg", "int_pos"]}
    ka = CFG["known_answers"]
    nxn = per["SGA_NxN.txt"]
    entries = {
        T + "test_costanzo_stringent_count_readme": {
            "species": "yeast", "expected": ka["stringent_all_files"]["value"], "tolerance_rel": 0.01,
            "observed": tot["stringent"], "rel_diff": round(tot["stringent"] / ka["stringent_all_files"]["value"] - 1, 4),
            "per_file": {f: per[f]["stringent"] for f in per},
            "note": ("Readme: HighestConfidence = raw (already reciprocal-filtered) + p<0.05 & (eps<-0.12 | eps>0.16) = "
                     "380,059 allele interactions. Applying that threshold to all four raw files gives more; which "
                     "screens/rows the readme count covers is not stated. Not used for any label rule; the label "
                     "thresholds are unaffected."),
        },
        T + "test_costanzo_intermediate_counts_readme": {
            "species": "yeast", "expected": {"neg": 550000, "pos": 350000}, "tolerance_rel": 0.05,
            "observed": {"neg": tot["int_neg"], "pos": tot["int_pos"]},
            "rel_diff": {"neg": round(tot["int_neg"] / 550000 - 1, 4), "pos": round(tot["int_pos"] / 350000 - 1, 4)},
            "per_file": {f: {"neg": per[f]["int_neg"], "pos": per[f]["int_pos"]} for f in per},
            "note": "Readme gives rounded '~550,000 negative and ~350,000 positive'; unit (allele pairs vs gene pairs) and file coverage not stated.",
        },
        T + "test_nxn_dmf_min": {
            "species": "yeast", "expected": ka["nxn_dmf_min"]["value"], "tolerance": "exact (4 dp)",
            "observed_all_nxn_rows": nxn["dmf_min_all_rows"], "observed_sn_dma_rows": nxn["dmf_min_sn_dma_rows"],
            "note": ("Pre-registered as the minimum over all NxN rows; the reviewer's value (data.md R2-2) is reproduced "
                     "on the sn x dma (null x null) rows only, so the reviewer's population was null x null. Test kept as "
                     "pre-registered (xfail strict); no label rule depends on it (DMF <= 0.5 is far from either floor)."),
        },
    }
    gaps = json.loads(OUT.read_text()) if OUT.exists() else {}
    if not isinstance(gaps, dict):
        raise SystemExit(f"{OUT} is not a JSON object; refusing to overwrite")
    gaps.update(entries)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(gaps, indent=1))
    print(json.dumps(entries, indent=1))


if __name__ == "__main__":
    main()
