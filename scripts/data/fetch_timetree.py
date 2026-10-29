#!/usr/bin/env python
"""Fetch TimeTree 5 pairwise divergence times between each query species and each panel species
(manifest entry ``timetree5_pairwise``; evaluation only: Spearman validation of Tier B d_s).

API: https://timetree.org/api/pairwise/{taxid_a}/{taxid_b} -> CSV with header
taxon_a_id,taxon_b_id,scientific_name_a,scientific_name_b,all_total,precomputed_age,
precomputed_ci_low,precomputed_ci_high,adjusted_age (Mya). One request per pair.

Writes data/raw/timetree/timetree_pairwise.tsv (one row per (query, panel code); failures kept
with an error message), SOURCE.json (fetch time, URL template, counts) and SHA256SUMS.txt.

Usage: .venv/bin/python scripts/data/fetch_timetree.py
"""

from __future__ import annotations

import csv
import datetime
import hashlib
import io
import json
import os
import sys
import time

import requests

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "data", "raw", "timetree")
URL = "https://timetree.org/api/pairwise/{a}/{b}"
UA = "sl_run2-data-acquisition/1.0 (academic research; python-requests)"
# Query (training / unseen-test) species of the project: NCBI taxids of the proteomes in
# data/interim/proteomes (UniProt reference proteomes; see DATA.md).
QUERIES = {"human": 9606, "yeast": 559292, "arabidopsis": 3702, "worm": 6239}


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(ROOT, "config", "panel_species.tsv")) as fh:
        panel = list(csv.DictReader(fh, delimiter="\t"))
    rows, n_fail = [], 0
    started = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    sess = requests.Session()
    for qname, qtax in QUERIES.items():
        for p in panel:
            rec = {"query": qname, "query_taxid": qtax, "code": p["code"], "taxid": int(p["taxid"]),
                   "adjusted_age": "", "precomputed_age": "", "ci_low": "", "ci_high": "",
                   "n_studies": "", "name_a": "", "name_b": "", "error": ""}
            if int(p["taxid"]) == qtax:
                rec["error"] = "same taxid"
                rows.append(rec)
                continue
            err = None
            for attempt in range(4):
                try:
                    r = sess.get(URL.format(a=qtax, b=p["taxid"]), timeout=60, headers={"User-Agent": UA})
                    if r.status_code != 200:
                        err = f"HTTP {r.status_code}"
                    else:
                        lines = [ln for ln in r.text.strip().splitlines() if ln.strip()]
                        rd = list(csv.DictReader(io.StringIO("\n".join(lines))))
                        if not rd:
                            err = "empty"
                        else:
                            d = rd[-1]
                            rec.update(adjusted_age=d.get("adjusted_age", ""),
                                       precomputed_age=d.get("precomputed_age", ""),
                                       ci_low=d.get("precomputed_ci_low", ""),
                                       ci_high=d.get("precomputed_ci_high", ""),
                                       n_studies=d.get("all_total", ""),
                                       name_a=d.get("scientific_name_a", ""),
                                       name_b=d.get("scientific_name_b", ""))
                            err = None if rec["adjusted_age"] not in ("", "0") else "no age"
                    if err is None or err == "no age":
                        break
                except requests.RequestException as e:  # retry on network errors
                    err = type(e).__name__
                time.sleep(2 * (attempt + 1))
            if err:
                rec["error"] = err
                n_fail += 1
            rows.append(rec)
            time.sleep(0.2)
        print(qname, "done", flush=True)
    path = os.path.join(OUT, "timetree_pairwise.tsv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    src = {"id": "timetree5_pairwise", "url_template": URL, "fetched_at_start": started,
           "fetched_at_end": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "n_rows": len(rows), "n_failed": n_fail, "queries": QUERIES,
           "file": "timetree_pairwise.tsv", "sha256": sha256(path),
           "field_used": "adjusted_age (Mya)"}
    # same layout as scripts/data/download.py SOURCE.json files (read by update_data_md.py)
    src_json = {"group": "timetree", "files": [{
        "id": "timetree5_pairwise", "filename": "timetree_pairwise.tsv", "url": URL,
        "date": src["fetched_at_end"], "sha256": src["sha256"], "bytes": os.path.getsize(path),
        "status": "ok", **{k: src[k] for k in ("fetched_at_start", "n_rows", "n_failed", "queries",
                                                "field_used")}}]}
    with open(os.path.join(OUT, "SOURCE.json"), "w") as fh:
        json.dump(src_json, fh, indent=1)
    with open(os.path.join(OUT, "SHA256SUMS.txt"), "w") as fh:
        fh.write(f"{src['sha256']}  timetree_pairwise.tsv\n")
    print(json.dumps({k: src[k] for k in ("n_rows", "n_failed")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
