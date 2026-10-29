"""Brief benchmark of the panel DB layout (Tier B): per-species MMseqs2 DBs vs one panel DB with
per-species taxonomy filtering (``--taxon-list``), same search settings (config/panel.yaml
tier_b.search, pinned scoring of config/caller.yaml).

Queries: a fixed random sample of the human proteome (seed = first entry of config/seeds.yaml).
Species: a few panel species of different sizes. Reports wall time per species search and the
per-(query, species) best-hit concordance at >= presence bits between the two layouts.

Usage: .venv/bin/python scripts/bench/panel_db_layout.py [--n 300] [--threads 8]
Writes results/bench/panel_db_layout.json.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from slparalog import fasta, panel, search  # noqa: E402
from slparalog.config import CONFIG_DIR, load_yaml  # noqa: E402


def lowpowermode():
    out = subprocess.run(["pmset", "-g"], capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "lowpowermode" in line:
            return line.split()[-1]
    return None


def best_hits(h: pd.DataFrame, presence: float) -> pd.DataFrame:
    h = h[h["bits"] >= presence].sort_values(["query", "bits", "target"], ascending=[True, False, True])
    return h.drop_duplicates("query")[["query", "target", "bits"]]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--species", nargs="*", default=["MOUSE", "DROME", "ARATH", "YEAST"])
    a = ap.parse_args()
    ccfg = load_yaml(CONFIG_DIR / "caller.yaml")
    pcfg = panel.load_panel_config()
    seed = int(load_yaml(CONFIG_DIR / "seeds.yaml")["seeds"][0])
    work = ROOT / "data/interim/bench_panel_layout"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    binary = search.resolve_binary(ccfg)
    mm = search.MMseqs(binary=binary, threads=a.threads, workdir=work)
    st = panel.build_dbs(pcfg, mm)

    seqs_by_id, _ = fasta.read_proteome(ROOT / "data/interim/proteomes/human.fasta")
    ids = fasta.canonical_order(seqs_by_id)
    rng = np.random.default_rng(seed)
    pick = np.sort(rng.choice(len(ids), size=a.n, replace=False))
    qseqs = [seqs_by_id[ids[k]] for k in pick]
    fasta.write_fasta(work / "q.fasta", ((str(k), s) for k, s in enumerate(qseqs)))
    search.createdb(mm, work / "q.fasta", "q")

    # layout B: one DB of the whole panel + taxonomy (minimal NCBI-style dump, one node per species)
    sp = panel.species_table(pcfg)
    tax = {c: int(t) for c, t in zip(sp["code"], sp["taxid"])}
    t0 = time.time()
    dump = work / "taxdump"
    dump.mkdir()
    (dump / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n" + "".join(
        f"{t}\t|\t1\t|\tspecies\t|\n" for t in sorted(set(tax.values())) if t != 1))
    (dump / "names.dmp").write_text("1\t|\troot\t|\t\t|\tscientific name\t|\n" + "".join(
        f"{t}\t|\tsp{t}\t|\t\t|\tscientific name\t|\n" for t in sorted(set(tax.values())) if t != 1))
    (dump / "merged.dmp").write_text("")
    (dump / "delnodes.dmp").write_text("")
    mm.run("createdb", str(panel._abs(pcfg["panel"]["fasta"])), "paneldb", "--shuffle", "0",
           "--dbtype", "1", threads=False)
    lk = pd.read_csv(work / "paneldb.lookup", sep="\t", header=None, names=["key", "acc", "file"])
    lk["taxid"] = lk["acc"].str.split("|").str[0].map(tax)
    lk[["acc", "taxid"]].to_csv(work / "taxmap.tsv", sep="\t", header=False, index=False)
    mm.run("createtaxdb", "paneldb", "tmp_tax", "--ncbi-tax-dump", str(dump),
           "--tax-mapping-file", str(work / "taxmap.tsv"))
    build_b = time.time() - t0

    presence = float(pcfg["tier_b"]["presence_min_bits"])
    scfg = pcfg["tier_b"]["search"]
    rows = []
    for code in a.species:
        t1 = time.time()
        ha, _ = panel.search_species(mm, ccfg, pcfg, "q", code, int(st["counts"][code]))
        ta = time.time() - t1
        t1 = time.time()
        mm.run("prefilter", "q", "paneldb", f"prefB_{code}", "--max-seqs", str(scfg["max_seqs"]),
               "-s", str(scfg["sensitivity"]), "--taxon-list", str(tax[code]),
               *panel.prefilter_scoring_args(ccfg))
        mm.run("align", "q", "paneldb", f"prefB_{code}", f"alnB_{code}", "-a", "--alignment-mode",
               str(ccfg["mmseqs"]["alignment_mode"]), "-e", str(scfg["evalue"]),
               *panel.align_scoring_args(ccfg))
        mm.run("convertalis", "q", "paneldb", f"alnB_{code}", f"alnB_{code}.tsv",
               "--format-output", ",".join(panel.PANEL_COLUMNS))
        tb = time.time() - t1
        hb = pd.read_csv(work / f"alnB_{code}.tsv", sep="\t", header=None, names=panel.PANEL_COLUMNS,
                         dtype={"target": str, "cigar": str}, keep_default_na=False)
        hb["target"] = hb["target"].str.split("|").str[1].astype(int)   # header '<code>|<index>' 
        ba, bb = best_hits(ha, presence), best_hits(hb, presence)
        m = ba.merge(bb, on="query", how="outer", suffixes=("_a", "_b"))
        both = m.dropna()
        rows.append({"species": code, "n_target": int(st["counts"][code]),
                     "per_species_db_s": round(ta, 2), "taxon_filter_s": round(tb, 2),
                     "cells_present_per_species": int(len(ba)), "cells_present_taxfilter": int(len(bb)),
                     "best_hit_concordance": round(float((both["target_a"] == both["target_b"]).mean()), 4)
                     if len(both) else None})
        print(rows[-1], flush=True)
    res = {"n_queries": a.n, "query_sample_seed": seed, "threads": a.threads,
           "settings": scfg, "taxonomy_db_build_s": round(build_b, 1), "rows": rows,
           "lowpowermode": lowpowermode(), "loadavg": os.getloadavg(),
           "decision_rule": "per-species DBs unless the taxonomy-filtered single DB is faster per species at equal results"}
    tot_a = sum(r["per_species_db_s"] for r in rows)
    tot_b = sum(r["taxon_filter_s"] for r in rows)
    res["total_per_species_db_s"], res["total_taxon_filter_s"] = round(tot_a, 1), round(tot_b, 1)
    res["chosen"] = "per_species" if tot_a <= tot_b else "taxon_filter"
    out = ROOT / "results/bench/panel_db_layout.json"
    out.write_text(json.dumps(res, indent=1))
    shutil.rmtree(work, ignore_errors=True)
    print(json.dumps({k: res[k] for k in ("total_per_species_db_s", "total_taxon_filter_s", "chosen")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
