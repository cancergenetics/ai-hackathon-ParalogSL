#!/usr/bin/env python
"""Parse Table 4 (synthetic RNAi phenotype pairs) of Tischler et al. 2006 (PMC1779603) into positives_table4.tsv.

Input : data/raw/tischler2006/PMC1779603_fulltext.xml (manifest id tischler2006_fulltext_xml; JATS XML via Europe PMC)
Output: data/raw/tischler2006/positives_table4.tsv, columns gpname_a, gpname_b (GPName = sequence names as in Additional file 2),
        the format documented in config/labels_unseen.yaml (worm.positive.positives_source).

Table 4 lists pairs by C. elegans gene names (some are former names), whereas the label builder keys on Additional file 2's
GPName. The map is derived from data only (no hand-entered gene knowledge): name -> WBGene via the Alliance/WormBase gene table
(symbol or synonym, case-insensitive; manifest id alliance_bgi_wb) -> the S2 duplicate pair (consecutive rows sharing an ortholog ID)
with exactly that WBGene pair -> its GPNames. A name may map to several WBGenes (reused names); the S2 pair must then be unique.
Fallback (counted separately, reported in results/labels/tischler_table4_extraction.json): if one name reaches no S2 pair but the
other name reaches a gene that sits in exactly one S2 pair, and the unreached name resolves to no gene present in S2 at all, the
S2 partner is taken (covers a former name that the gene table no longer carries). Anything else is unresolved: the script writes
nothing and exits non-zero.
Prints counts only, never gene names or pairs.

Run: .venv/bin/python -I scripts/data/extract_tischler_table4.py
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "tischler2006"
XML = RAW / "PMC1779603_fulltext.xml"
OUT = RAW / "positives_table4.tsv"
BGI = RAW / "BGI_WB.json.gz"
COUNTS = ROOT / "results" / "labels" / "tischler_table4_extraction.json"
EXPECTED = 16  # Table 4 rows, per config/labels_unseen.yaml (also checked against the parsed table, not assumed)


def load_builder():
    spec = importlib.util.spec_from_file_location("build_unseen_labels", ROOT / "scripts" / "labels" / "build_unseen_labels.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def text_of(el) -> str:
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip()


def parse_table4(xml_path: Path) -> list[tuple[str, str]]:
    root = ET.parse(xml_path).getroot()
    tw = next((t for t in root.iter("table-wrap") if t.get("id") == "T4"), None)
    if tw is None or text_of(tw.find("label")) != "Table 4":
        raise SystemExit("Table 4 (id T4) not found in full text")
    pairs = []
    for tr in tw.find(".//tbody").iter("tr"):
        cells = tr.findall("td")
        if not cells:
            continue
        first = text_of(cells[0])
        parts = [p.strip() for p in first.split("+")]
        if len(parts) != 2 or not all(parts):
            raise SystemExit(f"unparseable Table 4 row label ({len(parts)} parts)")
        pairs.append((parts[0], parts[1]))
    return pairs


def name_map(bgi_path: Path) -> dict[str, set[str]]:
    """lower-cased symbol / synonym -> set of WBGene IDs (only these two fields are read)."""
    with gzip.open(bgi_path, "rt") as fh:
        data = json.load(fh)["data"]
    nm: dict[str, set[str]] = {}
    for g in data:
        wb = g["basicGeneticEntity"]["primaryId"].split(":", 1)[1]
        for n in [g["symbol"]] + g["basicGeneticEntity"].get("synonyms", []):
            nm.setdefault(n.lower(), set()).add(wb)
    return nm


def main() -> int:
    B = load_builder()
    pairs = parse_table4(XML)

    s2 = B.textutil_to_txt(RAW / "gb-2006-7-8-r69-S2.doc", ROOT / "data" / "interim" / "tischler2006" / "S2.txt")
    s2_pairs = B.parse_tischler_s2(s2)  # ((gp, wb, lga), (gp, wb, lga))
    index = {}
    for a, b in s2_pairs:
        index.setdefault(frozenset((a[1], b[1])), []).append((a, b))

    nm = name_map(BGI)
    s2_genes = {g[1] for p in s2_pairs for g in p}
    out_rows = []
    cnt = {"table4_rows": len(pairs), "resolved_direct": 0, "resolved_anchor_fallback": 0,
           "unresolved": 0, "lga_flagged_in_s2": 0, "source_sha256": B.sha256_file(XML), "name_table_sha256": B.sha256_file(BGI)}
    for na, nb in pairs:
        wa, wb = nm.get(na.lower(), set()), nm.get(nb.lower(), set())
        direct = {frozenset((x, y)) for x in wa for y in wb if frozenset((x, y)) in index}
        cand = None
        if len(direct) == 1:
            cand, kind = index[next(iter(direct))], "resolved_direct"
        elif not direct:
            # anchor fallback: exactly one side reaches S2 (via a gene in exactly one S2 pair), the other side reaches no S2 gene
            for known, other in ((wa, wb), (wb, wa)):
                k = [g for g in known if g in s2_genes]
                partners = {frozenset(h) for g in k for h in index if g in h}
                if len(k) == 1 and len(partners) == 1 and not (other & s2_genes) and len(index[next(iter(partners))]) == 1:
                    cand, kind = index[next(iter(partners))], "resolved_anchor_fallback"
        if cand is None or len(cand) != 1:
            cnt["unresolved"] += 1
            continue
        a, b = cand[0]
        cnt[kind] += 1
        cnt["lga_flagged_in_s2"] += bool(a[2] or b[2])
        out_rows.append((a[0], b[0]))
    cnt["unique_pairs"] = len({frozenset(r) for r in out_rows})
    COUNTS.parent.mkdir(parents=True, exist_ok=True)
    COUNTS.write_text(json.dumps(cnt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in cnt.items() if not k.endswith("sha256")}, sort_keys=True))
    if cnt["unresolved"] or cnt["unique_pairs"] != EXPECTED or len(out_rows) != EXPECTED or len(pairs) != EXPECTED:
        print(f"FAIL: expected {EXPECTED} parsed rows and {EXPECTED} unique resolved pairs; {OUT.name} not written")
        return 1
    OUT.write_text("gpname_a\tgpname_b\n" + "".join(f"{a}\t{b}\n" for a, b in out_rows))
    print(f"wrote {OUT.relative_to(ROOT)}: {len(out_rows)} pairs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
