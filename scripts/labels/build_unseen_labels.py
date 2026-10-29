#!/usr/bin/env python
"""Build SEALED unseen-species labels (DESIGN v3 §2.3, §6.3; config/labels_unseen.yaml).

Reads raw unseen-species sources (this script is allowlisted to), writes ONLY
  data/labels/sealed/unseen_<species>.tsv, data/labels/sealed/MANIFEST.sha256 (entries for these files),
  results/labels/unseen_label_counts.json (counts only), data/interim/tischler2006/*.txt (converted copies).
Prints counts only, never pairs or labels. Run: .venv/bin/python scripts/labels/build_unseen_labels.py
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
CFG_PATH = ROOT / "config" / "labels_unseen.yaml"
AGI_RE = re.compile(r"^AT[1-5CM]G\d{5}$")
WB_RE = re.compile(r"^WBGene\d{8}$")


# ----------------------------------------------------------------------------- helpers
def pair_key(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def fasta_accessions(path: Path) -> set[str]:
    acc = set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                parts = line[1:].split("|")
                if len(parts) >= 2:
                    acc.add(parts[1])
    return acc


def proteome_gene_ids(upid_taxid: str, xref: str, uniprot_dir: Path) -> set[str]:
    """Gene IDs (AGI via TAIR / WBGene via WormBase xref) of canonical accessions of the proteome."""
    accs = fasta_accessions(uniprot_dir / f"{upid_taxid}.fasta.gz")
    ids = set()
    with gzip.open(uniprot_dir / f"{upid_taxid}.idmapping.gz", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[1] == xref and f[0] in accs:
                ids.add(f[2])
    return ids


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow(["NA" if v is None else v for v in r])


def update_manifest(sealed_dir: Path, names: list[str]) -> None:
    """sha256sum-format manifest (names relative to sealed_dir); replaces only our entries, keeps others."""
    mp = sealed_dir / "MANIFEST.sha256"
    entries: dict[str, str] = {}
    if mp.exists():
        for line in mp.read_text().splitlines():
            if line.strip():
                h, n = line.split(None, 1)
                entries[n.strip().lstrip("*")] = h
    for n in names:
        entries[n] = sha256_file(sealed_dir / n)
    mp.write_text("".join(f"{h}  {n}\n" for n, h in sorted(entries.items())))


# ----------------------------------------------------------------------------- Arabidopsis
def read_cusack(path: Path, cols: dict) -> dict[tuple[str, str], tuple[str, float]]:
    """pair -> (class, proxy). Reads ONLY the whitelisted columns."""
    csv.field_size_limit(10**9)
    out = {}
    with open(path, newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        h = next(r)
        ip, ic, il = h.index(cols["pair_id"]), h.index(cols["class"]), h.index(cols["single_lethal_proxy"])
        for row in r:
            a, b = row[ip].split("_")
            lv = row[il]
            out[pair_key(a, b)] = (row[ic], float("nan") if lv in ("", "NA") else float(lv))
    return out


def build_arabidopsis(cfg: dict, uniprot_dir: Path) -> tuple[list[list], dict]:
    c = cfg["arabidopsis"]
    raw = ROOT / c["raw_dir"]
    d1 = read_cusack(raw / "Dataset_1.txt", c["columns_read"])
    d2 = read_cusack(raw / "Dataset_2.txt", c["columns_read"])
    keys = sorted(set(d1) | set(d2))
    genes_ok = proteome_gene_ids(c["proteome"], "TAIR", uniprot_dir)
    cnt = {
        "raw_dataset1": {"pos": sum(v[0] == "positive" for v in d1.values()), "neg": sum(v[0] == "negative" for v in d1.values())},
        "raw_dataset2": {"pos": sum(v[0] == "positive" for v in d2.values()), "neg": sum(v[0] == "negative" for v in d2.values())},
        "unique_pairs": len(keys),
    }
    conflicts = [k for k in set(d1) & set(d2) if d1[k][0] != d2[k][0]]
    cnt["class_conflicts_dropped"] = len(conflicts)
    bad_id = [k for k in keys if not (AGI_RE.match(k[0]) and AGI_RE.match(k[1]))]
    unmapped = [k for k in keys if k not in set(bad_id) and not (k[0] in genes_ok and k[1] in genes_ok)]
    cnt["malformed_id_pairs_dropped"] = len(bad_id)
    cnt["mapping_loss_pairs"] = len(unmapped)
    all_genes = {g for k in keys for g in k}
    cnt["genes_total"] = len(all_genes)
    cnt["genes_unmapped"] = len(all_genes - genes_ok)
    drop = set(conflicts) | set(bad_id) | set(unmapped)
    rows = []
    for k in keys:
        if k in drop:
            continue
        cl1 = d1.get(k, (None, float("nan")))[0]
        cl2 = d2.get(k, (None, float("nan")))[0]
        proxy = [v[k][1] for v in (d1, d2) if k in v]
        n_sl = int(max(proxy)) if proxy and proxy[0] == proxy[0] else None
        # strict tier
        strict = None if cl1 is None else (1 if cl1 == "positive" else 0)
        broad = None if cl2 is None else (1 if cl2 == "positive" else 0)
        note = ""
        if (cl1 == "negative" or cl2 == "negative") and n_sl is not None and n_sl >= 1:
            strict = broad = None
            note = "excluded_negative_lethal_single"
        elif cl1 is None and cl2 == "positive":
            note = "broad_only_positive"
        rows.append([k[0], k[1], strict, broad, n_sl, note])
    st = lambda i, v: sum(1 for r in rows if r[i] == v)  # noqa: E731
    cnt["tiers_after_recheck"] = {
        "strict": {"pos": st(2, 1), "neg": st(2, 0)},
        "broad": {"pos": st(3, 1), "neg": st(3, 0)},
    }
    cnt["negatives_excluded_lethal_single"] = sum(1 for r in rows if r[5] == "excluded_negative_lethal_single")
    cnt["positives_with_lethal_single_retained"] = {
        "strict": sum(1 for r in rows if r[2] == 1 and (r[4] or 0) >= 1),
        "broad": sum(1 for r in rows if r[3] == 1 and (r[4] or 0) >= 1),
    }
    cnt["tiers_before_recheck_mapped"] = {
        "strict": {"pos": sum(1 for r in rows if r[2] == 1), "neg": sum(1 for r in rows if r[2] == 0) + cnt["negatives_excluded_lethal_single"]},
        "broad": {"pos": sum(1 for r in rows if r[3] == 1), "neg": sum(1 for r in rows if r[3] == 0) + cnt["negatives_excluded_lethal_single"]},
    }
    cnt["proxy"] = "continuous_lethality_binary_number_in_pair_noTF (Cusack per-gene lethality annotation); Lloyd&Meinke Table S6 not on disk"
    return rows, cnt


# ----------------------------------------------------------------------------- C. elegans
def textutil_to_txt(doc: Path, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    if not out.exists() or out.stat().st_mtime < doc.stat().st_mtime:
        subprocess.run(["textutil", "-convert", "txt", "-output", str(out), str(doc)], check=True)
    return out


def parse_tischler_s2(txt: Path):
    """Return list of pairs ((gp_a, wb_a, lga_a), (gp_b, wb_b, lga_b)); 5 lines per gene row after 'Comments'."""
    lines = [ln.rstrip("\n") for ln in txt.read_text().splitlines()]
    body = lines[lines.index("Comments") + 1:]
    rows = [body[i:i + 5] for i in range(0, len(body) - 4, 5)]
    rows = [r for r in rows if r[0]]
    genes = [(r[2], r[3], r[4] == "Lga", r[0]) for r in rows]
    pairs = []
    for i in range(0, len(genes) - 1, 2):
        a, b = genes[i], genes[i + 1]
        if a[3] != b[3]:
            raise ValueError("S2 row pairing broken (ortholog IDs differ)")
        pairs.append((a[:3], b[:3]))
    return pairs


def parse_tischler_s4(txt: Path) -> set[frozenset]:
    lines = [ln.rstrip("\n") for ln in txt.read_text().splitlines()]
    body = lines[lines.index("Ks") + 1:]
    rows = [body[i:i + 6] for i in range(0, len(body), 6)]
    rows = [r for r in rows if r and r[0].startswith("WBGene")]
    return {frozenset((rows[i][0], rows[i + 1][0])) for i in range(0, len(rows) - 1, 2)}


def build_worm(cfg: dict, uniprot_dir: Path) -> tuple[list[list] | None, dict]:
    c = cfg["worm"]
    raw = ROOT / c["raw_dir"]
    interim = ROOT / "data" / "interim" / "tischler2006"
    s2 = textutil_to_txt(raw / "gb-2006-7-8-r69-S2.doc", interim / "S2.txt")
    s4 = textutil_to_txt(raw / "gb-2006-7-8-r69-S4.doc", interim / "S4.txt")
    pairs = parse_tischler_s2(s2)
    s4_pairs = parse_tischler_s4(s4)
    tested = [(a, b) for a, b in pairs if not (a[2] or b[2])]
    cnt = {
        "s2_pairs": len(pairs),
        "s2_pairs_lga_excluded": len(pairs) - len(tested),
        "tested_pairs": len(tested),
        "tested_in_s4": sum(frozenset((a[1], b[1])) in s4_pairs for a, b in tested),
        "s4_pairs": len(s4_pairs),
    }
    genes_ok = proteome_gene_ids(c["proteome"], "WormBase", uniprot_dir)
    bad = [(a, b) for a, b in tested if not (WB_RE.match(a[1]) and WB_RE.match(b[1]))]
    unm = [(a, b) for a, b in tested if (a, b) not in bad and not (a[1] in genes_ok and b[1] in genes_ok)]
    cnt["malformed_id_pairs"] = len(bad)
    cnt["mapping_loss_pairs"] = len(unm)
    cnt["self_pairs_same_wbgene"] = sum(a[1] == b[1] for a, b in tested)
    pos_file = raw / "positives_table4.tsv"
    if not pos_file.exists():
        cnt["status"] = "blocked_positives"
        cnt["positives"] = None
        cnt["note"] = "synthetic-phenotype positives are not in Additional files 2/4; need Table 4 (PMC1779603); no worm file written"
        return None, cnt
    gp2wb = {}
    for a, b in tested:
        gp2wb[a[0]], gp2wb[b[0]] = a[1], b[1]
    pos = set()
    with open(pos_file) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            pos.add(pair_key(gp2wb[row["gpname_a"]], gp2wb[row["gpname_b"]]))
    drop = {pair_key(a[1], b[1]) for a, b in bad + unm + [t for t in tested if t[0][1] == t[1][1]]}
    rows = []
    for a, b in tested:
        k = pair_key(a[1], b[1])
        if k in drop:
            continue
        rows.append([k[0], k[1], 1 if k in pos else 0])
    cnt["status"] = "written"
    cnt["positives"] = sum(r[2] == 1 for r in rows)
    cnt["negatives"] = sum(r[2] == 0 for r in rows)
    return rows, cnt


# ----------------------------------------------------------------------------- pombe
def pombe_status(cfg: dict) -> dict:
    d = ROOT / "data" / "raw" / "ryan2012_pombe"
    present = [p.name for p in d.glob("NIHMS380848-supplement-*") if p.suffix in (".xls", ".txt") and p.name != "MANUAL_DOWNLOAD.txt"]
    if not present:
        return {"status": "unavailable", "reason": "full E-MAP S-score matrix not on disk (manual download pending); dropped per DESIGN §2.3"}
    return {"status": "files_present_parser_not_implemented", "n_files": len(present)}


# ----------------------------------------------------------------------------- main
def main() -> int:
    cfg = yaml.safe_load(CFG_PATH.read_text())
    uniprot_dir = ROOT / "data" / "raw" / "uniprot_2026_03"
    sealed = ROOT / cfg["outputs"]["sealed_dir"]
    counts = {"config_version": cfg["version"], "uniprot_release": cfg["release"]}
    written = []

    rows, counts["arabidopsis"] = build_arabidopsis(cfg, uniprot_dir)
    fn = cfg["outputs"]["files"]["arabidopsis"]
    write_tsv(sealed / fn, ["gene_a", "gene_b", "label_strict", "label_broad", "single_lethal_n", "note"], rows)
    written.append(fn)

    wrows, counts["worm"] = build_worm(cfg, uniprot_dir)
    if wrows is not None:
        fn = cfg["outputs"]["files"]["worm"]
        write_tsv(sealed / fn, ["gene_a", "gene_b", "label"], wrows)
        written.append(fn)
    counts["pombe"] = pombe_status(cfg)
    counts["insects"] = {"status": "no_test_set", "note": "stated limitation (DESIGN §2.3)"}

    update_manifest(sealed, written)
    out = ROOT / cfg["outputs"]["counts_json"]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(counts, indent=2, sort_keys=True) + "\n")
    print(json.dumps(counts, indent=2, sort_keys=True))  # counts only
    return 0


if __name__ == "__main__":
    sys.exit(main())
