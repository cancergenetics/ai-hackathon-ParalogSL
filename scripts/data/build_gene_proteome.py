#!/usr/bin/env python
"""Build per-gene proteomes (one FASTA record per gene) from the pinned UniProt reference proteomes (DESIGN v3 §2.1, R2-8).

Run:  .venv/bin/python scripts/data/build_gene_proteome.py [--species human yeast ...]

Inputs  data/raw/uniprot_2026_03/<UPID>_<taxid>.fasta.gz + .idmapping.gz   (UniProt 2026_03)
        human only: data/raw/hgnc/hgnc_complete_set_2026-10-06.tsv           (HGNC ID -> Ensembl gene)
                    data/raw/ensembl_116/Homo_sapiens.GRCh38.pep.all.fa.gz   (primary-assembly vs patch/haplotype ENSG)
Outputs data/interim/proteomes/<species>.fasta            one record per gene; header is the gene ID only
        data/interim/proteomes/<species>_gene2acc.tsv     gene -> accession with provenance
        results/data/proteome_<species>.json              QC (n genes, n accessions, n expanded, n dropped by reason, ...)

Rules (gene-ID convention: config/gene_ids.yaml; R2-8 adapted to the ENSG convention):
 1. Only accessions in the reference FASTA are used.
 2. Non-human species: gene IDs of an accession come from the release's .idmapping, first ID type in the species' priority
    list that yields an ID (yeast Gene_OrderedLocusName; arabidopsis TAIR > Araport > Gene_OrderedLocusName, upper-cased and
    restricted to the AGI pattern; pombe PomBase; worm WormBase; fly FlyBase). Worm/fly accessions without a MOD xref get
    the MOD ID that their Gene_ORFName carries on other accessions of the same idmapping file, if that is unique and
    unused ('orfname_xref').
 3. Human (genes are HGNC genes, IDs are their Ensembl gene IDs): for each HGNC ID of the accession (idmapping 'HGNC'), the
    gene ID is the HGNC complete set's ensembl_gene_id. UniProt's own Ensembl xrefs are NOT expanded directly, because they
    include patch/alt-haplotype copies, PAR-Y copies and read-through/novel overlapping models of the same gene (verified:
    1,159 accessions with 2 primary-or-patch ENSG but 1 HGNC). HGNC IDs without an Ensembl ID take an unused
    primary-assembly UniProt ENSG of the accession, if any. Accessions without HGNC ID take ONE UniProt ENSG
    (primary-assembly preferred, then smallest ID; logged).
 4. An accession carrying k >= 2 genes (identical-sequence genes collapsed by UniProt) is expanded into k records with the
    same sequence ('expanded').
 5. An accession with no gene ID keeps its sequence under the accession as record ID ('accession_fallback', logged).
 6. A gene reached from >= 2 accessions keeps one, chosen by: accession carries this gene only > Swiss-Prot > longer sequence >
    smaller accession (logged). Accessions left with no gene are dropped ('gene_assigned_to_other_accession').
"""
import argparse, collections, csv, gzip, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from manifest_io import ROOT, sha256_file  # noqa: E402

RAW_UNIPROT = os.path.join(ROOT, "data", "raw", "uniprot_2026_03")
ENSEMBL_PEP = os.path.join(ROOT, "data", "raw", "ensembl_116", "Homo_sapiens.GRCh38.pep.all.fa.gz")
HGNC_FILE = os.path.join(ROOT, "data", "raw", "hgnc", "hgnc_complete_set_2026-10-06.tsv")
OUT_DIR = os.path.join(ROOT, "data", "interim", "proteomes")
QC_DIR = os.path.join(ROOT, "results", "data")

# File-format plumbing only (which idmapping column holds the MOD gene ID); the ID convention itself is config/gene_ids.yaml.
SPECIES = {
    "human": {"mode": "hgnc_ensembl"},
    "yeast": {"id_types": ["Gene_OrderedLocusName"]},
    "arabidopsis": {"id_types": ["TAIR", "Araport", "Gene_OrderedLocusName"], "upper": True, "require_pattern": True},
    "pombe": {"id_types": ["PomBase"]},
    "worm": {"id_types": ["WormBase"], "orfname_xref": True},
    "fly": {"id_types": ["FlyBase"], "orfname_xref": True},
}
PRIMARY_SCAFFOLD = re.compile(r"^(GL|KI)\d+\.\d+$")


def load_gene_id_config(path=os.path.join(ROOT, "config", "gene_ids.yaml")):
    import yaml
    with open(path) as f:
        return yaml.safe_load(f)


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_uniprot_fasta(path):
    """-> list of (accession, 'sp'|'tr', sequence) in file order."""
    out, acc, db, buf = [], None, None, []
    with opener(path) as f:
        for line in f:
            if line.startswith(">"):
                if acc is not None:
                    out.append((acc, db, "".join(buf)))
                parts = line[1:].split("|")
                db, acc, buf = parts[0], parts[1], []
            else:
                buf.append(line.strip())
    if acc is not None:
        out.append((acc, db, "".join(buf)))
    return out


def read_idmapping(path, id_types, accessions=None):
    """-> {acc: {id_type: [unique ids in file order]}} for the given ID types (and accessions, if given)."""
    want = set(id_types)
    m = collections.defaultdict(lambda: collections.defaultdict(list))
    with opener(path) as f:
        for line in f:
            a, t, v = line.rstrip("\n").split("\t")
            if t in want and (accessions is None or a in accessions) and v not in m[a][t]:
                m[a][t].append(v)
    return m


def read_ensembl_pep_genes(path):
    """-> (primary-assembly ENSG set, patch/haplotype-only ENSG set) from an Ensembl pep.all FASTA."""
    prim, other = set(), set()
    with opener(path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            toks = line[1:].split()
            loc = toks[2].split(":")
            gene = next(t for t in toks if t.startswith("gene:"))[5:].split(".")[0]
            (prim if loc[0] == "chromosome" or PRIMARY_SCAFFOLD.match(loc[2]) else other).add(gene)
    return prim, other - prim


def read_hgnc_ensembl(path):
    """-> {HGNC:id: ENSG or ''}"""
    with opener(path) as f:
        return {r["hgnc_id"]: r["ensembl_gene_id"].strip() for r in csv.DictReader(f, delimiter="\t")}


def mod_gene_ids(acc, idmap, cfg, pattern=None):
    """Non-human: -> (sorted gene IDs, id_source)."""
    for t in cfg["id_types"]:
        ids = list(idmap.get(acc, {}).get(t, []))
        if cfg.get("upper"):
            ids = list(dict.fromkeys(i.upper() for i in ids))
        if cfg.get("require_pattern") and pattern:
            ids = [i for i in ids if re.match(pattern, i)]
        if ids:
            return sorted(set(ids)), t
    return [], "accession_fallback"


def human_gene_ids(acc, idmap, h2e, primary, nonprimary):
    """Human: -> (sorted ENSG IDs, id_source, note)."""
    hgnc = idmap.get(acc, {}).get("HGNC", [])
    uens = sorted(set(e.split(".")[0] for e in idmap.get(acc, {}).get("Ensembl", [])))
    uprim = [e for e in uens if e in primary]
    if hgnc:
        genes = sorted(set(h2e[h] for h in hgnc if h2e.get(h)))
        n_unmapped = sum(1 for h in hgnc if not h2e.get(h))
        if not n_unmapped:
            return genes, "HGNC>Ensembl", ""
        leftover = [e for e in uprim if e not in genes][:n_unmapped]
        genes = sorted(set(genes) | set(leftover))
        note = f"hgnc_without_ensembl_{n_unmapped}_filled_{len(leftover)}_from_uniprot_primary_ensg"
        return genes, ("HGNC>Ensembl+UniProt_Ensembl" if leftover else "HGNC>Ensembl") if genes else "accession_fallback", note
    if uens:
        pool = uprim or [e for e in uens if e not in nonprimary] or uens
        note = "no_hgnc_" + ("primary" if uprim else "nonprimary") + (f"_ensg_kept_1_of_{len(uens)}" if len(uens) > 1 else "_ensg")
        return [pool[0]], "UniProt_Ensembl", note
    return [], "accession_fallback", "no_hgnc_no_ensembl"


def build(species, fasta_path, idmapping_path, cfg, pattern=None, h2e=None, primary=None, nonprimary=None):
    """Core logic (pure, testable). -> (records [(gene_id, seq)], mapping rows [dict], qc dict)."""
    entries = read_uniprot_fasta(fasta_path)
    accs = [a for a, _, _ in entries]
    assert len(accs) == len(set(accs)), "duplicate accessions in FASTA"
    info = {a: (db, s) for a, db, s in entries}
    human = cfg.get("mode") == "hgnc_ensembl"
    types = ["HGNC", "Ensembl"] if human else list(cfg["id_types"])
    idmap = read_idmapping(idmapping_path, types + (["Gene_ORFName"] if cfg.get("orfname_xref") else []),
                           None if cfg.get("orfname_xref") else set(accs))

    acc_genes, acc_src, notes = {}, {}, collections.Counter()
    for a in accs:
        if human:
            ids, src, note = human_gene_ids(a, idmap, h2e, primary, nonprimary)
        else:
            (ids, src), note = mod_gene_ids(a, idmap, cfg, pattern), ""
        acc_genes[a], acc_src[a] = ids, [src, note]
        if note:
            notes[re.sub(r"_\d+", "_N", note)] += 1

    if cfg.get("orfname_xref"):  # MOD ID via Gene_ORFName shared with other accessions of the same idmapping file
        used = {g for a in accs for g in acc_genes[a]}
        key2mod = collections.defaultdict(set)
        for a, d in idmap.items():
            mods = set(d.get(cfg["id_types"][0], []))
            for k in d.get("Gene_ORFName", []):
                key2mod[k] |= mods
        for a in accs:
            if acc_genes[a]:
                continue
            cand = set().union(*[key2mod.get(k, set()) for k in idmap.get(a, {}).get("Gene_ORFName", [])] or [set()])
            if len(cand) == 1 and not (cand & used):
                acc_genes[a], acc_src[a] = sorted(cand), ["orfname_xref", ""]
                used |= cand
                notes["orfname_xref"] += 1

    for a in accs:
        if not acc_genes[a]:
            acc_genes[a], acc_src[a][0] = [a], "accession_fallback"

    gene_accs = collections.defaultdict(list)
    for a in accs:
        for g in acc_genes[a]:
            gene_accs[g].append(a)

    def rank(a):  # smaller is better
        db, s = info[a]
        return (len(acc_genes[a]) != 1, db != "sp", -len(s), a)

    chosen, multi_acc_genes = {}, []
    for g, al in gene_accs.items():
        best = min(al, key=rank)
        chosen[g] = best
        if len(al) > 1:
            multi_acc_genes.append({"gene": g, "kept": best, "others": sorted(x for x in al if x != best)})

    acc_final = collections.defaultdict(list)
    for g, a in chosen.items():
        acc_final[a].append(g)
    dropped = [a for a in accs if a not in acc_final]

    records, rows = [], []
    for g in sorted(chosen):
        a = chosen[g]
        db, s = info[a]
        records.append((g, s))
        rows.append({"gene_id": g, "accession": a, "db": db, "seq_len": len(s), "id_source": acc_src[a][0],
                     "n_genes_on_accession": len(acc_final[a]), "expanded": int(len(acc_final[a]) > 1),
                     "matches_id_pattern": int(bool(pattern and re.match(pattern, g))), "note": acc_src[a][1]})

    seq_count = collections.Counter(s for _, s in records)
    expanded = sorted(a for a, gl in acc_final.items() if len(gl) > 1)
    qc = {
        "species": species,
        "n_fasta_accessions": len(accs),
        "n_records_genes": len(records),
        "n_accessions_used": len(acc_final),
        "n_accessions_expanded": len(expanded),
        "n_extra_records_from_expansion": sum(len(acc_final[a]) - 1 for a in expanded),
        "n_accessions_dropped": len(dropped),
        "dropped_by_reason": {"gene_assigned_to_other_accession": len(dropped)} if dropped else {},
        "n_records_accession_fallback_id": sum(1 for r in rows if r["id_source"] == "accession_fallback"),
        "n_records_not_matching_id_pattern": sum(1 for r in rows if not r["matches_id_pattern"]),
        "n_genes_reached_from_multiple_accessions": len(multi_acc_genes),
        "mapping_notes": dict(notes),
        "id_source_counts": dict(collections.Counter(r["id_source"] for r in rows)),
        "n_records_with_identical_sequence_to_another": sum(c for c in seq_count.values() if c > 1),
        "n_identical_sequence_groups": sum(1 for c in seq_count.values() if c > 1),
        "expanded_accessions": {a: acc_final[a] for a in expanded[:300]},
        "multi_accession_genes": multi_acc_genes[:300],
        "dropped_accessions": dropped[:500],
    }
    return records, rows, qc


def write_fasta(records, path, width=60):
    with open(path, "w") as f:
        for g, s in records:
            f.write(f">{g}\n")
            for i in range(0, len(s), width):
                f.write(s[i:i + width] + "\n")


def write_tsv(rows, path):
    cols = ["gene_id", "accession", "db", "seq_len", "id_source", "n_genes_on_accession", "expanded", "matches_id_pattern", "note"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="*", default=list(SPECIES))
    args = ap.parse_args()
    gcfg = load_gene_id_config()
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(QC_DIR, exist_ok=True)
    for sp in args.species:
        cfg, gc = SPECIES[sp], gcfg[sp]
        upid = gc["uniprot_proteome"]
        fa = next(os.path.join(RAW_UNIPROT, x) for x in sorted(os.listdir(RAW_UNIPROT)) if x.startswith(upid + "_") and x.endswith(".fasta.gz"))
        im = fa.replace(".fasta.gz", ".idmapping.gz")
        inputs = [fa, im]
        kw = {}
        if cfg.get("mode") == "hgnc_ensembl":
            primary, nonprimary = read_ensembl_pep_genes(ENSEMBL_PEP)
            kw = {"h2e": read_hgnc_ensembl(HGNC_FILE), "primary": primary, "nonprimary": nonprimary}
            inputs += [HGNC_FILE, ENSEMBL_PEP]
        records, rows, qc = build(sp, fa, im, cfg, gc.get("example_pattern"), **kw)
        out_fa = os.path.join(OUT_DIR, f"{sp}.fasta")
        out_tsv = os.path.join(OUT_DIR, f"{sp}_gene2acc.tsv")
        write_fasta(records, out_fa)
        write_tsv(rows, out_tsv)
        qc.update({"uniprot_release": gcfg.get("release"), "uniprot_proteome": upid, "id_type": gc["id_type"],
                   "id_pattern": gc.get("example_pattern"),
                   "inputs_sha256": {os.path.relpath(p, ROOT): sha256_file(p) for p in inputs},
                   "outputs_sha256": {os.path.relpath(out_fa, ROOT): sha256_file(out_fa), os.path.relpath(out_tsv, ROOT): sha256_file(out_tsv)},
                   "script": "scripts/data/build_gene_proteome.py"})
        if kw:
            g = {r["gene_id"] for r in rows}
            qc["n_records_ensg_on_ensembl116_primary_with_peptide"] = len(g & kw["primary"])
            qc["n_records_ensg_patch_or_haplotype_only"] = len(g & kw["nonprimary"])
            qc["n_records_ensg_without_ensembl116_peptide"] = len({x for x in g if x.startswith("ENSG")} - kw["primary"] - kw["nonprimary"])
        with open(os.path.join(QC_DIR, f"proteome_{sp}.json"), "w") as f:
            json.dump(qc, f, indent=1)
        print(f"{sp:12s} accessions={qc['n_fasta_accessions']} genes={qc['n_records_genes']} expanded_acc={qc['n_accessions_expanded']} "
              f"(+{qc['n_extra_records_from_expansion']}) dropped={qc['n_accessions_dropped']} fallback={qc['n_records_accession_fallback_id']} "
              f"multi_acc_genes={qc['n_genes_reached_from_multiple_accessions']} notes={qc['mapping_notes']}")


if __name__ == "__main__":
    main()
