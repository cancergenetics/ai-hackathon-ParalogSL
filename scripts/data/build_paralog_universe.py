#!/usr/bin/env python
"""Label universe (DESIGN v3 §2.3, R2-4): within-species paralog pairs from ONE pinned Ensembl Compara release (116).

Run:  .venv/bin/python scripts/data/build_paralog_universe.py [--species human yeast]

Input   data/raw/ensembl_116/Compara.116.protein_default.homologies.<ensembl species>.tsv.gz  (config/data_manifest.yaml)
        data/interim/proteomes/<species>_gene2acc.tsv   (our per-gene proteome, scripts/data/build_gene_proteome.py)
Output  data/interim/universe/<species>_paralog_pairs.tsv  unordered gene pairs (gene_a < gene_b), our gene-ID convention
        results/data/universe_<species>.json                QC: n pairs, n genes, by homology type, coverage by our proteome

Rows kept: species == homology_species == the target species and homology_type in PARALOG_TYPES (every within-species
paralogue type Compara emits; 'gene_split' rows are not paralogues and are counted, not kept). Gene IDs are Compara gene
stable IDs, which already follow config/gene_ids.yaml (human ENSG without version, yeast ordered locus name).
The Ensembl identity / homology-type columns are kept ONLY for evaluation (stratification, caller recall); they must never
reach features (Ensembl Compara is an annotation, DESIGN §1 D1).
"""
import argparse, collections, csv, gzip, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from manifest_io import ROOT, sha256_file  # noqa: E402

RELEASE = "116"
SPECIES = {"human": "homo_sapiens", "yeast": "saccharomyces_cerevisiae"}
PARALOG_TYPES = ("within_species_paralog", "other_paralog")
RAW = os.path.join(ROOT, "data", "raw", f"ensembl_{RELEASE}")
OUT_DIR = os.path.join(ROOT, "data", "interim", "universe")
QC_DIR = os.path.join(ROOT, "results", "data")
COLS = ["gene_a", "gene_b", "homology_type", "ens_identity_a", "ens_identity_b", "homology_id", "a_in_proteome", "b_in_proteome"]


def parse_compara(path, ens_species):
    """-> (pairs {(a,b): row}, counters). a < b; one row per unordered pair (first homology_id kept, duplicates counted)."""
    pairs, c = {}, collections.Counter()
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            c["rows_total"] += 1
            if r["species"] != ens_species or r["homology_species"] != ens_species:
                continue
            t = r["homology_type"]
            c[f"within_species_rows_{t}"] += 1
            if t not in PARALOG_TYPES:
                continue
            a, b, ia, ib = r["gene_stable_id"], r["homology_gene_stable_id"], r["identity"], r["homology_identity"]
            if a == b:
                c["self_pairs_skipped"] += 1
                continue
            if b < a:
                a, b, ia, ib = b, a, ib, ia
            if (a, b) in pairs:
                c["duplicate_pair_rows"] += 1
                continue
            pairs[(a, b)] = {"gene_a": a, "gene_b": b, "homology_type": t, "ens_identity_a": ia, "ens_identity_b": ib,
                             "homology_id": r["homology_id"]}
    return pairs, c


def build(species, compara_path, proteome_genes):
    pairs, c = parse_compara(compara_path, SPECIES[species])
    rows = []
    for (a, b) in sorted(pairs):
        r = dict(pairs[(a, b)])
        r["a_in_proteome"], r["b_in_proteome"] = int(a in proteome_genes), int(b in proteome_genes)
        rows.append(r)
    genes = {g for ab in pairs for g in ab}
    qc = {"species": species, "ensembl_release": RELEASE, "paralog_types_kept": list(PARALOG_TYPES),
          "n_pairs": len(rows), "n_genes": len(genes),
          "n_pairs_by_type": dict(collections.Counter(r["homology_type"] for r in rows)),
          "parse_counts": dict(c),
          "n_genes_in_proteome": len(genes & proteome_genes),
          "frac_genes_in_proteome": round(len(genes & proteome_genes) / len(genes), 6) if genes else None,
          "n_pairs_both_in_proteome": sum(1 for r in rows if r["a_in_proteome"] and r["b_in_proteome"]),
          "frac_pairs_both_in_proteome": round(sum(1 for r in rows if r["a_in_proteome"] and r["b_in_proteome"]) / len(rows), 6) if rows else None,
          "n_proteome_genes": len(proteome_genes),
          "frac_proteome_genes_in_universe": round(len(genes & proteome_genes) / len(proteome_genes), 6) if proteome_genes else None,
          "genes_not_in_proteome_examples": sorted(genes - proteome_genes)[:30]}
    return rows, qc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="*", default=list(SPECIES))
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    for sp in args.species:
        src = os.path.join(RAW, f"Compara.{RELEASE}.protein_default.homologies.{SPECIES[sp]}.tsv.gz")
        g2a = os.path.join(ROOT, "data", "interim", "proteomes", f"{sp}_gene2acc.tsv")
        with open(g2a) as f:
            proteome_genes = {r["gene_id"] for r in csv.DictReader(f, delimiter="\t")}
        rows, qc = build(sp, src, proteome_genes)
        out = os.path.join(OUT_DIR, f"{sp}_paralog_pairs.tsv")
        with open(out, "w") as f:
            f.write("\t".join(COLS) + "\n")
            for r in rows:
                f.write("\t".join(str(r[k]) for k in COLS) + "\n")
        qc.update({"inputs_sha256": {os.path.relpath(src, ROOT): sha256_file(src), os.path.relpath(g2a, ROOT): sha256_file(g2a)},
                   "source_url": f"https://ftp.ensembl.org/pub/release-{RELEASE}/tsv/ensembl-compara/homologies/{SPECIES[sp]}/Compara.{RELEASE}.protein_default.homologies.tsv.gz",
                   "output_sha256": {os.path.relpath(out, ROOT): sha256_file(out)},
                   "script": "scripts/data/build_paralog_universe.py",
                   "note": "Ensembl identity/homology_type columns are for evaluation only, never features."})
        with open(os.path.join(QC_DIR, f"universe_{sp}.json"), "w") as f:
            json.dump(qc, f, indent=1)
        print(f"{sp}: {qc['n_pairs']} pairs, {qc['n_genes']} genes, by type {qc['n_pairs_by_type']}, "
              f"genes in proteome {qc['frac_genes_in_proteome']}, pairs both in proteome {qc['frac_pairs_both_in_proteome']}; {qc['parse_counts']}")


if __name__ == "__main__":
    main()
