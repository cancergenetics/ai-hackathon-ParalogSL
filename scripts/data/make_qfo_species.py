#!/usr/bin/env python
"""Write config/qfo_species.tsv (taxon id, UniProt proteome id, oscode, name, protein count) for the QfO panel in data/raw/qfo/.
Protein counts are re-counted from the FASTA headers (gzip); names/oscodes come from data/raw/qfo/species_list.tsv. Stdlib only."""
import csv, gzip, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
Q = os.path.join(ROOT, "data", "raw", "qfo")
rows, bad = [], []
for line in open(os.path.join(Q, "species_list.tsv")):
    upid, taxid, osc, n_listed, name = line.rstrip("\n").split("\t")
    fa = os.path.join(Q, f"{upid}_{taxid}.fasta.gz")
    with gzip.open(fa, "rt") as f:
        n = sum(1 for l in f if l.startswith(">"))
    if n != int(n_listed):
        bad.append((upid, n_listed, n))
    rows.append((int(taxid), upid, osc, name, n))
rows.sort()
out = os.path.join(ROOT, "config", "qfo_species.tsv")
with open(out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["taxon_id", "uniprot_proteome", "oscode", "name", "n_proteins"])
    w.writerows(rows)
print(f"{len(rows)} species, {sum(r[4] for r in rows)} proteins -> {out}")
if bad:
    print("WARNING count mismatches vs species_list.tsv:", bad); sys.exit(1)
