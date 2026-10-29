"""Known-answer tests for scripts/data/build_gene_proteome.py on tiny hand-made fixtures (no real data needed)."""
import gzip, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "data"))
import build_gene_proteome as bgp  # noqa: E402


def _fasta(path, entries):
    with gzip.open(path, "wt") as f:
        for db, acc, seq in entries:
            f.write(f">{db}|{acc}|{acc}_TEST Some protein OS=Test OX=1 GN=x PE=1 SV=1\n{seq}\n")


def _idmap(path, rows):
    with gzip.open(path, "wt") as f:
        for r in rows:
            f.write("\t".join(r) + "\n")


def test_yeast_like_expansion_fallback_and_duplicates(tmp_path):
    fa, im = str(tmp_path / "y.fasta.gz"), str(tmp_path / "y.idmapping.gz")
    _fasta(fa, [("sp", "P00001", "MAAAA"),      # two ORFs collapsed by UniProt -> 2 records
                ("tr", "P00002", "MCCCC"),      # one ORF
                ("tr", "P00003", "MDDDD"),      # no gene ID -> accession fallback
                ("tr", "P00004", "MEEEEEEEE"),  # same ORF as P00005, longer but TrEMBL
                ("sp", "P00005", "MEEEE")])     # Swiss-Prot wins
    _idmap(im, [("P00001", "Gene_OrderedLocusName", "YAL001W"), ("P00001", "Gene_OrderedLocusName", "YAL002C"),
                ("P00001", "SGD", "S000000001"),
                ("P00002", "Gene_OrderedLocusName", "YBR001W"),
                ("P00004", "Gene_OrderedLocusName", "YCR001W"), ("P00005", "Gene_OrderedLocusName", "YCR001W")])
    recs, rows, qc = bgp.build("yeast", fa, im, bgp.SPECIES["yeast"], pattern=r"^Y[A-P][LR]\d{3}[WC](-[A-Z])?$|^Q\d{4}$")
    d = dict(recs)
    assert set(d) == {"YAL001W", "YAL002C", "YBR001W", "P00003", "YCR001W"}
    assert d["YAL001W"] == d["YAL002C"] == "MAAAA"
    assert d["YCR001W"] == "MEEEE"
    assert qc["n_accessions_expanded"] == 1 and qc["n_extra_records_from_expansion"] == 1
    assert qc["n_accessions_dropped"] == 1 and qc["dropped_accessions"] == ["P00004"]
    assert qc["n_records_accession_fallback_id"] == 1
    assert qc["n_genes_reached_from_multiple_accessions"] == 1
    assert qc["n_records_not_matching_id_pattern"] == 1  # the accession fallback
    m = {r["gene_id"]: r for r in rows}
    assert m["YAL002C"]["accession"] == "P00001" and m["YAL002C"]["expanded"] == 1


def test_human_hgnc_driven_mapping(tmp_path):
    fa, im = str(tmp_path / "h.fasta.gz"), str(tmp_path / "h.idmapping.gz")
    _fasta(fa, [("sp", "H00001", "MHIST"), ("sp", "H00002", "MCDC"), ("tr", "H00003", "MNOHG"),
                ("tr", "H00004", "MNONE"), ("sp", "H00005", "MNOENS")])
    _idmap(im, [
        # 2 HGNC genes with identical protein; UniProt also lists a patch copy ENSG00000000009 -> exactly 2 records
        ("H00001", "HGNC", "HGNC:1"), ("H00001", "HGNC", "HGNC:2"),
        ("H00001", "Ensembl", "ENSG00000000001.3"), ("H00001", "Ensembl", "ENSG00000000002.1"), ("H00001", "Ensembl", "ENSG00000000009.1"),
        # 1 HGNC gene; UniProt also lists a primary-assembly read-through model ENSG00000000004 -> 1 record (ENSG3)
        ("H00002", "HGNC", "HGNC:3"), ("H00002", "Ensembl", "ENSG00000000003.2"), ("H00002", "Ensembl", "ENSG00000000004.1"),
        # no HGNC: one patch + one primary ENSG -> the primary one
        ("H00003", "Ensembl", "ENSG00000000007.1"), ("H00003", "Ensembl", "ENSG00000000005.1"),
        # HGNC without Ensembl ID in HGNC -> filled from UniProt primary ENSG
        ("H00005", "HGNC", "HGNC:6"), ("H00005", "Ensembl", "ENSG00000000006.1")])
    h2e = {"HGNC:1": "ENSG00000000001", "HGNC:2": "ENSG00000000002", "HGNC:3": "ENSG00000000003", "HGNC:6": ""}
    primary = {"ENSG00000000001", "ENSG00000000002", "ENSG00000000003", "ENSG00000000004", "ENSG00000000005", "ENSG00000000006"}
    nonprimary = {"ENSG00000000007", "ENSG00000000009"}
    recs, rows, qc = bgp.build("human", fa, im, bgp.SPECIES["human"], pattern=r"^ENSG\d{11}$", h2e=h2e, primary=primary, nonprimary=nonprimary)
    d = dict(recs)
    assert set(d) == {"ENSG00000000001", "ENSG00000000002", "ENSG00000000003", "ENSG00000000005", "ENSG00000000006", "H00004"}
    assert d["ENSG00000000001"] == d["ENSG00000000002"] == "MHIST"
    assert qc["n_accessions_expanded"] == 1 and qc["n_extra_records_from_expansion"] == 1
    assert qc["n_records_accession_fallback_id"] == 1 and qc["n_accessions_dropped"] == 0


def test_orfname_xref_fills_missing_mod_id(tmp_path):
    fa, im = str(tmp_path / "w.fasta.gz"), str(tmp_path / "w.idmapping.gz")
    _fasta(fa, [("tr", "W00001", "MWORM"), ("tr", "W00002", "MWORMB"), ("tr", "W00003", "MWORMC")])
    _idmap(im, [("W00001", "Gene_ORFName", "X.1"),
                ("W00009", "WormBase", "WBGene00000009"), ("W00009", "Gene_ORFName", "X.1"),  # not in FASTA
                ("W00002", "WormBase", "WBGene00000002"),
                ("W00003", "Gene_ORFName", "Y.1")])                                           # no xref anywhere
    recs, rows, qc = bgp.build("worm", fa, im, bgp.SPECIES["worm"], pattern=r"^WBGene\d{8}$")
    assert dict(recs).keys() == {"WBGene00000009", "WBGene00000002", "W00003"}
    assert qc["mapping_notes"] == {"orfname_xref": 1}


def test_ensembl_pep_primary_vs_patch(tmp_path):
    p = str(tmp_path / "pep.fa.gz")
    with gzip.open(p, "wt") as f:
        f.write(">ENSP1.1 pep chromosome:GRCh38:1:1:10:1 gene:ENSG00000000001.2 transcript:ENST1.1\nMA\n"
                ">ENSP2.1 pep scaffold:GRCh38:KI270711.1:1:10:1 gene:ENSG00000000002.1 transcript:ENST2.1\nMA\n"
                ">ENSP3.1 pep scaffold:GRCh38:HSCHR6_MHC_COX_CTG1:1:10:1 gene:ENSG00000000003.1 transcript:ENST3.1\nMA\n"
                ">ENSP4.1 pep scaffold:GRCh38:HG1_PATCH:1:10:1 gene:ENSG00000000001.2 transcript:ENST4.1\nMA\n")
    prim, nonprim = bgp.read_ensembl_pep_genes(p)
    assert prim == {"ENSG00000000001", "ENSG00000000002"}
    assert nonprim == {"ENSG00000000003"}
