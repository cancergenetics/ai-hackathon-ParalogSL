"""Known-answer tests for the panel gap-filling rule (scripts/data/build_panel.py) and the Compara paralog universe parser
(scripts/data/build_paralog_universe.py) on tiny hand-made fixtures."""
import gzip, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "data"))
import build_panel as bp  # noqa: E402
import build_paralog_universe as bu  # noqa: E402


def _rec(upid, taxid, genus, family, order, S, D, F, n=10000, extra=()):
    lineage = [("Eukaryota", "domain"), *[(x, "clade") for x in extra], (order, "order"), (family, "family"), (genus, "genus")]
    tot = 1000
    return {"id": upid, "taxonomy": {"scientificName": f"{genus} sp", "taxonId": taxid, "mnemonic": upid[-5:]}, "proteinCount": n,
            "proteomeCompletenessReport": {"buscoReport": {"total": tot, "completeSingle": int(S * 10), "completeDuplicated": int(D * 10),
                                                          "fragmented": int(F * 10)}},
            "taxonLineage": [{"scientificName": s, "rank": r} for s, r in lineage]}


CFG = {"filters": {"max_protein_count": 50000, "max_busco_duplicated_pct": 10, "ladder": [[85, 8], [80, 8], [75, 10]]},
       "clades": [{"name": "toy", "lineage_any": ["Toyida"], "unit": "family", "quota": 2, "exclude_units": ["Banned"]}]}


def test_select_rule_known_answer():
    recs = [_rec("UP01", 1, "G1", "Core", "O", 99, 1, 0, extra=["Toyida"]),        # core proteome -> its family is used
            _rec("UP02", 2, "G2", "Core", "O", 98, 1, 0, extra=["Toyida"]),        # same family as core -> excluded
            _rec("UP03", 3, "G3", "FamA", "O", 95, 1, 1, extra=["Toyida"]),        # best of FamA
            _rec("UP04", 4, "G4", "FamA", "O", 94, 1, 1, extra=["Toyida"]),        # second in FamA -> one per unit
            _rec("UP05", 5, "G5", "FamB", "O", 97, 20, 1, extra=["Toyida"]),       # D > 10 -> out
            _rec("UP06", 6, "G6", "FamC", "O", 96, 1, 1, n=60000, extra=["Toyida"]),  # too many proteins -> out
            _rec("UP07", 7, "G7", "Banned", "O", 99, 1, 1, extra=["Toyida"]),      # excluded unit
            _rec("UP08", 8, "G8", "FamD", "O", 90, 1, 1, extra=["Toyida"]),        # second pick
            _rec("UP09", 9, "G9", "FamE", "O", 99, 1, 1)]                          # not in clade
    allp = {r["id"]: bp.parse_record(r) for r in recs}
    sel, log = bp.select_rule(allp, {"UP01"}, {"1"}, CFG)
    assert [x["id"] for _, _, x in sel] == ["UP03", "UP08"]
    assert log[0]["ladder_step_used"] == [85, 8] and log[0]["quota_met"]


def test_select_rule_ladder_relaxes_only_when_quota_unmet():
    recs = [_rec("UP11", 11, "G1", "FamA", "O", 90, 1, 1, extra=["Toyida"]),
            _rec("UP12", 12, "G2", "FamB", "O", 78, 1, 9, extra=["Toyida"])]   # passes only the last step (S>=75, F<=10)
    allp = {r["id"]: bp.parse_record(r) for r in recs}
    sel, log = bp.select_rule(allp, set(), set(), CFG)
    assert sorted(x["id"] for _, _, x in sel) == ["UP11", "UP12"]
    assert log[0]["ladder_step_used"] == [75, 10]


def test_compara_universe_parser(tmp_path):
    p = str(tmp_path / "c.tsv.gz")
    hdr = ["gene_stable_id", "protein_stable_id", "species", "identity", "homology_type", "homology_gene_stable_id",
           "homology_protein_stable_id", "homology_species", "homology_identity", "dn", "ds", "goc_score", "wga_coverage",
           "is_high_confidence", "homology_id"]
    rows = [["YB", "p", "saccharomyces_cerevisiae", "40", "within_species_paralog", "YA", "q", "saccharomyces_cerevisiae", "50", "NULL", "NULL", "NULL", "NULL", "NULL", "1"],
            ["YA", "p", "saccharomyces_cerevisiae", "50", "within_species_paralog", "YB", "q", "saccharomyces_cerevisiae", "40", "NULL", "NULL", "NULL", "NULL", "NULL", "2"],  # duplicate pair
            ["YA", "p", "saccharomyces_cerevisiae", "20", "other_paralog", "YC", "q", "saccharomyces_cerevisiae", "25", "NULL", "NULL", "NULL", "NULL", "NULL", "3"],
            ["YA", "p", "saccharomyces_cerevisiae", "90", "gene_split", "YD", "q", "saccharomyces_cerevisiae", "90", "NULL", "NULL", "NULL", "NULL", "NULL", "4"],
            ["YA", "p", "saccharomyces_cerevisiae", "60", "ortholog_one2one", "ENSG1", "q", "homo_sapiens", "60", "NULL", "NULL", "NULL", "NULL", "NULL", "5"]]
    with gzip.open(p, "wt") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    out, qc = bu.build("yeast", p, {"YA", "YB"})
    assert [(r["gene_a"], r["gene_b"], r["homology_type"]) for r in out] == [("YA", "YB", "within_species_paralog"), ("YA", "YC", "other_paralog")]
    assert out[0]["ens_identity_a"] == "50" and out[0]["ens_identity_b"] == "40"   # identities follow the swapped order
    assert qc["n_pairs"] == 2 and qc["n_genes"] == 3 and qc["n_genes_in_proteome"] == 2
    assert qc["parse_counts"]["duplicate_pair_rows"] == 1 and qc["parse_counts"]["within_species_rows_gene_split"] == 1
