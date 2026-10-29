"""Known-answer tests for the human label parsers (config/labels_human.yaml).

Expected values and tolerances were written here BEFORE the parsers / build script were first run:
  Dede 2020 Table S2, zdLFC < -3: A549 17, HT29 21, OVCAR8 19 hits; 24 unique pairs; 403 pairs     -- exact
  Thompson 2021 SD1 D28 rra&t FDR<0.1 & filter: A375 44, MeWo 40, RPE-1 57, equal to SD5 lists       -- exact
  Thompson 2021: union 105 SL pairs, 72 of class Paralogous_gene_pair (paper)                          -- exact
  Thompson 2021 SD1 classes 645/447/95/4/1; filter 'No' 549/431/388                                    -- exact
  Harle 2025 S5: 117 pairs SL in >=1 line, 355 never, 882 pair x line hits, 4 pairs in >=25/27 lines  -- exact
  Harle 2025 S4 authors' rule (GI<-0.5, FDR<0.01, Bassik hit, no single depletion) == S5 cell-wise     -- exact
  Chan 2025 T6 authors' rule: 105 pairs SL in >=1 line; 23 pairs SL in >=6 lines                       -- exact
  CHyMErA S8 rows 692/692/691/691 (HAP1.T12, HAP1.T18, RPE1.T18, RPE1.T24); duplicated rows
      (rows - distinct unordered pairs) 4/4/3/3                                                        -- exact
  Ito 2021: 4,859 GEMINI pairs (Supp. table 10); 5,065 DKO rows in harmonised Supp_table1 Ito sheets   -- exact
  Ito 2021: GEMINI FDR < 0.05 hits per line within [325, 610] (R2-5 reviewer count)                   -- range
  in4mer/Chou 2025 T1: 4,540 pairs; prototype (K562 & A549) 2,070; Inzolia (MelJuso & A375) 4,414     -- exact
  Chou ZdLFC re-implementation on the gene-level LFC file reproduces T1 ZdLFC                         -- max |diff| <= 0.01
  in4mer Supp_table1: Parrish sheets 1,030 pairs; exp_mean == gene1_mean + gene2_mean in every sheet  -- |diff| <= 0.015
A count that cannot be reproduced is marked xfail(strict=True) with a note and listed in
results/labels/known_answer_gaps.json; thresholds in config/labels_human.yaml are never changed to hit a count.
"""
from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

pytestmark = pytest.mark.data

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((ROOT / "config" / "labels_human.yaml").read_text())
sys.path.insert(0, str(ROOT / "scripts" / "labels"))
hp = importlib.import_module("human_parsers")
ST1 = ROOT / "data/raw/in4mer2024/Supp_table1.xlsx"
needs = pytest.mark.skipif(not ST1.exists(), reason="raw human screen tables absent")


@pytest.fixture(scope="module")
def st1():
    return hp.read_supp_table1(ST1)


def _k(s, sep):
    a, b = str(s).split(sep)
    return hp.pair_key(a, b)


# ---------------------------------------------------------------------------------------------------------------- Dede
@needs
def test_dede_zdlfc_counts():
    raw = pd.read_csv(ROOT / CFG["studies"]["dede2020"]["significance_table"], sep="\t", index_col=0)
    assert len(raw) == 403
    hits = raw < -3
    assert hits.sum().to_dict() == {"A549": 17, "HT29": 21, "OVCAR8": 19}
    assert int(hits.any(axis=1).sum()) == 24


# ------------------------------------------------------------------------------------------------------------ Thompson
@needs
def test_thompson_authors_rule_reproduces_sd5():
    sc = CFG["studies"]["thompson2021"]
    sd1 = hp.read_thompson_sd1(ROOT, sc)
    sd5 = pd.read_excel(ROOT / "data/raw/thompson2021/41467_2021_21478_MOESM8_ESM.xlsx", sheet_name="5", header=3)
    expected = {"A375": 44, "MeWo": 40, "RPE1": 57}
    sd5_col = {"A375": "A375", "MeWo": "Mewo", "RPE1": "RPE"}
    union, union_par = set(), set()
    for lab, spec in sc["lines"].items():
        m = hp.thompson_strict(sd1, spec["sig_prefix"], spec["filter_col"])
        ours = {_k(p, "_") for p in sd1.loc[m, "GENE PAIR"]}
        theirs = {_k(p, "_") for p in sd5[sd5_col[lab]].dropna()}
        assert len(ours) == expected[lab]
        assert ours == theirs
        union |= ours
        union_par |= {_k(p, "_") for p in sd1.loc[m & (sd1["Gene pair class"] == "Paralogous_gene_pair"), "GENE PAIR"]}
    assert len(union) == 105
    assert len(union_par) == 72


@needs
def test_thompson_classes_and_filters():
    sd1 = hp.read_thompson_sd1(ROOT, CFG["studies"]["thompson2021"])
    vc = sd1["Gene pair class"].value_counts().to_dict()
    assert vc["Paralogous_gene_pair"] == 645 and vc["Bioinformatically_derived_pair"] == 447
    assert vc["SynLethDB_gene_pair"] == 95 and vc["Miscellaneous_pair"] == 4
    assert sum(v for k, v in vc.items() if str(k).startswith("Control")) == 1
    no = [int((sd1[c] == "No").sum()) for c in ("Passes A375 filter?", "Passes the Mewo filter?", "Passes RPE filter?")]
    assert no == [549, 431, 388]


# --------------------------------------------------------------------------------------------------------------- Harle
@needs
def test_harle_s5_counts():
    s5 = hp.read_harle_s5(ROOT, CFG["studies"]["harle2025"])
    lines = [c for c in s5.columns[1:28]]
    assert len(s5) == 472 and len(lines) == 27
    m = s5[lines].astype(int)
    per_pair = m.sum(axis=1)
    assert int((per_pair >= 1).sum()) == 117
    assert int((per_pair == 0).sum()) == 355
    assert int(m.values.sum()) == 882
    assert int((per_pair >= 25).sum()) == 4


@needs
def test_harle_s4_authors_rule_equals_s5():
    sc = CFG["studies"]["harle2025"]
    s4, s5 = hp.read_harle_s4(ROOT, sc), hp.read_harle_s5(ROOT, sc)
    call = hp.authors_call_harle(s4)
    ours = {(str(p), str(l)) for p, l, c in zip(s4["sorted_gene_pair"], s4["cell_line_label"], call) if c}
    lines = list(s5.columns[1:28])
    theirs = {(str(p), str(l)) for _, r in s5.iterrows() for l in lines if int(r[l]) == 1 for p in [r["sorted_gene_pair"]]}
    assert len(ours) == 882
    assert ours == theirs


# ---------------------------------------------------------------------------------------------------------------- Chan
@needs
def test_chan_authors_rule_counts():
    t6 = hp.read_chan_t6(ROOT, CFG["studies"]["chan2025"])
    assert t6["cell_line_label"].nunique() == 10
    call = hp.authors_call_harle(t6)
    per_pair = t6.assign(c=call.astype(int)).groupby("sorted_gene_pair")["c"].sum()
    assert int((per_pair >= 1).sum()) == 105
    assert int((per_pair >= 6).sum()) == 23


# ------------------------------------------------------------------------------------------------------------- CHyMErA
@needs
def test_chymera_sheet_rows_and_duplicates():
    sc = CFG["studies"]["chymera2020"]
    exp = {"HAP1.T12": (692, 4), "HAP1.T18": (692, 4), "RPE1.T18": (691, 3), "RPE1.T24": (691, 3)}
    for sheet, (n, dup) in exp.items():
        df = hp.read_chymera_sheet(ROOT, sc, sheet)
        assert len(df) == n, sheet
        assert len(df) - df["k"].nunique() == dup, sheet


# ----------------------------------------------------------------------------------------------------------------- Ito
@needs
def test_ito_table_sizes(st1):
    fdr = hp.read_ito_fdr(ROOT, CFG["studies"]["ito2021"])
    assert len(fdr) == 4859
    for spec in CFG["studies"]["ito2021"]["lines"].values():
        assert len(st1[spec["lfc_sheet"]]) == 5065


@needs
def test_ito_fdr_hits_per_line_range():
    fdr = hp.read_ito_fdr(ROOT, CFG["studies"]["ito2021"])
    for lab in CFG["studies"]["ito2021"]["lines"]:
        n = int((pd.to_numeric(fdr[lab], errors="coerce") < 0.05).sum())
        assert 325 <= n <= 610, (lab, n)


# -------------------------------------------------------------------------------------------------------- in4mer/Chou
@needs
def test_chou_t1_library_sizes():
    z = pd.read_excel(ROOT / CFG["studies"]["in4mer2024"]["significance_table"], sheet_name="ZdLFC")
    assert len(z) == 4540
    assert int(z[["K562_ZdLFC", "A549_ZdLFC"]].notna().all(axis=1).sum()) == 2070
    assert int(z[["MelJuso_ZdLFC", "A375_ZdLFC"]].notna().all(axis=1).sum()) == 4414


@needs
def test_chou_zscore_reimplementation_matches_t1():
    sc = CFG["studies"]["in4mer2024"]
    z = pd.read_excel(ROOT / sc["significance_table"], sheet_name="ZdLFC").set_index("Unnamed: 0")
    fc = hp.read_in4mer_lfc(ROOT, sc)
    pairs = [p for p in fc.index if len(str(p).split("_")) == 2]
    libs = {("K562_T8", "A549_T21"): ("K562", "A549"), ("MelJuso_T21", "A375_T21"): ("MelJuso", "A375")}
    for cols, labs in libs.items():
        lib = fc.loc[pairs, list(cols)].dropna().index
        # as in Chou's notebook (calculate_fcOBS_fcEXP_dLFC): pairs whose single-gene rows are missing are skipped
        lib = [p for p in lib if all(g in fc.index for g in p.split("_"))]
        for col, lab in zip(cols, labs):
            obs = fc.loc[lib, col].astype(float)
            exp = np.array([fc.loc[p.split("_"), col].sum() for p in lib])
            mine = pd.Series(hp.chou_zscore(obs.values - exp), index=lib)
            ref = z.loc[lib, f"{lab}_ZdLFC"].astype(float)
            assert float((mine - ref).abs().max()) <= 0.01, lab


# ------------------------------------------------------------------------------------------------------ Supp_table1
@needs
def test_supp_table1_consistency(st1):
    assert len(st1["Parrish_PC9"]) == 1030 and len(st1["Parrish_HeLa"]) == 1030
    for name, df in st1.items():
        d = (df["exp_mean"] - df["gene1_mean"] - df["gene2_mean"]).abs().max()
        assert d <= 0.015, name


# ----------------------------------------------------------------------------------------------- built label file
LAB = ROOT / "data/labels/human_labels_all.tsv"


@pytest.mark.skipif(not LAB.exists(), reason="labels not built")
def test_label_file_invariants():
    lab = pd.read_csv(LAB, sep="\t")
    assert lab["gene_a"].str.match(r"^ENSG\d{11}$").all() and lab["gene_b"].str.match(r"^ENSG\d{11}$").all()
    assert (lab["gene_a"] < lab["gene_b"]).all()
    assert not lab.duplicated(["gene_a", "gene_b"]).any()
    banned = re.compile(r"expr|tpm|essential|depmap|ppi|copy|split|fold", re.I)
    assert not [c for c in lab.columns if banned.search(c)]
    # nesting: primary positives/negatives are lenient positives/negatives; two-study positives are lenient positives
    assert (lab.loc[lab.label_primary == 1, "label_lenient"] == 1).all()
    assert (lab.loc[lab.label_primary == 0, "label_lenient"] == 0).all()
    assert (lab.loc[lab.label_two_studies == 1, "label_lenient"] == 1).all()
    pos = lab[lab.label_primary == 1]
    assert (pos.n_sl >= 2).all() and (pos.penetrance >= 0.10).all()
    neg = lab[lab.label_primary == 0]
    assert (neg.n_testable >= 2).all() and (neg.n_sl == 0).all()
    assert (neg.n_borderline <= np.floor(0.1 * neg.n_testable)).all()


@pytest.mark.skipif(not LAB.exists(), reason="labels not built")
def test_label_file_line_ids_and_d8_columns():
    lab = pd.read_csv(LAB, sep="\t")
    ids = {x.split("=")[0] for v in lab["merged_line_calls"].dropna() for x in v.split(",")}
    assert all(re.match(r"^(CVCL_\w+|name:[A-Z0-9]+)$", i) for i in ids), sorted(i for i in ids if "CVCL" not in i)[:5]
    assert "harle_chan_dependent" in lab.columns
    assert lab["harle_chan_dependent"].isin([True, False]).all()
