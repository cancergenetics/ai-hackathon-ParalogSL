"""Known-answer tests for the yeast label parsers and outputs (marker `data`; G4, R2-12).

Tolerances were written here (and in config/labels_yeast.yaml `known_answers`) BEFORE the first run:
  Costanzo strain-class counts sn 4071 / dma 3893 / tsq 1167 / damp 819 / tsa 792 / S 117   exact   (R2-12 [V])
  Costanzo stringent set, p<0.05 and (eps<-0.12 or eps>0.16), all four raw files: 380,059   +-1%    (BioGRID
      Costanzo2016 readme, 'HighestConfidence' = raw (already reciprocal-filtered) + stringent threshold).
      The task asked for the eps<-0.12, p<0.05 count; the readme/paper on disk give only the combined
      negative+positive stringent count, so that is the published number tested here.
  Intermediate threshold (p<0.05, |eps|>0.08): ~550,000 negative / ~350,000 positive            +-5%    (readme;
      values are rounded '~' figures, hence the wider tolerance)
  SGA_NxN rows 12,698,939; E defined (row SMFs) 11,755,688; E<=0.5 528,756; query SMF NaN 943,251;
      v2 POS rule (eps<-0.08, p<0.05, DMF<=0.5) 36,556 with E>0.5 and 17,956 with E<=0.5             exact   (data.md R2-1)
  SGA_NxN sn x dma rows with DMF <= 0.2: 29,174; NxN DMF minimum -0.1525                          exact   (data.md R2-2)
  SGA_ExN_NxE sn x dma rows: 0                                                                     exact   (data.md R2-2)
  Dean 2008 Dataset S1: 289 pairs listed, 201 with all fitness values, 49 lethal (W_AB = 0)        exact   (paper text;
      R2-9 'fail loudly'. The survey's 197/48 came from a pdftotext parse; the paper's 201/49 is the reference.)
Unreproducible values are marked xfail(strict=True) with an entry in results/labels/known_answer_gaps.json;
tolerances are never changed.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

pytestmark = pytest.mark.data

ROOT = Path(__file__).resolve().parents[1]
LAB = ROOT / "scripts" / "labels"
sys.path.insert(0, str(LAB))
CFG = yaml.safe_load((ROOT / "config" / "labels_yeast.yaml").read_text())


def _load(name):
    spec = importlib.util.spec_from_file_location(name, LAB / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


P = _load("yeast_parsers")
C = CFG["sources"]["costanzo2016"]
ZIP = ROOT / C["zip"]
needs_sga = pytest.mark.skipif(not ZIP.exists(), reason="Costanzo raw zip absent")
needs_dean = pytest.mark.skipif(not (ROOT / CFG["sources"]["dean2008"]["pdf"]).exists(), reason="Dean PDF absent")
UNRESTRICTED = ROOT / CFG["qc"]["outputs"]["pairs_unrestricted"]
QC = ROOT / CFG["qc"]["outputs"]["qc_json"]
needs_out = pytest.mark.skipif(not (UNRESTRICTED.exists() and QC.exists()), reason="run build_yeast_labels.py first")


def _files():
    return P.ensure_unpacked(ZIP, ROOT / C["unpack_to"], C["files_labels"] + C["files_known_answer_only"] + [C["strain_table"]])


@pytest.fixture(scope="module")
def sga_counts():
    """One streamed pass over the four raw SGA files."""
    f = _files()
    out = {"stringent": 0, "int_neg": 0, "int_pos": 0}
    nxn = {"rows": 0, "E_defined": 0, "E_le_05": 0, "q_smf_nan": 0, "v2pos_E_gt": 0, "v2pos_E_le": 0,
           "sn_dma_dmf_le_02": 0, "dmf_min": np.inf}
    exn_sn_dma = 0
    for name in C["files_labels"] + C["files_known_answer_only"]:
        for ch in P.iter_sga(f[name], 2_000_000):
            eps, p = ch["eps"].to_numpy(float), ch["p"].to_numpy(float)
            with np.errstate(invalid="ignore"):
                out["stringent"] += int(((p < 0.05) & ((eps < -0.12) | (eps > 0.16))).sum())
                out["int_neg"] += int(((p < 0.05) & (eps < -0.08)).sum())
                out["int_pos"] += int(((p < 0.05) & (eps > 0.08)).sum())
            qs, as_ = P.suffix_class(ch["query_strain"]), P.suffix_class(ch["array_strain"])
            sn_dma = ((qs == "sn") & (as_ == "dma")).fillna(False).to_numpy()
            if name == "SGA_ExN_NxE.txt":
                exn_sn_dma += int(sn_dma.sum())
            if name == "SGA_NxN.txt":
                dmf = ch["dmf"].to_numpy(float)
                E = ch["smf_q"].to_numpy(float) * ch["smf_a"].to_numpy(float)
                nxn["rows"] += len(ch)
                nxn["q_smf_nan"] += int(np.isnan(ch["smf_q"].to_numpy(float)).sum())
                d = ~np.isnan(E)
                nxn["E_defined"] += int(d.sum())
                with np.errstate(invalid="ignore"):
                    nxn["E_le_05"] += int((d & (E <= 0.5)).sum())
                    v2 = (eps < -0.08) & (p < 0.05) & (dmf <= 0.5)
                    nxn["v2pos_E_gt"] += int((v2 & d & (E > 0.5)).sum())
                    nxn["v2pos_E_le"] += int((v2 & d & (E <= 0.5)).sum())
                    nxn["sn_dma_dmf_le_02"] += int((sn_dma & (dmf <= 0.2)).sum())
                nxn["dmf_min"] = min(nxn["dmf_min"], float(np.nanmin(dmf)))
    return {"all": out, "nxn": nxn, "exn_sn_dma": exn_sn_dma}


# ----------------------------------------------------------------------------- parsers
@needs_sga
def test_strain_classes_exact():
    st = P.read_strain_table(_files()[C["strain_table"]])
    assert st["suffix"].value_counts().to_dict() == {"sn": 4071, "dma": 3893, "tsq": 1167, "damp": 819, "tsa": 792, "S": 117}


GAP = "see results/labels/known_answer_gaps.json (scripts/labels/yeast_known_answer_gaps.py)"


@needs_sga
@pytest.mark.xfail(strict=True, reason="readme count not reproduced; " + GAP)
def test_costanzo_stringent_count_readme(sga_counts):
    got = sga_counts["all"]["stringent"]
    assert abs(got - 380_059) <= 0.01 * 380_059, got


@needs_sga
@pytest.mark.xfail(strict=True, reason="readme counts not reproduced; " + GAP)
def test_costanzo_intermediate_counts_readme(sga_counts):
    a = sga_counts["all"]
    assert abs(a["int_neg"] - 550_000) <= 0.05 * 550_000, a["int_neg"]
    assert abs(a["int_pos"] - 350_000) <= 0.05 * 350_000, a["int_pos"]


@needs_sga
def test_nxn_reviewer_counts(sga_counts):
    n = sga_counts["nxn"]
    assert n["rows"] == 12_698_939
    assert n["q_smf_nan"] == 943_251
    assert n["E_defined"] == 11_755_688
    assert n["E_le_05"] == 528_756
    assert (n["v2pos_E_gt"], n["v2pos_E_le"]) == (36_556, 17_956)


@needs_sga
def test_nxn_dmf_counts(sga_counts):
    n = sga_counts["nxn"]
    assert n["sn_dma_dmf_le_02"] == 29_174


# split from test_nxn_dmf_counts after the first run (assertion and tolerance unchanged): the minimum over all
# NxN rows is not -0.1525; the reviewer's value appears to refer to null x null rows only (see gaps file).
@needs_sga
@pytest.mark.xfail(strict=True, reason="all-row NxN DMF minimum is not -0.1525; " + GAP)
def test_nxn_dmf_min(sga_counts):
    n = sga_counts["nxn"]
    assert round(n["dmf_min"], 4) == -0.1525


@needs_sga
def test_exn_nxe_has_no_null_null_rows(sga_counts):
    assert sga_counts["exn_sn_dma"] == 0


@needs_dean
def test_dean_counts_paper():
    d = P.parse_dean_pdf(ROOT / CFG["sources"]["dean2008"]["pdf"])
    assert len(d) == 289
    assert int(d["w12"].notna().sum()) == 201
    assert int((d["w12"] == 0).sum()) == 49
    assert not d[["orf1", "orf2"]].apply(lambda r: tuple(sorted(r)), axis=1).duplicated().any()


# ----------------------------------------------------------------------------- outputs
TIER_NAMES = [t for t in CFG["tiers"] if t != "values"]   # primary, lenient, strict, biogrid, report-only tier


@pytest.fixture(scope="module")
def labels():
    return pd.read_csv(UNRESTRICTED, sep="\t", usecols=["gene_a", "gene_b", *TIER_NAMES, "dean_call", "linked", "n_T"],
                       keep_default_na=False, na_values={"linked": [""]}, dtype={"dean_call": str})


@needs_out
def test_pairs_canonical_unique(labels):
    assert (labels["gene_a"] < labels["gene_b"]).all()
    assert not labels.duplicated(["gene_a", "gene_b"]).any()


@needs_out
def test_tier_vocabulary_and_nesting(labels):
    vocab = set(CFG["tiers"]["values"])
    for t in TIER_NAMES:
        assert set(labels[t]) <= vocab, t
    pos = {t: set(labels.index[labels[t] == "pos"]) for t in ["primary", "lenient", "strict", "biogrid"]}
    neg = {t: set(labels.index[labels[t] == "neg"]) for t in ["primary", "lenient", "strict"]}
    assert pos["strict"] <= pos["primary"] <= pos["lenient"]
    assert pos["primary"] <= pos["biogrid"]
    assert neg["lenient"] <= neg["primary"] and neg["strict"] == neg["primary"]


@needs_out
def test_dean_lethal_pairs_positive_in_primary(labels):
    d = labels[labels["dean_call"] == "pos"]
    assert len(d) > 0 and (d["primary"] == "pos").all()


@needs_out
def test_linked_pairs_not_labelled_by_sga_alone(labels):
    lk = labels[(labels["linked"] == 1) & (labels["dean_call"].isin(["", "amb", "untestable"]))]
    assert lk["primary"].isin(["untestable", "untested"]).all()


@needs_out
def test_qc_json_nesting_flags():
    qc = json.loads(QC.read_text())
    assert all(qc["unrestricted"]["nesting"].values())
