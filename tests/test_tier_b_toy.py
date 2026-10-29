"""Tier B on the planted toy panel (tests/fixtures/make_toy_panel.py).

Known answers come from the construction (substitution-only copies with disjoint mutation sets,
so identities and therefore d_s are exact), not from a previous run: self-exclusion of an
identical species, d_s, pre/post/undetermined/absent classification, depth bounds, censoring flags,
conflict resolution, relative rate, outgroup asymmetry, masking (--panel-dmin), symmetry, the
NaN-vs-0 bin semantics, order invariance of the core set, identical scoring settings in self- and
panel searches, and no proteome-constant registered Tier B feature; the row alignment / pos mapping of
scripts/features/run_tierB.py on a FASTA in non-canonical order; an end-to-end regression against
a stored table; proteome-level QC flags (metadata only, never feature columns).
"""

import copy
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from slparalog import fasta, panel, pipeline, search
from slparalog.config import CONFIG_DIR, PROJECT_ROOT, load_yaml
from slparalog.features import registry, tier_b

FIX = Path(__file__).parent / "fixtures"
sys.path.insert(0, str(FIX))
import make_toy_panel as mk  # noqa: E402

sys.path.pop(0)
CCFG = load_yaml(CONFIG_DIR / "caller.yaml")


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("toy_panel")
    query, pan, exp = mk.build()
    qf, pf, st = mk.write(tmp, query, pan)
    pc = copy.deepcopy(panel.load_panel_config())
    pc["panel"].update({"fasta": str(pf), "species_table": str(st)})
    pc["tier_b"].update({"db_dir": str(tmp / "db"), "cache_dir": str(tmp / "hits")})
    pc["tier_b"]["core"]["min_size"] = 10   # toy has 30 core genes (real panel: 200)
    pcy = tmp / "panel.yaml"
    pcy.write_text(yaml.safe_dump(pc))
    wd = tmp / "work"
    df, meta = pipeline.run(qf, tmp / "out.parquet", threads=2, features_only=True, tier_b=True,
                            panel_config_path=pcy, workdir=wd)
    seqs_by_id, _ = fasta.read_proteome(qf)
    ids = fasta.canonical_order(seqs_by_id)
    return {"df": df, "meta": meta, "exp": exp, "pc": pc, "wd": wd, "ids": ids,
            "seqs": [seqs_by_id[k] for k in ids], "tmp": tmp}


def _row(df, a, b):
    r = df[(df.gene_a == a) & (df.gene_b == b)]
    assert len(r) == 1, (a, b)
    return r.iloc[0]


def _stage(toy, dmin=None, swap=False, internal=False):
    """Re-run the Tier B stage from the cached panel hits (no new search)."""
    mm = search.MMseqs(binary=search.resolve_binary(CCFG), threads=1, workdir=toy["wd"])
    aln = search.read_alignments(toy["wd"] / "realn.tsv")
    pos = {g: k for k, g in enumerate(toy["ids"])}
    df = toy["df"]
    i = df.gene_a.map(pos).to_numpy()
    j = df.gene_b.map(pos).to_numpy()
    pairs = pd.DataFrame({"i": j, "j": i} if swap else {"i": i, "j": j})
    paths, _ = panel.search_panel(mm, CCFG, toy["pc"], "db", toy["seqs"])
    n = len(toy["seqs"])
    _, _, _, _, selfbits = tier_b.directional(aln, i[:0], j[:0], n)
    pdata = panel.load_panel_data(paths, n, selfbits, toy["pc"])
    cal = panel.calibrate(pdata, toy["pc"])
    mask = panel.apply_mask(cal, toy["pc"], dmin)
    res = tier_b.compute(pairs, aln, toy["seqs"], pdata, cal, mask, toy["pc"], return_internal=internal)
    out = res[0].assign(gene_a=df.gene_a.to_numpy(), gene_b=df.gene_b.to_numpy())
    return (out, pdata, cal, mask, res) if internal else out


def test_expected_pairs_called(toy):
    got = set(zip(toy["df"].gene_a, toy["df"].gene_b))
    assert {("E1", "E2"), ("P1", "P2"), ("Q1", "Q2"), ("X1", "X2"), ("Y1", "Y2"), ("Z1", "Z2"),
            ("F1", "F2")} <= got


def test_self_exclusion_and_d_s(toy):
    cal = toy["meta"]["tier_b"]["calibration"]
    assert cal["self_excluded"] == ["S0"]
    assert cal["per_species"]["S0"]["median_rbh_identity"] == pytest.approx(1.0)
    for s, d in toy["exp"]["d"].items():
        assert cal["per_species"][s]["d_s"] == pytest.approx(d, abs=1e-5), s
    assert cal["core_size"] >= 30
    assert toy["meta"]["tier_b"]["mask"]["n_species"] == 3


@pytest.mark.parametrize("fam", ["P", "Y", "E", "Q", "Z", "X"])
def test_planted_histories(toy, fam):
    e = toy["exp"][fam]
    r = _row(toy["df"], f"{fam}1", f"{fam}2")
    assert r.tb_depth == pytest.approx(e["depth"], abs=1e-5)
    if "younger" in e:
        assert r.tb_younger_than_nearest == e["younger"]
    if "older" in e:
        assert r.tb_older_than_panel == e["older"]
    for k, col in (("pre", "tb_frac_pre"), ("post", "tb_frac_post"), ("undet", "tb_frac_undet"),
                   ("absent", "tb_frac_absent")):
        if k in e:
            assert r[col] * 3 == pytest.approx(e[k]), (fam, k)
    if "disc" in e:
        assert r.tb_discordant_frac == pytest.approx(e["disc"])
        assert r.tb_depth_range < 0   # conflict: L >= U


def test_pre_post_bounds_and_rel_rate(toy):
    d = toy["exp"]["d"]
    p = _row(toy["df"], "P1", "P2")
    assert p.tb_depth_range == pytest.approx(d["S3"] - d["S2"], abs=1e-5)
    assert p.tb_rel_rate == pytest.approx(toy["exp"]["P_rel_rate"], abs=1e-4)
    assert p.tb_no_evidence == 0
    y = _row(toy["df"], "Y1", "Y2")
    assert y.tb_depth_range == pytest.approx(d["S1"], abs=1e-5)   # L := 0
    assert y.tb_n_post_single == 3                                  # 1:2 pattern in every species
    f = _row(toy["df"], "F1", "F2")
    assert f.tb_fragment_suspect == 1 and f.tb_no_evidence == 1 and np.isnan(f.tb_depth)


def test_outgroup_asymmetry_known_answer(toy):
    y = _row(toy["df"], "Y1", "Y2")
    vals = []
    for ia, ib in toy["exp"]["Y_outgroup"]:
        da, db = panel.kimura(ia), panel.kimura(ib)
        vals.append((abs(da - db) / (da + db), abs(da - db), max(ia, ib)))
    med = np.median(np.array(vals), axis=0)
    assert y.tb_n_outgroups == 3
    assert y.tb_asym_norm == pytest.approx(med[0], abs=1e-6)
    assert y.tb_asym_abs == pytest.approx(med[1], abs=1e-6)
    assert y.tb_outgroup_ident_max == pytest.approx(med[2], abs=1e-6)
    assert toy["meta"]["tier_b"]["features"]["n_cigar_nident_mismatch"] == 0


def test_symmetry(toy):
    a = _stage(toy)
    b = _stage(toy, swap=True)
    pd.testing.assert_frame_equal(a, b)
    base = toy["df"][registry.TIER_B_COLUMNS].reset_index(drop=True)
    pd.testing.assert_frame_equal(a[registry.TIER_B_COLUMNS].reset_index(drop=True), base)


def test_masking_recomputes_depth_and_flags(toy):
    d = toy["exp"]["d"]
    m1 = _stage(toy, dmin=0.2)   # removes S1
    y = _row(m1, "Y1", "Y2")
    assert y.tb_depth == pytest.approx(d["S2"] / 2, abs=1e-5) and y.tb_younger_than_nearest == 1
    q = _row(m1, "Q1", "Q2")     # S1 was its only pre species
    assert q.tb_younger_than_nearest == 1 and q.tb_depth == pytest.approx(d["S3"] / 2, abs=1e-5)
    assert q.tb_frac_undet == pytest.approx(0.5)
    p = _row(m1, "P1", "P2")
    assert p.tb_depth == pytest.approx((d["S2"] + d["S3"]) / 2, abs=1e-5)
    m2 = _stage(toy, dmin=0.4)   # only S3 left
    p = _row(m2, "P1", "P2")
    assert p.tb_younger_than_nearest == 1 and p.tb_depth == pytest.approx(d["S3"] / 2, abs=1e-5)
    e = _row(m2, "E1", "E2")
    assert e.tb_older_than_panel == 1 and e.tb_depth == pytest.approx(d["S3"], abs=1e-5)
    z = _row(m2, "Z1", "Z2")
    assert z.tb_no_evidence == 1 and z.tb_frac_absent == 1.0
    m3 = _stage(toy, dmin=1.0)   # nothing left: every pair without evidence, nothing crashes
    assert (m3.tb_no_evidence == 1).all() and m3.tb_depth.isna().all()


def test_bin_nan_vs_zero(toy):
    out, pdata, cal, mask, res = _stage(toy, internal=True)
    bf = res[2]["bin_fractions"]
    k = int(np.flatnonzero((out.gene_a == "Z1").to_numpy())[0])
    # S1 in B2, S2 in B3, S3 in B4; B1, B5, B6 empty -> NaN ("not observable")
    for b in (0, 4, 5):
        assert np.isnan(bf["pre"][k, b]) and np.isnan(bf["absent"][k, b])
    assert bf["pre"][k, 1] == 1.0 and bf["absent"][k, 1] == 0.0      # species present, observed 0
    assert bf["absent"][k, 2] == 1.0 and bf["pre"][k, 2] == 0.0


def test_core_set_order_invariant(toy):
    out, pdata, cal, mask, _ = _stage(toy, internal=True)
    perm = np.random.default_rng(0).permutation(pdata.best_bits.shape[0])
    pd2 = panel.PanelData(codes=pdata.codes, hit_paths=pdata.hit_paths,
                          **{k: getattr(pdata, k)[perm] for k in
                             ("best_bits", "best_target", "best_ident", "best_tlen", "second_bits",
                              "rbh", "dup_count")})
    cal2 = panel.calibrate(pd2, toy["pc"])
    np.testing.assert_array_equal(cal.d, cal2.d)
    assert cal.core_size == cal2.core_size


def test_scoring_identical_to_self_search():
    assert panel.align_scoring_args(CCFG) == search.scoring_args(CCFG)
    pre = panel.prefilter_scoring_args(CCFG)
    for flag in ("--sub-mat", "--comp-bias-corr", "--comp-bias-corr-scale", "--seed-sub-mat",
                 "--mask", "--mask-prob", "--max-seq-len"):
        assert flag in pre


def test_panel_settings_in_schema_hash(toy):
    h0 = toy["meta"]["schema_hash"]
    assert h0 != pipeline.schema_hash(CCFG)                      # Tier B changes the hash
    pc = copy.deepcopy(toy["pc"])
    sha = toy["meta"]["tier_b"]["search"]["db_stamp"]["panel_fasta_sha256"]
    base = panel.settings_hash(pc, CCFG, sha)
    pc["tier_b"]["presence_min_bits"] = 51
    assert panel.settings_hash(pc, CCFG, sha) != base
    c2 = copy.deepcopy(CCFG)
    c2["mmseqs"]["scoring"]["comp_bias_corr"] = 0
    assert panel.settings_hash(toy["pc"], c2, sha) != base


def test_kimura():
    assert panel.kimura(1.0) == 0.0
    assert panel.kimura(0.9) == pytest.approx(-np.log(1 - 0.1 - 0.002))
    assert panel.kimura(0.1) == 3.0 and panel.kimura(0.16) == 3.0   # capped, monotone
    assert np.isnan(panel.kimura(np.nan))


def _assert_no_constant(df):
    bad = []
    for f in registry.FEATURES:
        if f.tier != "B":
            continue
        v = df[f.name].to_numpy(np.float64)
        if not (np.nanvar(v) > 0 if np.isfinite(v).any() else False):
            bad.append(f.name)
    assert not bad, f"proteome-constant Tier B features: {bad}"


def test_no_proteome_constant_feature_toy(toy):
    _assert_no_constant(toy["df"])


YEAST_TB = PROJECT_ROOT / "data/interim/features/yeast/tierB.parquet"


@pytest.mark.data
@pytest.mark.skipif(not YEAST_TB.exists(), reason="run scripts/features/run_tierB.py yeast first")
def test_no_proteome_constant_feature_yeast():
    _assert_no_constant(pd.read_parquet(YEAST_TB))


# ----------------------------------------------------------------------------------------------
# run_tierB.py row alignment, regression, QC flags
# ----------------------------------------------------------------------------------------------

def _load_run_tierB():
    spec = importlib.util.spec_from_file_location("run_tierB_script",
                                                  PROJECT_ROOT / "scripts/features/run_tierB.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_run_tierB_alignment_noncanonical_fasta_order(toy):
    """Tier B rows built by run_tierB.py from a Tier A-like pair table (shuffled rows) and a FASTA in
    non-canonical order must equal the pipeline.run rows for the same gene pairs."""
    tmp = toy["tmp"] / "script_order"
    tmp.mkdir()
    seqs_by_id, _ = fasta.read_proteome(toy["tmp"] / "query.fasta")
    ids_canon = fasta.canonical_order(seqs_by_id)
    rng = np.random.default_rng(1)
    shuffled = [ids_canon[k] for k in rng.permutation(len(ids_canon))]
    assert shuffled != ids_canon and shuffled != sorted(ids_canon)   # FASTA order != canonical order
    fa = tmp / "shuffled.fasta"
    fasta.write_fasta(fa, ((g, seqs_by_id[g]) for g in shuffled))
    file_order, _ = fasta.read_proteome(fa)                           # dict in FASTA order
    assert list(file_order) != ids_canon
    # Tier A-like table: pipeline pairs in shuffled row order, half of the rows orientation-swapped
    df = toy["df"]
    perm = rng.permutation(len(df))
    ta = df.iloc[perm][["gene_a", "gene_b"]].reset_index(drop=True)
    sw = np.arange(len(ta)) % 2 == 1
    ta.loc[sw, ["gene_a", "gene_b"]] = ta.loc[sw, ["gene_b", "gene_a"]].to_numpy()

    mod = _load_run_tierB()
    wd = tmp / "work"
    wd.mkdir()
    out, pairbits, meta, tool, timings = mod.compute_tierB(
        file_order, ta, CCFG, toy["pc"], wd, threads=1, species_jobs=1)
    # row order and identity follow the given table
    assert (out.gene_a.to_numpy() == ta.gene_a.to_numpy()).all()
    assert (out.gene_b.to_numpy() == ta.gene_b.to_numpy()).all()
    assert list(out.columns) == ["gene_a", "gene_b"] + registry.TIER_B_COLUMNS
    # pos mapping: positions refer to the canonical order, not to FASTA order
    i, j, pairs = mod.pair_positions(ids_canon, ta.gene_a, ta.gene_b)
    assert [ids_canon[k] for k in i] == list(ta.gene_a) and [ids_canon[k] for k in j] == list(ta.gene_b)
    assert (pairs.i <= pairs.j).all()
    # values equal pipeline.run for the same unordered pair (feature symmetry => swap-invariant)
    ref = df.set_index(["gene_a", "gene_b"])[registry.TIER_B_COLUMNS]
    key = [(a, b) if a < b else (b, a) for a, b in zip(out.gene_a, out.gene_b)]
    exp = ref.loc[key].reset_index(drop=True)
    pd.testing.assert_frame_equal(out[registry.TIER_B_COLUMNS].reset_index(drop=True), exp,
                                  check_exact=False, rtol=1e-9, atol=1e-12)
    # pairbits sidecar follows the table orientation: swapping a row swaps the two directions
    assert (pairbits.gene_a.to_numpy() == ta.gene_a.to_numpy()).all()
    wd2 = tmp / "work2"
    wd2.mkdir()
    ta0 = ta.copy()
    ta0.loc[sw, ["gene_a", "gene_b"]] = ta.loc[sw, ["gene_b", "gene_a"]].to_numpy()   # all gene_a < gene_b
    _, pb0, _, _, _ = mod.compute_tierB(file_order, ta0, CCFG, toy["pc"], wd2, threads=1)
    assert (pb0.bits_a_to_b > 0).all()
    np.testing.assert_allclose(pairbits.bits_a_to_b.to_numpy()[sw], pb0.bits_b_to_a.to_numpy()[sw])
    np.testing.assert_allclose(pairbits.bits_a_to_b.to_numpy()[~sw], pb0.bits_a_to_b.to_numpy()[~sw])
    assert meta["threads"] == 1 and meta["species_jobs"] == 1


def test_pair_positions_unknown_gene_raises():
    mod = _load_run_tierB()
    with pytest.raises(ValueError):
        mod.pair_positions(["A", "B"], ["A"], ["ZZ"])


EXPECTED = FIX / "expected_tierB_toy.parquet"


def test_end_to_end_regression_stored_table(toy):
    """pipeline.run(tier_b=True) on the toy proteome equals the stored expected table."""
    cols = ["gene_a", "gene_b"] + registry.TIER_B_COLUMNS
    got = toy["df"][cols].reset_index(drop=True)
    assert EXPECTED.exists(), "regenerate with tests/fixtures/make_expected_tierB_toy.py"
    pd.testing.assert_frame_equal(got, pd.read_parquet(EXPECTED), check_exact=False,
                                  rtol=1e-9, atol=1e-12)


def test_meta_logs_threads_and_species_jobs(toy):
    tb = toy["meta"]["tier_b"]
    assert tb["threads"] == 2 and tb["species_jobs"] == 1
    assert toy["meta"]["threads"] == 2


def test_qc_flags_toy_clean_and_not_features(toy):
    qc = toy["meta"]["tier_b"]["qc"]
    assert toy["meta"]["qc_flags"] == qc
    assert qc["warnings"] == [] and not qc["any"]
    assert qc["flags"] == {"core_set_small": False, "nearest_species_far": False,
                           "no_usable_species": False}
    assert qc["nearest_species"] == "S1" and qc["core_size"] >= 30
    # QC flags are metadata: no QC/flag column in the output beyond registered features
    assert not [c for c in toy["df"].columns if "warn" in c or c.startswith("qc") or "core_set" in c
                or "nearest_species_far" in c]
    assert list(toy["df"].columns).count("tb_depth") == 1


def test_qc_flags_trigger(toy):
    out, pdata, cal, mask, _ = _stage(toy, internal=True)
    # |G*| < 200 after fallback (real-panel threshold)
    q = panel.qc_flags(cal, toy["pc"], min_core=200)
    assert q["flags"]["core_set_small"] and q["any"] and "core set" in q["warnings"][0]
    # nearest retained species d_s > threshold (S1: d = Kimura(0.90) ~ 0.11)
    assert not panel.qc_flags(cal, toy["pc"])["flags"]["nearest_species_far"]
    q = panel.qc_flags(cal, toy["pc"], nearest_d_max=0.05)
    assert q["flags"]["nearest_species_far"] and q["nearest_species"] == "S1"
    assert q["d_nearest"] == pytest.approx(toy["exp"]["d"]["S1"], abs=1e-5)
    # config threshold path: core.min_size above the toy core size
    pc = copy.deepcopy(toy["pc"])
    pc["tier_b"]["core"]["min_size"] = 1000
    assert panel.qc_flags(cal, pc)["flags"]["core_set_small"]
    # default far-threshold is the reviewed 0.3
    assert panel.QC_NEAREST_D_MAX == 0.3
