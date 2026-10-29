"""G6: frozen test/development split (DESIGN.md v3 §6.1, §6.3; reviews/round2/agent_failure_modes.md G6;
consistency M11). Allowlisted (D4/D5) to read sealed split files and counts; prints no pair IDs.

Toy tests (no marker): weighted AP == sklearn on duplicated rows, supernode contraction, Leiden fallback,
greedy stratified draw, power simulation separates a null from a strong effect.
Data tests (marker ``data``): byte-identical re-run, zero test<->dev gene overlap, no labelled pair (any tier)
spans the split, per-species test positives within 1 of the target (D9; yeast >= 30 unless the
D4 yeast fallback applies), sealed files match MANIFEST.sha256 and protected hashes match.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MS = _load("make_splits", "scripts/labels/make_splits.py")
PS = _load("power_sim", "scripts/eval/power_sim.py")
CFG = MS.load_cfg()
SEALED = ROOT / CFG["outputs"]["sealed_dir"]
DEV = ROOT / CFG["outputs"]["dev_dir"]
SUMMARY = ROOT / CFG["outputs"]["summary"]


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# ============================================================================ toy tests (known answers)


def test_weighted_ap_equals_sklearn_on_duplicated_rows():
    rng = np.random.default_rng(0)
    for _ in range(20):
        n = 60
        y = (rng.random(n) < 0.2).astype(float)
        y[0] = 1.0
        s = rng.standard_normal(n)
        w = rng.integers(0, 4, n)
        w[0] = max(w[0], 1)
        o = np.argsort(-s, kind="mergesort")
        got = float(PS.weighted_ap(y[o], w[o]))
        rep = np.repeat(np.arange(n), w)
        want = average_precision_score(y[rep], s[rep])
        assert got == pytest.approx(want, abs=1e-12)
    # unit weights == plain AP; known answer from G5
    y = np.array([1, 0, 1, 0.0])
    assert float(PS.weighted_ap(y, np.ones(4))) == pytest.approx(5 / 6)
    # 2-D weights, one row with no positives -> nan
    W = np.array([[1, 1, 1, 1], [0, 1, 0, 1]])
    out = PS.weighted_ap(y, W)
    assert out[0] == pytest.approx(5 / 6) and np.isnan(out[1])


def test_supernodes_contract_every_labelled_pair():
    labels = {
        "human": pd.DataFrame({"key_a": ["human:A", "human:B", "human:D"], "key_b": ["human:B", "human:C", "human:E"]}),
        "yeast": pd.DataFrame({"key_a": ["yeast:X"], "key_b": ["yeast:Y"]}),
    }
    sn = MS.supernodes(labels)
    assert sn["human:A"] == sn["human:B"] == sn["human:C"] == "human:A"
    assert sn["human:D"] == sn["human:E"] == "human:D"
    assert sn["yeast:X"] == sn["yeast:Y"]
    assert len(set(sn.values())) == 3


def test_leiden_fallback_splits_two_cliques_and_respects_cap():
    pytest.importorskip("leidenalg")
    a = [f"human:a{i}" for i in range(5)]
    b = [f"human:b{i}" for i in range(5)]
    edges = [(x, y, 0.9) for grp in (a, b) for i, x in enumerate(grp) for y in grp[i + 1:]]
    edges.append((a[0], b[0], 0.3))
    nodes = sorted(a + b)
    edges = [(min(x, y), max(x, y), w) for x, y, w in edges]
    comms, gamma, ok = MS.leiden_partition(nodes, edges, lambda c: len(c) <= 5, CFG, 11)
    assert ok
    assert sorted(map(sorted, comms)) == sorted([sorted(a), sorted(b)])
    # deterministic
    assert MS.leiden_partition(nodes, edges, lambda c: len(c) <= 5, CFG, 11) == (comms, gamma, ok)


def test_greedy_draw_hits_stratified_targets_on_toy_groups():
    rng = np.random.default_rng(3)
    G = 400
    gc = pd.DataFrame({"n_human": rng.integers(1, 6, G), "n_yeast": rng.integers(0, 4, G)})
    gc["pos_human"] = rng.binomial(gc["n_human"], 0.1)
    gc["pos_yeast"] = rng.binomial(gc["n_yeast"], 0.25)   # >= 75 yeast positives: no yeast fallback
    assert gc["pos_yeast"].sum() >= 75
    gc["neg_human"] = gc["n_human"] - gc["pos_human"]
    gc["neg_yeast"] = gc["n_yeast"] - gc["pos_yeast"]
    eligible = np.ones(G, bool)
    eligible[:5] = False
    for f in (0.30, 0.35, 0.40):
        t = MS.draw_test(CFG, gc, eligible, f, 11)
        assert not t[:5].any(), "forced-development groups must never be drawn"
        for sp in ("human", "yeast"):
            frac = gc.loc[t, f"pos_{sp}"].sum() / gc[f"pos_{sp}"].sum()
            assert abs(frac - f) < 0.03
        assert np.array_equal(t, MS.draw_test(CFG, gc, eligible, f, 11))
    # yeast fallback (< 75 yeast positives): the yeast target is ceil(0.40 * P) whatever f is
    gc2 = gc.copy()
    gc2["pos_yeast"] = (np.arange(G) % 8 == 0).astype(int) * (gc2["n_yeast"] > 0)
    gc2["neg_yeast"] = gc2["n_yeast"] - gc2["pos_yeast"]
    P = int(gc2["pos_yeast"].sum())
    assert P < 75
    t = MS.draw_test(CFG, gc2, eligible, 0.30, 11)
    assert abs(int(gc2.loc[t, "pos_yeast"].sum()) - math.ceil(0.40 * P)) <= 1


def test_yeast_fallback_target():
    t, fb = MS.positive_targets(CFG, 0.30, {"human": 200, "yeast": 60})
    assert fb and t["yeast"] == math.ceil(0.40 * 60) and t["human"] == pytest.approx(60.0)
    t, fb = MS.positive_targets(CFG, 0.30, {"human": 200, "yeast": 100})
    assert not fb and t["yeast"] == pytest.approx(30.0)


def test_power_sim_null_vs_strong_effect_on_toy_composition():
    rng = np.random.default_rng(5)
    G = 300
    g = np.repeat(np.arange(G), 3)
    sp = (np.arange(G) % 3 == 0).astype(int)[g]   # every third group is yeast
    y = (rng.random(len(g)) < 0.12).astype(float)
    has_pos = np.bincount(g, weights=y, minlength=G) > 0
    comp_sp = np.bincount(g, weights=sp, minlength=G) > 0
    stratum = np.where(comp_sp, 1, 0) * 2 + has_pos
    comp = {"y": y, "sp": sp, "g": g, "n_groups": G, "stratum": stratum}
    cfg = {"cluster_score_sd": 0.5, "bootstrap_resamples": 200, "ci_level": 0.95, "n_sim": 40}
    mu1 = np.array([PS.mu_of_auc(0.73, 0.5)] * 2)
    null = PS.simulate((comp, mu1, mu1.copy(), 0.8, cfg, np.random.SeedSequence(1), 0.2))
    strong = PS.simulate((comp, mu1, np.array([PS.mu_of_auc(0.93, 0.5)] * 2), 0.8, cfg, np.random.SeedSequence(1), 0.2))
    assert null["power_rule"] <= 0.1
    assert strong["power_rule"] >= 0.9
    assert abs(null["mean_T"]) < 0.02
    # deterministic given the seed
    again = PS.simulate((comp, mu1, mu1.copy(), 0.8, cfg, np.random.SeedSequence(1), 0.2))
    assert again == null


def test_bootstrap_counts_preserve_stratum_sizes():
    st = np.array([0, 0, 0, 1, 1, 2, 3, 3, 3, 3])
    W = PS.bootstrap_counts(np.random.default_rng(0), st, 50)
    for v in np.unique(st):
        assert (W[:, st == v].sum(axis=1) == (st == v).sum()).all()


# ============================================================================ data tests (real split)

data = pytest.mark.data
needs_split = pytest.mark.skipif(not (SEALED / CFG["outputs"]["groups"]).exists(), reason="split not built")


def _labels():
    return MS.load_all_labels(CFG)


def _split_tables():
    test = {sp: pd.read_csv(SEALED / f"{sp}_test.tsv", sep="\t", dtype=str, keep_default_na=False) for sp in MS.SPECIES}
    dev = {sp: pd.read_csv(DEV / f"{sp}_dev.tsv", sep="\t", dtype=str, keep_default_na=False) for sp in MS.SPECIES}
    return test, dev


@data
@needs_split
def test_rerun_is_byte_identical(tmp_path):
    """Re-running both stages on the frozen inputs (cached search, frozen results/mde.json) reproduces every
    split output byte for byte."""
    cache = ROOT / CFG["outputs"]["interim_dir"]
    for stage in ("candidates", "finalize"):
        subprocess.run([PY, str(ROOT / "scripts/labels/make_splits.py"), stage, "--out-root", str(tmp_path),
                        "--search-cache", str(cache), "--mde", str(ROOT / "results/mde.json")],
                       check=True, capture_output=True)
    o = CFG["outputs"]
    files = [Path(o["sealed_dir"]) / n for n in (o["candidates"], o["groups"], o["test_pairs"],
                                                 "human_test.tsv", "yeast_test.tsv")]
    files += [Path(o["dev_dir"]) / f"{sp}_dev.tsv" for sp in MS.SPECIES]
    for rel in files:
        assert (tmp_path / rel).read_bytes() == (ROOT / rel).read_bytes(), f"{rel} differs on re-run"
    a = json.loads((tmp_path / o["summary"]).read_text())
    b = json.loads(SUMMARY.read_text())
    a.pop("mde"), b.pop("mde")
    assert a == b


@data
@needs_split
def test_zero_gene_overlap_between_test_and_dev():
    test, dev = _split_tables()
    for sp in MS.SPECIES:
        tg = set(test[sp].gene_a) | set(test[sp].gene_b)
        dg = set(dev[sp].gene_a) | set(dev[sp].gene_b)
        assert not (tg & dg), f"{sp}: {len(tg & dg)} genes in both test and development"
    groups = pd.read_csv(SEALED / CFG["outputs"]["groups"], sep="\t")
    seen = {}
    for gid, a, gs in zip(groups.group_id, groups.assignment, groups.genes):
        for g in gs.split(","):
            assert g not in seen, "gene in two groups"
            seen[g] = a
    for sp in MS.SPECIES:
        assert {seen[f"{sp}:{g}"] for g in set(test[sp].gene_a) | set(test[sp].gene_b)} == {"test"}
        assert {seen[f"{sp}:{g}"] for g in set(dev[sp].gene_a) | set(dev[sp].gene_b)} == {"dev"}


@data
@needs_split
def test_no_labelled_pair_of_any_tier_spans_the_split():
    labels = _labels()
    test, dev = _split_tables()
    for sp in MS.SPECIES:
        lab = {tuple(x) for x in labels[sp][["gene_a", "gene_b"]].to_numpy()}
        t = {tuple(x) for x in test[sp][["gene_a", "gene_b"]].to_numpy()}
        d = {tuple(x) for x in dev[sp][["gene_a", "gene_b"]].to_numpy()}
        assert not (t & d)
        assert t | d == lab, f"{sp}: labelled pairs missing from or extra in the split tables"
        assert len(t) == len(test[sp]) and len(d) == len(dev[sp])
        # both genes of every pair carry the same group
        for tab in (test[sp], dev[sp]):
            assert tab["group"].ne("").all()
    # every pair's two genes are in one group (re-derived from the group table)
    groups = pd.read_csv(SEALED / CFG["outputs"]["groups"], sep="\t")
    g2g = MS.gene_to_group(groups)
    for sp in MS.SPECIES:
        assert (labels[sp]["key_a"].map(g2g) == labels[sp]["key_b"].map(g2g)).all()


@data
@needs_split
def test_cv_folds_are_grouped_and_cover_every_seed():
    _, dev = _split_tables()
    seeds = MS.load_seeds()
    both = pd.concat(dev.values(), ignore_index=True)
    for s in seeds:
        col = f"fold_s{s}"
        assert set(both[col].astype(int)) == set(range(int(CFG["cv"]["n_folds"])))
        assert (both.groupby("group")[col].nunique() == 1).all(), "a group spans two CV folds"


@data
@needs_split
def test_test_positive_counts_on_target():
    """D9: per species |test primary positives - target| <= positive_target_margin (1), target = f * P_s with
    f = the chosen fraction (yeast: ceil(0.40 * P) under the D4 fallback); yeast >= 30 test positives
    unless the fallback applies."""
    labels = _labels()
    test, _ = _split_tables()
    ta = CFG["test_assignment"]
    margin = float(ta["positive_target_margin"])
    mde = json.loads((ROOT / "results/mde.json").read_text())
    f = float(mde["chosen_fraction"])
    assert f in [float(x) for x in ta["candidate_fractions"]]
    assert mde["candidates_sha256"] == _sha(SEALED / CFG["outputs"]["candidates"])
    pc = {sp: CFG["labels"][sp]["primary"] for sp in MS.SPECIES}
    P = {sp: int((labels[sp].y_primary == 1).sum()) for sp in MS.SPECIES}
    tp = {sp: int((test[sp][pc[sp]["column"]] == str(pc[sp]["pos"])).sum()) for sp in MS.SPECIES}
    yeast_fallback = P["yeast"] < int(ta["yeast_fallback"]["min_primary_positives"])
    target = {sp: f * P[sp] for sp in MS.SPECIES}
    if yeast_fallback:
        target["yeast"] = math.ceil(float(ta["yeast_fallback"]["target_fraction"]) * P["yeast"])
    for sp in MS.SPECIES:
        assert abs(tp[sp] - target[sp]) <= margin, f"{sp}: {tp[sp]} test positives vs target {target[sp]:.1f}"
    if not yeast_fallback:
        assert tp["yeast"] >= int(ta["min_yeast_test_positives"])


@data
@needs_split
def test_groups_above_pair_cap_are_in_development():
    """D9: a group holding > cap_fraction of a species' labelled pairs (any tier) is forced into development."""
    labels = _labels()
    groups = pd.read_csv(SEALED / CFG["outputs"]["groups"], sep="\t")
    capf = float(CFG["size_cap"]["cap_fraction"])
    assert CFG["size_cap"]["unit"] == "labelled_pairs_any_tier" and not CFG["size_cap"]["leiden"]["enabled"]
    over = np.zeros(len(groups), bool)
    for sp in MS.SPECIES:
        over |= groups[f"n_{sp}"].to_numpy() > capf * len(labels[sp])
    assert (groups.loc[over, "assignment"] == "dev").all()
    assert (groups["oversized"].to_numpy().astype(bool) == over).all()
    assert (groups["leiden"] == 0).all()


@data
@needs_split
def test_sealed_manifest_and_protected_hashes():
    manifest = {}
    for line in (SEALED / "MANIFEST.sha256").read_text().splitlines():
        d, name = line.split(None, 1)
        manifest[name.strip()] = d
    o = CFG["outputs"]
    for name in (o["candidates"], o["groups"], o["test_pairs"], "human_test.tsv", "yeast_test.tsv"):
        assert manifest.get(name) == _sha(SEALED / name), f"{name} does not match MANIFEST.sha256"
    protected = {}
    for line in (ROOT / "config/protected.sha256").read_text().splitlines():
        d, rel = line.split(None, 1)
        protected[rel.strip()] = d
    for rel in ("config/splits.yaml", "config/mde.yaml", "scripts/labels/make_splits.py",
                "scripts/eval/power_sim.py", "tests/test_splits.py", "data/labels/sealed/MANIFEST.sha256"):
        assert rel in protected, f"{rel} is not protected"
        assert protected[rel] == _sha(ROOT / rel), f"{rel} changed since it was protected"
