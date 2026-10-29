"""Toy known-answer tests for the robustness-test evaluation (slparalog.eval.robustness)."""

import numpy as np
import pandas as pd
import pytest

from slparalog.eval import perfold
from slparalog.eval import robustness as R
from slparalog.eval.bootstrap import make_bootstrap
from slparalog.eval.metrics import ap


def _table(pairs, vals):
    """Feature table keyed by gene pair with one registered-style column ``f1``."""
    t = pd.DataFrame(pairs, columns=["gene_a", "gene_b"])
    t["f1"] = vals
    return t


def test_degraded_frame_lost_gene_removed_not_called_and_renamed():
    pairs = pd.DataFrame({"gene_a": ["a", "a", "b", "c", "d"], "gene_b": ["b", "c", "c", "d", "e"]})
    mapping = {"a": "a", "b": "b", "c": "c", "d": "d"}          # e has no degraded record: removed
    table = _table([("a", "b"), ("c", "b"), ("c", "d")], [1.0, 2.0, 3.0])  # (a,c) no longer called; reversed order ok
    fr, cnt = R.degraded_frame(pairs, mapping, table)
    assert fr["_lost"].tolist() == [False, True, False, False, True]
    assert fr["_gene_lost"].tolist() == [False, False, False, False, True]
    assert fr.loc[~fr["_lost"], "f1"].tolist() == [1.0, 2.0, 3.0]
    assert np.isnan(fr.loc[fr["_lost"], "f1"]).all()
    assert cnt == {"n_pairs": 5, "n_lost": 2, "n_lost_gene_removed": 1, "n_lost_not_called": 1}


def test_degraded_frame_follows_mapping_to_new_ids():
    pairs = pd.DataFrame({"gene_a": ["a"], "gene_b": ["b"]})
    table = _table([("a2", "b")], [7.0])
    fr, cnt = R.degraded_frame(pairs, {"a": "a2", "b": "b"}, table)
    assert fr["f1"].tolist() == [7.0] and cnt["n_lost"] == 0
    fr, cnt = R.degraded_frame(pairs, {"a": "a", "b": "b"}, table)          # old ID not in the table: no longer called
    assert cnt["n_lost_not_called"] == 1


def test_load_degraded_uses_mapping_table(tmp_path):
    table = pd.DataFrame({"gene_a": ["a", "b"], "gene_b": ["b", "c"], "ident_global": [0.5, 0.6]})
    table.to_parquet(tmp_path / "f.parquet")
    pd.DataFrame({"role": ["gene", "gene", "derived_haplotig", "gene"], "orig_id": ["a", "b", "a", "c"],
                  "degraded_id": ["a", "b", "a_hap", ""]}).to_csv(tmp_path / "m.tsv", sep="\t", index=False)
    pairs = pd.DataFrame({"gene_a": ["a", "b"], "gene_b": ["b", "c"]})
    fr, cnt = R.load_degraded(pairs, tmp_path / "f.parquet", tmp_path / "m.tsv")
    assert fr["_lost"].tolist() == [False, True] and fr["_gene_lost"].tolist() == [False, True]
    assert cnt["n_lost_gene_removed"] == 1


def test_fill_min_one_tie_block_below_everything():
    s = np.array([0.3, np.nan, 0.9, np.nan])
    lost = np.isnan(s)
    f = R.fill_min(s, lost)
    assert f[1] == f[3] == pytest.approx(-0.7) and f[0] == 0.3 and f[2] == 0.9
    m = R.fill_min(np.array([[0.2, np.nan, 0.5], [0.4, np.nan, 0.6]]), np.array([False, True, False]))
    assert m[:, 1].tolist() == pytest.approx([-0.8, -0.6])
    assert s[1] != s[1]   # input untouched


def _toy(n_groups=30, per=4, seed=0):
    rng = np.random.default_rng(seed)
    groups = [f"g{i}" for i in range(n_groups) for _ in range(per)]
    n = len(groups)
    sp = np.array(["human"] * (n // 2) + ["yeast"] * (n - n // 2))
    df = pd.DataFrame({"group": groups, "species": sp, "label": rng.integers(0, 2, n)})
    df.loc[[0, 1, per * n_groups // 2, per * n_groups // 2 + 1], "label"] = [1, 0, 1, 0]
    gid = df["group"].str[1:].astype(int).to_numpy()
    for s in (11, 23):
        perm = np.random.default_rng(s).permutation(n_groups)
        df[f"fold_s{s}"] = perm[gid] % 3
    return df


def test_perfold_ap_point_matches_harness_perfold_point():
    df = _toy()
    rng = np.random.default_rng(1)
    seeds = [11, 23]
    mat = rng.random((2, len(df))) + df["label"].to_numpy() * 0.3
    base = rng.random(len(df))
    mine = R.perfold_ap_point(df, {"m": mat, "b": base}, seeds, ("human", "yeast"))
    ref = perfold.point(df, {"m": mat}, base, seeds)["m"]["per_species"]
    for sp in ("human", "yeast"):
        assert mine["m"][sp]["ap"] == pytest.approx(ref[sp]["ap_model"])
        assert mine["b"][sp]["ap"] == pytest.approx(ref[sp]["ap_base"])
        assert mine["m"][sp]["n_folds"] == ref[sp]["n_folds"]


def test_bootstrap_delta_matches_harness_perfold_bootstrap_and_drop_summary_algebra():
    df = _toy()
    rng = np.random.default_rng(2)
    seeds = [11, 23]
    y = df["label"].to_numpy()
    nat = rng.random((2, len(df))) + y * 0.6
    deg = rng.random((2, len(df))) + y * 0.1            # degraded model is worse
    b_nat = rng.random(len(df)) + y * 0.2
    b_deg = b_nat.copy()                                 # B1 unaffected
    bs = make_bootstrap(df, 200, 7)
    sc = {"model|N": nat, "b1|N": b_nat, "model|k": deg, "b1|k": b_deg}
    bt = R.perfold_ap_bootstrap(bs, df, sc, seeds, ("human", "yeast"))
    ref = perfold.bootstrap(bs, df, {"n": nat, "d": deg}, b_nat, seeds)       # same base for both (B1 unchanged)
    for sp in ("human", "yeast"):
        assert np.allclose(bt["model|N"][sp] - bt["b1|N"][sp], ref["n"][sp], equal_nan=True)
        assert np.allclose(bt["model|k"][sp] - bt["b1|k"][sp], ref["d"][sp], equal_nan=True)
    pt = R.perfold_ap_point(df, sc, seeds, ("human", "yeast"))
    d = R.drop_summary(pt, bt, "k", ("human", "yeast"))
    for sp in ("human", "yeast"):
        x = d[sp]
        assert x["drop_model"] == pytest.approx(x["ap_model_natural"] - x["ap_model_degraded"])
        assert x["drop_b1"] == pytest.approx(0.0, abs=1e-12) and x["drop_b1_ci"]["lo"] == 0 == x["drop_b1_ci"]["hi"]
        assert x["drop_delta"] == pytest.approx(x["drop_model"] - x["drop_b1"])
        assert x["drop_model"] > 0 and x["drop_model_ci"]["lo"] <= x["drop_model"] <= x["drop_model_ci"]["hi"] + 0.05
    t = d["T_pf"]
    assert t["drop"] == pytest.approx(np.mean([d[sp]["drop_delta"] for sp in ("human", "yeast")]))
    assert t["natural"] == pytest.approx(np.mean([d[sp]["delta_natural"] for sp in ("human", "yeast")]))
    d1 = R.drop_summary(pt, bt, "k", ("yeast",))
    assert "T_pf" not in d1 and set(d1) == {"yeast"}


def test_identity_degradation_gives_zero_drop():
    df = _toy()
    rng = np.random.default_rng(3)
    seeds = [11, 23]
    nat = rng.random((2, len(df)))
    b = rng.random(len(df))
    bs = make_bootstrap(df, 100, 1)
    sc = {"model|N": nat, "b1|N": b, "model|k": nat.copy(), "b1|k": b.copy()}
    d = R.drop_summary(R.perfold_ap_point(df, sc, seeds, ("human", "yeast")),
                       R.perfold_ap_bootstrap(bs, df, sc, seeds, ("human", "yeast")), "k", ("human", "yeast"))
    assert d["T_pf"]["drop"] == 0 and d["T_pf"]["drop_ci"]["lo"] == 0 == d["T_pf"]["drop_ci"]["hi"]


def test_lost_pairs_at_minimum_lower_ap_by_known_amount():
    # 6 pairs in one fold; positives are pairs 0 and 3. Natural scores rank both on top.
    y = np.array([1, 0, 0, 1, 0, 0])
    nat = np.array([0.9, 0.1, 0.2, 0.8, 0.3, 0.05])
    assert ap(y, nat) == pytest.approx(1.0)
    lost = np.array([False, False, False, True, False, False])            # positive 3 is lost
    deg = R.fill_min(np.where(lost, np.nan, nat), lost)
    # ranking: 0.9(1), 0.3(0), 0.2(0), 0.1(0), 0.05(0), lost(1): AP = (1/1 + 2/6) / 2
    assert ap(y, deg) == pytest.approx((1.0 + 2 / 6) / 2)
