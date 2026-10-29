"""Toy known-answer tests for the panel-thinning test and masking augmentation (slparalog.eval.thinning)."""

import numpy as np
import pandas as pd
import pytest

from slparalog.eval import thinning as T
from slparalog.eval.cv import oof_scores
from slparalog.eval.models import GBDTModel


def test_augment_spreads_weight_and_dedups_identical_copies():
    nat = np.array([[1.0, np.nan], [2.0, 3.0], [4.0, 5.0]])
    thin = np.array([[1.0, np.nan], [2.0, 9.0], [4.0, 0.0]])   # pair 0 identical (NaN == NaN), pairs 1, 2 differ
    y = np.array([1, 0, 0])
    w = np.array([3.0, 2.0, 1.0])
    X, yy, ww, n, k = T.augment([nat, thin], y, w)
    assert n == 3 and k.tolist() == [1, 2, 2]
    assert X.shape == (5, 2) and len(yy) == 5 and len(ww) == 5
    # pair weight preserved: rows are [natural pairs 0..2] then [thinned copies of pairs 1, 2]
    assert ww.tolist() == [3.0, 1.0, 0.5, 1.0, 0.5]
    assert ww.sum() == pytest.approx(w.sum())
    assert yy.tolist() == [1, 0, 0, 0, 0]
    assert X[3].tolist() == [2.0, 9.0] and X[4].tolist() == [4.0, 0.0]


def test_augment_three_views_dedup_against_all_earlier_copies():
    a = np.array([[1.0], [1.0]])
    b = np.array([[2.0], [1.0]])    # pair 1 identical to a
    c = np.array([[2.0], [2.0]])    # pair 0 identical to b, pair 1 new
    X, _, ww, _, k = T.augment([a, b, c], np.array([1, 0]), np.array([1.0, 1.0]))
    assert k.tolist() == [2, 2]
    assert len(X) == 4 and ww.tolist() == [0.5, 0.5, 0.5, 0.5]


def test_augment_single_view_is_identity():
    X0 = np.arange(6.0).reshape(3, 2)
    X, y, w, n, k = T.augment([X0], np.array([0, 1, 0]), np.array([1.0, 2.0, 3.0]))
    assert np.array_equal(X, X0) and w.tolist() == [1.0, 2.0, 3.0] and n == 3 and k.tolist() == [1, 1, 1]


def test_leaf_in_pairs():
    assert T.leaf_in_pairs(50, 100, 100) == 50
    assert T.leaf_in_pairs(50, 200, 100) == 100
    assert T.leaf_in_pairs(50, 150, 100) == 75
    assert T.leaf_in_pairs(50, 100, 100) == 50 and T.leaf_in_pairs(7, 101, 100) == 8  # ceil


def test_condition_set_names():
    cs = T.condition_set([0.03, 0.35])
    assert list(cs) == ["N->N", "N->T0.03", "T0.03->T0.03", "A0.03->N", "A0.03->T0.03",
                        "N->T0.35", "T0.35->T0.35", "A0.35->N", "A0.35->T0.35"]
    assert cs["A0.35->T0.35"] == (["N", "T0.35"], "T0.35")


def test_decide_known_answers():
    rng = np.random.default_rng(1)
    base = rng.normal(0, 0.01, 400)
    boot = {"N->N": 0.14 + base, "N->T0.35": 0.08 + base * 0.5 + rng.normal(0, 0.005, 400),
            "A0.35->N": 0.135 + base, "A0.35->T0.35": 0.12 + base}
    Tpt = {k: float(np.mean(v)) for k, v in boot.items()}
    d = T.decide(Tpt, boot, 0.35)
    assert d["drop"] == pytest.approx(0.06, abs=0.003) and d["cond1_drop_ci_excludes_0"]
    assert d["recovered_fraction"] == pytest.approx(0.04 / 0.06, abs=0.1) and d["cond2_recovers_half"]
    assert d["cond3_natural_loss_le_0.01"] and d["adopt"]
    # recovery too small
    d2 = T.decide({**Tpt, "A0.35->T0.35": Tpt["N->T0.35"] + 0.01}, boot, 0.35)
    assert not d2["cond2_recovers_half"] and not d2["adopt"]
    # natural panel lowered by more than 0.01
    d3 = T.decide({**Tpt, "A0.35->N": Tpt["N->N"] - 0.011}, boot, 0.35)
    assert not d3["cond3_natural_loss_le_0.01"] and not d3["adopt"]
    # exactly 0.01 is allowed
    d4 = T.decide({**Tpt, "A0.35->N": Tpt["N->N"] - 0.01}, boot, 0.35)
    assert d4["cond3_natural_loss_le_0.01"] or abs(d4["natural_loss"] - 0.01) < 1e-12
    # no drop: CI includes 0 -> not adopted even with perfect recovery
    flat = {"N->N": 0.14 + base, "N->T0.35": 0.14 + base + rng.normal(0, 0.01, 400),
            "A0.35->N": 0.14 + base, "A0.35->T0.35": 0.15 + base}
    d5 = T.decide({k: float(np.mean(v)) for k, v in flat.items()}, flat, 0.35)
    assert not d5["cond1_drop_ci_excludes_0"] and not d5["adopt"]
    # thinning helps (negative drop): no drop, no adoption
    neg = {"N->N": np.full(50, 0.1), "N->T0.35": np.full(50, 0.12), "A0.35->N": np.full(50, 0.1),
           "A0.35->T0.35": np.full(50, 0.2)}
    d6 = T.decide({k: 0.1 for k in neg} | {"N->T0.35": 0.12, "A0.35->T0.35": 0.2}, neg, 0.35)
    assert d6["drop"] < 0 and not d6["adopt"] and d6["recovered_fraction"] is None


def _toy(n_groups=24, per=4, seed=0):
    rng = np.random.default_rng(seed)
    groups = [f"g{i}" for i in range(n_groups) for _ in range(per)]
    n = len(groups)
    sp = np.array(["human"] * (n // 2) + ["yeast"] * (n - n // 2))
    df = pd.DataFrame({"group": groups, "species": sp, "label": rng.integers(0, 2, n)})
    df.loc[[0, per * n_groups // 2], "label"] = [1, 0]
    gid = df["group"].str[1:].astype(int).to_numpy()
    for s in (11, 23):
        perm = np.random.default_rng(s).permutation(n_groups)
        df[f"fold_s{s}"] = perm[gid] % 4
    X = rng.random((n, 3))
    X[:, 0] += df["label"].to_numpy() * 0.8
    return df, X


class _Spy:
    """Records the training set it is fitted on; predicts the first feature."""
    log = []

    def __init__(self, params, seed):
        self.p = params

    def fit(self, X, y, sample_weight=None):
        _Spy.log.append({"n": len(y), "w": float(np.sum(sample_weight)), "msl": self.p.get("min_samples_leaf")})
        return self

    def predict(self, X):
        return X[:, 0]


def test_oof_conditions_shapes_weights_and_leaf_in_pairs():
    df, X = _toy()
    thin = X.copy()
    thin[:, 0] *= 0.5
    views = {"N": X, "T0.35": thin}
    conds = T.condition_set([0.35])
    _Spy.log = []
    res = T.oof_conditions(df, views, conds, lambda p, s: _Spy(p, s), {"min_samples_leaf": 10}, [11], n_jobs=1)
    assert set(res) == set(conds) and all(m.shape == (1, len(df)) for m in res.values())
    assert all(np.isfinite(m).all() for m in res.values())       # one score per pair
    # N->N: predictions are the natural first feature; N->T: the thinned one
    assert np.allclose(res["N->N"][0], X[:, 0]) and np.allclose(res["N->T0.35"][0], thin[:, 0])
    assert np.allclose(res["A0.35->T0.35"][0], thin[:, 0])
    # fits: per fold 3 distinct training sets (N, T, N+T); the augmented one has 2x rows, weights sum = pairs
    nf = df["fold_s11"].nunique()
    assert len(_Spy.log) == 3 * nf
    plain = [r for r in _Spy.log if r["msl"] == 10]
    aug = [r for r in _Spy.log if r["msl"] == 20]
    assert len(plain) == 2 * nf and len(aug) == nf
    for r in aug:
        assert r["w"] == pytest.approx(r["n"] / 2)               # sum of weights == number of pairs
    for r in plain:
        assert r["w"] == pytest.approx(r["n"])                    # species-balanced weights average 1


def test_n_to_n_equals_harness_oof_scores_for_gbdt():
    df, X = _toy()
    params = {"max_iter": 20, "max_leaf_nodes": 4, "min_samples_leaf": 5, "learning_rate": 0.1,
              "early_stopping": False}
    cols = ["a", "b", "c"]
    mk = lambda p, s: GBDTModel(cols, p, seed=s)   # noqa: E731
    res = T.oof_conditions(df, {"N": X}, {"N->N": (["N"], "N")}, mk, params, [11, 23], n_jobs=1)
    _, ref = oof_scores(df, X, lambda s: GBDTModel(cols, params, seed=s), [11, 23])
    assert np.allclose(res["N->N"], ref)


def test_augmented_training_never_sees_test_pairs():
    df, X = _toy()
    X = X.copy()
    X[:, 2] = np.arange(len(df))             # row id as a feature
    thin = X.copy()
    thin[:, 0] *= 0.5

    class Seen:
        def __init__(self, p, s):
            pass

        def fit(self, X, y, sample_weight=None):
            self.ids = set(X[:, 2].astype(int).tolist())
            return self

        def predict(self, X):
            return np.array([float(int(i) in self.ids) for i in X[:, 2]])

    res = T.oof_conditions(df, {"N": X, "T0.35": thin}, T.condition_set([0.35]), lambda p, s: Seen(p, s), {}, [11])
    assert all((m == 0).all() for m in res.values())
