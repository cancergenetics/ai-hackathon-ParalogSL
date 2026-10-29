"""G5: pinned metric, toy known answers, same-pair enforcement, config-driven constants."""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import sklearn
import yaml

from slparalog.config import load_success
from slparalog.eval import metrics
from slparalog.eval.metrics import PairSetMismatch, ap, ap_pair

ROOT = Path(__file__).resolve().parents[1]


def test_ap_toy_known_answer():
    assert ap([1, 0, 1, 0], [0.9, 0.8, 0.7, 0.1]) == pytest.approx(5 / 6)


def test_ap_perfect_and_worst():
    assert ap([1, 1, 0, 0], [0.9, 0.8, 0.2, 0.1]) == pytest.approx(1.0)
    # one positive ranked last of 4: AP = 1/4
    assert ap([0, 0, 0, 1], [0.9, 0.8, 0.7, 0.1]) == pytest.approx(0.25)


def test_constant_scores_give_prevalence():
    y = [1, 0, 0, 0, 1, 0, 0, 0, 0, 0]
    assert ap(y, [0.5] * len(y)) == pytest.approx(0.2)


def test_ap_invariant_to_row_order():
    y = np.array([1, 0, 1, 0, 0, 1, 0, 0])
    s = np.array([0.9, 0.8, 0.7, 0.1, 0.35, 0.6, 0.2, 0.05])
    base = ap(y, s)
    rng = np.random.default_rng(0)
    for _ in range(5):
        perm = rng.permutation(len(y))
        assert ap(y[perm], s[perm]) == pytest.approx(base)


def test_ap_matches_sklearn_directly():
    from sklearn.metrics import average_precision_score

    y = [1, 0, 1, 1, 0, 0, 1, 0]
    s = [0.2, 0.4, 0.9, 0.3, 0.8, 0.1, 0.5, 0.6]
    assert ap(y, s) == average_precision_score(y, s)


def test_ap_requires_positives_and_finite_scores():
    with pytest.raises(ValueError):
        ap([0, 0, 0], [0.1, 0.2, 0.3])
    with pytest.raises(ValueError):
        ap([1, 0], [np.nan, 0.1])


def _frames():
    ids = [f"p{i}" for i in range(6)]
    y = [1, 0, 0, 1, 0, 0]
    model = pd.DataFrame({"pair_id": ids, "label": y, "score": [0.9, 0.2, 0.1, 0.8, 0.3, 0.05]})
    base = pd.DataFrame({"pair_id": ids[::-1], "label": y[::-1], "score": [0.2, 0.5, 0.4, 0.3, 0.9, 0.1]})
    return model, base


def test_ap_pair_same_pairs_different_row_order_ok():
    model, base = _frames()
    out = ap_pair(model, base)
    assert out["n"] == 6 and out["n_pos"] == 2 and out["n_neg"] == 4
    assert out["prevalence"] == pytest.approx(2 / 6)
    assert out["delta"] == pytest.approx(out["ap_model"] - out["ap_baseline"])
    # baseline scores computed independently on the same pairs
    b = base.set_index("pair_id").loc[model["pair_id"]]
    assert out["ap_baseline"] == pytest.approx(ap(b["label"], b["score"]))


def test_ap_refuses_mismatched_pair_sets():
    model, base = _frames()
    with pytest.raises(PairSetMismatch):
        ap(model.label, model.score, pair_ids=model.pair_id, reference_pair_ids=base.pair_id[:-1])
    with pytest.raises(PairSetMismatch):
        ap_pair(model, base.iloc[:-1])
    swapped = base.copy()
    swapped.loc[swapped.index[0], "pair_id"] = "other"
    with pytest.raises(PairSetMismatch):
        ap_pair(model, swapped)
    dup = base.copy()
    dup.loc[dup.index[0], "pair_id"] = dup.loc[dup.index[1], "pair_id"]
    with pytest.raises(PairSetMismatch):
        ap_pair(model, dup)
    with pytest.raises(PairSetMismatch):
        ap(model.label, model.score, reference_pair_ids=base.pair_id)  # pair_ids missing


def test_ap_pair_refuses_label_disagreement():
    model, base = _frames()
    base.loc[base.index[0], "label"] = 1 - base.loc[base.index[0], "label"]
    with pytest.raises(PairSetMismatch):
        ap_pair(model, base)


def test_sklearn_version_matches_lock():
    lock = (ROOT / "requirements-lock.txt").read_text()
    m = re.search(r"^scikit-learn==(\S+)$", lock, flags=re.M)
    assert m, "scikit-learn not pinned in requirements-lock.txt"
    assert sklearn.__version__ == m.group(1)


def test_success_constants_come_from_config():
    cfg = load_success()
    raw = yaml.safe_load((ROOT / "config" / "success.yaml").read_text())
    assert cfg == raw
    for key in ("relative_gain", "bootstrap_resamples", "bootstrap_seed", "t_target_ap",
                "power_target", "test_fraction_range"):
        assert key in cfg
    lo, hi = cfg["test_fraction_range"]
    assert 0 < lo < hi < 1
    # the metric module must not carry the criterion constants as literals
    src = (ROOT / "src" / "slparalog" / "eval" / "metrics.py").read_text()
    for lit in ("10000", "10_000", "0.20", "0.04"):
        assert lit not in src
    assert not hasattr(metrics, "BOOTSTRAP_SEED")
