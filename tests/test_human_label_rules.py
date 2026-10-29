"""Toy-data tests with known answers for the human label rules (no raw data needed)."""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "labels"))
hp = importlib.import_module("human_parsers")
bh = importlib.import_module("build_human_labels")
CFG = yaml.safe_load((ROOT / "config" / "labels_human.yaml").read_text())


def test_robust_scale_known_answer():
    m, s, T = hp.robust_scale(pd.Series([0.0, 1.0, -1.0, 2.0, -2.0]), 2.0, 0.5)
    assert m == 0.0 and abs(s - 1.4826) < 1e-9 and abs(T - 2.9652) < 1e-9
    assert hp.robust_scale(pd.Series([0.0, 0.01, -0.01]), 2.0, 0.5)[2] == 0.5   # floor


def test_chou_zscore_ignores_outliers_in_fit():
    rng = np.random.default_rng(0)
    x = np.concatenate([rng.normal(0, 1, 5000), [-50.0]])
    z = hp.chou_zscore(x)
    assert abs(np.mean(z[:-1])) < 0.05 and z[-1] < -40


def _obs(**kw):
    base = dict(study="toy", line_label="L", line_key="L", sym1="A", sym2="B", unit="", gi_strict=False,
                gi_relaxed=False, flag_untestable=False, lfc1=0.0, lfc2=0.0, lfc12=0.0, score_name="s", score=0.0,
                author_pair_class="")
    base.update(kw)
    return base


def test_call_rows_r2_5_rule():
    # singles median 0, MAD tiny -> T = 0.5
    singles = {"L": pd.Series(np.zeros(100))}
    rows = pd.DataFrame([
        _obs(gi_strict=True, gi_relaxed=True, lfc12=-1.0),                    # SL
        _obs(gi_strict=True, gi_relaxed=True, lfc12=-0.4),                    # GI but DKO too shallow, ab<-T/2 -> B
        _obs(gi_strict=True, gi_relaxed=True, lfc12=-0.1),                    # GI, ab > -T/2 -> N
        _obs(gi_strict=False, gi_relaxed=False, lfc12=-2.0),                  # no GI -> N
        _obs(gi_strict=True, gi_relaxed=True, lfc1=-0.3, lfc2=-0.3, lfc12=-2),  # a+b=-0.6 <= -T -> U
        _obs(gi_strict=True, gi_relaxed=True, flag_untestable=True, lfc12=-2),  # study flag -> U
        _obs(gi_strict=True, gi_relaxed=True, lfc12=np.nan),                  # no LFC -> X
        _obs(gi_strict=True, gi_relaxed=True, lfc12=-2, score=np.nan),        # no significance score -> X
    ])
    out = bh.call_rows(rows, singles, CFG, {})
    assert out["call"].tolist() == ["S", "B", "N", "N", "U", "U", "X", "X"]


def test_call_rows_no_lfc_rule():
    rows = pd.DataFrame([_obs(gi_strict=True, gi_relaxed=True), _obs(gi_relaxed=True), _obs(),
                         _obs(gi_strict=True, gi_relaxed=True, flag_untestable=True)])
    assert bh.call_rows(rows, {}, CFG, {})["call"].tolist() == ["S", "B", "N", "U"]


def test_line_shift_uses_screen_median():
    # a screen whose singles sit at +0.3: raw DKO -0.3 is only 0.6 below m
    singles = {"L": pd.Series(np.full(100, 0.3))}
    rows = pd.DataFrame([_obs(gi_strict=True, gi_relaxed=True, lfc1=0.3, lfc2=0.3, lfc12=-0.3),
                         _obs(gi_strict=True, gi_relaxed=True, lfc1=0.3, lfc2=0.3, lfc12=0.0)])
    assert bh.call_rows(rows, singles, CFG, {})["call"].tolist() == ["S", "B"]


def test_collapse_rules():
    assert bh.collapse_dups(["S", "S"]) == "S"
    assert bh.collapse_dups(["S", "N"]) == "B"
    assert bh.collapse_dups(["U", "N"]) == "N"
    assert bh.collapse_dups(["U", "X"]) == "U"
    assert bh.collapse_dups(["X"]) == "X"
    assert bh.collapse_domains(["S", "S", "N"]) == "S"
    assert bh.collapse_domains(["S", "N"]) == "B"      # 50% is not > 50%
    assert bh.collapse_domains(["N", "B"]) == "B"
    assert bh.collapse_domains(["N", "N", "U"]) == "N"
    assert bh.collapse_domains(["U", "U"]) == "U"


def test_merge_line_rules():
    assert bh.merge_line(["S", "B"]) == "S"
    assert bh.merge_line(["S", "N"]) == "D"
    assert bh.merge_line(["N", "N", "U"]) == "N"
    assert bh.merge_line(["N", "B"]) == "D"
    assert bh.merge_line(["B"]) == "B"
    assert bh.merge_line(["U", "X"]) == "U"


def test_norm_line_id_strips_nbsp():
    assert hp.norm_line_id("\xa0CVCL_0313") == "CVCL_0313"
    assert hp.norm_line_id(" CVCL_1119 ") == "CVCL_1119"


def test_config_thresholds_match_rule_strings():
    """D8: numeric thresholds (config section 8) must be the values written in the section-3 rule strings."""
    st, thr = CFG["studies"], CFG["thresholds"]

    def has(rule: str, v) -> bool:
        return any(tok in rule for tok in {f"< {v}", f"<= {v}", f"> {v}", f"== {v}"})

    checks = [
        ("dede2020", "GI_strict", thr["dede2020"]["gi_strict_zdlfc_lt"]),
        ("dede2020", "GI_relaxed", thr["dede2020"]["gi_relaxed_zdlfc_lt"]),
        ("thompson2021", "GI_strict", thr["thompson2021"]["gi_strict_fdr_lt"]),
        ("thompson2021", "GI_relaxed", thr["thompson2021"]["gi_relaxed_fdr_lt"]),
        ("thompson2021", "GI_relaxed", thr["thompson2021"]["gi_relaxed_d14_fdr_lt"]),
        ("parrish2021", "GI_strict", thr["parrish2021"]["gi_strict_zdlfc_lt"]),
        ("parrish2021", "GI_relaxed", thr["parrish2021"]["gi_relaxed_zdlfc_lt"]),
        ("ito2021", "GI_strict", thr["ito2021"]["gi_strict_fdr_lt"]),
        ("ito2021", "GI_relaxed", thr["ito2021"]["gi_relaxed_fdr_lt"]),
        ("chymera2020", "GI_strict", thr["chymera2020"]["effect_size_lt"]),
        ("chymera2020", "GI_relaxed", thr["chymera2020"]["effect_size_lt"]),
        ("in4mer2024", "GI_strict", thr["in4mer2024"]["gi_strict_zdlfc_lt"]),
        ("in4mer2024", "GI_relaxed", thr["in4mer2024"]["gi_relaxed_zdlfc_lt"]),
        ("klingbeil2024", "GI_strict", thr["klingbeil2024"]["gi_strict_fdr_lt"]),
        ("klingbeil2024", "GI_relaxed", thr["klingbeil2024"]["gi_relaxed_fdr_lt"]),
        ("klingbeil2024", "domain_collapse", thr["klingbeil2024"]["domain_sl_fraction_gt"]),
        ("harle2025", "GI_strict", thr["harle_chan"]["gi_strict_gi_lt"]),
        ("harle2025", "GI_strict", thr["harle_chan"]["gi_strict_fdr_lt"]),
        ("harle2025", "GI_strict", thr["harle_chan"]["bassik_hit_eq"]),
        ("harle2025", "GI_relaxed", thr["harle_chan"]["gi_relaxed_gi_lt"]),
        ("harle2025", "GI_relaxed", thr["harle_chan"]["gi_relaxed_fdr_lt"]),
        ("harle2025", "GI_relaxed", thr["harle_chan"]["gi_relaxed_direction_gi_lt"]),
    ]
    bad = [(s, k, v) for s, k, v in checks if not has(str(st[s][k]), v)]
    assert not bad, bad
    assert hp.THR == thr
