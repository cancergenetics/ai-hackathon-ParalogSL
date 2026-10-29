"""Toy-data tests with hand-derived answers for the yeast label rules (scripts/labels/yeast_rules.py).

No data on disk needed. Also checks that the rule constants equal the pre-registered config values.
"""
import importlib.util
import math
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("yeast_rules", ROOT / "scripts" / "labels" / "yeast_rules.py")
R = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(R)
CFG = yaml.safe_load((ROOT / "config" / "labels_yeast.yaml").read_text())


def test_constants_match_config():
    m = CFG["measurement"]
    assert m["pos"] == {"eps_lt": R.POS_EPS, "p_lt": R.POS_P, "dmf_le": R.POS_DMF}
    assert m["testable"] == f"E_m > {R.E_FLOOR}"
    assert "eps >= -0.04" in m["neg"] and "eps >= -0.08" in m["neg"] and "DMF > 0.5" in m["neg"]
    assert "eps < -0.12" in m["lenient_pos"] and "DMF <= 0.2" in m["strict_pos"] and "eps < -0.35" in m[[k for k in m if k.endswith("comparable_pos")][0]]
    assert "floor(0.25 * n_T)" in CFG["pair_aggregation"]["negative"] and R.MID_FRAC == 0.25
    assert CFG["dean2008"]["testable"] == "W_A * W_B > 0.5" and R.DEAN_E_FLOOR == 0.5
    assert "-0.05" in CFG["dean2008"]["negative"] and R.DEAN_NEG_EPS == -0.05
    assert CFG["linked"]["window_bp"] == R.LINK_WINDOW_BP == 50000


def test_measurement_calls_toy():
    #            eps    p     dmf   E     -> primary, lenient, strict, cmp
    rows = [
        (-0.30, 0.01, 0.40, 0.90),   # POS, POS, not strict (dmf .4), not cmp
        (-0.70, 0.01, 0.10, 0.90),   # POS, POS, strict, cmp
        (-0.15, 0.01, 0.70, 0.90),   # MID (dmf>0.5 blocks POS; eps<-0.04, p<.05 blocks NEG), lenient POS
        (-0.10, 0.20, 0.70, 0.90),   # NEG (p>=.05 and eps>=-0.08? no: eps=-0.10 < -0.08) -> MID; lenient MID
        (-0.06, 0.20, 0.70, 0.90),   # NEG (p>=.05 and eps>=-0.08)
        (-0.02, 0.01, 0.70, 0.90),   # NEG (eps>=-0.04)
        (-0.02, 0.01, 0.45, 0.90),   # MID (DMF <= 0.5 blocks NEG)
        (-0.30, 0.01, 0.10, 0.50),   # UNTESTABLE (E == 0.5 is not > 0.5)
        (-0.30, 0.01, np.nan, 0.90), # MISSING
        (-0.30, 0.01, 0.10, np.nan), # UNTESTABLE (E undefined)
        (-0.40, 0.01, np.nan, np.nan),  # UNTESTABLE, but cmp pos (eps/p present)
    ]
    eps, p, dmf, E = (np.array(c, dtype=float) for c in zip(*rows))
    st, lst, spos, cmp_pos = R.measurement_calls(eps, p, dmf, E)
    P, N, M, U, X = R.POS, R.NEG, R.MID, R.UNTESTABLE, R.MISSING
    assert st.tolist() == [P, P, M, M, N, N, M, U, X, U, U]
    assert lst.tolist() == [P, P, P, M, N, N, M, U, X, U, U]
    assert spos.tolist() == [False, True, False, False, False, False, False, False, False, False, False]
    assert cmp_pos.tolist() == [False, True, False, False, False, False, False, False, False, False, True]


def test_aggregate_count_neutral():
    #           pos neg mid -> call, discordant
    cases = [
        (1, 0, 0, "pos", False),
        (1, 1, 0, "pos", False),     # n_POS >= n_NEG
        (1, 2, 0, "amb", True),      # discordant
        (0, 3, 1, "neg", False),     # n_MID 1 <= floor(0.25*4)=1
        (0, 2, 1, "amb", False),     # n_MID 1 > floor(0.25*3)=0
        (0, 0, 2, "amb", False),
        (0, 0, 0, "untestable", False),
        (2, 1, 5, "pos", False),
        (0, 8, 2, "neg", False),     # floor(2.5)=2
    ]
    call, disc = R.aggregate([c[0] for c in cases], [c[1] for c in cases], [c[2] for c in cases])
    assert call.tolist() == [c[3] for c in cases]
    assert disc.tolist() == [c[4] for c in cases]


def test_strict_sga():
    out = R.strict_sga(["pos", "pos", "pos", "neg", "amb"], [1, 0, 1, 0, 1], [0, 0, 2, 3, 0])
    assert out.tolist() == ["pos", "amb", "amb", "neg", "amb"]


def test_dean_rule():
    w1 = [0.9, 0.9, 0.9, 0.6, 0.9, np.nan]
    w2 = [1.0, 1.0, 1.0, 0.8, 1.0, np.nan]
    w12 = [0.0, 0.88, 0.80, 0.0, np.nan, np.nan]
    # exp: .9 -> pos ; .88-.9=-.02 -> neg ; .8-.9=-.10 -> amb ; exp .48 -> untestable ; missing -> untestable
    assert R.dean_call(w1, w2, w12).tolist() == ["pos", "neg", "amb", "untestable", "untestable", "untestable"]


def test_merge_dean_and_strict_and_biogrid():
    sga = np.array(["neg", "pos", "pos", "amb", "untested", "neg", "pos", "untestable"], dtype=object)
    dean = np.array(["pos", "neg", "pos", "neg", "pos", "", "", "amb"], dtype=object)
    call, ov, cf = R.merge_dean(sga, dean)
    assert call.tolist() == ["pos", "amb", "pos", "neg", "pos", "neg", "pos", "untestable"]
    assert ov.tolist() == [True, False, False, False, False, False, False, False]
    assert cf.tolist() == [False, True, False, False, False, False, False, False]
    sga_strict = np.array(["neg", "amb", "pos", "amb", "untested", "neg", "amb", "untestable"], dtype=object)
    st = R.strict_final(sga_strict, dean, call)
    assert st.tolist() == ["pos", "amb", "pos", "neg", "pos", "neg", "amb", "untestable"]
    # nested: strict pos subset of primary pos
    assert set(np.where(st == "pos")[0]) <= set(np.where(call == "pos")[0])
    bg, only, conf = R.biogrid_tier(call, [0, 1, 0, 2, 0, 1, 0, 0])
    assert bg.tolist() == ["pos", "pos", "pos", "pos", "pos", "pos", "pos", "untestable"]
    assert only.tolist() == [False, True, False, True, False, True, False, False]
    assert conf.tolist() == [False, False, False, True, False, True, False, False]


def test_lenient_superset_on_random_measurements():
    rng = np.random.default_rng(11)
    n = 5000
    eps = rng.normal(-0.05, 0.2, n); p = rng.uniform(0, 0.2, n); dmf = rng.uniform(-0.1, 1.1, n); E = rng.uniform(0.3, 1.1, n)
    st, lst, spos, _ = R.measurement_calls(eps, p, dmf, E)
    assert np.all(lst[st == R.POS] == R.POS)
    assert np.array_equal(st == R.NEG, lst == R.NEG)
    assert np.all(st[spos] == R.POS)


def test_linked_pairs():
    linked, gap = R.linked_pairs(["I", "I", "I", "II", None], [100, 100, 100, 100, 1], [200, 200, 200, 200, 2],
                                 ["I", "I", "I", "I", "I"], [150, 50_200, 50_201, 100, 1], [300, 50_300, 50_400, 200, 2])
    assert linked[:4].tolist() == [1.0, 1.0, 0.0, 0.0] and math.isnan(linked[4])
    assert gap[0] == 0 and gap[1] == 50_000 and gap[2] == 50_001 and math.isnan(gap[3])
