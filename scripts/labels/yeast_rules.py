"""Yeast label rules (config/labels_yeast.yaml v1; reviews/round2/data.md R2-1, R2-3, R2-9, R2-10, R2-11).

Pure, vectorised functions so that they can be unit-tested on toy data (tests/test_yeast_label_rules.py).
Thresholds are read from nowhere else: they are module constants that mirror the pre-registered config and a
test asserts that they equal the config values.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# measurement status codes
MISSING, UNTESTABLE, POS, NEG, MID = 0, 1, 2, 3, 4
STATUS_NAME = {MISSING: "MISSING", UNTESTABLE: "UNTESTABLE", POS: "POS", NEG: "NEG", MID: "MID"}

E_FLOOR = 0.5            # testable iff E > 0.5
POS_EPS, POS_P, POS_DMF = -0.08, 0.05, 0.5
NEG_EPS_STRONG, NEG_EPS_WEAK, NEG_DMF = -0.04, -0.08, 0.5
LEN_EPS = -0.12          # lenient addition: eps < -0.12 and p < 0.05
STRICT_DMF = 0.2         # strict: POS and DMF <= 0.2
DCMP_EPS = -0.35         # report-only comparable tier: eps < -0.35 and p < 0.05
MID_FRAC = 0.25          # negative needs n_MID <= floor(0.25 * n_T)
DEAN_E_FLOOR = 0.5       # Dean testable iff W_A * W_B > 0.5
DEAN_NEG_EPS = -0.05     # Dean negative iff W_AB - W_A*W_B >= -0.05
LINK_WINDOW_BP = 50_000


def measurement_calls(eps, p, dmf, E):
    """Return (primary_status, lenient_status, strict_pos, cmp_pos) arrays.

    Order: testability first (E > 0.5, E NaN -> UNTESTABLE), then MISSING if eps/p/DMF is NaN, then POS/NEG/MID.
    lenient_status replaces POS by POS or (eps < -0.12 and p < 0.05); NEG is identical (the two are disjoint).
    cmp_pos (report-only tier) = eps < -0.35 and p < 0.05 for rows with eps and p present, regardless of E.
    """
    eps = np.asarray(eps, dtype=float)
    p = np.asarray(p, dtype=float)
    dmf = np.asarray(dmf, dtype=float)
    E = np.asarray(E, dtype=float)
    with np.errstate(invalid="ignore"):
        testable = E > E_FLOOR
        missing = np.isnan(eps) | np.isnan(p) | np.isnan(dmf)
        pos = (eps < POS_EPS) & (p < POS_P) & (dmf <= POS_DMF)
        neg = ((eps >= NEG_EPS_STRONG) | ((p >= POS_P) & (eps >= NEG_EPS_WEAK))) & (dmf > NEG_DMF)
        len_extra = (eps < LEN_EPS) & (p < POS_P)
        strict = pos & (dmf <= STRICT_DMF)
        cmp_pos = (eps < DCMP_EPS) & (p < POS_P) & ~(np.isnan(eps) | np.isnan(p))
    st = np.full(eps.shape, MID, dtype=np.int8)
    st[neg] = NEG
    st[pos] = POS
    lst = st.copy()
    lst[len_extra] = POS
    for arr in (st, lst):
        arr[missing] = MISSING
        arr[~testable] = UNTESTABLE
    ok = testable & ~missing
    return st, lst, strict & ok, cmp_pos


def aggregate(n_pos, n_neg, n_mid):
    """R2-3 count-neutral pair rule. Returns (call array of str, discordant bool array).

    untestable if n_T == 0; pos if n_POS >= 1 and n_POS >= n_NEG; amb+discordant if n_POS >= 1 and n_NEG > n_POS;
    neg if n_POS == 0 and n_NEG >= 1 and n_MID <= floor(0.25 n_T); amb otherwise."""
    n_pos = np.asarray(n_pos, dtype=np.int64)
    n_neg = np.asarray(n_neg, dtype=np.int64)
    n_mid = np.asarray(n_mid, dtype=np.int64)
    n_t = n_pos + n_neg + n_mid
    out = np.full(n_t.shape, "amb", dtype=object)
    discordant = (n_pos >= 1) & (n_neg > n_pos)
    pos = (n_pos >= 1) & (n_pos >= n_neg)
    neg = (n_pos == 0) & (n_neg >= 1) & (n_mid <= np.floor(MID_FRAC * n_t))
    out[pos] = "pos"
    out[neg] = "neg"
    out[n_t == 0] = "untestable"
    return out, discordant


def strict_sga(primary_call, n_spos, n_neg):
    """Strict SGA call: pos iff primary pos and n_SPOS >= 1 and n_SPOS >= n_NEG; neg iff primary neg."""
    primary_call = np.asarray(primary_call, dtype=object)
    n_spos = np.asarray(n_spos)
    n_neg = np.asarray(n_neg)
    out = primary_call.copy()
    out[(primary_call == "pos") & ~((n_spos >= 1) & (n_spos >= n_neg))] = "amb"
    return out


def dean_call(w1, w2, w12):
    """R2-9 point rule. Returns array of pos/neg/amb/untestable."""
    w1 = np.asarray(w1, dtype=float)
    w2 = np.asarray(w2, dtype=float)
    w12 = np.asarray(w12, dtype=float)
    exp = w1 * w2
    out = np.full(w1.shape, "untestable", dtype=object)
    with np.errstate(invalid="ignore"):
        testable = (exp > DEAN_E_FLOOR) & ~np.isnan(w12)
        pos = testable & (w12 == 0)
        neg = testable & ~pos & (w12 - exp >= DEAN_NEG_EPS)
    out[testable] = "amb"
    out[pos] = "pos"
    out[neg] = "neg"
    return out


def merge_dean(sga, dean):
    """Per-tier Dean merge. sga in {pos,neg,amb,untestable,untested}; dean in {pos,neg,amb,untestable,None/''}.

    Returns (call, dean_override, dean_conflict)."""
    sga = np.asarray(sga, dtype=object)
    dean = np.asarray(dean, dtype=object)
    out = sga.copy()
    override = (dean == "pos") & (sga == "neg")
    conflict = (dean == "neg") & (sga == "pos")
    out[dean == "pos"] = "pos"
    out[(dean == "neg") & (sga != "pos")] = "neg"
    out[conflict] = "amb"
    return out, override, conflict


def strict_final(sga_strict, dean, primary_final):
    """strict pos = (SGA strict pos and Dean != neg) or Dean pos; neg = primary neg; untestable/untested kept; else amb."""
    sga_strict = np.asarray(sga_strict, dtype=object)
    dean = np.asarray(dean, dtype=object)
    primary_final = np.asarray(primary_final, dtype=object)
    out = np.full(primary_final.shape, "amb", dtype=object)
    keep = np.isin(primary_final, ["untestable", "untested", "neg"])
    out[keep] = primary_final[keep]
    pos = ((sga_strict == "pos") & (dean != "neg")) | (dean == "pos")
    out[pos] = "pos"
    return out


def biogrid_tier(primary_final, n_lt_sl):
    """primary, plus pairs with >= 1 qualifying BioGRID low-throughput SL record -> pos. Returns (call, only, conflict)."""
    primary_final = np.asarray(primary_final, dtype=object)
    has = np.asarray(n_lt_sl) > 0
    out = primary_final.copy()
    only = has & (primary_final != "pos")
    conflict = has & (primary_final == "neg")
    out[has] = "pos"
    return out, only, conflict


def linked_pairs(chrom_a, start_a, end_a, chrom_b, start_b, end_b, window=LINK_WINDOW_BP):
    """1 if same chromosome and the gap between the gene intervals (0 if overlapping) <= window; NaN if a gene
    lacks coordinates (chrom is None/NaN); else 0."""
    ca = pd.Series(chrom_a, dtype=object)
    cb = pd.Series(chrom_b, dtype=object)
    sa, ea = np.asarray(start_a, float), np.asarray(end_a, float)
    sb, eb = np.asarray(start_b, float), np.asarray(end_b, float)
    gap = np.maximum(0.0, np.maximum(sa, sb) - np.minimum(ea, eb))
    same = (ca.values == cb.values)
    out = np.where(same & (gap <= window), 1.0, 0.0)
    unknown = ca.isna().values | cb.isna().values
    out[unknown] = np.nan
    return out, np.where(same & ~unknown, gap, np.nan)
