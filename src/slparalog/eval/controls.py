"""Controls and pre-declared strata for the development harness (DESIGN v3 §6.5; RESEARCH_TREE D8, D10).

- shuffled-label negative control: labels permuted within species, the full CV re-run; AP should be
  close to prevalence.
- species-identification control: grouped-CV AUROC of predicting species from the model features.
- source-library control (human): grouped-CV AUROC of predicting a heavily screened library
  (n_studies_testable >= threshold) from the model features.
- identity-stratified comparison (helper for Tier C/D): model vs B1 within B1-score quantile strata.
- subset comparison: point AP (sklearn) and bootstrap ΔAP CI vs a baseline on any subset of pairs.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd

from slparalog.eval.bootstrap import GroupBootstrap, boot_metric, ci
from slparalog.eval.cv import oof_scores
from slparalog.eval.stats import auroc, summary


def permute_within_species(df: pd.DataFrame, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    y = df["label"].to_numpy().copy()
    for _, idx in df.groupby("species").indices.items():
        y[idx] = y[idx][rng.permutation(len(idx))]
    return y


def shuffled_label_control(df: pd.DataFrame, X: np.ndarray, make_model: Callable[[int], object],
                           seeds: list[int], species=("human", "yeast")) -> dict:
    """One permutation per seed (that seed's folds); AP of OOF scores against the permuted labels."""
    rows = []
    for s in seeds:
        yp = permute_within_species(df, s)
        sc, _ = oof_scores(df, X, make_model, [s], target=yp)
        for sp in species:
            m = (df["species"] == sp).to_numpy()
            r = summary(yp[m], sc[m])
            rows.append({"seed": s, "species": sp, "ap": r["ap"], "prevalence": r["prevalence"], "auroc": r["auroc"]})
    t = pd.DataFrame(rows)
    out = {"per_permutation": rows}
    for sp, d in t.groupby("species"):
        out[sp] = {"ap_mean": float(d["ap"].mean()), "ap_sd": float(d["ap"].std(ddof=1)),
                   "auroc_mean": float(d["auroc"].mean()), "prevalence": float(d["prevalence"].iloc[0]),
                   "ap_over_prevalence": float(d["ap"].mean() / d["prevalence"].iloc[0])}
    return out


def binary_target_control(df: pd.DataFrame, X: np.ndarray, target: np.ndarray, make_model, seeds: list[int],
                          bs: GroupBootstrap, level: float, mask: np.ndarray | None = None) -> dict:
    """Grouped-CV AUROC (and AP) of predicting a binary ``target`` from features (unweighted fits)."""
    m = np.ones(len(df), bool) if mask is None else np.asarray(mask, bool)
    sub = df.loc[m].reset_index(drop=True)
    t = np.asarray(target)[m].astype(int)
    sc, _ = oof_scores(sub, X[m], make_model, seeds, target=t, weights=None)
    boot = boot_metric(bs, sub["group"].to_numpy(), t, {"s": sc}, metric="auroc")["s"]
    s = summary(t, sc)
    return {"n": int(len(t)), "n_target_pos": int(t.sum()), "target_prevalence": float(t.mean()),
            "auroc": auroc(t, sc), "auroc_ci": ci(boot, level), "ap": s["ap"]}


def subset_compare(df: pd.DataFrame, scores: pd.DataFrame, mask: np.ndarray, methods: list[str], base: str,
                   bs: GroupBootstrap, level: float, min_pos_for_ci: int = 10, with_ci: list[str] | None = None) -> dict:
    """Point AP for ``methods`` on the subset; ΔAP vs ``base`` with bootstrap CI for ``with_ci`` methods."""
    m = np.asarray(mask, bool)
    sub = df.loc[m]
    y = sub["label"].to_numpy()
    out = {"n": int(m.sum()), "n_pos": int(y.sum()), "n_neg": int(len(y) - y.sum()),
           "prevalence": float(y.mean()) if len(y) else float("nan"), "ap": {}, "delta_vs_" + base: {}}
    if y.sum() == 0:
        return out
    for meth in methods:
        out["ap"][meth] = summary(y, scores.loc[m, meth].to_numpy())["ap"]
    if with_ci and y.sum() >= min_pos_for_ci:
        cols = list(dict.fromkeys(with_ci + [base]))
        boot = boot_metric(bs, sub["group"].to_numpy(), y, {c: scores.loc[m, c].to_numpy() for c in cols})
        for meth in with_ci:
            d = boot[meth] - boot[base]
            out["delta_vs_" + base][meth] = {"point": out["ap"][meth] - out["ap"][base], **ci(d, level)}
    elif with_ci:
        out["note"] = f"fewer than {min_pos_for_ci} positives: no CI"
    return out


def quantile_strata(x: pd.Series, q: int) -> pd.Series:
    """Quantile stratum index 0..q-1 (rank-based so ties do not create empty strata)."""
    r = x.rank(method="first")
    return pd.qcut(r, q, labels=False)


def strat_avg_ap(y: np.ndarray, s: np.ndarray, nst: np.ndarray, human: np.ndarray) -> dict:
    """Human AP within n_testable_studies 1 / 2 / >=3 and their mean over strata with a positive (D8(b), D12(b)).

    Shared by the development harness (pooled OOF scores) and scripts/eval/final_test.py (test scores).
    """
    out = {}
    for lab, m in (("1", nst == 1), ("2", nst == 2), (">=3", nst >= 3)):
        mm = human & m & np.isfinite(s)
        out[lab] = summary(y[mm], s[mm])["ap"] if y[mm].sum() > 0 else float("nan")
    vals = [v for v in out.values() if np.isfinite(v)]
    out["mean"] = float(np.mean(vals)) if vals else float("nan")
    return out


def end_to_end_ap(labelled: pd.DataFrame, scores: pd.DataFrame, methods: list[str], species) -> dict:
    """End-to-end AP (DESIGN 6.4): all labelled pairs of the tier; pairs the caller missed get the minimum
    score minus one (one tie block) for every method. ``scores`` is indexed by pair_id (called pairs only).
    Shared by dev_report.py and final_test.py."""
    out = {}
    for sp in species:
        L = labelled.loc[labelled["species"] == sp]
        out[sp] = {"n": int(len(L)), "n_pos": int(L["label"].sum()), "n_not_called": int((~L["called"]).sum()),
                   "n_pos_not_called": int(L.loc[~L["called"], "label"].sum()), "ap": {}}
        for meth in methods:
            s = L["pair_id"].map(scores[meth]).to_numpy(dtype=float)
            s = np.where(np.isfinite(s), s, np.nanmin(s) - 1.0)
            out[sp]["ap"][meth] = summary(L["label"].to_numpy(), s)["ap"]
    return out
