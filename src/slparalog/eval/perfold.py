"""Per-fold paired ΔAP (D12(a); independent review M1).

For every seed repeat s and fold k, the model fitted without fold k and B1 are scored on exactly the pairs
of fold k (one species at a time); ΔAP_{s,k,sp} = AP(model) − AP(B1). The per-species statistic is the
mean of ΔAP over the seed × fold test sets that hold ≥1 positive of that species; T_pf is the mean over
species. Unlike pooled OOF, each AP is computed within one fitted model's scores, as at the final test.

Uncertainty: the same stratified group bootstrap as everywhere (``GroupBootstrap``): in each resample the
fold-level APs are recomputed from the group multiplicities and averaged as above, so Δ is paired and the
SE reflects group sampling, not fold-to-fold scatter (which is anti-conservative, statistics R2 item 15).
Fold test sets whose resample has no positive weight are skipped for that resample.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from slparalog.eval.bootstrap import GroupBootstrap, boot_metric, ci
from slparalog.eval.cv import fold_col
from slparalog.eval.metrics import ap


def point(df: pd.DataFrame, per_seed: dict[str, np.ndarray], base: np.ndarray, seeds: list[int],
          species=("human", "yeast")) -> dict[str, dict]:
    """{name: {per_species: {sp: {delta, ap_model, ap_base, n_folds}}, T}} from sklearn AP per fold."""
    y = df["label"].to_numpy()
    sp_arr = df["species"].to_numpy()
    out = {}
    for name, mat in per_seed.items():
        res = {"per_species": {}}
        for sp in species:
            d, am, ab = [], [], []
            for i, s in enumerate(seeds):
                f = df[fold_col(s)].to_numpy()
                for k in np.unique(f):
                    m = (f == k) & (sp_arr == sp) & np.isfinite(mat[i])
                    if m.sum() == 0 or y[m].sum() == 0:
                        continue
                    a_m = ap(y[m], mat[i][m])
                    a_b = ap(y[m], base[m])
                    d.append(a_m - a_b)
                    am.append(a_m)
                    ab.append(a_b)
            if d:
                res["per_species"][sp] = {"delta": float(np.mean(d)), "ap_model": float(np.mean(am)),
                                          "ap_base": float(np.mean(ab)), "n_folds": len(d)}
        if len(res["per_species"]) == len(species):
            res["T"] = float(np.mean([res["per_species"][sp]["delta"] for sp in species]))
        out[name] = res
    return out


def _boot_seed(bs, df_small, mats_i: dict[str, np.ndarray], base, seed, species):
    y = df_small["label"].to_numpy()
    sp_arr = df_small["species"].to_numpy()
    g = df_small["group"].to_numpy()
    f = df_small[fold_col(seed)].to_numpy()
    sums = {sp: {n: np.zeros(bs.B) for n in mats_i} for sp in species}
    cnts = {sp: {n: np.zeros(bs.B) for n in mats_i} for sp in species}
    for sp in species:
        for k in np.unique(f):
            m = (f == k) & (sp_arr == sp)
            if m.sum() == 0 or y[m].sum() == 0:
                continue
            names = [n for n in mats_i if np.isfinite(mats_i[n][m]).all()]
            if not names:
                continue
            scores = {n: mats_i[n][m] for n in names}
            scores["__base"] = base[m]
            b = boot_metric(bs, g[m], y[m], scores, check=False)
            for n in names:
                d = b[n] - b["__base"]
                ok = np.isfinite(d)
                sums[sp][n][ok] += d[ok]
                cnts[sp][n][ok] += 1
    return sums, cnts


def bootstrap(bs: GroupBootstrap, df: pd.DataFrame, per_seed: dict[str, np.ndarray], base: np.ndarray,
              seeds: list[int], species=("human", "yeast"), n_jobs: int = 1) -> dict[str, dict[str, np.ndarray]]:
    """Per-resample per-fold mean ΔAP: {name: {sp: (B,), 'T': (B,)}} (NaN where no fold had positives)."""
    df_small = df[["species", "group", "label"] + [fold_col(s) for s in seeds]].reset_index(drop=True)
    args = [(bs, df_small, {n: m[i] for n, m in per_seed.items()}, np.asarray(base, float), s, species)
            for i, s in enumerate(seeds)]
    if n_jobs > 1:
        from joblib import Parallel, delayed

        parts = Parallel(n_jobs=min(n_jobs, len(args)), backend="loky")(delayed(_boot_seed)(*a) for a in args)
    else:
        parts = [_boot_seed(*a) for a in args]
    out: dict[str, dict[str, np.ndarray]] = {}
    for n in per_seed:
        out[n] = {}
        for sp in species:
            s_ = sum(p[0][sp][n] for p in parts)
            c_ = sum(p[1][sp][n] for p in parts)
            with np.errstate(invalid="ignore", divide="ignore"):
                out[n][sp] = np.where(c_ > 0, s_ / np.where(c_ > 0, c_, 1), np.nan)
        out[n]["T"] = np.mean(np.vstack([out[n][sp] for sp in species]), axis=0)
    return out


def summarize(pt: dict, bt: dict, level: float) -> dict:
    """Merge point estimates with bootstrap CIs."""
    out = {}
    for n, r in pt.items():
        o = {"per_species": {}}
        for sp, v in r["per_species"].items():
            o["per_species"][sp] = {**v, "delta_ci": ci(bt[n][sp], level)}
        if "T" in r:
            o["T"] = {"point": r["T"], **ci(bt[n]["T"], level)}
        out[n] = o
    return out
