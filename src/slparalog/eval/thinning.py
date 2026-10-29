"""Panel-thinning test and masking augmentation (DESIGN v3 §4 "Panel-density robustness"; D17(d)).

Feature "views" are the same model features of the same pairs (same row order) computed on the natural
panel and on panels thinned by ``--panel-dmin c`` (species with d_s < c removed). Conditions, for each c:

* ``N->N``        train natural, evaluate natural (the incumbent, reproduced)
* ``N->T<c>``     train natural, evaluate on the thinned-panel features (the thinning test)
* ``T<c>->T<c>``  train and evaluate on the thinned features
* ``A<c>->N``     masking augmentation: train on natural + c-masked copies of every training pair, evaluate natural
* ``A<c>->T<c>``  the same augmented model, evaluate on the thinned features

Every condition uses the same fold columns as the incumbent: all copies of a pair sit in the pair's fold
(the augmentation is applied inside each training fold; test pairs are never in training), evaluation uses
one row per pair, and the statistic is T_pf of ``slparalog.eval.perfold`` (per-fold paired ΔAP vs B1).

Augmentation details (DESIGN §4): species-balanced weights are computed on pairs and **spread** over the
pair's unique copies, so a pair's weights sum to its pair weight (identical copies are deduplicated);
``min_samples_leaf`` is counted in pairs, i.e. multiplied by rows/pairs of the augmented training set
(early stopping stays off, as for every GBDT of the harness).
"""

from __future__ import annotations

import math

import numpy as np

from slparalog.eval.bootstrap import ci
from slparalog.eval.cv import fold_col
from slparalog.eval.stats import species_balanced_weights


def row_equal(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Per row: all features equal, NaN equal to NaN."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return ((a == b) | (np.isnan(a) & np.isnan(b))).all(axis=1)


def augment(views: list[np.ndarray], y: np.ndarray, w: np.ndarray):
    """Stack the views of every training pair, drop copies identical to an earlier copy of the same pair and
    spread the pair weight ``w`` over the unique copies.

    Returns (X_aug, y_aug, w_aug, n_pairs, copies_per_pair). The first view is always kept. Within a pair the
    weights of the unique copies sum to ``w`` (so the total weight is unchanged and equals the pair-level
    weighting)."""
    n = len(y)
    views = [np.asarray(v, dtype=float) for v in views]
    keep = [np.ones(n, bool)]
    for k in range(1, len(views)):
        dup = np.zeros(n, bool)
        for j in range(k):
            dup |= row_equal(views[k], views[j])
        keep.append(~dup)
    n_unique = np.sum(keep, axis=0)
    X = np.vstack([v[m] for v, m in zip(views, keep)])
    yy = np.concatenate([np.asarray(y)[m] for m in keep])
    ww = np.concatenate([(np.asarray(w, dtype=float) / n_unique)[m] for m in keep])
    return X, yy, ww, n, n_unique


def leaf_in_pairs(min_samples_leaf: int, n_rows: int, n_pairs: int) -> int:
    """``min_samples_leaf`` counted in pairs -> rows (ceil, never below the pair-level value)."""
    return max(int(min_samples_leaf), int(math.ceil(min_samples_leaf * n_rows / n_pairs)))


def fit_condition(make_model, params: dict, views: list[np.ndarray], y, species, seed: int):
    """Fit ``make_model(params, seed)`` on the pair-weighted, deduplicated stack of ``views`` (one view = plain
    training on that view). Species-balanced pair weights (as the harness) are spread over the copies."""
    w = species_balanced_weights(species, y)
    X, yy, ww, n_pairs, _ = augment(views, y, w)
    p = dict(params)
    if len(views) > 1 and "min_samples_leaf" in p:
        p["min_samples_leaf"] = leaf_in_pairs(p["min_samples_leaf"], len(yy), n_pairs)
    m = make_model(p, seed)
    m.fit(X, yy, sample_weight=ww)
    return m


def _one_seed(views, conditions, y, species, folds, make_model, params, seed):
    """OOF scores of one seed for every condition: {cond: (n,)}. ``views`` {name: (n, p)}; ``conditions``
    {cond: (train_view_names, eval_view_name)}. Training fits are shared between conditions with the same
    training views."""
    n = len(y)
    out = {c: np.full(n, np.nan) for c in conditions}
    train_keys = {tuple(tv) for tv, _ in conditions.values()}
    for k in np.unique(folds):
        tr, te = folds != k, folds == k
        if te.sum() == 0 or tr.sum() == 0:
            continue
        for tk in train_keys:
            m = fit_condition(make_model, params, [views[v][tr] for v in tk], y[tr], species[tr], seed)
            for c, (tv, ev) in conditions.items():
                if tuple(tv) == tk:
                    out[c][te] = m.predict(views[ev][te])
    return out


def oof_conditions(df, views: dict[str, np.ndarray], conditions: dict[str, tuple[list[str], str]],
                   make_model, params: dict, seeds: list[int], n_jobs: int = 1) -> dict[str, np.ndarray]:
    """{cond: (n_seeds, n) per-seed OOF score matrix} (one score per pair per seed)."""
    y = df["label"].to_numpy()
    species = df["species"].to_numpy()
    args = [(views, conditions, y, species, df[fold_col(s)].to_numpy(), make_model, params, s) for s in seeds]
    if n_jobs > 1 and len(args) > 1:
        from joblib import Parallel, delayed

        rows = Parallel(n_jobs=min(n_jobs, len(args)), backend="loky")(delayed(_one_seed)(*a) for a in args)
    else:
        rows = [_one_seed(*a) for a in args]
    return {c: np.vstack([r[c] for r in rows]) for c in conditions}


def condition_set(cs: list[float]) -> dict[str, tuple[list[str], str]]:
    """The declared conditions for the ladder ``cs`` ({name: (training views, evaluation view)})."""
    out = {"N->N": (["N"], "N")}
    for c in cs:
        t = f"T{c:g}"
        out[f"N->{t}"] = (["N"], t)
        out[f"{t}->{t}"] = ([t], t)
        out[f"A{c:g}->N"] = (["N", t], "N")
        out[f"A{c:g}->{t}"] = (["N", t], t)
    return out


def decide(T: dict, boot: dict, c: float, level: float = 0.95, recover_frac: float = 0.5,
           max_natural_loss: float = 0.01) -> dict:
    """Adoption rule for masking augmentation at ``c`` (pre-declared; DESIGN §4 + PI 2026-10-28).

    ``T`` {cond: T_pf point}; ``boot`` {cond: (B,) bootstrap T_pf, same resamples}. Adopt iff
      (1) drop = T[N->N] - T[N->T<c>] > 0 and the paired percentile CI of the drop has lo > 0,
      (2) (T[A<c>->T<c>] - T[N->T<c>]) >= recover_frac * drop,
      (3) T[N->N] - T[A<c>->N] <= max_natural_loss.
    """
    t = f"T{c:g}"
    nn, nt, an, at = "N->N", f"N->{t}", f"A{c:g}->N", f"A{c:g}->{t}"
    drop = T[nn] - T[nt]
    drop_ci = ci(boot[nn] - boot[nt], level)
    recovered = T[at] - T[nt]
    nat_loss = T[nn] - T[an]
    c1 = bool(drop > 0 and drop_ci["lo"] > 0)
    c2 = bool(drop > 0 and recovered >= recover_frac * drop)
    c3 = bool(nat_loss <= max_natural_loss)
    return {"c": c, "drop": float(drop), "drop_ci": drop_ci,
            "recovered": float(recovered), "recovered_fraction": float(recovered / drop) if drop > 0 else None,
            "recovered_ci": ci(boot[at] - boot[nt], level),
            "natural_loss": float(nat_loss), "natural_loss_ci": ci(boot[nn] - boot[an], level),
            "cond1_drop_ci_excludes_0": c1, "cond2_recovers_half": c2, "cond3_natural_loss_le_0.01": c3,
            "adopt": bool(c1 and c2 and c3)}


# Text logged in results/experiments.tsv (declaration row, before any score of the run is computed).
RULE_TEXT = (
    "PANEL-THINNING TEST (DESIGN §4 panel-density robustness; D17(d)), GBDT_AB_17 on feature set A+B (tb_asym_* excluded), "
    "hyper-parameters fixed, no re-tuning, seeds config/seeds.yaml, same folds. Thinned-panel Tier B tables "
    "tierB_dmin<c>.parquet, c in {0.03, 0.12, 0.35}, recomputed from the cached panel hits (species with d_s < c removed). "
    "Statistic T_pf (per-fold paired dAP vs B1, mean over species) with per-species dAP vs B1. Conditions per c: "
    "N->N (incumbent), N->T<c> (train natural, evaluate thinned), T<c>->T<c> (train and evaluate thinned), "
    "A<c>->N and A<c>->T<c> (masking augmentation: train on natural plus c-masked copies of each training pair; all copies of "
    "a pair in the pair's fold; pair weights (species-balanced, computed on pairs) spread over the pair's unique copies so each "
    "pair sums to its pair weight; identical copies deduplicated; min_samples_leaf counted in pairs (x rows/pairs); no early "
    "stopping; evaluation on one row per pair). ADOPTION RULE for masking augmentation, decided at c=0.35: adopt iff "
    "(1) T_pf(N->N) - T_pf(N->T0.35) > 0 with the 95% percentile paired-bootstrap CI of the drop excluding 0, AND "
    "(2) T_pf(A0.35->T0.35) - T_pf(N->T0.35) >= 0.5 x the drop (recovers at least half), AND "
    "(3) T_pf(N->N) - T_pf(A0.35->N) <= 0.01 (natural-panel T_pf not lowered by more than 0.01). Point estimates for (2) and (3). "
    "c = 0.03 and 0.12 are reported, not decisive."
)
