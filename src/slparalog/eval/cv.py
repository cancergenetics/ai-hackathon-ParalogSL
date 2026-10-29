"""Grouped cross-validation on the stored fold columns (``fold_s<seed>`` for each seed in
``config/seeds.yaml``; DESIGN v3 §6.1). Folds were assigned per split group, so a group (and every
labelled pair in it) is always wholly in train or in test.

Out-of-fold (OOF) scores are averaged over the seed repeats (statistics R2 item 15). LOSO uses the same
folds: train on species X in folds ≠ k, test on species Y in fold k (statistics R2 item 14), so
orthologous groups that hold pairs of both species never cross train/test.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd

from slparalog.eval.stats import species_balanced_weights

# Number of worker processes for the seed repeats (set by the calling script; 1 = in-process).
N_JOBS = 1


def fold_col(seed: int) -> str:
    return f"fold_s{seed}"


def check_folds(df: pd.DataFrame, seeds: list[int], group_col: str = "group") -> None:
    """Every seed has a fold column and each group sits in exactly one fold."""
    for s in seeds:
        c = fold_col(s)
        if c not in df.columns:
            raise KeyError(f"missing fold column {c}")
        if (df.groupby(group_col)[c].nunique() > 1).any():
            raise ValueError(f"{c}: a group spans more than one fold")


def _fit_one_seed(X, y, species, folds, make_model, seed, tr_all, te_all, weights) -> np.ndarray:
    """OOF predictions of one seed repeat (module-level so joblib workers can run it)."""
    n = len(y)
    row = np.full(n, np.nan)
    for k in np.unique(folds):
        tr = tr_all & (folds != k)
        te = te_all & (folds == k)
        if te.sum() == 0 or tr.sum() == 0:
            continue
        w = species_balanced_weights(species[tr], y[tr]) if weights == "species_balanced" else None
        m = make_model(seed)
        m.fit(X[tr], y[tr], sample_weight=w)
        row[te] = m.predict(X[te])
    return row


def _prep(df, X, seeds, target, train_mask, test_mask):
    n = len(df)
    y = df["label"].to_numpy() if target is None else np.asarray(target)
    tr_all = np.ones(n, bool) if train_mask is None else np.asarray(train_mask, bool)
    te_all = np.ones(n, bool) if test_mask is None else np.asarray(test_mask, bool)
    species = df["species"].to_numpy()
    folds = {s: df[fold_col(s)].to_numpy() for s in seeds}
    return np.asarray(X, dtype=float), y, species, folds, tr_all, te_all


def oof_scores(
    df: pd.DataFrame,
    X: np.ndarray,
    make_model: Callable[[int], object],
    seeds: list[int],
    target: np.ndarray | None = None,
    weights: str | None = "species_balanced",
    train_mask: np.ndarray | None = None,
    test_mask: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """OOF scores for every row in ``test_mask`` (NaN elsewhere), averaged over seeds.

    ``make_model(seed)`` returns an object with ``fit(X, y, sample_weight)`` and ``predict(X)``.
    ``weights`` = "species_balanced" (equal species weight and weighted prevalence, computed on the
    training rows of each fold) or None (unweighted). Returns (mean over seeds, per-seed matrix).
    """
    res = oof_many(df, {"_": dict(make_model=make_model, X=X, target=target, weights=weights,
                                  train_mask=train_mask, test_mask=test_mask)}, seeds)
    return res["_"]


def oof_many(df: pd.DataFrame, specs: dict[str, dict], seeds: list[int], n_jobs: int | None = None
             ) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Run several OOF specifications, parallel over (spec, seed) with ``n_jobs`` (default N_JOBS) workers.

    ``specs[name]`` holds make_model, X and optionally target, weights ("species_balanced"), train_mask,
    test_mask (as in ``oof_scores``). Returns {name: (mean over seeds, per-seed matrix)}.
    """
    n_jobs = N_JOBS if n_jobs is None else n_jobs
    jobs = []
    for name, sp in specs.items():
        X, y, species, folds, tr_all, te_all = _prep(df, sp["X"], seeds, sp.get("target"),
                                                     sp.get("train_mask"), sp.get("test_mask"))
        w = sp.get("weights", "species_balanced")
        for s in seeds:
            jobs.append((name, (X, y, species, folds[s], sp["make_model"], s, tr_all, te_all, w)))
    if n_jobs > 1 and len(jobs) > 1:
        from joblib import Parallel, delayed

        rows = Parallel(n_jobs=min(n_jobs, len(jobs)), backend="loky")(
            delayed(_fit_one_seed)(*a) for _, a in jobs)
    else:
        rows = [_fit_one_seed(*a) for _, a in jobs]
    out: dict[str, list] = {}
    for (name, _), row in zip(jobs, rows):
        out.setdefault(name, []).append(row)
    return {k: (np.vstack(v).mean(axis=0), np.vstack(v)) for k, v in out.items()}
