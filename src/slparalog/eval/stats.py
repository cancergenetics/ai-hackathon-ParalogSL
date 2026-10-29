"""Secondary metrics, fast weighted AP/AUROC for the bootstrap, and species-balanced sample weights.

Point estimates of AUPRC always come from ``slparalog.eval.metrics.ap`` (sklearn
``average_precision_score``; AGENT_RULES 5). The vectorised ``ap_weighted_batch`` /
``auroc_weighted_batch`` below are used only inside the bootstrap, where each resample is a vector of
integer multiplicities; they reproduce sklearn's ``average_precision_score`` / ``roc_auc_score`` with
``sample_weight`` exactly (tie blocks handled as one threshold) — verified in
``tests/test_eval_harness.py`` and re-checked at run time on a few resamples (``check_against_sklearn``).
"""

from __future__ import annotations

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

from slparalog.eval.metrics import ap


def prevalence(y) -> float:
    y = np.asarray(y)
    return float(y.mean()) if len(y) else float("nan")


def auroc(y, s) -> float:
    y = np.asarray(y)
    if y.min() == y.max():
        return float("nan")
    return float(roc_auc_score(y, np.asarray(s, dtype=float)))


def precision_at_k(y, s, k: int) -> float:
    """Precision among the top-k scores; ties at the k-th score are resolved by averaging over the tie
    block (expected precision under random tie-breaking)."""
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    k = int(min(k, len(y)))
    if k <= 0:
        return float("nan")
    order = np.argsort(-s, kind="mergesort")
    s_sorted = s[order]
    kth = s_sorted[k - 1]
    above = s > kth
    tie = s == kth
    n_above = int(above.sum())
    n_tie_needed = k - n_above
    tp = y[above].sum() + (y[tie].mean() * n_tie_needed if n_tie_needed > 0 else 0.0)
    return float(tp / k)


def summary(y, s, ks=(10, 50)) -> dict:
    """AP (sklearn), AUROC, prevalence, lift, precision@k and counts for one score vector."""
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=float)
    n_pos = int(y.sum())
    prev = prevalence(y)
    out = {"n": int(len(y)), "n_pos": n_pos, "n_neg": int(len(y) - n_pos), "prevalence": prev}
    if n_pos == 0:
        out.update({"ap": float("nan"), "auroc": float("nan"), "lift": float("nan")})
        return out
    a = ap(y, s)
    out["ap"] = a
    out["auroc"] = auroc(y, s)
    out["lift"] = a / prev
    out["ap_above_chance"] = (a - prev) / (1 - prev) if prev < 1 else float("nan")  # all-positive subset
    for k in ks:
        out[f"p_at_{k}"] = precision_at_k(y, s, k)
    out["p_at_npos"] = precision_at_k(y, s, n_pos)
    return out


# ---------------------------------------------------------------- fast weighted AP / AUROC


def _blocks_desc(s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Order sorting ``s`` descending and start indices of tie blocks in that order."""
    order = np.argsort(-s, kind="mergesort")
    ss = s[order]
    starts = np.r_[0, np.flatnonzero(np.diff(ss) != 0) + 1]
    return order, starts


def ap_weighted_batch(y: np.ndarray, s: np.ndarray, W: np.ndarray) -> np.ndarray:
    """AP for each row of the weight matrix ``W`` (shape B x n), same definition as sklearn's
    ``average_precision_score(y, s, sample_weight=W[b])``. Rows with zero positive weight -> NaN."""
    y = np.asarray(y).astype(float)
    s = np.asarray(s, dtype=float)
    W = np.atleast_2d(np.asarray(W, dtype=float))
    order, starts = _blocks_desc(s)
    Wo = W[:, order]
    yo = y[order]
    tp_blk = np.add.reduceat(Wo * yo, starts, axis=1)
    fp_blk = np.add.reduceat(Wo * (1.0 - yo), starts, axis=1)
    tps = np.cumsum(tp_blk, axis=1)
    fps = np.cumsum(fp_blk, axis=1)
    P = tps[:, -1:]
    denom = tps + fps
    with np.errstate(invalid="ignore", divide="ignore"):
        prec = np.where(denom > 0, tps / np.where(denom > 0, denom, 1.0), 0.0)
        rec = np.where(P > 0, tps / np.where(P > 0, P, 1.0), 0.0)
    drec = np.diff(np.concatenate([np.zeros((rec.shape[0], 1)), rec], axis=1), axis=1)
    out = (drec * prec).sum(axis=1)
    out[P[:, 0] <= 0] = np.nan
    return out


def auroc_weighted_batch(y: np.ndarray, s: np.ndarray, W: np.ndarray) -> np.ndarray:
    """AUROC for each row of ``W`` (ties count 1/2), as sklearn ``roc_auc_score`` with sample_weight."""
    y = np.asarray(y).astype(float)
    s = np.asarray(s, dtype=float)
    W = np.atleast_2d(np.asarray(W, dtype=float))
    order = np.argsort(s, kind="mergesort")
    ss = s[order]
    starts = np.r_[0, np.flatnonzero(np.diff(ss) != 0) + 1]
    Wo = W[:, order]
    yo = y[order]
    pos_blk = np.add.reduceat(Wo * yo, starts, axis=1)
    neg_blk = np.add.reduceat(Wo * (1.0 - yo), starts, axis=1)
    neg_below = np.cumsum(neg_blk, axis=1) - neg_blk
    num = (pos_blk * (neg_below + 0.5 * neg_blk)).sum(axis=1)
    P = pos_blk.sum(axis=1)
    N = neg_blk.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / (P * N)
    out[(P <= 0) | (N <= 0)] = np.nan
    return out


def check_against_sklearn(y, s, W, n_check: int = 3, tol: float = 1e-9) -> None:
    """Assert that the fast batch AP/AUROC agree with sklearn on the first ``n_check`` rows of W."""
    W = np.atleast_2d(W)
    fast_ap = ap_weighted_batch(y, s, W[:n_check])
    fast_auc = auroc_weighted_batch(y, s, W[:n_check])
    for i in range(min(n_check, W.shape[0])):
        w = W[i]
        if (w * np.asarray(y)).sum() == 0:
            continue
        ref = average_precision_score(y, s, sample_weight=w)
        if abs(ref - fast_ap[i]) > tol:
            raise AssertionError(f"fast AP {fast_ap[i]} != sklearn {ref}")
        if ((1 - np.asarray(y)) * w).sum() > 0:
            ref2 = roc_auc_score(y, s, sample_weight=w)
            if abs(ref2 - fast_auc[i]) > tol:
                raise AssertionError(f"fast AUROC {fast_auc[i]} != sklearn {ref2}")


# ---------------------------------------------------------------- sample weights


def species_balanced_weights(species, y, target_prevalence: float | None = None) -> np.ndarray:
    """Sample weights giving each species equal total weight **and** equal weighted prevalence
    (DESIGN v3 §5).

    With S species and N rows, every species gets total weight N/S; within species s, positives share
    (N/S)·π* and negatives (N/S)·(1−π*), where π* = ``target_prevalence`` or, by default, the mean of
    the species' prevalences. The weights sum to N (mean weight 1). With one species and the default
    target, all weights are 1.
    """
    species = np.asarray(species)
    y = np.asarray(y).astype(int)
    sps = list(dict.fromkeys(species.tolist()))
    prevs = {s: y[species == s].mean() for s in sps}
    for s, p in prevs.items():
        if p <= 0 or p >= 1:
            raise ValueError(f"species {s} has prevalence {p}; weights undefined")
    pi = float(np.mean(list(prevs.values()))) if target_prevalence is None else float(target_prevalence)
    N = len(y)
    S = len(sps)
    w = np.empty(N, dtype=float)
    for s in sps:
        m = species == s
        npos = (y[m] == 1).sum()
        nneg = (y[m] == 0).sum()
        w[m & (y == 1)] = (N / S) * pi / npos
        w[m & (y == 0)] = (N / S) * (1 - pi) / nneg
    return w
