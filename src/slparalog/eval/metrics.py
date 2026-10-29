"""Pinned evaluation metric: AUPRC = sklearn.metrics.average_precision_score.

Model and baseline are always scored on exactly the same pairs; ``ap`` refuses
mismatched pair-ID sets (DESIGN.md section 6.2, gate G5).
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
from sklearn.metrics import average_precision_score


class PairSetMismatch(ValueError):
    """Model and baseline were not scored on the same set of pair IDs."""


def _check_unique(ids: Sequence, name: str) -> set:
    id_list = list(ids)
    id_set = set(id_list)
    if len(id_set) != len(id_list):
        raise PairSetMismatch(f"{name} pair IDs contain duplicates")
    return id_set


def ap(
    y_true: Iterable,
    y_score: Iterable,
    *,
    pair_ids: Sequence | None = None,
    reference_pair_ids: Sequence | None = None,
) -> float:
    """Average precision (AUPRC) via sklearn.

    If ``reference_pair_ids`` is given (the pair IDs the other side of a model-vs-baseline
    comparison was scored on), ``pair_ids`` must be given too and the two sets must be
    identical, otherwise ``PairSetMismatch`` is raised.
    """
    y = np.asarray(list(y_true))
    s = np.asarray(list(y_score), dtype=float)
    if y.shape != s.shape or y.ndim != 1:
        raise ValueError("y_true and y_score must be 1-D and the same length")
    if reference_pair_ids is not None:
        if pair_ids is None:
            raise PairSetMismatch("reference_pair_ids given without pair_ids")
        ids = _check_unique(pair_ids, "scored")
        ref = _check_unique(reference_pair_ids, "reference")
        if ids != ref:
            raise PairSetMismatch(
                f"pair-ID sets differ: {len(ids - ref)} only in scored, {len(ref - ids)} only in reference"
            )
    if pair_ids is not None and len(list(pair_ids)) != len(y):
        raise ValueError("pair_ids length must match y_true")
    if not np.isfinite(s).all():
        raise ValueError("y_score contains non-finite values")
    if not (y == 1).any():
        raise ValueError("no positives: AP undefined")
    return float(average_precision_score(y, s))


def ap_pair(
    model,
    baseline,
    *,
    id_col: str = "pair_id",
    label_col: str = "label",
    score_col: str = "score",
) -> dict:
    """AP of a model frame and a baseline frame on the same pairs (pandas DataFrames).

    Both frames must hold the same pair IDs (same set, unique) and agree on the labels.
    Returns AP for each, their difference and the pair counts.
    """
    m_ids = list(model[id_col])
    b_ids = list(baseline[id_col])
    # Raises PairSetMismatch on differing or duplicated IDs.
    ap_model = ap(model[label_col], model[score_col], pair_ids=m_ids, reference_pair_ids=b_ids)
    m = model.set_index(id_col)[label_col]
    b = baseline.set_index(id_col)
    b_aligned = b.loc[m.index]
    if not (m.to_numpy() == b_aligned[label_col].to_numpy()).all():
        raise PairSetMismatch("model and baseline frames disagree on labels for the same pair IDs")
    ap_base = ap(b_aligned[label_col], b_aligned[score_col])
    n_pos = int((m == 1).sum())
    n = int(len(m))
    return {
        "ap_model": ap_model,
        "ap_baseline": ap_base,
        "delta": ap_model - ap_base,
        "n": n,
        "n_pos": n_pos,
        "n_neg": n - n_pos,
        "prevalence": n_pos / n,
    }
