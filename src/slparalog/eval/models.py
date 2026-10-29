"""Model wrappers for the development harness: standardised logistic regression and sklearn
HistGradientBoosting (early_stopping=False), both fitted with sample weights (DESIGN v3 §5).

Grids, fixed parameters and the simplicity order used by the one-SE rule are read from
``config/models.yaml``; nothing here is tuned.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from slparalog.config import CONFIG_DIR, load_yaml
from slparalog.features import registry


def load_models_config(path=None) -> dict:
    return load_yaml(path or CONFIG_DIR / "models.yaml")


def log1p_columns(feature_cols: list[str], extra: tuple[str, ...] = ("bits_pair",)) -> list[str]:
    """Columns log1p-transformed before standardisation in LR: registered integer counts/lengths/ranks
    (registry dtype int32) plus ``extra`` (heavy-tailed, non-negative)."""
    out = []
    for c in feature_cols:
        f = registry.BY_NAME.get(c)
        if (f is not None and f.dtype == "int32") or c in extra:
            out.append(c)
    return out


@dataclass
class LRModel:
    """log1p (count columns) -> median imputation -> StandardScaler -> L2 logistic regression."""

    feature_cols: list[str]
    C: float = 1.0
    max_iter: int = 5000
    seed: int = 0
    log_cols: list[str] = field(default_factory=list)

    def _transform_raw(self, X: np.ndarray) -> np.ndarray:
        X = np.array(X, dtype=float, copy=True)
        for j, c in enumerate(self.feature_cols):
            if c in self.log_cols:
                X[:, j] = np.log1p(np.clip(X[:, j], 0, None))
        return X

    def fit(self, X, y, sample_weight=None):
        Xt = self._transform_raw(X)
        self.imp_ = SimpleImputer(strategy="median", keep_empty_features=True).fit(Xt)
        Xt = self.imp_.transform(Xt)
        self.sc_ = StandardScaler().fit(Xt, sample_weight=sample_weight)
        Xt = self.sc_.transform(Xt)
        self.lr_ = LogisticRegression(C=self.C, max_iter=self.max_iter, random_state=self.seed)
        self.lr_.fit(Xt, y, sample_weight=sample_weight)
        return self

    def predict(self, X) -> np.ndarray:
        Xt = self.sc_.transform(self.imp_.transform(self._transform_raw(X)))
        return self.lr_.decision_function(Xt)

    @property
    def coef(self) -> dict:
        return dict(zip(self.feature_cols, self.lr_.coef_[0].tolist()))


@dataclass
class GBDTModel:
    """HistGradientBoosting with bin edges from the *unweighted* training quantiles.

    Deviation from plain ``HistGradientBoostingClassifier.fit(X, y, sample_weight=w)``: sklearn 1.9 computes
    *weighted* bin quantiles when ``sample_weight`` is given, which costs ~8 CPU-s per fit here (15x the
    boosting). We pre-bin with sklearn's own ``_BinMapper`` fitted without weights and pass bin indices to
    the classifier; sample weights still enter every gradient/hessian. Unweighted fits are reproduced
    exactly (test_gbdt_prebinning_matches_unweighted_sklearn); weighted fits differ from sklearn's native
    weighted binning (independent review, reviews/results/eval_harness_review.md m3: |ΔT| <= 0.004, prediction
    Spearman 0.95-0.99, far below SE(T) ~0.02).
    """

    feature_cols: list[str]
    params: dict = field(default_factory=dict)
    seed: int = 0
    prebin: bool = True

    def _bin(self, X) -> np.ndarray:
        X = np.asarray(X, dtype=float)
        if not self.prebin:
            return X
        Xb = self.bm_.transform(X).astype(float)
        Xb[np.isnan(X)] = np.nan
        return Xb

    def fit(self, X, y, sample_weight=None):
        p = dict(self.params)
        p.setdefault("early_stopping", False)
        if p["early_stopping"]:
            raise ValueError("early_stopping must be False (DESIGN §5: no row-random early stopping)")
        X = np.asarray(X, dtype=float)
        if self.prebin:
            from sklearn.ensemble._hist_gradient_boosting.binning import _BinMapper

            self.bm_ = _BinMapper(n_bins=256, random_state=self.seed).fit(X)
        self.m_ = HistGradientBoostingClassifier(random_state=self.seed, **p)
        self.m_.fit(self._bin(X), y, sample_weight=sample_weight)
        return self

    def predict(self, X) -> np.ndarray:
        return self.m_.predict_proba(self._bin(X))[:, 1]


DC_LEVELS = (1, 2, 3)  # n_testable_studies coded 1, 2, >=3


def deconfound_covariates(n_studies, species, yeast_level: int = 3) -> np.ndarray:
    """Training-only nuisance covariates for D12(c): indicators of human n_testable_studies == 1 and == 2
    (reference level >= 3 = both zero). Yeast rows get the reference level (both zero): SGA screens test
    pairs uniformly, closest to the most thoroughly screened human stratum. Returns (n, 2)."""
    lvl = np.where(np.asarray(species) == "human", np.minimum(np.nan_to_num(np.asarray(n_studies, float), nan=3), 3),
                   yeast_level).astype(int)
    return np.column_stack([(lvl == 1).astype(float), (lvl == 2).astype(float)])


def reference_covariates(n: int, ref: int) -> np.ndarray:
    """Covariates for prediction: every pair set to reference level ``ref`` in {1, 2, 3 (= >=3)}."""
    if ref not in DC_LEVELS:
        raise ValueError(f"reference level must be one of {DC_LEVELS}")
    return np.tile([float(ref == 1), float(ref == 2)], (n, 1))


@dataclass
class DeconfoundModel:
    """Wraps a model trained on [features | nuisance covariates] (last ``n_cov`` columns of X, built by
    ``deconfound_covariates``). At prediction the covariates are replaced by the fixed reference level
    ``ref`` for every pair, so scores never depend on the pair's ascertainment (D12(c)). The covariates are
    label-process metadata, never registered features, and never an input of the shipped tool."""

    base: object
    ref: int = 3
    n_cov: int = 2

    def fit(self, X, y, sample_weight=None):
        self.base.fit(X, y, sample_weight=sample_weight)
        return self

    def predict(self, X) -> np.ndarray:
        X = np.array(X, dtype=float, copy=True)
        X[:, -self.n_cov:] = reference_covariates(len(X), self.ref)
        return self.base.predict(X)


def grid(spec: dict) -> list[dict]:
    """Cartesian product of a {param: [values]} grid, in a fixed (sorted-key) order."""
    keys = sorted(spec)
    return [dict(zip(keys, vals)) for vals in itertools.product(*(spec[k] for k in keys))]


def gbdt_complexity(p: dict) -> tuple:
    """Simplicity order for the one-SE rule (smaller = simpler); see config/models.yaml."""
    return (p["max_leaf_nodes"], p["max_iter"] * p["learning_rate"], -p["min_samples_leaf"], -p["l2_regularization"])


def lr_complexity(p: dict) -> tuple:
    return (p["C"],)


def one_se_select(results: list[dict], complexity) -> dict:
    """One-SE rule: among configs with T >= T_best - SE(T_best), the simplest (lowest ``complexity``).

    Each result needs keys ``params``, ``T`` and ``T_se``. Returns {selected, best, threshold, eligible}.
    """
    best = max(results, key=lambda r: r["T"])
    thr = best["T"] - best["T_se"]
    eligible = [r for r in results if r["T"] >= thr]
    sel = min(eligible, key=lambda r: (complexity(r["params"]), -r["T"]))
    return {"selected": sel, "best": best, "threshold": thr, "n_eligible": len(eligible)}


def guarded_select(results: list[dict], complexity) -> dict:
    """D12 selection: one-SE rule (on ``T``/``T_se``, the per-fold paired statistic) among configs whose
    ``guard_pass`` is True (human stratum-averaged AP not below the incumbent's). Returns the one-SE result
    plus ``n_guard_pass``; ``selected`` is None if no config passes the guard."""
    ok = [r for r in results if r.get("guard_pass")]
    if not ok:
        return {"selected": None, "best": None, "threshold": None, "n_eligible": 0, "n_guard_pass": 0}
    out = one_se_select(ok, complexity)
    out["n_guard_pass"] = len(ok)
    return out


def _species_pct_rank(s: np.ndarray, sp: np.ndarray) -> np.ndarray:
    """Percentile rank (average ties, in (0, 1]) of ``s`` within each species code."""
    from scipy.stats import rankdata

    out = np.empty(len(s), dtype=float)
    for c in np.unique(sp):
        m = sp == c
        out[m] = rankdata(s[m], method="average") / m.sum()
    return out


@dataclass
class RankEnsembleModel:
    """Identity-anchored rank-average ensemble (D14 candidate (a)).

    X = [features | B1 score | group code | species code]. score = w * rank(model) + (1 - w) * rank(B1), with
    percentile ranks computed within the scored set per species (as a proteome would be ranked at inference).
    ``w`` is chosen **inside the training fold only**: inner grouped K-fold (GroupKFold on the group code)
    produces inner out-of-fold model scores; ranks are taken within each inner test fold and species, and
    the w in ``w_grid`` maximising the species-mean AP of the combined inner-OOF ranks is kept (ties -> the
    smaller w, i.e. more identity). The base model is then refitted on the whole training fold.
    """

    make_base: object  # callable () -> model with fit(X, y, sample_weight) / predict(X)
    w_grid: tuple = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    inner_folds: int = 3

    def fit(self, X, y, sample_weight=None):
        from sklearn.metrics import average_precision_score
        from sklearn.model_selection import GroupKFold

        X = np.asarray(X, dtype=float)
        y = np.asarray(y)
        F, b1, g, sp = X[:, :-3], X[:, -3], X[:, -2], X[:, -1]
        sw = None if sample_weight is None else np.asarray(sample_weight)
        oof = np.full(len(y), np.nan)
        rm = np.full(len(y), np.nan)
        rb = np.full(len(y), np.nan)
        for tr, te in GroupKFold(n_splits=self.inner_folds).split(F, y, g):
            m = self.make_base()
            m.fit(F[tr], y[tr], sample_weight=None if sw is None else sw[tr])
            oof[te] = m.predict(F[te])
            rm[te] = _species_pct_rank(oof[te], sp[te])
            rb[te] = _species_pct_rank(b1[te], sp[te])
        best = None
        self.w_scores_ = {}
        for w in self.w_grid:
            s = w * rm + (1 - w) * rb
            aps = [average_precision_score(y[sp == c], s[sp == c]) for c in np.unique(sp) if y[sp == c].sum() > 0]
            val = float(np.mean(aps))
            self.w_scores_[w] = val
            if best is None or val > best[1] + 1e-12:
                best = (w, val)
        self.w_ = best[0]
        self.base_ = self.make_base()
        self.base_.fit(F, y, sample_weight=sw)
        return self

    def predict(self, X) -> np.ndarray:
        X = np.asarray(X, dtype=float)
        F, b1, sp = X[:, :-3], X[:, -3], X[:, -1]
        return self.w_ * _species_pct_rank(self.base_.predict(F), sp) + (1 - self.w_) * _species_pct_rank(b1, sp)


def resolve_feature_set(spec: str, available: list[str], groups: dict[str, list[str]]) -> list[str]:
    """Columns of a feature-set spec 'tok+tok+...': a token is a registry tier letter (A, B, C, ...) or a
    Tier A feature-group name from config/models.yaml ``feature_groups`` (prefix lists; *_pct follow their base)."""
    cols: list[str] = []
    for tok in spec.split("+"):
        if tok in groups:
            pref = tuple(groups[tok])
            sel = [f for f in available if registry.BY_NAME[f].tier == "A" and f.startswith(pref)]
        else:
            sel = [f for f in available if registry.BY_NAME[f].tier == tok]
        if not sel:
            raise ValueError(f"feature-set token {tok!r} matches no available registered feature")
        cols += [f for f in sel if f not in cols]
    return [f for f in available if f in cols]  # registry order
