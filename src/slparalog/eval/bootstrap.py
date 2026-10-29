"""Stratified group bootstrap (DESIGN v3 §6.2; statistics R2 items 6, 15).

Resampling unit: the split group (labelled pairs are contracted into groups, so a group never spans
folds). Strata: species composition of the group (human-only, yeast-only, mixed) × whether the group
holds ≥1 positive (``config/success.yaml: bootstrap_strata``). Within each stratum, as many groups as
the stratum holds are drawn with replacement. A resample is a vector of group multiplicities; a pair's
weight in the resample is the multiplicity of its group. Every statistic in a resample (AP per species,
ΔAP, T = species mean of ΔAP, stratum averages, ...) is computed from the same multiplicities, so all
differences are paired. Resamples in which a species (or subset) has no positive weight give NaN for
that statistic; they are dropped from that statistic and counted.

CIs are percentile intervals at ``ci_level``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from slparalog.eval.stats import ap_weighted_batch, auroc_weighted_batch, check_against_sklearn


@dataclass
class GroupBootstrap:
    """Group multiplicity matrix for a fixed population of pairs."""

    groups: np.ndarray  # unique group IDs (G,)
    counts: np.ndarray  # (B, G) multiplicities, int16
    strata: pd.DataFrame  # per group: species composition, has_pos, stratum
    seed: int

    @property
    def B(self) -> int:
        return self.counts.shape[0]

    def pair_index(self, pair_groups) -> np.ndarray:
        idx = pd.Index(self.groups).get_indexer(np.asarray(pair_groups).astype(str))
        if (idx < 0).any():
            raise KeyError("pairs whose group is not in the bootstrap population")
        return idx


def group_strata(df: pd.DataFrame, group_col: str = "group", species_col: str = "species",
                 label_col: str = "label") -> pd.DataFrame:
    """Per-group stratum: composition in {<species>-only, mixed} × has_pos."""
    g = df.groupby(group_col).agg(
        n_species=(species_col, "nunique"),
        first_species=(species_col, "first"),
        has_pos=(label_col, "max"),
    )
    comp = np.where(g["n_species"] > 1, "mixed", g["first_species"].astype(str) + "_only")
    g["composition"] = comp
    g["has_pos"] = g["has_pos"].astype(int)
    g["stratum"] = g["composition"] + ":" + g["has_pos"].astype(str)
    return g


def make_bootstrap(df: pd.DataFrame, B: int, seed: int, group_col: str = "group") -> GroupBootstrap:
    """Draw B stratified group resamples for the population ``df`` (columns group, species, label)."""
    st = group_strata(df, group_col=group_col)
    groups = st.index.to_numpy().astype(str)
    rng = np.random.default_rng(seed)
    G = len(groups)
    counts = np.zeros((B, G), dtype=np.int64)
    # strata in sorted order so the draw sequence is deterministic given the seed
    strata_idx = st.reset_index().groupby("stratum").indices
    for name in sorted(strata_idx):
        idx = np.asarray(strata_idx[name])
        n = len(idx)
        draws = idx[rng.integers(0, n, size=(B, n))]  # (B, n) group positions
        lin = (np.arange(B)[:, None] * G + draws).ravel()
        counts += np.bincount(lin, minlength=B * G).reshape(B, G)
    return GroupBootstrap(groups=groups, counts=counts.astype(np.int16), strata=st, seed=seed)


def boot_metric(bs: GroupBootstrap, pair_groups, y, scores: dict[str, np.ndarray], metric: str = "ap",
                chunk: int = 1000, check: bool = True) -> dict[str, np.ndarray]:
    """Per-resample AP (or AUROC) for several score vectors on one subset of pairs.

    Returns {name: array(B,)} with NaN where the subset has no positive (or no negative) weight.
    """
    gidx = bs.pair_index(pair_groups)
    y = np.asarray(y).astype(int)
    fn = ap_weighted_batch if metric == "ap" else auroc_weighted_batch
    out = {k: np.empty(bs.B) for k in scores}
    for start in range(0, bs.B, chunk):
        W = bs.counts[start:start + chunk][:, gidx].astype(float)
        for k, s in scores.items():
            s = np.asarray(s, dtype=float)
            if check and start == 0:
                check_against_sklearn(y, s, W, n_check=2)
            out[k][start:start + chunk] = fn(y, s, W)
    return out


def ci(x: np.ndarray, level: float = 0.95) -> dict:
    """Percentile CI, SD and the number of finite resamples."""
    x = np.asarray(x, dtype=float)
    f = x[np.isfinite(x)]
    if len(f) == 0:
        return {"lo": float("nan"), "hi": float("nan"), "se": float("nan"), "n_ok": 0, "n_dropped": int(len(x))}
    a = (1 - level) / 2
    return {
        "lo": float(np.quantile(f, a)),
        "hi": float(np.quantile(f, 1 - a)),
        "se": float(np.std(f, ddof=1)) if len(f) > 1 else float("nan"),
        "n_ok": int(len(f)),
        "n_dropped": int(len(x) - len(f)),
    }


def species_boot(bs: GroupBootstrap, df: pd.DataFrame, score_cols: list[str], metric: str = "ap",
                 mask: np.ndarray | None = None, species=("human", "yeast"), n_jobs: int = 1
                 ) -> dict[str, dict[str, np.ndarray]]:
    """{species: {score_col: per-resample metric}} for the pairs of ``df`` (optionally masked).
    With ``n_jobs`` > 1 the score columns are split across joblib workers (identical results)."""
    out: dict[str, dict[str, np.ndarray]] = {}
    for sp in species:
        m = (df["species"] == sp).to_numpy()
        if mask is not None:
            m = m & np.asarray(mask, bool)
        if m.sum() == 0:
            continue
        sub = df.loc[m]
        g, y = sub["group"].to_numpy(), sub["label"].to_numpy()
        scores = {c: sub[c].to_numpy() for c in score_cols}
        if n_jobs > 1 and len(score_cols) > 1:
            from joblib import Parallel, delayed

            parts = [score_cols[i::n_jobs] for i in range(min(n_jobs, len(score_cols)))]
            res = Parallel(n_jobs=len(parts), backend="loky")(
                delayed(boot_metric)(bs, g, y, {c: scores[c] for c in part}, metric) for part in parts)
            out[sp] = {k: v for r in res for k, v in r.items()}
        else:
            out[sp] = boot_metric(bs, g, y, scores, metric=metric)
    return out


def species_mean_delta(boot: dict[str, dict[str, np.ndarray]], model: str, base: str) -> np.ndarray:
    """Per-resample T = mean over species of (metric_model − metric_base); NaN if any species is NaN."""
    deltas = [boot[sp][model] - boot[sp][base] for sp in boot]
    return np.mean(np.vstack(deltas), axis=0)
