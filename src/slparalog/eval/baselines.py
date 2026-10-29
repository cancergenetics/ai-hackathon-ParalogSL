"""Baselines B0-B4 (DESIGN v3 §5; RESEARCH_TREE D8).

- B0 random: uniform scores from a seeded generator.
- B1 identity alone: one of the pre-registered identity variants (``B1_VARIANTS``), chosen on
  development data by species-mean AP. Identity needs no fitting, so its OOF score is the variant
  itself; every variant is logged.
- B2 identity + family size + closest-pair flag, logistic regression (columns in config/models.yaml).
- B3 gene-prior leakage baseline: per fold, score(a, b) = max over the two genes of the mean training
  label of the training pairs containing that gene; genes without training pairs get the training
  prevalence of their species.
- B4 (human, diagnostic only, never a tool feature): number of testable cell lines from the labels.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from slparalog.eval.cv import fold_col

# Pre-registered B1 candidate set (DESIGN v3 §5; reviews/round1/statistics.md item 8).
# name -> (description, function of the feature frame)
B1_VARIANTS: dict[str, tuple[str, object]] = {
    "ident_local_min": ("local identity, min over orientations", lambda d: d["ident_local_min"]),
    "ident_local_max": ("local identity, max over orientations", lambda d: d["ident_local_max"]),
    "ident_ens_min": ("Ensembl-style identity, min (identical / longer)", lambda d: d["ident_ens_min"]),
    "ident_ens_max": ("Ensembl-style identity, max (identical / shorter)", lambda d: d["ident_ens_max"]),
    "ident_ens_mean": ("Ensembl-style identity, mean of min and max", lambda d: (d["ident_ens_min"] + d["ident_ens_max"]) / 2),
    "ident_global": ("global Needleman-Wunsch identity", lambda d: d["ident_global"]),
    "ident_x_cov": ("local identity (max over orientations) x min coverage of the two proteins",
                    lambda d: d["ident_local_max"] * d["cov_min"]),
}


def b0_random(n: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).random(n)


def b1_variant_scores(df: pd.DataFrame) -> pd.DataFrame:
    """All B1 candidate scores for the rows of ``df`` (columns ``b1_<variant>``)."""
    return pd.DataFrame({f"b1_{k}": np.asarray(fn(df), dtype=float) for k, (_, fn) in B1_VARIANTS.items()},
                        index=df.index)


def b3_gene_prior(eval_df: pd.DataFrame, train_pool: pd.DataFrame, seeds: list[int]) -> tuple[np.ndarray, dict]:
    """Gene-prior leakage baseline, OOF over the stored folds, averaged over seeds.

    ``train_pool``: all labelled pairs of the tier (called or not) with species, gene_a, gene_b, label,
    fold columns. For each seed and fold k, gene priors come from pool pairs in folds ≠ k; each
    evaluation pair in fold k gets max(prior(a), prior(b)), with the species' training prevalence for
    genes absent from the training pairs. Returns (scores, diagnostics: fraction of evaluation pairs
    with at least one gene seen in training).
    """
    n = len(eval_df)
    per_seed = np.full((len(seeds), n), np.nan)
    seen_frac = []
    for i, seed in enumerate(seeds):
        fc = fold_col(seed)
        folds_e = eval_df[fc].to_numpy()
        n_seen = 0
        for k in np.unique(folds_e):
            tr = train_pool.loc[train_pool[fc] != k]
            long = pd.concat([
                tr[["species", "gene_a", "label"]].rename(columns={"gene_a": "gene"}),
                tr[["species", "gene_b", "label"]].rename(columns={"gene_b": "gene"}),
            ])
            prior = long.groupby(["species", "gene"])["label"].mean()
            prev = tr.groupby("species")["label"].mean()
            te = folds_e == k
            sub = eval_df.loc[te]
            ka = pd.MultiIndex.from_arrays([sub["species"], sub["gene_a"]])
            kb = pd.MultiIndex.from_arrays([sub["species"], sub["gene_b"]])
            pa = prior.reindex(ka).to_numpy()
            pb = prior.reindex(kb).to_numpy()
            fb = sub["species"].map(prev).to_numpy()
            n_seen += int((np.isfinite(pa) | np.isfinite(pb)).sum())
            pa = np.where(np.isfinite(pa), pa, fb)
            pb = np.where(np.isfinite(pb), pb, fb)
            per_seed[i, te] = np.maximum(pa, pb)
        seen_frac.append(n_seen / n)
    return per_seed.mean(axis=0), {"frac_pairs_with_gene_seen_in_training": float(np.mean(seen_frac))}
