"""Robustness test: model trained on natural features, evaluated on features of degraded proteomes
(DESIGN §6.5 "Robustness test"; exploratory, reporting level 3; template: the thinning N->T design of
``slparalog.eval.thinning``).

Views are the model features of the SAME labelled & called dev pairs (rows of the natural primary population, same
order) taken (N) from the natural tables and (D_<kind>) from the feature table computed on the degraded proteome of
``kind`` (``scripts/robustness/degrade_proteome.py``; gene-ID mapping table). A pair is **lost** in a degraded view if
(i) one of its genes has no degraded record (removed gene) or (ii) both genes survive but the caller no longer calls
the pair. Lost pairs are scored at the minimum (one tie block: the smallest finite score of that score vector minus 1,
as in ``controls.end_to_end_ap``), for the model and for B1 alike. The population (natural primary pairs) and the
folds are never changed, so natural and degraded APs are paired on identical pairs.

Conditions: ``N->N`` (the incumbent, reproduced) and ``N-><kind>`` (train on the natural features of the training
folds, evaluate the held-out fold on the degraded features). Statistic: per-fold AP of the model and of B1 (the same
seed x fold test sets as ``slparalog.eval.perfold``), their paired difference, and the natural -> degraded drops, each
with a stratified group-bootstrap CI (same resamples throughout, hence paired).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from slparalog.eval.bootstrap import boot_metric, ci
from slparalog.eval.cv import fold_col
from slparalog.eval.dataset import load_feature_table, pair_key
from slparalog.eval.metrics import ap

RULE_TEXT = (
    "ROBUSTNESS TEST (DESIGN §6.5 'Robustness test'; reporting level 3 = exploratory; NO adoption decision depends on it). "
    "Model fixed: GBDT_AB_17 on feature set A+B (tb_asym_* excluded), hyper-parameters fixed, no re-tuning, seeds "
    "config/seeds.yaml, same dev folds as the incumbent. Degraded proteomes (scripts/robustness/degrade_proteome.py, seed 11 = "
    "seeds[0]): fragments (20% of proteins truncated to a random 50-80% of their length at a random N- or C-terminus), split (5% "
    "of genes split at a random internal point; longer piece keeps the gene ID, shorter gets a new ID), fused (2% of genes get "
    "another random non-labelled gene's sequence appended; donor removed), missing (15% of genes removed), haplotig (5% of genes "
    "get an extra near-identical copy, 1% random substitutions, new ID), combined (all of them: missing -> fused -> split -> "
    "fragments -> haplotig). Yeast: all six kinds (+ an identity control = natural proteome through the same route); human: "
    "combined only. Tier A+B features of each degraded proteome from `slparalog predict --features-only --tier-b` (existing panel "
    "databases; <=8 threads; runtime recorded). Evaluation (thinning N->T design): train on the NATURAL features of the training "
    "folds, evaluate the held-out fold on the DEGRADED features of the same labelled pairs; pairs map to degraded gene IDs "
    "through the mapping table; pairs that lose a gene or are no longer called are scored at the minimum score (one tie block), "
    "B1 (ident_global) identically; the population (natural labelled & called dev pairs) is unchanged. Reported per species x "
    "kind: per-fold mean AP of the model and of B1, natural and degraded, dAP vs B1 per species (T_pf for the combined kind, where "
    "both species exist), and the natural-to-degraded drops (model, B1, model-minus-B1) with 95% stratified group-bootstrap "
    "percentile CIs (same resamples, paired); pooled-OOF APs, prevalence and counts, and the number of lost pairs. No decision "
    "rule: no thresholds, nothing is adopted or rejected on this result."
)


# --------------------------------------------------------------------------- views
def map_pairs(gene_a: np.ndarray, gene_b: np.ndarray, mapping: dict[str, str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map pair genes through ``mapping`` ({original ID: degraded ID}; absent = removed).
    Returns (mapped a, mapped b, gene_lost mask); mapped IDs are '' where lost."""
    ma = np.array([mapping.get(g, "") for g in gene_a], dtype=object)
    mb = np.array([mapping.get(g, "") for g in gene_b], dtype=object)
    return ma, mb, (ma == "") | (mb == "")


def degraded_frame(pairs: pd.DataFrame, mapping: dict[str, str], table: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Rows of the degraded feature ``table`` (keyed by unordered gene pair, registered columns) for the pairs of
    ``pairs`` (columns gene_a, gene_b). Lost pairs get all-NaN rows. Returns (frame aligned to ``pairs``, counts)."""
    ma, mb, gene_lost = map_pairs(pairs["gene_a"].to_numpy(), pairs["gene_b"].to_numpy(), mapping)
    lo, hi = pair_key(pd.Series(ma), pd.Series(mb))
    tlo, thi = pair_key(table["gene_a"], table["gene_b"])  # tables from load_feature_table are already ordered
    idx = pd.MultiIndex.from_arrays([tlo.to_numpy(), thi.to_numpy()])
    pos = idx.get_indexer(pd.MultiIndex.from_arrays([lo.to_numpy(), hi.to_numpy()]))
    pos = np.where(gene_lost, -1, pos)
    feat = table.drop(columns=["gene_a", "gene_b"])
    frame = feat.iloc[np.where(pos >= 0, pos, 0)].reset_index(drop=True).astype(float)
    lost = pos < 0
    frame.loc[lost, :] = np.nan
    counts = {"n_pairs": int(len(pairs)), "n_lost": int(lost.sum()), "n_lost_gene_removed": int(gene_lost.sum()),
              "n_lost_not_called": int((lost & ~gene_lost).sum())}
    frame["_lost"] = lost
    frame["_gene_lost"] = gene_lost
    return frame, counts


def load_degraded(pairs: pd.DataFrame, features_path, mapping_path) -> tuple[pd.DataFrame, dict]:
    m = pd.read_csv(mapping_path, sep="\t", dtype=str, keep_default_na=False)
    g = m.loc[(m["role"] == "gene") & (m["degraded_id"] != "")]
    mapping = dict(zip(g["orig_id"], g["degraded_id"]))
    return degraded_frame(pairs, mapping, load_feature_table(features_path))


def fill_min(score: np.ndarray, lost: np.ndarray) -> np.ndarray:
    """Lost pairs at the minimum score minus one (one tie block); ``score`` 1-D or (n_seeds, n)."""
    s = np.array(score, dtype=float, copy=True)
    if s.ndim == 1:
        s[lost] = np.nanmin(s[~lost]) - 1.0
    else:
        for i in range(s.shape[0]):
            s[i, lost] = np.nanmin(s[i, ~lost]) - 1.0
    return s


# --------------------------------------------------------------------------- per-fold AP
def _folds(df: pd.DataFrame, seeds):
    sp_arr = df["species"].to_numpy()
    return [(df[fold_col(s)].to_numpy(), sp_arr) for s in seeds]


def perfold_ap_point(df: pd.DataFrame, scores: dict[str, np.ndarray], seeds: list[int], species) -> dict:
    """{name: {sp: {ap, n_folds}}}: mean over the seed x fold test sets (>=1 positive of the species) of sklearn AP.
    ``scores[name]``: (n_seeds, n) or (n,) (the same for every seed, e.g. B1)."""
    y = df["label"].to_numpy()
    out = {n: {} for n in scores}
    for n, mat in scores.items():
        mat = np.atleast_2d(mat)
        for sp in species:
            vals = []
            for i, (f, sp_arr) in enumerate(_folds(df, seeds)):
                row = mat[i if mat.shape[0] > 1 else 0]
                for k in np.unique(f):
                    m = (f == k) & (sp_arr == sp)
                    if m.sum() == 0 or y[m].sum() == 0:
                        continue
                    vals.append(ap(y[m], row[m]))
            if vals:
                out[n][sp] = {"ap": float(np.mean(vals)), "n_folds": len(vals)}
    return out


def _boot_seed(bs, df_small, rows: dict[str, np.ndarray], seed, species):
    y = df_small["label"].to_numpy()
    sp_arr = df_small["species"].to_numpy()
    g = df_small["group"].to_numpy()
    f = df_small[fold_col(seed)].to_numpy()
    sums = {sp: {n: np.zeros(bs.B) for n in rows} for sp in species}
    cnts = {sp: np.zeros(bs.B) for sp in species}
    for sp in species:
        for k in np.unique(f):
            m = (f == k) & (sp_arr == sp)
            if m.sum() == 0 or y[m].sum() == 0:
                continue
            b = boot_metric(bs, g[m], y[m], {n: r[m] for n, r in rows.items()}, check=False)
            ok = np.isfinite(next(iter(b.values())))  # finiteness depends on labels and weights only
            for n in rows:
                sums[sp][n][ok] += b[n][ok]
            cnts[sp][ok] += 1
    return sums, cnts


def perfold_ap_bootstrap(bs, df: pd.DataFrame, scores: dict[str, np.ndarray], seeds: list[int], species,
                         n_jobs: int = 1) -> dict[str, dict[str, np.ndarray]]:
    """{name: {sp: (B,)}} per-resample mean over seed x fold test sets of AP (NaN: no fold with positive weight)."""
    df_small = df[["species", "group", "label"] + [fold_col(s) for s in seeds]].reset_index(drop=True)
    mats = {n: np.atleast_2d(np.asarray(m, float)) for n, m in scores.items()}
    args = [(bs, df_small, {n: (m[i] if m.shape[0] > 1 else m[0]) for n, m in mats.items()}, s, species)
            for i, s in enumerate(seeds)]
    if n_jobs > 1:
        from joblib import Parallel, delayed

        parts = Parallel(n_jobs=min(n_jobs, len(args)), backend="loky")(delayed(_boot_seed)(*a) for a in args)
    else:
        parts = [_boot_seed(*a) for a in args]
    out: dict[str, dict[str, np.ndarray]] = {n: {} for n in scores}
    for sp in species:
        c_ = sum(p[1][sp] for p in parts)
        for n in scores:
            s_ = sum(p[0][sp][n] for p in parts)
            with np.errstate(invalid="ignore", divide="ignore"):
                out[n][sp] = np.where(c_ > 0, s_ / np.where(c_ > 0, c_, 1), np.nan)
    return out


# --------------------------------------------------------------------------- summaries
def drop_summary(point: dict, boot: dict, kind_key: str, species, level: float = 0.95) -> dict:
    """Natural -> degraded drops per species (and the species mean 'T' if all of ``species`` are given).

    ``point[name][sp]['ap']`` / ``boot[name][sp]`` with names ``model|N``, ``b1|N``, ``model|<kind>``, ``b1|<kind>``.
    Returns {sp: {ap_model_nat, ap_model_deg, ap_b1_nat, ap_b1_deg, delta_nat, delta_deg, drop_model, drop_b1,
    drop_delta (each with CI)}}; delta = AP(model) - AP(B1)."""
    out = {}
    for sp in species:
        mn, md = point["model|N"][sp]["ap"], point[f"model|{kind_key}"][sp]["ap"]
        bn, bd = point["b1|N"][sp]["ap"], point[f"b1|{kind_key}"][sp]["ap"]
        Bm_n, Bm_d = boot["model|N"][sp], boot[f"model|{kind_key}"][sp]
        Bb_n, Bb_d = boot["b1|N"][sp], boot[f"b1|{kind_key}"][sp]
        out[sp] = {
            "n_folds": point["model|N"][sp]["n_folds"],
            "ap_model_natural": mn, "ap_model_degraded": md, "ap_b1_natural": bn, "ap_b1_degraded": bd,
            "delta_natural": mn - bn, "delta_degraded": md - bd,
            "delta_natural_ci": ci(Bm_n - Bb_n, level), "delta_degraded_ci": ci(Bm_d - Bb_d, level),
            "drop_model": mn - md, "drop_model_ci": ci(Bm_n - Bm_d, level),
            "drop_b1": bn - bd, "drop_b1_ci": ci(Bb_n - Bb_d, level),
            "drop_delta": (mn - bn) - (md - bd), "drop_delta_ci": ci((Bm_n - Bb_n) - (Bm_d - Bb_d), level),
        }
    if len(species) > 1:  # T_pf = species mean of dAP vs B1
        Tn = np.mean(np.vstack([boot["model|N"][sp] - boot["b1|N"][sp] for sp in species]), axis=0)
        Td = np.mean(np.vstack([boot[f"model|{kind_key}"][sp] - boot[f"b1|{kind_key}"][sp] for sp in species]), axis=0)
        tn = float(np.mean([out[sp]["delta_natural"] for sp in species]))
        td = float(np.mean([out[sp]["delta_degraded"] for sp in species]))
        out["T_pf"] = {"natural": tn, "natural_ci": ci(Tn, level), "degraded": td, "degraded_ci": ci(Td, level),
                       "drop": tn - td, "drop_ci": ci(Tn - Td, level)}
    return out
