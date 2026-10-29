"""Toy known-answer tests for the M2 development harness (metrics, bootstrap, weights, CV, baselines)."""

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score, roc_auc_score

from slparalog.eval import baselines, controls
from slparalog.eval.bootstrap import boot_metric, ci, group_strata, make_bootstrap, species_mean_delta
from slparalog.eval.cv import check_folds, oof_scores
from slparalog.eval.dataset import assemble, pair_key
from slparalog.eval.models import gbdt_complexity, grid, one_se_select
from slparalog.eval.stats import (ap_weighted_batch, auroc_weighted_batch, precision_at_k,
                                  species_balanced_weights, summary)

# ------------------------------------------------------------------ fast weighted AP / AUROC


def test_fast_ap_toy_known_answer():
    y = np.array([1, 0, 1, 0])
    s = np.array([0.9, 0.8, 0.7, 0.1])
    assert ap_weighted_batch(y, s, np.ones((1, 4)))[0] == pytest.approx(5 / 6)
    assert auroc_weighted_batch(y, s, np.ones((1, 4)))[0] == pytest.approx(3 / 4)


def test_fast_ap_ties_known_answer():
    # all scores tied: AP = prevalence, AUROC = 0.5
    y = np.array([1, 0, 0, 0])
    s = np.zeros(4)
    assert ap_weighted_batch(y, s, np.ones((1, 4)))[0] == pytest.approx(0.25)
    assert auroc_weighted_batch(y, s, np.ones((1, 4)))[0] == pytest.approx(0.5)


def test_fast_ap_matches_sklearn_with_weights_and_ties():
    rng = np.random.default_rng(0)
    for _ in range(30):
        n = int(rng.integers(5, 80))
        y = rng.integers(0, 2, n)
        y[0], y[1] = 1, 0
        s = rng.integers(0, 6, n).astype(float)  # many ties
        W = rng.integers(0, 4, (5, n)).astype(float)
        W[:, 0] = np.maximum(W[:, 0], 1)
        W[:, 1] = np.maximum(W[:, 1], 1)
        fa = ap_weighted_batch(y, s, W)
        fu = auroc_weighted_batch(y, s, W)
        for b in range(5):
            assert fa[b] == pytest.approx(average_precision_score(y, s, sample_weight=W[b]), abs=1e-12)
            assert fu[b] == pytest.approx(roc_auc_score(y, s, sample_weight=W[b]), abs=1e-12)


def test_multiplicity_weights_equal_replication():
    y = np.array([1, 0, 1, 0, 0, 1])
    s = np.array([0.3, 0.9, 0.3, 0.2, 0.5, 0.8])
    w = np.array([2, 1, 0, 3, 1, 1])
    rep_y, rep_s = np.repeat(y, w), np.repeat(s, w)
    assert ap_weighted_batch(y, s, w[None])[0] == pytest.approx(average_precision_score(rep_y, rep_s))


def test_fast_ap_no_positive_weight_is_nan():
    y = np.array([1, 0, 0])
    assert np.isnan(ap_weighted_batch(y, np.array([0.1, 0.2, 0.3]), np.array([[0, 1, 1]]))[0])


# ------------------------------------------------------------------ secondary metrics


def test_precision_at_k_known_answers():
    y = [1, 0, 1, 0, 0]
    s = [0.9, 0.8, 0.7, 0.1, 0.0]
    assert precision_at_k(y, s, 1) == 1.0
    assert precision_at_k(y, s, 2) == 0.5
    assert precision_at_k(y, s, 3) == pytest.approx(2 / 3)
    # tie at the k-th score: expected precision over the tie block
    assert precision_at_k([1, 0, 1, 0], [0.9, 0.5, 0.5, 0.5], 2) == pytest.approx((1 + 1 / 3) / 2)


def test_summary_lift_and_counts():
    r = summary([1, 0, 1, 0], [0.9, 0.8, 0.7, 0.1])
    assert (r["n"], r["n_pos"], r["n_neg"]) == (4, 2, 2)
    assert r["ap"] == pytest.approx(5 / 6)
    assert r["lift"] == pytest.approx((5 / 6) / 0.5)


# ------------------------------------------------------------------ weights


def test_species_balanced_weights_known_answer():
    sp = np.array(["a", "a", "a", "a", "b", "b"])
    y = np.array([1, 0, 0, 0, 1, 0])
    w = species_balanced_weights(sp, y)
    # N=6, S=2, pi* = mean(0.25, 0.5) = 0.375
    assert w.tolist() == pytest.approx([1.125, 0.625, 0.625, 0.625, 1.125, 1.875])
    assert w.sum() == pytest.approx(6)
    for s in "ab":
        m = sp == s
        assert w[m].sum() == pytest.approx(3)
        assert (w[m] * y[m]).sum() / w[m].sum() == pytest.approx(0.375)


def test_single_species_weights_are_uniform():
    w = species_balanced_weights(np.array(["a"] * 5), np.array([1, 0, 0, 1, 0]))
    assert w == pytest.approx(np.ones(5))


# ------------------------------------------------------------------ bootstrap


def _toy_pop():
    # groups: g1 human pos, g2 human neg, g3 human neg, g4 yeast pos, g5 mixed pos, g6 yeast neg
    rows = [
        ("human", "g1", 1, 0.9), ("human", "g1", 0, 0.2),
        ("human", "g2", 0, 0.5), ("human", "g3", 0, 0.1),
        ("yeast", "g4", 1, 0.8), ("yeast", "g4", 0, 0.3),
        ("human", "g5", 1, 0.7), ("yeast", "g5", 0, 0.6), ("yeast", "g5", 1, 0.4),
        ("yeast", "g6", 0, 0.05),
    ]
    return pd.DataFrame(rows, columns=["species", "group", "label", "score"])


def test_group_strata_composition():
    st = group_strata(_toy_pop())
    assert st.loc["g5", "composition"] == "mixed"
    assert st.loc["g1", "stratum"] == "human_only:1"
    assert st.loc["g2", "stratum"] == "human_only:0"
    assert st.loc["g6", "stratum"] == "yeast_only:0"


def test_bootstrap_preserves_stratum_sizes_and_is_deterministic():
    df = _toy_pop()
    bs = make_bootstrap(df, B=500, seed=11)
    bs2 = make_bootstrap(df, B=500, seed=11)
    assert np.array_equal(bs.counts, bs2.counts)
    st = bs.strata.loc[bs.groups]
    for name, sub in st.groupby("stratum"):
        cols = np.flatnonzero(st["stratum"].to_numpy() == name)
        assert (bs.counts[:, cols].sum(axis=1) == len(sub)).all()
    # every group is sometimes drawn more than once (with replacement) in the 2-group stratum
    assert bs.counts.max() >= 2


def test_paired_delta_identical_scores_is_zero():
    df = _toy_pop()
    bs = make_bootstrap(df, B=300, seed=1)
    out = boot_metric(bs, df["group"], df["label"], {"a": df["score"], "b": df["score"].copy()})
    d = out["a"] - out["b"]
    assert np.nanmax(np.abs(d)) == 0
    c = ci(d)
    assert c["lo"] == 0 and c["hi"] == 0


def test_bootstrap_point_mass_population():
    # one group per stratum -> every resample equals the full data: AP CI collapses to the point value
    df = pd.DataFrame({"species": ["human"] * 4, "group": ["g1", "g1", "g2", "g2"],
                       "label": [1, 0, 0, 0], "score": [0.9, 0.1, 0.5, 0.2]})
    bs = make_bootstrap(df, B=50, seed=3)
    out = boot_metric(bs, df["group"], df["label"], {"s": df["score"]})["s"]
    assert np.allclose(out, average_precision_score(df["label"], df["score"]))


def test_species_mean_delta():
    boot = {"human": {"m": np.array([0.5, 0.6]), "b": np.array([0.4, 0.4])},
            "yeast": {"m": np.array([0.3, 0.1]), "b": np.array([0.2, 0.2])}}
    assert species_mean_delta(boot, "m", "b") == pytest.approx([0.1, 0.05])


def test_ci_percentile_known_answer():
    x = np.arange(101, dtype=float)
    c = ci(np.r_[x, np.nan], 0.95)
    assert c["lo"] == pytest.approx(2.5) and c["hi"] == pytest.approx(97.5)
    assert c["n_dropped"] == 1 and c["n_ok"] == 101


# ------------------------------------------------------------------ CV and baselines


class _Recorder:
    """Toy model: remembers training groups; predicts 1 for rows whose group it saw in training."""

    def __init__(self, groups):
        self.groups = groups

    def fit(self, X, y, sample_weight=None):
        self.seen = set(X[:, 0].tolist())
        return self

    def predict(self, X):
        return np.array([1.0 if g in self.seen else 0.0 for g in X[:, 0]])


def _toy_cv_frame():
    rng = np.random.default_rng(0)
    groups = [f"g{i}" for i in range(20) for _ in range(3)]
    df = pd.DataFrame({"group": groups, "species": ["human"] * 30 + ["yeast"] * 30,
                       "label": rng.integers(0, 2, 60)})
    gid = df["group"].str[1:].astype(int)
    for s in (11, 23):
        perm = np.random.default_rng(s).permutation(20)
        df[f"fold_s{s}"] = (perm[gid] % 5)
    return df, gid.to_numpy().astype(float)[:, None]


def test_oof_never_scores_a_training_group():
    df, X = _toy_cv_frame()
    check_folds(df, [11, 23])
    sc, per_seed = oof_scores(df, X, lambda s: _Recorder(None), [11, 23], weights=None)
    assert np.all(per_seed == 0)  # no test row's group was seen in training
    assert np.isfinite(sc).all()


def test_loso_masks_train_one_species_only():
    df, X = _toy_cv_frame()
    sp = df["species"].to_numpy()

    class M:
        def fit(self, X, y, sample_weight=None):
            assert set(np.unique(X[:, 1])) == {0.0}
            return self

        def predict(self, X):
            return np.zeros(len(X))

    X2 = np.column_stack([X[:, 0], (sp == "yeast").astype(float)])
    sc, _ = oof_scores(df, X2, lambda s: M(), [11], train_mask=sp == "human", test_mask=sp == "yeast")
    assert np.isnan(sc[sp == "human"]).all() and np.isfinite(sc[sp == "yeast"]).all()


def test_b3_gene_prior_known_answer():
    pool = pd.DataFrame({
        "species": ["h"] * 5, "gene_a": ["a", "a", "c", "x", "d"], "gene_b": ["b", "c", "d", "y", "e"],
        "label": [1, 0, 1, 0, 0], "fold_s11": [0, 0, 1, 1, 2],
    })
    ev = pool.copy()
    sc, diag = baselines.b3_gene_prior(ev, pool, [11])
    # pair (a,b) in fold 0: training folds 1,2 -> c:1 (c,d), d: mean(1,0)=0.5, e:0, x,y:0; a,b unseen
    # -> fallback = training prevalence of folds {1,2} = 1/3
    assert sc[0] == pytest.approx(1 / 3)
    # pair (a,c) fold 0: c prior = 1.0 -> max(1/3, 1) = 1
    assert sc[1] == pytest.approx(1.0)
    # pair (c,d) fold 1: training folds 0,2 -> c: 0 (a,c), d: 0 (d,e) -> max = 0
    assert sc[2] == pytest.approx(0.0)
    # pair (d,e) fold 2: training folds 0,1 -> d: 1 (c,d); e unseen -> max(1, prev) = 1
    assert sc[4] == pytest.approx(1.0)
    assert 0 < diag["frac_pairs_with_gene_seen_in_training"] <= 1


def test_b1_variants_cover_preregistered_set():
    assert set(baselines.B1_VARIANTS) == {"ident_local_min", "ident_local_max", "ident_ens_min", "ident_ens_max",
                                          "ident_ens_mean", "ident_global", "ident_x_cov"}
    d = pd.DataFrame({"ident_local_min": [0.5], "ident_local_max": [0.6], "ident_ens_min": [0.2],
                      "ident_ens_max": [0.4], "ident_global": [0.3], "cov_min": [0.5]})
    v = baselines.b1_variant_scores(d).iloc[0]
    assert v["b1_ident_ens_mean"] == pytest.approx(0.3)
    assert v["b1_ident_x_cov"] == pytest.approx(0.3)


# ------------------------------------------------------------------ selection and controls


def test_one_se_rule_picks_simplest_within_one_se():
    rs = [
        {"name": "big", "params": {"max_leaf_nodes": 15, "max_iter": 300, "learning_rate": 0.1,
                                   "min_samples_leaf": 20, "l2_regularization": 0.0}, "T": 0.10, "T_se": 0.03},
        {"name": "small", "params": {"max_leaf_nodes": 7, "max_iter": 100, "learning_rate": 0.03,
                                     "min_samples_leaf": 50, "l2_regularization": 1.0}, "T": 0.08, "T_se": 0.03},
        {"name": "bad", "params": {"max_leaf_nodes": 7, "max_iter": 100, "learning_rate": 0.03,
                                   "min_samples_leaf": 20, "l2_regularization": 1.0}, "T": 0.01, "T_se": 0.03},
    ]
    r = one_se_select(rs, gbdt_complexity)
    assert r["best"]["name"] == "big" and r["selected"]["name"] == "small" and r["n_eligible"] == 2


def test_grid_order_is_deterministic():
    g = grid({"b": [1, 2], "a": [3]})
    assert g == [{"a": 3, "b": 1}, {"a": 3, "b": 2}]


def test_permute_within_species_keeps_prevalence():
    df = pd.DataFrame({"species": ["h"] * 10 + ["y"] * 10, "label": [1] * 3 + [0] * 7 + [1] * 5 + [0] * 5})
    yp = controls.permute_within_species(df, 11)
    assert yp[:10].sum() == 3 and yp[10:].sum() == 5


# ------------------------------------------------------------------ dataset assembly


def test_pair_key_is_order_invariant():
    lo, hi = pair_key(pd.Series(["b", "a"]), pd.Series(["a", "b"]))
    assert lo.tolist() == ["a", "a"] and hi.tolist() == ["b", "b"]


def test_assemble_joins_unordered_and_counts_uncalled(tmp_path):
    folds = {f"fold_s{s}": [0, 1, 2, 3] for s in (11, 23, 37, 41, 59)}
    hum = pd.DataFrame({"gene_a": ["hb", "hc", "hd", "he"], "gene_b": ["ha", "ha", "ha", "hf"],
                        "label_primary": [1, 0, -1, 1], "group": ["G1", "G2", "G3", "G4"],
                        "n_testable": [3, 4, 5, 6], "n_studies_testable": [1, 2, 3, 1],
                        "harle_chan_dependent": [False] * 4, **folds})
    yea = pd.DataFrame({"gene_a": ["ya", "yc"], "gene_b": ["yb", "yd"], "primary": ["pos", "neg"],
                        "group": ["G5", "G6"], "orientation_conflict": [0, 0],
                        **{k: v[:2] for k, v in folds.items()}})
    hp, yp = tmp_path / "h.tsv", tmp_path / "y.tsv"
    hum.to_csv(hp, sep="\t", index=False)
    yea.to_csv(yp, sep="\t", index=False)
    # caller table: (ha, hb) stored reversed, (ha, hc) present, (he, hf) missing -> not called
    fh = pd.DataFrame({"gene_a": ["ha", "ha", "hx"], "gene_b": ["hb", "hc", "hy"],
                       "ident_local_max": [0.9, 0.5, 0.1], "not_registered_col": [1, 2, 3]})
    fy = pd.DataFrame({"gene_a": ["yb", "yc"], "gene_b": ["ya", "yd"], "ident_local_max": [0.8, 0.2],
                       "not_registered_col": [1, 2]})
    fhp, fyp = tmp_path / "fh.parquet", tmp_path / "fy.parquet"
    fh.to_parquet(fhp)
    fy.to_parquet(fyp)
    A = assemble({"human": {"tierA": fhp}, "yeast": {"tierA": fyp}}, label_paths={"human": hp, "yeast": yp})
    c = A.counts()
    assert c["human"]["labelled"] == 3  # ambiguous (-1) dropped
    assert c["human"]["called"] == 2 and c["human"]["not_called"] == 1 and c["human"]["not_called_pos"] == 1
    assert c["yeast"]["called"] == 2
    assert "not_registered_col" not in A.primary.columns
    row = A.primary.set_index("pair_id").loc["human:ha|hb"]
    assert row["ident_local_max"] == 0.9 and row["label"] == 1
    assert A.feature_columns == ["ident_local_max"]


def test_gbdt_prebinning_matches_unweighted_sklearn():
    from sklearn.ensemble import HistGradientBoostingClassifier

    from slparalog.eval.models import GBDTModel

    rng = np.random.default_rng(5)
    X = rng.random((400, 4))
    X[:7, 1] = np.nan
    y = (X[:, 0] + 0.3 * rng.random(400) > 0.8).astype(int)
    p = {"max_iter": 20, "max_leaf_nodes": 7}
    ours = GBDTModel(["a", "b", "c", "d"], p, seed=0).fit(X, y)
    ref = HistGradientBoostingClassifier(early_stopping=False, random_state=0, **p).fit(X, y)
    assert np.allclose(ours.predict(X), ref.predict_proba(X)[:, 1])
    with pytest.raises(ValueError):
        GBDTModel(["a"], {"early_stopping": True}).fit(X[:, :1], y)


# ------------------------------------------------------------------ D12 additions


def _pf_frame():
    df = pd.DataFrame({
        "species": ["human"] * 8,
        "group": ["a", "a", "b", "b", "c", "c", "d", "d"],
        "label": [1, 0, 0, 1, 1, 0, 0, 0],
        "fold_s11": [0, 0, 0, 0, 1, 1, 1, 1],
    })
    return df


def test_perfold_point_known_answer():
    from slparalog.eval import perfold

    df = _pf_frame()
    base = np.array([0.9, 0.1, 0.8, 0.2, 0.3, 0.6, 0.5, 0.4])
    model = np.array([0.9, 0.1, 0.2, 0.8, 0.9, 0.1, 0.2, 0.3])
    pt = perfold.point(df, {"m": model[None]}, base, [11], species=("human",))
    # fold 0: base AP of y=[1,0,0,1] with s=[.9,.1,.8,.2]: ranks 1(+),0.8(-),0.2(+),0.1(-) -> (1 + 2/3)/2
    # model fold 0: perfect -> 1. fold 1: y=[1,0,0,0]; base s=[.3,.6,.5,.4] ranks + last -> 1/4; model perfect -> 1
    d0 = 1 - (1 + 2 / 3) / 2
    d1 = 1 - 1 / 4
    r = pt["m"]["per_species"]["human"]
    assert r["delta"] == pytest.approx((d0 + d1) / 2)
    assert r["n_folds"] == 2
    assert pt["m"]["T"] == pytest.approx((d0 + d1) / 2)


def test_perfold_bootstrap_identical_scores_zero_and_point_mass():
    from slparalog.eval import perfold

    df = _pf_frame()
    base = np.array([0.9, 0.1, 0.8, 0.2, 0.3, 0.6, 0.5, 0.4])
    bs = make_bootstrap(df, B=200, seed=2)
    bt = perfold.bootstrap(bs, df, {"same": base[None]}, base, [11], species=("human",))
    assert np.nanmax(np.abs(bt["same"]["T"])) == 0


def test_deconfound_covariates_and_reference():
    from slparalog.eval.models import DeconfoundModel, deconfound_covariates, reference_covariates

    c = deconfound_covariates([1, 2, 3, 7, np.nan, 1], ["human"] * 5 + ["yeast"])
    assert c.tolist() == [[1, 0], [0, 1], [0, 0], [0, 0], [0, 0], [0, 0]]
    assert reference_covariates(2, 1).tolist() == [[1, 0], [1, 0]]
    with pytest.raises(ValueError):
        reference_covariates(1, 4)

    class Echo:
        def fit(self, X, y, sample_weight=None):
            self.X = X
            return self

        def predict(self, X):
            return X[:, -2] * 10 + X[:, -1]

    X = np.column_stack([np.zeros(3), c[:3]])
    m = DeconfoundModel(Echo(), ref=2).fit(X, np.array([0, 1, 0]))
    assert np.array_equal(m.base.X, X)  # training sees the real covariates
    assert m.predict(X).tolist() == [1, 1, 1]  # prediction sees the reference level only
    assert X[:, 1:].tolist() == c[:3].tolist()  # caller's matrix untouched


def test_guarded_select():
    from slparalog.eval.models import guarded_select, lr_complexity

    rs = [{"name": "a", "params": {"C": 1.0}, "T": 0.05, "T_se": 0.01, "guard_pass": False},
          {"name": "b", "params": {"C": 0.1}, "T": 0.03, "T_se": 0.01, "guard_pass": True},
          {"name": "c", "params": {"C": 0.01}, "T": 0.025, "T_se": 0.01, "guard_pass": True}]
    r = guarded_select(rs, lr_complexity)
    assert r["best"]["name"] == "b" and r["selected"]["name"] == "c" and r["n_guard_pass"] == 2
    assert guarded_select([{**x, "guard_pass": False} for x in rs], lr_complexity)["selected"] is None


def test_explog_schema_growth(tmp_path):
    from slparalog.eval import explog

    p = tmp_path / "exp.tsv"
    pd.DataFrame({c: ["x"] for c in explog.COLUMNS[:5]}).to_csv(p, sep="\t", index=False)
    explog.append([{"exp_id": "new", "params": {"a": 1}}], path=p)
    t = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    assert list(t.columns) == explog.COLUMNS and t["exp_id"].tolist() == ["x", "new"]


# ------------------------------------------------------------------ D14 additions


def test_rank_ensemble_picks_identity_when_model_is_noise():
    from slparalog.eval.models import RankEnsembleModel

    rng = np.random.default_rng(0)
    n = 600
    y = (rng.random(n) < 0.15).astype(int)
    b1 = y + rng.normal(0, 0.6, n)  # informative "identity"
    noise = rng.normal(size=(n, 2))  # uninformative features

    class Mean:
        def fit(self, X, y, sample_weight=None):
            return self

        def predict(self, X):
            return X[:, 0]

    X = np.column_stack([noise, b1, np.arange(n) % 60, (np.arange(n) < 300).astype(float)])
    m = RankEnsembleModel(lambda: Mean()).fit(X, y)
    assert m.w_ == 0.0  # pure identity wins over noise
    s = m.predict(X)
    # w = 0: prediction is the per-species percentile rank of b1
    for c in (0.0, 1.0):
        mm = X[:, -1] == c
        assert np.array_equal(np.argsort(s[mm], kind="stable"), np.argsort(b1[mm], kind="stable"))


def test_rank_ensemble_picks_model_when_model_is_better():
    from slparalog.eval.models import RankEnsembleModel

    rng = np.random.default_rng(1)
    n = 600
    y = (rng.random(n) < 0.15).astype(int)
    good = y + rng.normal(0, 0.1, n)
    b1 = rng.normal(size=n)

    class First:
        def fit(self, X, y, sample_weight=None):
            return self

        def predict(self, X):
            return X[:, 0]

    X = np.column_stack([good, b1, np.arange(n) % 60, np.zeros(n)])
    assert RankEnsembleModel(lambda: First()).fit(X, y).w_ >= 0.9


def test_resolve_feature_set_groups_and_tiers():
    from slparalog.eval.models import resolve_feature_set

    avail = ["ident_global", "cov_min", "hits_all_min", "len_min", "len_ratio"]
    groups = {"sim": ["ident_", "cov_", "len_ratio"], "gene": ["len_min"]}
    assert resolve_feature_set("sim", avail, groups) == ["ident_global", "cov_min", "len_ratio"]
    assert resolve_feature_set("sim+gene", avail, groups) == ["ident_global", "cov_min", "len_min", "len_ratio"]
    assert resolve_feature_set("A", avail, groups) == avail
    with pytest.raises(ValueError):
        resolve_feature_set("Z", avail, groups)
