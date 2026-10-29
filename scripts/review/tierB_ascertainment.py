"""Review step 3: is the Tier B gain learned screening ascertainment?

(a) within each human source library (pairs whose `studies` lists the library): pooled seed-averaged OOF AP
    of B1, GBDT-A, GBDT-A+B against (i) the primary label and (ii) the library's own SL call
    (library in `studies_sl`); dAP with stratified group bootstrap (B=2000, seed 11).
(b) within n_studies_testable strata (1, 2, >=3): per-fold paired dAP and pooled dAP with CI.
(c) regressed out:
    c1 score level: grouped-CV LR (human only, stored folds) of label on ascertainment covariates
        (n_studies_testable indicators + 9 library indicators) + B1, then + GBDT-A score, then + GBDT-A+B score.
    c2 feature level: every Tier B feature residualised (OLS, human rows only) on the same covariates, A+B refit.
    c3 harness D12(c) deconfounded training (n_studies_testable covariate, predicted at reference >=3).
(d) label-free: does Tier B predict "this Compara paralog pair was screened"? Universe U1 = all human Compara
    paralog pairs called by the caller; U2 = those touching >=1 dev gene (which excludes every held-out test pair,
    since split groups are gene-connected). Target = pair is a dev-labelled pair (any tier). GBDT (same params,
    unweighted), 5-fold GroupKFold on gene-connected components. Plus the cross-fitted screening propensity
    used as an SL score on the dev folds.
Output results/review/tierB_gain/ascertainment.json"""

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

from tierB_common import (OUT, ROOT, SEEDS, SPECIES, boot_delta, dump, fit, load, mk_gb, perfold, pooled_ap)
from slparalog.eval.cv import fold_col, oof_many
from slparalog.eval.models import DeconfoundModel, GBDTModel, LRModel, deconfound_covariates
from slparalog.eval.stats import auroc
from slparalog.io.labels import load_labels

LIBS = ["chymera2020", "parrish2021", "in4mer2024", "thompson2021", "ito2021", "klingbeil2024", "dede2020",
        "chan2025", "harle2025"]


def human_meta(df):
    raw = load_labels(ROOT / "data" / "labels" / "dev" / "human_dev.tsv", low_memory=False)
    a, b = raw["gene_a"].astype(str), raw["gene_b"].astype(str)
    raw["pair_id"] = "human:" + np.where(a <= b, a, b) + "|" + np.where(a <= b, b, a)
    m = df[["pair_id"]].merge(raw[["pair_id", "studies", "studies_sl"]], on="pair_id", how="left")
    st = m["studies"].fillna("").astype(str)
    ssl = m["studies_sl"].fillna("").astype(str)
    lib = pd.DataFrame({L: st.str.split(",").apply(lambda x: L in x).to_numpy() for L in LIBS})
    libsl = pd.DataFrame({L: ssl.str.split(",").apply(lambda x: L in x).to_numpy() for L in LIBS})
    return lib, libsl


def main():
    df, fa, fb, _ = load()
    y = df["label"].to_numpy()
    sp = df["species"].to_numpy()
    hm = sp == "human"
    base = df["B1"].to_numpy()
    core = np.load(OUT / "oof_core.npz")
    A_m, AB_m = core["A"], core["AB"]
    S = {"B1": base, "A": A_m.mean(0), "AB": AB_m.mean(0), "B": core["B"].mean(0)}
    lib, libsl = human_meta(df)
    nst = pd.to_numeric(df["meta_n_studies_testable"], errors="coerce").to_numpy()
    out = {"note": "pooled = seed-averaged OOF scores; perfold = harness per-fold paired statistic"}

    # (a) within library
    within = {}
    for L in LIBS:
        m = hm & lib[L].to_numpy()
        r = boot_delta(df, m, S, "B1")
        r["delta_AB_vs_A"] = r["ap"].get("AB", np.nan) - r["ap"].get("A", np.nan)
        rA = boot_delta(df, m, {"A": S["A"], "AB": S["AB"]}, "A")
        r["delta_AB_vs_A_ci"] = rA["delta_vs_A"].get("AB")
        # library's own SL call
        dfl = df.copy()
        dfl["label"] = libsl[L].to_numpy().astype(int)
        rl = boot_delta(dfl, m, S, "B1")
        rl2 = boot_delta(dfl, m, {"A": S["A"], "AB": S["AB"]}, "A")
        rl["delta_AB_vs_A_ci"] = rl2["delta_vs_A"].get("AB")
        within[L] = {"primary_label": r, "library_own_call": rl}
        print(L, r["n"], r["n_pos"], {k: round(v, 3) for k, v in r["ap"].items()}, "own:",
              rl["n_pos"], {k: round(v, 3) for k, v in rl["ap"].items()})
    out["within_library"] = within

    # (b) within n_studies_testable strata
    strata = {}
    for name, m in (("1", nst == 1), ("2", nst == 2), (">=3", nst >= 3)):
        mm = hm & m
        r = boot_delta(df, mm, S, "B1")
        r["delta_AB_vs_A_ci"] = boot_delta(df, mm, {"A": S["A"], "AB": S["AB"]}, "A")["delta_vs_A"].get("AB")
        r["perfold_AB"] = perfold(df, AB_m, base, species=("human",), mask=mm)
        r["perfold_A"] = perfold(df, A_m, base, species=("human",), mask=mm)
        strata[name] = r
        print("stratum", name, r["n_pos"], {k: round(v, 3) for k, v in r["ap"].items()},
              r["perfold_AB"].get("human"), r["perfold_A"].get("human"))
    out["within_n_studies_strata"] = strata

    # (c1) score-level regression with ascertainment covariates (human only)
    cov = np.column_stack([(nst == 1), (nst == 2)] + [lib[L].to_numpy() for L in LIBS]).astype(float)
    hdf = df.loc[hm].reset_index(drop=True)
    ch = cov[hm]

    def lgt(p):
        p = np.clip(p, 1e-6, 1 - 1e-6)
        return np.log(p / (1 - p))
    sets = {"ASC": ch, "ASC+B1": np.column_stack([ch, base[hm]]),
            "ASC+B1+A": np.column_stack([ch, base[hm], lgt(S["A"][hm])]),
            "ASC+B1+AB": np.column_stack([ch, base[hm], lgt(S["AB"][hm])]),
            "ASC+B1+A+AB": np.column_stack([ch, base[hm], lgt(S["A"][hm]), lgt(S["AB"][hm])]),
            "B1+AB": np.column_stack([base[hm], lgt(S["AB"][hm])])}
    specs = {k: dict(make_model=(lambda cols: (lambda s: LRModel(cols, C=1.0, seed=s)))
                     ([f"x{i}" for i in range(v.shape[1])]), X=v, weights=None) for k, v in sets.items()}
    r = oof_many(hdf, specs, SEEDS, n_jobs=6)
    c1 = {}
    for k in sets:
        pf = perfold(hdf, r[k][1], hdf["B1"].to_numpy(), species=("human",))
        c1[k] = {"pooled_ap": pooled_ap(hdf["label"].to_numpy(), r[k][0]), "perfold_dAP_vs_B1": pf["human"]["delta"],
                 "perfold_ap": pf["human"]["ap_model"]}
    # in-sample LR on all human rows: coefficient / LR test of the AB score given ASC + B1 + A
    from scipy.stats import chi2
    from sklearn.linear_model import LogisticRegression
    yh = hdf["label"].to_numpy()

    def llf(X):
        m = LogisticRegression(C=1e6, max_iter=10000).fit(X, yh)
        pr = np.clip(m.predict_proba(X)[:, 1], 1e-12, 1 - 1e-12)
        return float(np.sum(yh * np.log(pr) + (1 - yh) * np.log(1 - pr))), m
    l0, _ = llf(sets["ASC+B1+A"])
    l1, m1 = llf(sets["ASC+B1+A+AB"])
    llr = 2 * (l1 - l0)
    c1["in_sample_LR_test_AB_given_ASC_B1_A"] = {"llr": llr, "p_chi2_df1": float(chi2.sf(llr, 1)),
                                                 "coef_AB_logit": float(m1.coef_[0][-1]),
                                                 "note": "in-sample (OOF scores as covariates), unpenalised"}
    out["regress_out_score_level"] = c1
    print("c1", {k: (round(v["pooled_ap"], 3), round(v["perfold_dAP_vs_B1"], 3)) for k, v in c1.items() if "pooled_ap" in v},
          c1["in_sample_LR_test_AB_given_ASC_B1_A"])

    # (c2) Tier B residualised on covariates in human rows
    XB = df[fb].to_numpy(float).copy()
    D = np.column_stack([np.ones(hm.sum()), ch])
    r2 = {}
    for j in range(XB.shape[1]):
        col = XB[hm, j]
        ok = np.isfinite(col)
        beta, *_ = np.linalg.lstsq(D[ok], col[ok], rcond=None)
        pred = D @ beta
        r2[fb[j]] = float(1 - np.nanvar(col - pred) / np.nanvar(col))
        newc = col - pred + np.nanmean(col)
        XB[np.flatnonzero(hm), j] = np.where(ok, newc, np.nan)
    XA = df[fa].to_numpy(float)
    rr = fit(df, {"AB_resid": dict(make_model=mk_gb(fa + fb), X=np.column_stack([XA, XB]))})
    out["regress_out_feature_level"] = {"R2_tierB_on_covariates_human": r2,
                                        "AB_resid": perfold(df, rr["AB_resid"][1], base)}
    print("c2", out["regress_out_feature_level"]["AB_resid"])

    # (c3) deconfounded training (D12(c))
    Xdc = np.column_stack([XA, df[fb].to_numpy(float), deconfound_covariates(nst, sp)])
    cdc = fa + fb + ["dc1", "dc2"]
    Xdca = np.column_stack([XA, deconfound_covariates(nst, sp)])
    p = {"early_stopping": False, "l2_regularization": 1.0, "learning_rate": 0.03, "max_iter": 100,
         "max_leaf_nodes": 7, "min_samples_leaf": 50}
    rdc = fit(df, {"AB_dc": dict(make_model=lambda s: DeconfoundModel(GBDTModel(cdc, p, seed=s), ref=3), X=Xdc),
                   "A_dc": dict(make_model=lambda s: DeconfoundModel(GBDTModel(fa + ["dc1", "dc2"], p, seed=s), ref=3),
                                X=Xdca)})
    out["deconfounded_training"] = {k: perfold(df, v[1], base) for k, v in rdc.items()}
    for k, v in rdc.items():
        mm = v[0]
        out["deconfounded_training"][k]["strata_pooled_ap"] = {
            n: pooled_ap(y[hm & m], mm[hm & m]) for n, m in (("1", nst == 1), ("2", nst == 2), (">=3", nst >= 3))}
    print("c3", {k: (v["T"], v["human"]["delta"]) for k, v in out["deconfounded_training"].items()})

    # (d) universe: predicting screening from features (label-free)
    u = pd.read_csv(ROOT / "data" / "interim" / "universe" / "human_paralog_pairs.tsv", sep="\t")
    a, b = u["gene_a"].astype(str), u["gene_b"].astype(str)
    u["ga"], u["gb"] = np.where(a <= b, a, b), np.where(a <= b, b, a)
    u = u.drop_duplicates(["ga", "gb"])
    ta = pd.read_parquet(ROOT / "data/interim/features/human/tierA.parquet")
    tb = pd.read_parquet(ROOT / "data/interim/features/human/tierB.parquet")
    a, b = ta["gene_a"].astype(str), ta["gene_b"].astype(str)
    ta["ga"], ta["gb"] = np.where(a <= b, a, b), np.where(a <= b, b, a)
    tb = tb.drop(columns=["gene_a", "gene_b"])
    feat = pd.concat([ta[["ga", "gb"] + fa], tb[fb]], axis=1)
    U = u[["ga", "gb"]].merge(feat, on=["ga", "gb"], how="inner")
    raw = load_labels(ROOT / "data" / "labels" / "dev" / "human_dev.tsv", low_memory=False)
    a, b = raw["gene_a"].astype(str), raw["gene_b"].astype(str)
    dev = pd.DataFrame({"ga": np.where(a <= b, a, b), "gb": np.where(a <= b, b, a)})
    dev["screened"] = 1
    U = U.merge(dev.drop_duplicates(), on=["ga", "gb"], how="left")
    U["screened"] = U["screened"].fillna(0).astype(int)
    devgenes = set(dev["ga"]) | set(dev["gb"])
    # gene-connected components for grouping
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for x1, x2 in zip(U["ga"], U["gb"]):
        r1, r2_ = find(x1), find(x2)
        if r1 != r2_:
            parent[r1] = r2_
    U["comp"] = [find(x) for x in U["ga"]]
    U["ident_global_b1"] = U["ident_global"]
    uni = {}
    for uname, um in (("U1_all_called_compara", np.ones(len(U), bool)),
                      ("U2_touching_dev_gene", (U["ga"].isin(devgenes) | U["gb"].isin(devgenes)).to_numpy())):
        Uu = U.loc[um].reset_index(drop=True)
        comps = Uu["comp"].value_counts()
        gkf = GroupKFold(n_splits=5)
        res = {"n_pairs": int(len(Uu)), "n_screened": int(Uu["screened"].sum()),
               "n_components": int(len(comps)), "largest_component_pairs": int(comps.iloc[0])}
        for fsname, cols in (("B1_identity", ["ident_global"]), ("A", fa), ("B", fb), ("A+B", fa + fb)):
            Xu = Uu[cols].to_numpy(float)
            sc = np.full(len(Uu), np.nan)
            for tr, te in gkf.split(Xu, Uu["screened"], Uu["comp"]):
                if fsname == "B1_identity":
                    sc[te] = Xu[te, 0]
                    continue
                m = GBDTModel(cols, p, seed=11).fit(Xu[tr], Uu["screened"].to_numpy()[tr])
                sc[te] = m.predict(Xu[te])
            res[fsname] = {"auroc": auroc(Uu["screened"].to_numpy(), sc),
                           "ap": pooled_ap(Uu["screened"].to_numpy(), sc)}
        res["prevalence"] = float(Uu["screened"].mean())
        uni[uname] = res
        print(uname, res)
    out["universe_screening_prediction"] = uni
    out["universe_note"] = ("U1 counts held-out test screened pairs as unscreened (unknown to the reviewer), which "
                            "attenuates U1; U2 contains no held-out test pair.")

    # (d2) cross-fitted screening propensity (A+B and B features) as an SL score on dev folds (human)
    hdf2 = df.loc[hm].reset_index(drop=True)
    keyh = hdf2["pair_id"].str.replace("human:", "", regex=False).str.split("|", expand=True)
    hdf2["ga"], hdf2["gb"] = keyh[0], keyh[1]
    prop = {k: np.full((len(SEEDS), len(hdf2)), np.nan) for k in ("A+B", "B")}
    for i, s in enumerate(SEEDS):
        f = hdf2[fold_col(s)].to_numpy()
        for k in np.unique(f):
            te = f == k
            gk = set(hdf2.loc[te, "ga"]) | set(hdf2.loc[te, "gb"])
            trm = ~(U["ga"].isin(gk) | U["gb"].isin(gk))
            Ut = U.loc[trm]
            for name, cols in (("A+B", fa + fb), ("B", fb)):
                m = GBDTModel(cols, p, seed=s).fit(Ut[cols].to_numpy(float), Ut["screened"].to_numpy())
                prop[name][i, te] = m.predict(hdf2.loc[te, cols].to_numpy(float))
    bh = hdf2["B1"].to_numpy()
    out["screening_propensity_as_SL_score"] = {
        k: {"perfold": perfold(hdf2, v, bh, species=("human",)),
            "pooled_ap": pooled_ap(hdf2["label"].to_numpy(), v.mean(0)),
            "auroc_vs_n_studies_ge3": auroc((pd.to_numeric(hdf2["meta_n_studies_testable"]) >= 3).astype(int).to_numpy(),
                                            v.mean(0))} for k, v in prop.items()}
    # does AB add beyond B1 + propensity? grouped-CV LR on human
    sets2 = {"B1+prop": np.column_stack([bh, lgt(prop["A+B"].mean(0))]),
             "B1+prop+AB": np.column_stack([bh, lgt(prop["A+B"].mean(0)), lgt(S["AB"][hm])]),
             "B1+AB": np.column_stack([bh, lgt(S["AB"][hm])])}
    specs2 = {k: dict(make_model=(lambda cols: (lambda s: LRModel(cols, C=1.0, seed=s)))
                      ([f"x{i}" for i in range(v.shape[1])]), X=v, weights=None) for k, v in sets2.items()}
    r3 = oof_many(hdf2, specs2, SEEDS, n_jobs=6)
    out["propensity_stack"] = {k: perfold(hdf2, r3[k][1], bh, species=("human",))["human"] for k in sets2}
    print("prop", {k: v["perfold"]["human"] for k, v in out["screening_propensity_as_SL_score"].items()})
    print("stack", out["propensity_stack"])
    dump(out, "ascertainment.json")


if __name__ == "__main__":
    main()
