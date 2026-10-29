"""Review step 4: family concentration and cross-species transfer of the Tier B gain.

(a) The 5 largest positive-holding split groups (by primary positives, both species) are dropped in turn
    (and all 5 together) from training AND evaluation; GBDT A and A+B refitted; per-fold dAP vs B1.
(b) LOSO on the stored folds: train on one species (folds != k), test on the other species (fold k);
    A vs A+B (uniform weights, single species).
(c) Within-ascertainment cells: human pairs grouped by their exact set of testing studies (`studies`);
    positive-weighted mean of the within-cell AP over cells with >=1 positive and >=1 negative.
Output results/review/tierB_gain/concentration.json"""

import numpy as np
import pandas as pd

from tierB_common import OUT, ROOT, SEEDS, SPECIES, dump, fit, load, mk_gb, perfold, pooled_ap
from slparalog.eval.bootstrap import make_bootstrap, ci
from slparalog.io.labels import load_labels


def main():
    df, fa, fb, _ = load()
    y = df["label"].to_numpy()
    sp = df["species"].to_numpy()
    out = {}
    g = df.loc[y == 1].groupby("group").agg(n_pos=("label", "size"),
                                              human_pos=("species", lambda s: int((s == "human").sum())),
                                              yeast_pos=("species", lambda s: int((s == "yeast").sum())))
    g["n_pairs"] = df.groupby("group").size().reindex(g.index)
    top = g.sort_values(["n_pos", "n_pairs"], ascending=False).head(5)
    out["top5_groups"] = top.reset_index().to_dict(orient="records")
    out["total_pos"] = {c: int(y[sp == c].sum()) for c in SPECIES}
    drops = {f"drop_{gid}": [gid] for gid in top.index}
    drops["drop_all5"] = list(top.index)
    drops["none"] = []
    res = {}
    for name, gl in drops.items():
        keep = ~df["group"].isin(gl).to_numpy()
        d2 = df.loc[keep].reset_index(drop=True)
        r = fit(d2, {"A": dict(make_model=mk_gb(fa), X=d2[fa].to_numpy(float)),
                     "AB": dict(make_model=mk_gb(fa + fb), X=d2[fa + fb].to_numpy(float))})
        b = d2["B1"].to_numpy()
        pa, pab = perfold(d2, r["A"][1], b), perfold(d2, r["AB"][1], b)
        res[name] = {"n_pos_removed": {c: int(y[~keep & (sp == c)].sum()) for c in SPECIES},
                     "A": pa, "AB": pab,
                     "human_dAP_AB_vs_B1": pab["human"]["delta"], "human_dAP_AB_minus_A": pab["human"]["delta"] - pa["human"]["delta"],
                     "T_AB": pab["T"]}
        print(name, res[name]["n_pos_removed"], round(pab["human"]["delta"], 3), round(pa["human"]["delta"], 3), round(pab["T"], 3))
    out["drop_groups"] = res

    # (b) LOSO
    lo = {}
    b = df["B1"].to_numpy()
    for tr_sp, te_sp in (("yeast", "human"), ("human", "yeast")):
        specs = {k: dict(make_model=mk_gb(c), X=df[c].to_numpy(float), train_mask=sp == tr_sp, test_mask=sp == te_sp,
                         weights=None) for k, c in (("A", fa), ("AB", fa + fb), ("B", fb))}
        r = fit(df, specs)
        lo[f"{tr_sp}_to_{te_sp}"] = {k: perfold(df, v[1], b, species=(te_sp,))[te_sp] for k, v in r.items()}
        lo[f"{tr_sp}_to_{te_sp}"]["pooled_ap"] = {k: pooled_ap(y[sp == te_sp], v[0][sp == te_sp]) for k, v in r.items()}
        lo[f"{tr_sp}_to_{te_sp}"]["pooled_ap"]["B1"] = pooled_ap(y[sp == te_sp], b[sp == te_sp])
        print("LOSO", tr_sp, te_sp, {k: round(v[te_sp]["delta"] if False else lo[f"{tr_sp}_to_{te_sp}"][k]["delta"], 3) for k in r})
    out["loso"] = lo

    # (c) within-ascertainment cells (human)
    raw = load_labels(ROOT / "data" / "labels" / "dev" / "human_dev.tsv", low_memory=False)
    a_, b_ = raw["gene_a"].astype(str), raw["gene_b"].astype(str)
    raw["pair_id"] = "human:" + np.where(a_ <= b_, a_, b_) + "|" + np.where(a_ <= b_, b_, a_)
    hm = sp == "human"
    h = df.loc[hm, ["pair_id", "group", "label", "species"]].merge(raw[["pair_id", "studies"]], on="pair_id", how="left")
    core = np.load(OUT / "oof_core.npz")
    scores = {"B1": b[hm], "A": core["A"].mean(0)[hm], "AB": core["AB"].mean(0)[hm], "B": core["B"].mean(0)[hm]}
    h["cell"] = h["studies"].fillna("").astype(str).apply(lambda s: ",".join(sorted(s.split(","))))
    cells = h.groupby("cell")["label"].agg(["sum", "size"])
    ok = cells[(cells["sum"] > 0) & (cells["sum"] < cells["size"])].index
    hh = h[h["cell"].isin(ok)].reset_index(drop=True)
    idx = h["cell"].isin(ok).to_numpy()

    def cell_ap(w=None):
        res_ = {}
        for k, s in scores.items():
            s = s[idx]
            num, den = 0.0, 0.0
            for c, ii in hh.groupby("cell").indices.items():
                yy = hh["label"].to_numpy()[ii]
                ww = np.ones(len(ii)) if w is None else w[ii]
                if (yy * ww).sum() == 0 or ((1 - yy) * ww).sum() == 0:
                    continue
                from sklearn.metrics import average_precision_score
                apv = average_precision_score(yy, s[ii], sample_weight=ww)
                num += apv * (yy * ww).sum()
                den += (yy * ww).sum()
            res_[k] = num / den if den else np.nan
        return res_
    pt = cell_ap()
    bs = make_bootstrap(hh, 1000, 11)
    gi = bs.pair_index(hh["group"].to_numpy())
    boots = {k: [] for k in scores}
    for r_ in range(bs.B):
        w = bs.counts[r_][gi].astype(float)
        v = cell_ap(w)
        for k in scores:
            boots[k].append(v[k])
    boots = {k: np.array(v) for k, v in boots.items()}
    out["within_study_set_cells"] = {
        "n_cells_used": int(len(ok)), "n_pairs": int(len(hh)), "n_pos": int(hh["label"].sum()),
        "n_pos_excluded_cells": int(h.loc[~h["cell"].isin(ok), "label"].sum()),
        "pos_weighted_ap": pt,
        "delta_AB_vs_B1": {"point": pt["AB"] - pt["B1"], **ci(boots["AB"] - boots["B1"])},
        "delta_AB_vs_A": {"point": pt["AB"] - pt["A"], **ci(boots["AB"] - boots["A"])},
        "delta_A_vs_B1": {"point": pt["A"] - pt["B1"], **ci(boots["A"] - boots["B1"])},
        "delta_B_vs_B1": {"point": pt["B"] - pt["B1"], **ci(boots["B"] - boots["B1"])},
        "bootstrap_B": 1000}
    print("cells", out["within_study_set_cells"])
    dump(out, "concentration.json")


if __name__ == "__main__":
    main()
