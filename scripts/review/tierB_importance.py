"""Review step 5: which Tier B features drive the A+B gain, and does their direction match the priors?

(a) Held-out permutation importance: for each seed x fold, GBDT A+B (GBDT_AB_17) is fitted on the training folds;
    each feature (and the blocks 'all Tier B', 'all Tier A') is permuted within the test fold within species
    (one permutation per seed x fold, rng seed = 100*seed + fold); importance = drop in per-species AP.
(b) Direction: for the top Tier B features, standardized mean difference positives vs negatives (human all,
    human n_studies_testable>=3, yeast) and Spearman with human n_studies_testable (ascertainment link).
(c) Supplementary ascertainment: within the >=3 stratum, grouped-CV LR (stored folds) label ~ ASC + B1 vs
    ASC + B1 + logit(AB); within study-set cells with >= 20 pairs and >= 2 positives, AP of B1 / A / AB and of the
    yeast-trained (LOSO) A / AB scores.
Output results/review/tierB_gain/importance.json"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score as aps

from tierB_common import GB_PARAMS, OUT, ROOT, SEEDS, SPECIES, dump, fit, load, mk_gb, perfold
from slparalog.eval.cv import fold_col, oof_many
from slparalog.eval.models import GBDTModel, LRModel
from slparalog.eval.stats import species_balanced_weights
from slparalog.io.labels import load_labels
from joblib import Parallel, delayed

LIBS = ["chymera2020", "parrish2021", "in4mer2024", "thompson2021", "ito2021", "klingbeil2024", "dede2020",
        "chan2025", "harle2025"]


def one(X, y, sp, f, k, s, cols, blocks):
    tr, te = f != k, f == k
    w = species_balanced_weights(sp[tr], y[tr])
    m = GBDTModel(cols, dict(GB_PARAMS), seed=s).fit(X[tr], y[tr], sample_weight=w)
    Xt = X[te]
    yt, spt = y[te], sp[te]
    base = m.predict(Xt)
    rng = np.random.default_rng(100 * s + int(k))
    res = {}
    for c in SPECIES:
        mc = spt == c
        if yt[mc].sum() == 0:
            continue
        b_ap = aps(yt[mc], base[mc])
        idx = np.flatnonzero(mc)
        perm = idx[rng.permutation(len(idx))]
        for name, js in blocks.items():
            Xp = Xt.copy()
            Xp[np.ix_(idx, js)] = Xt[np.ix_(perm, js)]
            res[(c, name)] = b_ap - aps(yt[mc], m.predict(Xp)[mc])
    return res


def main():
    df, fa, fb, _ = load()
    y = df["label"].to_numpy()
    sp = df["species"].to_numpy()
    cols = fa + fb
    X = df[cols].to_numpy(float)
    blocks = {c: [j] for j, c in enumerate(cols)}
    blocks["ALL_TIER_B"] = list(range(len(fa), len(cols)))
    blocks["ALL_TIER_A"] = list(range(len(fa)))
    jobs = [(X, y, sp, df[fold_col(s)].to_numpy(), k, s, cols, blocks) for s in SEEDS
            for k in np.unique(df[fold_col(s)])]
    rs = Parallel(n_jobs=6, backend="loky")(delayed(one)(*j) for j in jobs)
    imp = {}
    for c in SPECIES:
        for name in blocks:
            v = [r[(c, name)] for r in rs if (c, name) in r]
            imp.setdefault(c, {})[name] = {"mean_drop_ap": float(np.mean(v)), "sd": float(np.std(v, ddof=1)),
                                           "n": len(v)}
    out = {"importance": imp}
    for c in SPECIES:
        top = sorted(((n, v["mean_drop_ap"]) for n, v in imp[c].items()), key=lambda t: -t[1])[:15]
        print(c, [(n, round(v, 4)) for n, v in top])

    # (b) direction
    nst = pd.to_numeric(df["meta_n_studies_testable"], errors="coerce").to_numpy()
    hm = sp == "human"
    direc = {}
    for c in fb:
        x = df[c].to_numpy(float)

        def smd(m):
            a_, b_ = x[m & (y == 1)], x[m & (y == 0)]
            a_, b_ = a_[np.isfinite(a_)], b_[np.isfinite(b_)]
            sd = np.sqrt((a_.var() + b_.var()) / 2)
            return float((a_.mean() - b_.mean()) / sd) if sd > 0 else float("nan")
        okh = hm & np.isfinite(x)
        direc[c] = {"smd_human": smd(hm), "smd_human_ge3": smd(hm & (nst >= 3)), "smd_human_1": smd(hm & (nst == 1)),
                    "smd_yeast": smd(sp == "yeast"),
                    "spearman_with_n_studies_human": float(spearmanr(x[okh], nst[okh]).statistic),
                    "drop_ap_human": imp["human"][c]["mean_drop_ap"], "drop_ap_yeast": imp["yeast"][c]["mean_drop_ap"]}
    out["direction"] = direc
    print(pd.DataFrame(direc).T.sort_values("drop_ap_human", ascending=False).round(3).head(12).to_string())

    # (c) supplementary ascertainment
    core = np.load(OUT / "oof_core.npz")
    raw = load_labels(ROOT / "data" / "labels" / "dev" / "human_dev.tsv", low_memory=False)
    a_, b_ = raw["gene_a"].astype(str), raw["gene_b"].astype(str)
    raw["pair_id"] = "human:" + np.where(a_ <= b_, a_, b_) + "|" + np.where(a_ <= b_, b_, a_)
    h = df.loc[hm].reset_index(drop=True)
    h = h.merge(raw[["pair_id", "studies"]], on="pair_id", how="left")
    lib = np.column_stack([h["studies"].fillna("").str.split(",").apply(lambda z: L in z).to_numpy() for L in LIBS])
    nh = nst[hm]
    ASC = np.column_stack([nh == 1, nh == 2, np.minimum(nh, 9), lib]).astype(float)

    def lgt(p):
        p = np.clip(p, 1e-6, 1 - 1e-6)
        return np.log(p / (1 - p))
    ab_h, a_h = core["AB"].mean(0)[hm], core["A"].mean(0)[hm]
    # LOSO yeast-trained scores for human
    r = fit(df, {k: dict(make_model=mk_gb(cc), X=df[cc].to_numpy(float), train_mask=sp == "yeast",
                         test_mask=hm, weights=None) for k, cc in (("A_y", fa), ("AB_y", fa + fb))})
    ay_h, aby_h = r["A_y"][0][hm], r["AB_y"][0][hm]
    g3 = nh >= 3
    h3 = h.loc[g3].reset_index(drop=True)
    sets = {"ASC+B1": np.column_stack([ASC, h["B1"]])[g3],
            "ASC+B1+AB": np.column_stack([ASC, h["B1"], lgt(ab_h)])[g3],
            "ASC+B1+AB_yeast_trained": np.column_stack([ASC, h["B1"], lgt(aby_h)])[g3]}
    specs = {k: dict(make_model=(lambda cc: (lambda s: LRModel(cc, C=1.0, seed=s)))([f"x{i}" for i in range(v.shape[1])]),
                     X=v, weights=None) for k, v in sets.items()}
    rr = oof_many(h3, specs, SEEDS, n_jobs=6)
    out["ge3_stratum_asc_stack"] = {k: perfold(h3, rr[k][1], h3["B1"].to_numpy(), species=("human",))["human"]
                                    for k in sets}
    out["ge3_stratum_asc_stack"]["n"] = int(len(h3))
    out["ge3_stratum_asc_stack"]["n_pos"] = int(h3["label"].sum())
    print("ge3 stack", {k: v for k, v in out["ge3_stratum_asc_stack"].items()})
    # all-human stack for the yeast-trained AB score
    hsets = {"ASC+B1": np.column_stack([ASC, h["B1"]]), "ASC+B1+AB_yeast_trained": np.column_stack([ASC, h["B1"], lgt(aby_h)]),
             "ASC+B1+A_yeast_trained": np.column_stack([ASC, h["B1"], lgt(ay_h)])}
    specs = {k: dict(make_model=(lambda cc: (lambda s: LRModel(cc, C=1.0, seed=s)))([f"x{i}" for i in range(v.shape[1])]),
                     X=v, weights=None) for k, v in hsets.items()}
    rr = oof_many(h, specs, SEEDS, n_jobs=6)
    out["human_asc_stack_yeast_trained"] = {k: perfold(h, rr[k][1], h["B1"].to_numpy(), species=("human",))["human"]
                                            for k in hsets}
    print("human stack yeast-trained", out["human_asc_stack_yeast_trained"])
    # large cells
    h["cell"] = h["studies"].fillna("").apply(lambda z: ",".join(sorted(z.split(","))))
    sc = {"B1": h["B1"].to_numpy(), "A": a_h, "AB": ab_h, "A_yeast_trained": ay_h, "AB_yeast_trained": aby_h}
    cells = []
    for cname, ii in h.groupby("cell").indices.items():
        yy = h["label"].to_numpy()[ii]
        if len(ii) >= 20 and yy.sum() >= 2:
            cells.append({"cell": cname, "n": int(len(ii)), "pos": int(yy.sum()),
                          **{k: float(aps(yy, v[ii])) for k, v in sc.items()}})
    ct = pd.DataFrame(cells)
    wpos = ct["pos"] / ct["pos"].sum()
    out["large_cells"] = {"cells": cells, "pos_weighted_ap": {k: float((ct[k] * wpos).sum()) for k in sc},
                          "n_cells": len(cells), "n_pos": int(ct["pos"].sum()),
                          "n_cells_AB_gt_B1": int((ct["AB"] > ct["B1"]).sum())}
    print(ct.round(3).to_string())
    print("large cells", out["large_cells"]["pos_weighted_ap"])
    dump(out, "importance.json")


if __name__ == "__main__":
    main()
