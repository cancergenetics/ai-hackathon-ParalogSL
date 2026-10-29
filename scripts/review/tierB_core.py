"""Review step 1: reproduce GBDT A / A+B T_pf, then controls (Tier B rows shuffled within species, labels
shuffled within species, Tier B alone) and the per-seed spread. Output results/review/tierB_gain/core.json and
oof_core.npz (per-seed OOF matrices, reused by later review scripts)."""

from tierB_common import (SEEDS, SPECIES, OUT, dump, fit, load, mk_gb, mk_lr, perfold, perfold_seed, pooled_ap)

import numpy as np


def main():
    df, fa, fb, _ = load()
    y = df["label"].to_numpy()
    sp = df["species"].to_numpy()
    base = df["B1"].to_numpy()
    XA, XB = df[fa].to_numpy(float), df[fb].to_numpy(float)
    XAB = np.column_stack([XA, XB])
    cab = fa + fb
    specs = {"A": dict(make_model=mk_gb(fa), X=XA), "AB": dict(make_model=mk_gb(cab), X=XAB),
             "B": dict(make_model=mk_gb(fb), X=XB),
             "LR_A": dict(make_model=mk_lr(fa), X=XA), "LR_AB": dict(make_model=mk_lr(cab), X=XAB)}
    res = fit(df, specs)
    out = {"n_features": {"A": len(fa), "B": len(fb)}, "tierB_features": fb,
           "counts": {s: {"n": int((sp == s).sum()), "pos": int(y[sp == s].sum()),
                          "prevalence": float(y[sp == s].mean())} for s in SPECIES},
           "params": "GBDT_17 == GBDT_AB_17 {l2 1.0, lr 0.03, iter 100, leaves 7, min_leaf 50}; LR C=0.01"}
    out["perfold"] = {k: perfold(df, v[1], base) for k, v in res.items()}
    out["perfold_per_seed"] = {k: perfold_seed(df, res[k][1], base) for k in ("A", "AB", "B")}
    out["pooled_ap"] = {k: {s: pooled_ap(y[sp == s], v[0][sp == s]) for s in SPECIES} for k, v in res.items()}
    out["pooled_ap"]["B1"] = {s: pooled_ap(y[sp == s], base[sp == s]) for s in SPECIES}

    # --- control: Tier B rows shuffled within species (one permutation per seed, that seed's folds)
    shuf_rows, lab_rows = [], []
    for s in SEEDS:
        rng = np.random.default_rng(1000 + s)
        XBs = XB.copy()
        for c in SPECIES:
            idx = np.flatnonzero(sp == c)
            XBs[idx] = XB[idx[rng.permutation(len(idx))]]
        r = fit(df, {"ABshufB": dict(make_model=mk_gb(cab), X=np.column_stack([XA, XBs]))}, seeds=[s])
        shuf_rows.append(r["ABshufB"][1][0])
        # --- control: labels shuffled within species
        yp = y.copy()
        rng2 = np.random.default_rng(2000 + s)
        for c in SPECIES:
            idx = np.flatnonzero(sp == c)
            yp[idx] = y[idx[rng2.permutation(len(idx))]]
        r2 = fit(df, {"ABshufY": dict(make_model=mk_gb(cab), X=XAB, target=yp)}, seeds=[s])
        sc = r2["ABshufY"][0]
        lab_rows.append({"seed": s, **{f"ap_{c}": pooled_ap(yp[sp == c], sc[sp == c]) for c in SPECIES},
                         "perfold_vs_B1_on_permuted": perfold(df, r2["ABshufY"][1], base, [s], y=yp),
                         "B1_ap_on_permuted": {c: pooled_ap(yp[sp == c], base[sp == c]) for c in SPECIES}})
    shuf = np.vstack(shuf_rows)
    out["control_shuffle_tierB_rows"] = {"perfold": perfold(df, shuf, base),
                                         "per_seed": perfold_seed(df, shuf, base),
                                         "pooled_ap": {c: pooled_ap(y[sp == c], shuf.mean(0)[sp == c]) for c in SPECIES}}
    out["control_shuffle_labels"] = lab_rows
    d_ab_a = {s: out["perfold_per_seed"]["AB"][s]["T"] - out["perfold_per_seed"]["A"][s]["T"] for s in SEEDS}
    out["seed_spread"] = {"T_AB": [out["perfold_per_seed"]["AB"][s]["T"] for s in SEEDS],
                          "T_A": [out["perfold_per_seed"]["A"][s]["T"] for s in SEEDS],
                          "dT_AB_minus_A": list(d_ab_a.values()),
                          "human_delta_AB": [out["perfold_per_seed"]["AB"][s]["human"]["delta"] for s in SEEDS],
                          "yeast_delta_AB": [out["perfold_per_seed"]["AB"][s]["yeast"]["delta"] for s in SEEDS]}
    np.savez_compressed(OUT / "oof_core.npz", pair_id=df["pair_id"].to_numpy().astype(str),
                        A=res["A"][1], AB=res["AB"][1], B=res["B"][1], ABshufB=shuf)
    dump(out, "core.json")
    for k in ("A", "AB", "B", "LR_A", "LR_AB"):
        p = out["perfold"][k]
        print(k, round(p["T"], 4), round(p["human"]["delta"], 4), round(p["yeast"]["delta"], 4))
    p = out["control_shuffle_tierB_rows"]["perfold"]
    print("ABshufB", round(p["T"], 4), round(p["human"]["delta"], 4), round(p["yeast"]["delta"], 4))
    print("seed spread", out["seed_spread"])
    print("label shuffle", [(r["ap_human"], r["ap_yeast"]) for r in lab_rows])


if __name__ == "__main__":
    main()
