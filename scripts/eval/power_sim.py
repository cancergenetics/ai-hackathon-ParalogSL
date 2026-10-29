#!/usr/bin/env python
"""M1 power simulation of the full primary decision rule (DESIGN.md v3 §6.1 "Power", §6.2; config/mde.yaml).

Allowlisted (D4/D5) to read sealed split candidates and label counts. For each candidate test fraction, the
real candidate test set (one deterministic draw per fraction, data/labels/sealed/split_candidates.tsv) gives
the primary labels, the split-group membership and the group sizes of every primary-labelled test pair.
Only the scores are simulated (binormal, B1 vs model, within-class correlation rho, one shared random offset
per split group). **No score of any kind (identity included) is computed on any labelled pair.**

Decision rule simulated (config/success.yaml, DESIGN §6.2): T = mean over {human, yeast} of
dAP_s = AP(model) - AP(B1); pass iff the percentile 95% CI of T from the stratified paired group bootstrap
(strata: species composition x group has a primary positive) has lower bound > 0 AND dAP_s > 0 in both
species. The model AUROC in each species is calibrated so that E[dAP_s] = t_target_ap at the real composition.

Chosen fraction = smallest candidate with power >= power_target at the deciding rho; otherwise 0.40 and the
primary comparison is pre-declared underpowered. Writes results/mde.json.

    .venv/bin/python scripts/eval/power_sim.py [--workers 8] [--out-root DIR]
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import yaml
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
MDE_CFG = ROOT / "config" / "mde.yaml"
SPECIES = ("human", "yeast")


def _load_make_splits():
    spec = importlib.util.spec_from_file_location("make_splits", ROOT / "scripts" / "labels" / "make_splits.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ----------------------------------------------------------------------------- weighted AP


def weighted_ap(y_sorted: np.ndarray, w_sorted: np.ndarray) -> np.ndarray:
    """AP of pairs already sorted by decreasing score, each pair carrying weight w (multiplicity).

    Equal to sklearn.metrics.average_precision_score on the data with every pair repeated w times when scores
    have no ties across distinct pairs (tests/test_splits.py checks this). ``w_sorted`` may be 2-D (B x n);
    returns nan where the weighted positive count is 0."""
    w = np.asarray(w_sorted, dtype=np.float64)
    wy = w * y_sorted
    cw = np.cumsum(w, axis=-1)
    ctp = np.cumsum(wy, axis=-1)
    prec = np.divide(ctp, cw, out=np.zeros_like(ctp), where=cw > 0)
    num = np.sum(wy * prec, axis=-1)
    den = ctp[..., -1]
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)


# ----------------------------------------------------------------------------- composition


def composition(cand, labels, f: float, ms) -> dict:
    """Primary-labelled test pairs of the candidate draw at fraction f: y, species code, group index, and the
    bootstrap stratum of every test group."""
    g2g = ms.gene_to_group(cand)
    test_groups = set(cand.loc[cand[ms.frac_key(f)] == 1, "group_id"])
    ys, sps, grs = [], [], []
    for k, sp in enumerate(SPECIES):
        df = labels[sp]
        df = df[df.y_primary >= 0]
        grp = df["key_a"].map(g2g)
        m = grp.isin(test_groups).to_numpy()
        ys.append(df["y_primary"].to_numpy()[m])
        sps.append(np.full(int(m.sum()), k))
        grs.append(grp.to_numpy()[m])
    y = np.concatenate(ys).astype(np.float64)
    sp = np.concatenate(sps)
    gnames = np.concatenate(grs)
    uniq, gidx = np.unique(gnames, return_inverse=True)
    G = len(uniq)
    has_h = np.bincount(gidx, weights=(sp == 0), minlength=G) > 0
    has_y = np.bincount(gidx, weights=(sp == 1), minlength=G) > 0
    has_pos = np.bincount(gidx, weights=y, minlength=G) > 0
    comp = np.where(has_h & has_y, 2, np.where(has_h, 0, 1))
    stratum = comp * 2 + has_pos.astype(int)
    return {"y": y, "sp": sp, "g": gidx, "n_groups": G, "stratum": stratum,
            "counts": {s: {"pos": int(y[sp == k].sum()), "neg": int((y[sp == k] == 0).sum())}
                       for k, s in enumerate(SPECIES)},
            "strata_counts": {name: int((stratum == v).sum()) for v, name in enumerate(
                ["human_only_nopos", "human_only_pos", "yeast_only_nopos", "yeast_only_pos",
                 "mixed_nopos", "mixed_pos"])}}


# ----------------------------------------------------------------------------- simulation


def mu_of_auc(auc: float, sd_c: float) -> float:
    return float(np.sqrt(2.0 * (1.0 + sd_c ** 2)) * norm.ppf(auc))


def auc_of_mu(mu: float, sd_c: float) -> float:
    return float(norm.cdf(mu / np.sqrt(2.0 * (1.0 + sd_c ** 2))))


def draw_scores(rng, comp, mu1: np.ndarray, mu2: np.ndarray, rho: float, sd_c: float):
    n = len(comp["y"])
    z1 = rng.standard_normal(n)
    z2 = rho * z1 + np.sqrt(1.0 - rho ** 2) * rng.standard_normal(n)
    u = rng.normal(0.0, sd_c, comp["n_groups"])[comp["g"]]
    y, sp = comp["y"], comp["sp"]
    s1 = z1 + mu1[sp] * y + u
    s2 = z2 + mu2[sp] * y + u
    return s1, s2


def point_ap(y, s) -> float:
    o = np.argsort(-s, kind="mergesort")
    return float(weighted_ap(y[o], np.ones(len(y))))


def calibrate(comp, mu1: np.ndarray, cfg: dict, seed_seq) -> tuple[np.ndarray, dict]:
    """Per species, the model mean shift mu2 with E[dAP_s] = t_target_ap (bisection, common random numbers)."""
    c = cfg["calibration"]
    sd_c = float(cfg["cluster_score_sd"])
    T = float(cfg["t_target_ap"])
    rng = np.random.default_rng(seed_seq)
    n = len(comp["y"])
    Z1 = rng.standard_normal((int(c["n_sim"]), n))
    Z2 = rng.standard_normal((int(c["n_sim"]), n))
    U = rng.normal(0.0, sd_c, (int(c["n_sim"]), comp["n_groups"]))[:, comp["g"]]
    y, sp = comp["y"], comp["sp"]
    mu2 = mu1.copy()
    info = {}
    for k, s in enumerate(SPECIES):
        m = sp == k
        yk = y[m]
        base = np.array([point_ap(yk, Z1[i, m] + mu1[k] * yk + U[i, m]) for i in range(len(Z1))])

        def mean_dap(mu):
            return float(np.mean([point_ap(yk, Z2[i, m] + mu * yk + U[i, m]) for i in range(len(Z2))]) - base.mean())

        lo, hi = mu1[k], mu_of_auc(float(c["max_auroc"]), sd_c)
        reachable = mean_dap(hi) >= T
        if reachable:
            for _ in range(int(c["iterations"])):
                mid = 0.5 * (lo + hi)
                if mean_dap(mid) < T:
                    lo = mid
                else:
                    hi = mid
            mu2[k] = 0.5 * (lo + hi)
        else:
            mu2[k] = hi
        info[s] = {"b1_auroc": auc_of_mu(mu1[k], sd_c), "model_auroc": auc_of_mu(mu2[k], sd_c),
                   "delta_auroc": auc_of_mu(mu2[k], sd_c) - auc_of_mu(mu1[k], sd_c),
                   "expected_b1_ap": float(base.mean()), "expected_delta_ap": mean_dap(mu2[k]),
                   "target_reachable": bool(reachable)}
    return mu2, info


def bootstrap_counts(rng, stratum: np.ndarray, B: int) -> np.ndarray:
    """B x n_groups multiplicities; groups resampled with replacement within strata."""
    G = len(stratum)
    W = np.zeros((B, G), dtype=np.int32)
    for v in np.unique(stratum):
        members = np.where(stratum == v)[0]
        k = len(members)
        idx = rng.integers(0, k, (B, k))
        cnt = np.bincount((idx + np.arange(B)[:, None] * k).ravel(), minlength=B * k).reshape(B, k)
        W[:, members] = cnt
    return W


def simulate(args) -> dict:
    comp, mu1, mu2, rho, cfg, seed_seq, rel_gain = args
    sd_c = float(cfg["cluster_score_sd"])
    B = int(cfg["bootstrap_resamples"])
    alpha = 1.0 - float(cfg["ci_level"])
    rng = np.random.default_rng(seed_seq)
    y, sp = comp["y"], comp["sp"]
    rows = []
    dropped = 0
    for _ in range(int(cfg["n_sim"])):
        s1, s2 = draw_scores(rng, comp, mu1, mu2, rho, sd_c)
        Wg = bootstrap_counts(rng, comp["stratum"], B)
        d_point, rel, dboot = [], [], []
        for k in range(len(SPECIES)):
            m = np.where(sp == k)[0]
            yk = y[m]
            o1 = np.argsort(-s1[m], kind="mergesort")
            o2 = np.argsort(-s2[m], kind="mergesort")
            a1 = float(weighted_ap(yk[o1], np.ones(len(m))))
            a2 = float(weighted_ap(yk[o2], np.ones(len(m))))
            d_point.append(a2 - a1)
            rel.append(a2 / a1 - 1.0)
            Wk = Wg[:, comp["g"][m]]
            dboot.append(weighted_ap(yk[o2], Wk[:, o2]) - weighted_ap(yk[o1], Wk[:, o1]))
        Tb = np.mean(np.vstack(dboot), axis=0)
        ok = np.isfinite(Tb)
        dropped += int((~ok).sum())
        lo, hi = np.percentile(Tb[ok], [100 * alpha / 2, 100 * (1 - alpha / 2)])
        rule = bool(lo > 0 and min(d_point) > 0)
        rows.append((d_point[0], d_point[1], float(np.mean(d_point)), bool(lo > 0), rule,
                     bool(rule and np.mean(rel) >= rel_gain), float(hi - lo)))
    r = np.array(rows, dtype=float)
    return {"rho": rho, "n_sim": int(cfg["n_sim"]), "bootstrap_resamples": B,
            "mean_delta_ap_human": float(r[:, 0].mean()), "mean_delta_ap_yeast": float(r[:, 1].mean()),
            "mean_T": float(r[:, 2].mean()), "sd_T": float(r[:, 2].std(ddof=1)),
            "p_ci_excludes_0": float(r[:, 3].mean()), "power_rule": float(r[:, 4].mean()),
            "power_rule_and_relative_gain": float(r[:, 5].mean()),
            "mean_ci_width": float(r[:, 6].mean()), "mde80_approx_2.8sd": float(2.8 * r[:, 2].std(ddof=1)),
            "dropped_resamples_total": dropped}


# ----------------------------------------------------------------------------- main


def run(out_root: Path, workers: int) -> dict:
    ms = _load_make_splits()
    cfg = yaml.safe_load(MDE_CFG.read_text())
    succ = yaml.safe_load((ROOT / "config" / "success.yaml").read_text())
    scfg = ms.load_cfg()
    # constants must agree with config/success.yaml and config/splits.yaml
    assert float(cfg["t_target_ap"]) == float(succ["t_target_ap"]), "t_target_ap differs from success.yaml"
    assert float(cfg["power_target"]) == float(succ["power_target"]), "power_target differs from success.yaml"
    fr = [float(x) for x in cfg["candidate_fractions"]]
    assert fr == [float(x) for x in scfg["test_assignment"]["candidate_fractions"]], "fractions differ from splits.yaml"
    lo_r, hi_r = succ["test_fraction_range"]
    assert all(lo_r <= f <= hi_r for f in fr)
    rel_gain = float(succ["relative_gain"])

    P = ms.paths(scfg, out_root)
    cand = ms.read_candidates(scfg, out_root)
    labels = ms.load_all_labels(scfg)
    sd_c = float(cfg["cluster_score_sd"])
    n_pos_yeast = int((labels["yeast"].y_primary == 1).sum())
    yeast_fallback = n_pos_yeast < int(cfg["yeast_fallback_min_primary_positives"])

    root_ss = np.random.SeedSequence(int(cfg["seed"]))
    scen_specs = [("decisive", float(cfg["rho"]), cfg["b1_auroc"])]
    scen_specs += [(f"rho_{r}", float(r), cfg["b1_auroc"]) for r in cfg["rho_report"]]
    scen_specs += [("b1_auroc_report", float(cfg["rho"]), cfg["b1_auroc_report"])]
    out = {"script": "scripts/eval/power_sim.py", "config": "config/mde.yaml",
           "config_sha256": sha256_file(MDE_CFG),
           "candidates_sha256": sha256_file(P["candidates"]),
           "inputs": {f"labels_{sp}_sha256": sha256_file(ROOT / scfg["labels"][sp]["file"]) for sp in SPECIES},
           "note": "scores are simulated; no score of any kind was computed on any labelled pair",
           "recall_assumption": "every primary-labelled test pair is scored (caller recall 1.0)",
           "yeast_primary_positives": n_pos_yeast, "yeast_fallback": yeast_fallback,
           "fractions": {}}
    jobs, keys = [], []
    calib = {}
    for fi, f in enumerate(fr):
        comp = composition(cand, labels, f, ms)
        fkey = f"{f:.2f}"
        out["fractions"][fkey] = {"test_counts": comp["counts"], "n_test_groups_with_primary_pairs": comp["n_groups"],
                                  "strata_counts": comp["strata_counts"], "calibration": {}, "scenarios": {}}
        for si, (name, rho, aucs) in enumerate(scen_specs):
            akey = f"{aucs['human']}_{aucs['yeast']}"
            mu1 = np.array([mu_of_auc(float(aucs[s]), sd_c) for s in SPECIES])
            if (fkey, akey) not in calib:
                ss = np.random.SeedSequence([int(cfg["seed"]), fi, 1000 + len(calib)])
                calib[(fkey, akey)] = calibrate(comp, mu1, cfg, ss)
                out["fractions"][fkey]["calibration"][akey] = calib[(fkey, akey)][1]
            mu2 = calib[(fkey, akey)][0]
            ss = np.random.SeedSequence([int(cfg["seed"]), fi, si])
            jobs.append((comp, mu1, mu2, rho, cfg, ss, rel_gain))
            keys.append((fkey, name, akey))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        results = list(ex.map(simulate, jobs))
    for (fkey, name, akey), res in zip(keys, results):
        res["b1_auroc"] = akey
        out["fractions"][fkey]["scenarios"][name] = res
    chosen, underpowered = None, True
    for f in fr:
        if out["fractions"][f"{f:.2f}"]["scenarios"]["decisive"]["power_rule"] >= float(cfg["power_target"]):
            chosen, underpowered = f, False
            break
    if chosen is None:
        chosen = max(fr)
    out["decision"] = {"rule": "smallest fraction with power_rule >= power_target at rho (decisive scenario); "
                               "else the largest fraction, pre-declared underpowered",
                       "power_target": float(cfg["power_target"]), "rho": float(cfg["rho"])}
    out["chosen_fraction"] = chosen
    out["power_at_chosen"] = out["fractions"][f"{chosen:.2f}"]["scenarios"]["decisive"]["power_rule"]
    out["underpowered"] = underpowered
    path = out_root / cfg["output"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2) + "\n")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-root", type=Path, default=ROOT)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args(argv)
    res = run(args.out_root.resolve(), args.workers)
    print(json.dumps({"chosen_fraction": res["chosen_fraction"], "power_at_chosen": res["power_at_chosen"],
                      "underpowered": res["underpowered"], "yeast_fallback": res["yeast_fallback"]}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
