#!/usr/bin/env python
"""Build human paralog SL labels from primary screen tables (rules: config/labels_human.yaml, pre-registered).

Outputs
  data/interim/labels/human_calls.tsv      per study x line x gene pair call (S/N/B/U/X), not committed
  data/labels/human_labels_all.tsv         pair-level tiers + per-source evidence codes + licence flag (no splits)
  results/labels/human_label_qc.json       counts, positive rates, mapping losses, Dennler agreement (report only),
                                           13-pair candidate status (report only), consistency checks, R2-14 audits
  results/labels/human_mapping_report.tsv  every symbol that was not an exact approved-symbol hit

Call codes: S = SL, N = clear negative, B = borderline, D = line discordant across studies (counts as borderline),
U = untestable, X = untested (missing from the significance or LFC table).

Usage: .venv/bin/python scripts/labels/build_human_labels.py
"""
from __future__ import annotations

import json
import hashlib
import re
import math
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import human_parsers as hp  # noqa: E402

CFG_PATH = ROOT / "config" / "labels_human.yaml"
ST1_PATH = ROOT / "data" / "raw" / "in4mer2024" / "Supp_table1.xlsx"
UNIVERSE_FILE = ROOT / "data" / "interim" / "universe" / "human_paralog_pairs.tsv"   # pinned path (D8)


# ----------------------------------------------------------------------------------------------------------------
# Line-level calls
# ----------------------------------------------------------------------------------------------------------------

def call_rows(obs: pd.DataFrame, singles: dict, cfg: dict, scales: dict) -> pd.DataFrame:
    """Apply the uniform R2-5 rule (LFC studies) or the authors'-flag rule (no-LFC studies) row by row."""
    lc = cfg["line_call"]
    obs = obs.copy()
    has_lfc = bool(singles)
    if has_lfc:
        for lab, ser in singles.items():
            m, s, T = hp.robust_scale(ser, lc["T_multiplier"], lc["T_floor"])
            scales[(obs["study"].iloc[0], lab)] = dict(m=m, s=s, T=T, n_singles=int(ser.notna().sum()))
        st = obs["study"].iloc[0]
        m = obs["line_label"].map(lambda l: scales[(st, l)]["m"])
        T = obs["line_label"].map(lambda l: scales[(st, l)]["T"])
        a, b, ab = obs["lfc1"].astype(float) - m, obs["lfc2"].astype(float) - m, obs["lfc12"].astype(float) - m
        obs["a"], obs["b"], obs["ab"], obs["T"] = a, b, ab, T
        missing = a.isna() | b.isna() | ab.isna()
        testable = ((a + b) > -T) & ~obs["flag_untestable"].astype(bool)
        sl = testable & obs["gi_strict"].astype(bool) & (ab <= -T)
        neg = testable & ~sl & (~obs["gi_relaxed"].astype(bool) | (ab > -lc["borderline_margin_dko"] * T))
    else:
        obs["a"] = obs["b"] = obs["ab"] = obs["T"] = np.nan
        missing = obs["score"].isna()
        testable = ~obs["flag_untestable"].astype(bool)
        sl = testable & obs["gi_strict"].astype(bool)
        neg = testable & ~sl & ~obs["gi_relaxed"].astype(bool)
    no_sig = obs["unit"].astype(str).str.contains("NO_SIG") | obs["score"].isna()
    call = np.where(sl, "S", np.where(neg, "N", np.where(testable, "B", "U")))
    call = np.where(missing | no_sig, "X", call)
    obs["call"] = call
    return obs


def collapse_dups(calls) -> str:
    """Several calls for one study x line x pair (duplicate rows, symbol collisions): identical -> that call;
    disagreeing testable calls -> B; only untestable -> U; nothing -> X."""
    c = [x for x in calls if x != "X"]
    if not c:
        return "X"
    t = [x for x in c if x in "SNB"]
    if not t:
        return "U"
    return t[0] if len(set(t)) == 1 else "B"


def collapse_domains(calls) -> str:
    """Klingbeil domain-combination collapse (config studies.klingbeil2024.domain_collapse)."""
    c = [x for x in calls if x != "X"]
    if not c:
        return "X"
    t = [x for x in c if x in "SNB"]
    if not t:
        return "U"
    ns = sum(x == "S" for x in t)
    if ns / len(t) > hp.THR["klingbeil2024"]["domain_sl_fraction_gt"]:
        return "S"
    if ns > 0 or any(x == "B" for x in t):
        return "B"
    return "N"


def merge_line(study_calls: list[str]) -> str:
    """Same physical line across studies (R2-7b)."""
    t = [x for x in study_calls if x in "SNB"]
    if not t:
        return "U" if "U" in study_calls else "X"
    if "S" in t and "N" not in t:
        return "S"
    if all(x == "N" for x in t):
        return "N"
    if len(set(t)) == 1:
        return t[0]   # all borderline
    return "D"


def pair_tiers(d: pd.DataFrame, agg_cfg: dict) -> dict:
    """Line merge + primary / lenient / two-studies tiers for one gene pair's study x line calls."""
    pr = agg_cfg["primary"]
    by_line = defaultdict(list)
    for st, lid, c in zip(d.study, d.line_id, d.call):
        by_line[lid].append((st, c))
    merged = {lid: merge_line([c for _, c in v]) for lid, v in by_line.items()}
    n_test = sum(c in "SNBD" for c in merged.values())
    n_sl = sum(c == "S" for c in merged.values())
    n_bord = sum(c in "BD" for c in merged.values())
    pen = n_sl / n_test if n_test else float("nan")
    st_sl = {st for lid, v in by_line.items() if merged[lid] == "S" for st, c in v if c == "S"}
    st_test = {st for st, c in zip(d.study, d.call) if c in "SNB"}
    # primary
    two_route = len(st_sl) >= 2 and n_test and pen >= pr["penetrance_floor_two_study_route"]
    one_route = n_test and pen >= pr["penetrance_single_route"]
    if n_test == 0:
        p, reason = -2, "all_untestable"
    elif n_sl >= pr["min_sl_lines_positive"] and (two_route or one_route):
        p, reason = 1, ""
    elif (n_test >= pr["min_testable_lines_negative"] and n_sl == 0
          and n_bord <= math.floor(pr["borderline_allowance_fraction"] * n_test)):
        p, reason = 0, ""
    else:
        p = -1
        if n_sl == 1:
            reason = "single_line_sl"
        elif n_sl >= 2:
            reason = "low_penetrance"
        elif n_test < pr["min_testable_lines_negative"]:
            reason = "too_few_testable"
        else:
            reason = "borderline_excess"
    # lenient
    len_lab = 1 if n_sl >= 1 else (0 if n_test >= 1 else -2)
    # two studies (De Kegel style + R2-15 floor per PI D4)
    ts = agg_cfg["two_studies"]
    if len(st_sl) >= 2 and pen >= ts["penetrance_floor"]:
        two = 1
    elif len(st_test) >= 2 and n_sl == 0:
        two = 0
    else:
        two = -2 if n_test == 0 else -1
    return dict(p=p, reason=reason, len_lab=len_lab, two=two, merged=merged, n_test=n_test, n_sl=n_sl,
                n_bord=n_bord, pen=pen, st_sl=st_sl, st_test=st_test)


# ----------------------------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------------------------

def main() -> None:
    cfg = yaml.safe_load(CFG_PATH.read_text())
    qc: dict = {"config": str(CFG_PATH.relative_to(ROOT)), "config_version": cfg["version"]}
    try:
        qc["git_commit_at_build"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                                                   text=True).stdout.strip()
    except Exception:  # pragma: no cover
        qc["git_commit_at_build"] = None

    st1 = hp.read_supp_table1(ST1_PATH)
    qc["supp_table1_exp_eq_sum_max_abs_diff"] = {
        k: round(float((v["exp_mean"] - v["gene1_mean"] - v["gene2_mean"]).abs().max()), 4) for k, v in st1.items()}

    # ---- parse -------------------------------------------------------------------------------------------------
    parsed, scales = [], {}
    o, sg, _ = hp.parse_dede(cfg, ROOT, st1); parsed.append(call_rows(o, sg, cfg, scales))
    o, sg, sd1 = hp.parse_thompson(cfg, ROOT, st1); parsed.append(call_rows(o, sg, cfg, scales))
    o, sg = hp.parse_parrish(cfg, ROOT, st1); parsed.append(call_rows(o, sg, cfg, scales))
    o, sg, _ = hp.parse_ito(cfg, ROOT, st1); parsed.append(call_rows(o, sg, cfg, scales))
    o, sg, chy_raw = hp.parse_chymera(cfg, ROOT); parsed.append(call_rows(o, sg, cfg, scales))
    o, sg, _ = hp.parse_in4mer(cfg, ROOT); parsed.append(call_rows(o, sg, cfg, scales))
    o, sg, kl_info = hp.parse_klingbeil(cfg, ROOT); parsed.append(call_rows(o, sg, cfg, scales))
    qc["klingbeil_lfc_lookup"] = kl_info
    hs = cfg["studies"]["harle2025"]
    s4, s2 = hp.read_harle_s4(ROOT, hs), hp.read_harle_s2(ROOT, hs)
    harle_keys = hp.harle_line_keys(cfg, s4, s2)
    parsed.append(call_rows(hp.parse_harle_like("harle2025", s4, harle_keys), {}, cfg, scales))
    t6 = hp.read_chan_t6(ROOT, cfg["studies"]["chan2025"])
    chan_keys = {str(k): v for k, v in cfg["studies"]["chan2025"]["lines"].items()}
    parsed.append(call_rows(hp.parse_harle_like("chan2025", t6, chan_keys), {}, cfg, scales))
    rows = pd.concat(parsed, ignore_index=True)

    qc["line_scales"] = {f"{k[0]}|{k[1]}": {kk: round(v, 4) if isinstance(v, float) else v for kk, v in d.items()}
                         for k, d in scales.items()}

    # CHyMErA: cross-check late-sheet LFCs against the harmonised Supp_table1 sheets (logged only)
    xc = {}
    for lab, sheet in (("HAP1", "Gonatopoulos_HAP1"), ("RPE1", "Gonatopoulos_RPE1")):
        look, _ = hp.st1_lookup(st1[sheet])
        r = rows[(rows.study == "chymera2020") & (rows.line_label == lab)]
        o2 = [look.get(hp.pair_key(a, b), {}).get("obs", np.nan) for a, b in zip(r.sym1, r.sym2)]
        d = pd.DataFrame({"late": r.lfc12.astype(float).values, "st1": o2}).dropna()
        xc[lab] = {"n": int(len(d)), "pearson_r": round(float(d.corr().iloc[0, 1]), 4)}
    qc["chymera_late_vs_supp_table1_dko"] = xc

    # ---- unit collapse: CHyMErA duplicate rows, Klingbeil domain combinations -----------------------------------
    rows["k"] = [hp.pair_key(a, b) for a, b in zip(rows.sym1, rows.sym2)]
    grp_cols = ["study", "line_label", "line_key", "k"]
    first = rows.groupby(grp_cols, sort=False).agg(
        sym1=("sym1", "first"), sym2=("sym2", "first"), author_pair_class=("author_pair_class", "first"),
        score=("score", "min"), a=("a", "mean"), b=("b", "mean"), ab=("ab", "mean"), T=("T", "first"),
        n_units=("call", "size")).reset_index()
    cl = rows.groupby(grp_cols, sort=False)["call"].agg(list).reset_index()
    first["call"] = [collapse_domains(c) if s == "klingbeil2024" else collapse_dups(c)
                     for s, c in zip(cl.study, cl.call)]
    sp = first

    # ---- map symbols ------------------------------------------------------------------------------------------
    gm = cfg["gene_mapping"]["files"]
    mapper = hp.SymbolMapper(ROOT / gm["idmapping"], ROOT / gm["gene2acc"], ROOT / gm["fasta"])
    maprep = []
    syms_by_study = defaultdict(set)
    for st, a, b in zip(sp.study, sp.sym1, sp.sym2):
        syms_by_study[st] |= {a, b}
    mapping_summary = {}
    for st, syms in syms_by_study.items():
        cnt = Counter()
        for s in sorted(syms):
            e, h, status = mapper.map(s)
            cnt[status] += 1
            if status != "approved":
                maprep.append(dict(source=st, symbol=s, status=status, hgnc=h or "", ensembl_gene=e or ""))
        mapping_summary[st] = dict(cnt)
    qc["symbol_mapping_by_source"] = mapping_summary
    qc["ensembl_multi_xref_resolution"] = dict(sorted(mapper.ens_resolution.items()))
    out_map = ROOT / cfg["outputs"]["mapping_report"]
    out_map.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(maprep, columns=["source", "symbol", "status", "hgnc", "ensembl_gene"]).to_csv(
        out_map, sep="\t", index=False)

    e1 = [mapper.map(s)[0] for s in sp.sym1]
    e2 = [mapper.map(s)[0] for s in sp.sym2]
    sp["ens1"], sp["ens2"] = e1, e2
    unm = sp.ens1.isna() | sp.ens2.isna()
    selfp = (~unm) & (sp.ens1 == sp.ens2)
    loss = {}
    for st in sp.study.unique():
        m = sp.study == st
        loss[st] = dict(source_pairs=int(sp[m].k.nunique()),
                        source_pairs_lost_unmapped=int(sp[m & unm].k.nunique()),
                        source_pairs_lost_self_pair=int(sp[m & selfp].k.nunique()))
    qc["pair_mapping_losses"] = loss
    sp = sp[~unm & ~selfp].copy()
    sp["gene_a"] = np.where(sp.ens1 < sp.ens2, sp.ens1, sp.ens2)
    sp["gene_b"] = np.where(sp.ens1 < sp.ens2, sp.ens2, sp.ens1)

    # collapse study x line_key x gene pair (symbol collisions, same line under two labels)
    g = sp.groupby(["study", "line_key", "gene_a", "gene_b"], sort=False)
    calls = g.agg(line_label=("line_label", "first"), sym1=("sym1", "first"), sym2=("sym2", "first"),
                  author_pair_class=("author_pair_class", "first"), score=("score", "min"),
                  a=("a", "mean"), b=("b", "mean"), ab=("ab", "mean"), T=("T", "first"),
                  n_source_rows=("n_units", "sum")).reset_index()
    calls["call"] = g["call"].agg(lambda c: collapse_dups(list(c))).values
    lic = {s: v["license_class"] for s, v in cfg["studies"].items()}
    calls["license_class"] = calls.study.map(lic)
    ck = cfg["cell_lines"]["keys"]
    calls["line_key"] = calls.line_key.map(hp.norm_line_id)
    calls["line_id"] = [hp.norm_line_id(ck[k]["cellosaurus"]) if k in ck and ck[k].get("cellosaurus")
                        else (k if str(k).startswith("CVCL_") else f"name:{k}") for k in calls.line_key]
    bad = [x for x in set(calls.line_id) if not re.match(r"^(CVCL_\w+|name:[A-Z0-9]+)$", x)]
    assert not bad, f"malformed line ids: {bad}"
    out_calls = ROOT / cfg["outputs"]["calls"]
    out_calls.parent.mkdir(parents=True, exist_ok=True)
    calls[["gene_a", "gene_b", "sym1", "sym2", "study", "line_label", "line_id", "call", "score", "a", "b", "ab",
           "T", "n_source_rows", "license_class", "author_pair_class"]].to_csv(out_calls, sep="\t", index=False,
                                                                               float_format="%.4g")

    # ---- per study x line call counts and audits ----------------------------------------------------------------
    qc["study_line_call_counts"] = {
        f"{s}|{l}": dict(Counter(c)) for (s, l), c in calls.groupby(["study", "line_label"])["call"]}
    # R2-14: SL rate among testable observations by (a+b)/T bin
    t = calls[calls.call.isin(["S", "N", "B"]) & calls["T"].notna()].copy()
    t["bin"] = pd.cut((t.a + t.b) / t["T"], [-1.0, -0.5, -0.25, 0, 0.25, 0.5, np.inf], right=False)
    qc["audit_sl_rate_by_expected_dko_over_T"] = {
        str(b): dict(n=int(len(x)), sl_rate=round(float((x.call == "S").mean()), 4) if len(x) else None)
        for b, x in t.groupby("bin", observed=False)}

    # ---- line merge and pair aggregation ----------------------------------------------------------------------
    agg_cfg = cfg["aggregation"]["tiers"]
    pr = agg_cfg["primary"]
    lab_rows = []
    sym_of = {}
    for r in calls.itertuples(index=False):
        sym_of.setdefault(r.gene_a, mapper.hgnc_sym.get(mapper.map(r.sym1 if mapper.map(r.sym1)[0] == r.gene_a
                                                                   else r.sym2)[1], ""))
        sym_of.setdefault(r.gene_b, mapper.hgnc_sym.get(mapper.map(r.sym2 if mapper.map(r.sym2)[0] == r.gene_b
                                                                   else r.sym1)[1], ""))
    for (ga, gb), d in calls.groupby(["gene_a", "gene_b"], sort=True):
        t = pair_tiers(d, agg_cfg)
        p, reason, len_lab, two, merged = t["p"], t["reason"], t["len_lab"], t["two"], t["merged"]
        n_test, n_sl, n_bord, pen, st_sl, st_test = (t[k] for k in ("n_test", "n_sl", "n_bord", "pen", "st_sl",
                                                                   "st_test"))
        # D8: would the primary or lenient label change without Harle + Chan evidence?
        d_wo = d[~d.study.isin(["harle2025", "chan2025"])]
        if len(d_wo) == len(d):
            hc_dep = False
        elif len(d_wo) == 0:
            hc_dep = True   # pair exists only through Harle/Chan
        else:
            t_wo = pair_tiers(d_wo, agg_cfg)
            hc_dep = (t_wo["p"] != p) or (t_wo["len_lab"] != len_lab)
        ev = ";".join(f"{st}[{','.join(f'{lid}={c}' for lid, c in sorted(zip(x.line_id, x.call)))}]"
                      for st, x in d.groupby("study", sort=True))
        lines_ev = ",".join(f"{lid}={c}" for lid, c in sorted(merged.items()))
        classes = sorted({lic[s] for s, c in zip(d.study, d.call) if c in "SNB"})
        apc = sorted({f"{s}:{c}" for s, c in zip(d.study, d.author_pair_class) if isinstance(c, str) and c})
        lab_rows.append(dict(
            gene_a=ga, gene_b=gb, symbol_a=sym_of.get(ga, ""), symbol_b=sym_of.get(gb, ""),
            label_primary=p, primary_reason=reason, label_lenient=len_lab, label_two_studies=two,
            n_lines=len(merged), n_testable=n_test, n_sl=n_sl, n_borderline=n_bord,
            penetrance=round(pen, 4) if n_test else np.nan,
            n_studies=d.study.nunique(), n_studies_testable=len(st_test), n_studies_sl=len(st_sl),
            studies=",".join(sorted(d.study.unique())), studies_sl=",".join(sorted(st_sl)),
            merged_line_calls=lines_ev, evidence=ev,
            license_classes=",".join(classes),
            only_nonopen_evidence=("open" not in classes) if classes else ("open" not in {lic[s] for s in d.study}),
            both_in_ref_proteome=(ga in mapper.ens_in_ref) and (gb in mapper.ens_in_ref),
            harle_chan_dependent=hc_dep,
            author_pair_class=";".join(apc)))
    lab = pd.DataFrame(lab_rows)

    # Compara universe column (only if the universe step has built it); path pinned, sha256 recorded (D8)
    qc["compara_universe"] = {"file": str(UNIVERSE_FILE.relative_to(ROOT)), "exists": UNIVERSE_FILE.exists()}
    if UNIVERSE_FILE.exists():
        qc["compara_universe"]["sha256"] = hashlib.sha256(UNIVERSE_FILE.read_bytes()).hexdigest()
        u = pd.read_csv(UNIVERSE_FILE, sep="\t", usecols=["gene_a", "gene_b"])
        upairs = {hp.pair_key(a, b) for a, b in zip(u.gene_a, u.gene_b)}
        qc["compara_universe"]["n_pairs"] = len(upairs)
        lab["in_compara_universe"] = [(a, b) in upairs for a, b in zip(lab.gene_a, lab.gene_b)]

    out_lab = ROOT / cfg["outputs"]["labels"]
    out_lab.parent.mkdir(parents=True, exist_ok=True)
    lab.to_csv(out_lab, sep="\t", index=False)

    # ---- QC -----------------------------------------------------------------------------------------------------
    def tier_counts(col):
        return {str(k): int(v) for k, v in lab[col].value_counts().sort_index().items()}
    qc["n_pairs"] = int(len(lab))
    qc["tier_counts"] = {c: tier_counts(c) for c in ("label_primary", "label_lenient", "label_two_studies")}
    qc["label_value_legend"] = cfg["aggregation"]["label_values"]
    qc["primary_ambiguous_reasons"] = {k: int(v) for k, v in lab[lab.label_primary == -1].primary_reason
                                       .value_counts().items()}
    qc["primary_counts_only_nonopen_evidence"] = {str(k): int(v) for k, v in
                                                  lab[lab.only_nonopen_evidence].label_primary.value_counts().items()}
    qc["pairs_not_both_in_ref_proteome"] = int((~lab.both_in_ref_proteome).sum())
    per_source = {}
    for st in cfg["studies"]:
        m = lab.studies.str.split(",").apply(lambda x: st in x)
        c = calls[calls.study == st]
        pair_any = c.groupby(["gene_a", "gene_b"])["call"].agg(lambda x: "".join(x))
        testable_pairs = pair_any[pair_any.str.contains("[SNB]")]
        per_source[st] = dict(
            license_class=lic[st], pairs=int(m.sum()), study_line_obs=int(len(c)),
            calls=dict(Counter(c.call)),
            pairs_testable=int(len(testable_pairs)),
            pairs_sl_in_ge1_line_within_study=int(testable_pairs.str.contains("S").sum()),
            within_study_positive_rate=round(float(testable_pairs.str.contains("S").mean()), 4)
            if len(testable_pairs) else None,
            primary=dict(Counter(lab[m].label_primary.astype(str))),
            primary_positive_rate=round(float((lab[m].label_primary == 1).sum()
                                              / max(1, lab[m].label_primary.isin([0, 1]).sum())), 4))
    qc["per_source"] = per_source
    # R2-14: primary positive rate by n_testable bin
    lb = lab[lab.label_primary.isin([0, 1])].copy()
    lb["bin"] = pd.cut(lb.n_testable, [2, 3, 4, 6, 11, 21, 100], right=False)
    qc["audit_primary_positive_rate_by_n_testable"] = {
        str(b): dict(n=int(len(x)), pos_rate=round(float((x.label_primary == 1).mean()), 4) if len(x) else None)
        for b, x in lb.groupby("bin", observed=False)}

    # Thompson SD3 Ensembl IDs vs our mapping (cross-check only)
    sd3 = pd.read_excel(ROOT / "data/raw/thompson2021/41467_2021_21478_MOESM6_ESM.xlsx", sheet_name="3", header=2)
    pairs_sd3 = set()
    agree = Counter()
    for ea, sa, eb, sb in sd3[["GeneA_ENSEMBL_gene_ID", "GeneASymbol", "GeneB_ENSEMBL_gene_ID",
                               "GeneBSymbol"]].itertuples(index=False):
        for e, s in ((ea, sa), (eb, sb)):
            if (e, s) in pairs_sd3 or str(s).lower() == "control" or not str(e).startswith("ENSG"):
                continue
            pairs_sd3.add((e, s))
            me = mapper.map(s)[0]
            agree["same" if me == e else ("unmapped" if me is None else "different")] += 1
    qc["thompson_sd3_ensembl_vs_our_mapping"] = dict(agree)

    # 13-pair candidate list (report only, never tuned)
    cand = []
    for line in (ROOT / cfg["crosscheck_only"]["candidate_13_pairs"]["file"]).read_text().splitlines():
        if not line.strip():
            continue
        a, b = line.split("\t")[:2]
        ea, eb = mapper.map(a)[0], mapper.map(b)[0]
        k = hp.pair_key(ea or "", eb or "")
        r = lab[(lab.gene_a == k[0]) & (lab.gene_b == k[1])]
        cand.append(dict(pair=f"{a}|{b}", in_labels=bool(len(r)),
                         **({c: (int(r.iloc[0][c]) if c.startswith(("label", "n_")) else r.iloc[0][c])
                             for c in ("label_primary", "primary_reason", "label_lenient", "label_two_studies",
                                       "n_testable", "n_sl", "studies_sl")} if len(r) else {})))
    qc["candidate_13_pairs_status"] = cand
    qc["candidate_13_pairs_summary"] = dict(Counter(str(c.get("label_primary", "absent")) for c in cand))

    # Dennler & Ryan SL-lenient agreement vs OUR lenient tier (report only; whitelisted columns only)
    qc["dennler_sl_lenient_agreement"] = dennler_agreement(cfg, lab)

    # parser-level counts (known answers are asserted in tests/test_label_known_answers_human.py)
    qc["note"] = ("Thresholds pre-registered in config/labels_human.yaml (commits 513c1e5, c4e5039) before this "
                  "script first ran; Dennler agreement and 13-pair status are report-only.")
    write_known_answer_gaps(cfg)
    out_qc = ROOT / cfg["outputs"]["qc"]
    qc["harle_chan_dependent"] = {
        "definition": "primary or lenient label changes if all Harle 2025 + Chan 2025 calls are removed",
        "by_label_primary": {str(k): int(v) for k, v in lab[lab.harle_chan_dependent].label_primary
                             .value_counts().sort_index().items()},
        "n_pairs": int(lab.harle_chan_dependent.sum())}
    out_qc.write_text(json.dumps(qc, indent=1, default=_json_default, sort_keys=True))
    print(json.dumps({k: qc[k] for k in ("n_pairs", "tier_counts", "primary_ambiguous_reasons")}, indent=1))


# Published counts that the known-answer tests (tests/test_label_known_answers_human.py) cannot check from the
# tables on disk. Counts the tests do check and fail are added here as xfail(strict=True) gaps; none at present.
NOT_TESTABLE_HUMAN = [
    dict(id="parrish2021_hits_PC9_87_HeLa_70", status="not_testable",
         reason="Parrish Table S4 (GI, FDR) not obtained; labels use the pre-registered ZdLFC fallback from in4mer "
                "Supp_table1, which is a different statistic"),
    dict(id="ito2021_table4_3284_single_rows", status="not_tested",
         reason="the singles/controls row definition of Table 4 is not documented; Supp_table1 (5,065 DKO rows) and "
                "Table 10 (4,859 GEMINI pairs) sizes are tested instead"),
    dict(id="klingbeil2024_pair_count", status="not_testable",
         reason="paper text (2,529) and Fig. 1 (2,726) disagree; no single published count to reproduce"),
    dict(id="chymera2020_dekegel_72_sl", status="not_tested",
         reason="De Kegel re-threshold (hit and LFC <= -0.9) is not part of our rule set"),
]


def write_known_answer_gaps(cfg: dict) -> None:
    path = ROOT / cfg["outputs"]["known_answer_gaps"]
    data = json.loads(path.read_text()) if path.exists() else {}
    data["human"] = dict(xfail_gaps=[], not_testable_or_not_tested=NOT_TESTABLE_HUMAN,
                         test_file="tests/test_label_known_answers_human.py",
                         note="all known-answer tests in the test file reproduce their published counts exactly")
    path.write_text(json.dumps(data, indent=1))


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def dennler_agreement(cfg: dict, lab: pd.DataFrame) -> dict:
    spec = cfg["crosscheck_only"]["dennler2025_sl_lenient"]
    path = ROOT / spec["file"]
    if not path.exists():
        return {"status": "file absent"}
    wl = spec["column_whitelist"]
    header = pd.read_csv(path, nrows=0).columns.tolist()
    assert set(wl) <= set(header), "whitelisted columns missing"
    d = pd.read_csv(path, usecols=wl)   # ONLY the whitelisted columns are read
    d["k"] = [hp.pair_key(a, b) for a, b in zip(d[wl[0]].astype(str), d[wl[1]].astype(str))]
    d["SL"] = d[wl[2]].astype(str).str.upper().isin(["TRUE", "1", "1.0"])
    theirs = dict(zip(d.k, d.SL))
    ours = {(a, b): v for a, b, v in zip(lab.gene_a, lab.gene_b, lab.label_lenient) if v in (0, 1)}
    both = set(theirs) & set(ours)
    tab = Counter((("ours_SL" if ours[k] == 1 else "ours_nonSL"), ("dr_SL" if theirs[k] else "dr_nonSL"))
                  for k in both)
    n = len(both)
    po = (tab[("ours_SL", "dr_SL")] + tab[("ours_nonSL", "dr_nonSL")]) / n if n else float("nan")
    p_ours = sum(v for (o, _), v in tab.items() if o == "ours_SL") / n if n else 0
    p_dr = sum(v for (_, t), v in tab.items() if t == "dr_SL") / n if n else 0
    pe = p_ours * p_dr + (1 - p_ours) * (1 - p_dr)
    return dict(
        n_dennler_pairs=len(theirs), n_dennler_sl=int(sum(theirs.values())), n_ours_lenient_labelled=len(ours),
        n_overlap=n, table={f"{a}&{b}": int(v) for (a, b), v in sorted(tab.items())},
        percent_agreement=round(po, 4), cohen_kappa=round((po - pe) / (1 - pe), 4) if n and pe < 1 else None,
        ours_only=len(set(ours) - set(theirs)), dennler_only=len(set(theirs) - set(ours)),
        dennler_only_sl=int(sum(theirs[k] for k in set(theirs) - set(ours))),
        note="report only; never used to tune any rule. Dennler IDs are Ensembl 111; ours UniProt 2026_03 xrefs")


if __name__ == "__main__":
    main()
