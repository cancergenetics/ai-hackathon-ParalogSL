#!/usr/bin/env python
"""Build yeast paralog SL labels (DESIGN v3 §2.3; config/labels_yeast.yaml v1; data.md R2-1/2/3/9/10/11; D4).

Outputs
  data/interim/labels/yeast_calls.tsv                     measurement-level calls (null x null SGA rows; not committed)
  data/interim/labels/yeast_pair_labels_unrestricted.tsv  pair-level labels for every tested pair (not committed)
  results/labels/yeast_label_qc.json                      QC counts (committed)
  data/labels/yeast_labels_all.tsv                        written by scripts/labels/apply_universe.py (called at the
                                                          end of this script) once the paralog universe file exists.
Run:  .venv/bin/python scripts/labels/build_yeast_labels.py
Memory: SGA files are streamed in chunks; per-measurement state is kept as compact numpy arrays (~0.3 GB).
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import yeast_parsers as P  # noqa: E402
import yeast_rules as R  # noqa: E402

ROOT = HERE.parents[1]
CFG_PATH = ROOT / "config" / "labels_yeast.yaml"
TIERS = ["primary", "lenient", "strict", "biogrid", "dennler_comparable"]
CALLS = ["pos", "neg", "amb", "untestable", "untested"]


def log(msg: str) -> None:
    import resource
    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9  # bytes on macOS
    print(f"[{time.strftime('%H:%M:%S')}] {msg} (peak RSS so far {peak_gb:.2f} GB)", flush=True)


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_head() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    except Exception:
        return ""


def git_dirty() -> list[str]:
    """Yeast label code/config files that differ from HEAD at run time (should be empty)."""
    paths = ["config/labels_yeast.yaml", "scripts/labels/build_yeast_labels.py", "scripts/labels/yeast_parsers.py",
             "scripts/labels/yeast_rules.py", "scripts/labels/apply_universe.py"]
    r = subprocess.run(["git", "status", "--porcelain", "--", *paths], cwd=ROOT, capture_output=True, text=True)
    return [ln for ln in r.stdout.splitlines() if ln.strip()]


def vc(x) -> dict:
    s = pd.Series(x).value_counts(dropna=False)
    return {str(k): int(v) for k, v in s.items() if v > 0}


# ============================================================================ SGA measurements
def build_sga(cfg: dict, qc: dict):
    c = cfg["sources"]["costanzo2016"]
    wl = cfg["allele_whitelist"]["primary"]
    files = P.ensure_unpacked(ROOT / c["zip"], ROOT / c["unpack_to"],
                              c["files_labels"] + c["files_known_answer_only"] + [c["strain_table"]])
    strains = P.read_strain_table(files[c["strain_table"]])
    qc["strain_table"] = {"n_strains": int(len(strains)), "suffix_counts": vc(strains["suffix"])}
    genes = sorted(set(strains["gene"]))
    gidx = {g: i for i, g in enumerate(genes)}
    s_gene = dict(zip(strains["strain_id"], strains["gene"]))
    s_smf26 = dict(zip(strains["strain_id"], strains["smf26"]))
    s_smf30 = dict(zip(strains["strain_id"], strains["smf30"]))

    calls_path = ROOT / cfg["qc"]["outputs"]["calls"]
    calls_path.parent.mkdir(parents=True, exist_ok=True)
    if calls_path.exists():
        calls_path.unlink()
    first = True
    acc = {k: [] for k in ["ia", "ib", "aq", "st", "lst", "spos", "cmp", "ep", "minsmf", "E", "posraw"]}
    mq: dict = {"rows_read": {}, "suffix_combos": {}, "null_null_rows": 0, "unmapped_strain_rows": 0,
                "self_pair_rows": 0, "smf_source": {"row": 0, "strain_table": 0, "none": 0}}
    for fname in c["files_labels"]:
        n_rows = 0
        combos: dict = {}
        for ch in P.iter_sga(files[fname], int(c["chunksize"])):
            n_rows += len(ch)
            qs = P.suffix_class(ch["query_strain"])
            as_ = P.suffix_class(ch["array_strain"])
            cc = pd.DataFrame({"q": qs, "a": as_, "t": ch["arraytype_temp"]}).value_counts(dropna=False)
            for k, v in cc.items():
                key = "|".join(str(x) for x in k)
                combos[key] = combos.get(key, 0) + int(v)
            keep = qs.isin(wl["query_suffix"]) & as_.isin(wl["array_suffix"]) & ch["arraytype_temp"].isin(wl["arraytype_temp"])
            ch = ch[keep.fillna(False).to_numpy()].copy()
            if ch.empty:
                continue
            mq["null_null_rows"] += len(ch)
            gq = ch["query_strain"].map(s_gene)
            ga = ch["array_strain"].map(s_gene)
            unm = gq.isna() | ga.isna()
            mq["unmapped_strain_rows"] += int(unm.sum())
            ch, gq, ga = ch[~unm.to_numpy()], gq[~unm], ga[~unm]
            self_ = (gq == ga).to_numpy()
            mq["self_pair_rows"] += int(self_.sum())
            ch, gq, ga = ch[~self_], gq[~self_], ga[~self_]
            # SMF: row values first, strain-table value at the row's temperature as fallback (R2-1)
            t30 = ch["arraytype_temp"].str.endswith("30").to_numpy()
            fq = np.where(t30, ch["query_strain"].map(s_smf30), ch["query_strain"].map(s_smf26)).astype(float)
            fa = np.where(t30, ch["array_strain"].map(s_smf30), ch["array_strain"].map(s_smf26)).astype(float)
            rq, ra = ch["smf_q"].to_numpy(float), ch["smf_a"].to_numpy(float)
            row_ok = ~np.isnan(rq) & ~np.isnan(ra)
            smf_q = np.where(row_ok, rq, np.where(np.isnan(rq), fq, rq))
            smf_a = np.where(row_ok, ra, np.where(np.isnan(ra), fa, ra))
            E = smf_q * smf_a
            src = np.where(row_ok, "row", np.where(np.isnan(E), "none", "strain_table"))
            for k, v in pd.Series(src).value_counts().items():
                mq["smf_source"][k] += int(v)
            eps, p, dmf = ch["eps"].to_numpy(float), ch["p"].to_numpy(float), ch["dmf"].to_numpy(float)
            st, lst, spos, cmp_pos = R.measurement_calls(eps, p, dmf, E)
            with np.errstate(invalid="ignore"):
                posraw = (eps < R.POS_EPS) & (p < R.POS_P) & (dmf <= R.POS_DMF)
            iq = gq.map(gidx).to_numpy(np.int32)
            ia_ = ga.map(gidx).to_numpy(np.int32)
            a_is_q = iq < ia_
            lo, hi = np.minimum(iq, ia_), np.maximum(iq, ia_)
            acc["ia"].append(lo); acc["ib"].append(hi); acc["aq"].append(a_is_q)
            acc["st"].append(st); acc["lst"].append(lst); acc["spos"].append(spos); acc["cmp"].append(cmp_pos)
            acc["minsmf"].append(np.fmin(smf_q, smf_a).astype(np.float32)); acc["E"].append(E.astype(np.float32))
            acc["posraw"].append(posraw)
            acc["ep"].append(~np.isnan(eps) & ~np.isnan(p))
            out = pd.DataFrame({
                "source_file": fname, "query_strain": ch["query_strain"].to_numpy(), "array_strain": ch["array_strain"].to_numpy(),
                "arraytype_temp": ch["arraytype_temp"].to_numpy(),
                "gene_a": np.array(genes, dtype=object)[lo], "gene_b": np.array(genes, dtype=object)[hi],
                "gene_a_is_query": a_is_q.astype(np.int8),
                "eps": eps, "p": p, "smf_query": smf_q, "smf_array": smf_a, "smf_source": src, "E": E, "dmf": dmf,
                "call_primary": pd.Series(st).map(R.STATUS_NAME).to_numpy(),
                "call_lenient": pd.Series(lst).map(R.STATUS_NAME).to_numpy(),
                "strict_pos": spos.astype(np.int8), "eps_lt_m035_p_lt_005": cmp_pos.astype(np.int8),
            })
            out.to_csv(calls_path, sep="\t", index=False, header=first, mode="w" if first else "a")
            first = False
        mq["rows_read"][fname] = n_rows
        mq["suffix_combos"][fname] = combos
        log(f"{fname}: {n_rows} rows read; null x null so far {mq['null_null_rows']}")
    M = {k: np.concatenate(v) for k, v in acc.items()}
    mq["measurements_kept"] = int(len(M["st"]))
    mq["status_primary"] = {R.STATUS_NAME[k]: int((M["st"] == k).sum()) for k in R.STATUS_NAME}
    mq["status_lenient"] = {R.STATUS_NAME[k]: int((M["lst"] == k).sum()) for k in R.STATUS_NAME}
    mq["strict_pos"] = int(M["spos"].sum())
    mq["eps_lt_m035_p_lt_005"] = int(M["cmp"].sum())
    # positive rate by E bin (v2 POS rule without the floor; R2-1 / R2-14 audit)
    edges = [-np.inf] + list(cfg["qc"]["e_bins"]) + [np.inf]
    labels = ["E<=0.5", "0.5<E<=0.6", "0.6<E<=0.7", "0.7<E<=0.85", "E>0.85"]
    Eb = pd.cut(pd.Series(M["E"]), bins=edges, labels=labels, right=True)
    tab = pd.DataFrame({"bin": Eb, "posraw": M["posraw"]}).groupby("bin", observed=False)["posraw"].agg(["size", "sum"])
    mq["pos_rule_rate_by_E_bin"] = {str(k): {"n": int(r["size"]), "n_pos_rule": int(r["sum"]),
                                             "rate": (float(r["sum"]) / r["size"]) if r["size"] else None}
                                    for k, r in tab.iterrows()}
    mq["E_undefined"] = int(np.isnan(M["E"]).sum())
    unt = M["st"] == R.UNTESTABLE
    mq["untestable_reasons"] = {"no_smf": int((unt & np.isnan(M["E"])).sum()),
                                "sick_singles_E_le_0.5": int((unt & ~np.isnan(M["E"])).sum())}
    mq["smf_fallback_note"] = ("strain-table fallback filled no values: every strain whose row SMF is NaN also has NaN "
                               "SMF in the strain table, so those measurements are untestable (no_smf)")
    qc["measurements"] = mq
    return genes, strains, M


def aggregate_sga(genes, M, cfg, qc):
    """Pair-level counts with numpy bincount on factorised pair keys (memory-light). Index = pair key into `genes`."""
    G = len(genes)
    key = M["ia"].astype(np.int64) * G + M["ib"].astype(np.int64)
    uk, inv = np.unique(key, return_inverse=True)
    del key
    n = len(uk)

    def cnt(mask):
        return np.bincount(inv[mask], minlength=n).astype(np.int32)

    st = M["st"]
    pr = pd.DataFrame({
        "ia": (uk // G).astype(np.int32), "ib": (uk % G).astype(np.int32),
        "n_meas": np.bincount(inv, minlength=n).astype(np.int32),
        "n_pos": cnt(st == R.POS), "n_neg": cnt(st == R.NEG), "n_mid": cnt(st == R.MID),
        "n_untestable": cnt(st == R.UNTESTABLE), "n_missing": cnt(st == R.MISSING),
        "n_pos_len": cnt(M["lst"] == R.POS), "n_spos": cnt(M["spos"]), "n_cmp_pos": cnt(M["cmp"]),
        "n_ep": cnt(M["ep"]), "n_a_query": cnt(M["aq"]), "n_nosmf": cnt(np.isnan(M["E"])),
    })
    mins = np.full(n, np.inf, dtype=np.float32)
    np.fmin.at(mins, inv, M["minsmf"])
    pr["min_smf"] = np.where(np.isinf(mins), np.nan, mins)
    pr["n_b_query"] = pr["n_meas"] - pr["n_a_query"]
    pr["n_T"] = pr["n_pos"] + pr["n_neg"] + pr["n_mid"]
    pr["n_mid_len"] = pr["n_T"] - pr["n_pos_len"] - pr["n_neg"]
    pr["sga_primary"], pr["discordant"] = R.aggregate(pr["n_pos"], pr["n_neg"], pr["n_mid"])
    pr["sga_lenient"], pr["discordant_lenient"] = R.aggregate(pr["n_pos_len"], pr["n_neg"], pr["n_mid_len"])
    pr["sga_strict"] = R.strict_sga(pr["sga_primary"], pr["n_spos"], pr["n_neg"])
    # report-only comparable tier: any eps<-0.35 & p<0.05 among rows with eps/p present
    pr["sga_cmp"] = np.where(pr["n_cmp_pos"] > 0, "pos", np.where(pr["n_ep"] > 0, "neg", "untested"))
    # orientation-level calls (flag only, D4): orientation 1 = gene_a is the query
    oinv = inv.astype(np.int64) * 2 + M["aq"].astype(np.int64)
    oc = {}
    for nm, code in [("p", R.POS), ("n", R.NEG), ("m", R.MID)]:
        oc[nm] = np.bincount(oinv[st == code], minlength=2 * n).reshape(n, 2)
    calls = [R.aggregate(oc["p"][:, k], oc["n"][:, k], oc["m"][:, k])[0] for k in (0, 1)]
    pr["orientation_conflict"] = (((calls[0] == "pos") & (calls[1] == "neg"))
                                  | ((calls[0] == "neg") & (calls[1] == "pos"))).astype(np.int8)
    del oinv, oc, calls, inv
    for col in ["sga_primary", "sga_lenient", "sga_strict", "sga_cmp"]:
        pr[col] = pd.Categorical(pr[col], categories=CALLS)
    for col in ["discordant", "discordant_lenient"]:
        pr[col] = pr[col].astype(np.int8)
    garr = np.array(genes, dtype=object)
    qc["sga_pairs"] = {
        "n_pairs_measured": int(len(pr)),
        "sga_primary": vc(pr["sga_primary"]), "sga_lenient": vc(pr["sga_lenient"]), "sga_strict": vc(pr["sga_strict"]),
        "untestable_sick_singles": int(((pr["n_T"] == 0) & (pr["n_untestable"] > pr["n_nosmf"])).sum()),
        "untestable_no_smf_only": int(((pr["n_T"] == 0) & (pr["n_untestable"] > 0) & (pr["n_untestable"] == pr["n_nosmf"])).sum()),
        "untestable_note": ("pairs with n_T == 0 split by reason: sick_singles = at least one measurement with E <= 0.5; "
                            "no_smf_only = every measurement lacks SMF (E undefined). Both are label 'untestable' (rule unchanged)"),
        "untestable_missing_only": int(((pr["n_T"] == 0) & (pr["n_untestable"] == 0)).sum()),
        "discordant": int(pr["discordant"].sum()),
        "orientation_conflict_by_sga_primary": vc(pr.loc[pr["orientation_conflict"] == 1, "sga_primary"]),
        "both_orientations_measured": int(((pr["n_a_query"] > 0) & (pr["n_b_query"] > 0)).sum()),
        "max_n_T": int(pr["n_T"].max()),
    }
    b = pd.cut(pr["n_T"], bins=[-1, 0, 1, 2, 4, np.inf], labels=["0", "1", "2", "3-4", ">=5"])
    qc["sga_pairs"]["sga_primary_by_n_T"] = {str(k): vc(v) for k, v in pr.groupby(b, observed=True)["sga_primary"]}
    pr["gene_a"] = garr[pr["ia"].to_numpy()]
    pr["gene_b"] = garr[pr["ib"].to_numpy()]
    return pr


# ============================================================================ Dean, BioGRID, coordinates
def build_dean(cfg, qc):
    d = P.parse_dean_pdf(ROOT / cfg["sources"]["dean2008"]["pdf"])
    pa = cfg["dean2008"]["parser_assert"]
    measured = d["w12"].notna()
    got = {"pairs_listed": int(len(d)), "measured": int(measured.sum()), "lethal": int((d["w12"] == 0).sum())}
    if got != {k: pa[k] for k in got}:
        raise AssertionError(f"Dean parser assert failed: got {got}, expected {pa}")
    d["gene_a"] = np.where(d["orf1"] < d["orf2"], d["orf1"], d["orf2"])
    d["gene_b"] = np.where(d["orf1"] < d["orf2"], d["orf2"], d["orf1"])
    if d.duplicated(["gene_a", "gene_b"]).any():
        raise ValueError("duplicate Dean pairs")
    d["dean_call"] = R.dean_call(d["w1"], d["w2"], d["w12"])
    d["dean_evidence"] = ("DEAN[W=" + d["w1"].map("{:.2f}".format) + "," + d["w2"].map("{:.2f}".format) + ","
                          + d["w12"].map("{:.2f}".format) + ";" + d["dean_call"] + "]")
    qc["dean2008"] = {"parse": got, "calls": vc(d["dean_call"])}
    return d[["gene_a", "gene_b", "dean_call", "dean_evidence", "w1", "w2", "w12"]]


def build_biogrid(cfg, strains, qc):
    b = cfg["sources"]["biogrid"]
    org = str(b["organism_id"])
    null_genes = set(strains.loc[strains["suffix"].isin(["sn", "dma"]), "gene"])
    n_rows = n_gen = 0
    recs = []
    for ch in P.iter_biogrid(ROOT / b["zip"], b["member"]):
        n_rows += len(ch)
        ch = ch[ch["system_type"] == "genetic"]
        n_gen += len(ch)
        lt = ch["throughput"].str.contains("Low Throughput", regex=False)
        sel = ch[lt & ch["system"].isin(["Synthetic Lethality", "Synthetic Growth Defect"])
                 & (ch["org_a"] == org) & (ch["org_b"] == org)]
        recs.append(sel[["sys_a", "sys_b", "system", "publication"]])
    r = pd.concat(recs, ignore_index=True)
    st = {"rows_read": n_rows, "genetic_rows": n_gen, "lt_sl_or_sgd_records_559292": int(len(r))}
    bad = (r["sys_a"].isin(["", "-"])) | (r["sys_b"].isin(["", "-"]))
    st["dropped_no_systematic_name"] = int(bad.sum())
    r = r[~bad]
    selfp = r["sys_a"] == r["sys_b"]
    st["dropped_self"] = int(selfp.sum())
    r = r[~selfp].copy()
    ess = ~(r["sys_a"].isin(null_genes) & r["sys_b"].isin(null_genes))
    r["essential_member"] = ess
    st["sl_records_essential_member"] = int((ess & (r["system"] == "Synthetic Lethality")).sum())
    r["gene_a"] = np.where(r["sys_a"] < r["sys_b"], r["sys_a"], r["sys_b"])
    r["gene_b"] = np.where(r["sys_a"] < r["sys_b"], r["sys_b"], r["sys_a"])
    sl_q = r[(r["system"] == "Synthetic Lethality") & ~r["essential_member"]]
    sgd = r[(r["system"] == "Synthetic Growth Defect") & ~r["essential_member"]]
    a = sl_q.groupby(["gene_a", "gene_b"]).size().rename("n_biogrid_lt_sl")
    s = sgd.groupby(["gene_a", "gene_b"]).size().rename("n_biogrid_lt_sgd")
    e = r[(r["system"] == "Synthetic Lethality") & r["essential_member"]].groupby(["gene_a", "gene_b"]).size().rename("n_biogrid_lt_sl_essential")
    out = pd.concat([a, s, e], axis=1).fillna(0).astype(int).reset_index()
    st["qualifying_lt_sl_records"] = int(len(sl_q))
    st["qualifying_lt_sl_pairs"] = int((out["n_biogrid_lt_sl"] > 0).sum())
    st["lt_sgd_pairs_evidence_only"] = int((out["n_biogrid_lt_sgd"] > 0).sum())
    qc["biogrid"] = st
    return out


def kuzmin_status(cfg, qc):
    z = ROOT / cfg["sources"]["kuzmin2020"]["zip"]
    qc["kuzmin2020"] = {"status": "present_parser_not_implemented" if z.exists() else "not_obtained_dropped (R2-13 fallback)"}


# ============================================================================ assemble
def assemble(cfg, genes, pr, dean, bg, qc):
    """Align SGA pairs, Dean and BioGRID on integer pair keys over one sorted gene list (no string outer merges)."""
    all_genes = np.array(sorted(set(genes) | set(dean["gene_a"]) | set(dean["gene_b"]) | set(bg["gene_a"]) | set(bg["gene_b"])),
                         dtype=object)
    G2 = len(all_genes)
    gidx = pd.Index(all_genes)
    remap = gidx.get_indexer(np.array(genes, dtype=object)).astype(np.int64)
    k_sga = remap[pr["ia"].to_numpy()] * G2 + remap[pr["ib"].to_numpy()]
    k_dean = gidx.get_indexer(dean["gene_a"]).astype(np.int64) * G2 + gidx.get_indexer(dean["gene_b"])
    k_bg = gidx.get_indexer(bg["gene_a"]).astype(np.int64) * G2 + gidx.get_indexer(bg["gene_b"])
    keys = np.union1d(np.union1d(k_sga, k_dean), k_bg)
    N = len(keys)
    df = pd.DataFrame({"gene_a": all_genes[keys // G2], "gene_b": all_genes[keys % G2]})
    if not (df["gene_a"] < df["gene_b"]).all():
        raise AssertionError("pair keys not canonical")
    ps, pdn, pbg = np.searchsorted(keys, k_sga), np.searchsorted(keys, k_dean), np.searchsorted(keys, k_bg)
    in_sga = np.zeros(N, dtype=bool); in_sga[ps] = True
    in_sga = pd.Series(in_sga)
    for col in ["n_meas", "n_pos", "n_neg", "n_mid", "n_untestable", "n_missing", "n_pos_len", "n_mid_len", "n_spos",
                "n_cmp_pos", "n_a_query", "n_b_query", "n_T", "discordant", "discordant_lenient", "orientation_conflict"]:
        a = np.zeros(N, dtype=pr[col].dtype); a[ps] = pr[col].to_numpy(); df[col] = a
    a = np.full(N, np.nan, dtype=np.float32); a[ps] = pr["min_smf"].to_numpy(); df["min_smf"] = a
    for col in ["sga_primary", "sga_lenient", "sga_strict", "sga_cmp"]:
        codes = np.full(N, CALLS.index("untested"), dtype=np.int8); codes[ps] = pr[col].cat.codes.to_numpy()
        df[col] = pd.Categorical.from_codes(codes, categories=CALLS)
    del pr
    dc = np.full(N, "", dtype=object); dc[pdn] = dean["dean_call"].to_numpy(); df["dean_call"] = dc
    de = np.full(N, "", dtype=object); de[pdn] = dean["dean_evidence"].to_numpy(); df["dean_evidence"] = de
    for col in ["n_biogrid_lt_sl", "n_biogrid_lt_sgd", "n_biogrid_lt_sl_essential"]:
        a = np.zeros(N, dtype=np.int32); a[pbg] = bg[col].to_numpy(); df[col] = a
    # linked pairs (R2-11) from Ensembl 116 GTF coordinates
    coords = P.read_gtf_genes(ROOT / cfg["sources"]["gene_coordinates"]).set_index("gene")
    ca = df["gene_a"].map(coords["chrom"]); cb = df["gene_b"].map(coords["chrom"])
    linked, gap = R.linked_pairs(ca, df["gene_a"].map(coords["start"]), df["gene_a"].map(coords["end"]),
                                 cb, df["gene_b"].map(coords["start"]), df["gene_b"].map(coords["end"]),
                                 window=int(cfg["linked"]["window_bp"]))
    df["linked"] = linked
    df["same_chrom_gap_bp"] = gap
    lk = df["linked"] == 1
    sga_t = {}
    for t, col in [("primary", "sga_primary"), ("lenient", "sga_lenient"), ("strict", "sga_strict")]:
        v = df[col].astype(object).to_numpy()
        sga_t[t] = np.where(lk.to_numpy() & (v != "untested"), "untestable", v)
    df["primary"], df["dean_override"], df["dean_conflict"] = R.merge_dean(sga_t["primary"], df["dean_call"])
    df["lenient"], _, _ = R.merge_dean(sga_t["lenient"], df["dean_call"])
    df["strict"] = R.strict_final(sga_t["strict"], df["dean_call"], df["primary"])
    df["biogrid"], df["biogrid_only"], df["biogrid_conflict"] = R.biogrid_tier(df["primary"], df["n_biogrid_lt_sl"])
    df["dennler_comparable"] = df["sga_cmp"].astype(object)
    for col in ["dean_override", "dean_conflict", "biogrid_only", "biogrid_conflict"]:
        df[col] = df[col].astype(np.int8)
    # source as a categorical built from codes (memory); evidence strings are built per chunk at write time
    code = (in_sga.to_numpy().astype(int) + 2 * (df["dean_call"] != "").to_numpy().astype(int)
            + 4 * ((df["n_biogrid_lt_sl"] + df["n_biogrid_lt_sgd"] + df["n_biogrid_lt_sl_essential"]) > 0).to_numpy())
    names = [";".join(n for bit, n in zip((1, 2, 4), ("costanzo2016", "dean2008", "biogrid_lt")) if c & bit) for c in range(8)]
    df["source"] = pd.Categorical.from_codes(code, categories=names)
    orfs = P.proteome_orfs(ROOT / cfg["sources"]["uniprot_idmapping"])
    df["in_proteome"] = (df["gene_a"].isin(orfs) & df["gene_b"].isin(orfs)).astype(np.int8)
    df["license_note"] = pd.Categorical(np.where(in_sga.to_numpy(), "costanzo2016_no_explicit_license_derived_calls_only", ""))
    for col in TIERS + ["sga_primary", "sga_lenient", "sga_strict", "dean_call"]:
        df[col] = df[col].astype("category")
    # rows are already in (gene_a, gene_b) order: keys index a sorted gene list
    cols = ["gene_a", "gene_b", "primary", "lenient", "strict", "biogrid", "dennler_comparable",
            "sga_primary", "sga_lenient", "sga_strict", "dean_call",
            "n_meas", "n_T", "n_pos", "n_neg", "n_mid", "n_pos_len", "n_spos", "n_untestable", "n_missing", "n_cmp_pos",
            "n_a_query", "n_b_query", "min_smf", "n_biogrid_lt_sl", "n_biogrid_lt_sgd", "n_biogrid_lt_sl_essential",
            "orientation_conflict", "discordant", "linked", "same_chrom_gap_bp", "dean_override", "dean_conflict",
            "biogrid_only", "biogrid_conflict", "in_proteome", "source", "license_note", "dean_evidence"]
    return df[cols]


def evidence_strings(df: pd.DataFrame) -> pd.Series:
    """Compact per-pair evidence codes: SGA[...]|DEAN[...]|BG[...]."""
    in_sga = df["n_meas"] > 0
    sga_ev = ("SGA[m=" + df["n_meas"].astype(str) + ";T=" + df["n_T"].astype(str) + ";P=" + df["n_pos"].astype(str)
              + ";N=" + df["n_neg"].astype(str) + ";M=" + df["n_mid"].astype(str) + ";PL=" + df["n_pos_len"].astype(str)
              + ";SP=" + df["n_spos"].astype(str) + ";U=" + df["n_untestable"].astype(str)
              + ";o=" + df["n_a_query"].astype(str) + "/" + df["n_b_query"].astype(str) + "]")
    ev = pd.Series(np.where(in_sga, sga_ev, ""), index=df.index)
    ev = ev.str.cat(df["dean_evidence"].fillna("").astype(str), sep="|")
    bg = (df["n_biogrid_lt_sl"] + df["n_biogrid_lt_sgd"] + df["n_biogrid_lt_sl_essential"]) > 0
    bg_ev = ("BG[LT_SL=" + df["n_biogrid_lt_sl"].astype(str) + ";LT_SGD=" + df["n_biogrid_lt_sgd"].astype(str)
             + ";LT_SL_ess=" + df["n_biogrid_lt_sl_essential"].astype(str) + "]")
    ev = ev.str.cat(np.where(bg, bg_ev, ""), sep="|")
    return ev.str.replace(r"\|+", "|", regex=True).str.strip("|")


def write_pairs(df: pd.DataFrame, path: Path, chunk: int = 500_000) -> None:
    """Write the pair table in chunks, adding the evidence column per chunk (keeps peak memory low)."""
    cols = [c for c in df.columns if c not in ("dean_evidence", "license_note")] + ["evidence", "license_note"]
    for i in range(0, len(df), chunk):
        part = df.iloc[i:i + chunk].copy()
        part["evidence"] = evidence_strings(part)
        part[cols].to_csv(path, sep="\t", index=False, header=(i == 0), mode="w" if i == 0 else "a")


# ============================================================================ QC helpers (also used by apply_universe)
def tier_qc(df: pd.DataFrame) -> dict:
    out: dict = {"n_pairs": int(len(df))}
    for t in TIERS:
        out[t] = vc(df[t])
    out["primary_by_source"] = {s: vc(v) for s, v in df.groupby("source", observed=True)["primary"]}
    flags = ["orientation_conflict", "discordant", "dean_override", "dean_conflict", "biogrid_only", "biogrid_conflict"]
    out["flags"] = {f: int(df[f].sum()) for f in flags}
    out["orientation_conflict_by_primary"] = vc(df.loc[df["orientation_conflict"] == 1, "primary"])
    out["linked"] = {"n_linked": int((df["linked"] == 1).sum()), "n_unknown": int(df["linked"].isna().sum()),
                     "primary_of_linked": vc(df.loc[df["linked"] == 1, "primary"]),
                     "sga_primary_of_linked_before_rule": vc(df.loc[df["linked"] == 1, "sga_primary"])}
    pos = {t: (df[t] == "pos").to_numpy() for t in ["primary", "lenient", "strict", "biogrid"]}
    neg = {t: (df[t] == "neg").to_numpy() for t in ["primary", "lenient", "strict"]}
    out["nesting"] = {"strict_pos_subset_primary_pos": bool(not (pos["strict"] & ~pos["primary"]).any()),
                      "primary_pos_subset_lenient_pos": bool(not (pos["primary"] & ~pos["lenient"]).any()),
                      "primary_pos_subset_biogrid_pos": bool(not (pos["primary"] & ~pos["biogrid"]).any()),
                      "lenient_neg_subset_primary_neg": bool(not (neg["lenient"] & ~neg["primary"]).any()),
                      "strict_neg_equals_primary_neg": bool((neg["strict"] == neg["primary"]).all())}
    labelled = df["primary"].isin(["pos", "neg"])
    out["primary_prevalence"] = (float((df["primary"] == "pos").sum()) / labelled.sum()) if labelled.sum() else None
    smf_b = pd.cut(df["min_smf"], bins=[-np.inf, 0.7, 0.85, 0.95, np.inf], labels=["<=0.7", "0.7-0.85", "0.85-0.95", ">0.95"])
    out["primary_by_min_smf"] = {str(k): vc(v) for k, v in df.groupby(smf_b, observed=True)["primary"]}
    return out


def dennler_qc(cfg, df) -> dict:
    """Report-only agreement with Dennler & Ryan yeast SL files. Reads ONLY the whitelisted columns."""
    dq = cfg["qc"]["dennler"]
    wl = list(dq["column_whitelist"])
    res: dict = {"column_whitelist": wl, "note": "report only; never used to tune any rule"}
    ours = df[["gene_a", "gene_b", *dq["compare_to"]]]
    for name in ["file", "file_deadly"]:
        path = ROOT / dq[name]
        if not path.exists():
            res[name] = {"status": "missing"}
            continue
        hdr = pd.read_csv(path, nrows=0).columns.tolist()
        if not set(wl) <= set(hdr):
            res[name] = {"status": f"whitelisted columns absent: {sorted(set(wl) - set(hdr))}"}
            continue
        d = pd.read_csv(path, usecols=wl, dtype=str)
        a, b = d[wl[0]].str.strip(), d[wl[1]].str.strip()
        d = pd.DataFrame({"gene_a": np.where(a < b, a, b), "gene_b": np.where(a < b, b, a), "SL": d[wl[2]].str.strip()})
        d = d.drop_duplicates(["gene_a", "gene_b"])
        r = {"n_pairs": int(len(d)), "SL_values": vc(d["SL"])}
        j = d.merge(ours, on=["gene_a", "gene_b"], how="left")
        r["n_in_our_table"] = int(j[dq["compare_to"][0]].notna().sum())
        for t in dq["compare_to"]:
            ct = pd.crosstab(j["SL"], j[t].astype(object).fillna("absent"))
            r[f"crosstab_SL_vs_{t}"] = {str(i): {str(c): int(ct.loc[i, c]) for c in ct.columns} for i in ct.index}
        res[Path(dq[name]).name] = r
    return res


def restricted_report(cfg, df) -> dict:
    """Tier QC plus the report-only external agreement, for a (universe-restricted) label table."""
    return {"tiers": tier_qc(df), "external_agreement": dennler_qc(cfg, df)}


def coverage_by_distance(cfg, strains, pr, qc):
    """QC only: fraction of possible (sn query gene, dma array gene) same-chromosome combinations that were measured,
    by inter-gene gap (SGA linkage depletion; never used to change the window)."""
    coords = P.read_gtf_genes(ROOT / cfg["sources"]["gene_coordinates"]).set_index("gene")
    q = sorted(set(strains.loc[strains["suffix"] == "sn", "gene"]) & set(coords.index))
    a = sorted(set(strains.loc[strains["suffix"] == "dma", "gene"]) & set(coords.index))
    bins = [0, 10_000, 25_000, 50_000, 100_000, 200_000, 500_000]
    lab = ["0-10k", "10-25k", "25-50k", "50-100k", "100-200k", "200-500k"]
    gi = {g: i for i, g in enumerate(sorted(set(pr["gene_a"]) | set(pr["gene_b"]) | set(q) | set(a)))}
    G = len(gi)
    measured = np.sort(pr["gene_a"].map(gi).to_numpy(np.int64) * G + pr["gene_b"].map(gi).to_numpy(np.int64))
    cand_key, cand_bin = [], []
    cq = coords.loc[q]; cA = coords.loc[a]
    for chrom, gq in cq.groupby("chrom"):
        ga = cA[cA["chrom"] == chrom]
        ia = np.array([gi[x] for x in ga.index], dtype=np.int64)
        for gname, rq in gq.iterrows():
            gap = np.maximum(0, np.maximum(rq["start"], ga["start"].to_numpy()) - np.minimum(rq["end"], ga["end"].to_numpy()))
            sel = (ga.index != gname) & (gap < bins[-1])
            iq = gi[gname]
            lo, hi = np.minimum(iq, ia[sel]), np.maximum(iq, ia[sel])
            cand_key.append(lo * G + hi)
            cand_bin.append(np.searchsorted(bins, gap[sel], side="right") - 1)
    ck, cb = np.concatenate(cand_key), np.concatenate(cand_bin)
    pos_ = np.searchsorted(measured, ck)
    hit = (pos_ < len(measured)) & (measured[np.minimum(pos_, len(measured) - 1)] == ck)
    poss = {lab[i]: int((cb == i).sum()) for i in range(len(lab))}
    meas = {lab[i]: int(((cb == i) & hit).sum()) for i in range(len(lab))}
    qc["linked_coverage_qc"] = {k: {"ordered_sn_query_dma_array_combos": poss[k], "of_which_pair_measured_any_orientation": meas[k],
                                    "frac": (meas[k] / poss[k]) if poss[k] else None} for k in lab}


# ============================================================================ main
def main() -> None:
    t0 = time.time()
    cfg = yaml.safe_load(CFG_PATH.read_text())
    qc: dict = {"config": str(CFG_PATH.relative_to(ROOT)), "config_version": cfg["version"],
                "config_sha256": sha256_file(CFG_PATH), "git_head_at_run": git_head(),
                "uncommitted_yeast_code_or_config_at_run": git_dirty(),
                "inputs_sha256": {}}
    for k in ["costanzo2016", "dean2008", "biogrid"]:
        s = cfg["sources"][k]
        p = ROOT / (s.get("zip") or s.get("pdf"))
        qc["inputs_sha256"][k] = sha256_file(p)
    qc["inputs_sha256"]["gene_coordinates"] = sha256_file(ROOT / cfg["sources"]["gene_coordinates"])
    genes, strains, M = build_sga(cfg, qc)
    pr = aggregate_sga(genes, M, cfg, qc)
    M.clear()
    del M
    log(f"SGA pairs: {len(pr)}")
    dean = build_dean(cfg, qc)
    # Dean-independent QC (R2-9): SGA-only primary call vs Dean call, before merging
    j = dean.merge(pr[["gene_a", "gene_b", "sga_primary"]], on=["gene_a", "gene_b"], how="left")
    j["sga_primary"] = j["sga_primary"].fillna("untested")
    ct = pd.crosstab(j["dean_call"], j["sga_primary"])
    qc["dean2008"]["sga_only_vs_dean"] = {str(i): {str(c): int(ct.loc[i, c]) for c in ct.columns} for i in ct.index}
    bg = build_biogrid(cfg, strains, qc)
    kuzmin_status(cfg, qc)
    coverage_by_distance(cfg, strains, pr, qc)
    df = assemble(cfg, genes, pr, dean, bg, qc)
    del pr
    qc["mapping"] = {"pairs_both_genes_in_uniprot_proteome": int(df["in_proteome"].sum()), "pairs_total": int(len(df))}
    qc["linked_rule"] = {"window_bp": cfg["linked"]["window_bp"], "coordinates": cfg["sources"]["gene_coordinates"]}
    qc["unrestricted"] = tier_qc(df)
    qc["dennler_agreement_unrestricted"] = dennler_qc(cfg, df)
    out = ROOT / cfg["qc"]["outputs"]["pairs_unrestricted"]
    out.parent.mkdir(parents=True, exist_ok=True)
    write_pairs(df, out)
    log(f"wrote {out} ({len(df)} pairs)")
    qc_path = ROOT / cfg["qc"]["outputs"]["qc_json"]
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    qc["runtime_s"] = round(time.time() - t0, 1)
    qc_path.write_text(json.dumps(qc, indent=1, default=str))
    log(f"wrote {qc_path}")
    # universe restriction (separate, re-runnable step)
    subprocess.run([sys.executable, str(HERE / "apply_universe.py")], check=True)


if __name__ == "__main__":
    main()
