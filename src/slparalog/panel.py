"""Reference-panel search and divergence calibration for Tier B (DESIGN v3 §2.2, §4 Tier B).

Normative spec: reviews/round2/algorithms_comparative.md changes 2-7, 9-11 (D4). Settings live in
``config/panel.yaml`` (section ``tier_b``); the pinned MMseqs2 scoring settings are read from
``config/caller.yaml`` (``mmseqs.scoring``) so that self- and panel-search bitscores are comparable
(change 6, risk 5). Both enter the schema hash (``settings_hash``).

Steps
-----
1. ``build_dbs``: one MMseqs2 database per panel species (``<db_dir>/<code>/db``, headers = the
   0-based within-species index, so key == header == index), built once and stamped. Choice of
   per-species DBs over one DB with taxonomy filtering: results/bench/panel_db_layout.json.
2. ``search_panel``: for every species, ``prefilter`` (-s, --max-seqs) + ``align -a
   --alignment-mode 3`` + ``convertalis`` (with ``cigar``) of the query DB against the species DB,
   all hits with bits >= keep_min_bits kept. Saturated queries (>= saturation_fraction * cap hits at
   >= presence bits) are re-searched with a larger cap. Results are cached per (query sequences,
   settings) in ``<cache_dir>/<key>/<code>.parquet``.
3. ``reduce_species``: per query and species, the best hit (presence: bits >= 50), the second-best
   bits, the pseudo-RBH flag (change 5.1: x is q's unique top hit, q is x's unique top query, both
   coverages >= 0.5), the RBH identity (nident/alnlen), and the number of hits >= 0.5 x self-bits.
4. ``calibrate``: self-exclusion (median RBH identity >= 0.98, change 9), core set G* (RBH in >= 80%
   of retained species, fallback in 10% steps while |G*| < 200), d_s = Kimura distance of the median
   RBH identity over G* (change 5.3). D_s is computed once on the natural panel.
5. ``apply_mask``: ``--panel-dmin c`` removes species with d_s < c (change 7); bins (frozen edges)
   are recomputed on the remaining species. Everything downstream (tier_b) is recomputed per mask.
6. ``qc_flags``: proteome-level QC warnings (|G*| below ``core.min_size`` after the fallback; nearest
   retained species with d_s > ``QC_NEAREST_D_MAX``). They are output metadata only (Tier B meta
   ``qc``, pipeline meta ``qc_flags``), never a feature column: a proteome-level flag is constant
   within a proteome.

Bin edges (``tier_b.bins.edges`` in ``config/panel.yaml``): the edges are HUMAN-DERIVED (taken from
the reviewer's human core-set d values and frozen before any other species' d_s was computed). They
are used ONLY to form the 1/n_bin weights of the discordance count in ``tb_depth`` conflict
resolution (``apply_mask`` -> ``Mask.weight``); no per-bin quantity is a registered feature. On a
new proteome the same edges are applied to that proteome's own d_s scale.

Ties: a best hit (or top query) that ties on bitscore with the runner-up is not unique, so it never
forms a pseudo-RBH; best hits for presence/classification take the lowest target index among ties.
Query indices are the canonical (sequence-sorted) indices of the pipeline, so nothing here depends
on FASTA order.
"""

from __future__ import annotations

import hashlib
import json
import logging
import math
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .config import CONFIG_DIR, PROJECT_ROOT, load_yaml

log = logging.getLogger(__name__)

DEFAULT_PANEL_CONFIG = CONFIG_DIR / "panel.yaml"
PANEL_COLUMNS = ["query", "target", "bits", "nident", "alnlen", "qstart", "qend", "tstart",
                 "tend", "qlen", "tlen", "cigar"]
_DTYPES = {"query": np.int32, "target": np.int32, "bits": np.float32, "nident": np.int32,
           "alnlen": np.int32, "qstart": np.int32, "qend": np.int32, "tstart": np.int32,
           "tend": np.int32, "qlen": np.int32, "tlen": np.int32, "cigar": str}
DB_STAMP_VERSION = 1


class PanelError(RuntimeError):
    pass


# ----------------------------------------------------------------------------------------------
# configuration
# ----------------------------------------------------------------------------------------------

def _abs(p: str | Path) -> Path:
    p = Path(p)
    return p if p.is_absolute() else PROJECT_ROOT / p


def load_panel_config(path: str | Path | None = None) -> dict:
    """Full panel.yaml; the Tier B section is ``cfg['tier_b']`` and file paths ``cfg['panel']``."""
    return load_yaml(path or DEFAULT_PANEL_CONFIG)


def species_table(pcfg: dict) -> pd.DataFrame:
    df = pd.read_csv(_abs(pcfg["panel"]["species_table"]), sep="\t", dtype={"code": str})
    if df["code"].duplicated().any():
        raise PanelError("duplicate species codes in the species table")
    return df


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def prefilter_scoring_args(ccfg: dict) -> list[str]:
    """Prefilter options from the pinned scoring block of config/caller.yaml (same keys and values
    as the self-search prefilter)."""
    s = ccfg["mmseqs"]["scoring"]
    return ["--seed-sub-mat", s["seed_sub_mat"], "--mask", str(s["mask"]),
            "--mask-prob", str(s["mask_prob"]), "--mask-lower-case", str(s["mask_lower_case"]),
            "--sub-mat", s["sub_mat"], "--comp-bias-corr", str(s["comp_bias_corr"]),
            "--comp-bias-corr-scale", str(s["comp_bias_corr_scale"]),
            "--max-seq-len", str(s["max_seq_len"])]


def align_scoring_args(ccfg: dict) -> list[str]:
    """Alignment scoring options; identical to ``search.scoring_args`` (asserted by a test)."""
    s = ccfg["mmseqs"]["scoring"]
    return ["--sub-mat", s["sub_mat"], "--gap-open", s["gap_open"], "--gap-extend", s["gap_extend"],
            "--comp-bias-corr", str(s["comp_bias_corr"]),
            "--comp-bias-corr-scale", str(s["comp_bias_corr_scale"]),
            "--max-seq-len", str(s["max_seq_len"])]


def settings_payload(pcfg: dict, ccfg: dict, panel_fasta_sha256: str | None = None) -> dict:
    """Everything that changes panel hits or Tier B features (enters the schema hash): the Tier B
    settings (without paths), the pinned scoring settings shared with the self-search, the species table
    and the panel FASTA checksum (from the DB stamp)."""
    return {
        # file-system locations (db_dir, cache_dir) are provenance, not content: the hash must not depend on them
        "tier_b": {k: v for k, v in pcfg["tier_b"].items() if k not in ("db_dir", "cache_dir")},
        "scoring": ccfg["mmseqs"]["scoring"],
        "alignment_mode": ccfg["mmseqs"]["alignment_mode"],
        "species_table_sha256": _sha256_file(_abs(pcfg["panel"]["species_table"])),
        "panel_fasta_sha256": panel_fasta_sha256,
    }


def settings_hash(pcfg: dict, ccfg: dict, panel_fasta_sha256: str | None = None) -> str:
    payload = json.dumps(settings_payload(pcfg, ccfg, panel_fasta_sha256), sort_keys=True,
                         default=str).encode()
    return hashlib.sha256(payload).hexdigest()


def search_settings_hash(pcfg: dict, ccfg: dict, tool_sha: str) -> str:
    """Hash of what determines the raw hit tables (cache key part)."""
    t = pcfg["tier_b"]
    payload = {"search": t["search"], "scoring": ccfg["mmseqs"]["scoring"],
               "alignment_mode": ccfg["mmseqs"]["alignment_mode"], "tool": tool_sha,
               "db_layout": t.get("layout", "per_species"), "db_stamp_version": DB_STAMP_VERSION,
               "species_table_sha256": _sha256_file(_abs(pcfg["panel"]["species_table"]))}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def sequences_hash(seqs: list[str]) -> str:
    h = hashlib.sha256()
    for s in seqs:
        h.update(s.encode())
        h.update(b"\n")
    return h.hexdigest()


# ----------------------------------------------------------------------------------------------
# per-species databases
# ----------------------------------------------------------------------------------------------

def _split_panel_fasta(fasta: Path, codes: list[str], outdir: Path) -> dict[str, int]:
    """Write <outdir>/<code>/seqs.fasta with headers = within-species index. Returns counts."""
    handles, counts, expect = {}, {c: 0 for c in codes}, {c: 0 for c in codes}
    for c in codes:
        (outdir / c).mkdir(parents=True, exist_ok=True)
        handles[c] = open(outdir / c / "seqs.fasta.part", "w")
    cur = None
    try:
        with open(fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    code, idx = line[1:].strip().split("|")
                    if code not in handles:
                        raise PanelError(f"panel FASTA species {code!r} not in the species table")
                    if int(idx) != expect[code]:
                        raise PanelError(f"{code}: indices are not 0..n-1 in file order")
                    expect[code] += 1
                    counts[code] += 1
                    cur = handles[code]
                    cur.write(f">{idx}\n")
                else:
                    cur.write(line)
    finally:
        for h in handles.values():
            h.close()
    for c in codes:
        (outdir / c / "seqs.fasta.part").replace(outdir / c / "seqs.fasta")
    return counts


def build_dbs(pcfg: dict, mm, force: bool = False) -> dict:
    """Build (once) one MMseqs2 DB per panel species. ``mm`` is a ``search.MMseqs`` runner.
    Returns the stamp dict. The stamp records the panel FASTA sha256, the species list, the
    MMseqs2 binary sha256 and per-species protein counts; a mismatch triggers a rebuild."""
    db_dir = _abs(pcfg["tier_b"]["db_dir"])
    fasta = _abs(pcfg["panel"]["fasta"])
    codes = species_table(pcfg)["code"].tolist()
    stamp_path = db_dir / "stamp.json"
    tool_sha = _sha256_file(Path(mm.binary))
    if stamp_path.exists() and not force:
        st = json.loads(stamp_path.read_text())
        ok = (st.get("version") == DB_STAMP_VERSION and st.get("codes") == codes
              and st.get("panel_fasta_size") == fasta.stat().st_size
              and st.get("mmseqs_sha256") == tool_sha
              and all((db_dir / c / "db.dbtype").exists() for c in codes))
        if ok:
            return st
    t0 = time.time()
    db_dir.mkdir(parents=True, exist_ok=True)
    fasta_sha = _sha256_file(fasta)
    counts = _split_panel_fasta(fasta, codes, db_dir)
    for c in codes:
        d = db_dir / c
        for old in d.glob("db*"):
            old.unlink()
        mm.run("createdb", str(d / "seqs.fasta"), str(d / "db"), "--shuffle", "0",
               "--dbtype", "1", threads=False)
        lk = pd.read_csv(d / "db.lookup", sep="\t", header=None, dtype=str)
        if not (lk[0].astype(np.int64).to_numpy() == lk[1].astype(np.int64).to_numpy()).all():
            raise PanelError(f"{c}: database keys differ from header indices")
        if len(lk) != counts[c]:
            raise PanelError(f"{c}: createdb kept {len(lk)} of {counts[c]} sequences")
    st = {"version": DB_STAMP_VERSION, "codes": codes, "counts": counts,
          "panel_fasta_sha256": fasta_sha, "panel_fasta_size": fasta.stat().st_size,
          "mmseqs_sha256": tool_sha, "build_s": round(time.time() - t0, 1)}
    stamp_path.write_text(json.dumps(st, indent=1))
    return st


def species_db(pcfg: dict, code: str) -> Path:
    return _abs(pcfg["tier_b"]["db_dir"]) / code / "db"


def read_species_seqs(pcfg: dict, code: str, wanted: np.ndarray | None = None) -> dict[int, str]:
    """Sequences of species ``code`` by within-species index (only ``wanted`` if given)."""
    want = None if wanted is None else set(int(x) for x in wanted)
    out: dict[int, str] = {}
    cur, chunks = None, []
    with open(_abs(pcfg["tier_b"]["db_dir"]) / code / "seqs.fasta") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(chunks)
                k = int(line[1:])
                cur, chunks = (k if want is None or k in want else None), []
            elif cur is not None:
                chunks.append(line.strip())
    if cur is not None:
        out[cur] = "".join(chunks)
    return out


# ----------------------------------------------------------------------------------------------
# search
# ----------------------------------------------------------------------------------------------

def _read_hits(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame({c: pd.Series(dtype=(object if t is str else t)) for c, t in _DTYPES.items()})
    return pd.read_csv(path, sep="\t", header=None, names=PANEL_COLUMNS, dtype=_DTYPES,
                       keep_default_na=False)


def _search_round(mm, ccfg: dict, scfg: dict, qdb: str, tdb: str, tag: str, cap: int):
    pref, aln = f"pref_{tag}", f"aln_{tag}"
    mm.run("prefilter", qdb, tdb, pref, "--max-seqs", str(cap), "-s", str(scfg["sensitivity"]),
           *prefilter_scoring_args(ccfg))
    mm.run("result2stats", qdb, tdb, pref, f"{pref}_count", "--stat", "linecount")
    mm.run("createtsv", qdb, qdb, f"{pref}_count", f"{pref}_count.tsv")
    pc = pd.read_csv(mm.workdir / f"{pref}_count.tsv", sep="\t", header=None, names=["q", "n"])
    mm.run("align", qdb, tdb, pref, aln, "-a", "--alignment-mode",
           str(ccfg["mmseqs"]["alignment_mode"]), "-e", str(scfg["evalue"]),
           *align_scoring_args(ccfg))
    mm.run("convertalis", qdb, tdb, aln, f"{aln}.tsv", "--format-output", ",".join(PANEL_COLUMNS))
    hits = _read_hits(mm.workdir / f"{aln}.tsv")
    hits = hits[hits["bits"] >= float(scfg["keep_min_bits"])]
    for f in mm.workdir.glob(f"{pref}*"):
        f.unlink()
    for f in mm.workdir.glob(f"{aln}*"):
        f.unlink()
    return pc, hits


def _capped(hits: pd.DataFrame, presence: float, cap: int, frac: float) -> set[int]:
    q = hits.loc[hits["bits"] >= presence, "query"].to_numpy()
    u, c = np.unique(q, return_counts=True)
    return set(u[c >= frac * cap].tolist())


def search_species(mm, ccfg: dict, pcfg: dict, qdb: str, code: str, n_target: int):
    """One per-species search with saturation re-search. Returns (hits, log)."""
    t = pcfg["tier_b"]
    scfg = t["search"]
    presence = float(t["presence_min_bits"])
    frac = float(scfg["saturation_fraction"])
    cap = int(scfg["max_seqs"])
    tdb = str(species_db(pcfg, code))
    pc, hits = _search_round(mm, ccfg, scfg, qdb, tdb, f"{code}_r0", cap)
    capped = _capped(hits, presence, cap, frac)
    rounds = [{"round": 0, "cap": cap, "n_capped": len(capped),
               "n_prefilter_saturated": int((pc["n"] >= frac * cap).sum())}]
    for r in range(1, int(scfg["max_research_rounds"]) + 1):
        if not capped or cap >= n_target:
            break
        cap = min(cap * int(scfg["research_cap_factor"]), n_target)
        keyfile = mm.workdir / f"capped_{code}_r{r}.txt"
        keyfile.write_text("".join(f"{k}\n" for k in sorted(capped)))
        sub = str(mm.workdir / f"{Path(qdb).name}_{code}_capped_r{r}")
        mm.run("createsubdb", str(keyfile), qdb, sub, "--subdb-mode", "1", threads=False)
        spc, sub_hits = _search_round(mm, ccfg, scfg, sub, tdb, f"{code}_r{r}", cap)
        hits = pd.concat([hits[~hits["query"].isin(capped)], sub_hits], ignore_index=True)
        capped = _capped(sub_hits, presence, cap, frac) if cap < n_target else set()
        rounds.append({"round": r, "cap": cap, "n_queries": int(len(spc)),
                       "n_capped": len(capped),
                       "n_prefilter_saturated": int((spc["n"] >= frac * cap).sum())})
        for f in mm.workdir.glob(f"{Path(sub).name}*"):
            f.unlink()
    hits = hits.sort_values(["query", "bits", "target"], ascending=[True, False, True],
                            kind="mergesort").reset_index(drop=True)
    return hits, {"rounds": rounds, "n_capped_final": len(capped), "n_hits": int(len(hits))}


def search_panel(mm, ccfg: dict, pcfg: dict, qdb: str, seqs: list[str],
                 cache_root: Path | None = None, species_jobs: int = 1) -> tuple[dict[str, Path], dict]:
    """Search the query DB (``qdb``, keys 0..n-1 in canonical order) against every panel species.
    Returns ({code: hit parquet}, log). Cached per (query sequences, search settings).
    ``species_jobs`` > 1 runs that many species searches concurrently (each with ``mm.threads``
    threads, in its own sub-directory); results do not depend on it."""
    st = build_dbs(pcfg, mm)
    codes = st["codes"]
    tool_sha = st["mmseqs_sha256"]
    key = (sequences_hash(seqs)[:16] + "_" + search_settings_hash(pcfg, ccfg, tool_sha)[:12])
    root = Path(cache_root) if cache_root is not None else _abs(pcfg["tier_b"]["cache_dir"])
    cdir = root / key
    cdir.mkdir(parents=True, exist_ok=True)
    paths, slog = {}, {}
    t_all = time.time()
    n_cached = 0
    qdb_abs = str(Path(mm.workdir) / qdb) if not Path(qdb).is_absolute() else qdb
    todo = []
    for code in codes:
        out = cdir / f"{code}.parquet"
        lg = cdir / f"{code}.json"
        if out.exists() and lg.exists():
            paths[code] = out
            slog[code] = json.loads(lg.read_text())
            n_cached += 1
        else:
            todo.append(code)

    def one(code):
        out = cdir / f"{code}.parquet"
        sub = Path(mm.workdir) / f"sp_{code}"
        sub.mkdir(parents=True, exist_ok=True)
        smm = type(mm)(binary=mm.binary, threads=mm.threads, workdir=sub)
        t0 = time.time()
        hits, info = search_species(smm, ccfg, pcfg, qdb_abs, code, int(st["counts"][code]))
        info["search_s"] = round(time.time() - t0, 2)
        hits.to_parquet(out.with_suffix(".parquet.part"), index=False)
        out.with_suffix(".parquet.part").replace(out)
        (cdir / f"{code}.json").write_text(json.dumps(info))
        shutil.rmtree(sub, ignore_errors=True)
        log.info("panel %s: %d hits in %.1fs", code, info["n_hits"], info["search_s"])
        return code, out, info

    if species_jobs > 1 and len(todo) > 1:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=int(species_jobs)) as ex:
            results = list(ex.map(one, todo))
    else:
        results = [one(c) for c in todo]
    for code, out, info in results:
        paths[code] = out
        slog[code] = info
    paths = {c: paths[c] for c in codes}
    slog = {c: slog[c] for c in codes}
    meta = {"cache_key": key, "cache_dir": str(cdir), "n_species_cached": n_cached,
            "search_wall_s": round(time.time() - t_all, 1), "per_species": slog,
            "n_capped_first_round_total": int(sum(v["rounds"][0]["n_capped"] for v in slog.values())),
            "n_capped_final_total": int(sum(v["n_capped_final"] for v in slog.values())),
            "db_stamp": {k: st[k] for k in ("panel_fasta_sha256", "mmseqs_sha256", "version")}}
    return paths, meta


def load_hits(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    return pd.read_parquet(path, columns=columns)


# ----------------------------------------------------------------------------------------------
# reduction to per-(query, species) summaries
# ----------------------------------------------------------------------------------------------

@dataclass
class PanelData:
    """Per-(query, species) summaries; matrices are (n_query, n_species)."""
    codes: list[str]
    best_bits: np.ndarray      # float32, 0 where no hit >= presence bits
    best_target: np.ndarray    # int32, -1 where absent
    best_ident: np.ndarray     # float32 nident/alnlen of the best hit, NaN where absent
    best_tlen: np.ndarray      # int32 length of the best hit, 0 where absent
    second_bits: np.ndarray    # float32, second-best hit bits (>= presence) or 0
    rbh: np.ndarray            # bool pseudo-RBH (change 5.1)
    dup_count: np.ndarray      # int16 hits with bits >= dupability_bsr x self-bits (>= presence)
    hit_paths: dict[str, Path] = field(default_factory=dict)


def reduce_species(hits: pd.DataFrame, n: int, selfbits: np.ndarray, tcfg: dict) -> dict:
    presence = float(tcfg["presence_min_bits"])
    mincov = float(tcfg["rbh_min_cov"])
    bsr = float(tcfg["dupability_bsr"])
    P = hits[hits["bits"] >= presence]
    out = {"best_bits": np.zeros(n, np.float32), "best_target": np.full(n, -1, np.int32),
           "best_ident": np.full(n, np.nan, np.float32), "best_tlen": np.zeros(n, np.int32),
           "second_bits": np.zeros(n, np.float32), "rbh": np.zeros(n, bool),
           "dup_count": np.zeros(n, np.int16)}
    if len(P) == 0:
        return out
    P = P.sort_values(["query", "bits", "target"], ascending=[True, False, True], kind="mergesort")
    q = P["query"].to_numpy(np.int64)
    first = np.r_[True, q[1:] != q[:-1]]
    pos = np.arange(len(P)) - np.maximum.accumulate(np.where(first, np.arange(len(P)), 0))
    B = P[first]
    bq = B["query"].to_numpy(np.int64)
    out["best_bits"][bq] = B["bits"].to_numpy(np.float32)
    out["best_target"][bq] = B["target"].to_numpy(np.int32)
    out["best_ident"][bq] = (B["nident"].to_numpy(np.float64) / B["alnlen"].to_numpy(np.float64))
    out["best_tlen"][bq] = B["tlen"].to_numpy(np.int32)
    S2 = P[pos == 1]
    out["second_bits"][S2["query"].to_numpy(np.int64)] = S2["bits"].to_numpy(np.float32)
    # dup count: hits with bits >= bsr x self-bits of the query
    sb = selfbits[q]
    ok = P["bits"].to_numpy(np.float64) >= bsr * sb
    out["dup_count"][:] = np.minimum(np.bincount(q[ok], minlength=n), np.iinfo(np.int16).max)
    # pseudo-RBH
    unique_q = out["best_bits"][bq] > out["second_bits"][bq]
    T = P.sort_values(["target", "bits", "query"], ascending=[True, False, True], kind="mergesort")
    tt = T["target"].to_numpy(np.int64)
    tfirst = np.r_[True, tt[1:] != tt[:-1]]
    tpos = np.arange(len(T)) - np.maximum.accumulate(np.where(tfirst, np.arange(len(T)), 0))
    top = T[tfirst]
    sec = T[tpos == 1]
    top_q = pd.Series(top["query"].to_numpy(), index=top["target"].to_numpy())
    top_b = pd.Series(top["bits"].to_numpy(), index=top["target"].to_numpy())
    sec_b = pd.Series(sec["bits"].to_numpy(), index=sec["target"].to_numpy())
    bt = B["target"].to_numpy(np.int64)
    tq = top_q.reindex(bt).to_numpy()
    unique_t = top_b.reindex(bt).to_numpy() > sec_b.reindex(bt).fillna(-1.0).to_numpy()
    qcov = (B["qend"].to_numpy() - B["qstart"].to_numpy() + 1) / B["qlen"].to_numpy()
    tcov = (B["tend"].to_numpy() - B["tstart"].to_numpy() + 1) / B["tlen"].to_numpy()
    rbh = unique_q & (tq == bq) & unique_t & (qcov >= mincov) & (tcov >= mincov)
    out["rbh"][bq] = rbh
    return out


def load_panel_data(hit_paths: dict[str, Path], n: int, selfbits: np.ndarray, pcfg: dict) -> PanelData:
    codes = list(hit_paths)
    S = len(codes)
    mats = {"best_bits": np.zeros((n, S), np.float32), "best_target": np.full((n, S), -1, np.int32),
            "best_ident": np.full((n, S), np.nan, np.float32), "best_tlen": np.zeros((n, S), np.int32),
            "second_bits": np.zeros((n, S), np.float32), "rbh": np.zeros((n, S), bool),
            "dup_count": np.zeros((n, S), np.int16)}
    cols = [c for c in PANEL_COLUMNS if c != "cigar"]
    for k, code in enumerate(codes):
        r = reduce_species(load_hits(hit_paths[code], cols), n, selfbits, pcfg["tier_b"])
        for name, v in r.items():
            mats[name][:, k] = v
    return PanelData(codes=codes, hit_paths=dict(hit_paths), **mats)


# ----------------------------------------------------------------------------------------------
# calibration and masks
# ----------------------------------------------------------------------------------------------

def kimura(identity, cap: float = 3.0):
    """Kimura protein distance -ln(1 - p - 0.2 p^2), p = 1 - identity; cap where the argument is
    <= 0, and clipped at ``cap`` (keeps d monotone in p). NaN identity -> NaN."""
    p = 1.0 - np.asarray(identity, dtype=np.float64)
    arg = 1.0 - p - 0.2 * p * p
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.where(arg > 0, -np.log(np.where(arg > 0, arg, 1.0)), cap)
    d = np.minimum(d, cap)
    return np.where(np.isnan(p), np.nan, d)


@dataclass
class Calibration:
    codes: list[str]
    d: np.ndarray               # d_s per species (NaN: excluded or uncalibrated)
    usable: np.ndarray          # bool: retained after self-exclusion and calibrated
    self_excluded: np.ndarray   # bool
    median_rbh_ident: np.ndarray
    n_rbh: np.ndarray
    n_core_rbh: np.ndarray
    core_size: int
    core_frac_used: float
    log: dict


def calibrate(pdata: PanelData, pcfg: dict) -> Calibration:
    t = pcfg["tier_b"]
    S = len(pdata.codes)
    rbh = pdata.rbh
    ident = pdata.best_ident.astype(np.float64)
    n_rbh = rbh.sum(0)
    med = np.full(S, np.nan)
    for s in range(S):
        if n_rbh[s]:
            med[s] = float(np.median(ident[rbh[:, s], s]))
    thr = float(t["self_exclusion_median_identity"])
    self_ex = np.nan_to_num(med, nan=0.0) >= thr
    retained = (~self_ex) & (n_rbh > 0)
    c = t["core"]
    frac = float(c["min_species_frac"])
    step = float(c["fallback_step"])
    min_size = int(c["min_size"])
    n_ret = int(retained.sum())
    cnt = rbh[:, retained].sum(1)
    tried = []
    while True:
        core = cnt >= frac * n_ret - 1e-9 if n_ret else np.zeros(len(cnt), bool)
        tried.append({"frac": round(frac, 4), "size": int(core.sum())})
        if core.sum() >= min_size or frac - step < step - 1e-9:
            break
        frac = round(frac - step, 10)
    d = np.full(S, np.nan)
    n_core = np.zeros(S, np.int64)
    for s in np.flatnonzero(retained):
        sel = core & rbh[:, s]
        n_core[s] = int(sel.sum())
        if n_core[s]:
            d[s] = float(kimura(np.median(ident[sel, s]), float(t["kimura"]["cap"])))
    usable = retained & np.isfinite(d)
    logd = {
        "self_excluded": [pdata.codes[s] for s in np.flatnonzero(self_ex)],
        "no_rbh": [pdata.codes[s] for s in np.flatnonzero((~self_ex) & (n_rbh == 0))],
        "uncalibrated": [pdata.codes[s] for s in np.flatnonzero(retained & ~np.isfinite(d))],
        "n_retained": n_ret, "n_usable": int(usable.sum()),
        "core_size": int(core.sum()), "core_frac_used": frac, "core_thresholds_tried": tried,
        "per_species": {
            pdata.codes[s]: {"d_s": None if not np.isfinite(d[s]) else round(float(d[s]), 6),
                             "median_rbh_identity": None if not np.isfinite(med[s]) else round(float(med[s]), 6),
                             "n_rbh": int(n_rbh[s]), "n_core_rbh": int(n_core[s]),
                             "rbh_fraction": round(float(n_rbh[s]) / rbh.shape[0], 6) if rbh.shape[0] else None,
                             "core_rbh_fraction": (round(float(n_core[s]) / core.sum(), 6)
                                                   if core.sum() else None),
                             "self_excluded": bool(self_ex[s])}
            for s in range(S)},
    }
    return Calibration(codes=list(pdata.codes), d=d, usable=usable, self_excluded=self_ex,
                       median_rbh_ident=med, n_rbh=n_rbh, n_core_rbh=n_core,
                       core_size=int(core.sum()), core_frac_used=frac, log=logd)


QC_NEAREST_D_MAX = 0.3   # nearest retained species farther than this -> "nearest_species_far" (worm-like)


def qc_flags(cal: Calibration, pcfg: dict, nearest_d_max: float = QC_NEAREST_D_MAX,
             min_core: int | None = None) -> dict:
    """Proteome-level QC for Tier B: output metadata only, NEVER a feature column (constant within
    a proteome). Evaluated on the natural panel (``cal``, after self-exclusion and the core-set
    fallback), independent of any ``--panel-dmin`` mask.

    * ``core_set_small``: |G*| < ``min_core`` (default ``tier_b.core.min_size``) after the fallback;
    * ``nearest_species_far``: smallest d_s over usable species > ``nearest_d_max``;
    * ``no_usable_species``: no calibrated species at all (then the nearest-species flag is False).
    """
    if min_core is None:
        min_core = int(pcfg["tier_b"]["core"]["min_size"])
    d = cal.d[cal.usable]
    k = int(np.argmin(np.where(cal.usable, cal.d, np.inf))) if d.size else None
    d_nearest = None if k is None else float(cal.d[k])
    flags = {
        "core_set_small": bool(cal.core_size < min_core),
        "nearest_species_far": bool(d_nearest is not None and d_nearest > nearest_d_max),
        "no_usable_species": bool(d.size == 0),
    }
    warnings = []
    if flags["core_set_small"]:
        warnings.append(f"core set |G*|={cal.core_size} < {min_core} after fallback "
                        f"(fraction used {cal.core_frac_used}): d_s unreliable")
    if flags["nearest_species_far"]:
        warnings.append(f"nearest retained panel species {cal.codes[k]} has d_s={d_nearest:.3f} > "
                        f"{nearest_d_max}: young duplications unresolvable, depth often NaN")
    if flags["no_usable_species"]:
        warnings.append("no usable panel species after self-exclusion/calibration")
    return {"flags": flags, "warnings": warnings, "any": bool(warnings),
            "core_size": int(cal.core_size), "core_min_size": int(min_core),
            "d_nearest": None if d_nearest is None else round(d_nearest, 6),
            "nearest_species": None if k is None else cal.codes[k],
            "nearest_d_max": float(nearest_d_max)}


@dataclass
class Mask:
    dmin: float | None
    cols: np.ndarray       # species columns used (indices into PanelData.codes), ordered by (d, column)
    d: np.ndarray          # d_s of those columns
    bins: np.ndarray       # 0-based bin index per used column
    weight: np.ndarray     # 1 / n_bin per used column (bin counts over used columns)
    log: dict


def assign_bins(d: np.ndarray, edges: list[float]) -> np.ndarray:
    """0-based bin index of each d. The edges are human-derived and used only for the discordance
    weights (see module docstring); they carry no other meaning."""
    e = np.asarray(edges, dtype=np.float64)
    return np.searchsorted(e, np.asarray(d, dtype=np.float64), side="right") - 1


def apply_mask(cal: Calibration, pcfg: dict, dmin: float | None = None) -> Mask:
    """Species used for the features (optionally d_s >= dmin) and their bins/weights. The bin edges
    are human-derived and used only for the discordance weights 1/n_bin (module docstring)."""
    t = pcfg["tier_b"]
    use = cal.usable.copy()
    if dmin is not None:
        use &= cal.d >= float(dmin)
    cols = np.flatnonzero(use)
    order = np.lexsort((cols, cal.d[cols]))
    cols = cols[order]
    d = cal.d[cols]
    bins = assign_bins(d, t["bins"]["edges"])
    nb = np.bincount(bins, minlength=len(t["bins"]["edges"]))
    w = 1.0 / nb[bins] if len(bins) else np.zeros(0)
    labels = t["bins"]["labels"]
    lg = {"dmin": dmin, "n_species": int(len(cols)),
          "removed_by_mask": [cal.codes[s] for s in np.flatnonzero(cal.usable & ~use)],
          "bin_counts": {labels[i]: int(nb[i]) for i in range(len(labels))},
          "d_nearest": None if not len(d) else round(float(d.min()), 6),
          "nearest_species": None if not len(d) else cal.codes[cols[0]]}
    return Mask(dmin=dmin, cols=cols, d=d, bins=bins, weight=w, log=lg)
