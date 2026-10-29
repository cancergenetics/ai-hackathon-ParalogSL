"""MMseqs2 wrappers for the self-search (DESIGN v3 §3; algorithms_comparative §2.2).

All settings come from ``config/caller.yaml``. The input database is written by the pipeline with
headers ``0..n-1`` (canonical order, see ``fasta.canonical_order``) and ``createdb --shuffle 0``, so
database key == header name == internal protein index; ``_check_lookup`` asserts this.

Steps:
* ``all_vs_all``: ``prefilter`` (``-s``, ``--max-seqs``) then ``align -a --alignment-mode 3``.
  Running the two modules explicitly (what ``mmseqs search`` does for a single sensitivity step)
  gives direct access to the prefilter database; per-query candidate counts are logged
  (``result2stats --stat linecount``). Queries whose aligned hits passing the calling criterion
  reach saturation_fraction * cap are re-searched (``createsubdb``) with a larger cap and their
  rows replaced (rule of PI decision D11; see ``all_vs_all``).
* ``realign``: ``align`` on an explicit pair list (``tsv2db`` -> prefilter-format DB, type 7) with
  ``-e inf``, so every requested ordered pair (both orientations) and every self pair (explicit
  self-scores) gets backtraced statistics.
"""

from __future__ import annotations

import hashlib
import logging
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .config import PROJECT_ROOT

log = logging.getLogger(__name__)

ALN_COLUMNS = [
    "query", "target", "bits", "nident", "alnlen", "qstart", "qend",
    "tstart", "tend", "qlen", "tlen", "fident",
]
_ALN_DTYPES = {
    "query": np.int64, "target": np.int64, "bits": np.float64, "nident": np.int64,
    "alnlen": np.int64, "qstart": np.int64, "qend": np.int64, "tstart": np.int64,
    "tend": np.int64, "qlen": np.int64, "tlen": np.int64, "fident": np.float64,
}


class MMseqsError(RuntimeError):
    pass


def resolve_binary(cfg: dict) -> Path:
    """Configured MMseqs2 binary (relative to the project root); falls back to `mmseqs` on PATH.

    The version actually used is recorded by tool_info() in every output's metadata.
    """
    p = Path(cfg["mmseqs"]["binary"])
    p = p if p.is_absolute() else PROJECT_ROOT / p
    if not p.exists() and (on_path := shutil.which("mmseqs")):
        return Path(on_path)
    return p


def tool_info(binary: Path) -> dict:
    out = subprocess.run([str(binary), "version"], capture_output=True, text=True, check=True)
    h = hashlib.sha256(Path(binary).read_bytes()).hexdigest()
    return {"mmseqs_version": out.stdout.strip(), "mmseqs_sha256": h, "mmseqs_path": str(binary)}


def scoring_args(cfg: dict) -> list[str]:
    s = cfg["mmseqs"]["scoring"]
    return [
        "--sub-mat", s["sub_mat"],
        "--gap-open", s["gap_open"],
        "--gap-extend", s["gap_extend"],
        "--comp-bias-corr", str(s["comp_bias_corr"]),
        "--comp-bias-corr-scale", str(s["comp_bias_corr_scale"]),
        "--max-seq-len", str(s["max_seq_len"]),
    ]


def prefilter_args(cfg: dict) -> list[str]:
    s = cfg["mmseqs"]["scoring"]
    return [
        "-s", str(cfg["mmseqs"]["sensitivity"]),
        "--seed-sub-mat", s["seed_sub_mat"],
        "--mask", str(s["mask"]),
        "--mask-prob", str(s["mask_prob"]),
        "--mask-lower-case", str(s["mask_lower_case"]),
    ]


@dataclass
class MMseqs:
    binary: Path
    threads: int
    workdir: Path
    commands: list[list[str]] = field(default_factory=list)

    def run(self, *args: str, threads: bool = True) -> None:
        cmd = [str(self.binary), *args, "-v", "1"]
        if threads:
            cmd += ["--threads", str(self.threads)]
        self.commands.append(cmd)
        log.info("mmseqs %s", " ".join(args[:1]))
        res = subprocess.run(cmd, cwd=self.workdir, capture_output=True, text=True)
        if res.returncode != 0:
            raise MMseqsError(f"{' '.join(cmd)}\n{res.stdout[-2000:]}\n{res.stderr[-2000:]}")


def createdb(mm: MMseqs, fasta: Path, db: str) -> None:
    mm.run("createdb", str(fasta), db, "--shuffle", "0", "--dbtype", "1", threads=False)
    _check_lookup(mm.workdir / f"{db}.lookup")


def _check_lookup(path: Path) -> None:
    lk = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if not (lk[0].astype(np.int64).to_numpy() == lk[1].astype(np.int64).to_numpy()).all():
        raise MMseqsError(f"{path}: database keys differ from header indices")


def lowcomp_fractions(mm: MMseqs, cfg: dict, db: str, n: int) -> np.ndarray:
    """Fraction of residues soft-masked by MMseqs2 tantan (``masksequence``) per protein.
    Same masker and probability threshold as the search prefilter."""
    s = cfg["mmseqs"]["scoring"]
    mm.run("masksequence", db, f"{db}_masked", "--mask-prob", str(s["mask_prob"]))
    mm.run("convert2fasta", f"{db}_masked", f"{db}_masked.fasta", threads=False)
    frac = np.full(n, np.nan)
    cur, chunks = None, []

    def _flush():
        if cur is not None:
            seq = "".join(chunks)
            frac[cur] = sum(c.islower() for c in seq) / len(seq) if seq else np.nan

    with open(mm.workdir / f"{db}_masked.fasta") as fh:
        for line in fh:
            if line.startswith(">"):
                _flush()
                cur, chunks = int(line[1:].split()[0]), []
            else:
                chunks.append(line.strip())
    _flush()
    if np.isnan(frac).any():
        raise MMseqsError("masksequence output is missing proteins")
    return frac


def _read_tsv_counts(path: Path) -> dict[int, int]:
    df = pd.read_csv(path, sep="\t", header=None, names=["q", "n"], dtype=np.int64)
    return dict(zip(df.q.tolist(), df.n.tolist()))


def read_alignments(path: Path) -> pd.DataFrame:
    if path.stat().st_size == 0:
        return pd.DataFrame({c: pd.Series(dtype=t) for c, t in _ALN_DTYPES.items()})
    return pd.read_csv(path, sep="\t", header=None, names=ALN_COLUMNS, dtype=_ALN_DTYPES)


def _search_round(mm: MMseqs, cfg: dict, qdb: str, tdb: str, tag: str, cap: int):
    pref, aln = f"pref_{tag}", f"aln_{tag}"
    mm.run("prefilter", qdb, tdb, pref, "--max-seqs", str(cap),
           *prefilter_args(cfg), *_prefilter_scoring(cfg))
    mm.run("result2stats", qdb, tdb, pref, f"{pref}_count", "--stat", "linecount")
    mm.run("createtsv", qdb, qdb, f"{pref}_count", f"{pref}_count.tsv")
    counts = _read_tsv_counts(mm.workdir / f"{pref}_count.tsv")
    mm.run("align", qdb, tdb, pref, aln, "-a", "--alignment-mode",
           str(cfg["mmseqs"]["alignment_mode"]), "-e", str(cfg["mmseqs"]["search_evalue"]),
           *scoring_args(cfg))
    mm.run("convertalis", qdb, tdb, aln, f"{aln}.tsv", "--format-output", ",".join(ALN_COLUMNS))
    return counts, read_alignments(mm.workdir / f"{aln}.tsv")


def _prefilter_scoring(cfg: dict) -> list[str]:
    s = cfg["mmseqs"]["scoring"]
    return [
        "--sub-mat", s["sub_mat"],
        "--comp-bias-corr", str(s["comp_bias_corr"]),
        "--comp-bias-corr-scale", str(s["comp_bias_corr_scale"]),
        "--max-seq-len", str(s["max_seq_len"]),
    ]


def all_vs_all(mm: MMseqs, cfg: dict, db: str, n: int) -> tuple[pd.DataFrame, dict]:
    """Self-search with saturation detection and re-search. Returns (hits, saturation log).

    Saturation rule (PI decision D11): a query is capped when its number of aligned hits passing
    the calling criterion (E' <= evalue_prime_max, self hit included, as the cap also counts it)
    is >= saturation_fraction * cap. Only those queries are re-searched with a larger cap.
    Prefilter candidate counts are still logged (``n_prefilter_saturated``) but not used.
    """
    m = cfg["mmseqs"]
    if m.get("saturation_rule") != "called_hits":
        raise ValueError("only saturation_rule: called_hits is implemented (D11)")
    frac = float(m["saturation_fraction"])
    cap = int(m["max_seqs"])
    pcounts, hits = _search_round(mm, cfg, db, db, "r0", cap)
    capped = capped_queries(hits, cfg, cap)
    rounds = [{"round": 0, "cap": cap, "n_queries": n, "n_capped": len(capped),
               "n_prefilter_saturated": len(_capped(pcounts, cap, frac))}]
    first_capped = sorted(capped)
    for r in range(1, int(m["max_research_rounds"]) + 1):
        if not capped or cap >= n:
            break
        cap = min(cap * int(m["research_cap_factor"]), n)
        keyfile = mm.workdir / f"capped_r{r}.txt"
        keyfile.write_text("".join(f"{k}\n" for k in sorted(capped)))
        sub = f"{db}_capped_r{r}"
        mm.run("createsubdb", keyfile.name, db, sub, "--subdb-mode", "1", threads=False)
        sub_pcounts, sub_hits = _search_round(mm, cfg, sub, db, f"r{r}", cap)
        hits = pd.concat([hits[~hits["query"].isin(capped)], sub_hits], ignore_index=True)
        capped = capped_queries(sub_hits, cfg, cap) if cap < n else set()
        rounds.append({"round": r, "cap": cap, "n_queries": len(sub_pcounts),
                       "n_capped": len(capped),
                       "n_prefilter_saturated": len(_capped(sub_pcounts, cap, frac))})
    hits = hits.sort_values(["query", "target"], kind="mergesort").reset_index(drop=True)
    sat = {"rule": "called_hits", "rounds": rounds, "n_capped_first_round": len(first_capped),
           "n_capped_final": len(capped), "capped_first_round": first_capped}
    return hits, sat


def capped_queries(hits: pd.DataFrame, cfg: dict, cap: int) -> set[int]:
    """Queries whose aligned hits with E' <= threshold (self included) reach the saturation level."""
    from .caller import evalue_prime  # local import: caller does not depend on search

    c = cfg["caller"]
    ep = evalue_prime(hits["qlen"].to_numpy(), hits["bits"].to_numpy(), float(c["N_ref"]))
    ok = hits["query"].to_numpy()[ep <= float(c["evalue_prime_max"])]
    q, cnt = np.unique(ok, return_counts=True)
    return set(q[cnt >= float(cfg["mmseqs"]["saturation_fraction"]) * cap].tolist())


def _capped(counts: dict[int, int], cap: int, frac: float) -> set[int]:
    return {q for q, c in counts.items() if c >= frac * cap}


def realign(mm: MMseqs, cfg: dict, db: str, q: np.ndarray, t: np.ndarray) -> pd.DataFrame:
    """Backtraced alignment of each ordered pair (q[i] -> t[i]); include (i, i) for self-scores."""
    order = np.lexsort((t, q))
    q, t = np.asarray(q)[order], np.asarray(t)[order]
    tsv = mm.workdir / "pairs.tsv"
    pd.DataFrame({"q": q, "t": t, "s": 0, "d": 0}).to_csv(tsv, sep="\t", header=False, index=False)
    mm.run("tsv2db", tsv.name, "pairs_db", "--output-dbtype", "7", threads=False)
    mm.run("align", db, db, "pairs_db", "realn", "-a", "--alignment-mode",
           str(cfg["mmseqs"]["alignment_mode"]), "-e", str(cfg["mmseqs"]["realign_evalue"]),
           *scoring_args(cfg))
    mm.run("convertalis", db, db, "realn", "realn.tsv", "--format-output", ",".join(ALN_COLUMNS))
    out = read_alignments(mm.workdir / "realn.tsv")
    # MMseqs2 may report >1 HSP per pair only with --alt-ali > 0 (default 0); keep the best anyway.
    out = (out.sort_values(["query", "target", "bits"], ascending=[True, True, False],
                           kind="mergesort")
              .drop_duplicates(["query", "target"], keep="first").reset_index(drop=True))
    return out


def cleanup(workdir: Path) -> None:
    shutil.rmtree(workdir, ignore_errors=True)
