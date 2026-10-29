"""Proteome FASTA -> called paralog pairs + Tier A features (+ placeholder score).

Output: one row per unordered pair, ``gene_a < gene_b`` lexicographically, registry columns in
registry order, optional ``sl_score``/``sl_rank`` placeholders (NaN until the model step exists).
Written as parquet (``.parquet``) or TSV (anything else). A sidecar ``<output>.meta.json`` holds the
schema hash, tool versions, all caller/search settings, input QC and run statistics; for parquet
the same metadata is also stored in the file's schema metadata (key ``slparalog``).

Schema hash = sha256 over (feature registry, config/caller.yaml content). It changes whenever a
feature definition, a search/scoring setting, N_ref or a feature threshold changes.
"""

from __future__ import annotations

import hashlib
import json
import logging
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__, caller, fasta, search
from .config import CONFIG_DIR, load_yaml
from .features import registry, tier_a

log = logging.getLogger(__name__)

DEFAULT_CALLER_CONFIG = CONFIG_DIR / "caller.yaml"


def schema_hash(cfg: dict, panel_settings_hash: str | None = None) -> str:
    d = {"registry": registry.registry_hash(), "caller_config": cfg}
    if panel_settings_hash is not None:  # Tier B on: panel settings + shared scoring + panel checksum
        d["panel_settings"] = panel_settings_hash
    payload = json.dumps(d, sort_keys=True, default=str).encode()
    return hashlib.sha256(payload).hexdigest()


def run(proteome: str | Path, output: str | Path | None = None, threads: int = 1,
        workdir: str | Path | None = None, keep_work: bool = False,
        features_only: bool = False, config_path: str | Path | None = None,
        tier_b: bool = False, panel_dmin: float | None = None, species_jobs: int = 1,
        panel_config_path: str | Path | None = None, panel_cache: str | Path | None = None,
        tier_c: bool = False, plm_cache: str | Path | None = None, plm_model: str | None = None,
        ) -> tuple[pd.DataFrame, dict]:
    t0 = time.time()
    if panel_dmin is not None and not tier_b:
        raise ValueError("panel_dmin (--panel-dmin) needs tier_b (--tier-b)")
    cfg = load_yaml(config_path or DEFAULT_CALLER_CONFIG)
    timings: dict[str, float] = {}
    seqs_by_id, qc = fasta.read_proteome(proteome)
    ids = fasta.canonical_order(seqs_by_id)
    seqs = [seqs_by_id[k] for k in ids]
    n = len(ids)

    # A user-given workdir is never deleted (it may hold other files); a temporary one is removed
    # unless keep_work is set.
    if workdir is None:
        wd = Path(tempfile.mkdtemp(prefix="slparalog_"))
        remove = not keep_work
    else:
        wd = Path(workdir).resolve()  # mmseqs runs with the workdir as cwd: relative paths would break
        wd.mkdir(parents=True, exist_ok=True)
        remove = False
    try:
        binary = search.resolve_binary(cfg)
        mm = search.MMseqs(binary=binary, threads=int(threads), workdir=wd)
        fasta.write_fasta(wd / "input.fasta", ((str(k), s) for k, s in enumerate(seqs)))
        t = time.time()
        search.createdb(mm, wd / "input.fasta", "db")
        lowcomp = search.lowcomp_fractions(mm, cfg, "db", n)
        timings["createdb_mask_s"] = time.time() - t
        t = time.time()
        hits, sat = search.all_vs_all(mm, cfg, "db", n)
        timings["search_s"] = time.time() - t
        pairs, call_stats = caller.call_pairs(hits, cfg)
        t = time.time()
        q = np.concatenate([pairs["i"], pairs["j"], np.arange(n)])
        tt = np.concatenate([pairs["j"], pairs["i"], np.arange(n)])
        aln = search.realign(mm, cfg, "db", q, tt)
        timings["realign_s"] = time.time() - t
        t = time.time()
        feats, fsum = tier_a.compute(pairs, aln, seqs, lowcomp, cfg, threads=threads)
        timings["features_s"] = time.time() - t
        tb_meta = None
        if tier_b:  # Tier B hook (src/slparalog/features/tier_b.py); --panel-dmin is a hidden option
            from . import panel as panel_mod
            from .features import tier_b as tier_b_mod

            t = time.time()
            pcfg = panel_mod.load_panel_config(panel_config_path)
            tb, tb_meta = tier_b_mod.run_stage(mm, cfg, pcfg, "db", seqs, pairs, aln, threads=threads,
                                               panel_dmin=panel_dmin, cache_root=panel_cache,
                                               species_jobs=species_jobs)
            feats = pd.concat([feats, tb], axis=1)
            timings["tier_b_s"] = time.time() - t
        tc_meta = None
        if tier_c:  # Tier C: PLM embedding in its own process (no torch here), then numpy features
            from .features import tier_c as tier_c_mod

            t = time.time()
            emb = tier_c_mod.embed_subprocess(seqs, wd / "plm", plm_cache or (wd / "plm_cache"), plm_model)
            timings["plm_embed_s"] = time.time() - t
            t = time.time()
            tc, tc_meta = tier_c_mod.compute(pairs, emb, threads=threads)
            tc_meta["plm_model"] = plm_model or "default"
            feats = pd.concat([feats, tc], axis=1)
            timings["tier_c_s"] = time.time() - t
        tool = search.tool_info(binary)
    finally:
        if remove:
            search.cleanup(wd)

    id_arr = np.array(ids, dtype=object)
    ga, gb = caller.canonical_gene_order(id_arr[pairs["i"].to_numpy()], id_arr[pairs["j"].to_numpy()])
    out = pd.concat([pd.DataFrame({"gene_a": ga, "gene_b": gb}), feats.reset_index(drop=True)], axis=1)
    if not features_only:
        out["sl_score"] = np.nan  # model step not implemented yet (M3)
        out["sl_rank"] = np.nan
    out = out.sort_values(["gene_a", "gene_b"], kind="mergesort").reset_index(drop=True)

    tools_cfg = load_yaml(CONFIG_DIR / "tools.yaml")
    meta = {
        "slparalog_version": __version__,
        "schema_hash": schema_hash(cfg, tb_meta["settings_hash"] if tb_meta else None),
        "registry_hash": registry.registry_hash(),  # provenance only; compatibility is decided by columns_hash
        "columns_hash": registry.table_columns_hash(out.columns),
        "caller_config": cfg,
        "tool": tool,
        "tool_matches_config": (tool["mmseqs_version"] == str(tools_cfg["mmseqs"]["version"])
                                and tool["mmseqs_sha256"] == tools_cfg["mmseqs"]["sha256"]),
        "input": {"path": str(proteome), "sha256": _sha256(proteome), **qc.as_dict()},
        "n_proteins": n,
        "threads": int(threads),
        "features_only": bool(features_only),
        "caller": call_stats,
        "saturation": {k: v for k, v in sat.items() if k != "capped_first_round"},
        "n_search_hits_nonself": int((hits["query"] != hits["target"]).sum()),
        "features": fsum,
        "tier_b": tb_meta,
        # proteome-level QC (metadata only; never a feature column)
        "qc_flags": tb_meta["qc"] if tb_meta else None,
        "tier_c": tc_meta,
        "timings": {**timings, "total_s": time.time() - t0},
        "columns": list(out.columns),
    }
    if output is not None:
        write(out, output, meta)
    return out, meta


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write(df: pd.DataFrame, output: str | Path, meta: dict) -> None:
    output = Path(output)
    meta = {**meta}
    meta.setdefault("columns_hash", registry.table_columns_hash(df.columns))  # definitions of this table's own columns
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix == ".parquet":
        import pyarrow as pa
        import pyarrow.parquet as pq

        table = pa.Table.from_pandas(df, preserve_index=False)
        md = dict(table.schema.metadata or {})
        md[b"slparalog"] = json.dumps(meta, default=str).encode()
        pq.write_table(table.replace_schema_metadata(md), output)
    else:
        df.to_csv(output, sep="\t", index=False, float_format="%.6g")
    Path(str(output) + ".meta.json").write_text(json.dumps(meta, indent=1, default=str))
