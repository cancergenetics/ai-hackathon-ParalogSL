"""``slparalog predict``: proteome FASTA + model artifact -> every called paralog pair with an SL score.

The tiers the model needs (``model.json: required_tiers``) are computed by ``pipeline.run`` (Tier A always,
Tier B / Tier C only if the model uses a feature of that tier). Process isolation as in the pipeline: the
PLM embedding runs in a child process; this process uses numpy/sklearn only and never imports torch.

Output (TSV; ``pandas.read_csv(path, sep="\\t", comment="#")``): a ``#`` header block (tool version,
model id, input checksum, run time, QC warnings), then one row per called pair sorted by ``sl_rank``:
gene_a, gene_b, sl_score, sl_rank (1 = highest score; ties share the best rank), ident_global and the QC
flag columns of the tiers that were computed (registry ``kind == "qc"``). ``sl_score`` is a ranking score,
not a calibrated probability.
"""

from __future__ import annotations

import contextlib
import datetime as dt
import json
import os
from pathlib import Path

import pandas as pd
from scipy.stats import rankdata

from . import __version__, pipeline, setup_data
from .features import registry
from .model import ModelError, SLModel, load_model

CORE_COLUMNS = ["gene_a", "gene_b", "sl_score", "sl_rank", "ident_global"]


def qc_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c in registry.BY_NAME and registry.BY_NAME[c].kind == "qc"]


def qc_warnings(meta: dict, n_pairs: int) -> list[str]:
    """Proteome-level QC warnings from the pipeline / panel metadata (never model inputs)."""
    w: list[str] = []
    if n_pairs == 0:
        w.append("no paralog pairs were called in this proteome")
    qf = meta.get("qc_flags")
    if qf:
        w += [f"Tier B: {x}" for x in qf.get("warnings", [])]
    if meta.get("tool_matches_config") is False:
        w.append("MMseqs2 binary differs from the pinned version/checksum in config/tools.yaml")
    inp = meta.get("input", {})
    if inp.get("n_dropped_empty"):
        w.append(f"{inp['n_dropped_empty']} empty protein records were dropped from the input")
    if inp.get("n_internal_stops_to_x") or inp.get("n_nonstandard_to_x"):
        w.append(f"input residues replaced by X: {inp.get('n_internal_stops_to_x', 0)} internal stops, "
                 f"{inp.get('n_nonstandard_to_x', 0)} non-standard letters")
    return w


def build_table(feats: pd.DataFrame, model: SLModel) -> pd.DataFrame:
    """Score table from a pipeline feature table: core columns + QC flags, sorted by ``sl_rank``."""
    score = model.score(feats)
    if len(score) and not pd.Series(score).notna().all():
        raise ModelError("model produced NaN scores")
    out = pd.DataFrame({"gene_a": feats["gene_a"].to_numpy(), "gene_b": feats["gene_b"].to_numpy(),
                        "sl_score": score})
    out["sl_rank"] = rankdata(-score, method="min").astype("int64") if len(score) else pd.Series([], dtype="int64")
    out["ident_global"] = feats["ident_global"].to_numpy()
    for c in qc_columns(feats):
        out[c] = feats[c].to_numpy()
    return out.sort_values(["sl_rank", "gene_a", "gene_b"], kind="mergesort").reset_index(drop=True)


def header_lines(meta: dict, model: SLModel, n_pairs: int, warnings: list[str], started: dt.datetime,
                 tiers: list[str]) -> list[str]:
    m = model.meta
    lines = [
        f"slparalog {__version__}",
        f"model_id: {model.model_id}",
        f"model: type={m['model_type']} feature_set={m['feature_set']} n_features={len(model.features)} "
        f"registry_hash={m['registry_hash'][:16]} sklearn={m['sklearn_version']} git={m.get('git', {}).get('commit')}",
        f"input: {meta['input']['path']} sha256={meta['input']['sha256']}",
        f"n_proteins: {meta['n_proteins']}  n_pairs: {n_pairs}  tiers_computed: {','.join(tiers)}",
        f"run_utc: {started.strftime('%Y-%m-%dT%H:%M:%SZ')}  run_time_s: {meta['timings']['total_s']:.1f}",
        "score: sl_score is a ranking score in [0, 1], not a calibrated probability; sl_rank 1 = highest score",
    ]
    if warnings:
        lines += [f"qc_warning: {w}" for w in warnings]
    else:
        lines.append("qc_warnings: none")
    return lines


def write_tsv(table: pd.DataFrame, path: str | Path, header: list[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for ln in header:
            fh.write(f"# {ln}\n")
        table.to_csv(fh, sep="\t", index=False, float_format="%.6g")


@contextlib.contextmanager
def _env(updates: dict[str, str]):
    old = {k: os.environ.get(k) for k in updates}
    os.environ.update(updates)
    try:
        yield
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def run(proteome: str | Path, output: str | Path, model: str | Path | None = None, *, threads: int = 1,
        data_dir: str | Path | None = None, workdir: str | Path | None = None, keep_work: bool = False,
        config_path: str | Path | None = None, panel_config: str | Path | None = None,
        plm_cache: str | Path | None = None, plm_model: str | None = None,
        allow_sklearn_mismatch: bool = False) -> tuple[pd.DataFrame, dict]:
    started = dt.datetime.now(dt.timezone.utc)
    dd = setup_data.data_dir(data_dir)
    mdir = setup_data.resolve_model(model, dd)
    sl = load_model(mdir, allow_sklearn_mismatch=allow_sklearn_mismatch)
    tiers = sl.tiers
    want_b, want_c = "B" in tiers, "C" in tiers
    pcfg = None
    if want_b:
        pcfg = Path(panel_config) if panel_config else setup_data.panel_config_path(dd)
        if not pcfg.exists():
            raise setup_data.SetupError(f"model {sl.model_id} needs the Tier B reference panel, which is not set up "
                                        f"under {dd}: run `slparalog setup --panel` first")
    env: dict[str, str] = {}
    if want_c:
        hf = setup_data.plm_home(dd)
        if hf is not None:
            env["HF_HOME"] = str(hf)
    with _env(env):
        feats, meta = pipeline.run(proteome, None, threads=threads, workdir=workdir, keep_work=keep_work,
                                   features_only=True, config_path=config_path, tier_b=want_b,
                                   panel_config_path=pcfg, tier_c=want_c, plm_cache=plm_cache, plm_model=plm_model)
    table = build_table(feats, sl)
    warnings = qc_warnings(meta, len(table))
    meta = {**meta, "features_only": False, "model": {k: sl.meta[k] for k in (
        "model_id", "model_type", "feature_set", "features", "required_tiers", "registry_hash", "model_sha256",
        "sklearn_version", "git", "seeds")}, "model_path": str(sl.path), "qc_warnings": warnings,
        "columns": list(table.columns), "columns_hash": registry.table_columns_hash(table.columns)}
    meta["model"]["columns_hash"] = sl.meta.get("columns_hash")
    out = Path(output)
    if out.suffix == ".parquet":
        pipeline.write(table, out, meta)
    else:
        write_tsv(table, out, header_lines(meta, sl, len(table), warnings, started, tiers))
        Path(str(out) + ".meta.json").write_text(json.dumps(meta, indent=1, default=str))
    return table, meta
