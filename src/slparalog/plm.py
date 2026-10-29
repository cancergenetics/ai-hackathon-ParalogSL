"""Tier C: mean-pooled per-protein embeddings from a protein language model (default ESM-2 150M).

* Model id and pinned HF revision come from ``config/plm.yaml``; weights live under ``models_cache/``.
* Proteins longer than the window are split into ``ceil(L/window)`` equal-length chunks, embedded
  independently and combined by a length-weighted average of the chunk means (no residue dropped).
* Chunks are batched by length (token budget). fp16 on MPS, fp32 on CPU (fallback).
* Results are cached per *unique sequence* (key = SHA-256 of the canonical sequence) as float16 numpy
  shards in a user-given directory, so identical sequences are computed once and reruns are free.
* ``embed_fasta`` returns an (n_records x d) matrix aligned to the FASTA's record order.

PROCESS ISOLATION: torch is imported lazily and only when something has to be computed. Never import
this module's compute path in a process that fits LightGBM (libomp conflict); use
``embed_fasta_subprocess`` / ``python -m slparalog.plm`` / ``slparalog embed`` instead. The cache
readers (``read_fasta``, ``load_cached_matrix``) are torch-free.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from slparalog.config import CONFIG_DIR, PROJECT_ROOT, load_yaml

CACHE_META = "meta.json"
CACHE_FORMAT = 1
DEFAULT_LOCK = PROJECT_ROOT / "logs" / ".plm_embed.lock"


# --------------------------------------------------------------------------- config / env
def load_plm_config(path: str | Path | None = None) -> dict:
    return load_yaml(path or CONFIG_DIR / "plm.yaml")


def model_spec(model_key: str | None = None, config: dict | None = None) -> dict:
    cfg = config or load_plm_config()
    key = model_key or cfg["default_model"]
    if key not in cfg["models"]:
        raise KeyError(f"unknown model {key!r}; known: {sorted(cfg['models'])}")
    return {"key": key, **cfg["models"][key]}


def setup_env(config: dict | None = None) -> None:
    """Point HF/torch caches at models_cache/ (must run before transformers is imported)."""
    cfg = config or load_plm_config()
    os.environ.setdefault("HF_HOME", str(PROJECT_ROOT / cfg["hf_home"]))
    os.environ.setdefault("TORCH_HOME", str(PROJECT_ROOT / cfg["torch_home"]))
    os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")


# --------------------------------------------------------------------------- sequences
def read_fasta(path: str | Path) -> tuple[list[str], list[str]]:
    """Return (ids, sequences) in record order; id = first whitespace-delimited header token."""
    ids: list[str] = []
    seqs: list[str] = []
    buf: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if ids:
                    seqs.append("".join(buf))
                ids.append(line[1:].split()[0] if line[1:].split() else "")
                buf = []
            elif line.strip():
                buf.append(line.strip())
    if ids:
        seqs.append("".join(buf))
    return ids, seqs


def canonical(seq: str) -> str:
    """Upper-case, whitespace and trailing/internal '*' removed."""
    return "".join(seq.split()).upper().replace("*", "")


def seq_sha256(seq: str) -> str:
    return hashlib.sha256(canonical(seq).encode("ascii", "replace")).hexdigest()


def chunk_sequence(seq: str, window: int) -> list[str]:
    """Split into ceil(L/window) equal-length chunks (last may be shorter by < n_chunks residues)."""
    n = max(1, math.ceil(len(seq) / window))
    size = math.ceil(len(seq) / n)
    return [seq[i : i + size] for i in range(0, len(seq), size)]


def make_batches(items: list[tuple[int, str]], tok_budget: int, max_bs: int) -> list[list[tuple[int, str]]]:
    """Sort (idx, seq) by length; greedy batches under a padded-token budget (+2 special tokens)."""
    items = sorted(items, key=lambda x: len(x[1]))
    batches: list[list[tuple[int, str]]] = []
    cur: list[tuple[int, str]] = []
    curmax = 0
    for it in items:
        ln = len(it[1]) + 2
        if cur and (max(curmax, ln) * (len(cur) + 1) > tok_budget or len(cur) >= max_bs):
            batches.append(cur)
            cur, curmax = [], 0
        cur.append(it)
        curmax = max(curmax, ln)
    if cur:
        batches.append(cur)
    return batches


# --------------------------------------------------------------------------- cache
class EmbeddingCache:
    """Float16 shards ``shard_XXXXXX.npy`` (n x d) + ``shard_XXXXXX.keys`` (one SHA-256 per row).

    A shard is valid only when its ``.keys`` file exists (written last, atomically), so an
    interrupted run never leaves a half-written shard that is trusted.
    """

    def __init__(self, directory: str | Path, meta: dict):
        self.dir = Path(directory)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.meta = {"format": CACHE_FORMAT, **meta}
        mp = self.dir / CACHE_META
        if mp.exists():
            old = json.loads(mp.read_text())
            if old != self.meta:
                raise ValueError(
                    f"cache {self.dir} was built with different settings:\n  cache: {old}\n  now:   {self.meta}\n"
                    "use a separate cache directory per model/revision/window"
                )
        else:
            mp.write_text(json.dumps(self.meta, indent=1, sort_keys=True))
        self._shards: list[np.ndarray] = []
        self._index: dict[str, tuple[int, int]] = {}
        self._next = 0
        self._load()

    def _load(self) -> None:
        for kp in sorted(self.dir.glob("shard_*.keys")):
            npy = kp.with_suffix(".npy")
            if not npy.exists():
                continue
            keys = kp.read_text().split()
            arr = np.load(npy, mmap_mode="r")
            if arr.shape[0] != len(keys) or arr.dtype != np.float16:
                raise ValueError(f"corrupt cache shard {npy}")
            si = len(self._shards)
            self._shards.append(arr)
            for ri, k in enumerate(keys):
                self._index.setdefault(k, (si, ri))
            self._next = max(self._next, int(kp.stem.split("_")[1]) + 1)

    def __contains__(self, key: str) -> bool:
        return key in self._index

    def __len__(self) -> int:
        return len(self._index)

    def get(self, key: str) -> np.ndarray:
        si, ri = self._index[key]
        return np.asarray(self._shards[si][ri])

    def add_shard(self, keys: list[str], emb: np.ndarray) -> None:
        if not keys:
            return
        emb16 = np.ascontiguousarray(emb.astype(np.float16))
        name = f"shard_{self._next:06d}"
        self._next += 1
        tmp_npy = self.dir / f"{name}.tmp.npy"
        np.save(tmp_npy, emb16)
        os.replace(tmp_npy, self.dir / f"{name}.npy")
        tmp_keys = self.dir / f"{name}.keys.tmp"
        tmp_keys.write_text("\n".join(keys) + "\n")
        os.replace(tmp_keys, self.dir / f"{name}.keys")
        si = len(self._shards)
        self._shards.append(np.load(self.dir / f"{name}.npy", mmap_mode="r"))
        for ri, k in enumerate(keys):
            self._index.setdefault(k, (si, ri))


def _cache_meta(spec: dict, pooling: str) -> dict:
    return {
        "hf_id": spec["hf_id"],
        "revision": spec["revision"],
        "window": int(spec["window"]),
        "dim": int(spec["dim"]),
        "pooling": pooling,
    }


def load_cached_matrix(seqs: list[str], cache_dir: str | Path, model_key: str | None = None) -> np.ndarray:
    """torch-free: (n x d) float32 matrix from the cache only; KeyError if any sequence is missing."""
    cfg = load_plm_config()
    spec = model_spec(model_key, cfg)
    cache = EmbeddingCache(cache_dir, _cache_meta(spec, cfg["pooling"]))
    out = np.empty((len(seqs), spec["dim"]), np.float32)
    for i, s in enumerate(seqs):
        out[i] = cache.get(seq_sha256(s))
    return out


# --------------------------------------------------------------------------- model
class Embedder:
    """Loads the pinned model once; ``embed_chunks`` returns mean-pooled float32 vectors per chunk."""

    def __init__(self, spec: dict, device: str | None = None, dtype: str | None = None, config: dict | None = None):
        cfg = config or load_plm_config()
        setup_env(cfg)
        import torch  # lazy: keep torch out of processes that never compute
        from transformers import AutoModel, AutoTokenizer

        self.torch = torch
        self.spec = spec
        if device is None:
            device = "mps" if "mps" in cfg["device_preference"] and torch.backends.mps.is_available() else "cpu"
        self.device = device
        self.dtype_name = dtype or cfg["dtype"][device]
        dt = {"fp32": torch.float32, "fp16": torch.float16, "bf16": torch.bfloat16}[self.dtype_name]
        kw = dict(revision=spec["revision"])
        try:
            self.tok = AutoTokenizer.from_pretrained(spec["hf_id"], local_files_only=True, **kw)
            model = AutoModel.from_pretrained(spec["hf_id"], dtype=dt, local_files_only=True, **kw)
        except OSError:  # not cached yet -> fetch the pinned revision
            self.tok = AutoTokenizer.from_pretrained(spec["hf_id"], **kw)
            model = AutoModel.from_pretrained(spec["hf_id"], dtype=dt, **kw)
        try:
            self.model = model.to(device).eval()
        except Exception as e:  # MPS unavailable at runtime -> CPU fallback
            if device == "cpu":
                raise
            print(f"[plm] {device} failed ({e!r}); falling back to cpu fp32", file=sys.stderr)
            self.device, self.dtype_name = "cpu", "fp32"
            self.model = model.float().to("cpu").eval()
        assert self.model.config.hidden_size == spec["dim"], "config dim does not match model"

    def _sync(self) -> None:
        if self.device == "mps":
            self.torch.mps.synchronize()

    def embed_batch(self, seqs: list[str]) -> np.ndarray:
        torch = self.torch
        with torch.inference_mode():
            enc = self.tok(seqs, return_tensors="pt", padding=True)
            enc = {k: v.to(self.device) for k, v in enc.items()}
            h = self.model(**enc).last_hidden_state.float()
            attn = enc["attention_mask"]
            mask = attn.clone()
            mask[:, 0] = 0  # <cls>
            last = attn.sum(1) - 1
            mask[torch.arange(mask.shape[0], device=mask.device), last] = 0  # <eos>
            m = mask.unsqueeze(-1).float()
            emb = (h * m).sum(1) / m.sum(1).clamp(min=1)
            out = emb.cpu().numpy()
        return out

    def empty_cache(self) -> None:
        if self.device == "mps":
            self.torch.mps.empty_cache()


# --------------------------------------------------------------------------- core API
def embed_sequences(
    seqs: list[str],
    cache_dir: str | Path,
    model_key: str | None = None,
    *,
    device: str | None = None,
    window: int | None = None,
    log=None,
    config: dict | None = None,
) -> tuple[np.ndarray, dict]:
    """Embed sequences (any order, duplicates allowed). Returns ((n x d) float32, stats).

    Values are the float16-cached embeddings upcast to float32, so a rerun returns identical bytes.
    ``window`` overrides the config window (tests); it is recorded in the cache meta.
    """
    cfg = config or load_plm_config()
    spec = model_spec(model_key, cfg)
    if window is not None:
        spec["window"] = int(window)
    log = log or (lambda msg: print(msg, file=sys.stderr, flush=True))
    t_start = time.time()
    cache = EmbeddingCache(cache_dir, _cache_meta(spec, cfg["pooling"]))

    canon = [canonical(s) for s in seqs]
    if any(len(s) == 0 for s in canon):
        raise ValueError("empty sequence in input")
    keys = [hashlib.sha256(s.encode("ascii", "replace")).hexdigest() for s in canon]
    uniq: dict[str, str] = {}
    for k, s in zip(keys, canon):
        uniq.setdefault(k, s)
    todo = [k for k in uniq if k not in cache]
    stats = {
        "n_input": len(seqs),
        "n_unique": len(uniq),
        "n_cache_hits": len(uniq) - len(todo),
        "n_computed": len(todo),
        "model": spec["key"],
        "hf_id": spec["hf_id"],
        "revision": spec["revision"],
        "window": spec["window"],
    }
    log(f"[plm] {len(seqs)} input, {len(uniq)} unique, {len(todo)} to compute ({stats['n_cache_hits']} cached)")

    embed_s = 0.0
    if todo:
        t_load = time.time()
        emb = Embedder(spec, device=device, config=cfg)
        stats.update(device=emb.device, dtype=emb.dtype_name, load_s=round(time.time() - t_load, 2))
        items: list[tuple[int, str]] = []
        owner: list[int] = []
        remaining: list[int] = []
        for ti, k in enumerate(todo):
            chunks = chunk_sequence(uniq[k], spec["window"])
            remaining.append(len(chunks))
            for c in chunks:
                items.append((len(owner), c))
                owner.append(ti)
        batches = make_batches(items, int(spec["tok_budget"]), int(spec["max_batch"]))
        stats.update(n_chunks=len(items), n_residues=sum(len(c) for _, c in items), n_batches=len(batches))
        acc = np.zeros((len(todo), spec["dim"]), np.float64)
        wsum = np.zeros(len(todo), np.float64)
        pending: list[int] = []
        flush_n = int(cfg.get("cache_shard_size", 2000))
        # warm-up (kernel compilation / allocation), not timed against throughput
        emb.embed_batch([items[0][1][:8]])
        emb._sync()
        t0 = time.time()
        last_log = t0
        done = 0

        def flush():
            if pending:
                cache.add_shard([todo[i] for i in pending], (acc[pending] / wsum[pending, None]))
                pending.clear()

        for bi, b in enumerate(batches):
            out = emb.embed_batch([s for _, s in b])
            for (gi, s), e in zip(b, out):
                o = owner[gi]
                acc[o] += e.astype(np.float64) * len(s)
                wsum[o] += len(s)
                remaining[o] -= 1
                if remaining[o] == 0:
                    pending.append(o)
                    done += 1
            if len(pending) >= flush_n:
                flush()
            if bi % 200 == 199:
                emb.empty_cache()
            if time.time() - last_log > 60:
                el = time.time() - t0
                log(f"[plm] {done}/{len(todo)} seqs, batch {bi + 1}/{len(batches)}, {el:.0f}s, {done / el:.2f} seq/s")
                last_log = time.time()
        emb._sync()
        flush()
        embed_s = time.time() - t0
        del emb
    stats["embed_s"] = round(embed_s, 2)
    stats["seq_per_s_computed"] = round(len(todo) / embed_s, 3) if todo else None

    out = np.empty((len(seqs), spec["dim"]), np.float32)
    for i, k in enumerate(keys):
        out[i] = cache.get(k)
    stats["total_s"] = round(time.time() - t_start, 2)
    return out, stats


def embed_fasta(
    fasta: str | Path,
    cache_dir: str | Path,
    out: str | Path | None = None,
    model_key: str | None = None,
    **kw,
) -> tuple[np.ndarray, list[str], dict]:
    """(n_records x d) float32 matrix aligned to FASTA record order, the ids, and stats.

    If ``out`` is given, writes ``out`` (.npy) and ``<out stem>.ids.txt`` (one id per line).
    """
    ids, seqs = read_fasta(fasta)
    if not ids:
        raise ValueError(f"no records in {fasta}")
    mat, stats = embed_sequences(seqs, cache_dir, model_key, **kw)
    stats["fasta"] = str(fasta)
    stats["n_records"] = len(ids)
    if out is not None:
        out = Path(out)
        out.parent.mkdir(parents=True, exist_ok=True)
        np.save(out, mat)
        out.with_suffix(".ids.txt").write_text("\n".join(ids) + "\n")
    return mat, ids, stats


def embed_fasta_subprocess(
    fasta: str | Path, cache_dir: str | Path, out: str | Path, model_key: str | None = None
) -> np.ndarray:
    """Run the embedding in a fresh process (safe from a process that uses LightGBM); load result."""
    cmd = [sys.executable, "-m", "slparalog.plm", str(fasta), "--cache", str(cache_dir), "--out", str(out)]
    if model_key:
        cmd += ["--model", model_key]
    subprocess.run(cmd, check=True)
    return np.load(out)


# --------------------------------------------------------------------------- CLI
def lowpower_state() -> str:
    try:
        r = subprocess.run(["pmset", "-g"], capture_output=True, text=True, timeout=10).stdout
        return " | ".join(ln.strip() for ln in r.splitlines() if "lowpowermode" in ln.lower()) or "unknown"
    except Exception:
        return "unknown"


def peak_rss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round((r if sys.platform == "darwin" else r * 1024) / 1e9, 3)


def build_parser(parser: argparse.ArgumentParser | None = None) -> argparse.ArgumentParser:
    p = parser or argparse.ArgumentParser(prog="python -m slparalog.plm", description=__doc__.split("\n")[0])
    p.add_argument("fasta", help="proteome FASTA (record order defines the output rows)")
    p.add_argument("--cache", required=True, help="cache directory (per model; float16 shards)")
    p.add_argument("--out", required=True, help="output .npy; ids go to <stem>.ids.txt")
    p.add_argument("--model", default=None, help="key in config/plm.yaml (default: config default_model)")
    p.add_argument("--device", default=None, choices=["mps", "cpu"])
    p.add_argument("--bench-json", default=None, help="write wall time / seq/s / peak RSS / power state here")
    p.add_argument("--lock", default=str(DEFAULT_LOCK), help="lock file: only one embedding job at a time")
    return p


def run(args: argparse.Namespace) -> int:
    DEFAULT_LOCK.parent.mkdir(parents=True, exist_ok=True)
    lock_path = Path(args.lock)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lockfh = open(lock_path, "w")
    try:
        fcntl.flock(lockfh, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print(f"another embedding job holds {lock_path}; refusing to start", file=sys.stderr)
        return 3
    setup_env()
    t0 = time.time()
    _, ids, stats = embed_fasta(args.fasta, args.cache, args.out, args.model, device=args.device)
    wall = time.time() - t0
    stats.update(
        wall_s=round(wall, 2),
        seq_per_s_wall=round(len(ids) / wall, 3),
        peak_rss_gb=peak_rss_gb(),
        lowpowermode=lowpower_state(),
        out=str(args.out),
        cache=str(args.cache),
        cpu_count=os.cpu_count(),
        date=time.strftime("%Y-%m-%d %H:%M:%S"),
    )
    try:
        import torch
        import transformers

        stats.update(torch=torch.__version__, transformers=transformers.__version__)
    except ImportError:
        pass
    print(json.dumps(stats), flush=True)
    if args.bench_json:
        Path(args.bench_json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.bench_json).write_text(json.dumps(stats, indent=1, sort_keys=True) + "\n")
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
