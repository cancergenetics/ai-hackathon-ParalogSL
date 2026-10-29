"""``slparalog embed`` subcommand. Always runs the embedding in a fresh child process so that torch is
never imported into the (possibly LightGBM-using) CLI process."""

from __future__ import annotations

import argparse
import subprocess
import sys


def _run(args: argparse.Namespace) -> int:
    cmd = [sys.executable, "-m", "slparalog.plm", args.fasta, "--cache", args.cache, "--out", args.out]
    for flag, val in (("--model", args.model), ("--device", args.device), ("--bench-json", args.bench_json)):
        if val:
            cmd += [flag, val]
    return subprocess.run(cmd).returncode


def register(sub: argparse._SubParsersAction) -> None:
    from slparalog.plm import build_parser as _plm_parser  # torch-free import

    p = sub.add_parser("embed", help="mean-pooled PLM embeddings for a proteome FASTA (own process)")
    _plm_parser(p)
    p.set_defaults(func=_run)
