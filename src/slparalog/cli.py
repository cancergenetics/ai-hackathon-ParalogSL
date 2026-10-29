"""Console entry point: ``slparalog <subcommand>``.

``predict --model DIR`` runs the paralog caller and the tiers the model needs and writes the scored pairs
(``slparalog.predict``); ``predict --features-only`` writes the caller + feature table (no model). Without
a model (none given, none installed) ``predict`` keeps its legacy output: features plus NaN ``sl_score``/``sl_rank``
placeholders (pinned by tests/test_tier_a_toy.py::test_cli_smoke).
``setup`` fetches the panel / ESM-2 weights / model (``slparalog.setup_data``); ``train`` fits a model
artifact from feature tables and development labels (``slparalog.train``).
"""

from __future__ import annotations

import argparse
import json
import logging
import sys


def _model_available(args: argparse.Namespace) -> bool:
    """True if ``--model`` was given or a default model is installed under the data dir."""
    from . import setup_data

    if args.model is not None:
        return True
    try:
        setup_data.resolve_model(None, setup_data.data_dir(args.data_dir))
        return True
    except setup_data.SetupError:
        return False


def _predict(args: argparse.Namespace) -> int:
    from . import pipeline

    if not args.proteome or not args.output:
        print("predict needs a proteome FASTA and -o/--output", file=sys.stderr)
        return 1
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(asctime)s %(name)s %(message)s")
    if args.features_only and args.model is not None:
        print("--features-only cannot be combined with --model", file=sys.stderr)
        return 1
    if not args.features_only and _model_available(args):
        return _predict_scored(args)
    df, meta = pipeline.run(args.proteome, args.output, threads=args.threads,
                            workdir=args.workdir, keep_work=args.keep_work,
                            features_only=args.features_only, config_path=args.config,
                            tier_b=args.tier_b, panel_dmin=args.panel_dmin,
                            tier_c=args.tier_c, plm_cache=args.plm_cache, plm_model=args.plm_model)
    c = meta["caller"]
    print(f"proteins={meta['n_proteins']} pairs={c['n_pairs']} "
          f"capped_queries={meta['saturation']['n_capped_first_round']} "
          f"one_direction={c['frac_called_one_direction']:.4f} "
          f"time={meta['timings']['total_s']:.1f}s -> {args.output}")
    if not args.features_only:
        print("note: no model (--model, or one installed by `setup --model`): sl_score and sl_rank are NaN "
              "placeholders", file=sys.stderr)
    return 0


def _predict_scored(args: argparse.Namespace) -> int:
    from . import predict as predict_mod
    from .model import ModelError
    from .setup_data import SetupError

    if args.tier_b or args.tier_c or args.panel_dmin is not None:
        print("note: --tier-b/--tier-c/--panel-dmin are ignored with a model; the model decides which tiers are computed",
              file=sys.stderr)
    try:
        table, meta = predict_mod.run(
            args.proteome, args.output, args.model, threads=args.threads, data_dir=args.data_dir,
            workdir=args.workdir, keep_work=args.keep_work, config_path=args.config,
            panel_config=args.panel_config, plm_cache=args.plm_cache, plm_model=args.plm_model,
            allow_sklearn_mismatch=args.allow_sklearn_mismatch)
    except (ModelError, SetupError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    for w in meta["qc_warnings"]:
        print(f"QC warning: {w}", file=sys.stderr)
    print(f"proteins={meta['n_proteins']} pairs={len(table)} model={meta['model']['model_id']} "
          f"time={meta['timings']['total_s']:.1f}s -> {args.output}")
    return 0


def _setup(args: argparse.Namespace) -> int:
    from . import setup_data as sd

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(asctime)s %(name)s %(message)s")
    dd = sd.data_dir(args.data_dir)
    do_panel, do_plm = args.panel, args.plm
    if not (args.panel or args.plm or args.model):
        do_panel = do_plm = True  # default: everything that can be fetched
    try:
        if do_panel:
            r = sd.setup_panel(dd, threads=args.threads, keep_downloads=args.keep_downloads)
            print(f"panel: {r['panel_fasta']}, {r['n_species']} species, "
                  f"{len(r['downloaded'])} downloaded, {len(r['present'])} already present -> {dd / 'panel'}")
        if do_plm:
            r = sd.setup_plm(dd)
            print(f"plm: {r['model']} revision {r['revision'][:12]} {r['status']} -> {r['hf_home']}")
        if args.model:
            r = sd.setup_model(dd, args.model, allow_sklearn_mismatch=args.allow_sklearn_mismatch)
            print(f"model: {r['model_id']} {r['status']} -> {r['path']}")
    except sd.SetupError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    return 0


def _kv(items: list[str] | None, what: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for it in items or []:
        if "=" not in it:
            raise SystemExit(f"{what} expects SPECIES=PATH, got {it!r}")
        sp, path = it.split("=", 1)
        out.setdefault(sp, []).append(path)
    return out


def _train(args: argparse.Namespace) -> int:
    from . import train as train_mod
    from .model import ModelError

    feats = _kv(args.features, "--features")
    labs = _kv(args.labels, "--labels")
    if any(len(v) != 1 for v in labs.values()):
        print("--labels takes one table per species", file=sys.stderr)
        return 1
    params = {}
    for it in args.param or []:
        k, _, v = it.partition("=")
        params[k] = json.loads(v)
    try:
        meta = train_mod.train(feats, {k: v[0] for k, v in labs.items()}, args.feature_set, args.model_type,
                               args.output, tier=args.tier, params=params, seed=args.seed, overwrite=args.overwrite)
    except (PermissionError, ModelError, ValueError) as e:  # PermissionError: load_labels refuses held-out paths
        print(f"error: {e}", file=sys.stderr)
        return 1
    n = sum(v["n_pairs"] for v in meta["training"]["per_species"].values())
    print(f"model {meta['model_id']}: {len(meta['features'])} features, {n} training pairs -> {args.output}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="slparalog",
        description="Predict synthetic-lethal paralog pairs from a proteome.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_predict = sub.add_parser("predict", help="call paralog pairs in a proteome FASTA and score them with a model")
    p_predict.add_argument("proteome", nargs="?", help="proteome FASTA (one protein per gene)")
    p_predict.add_argument("-o", "--output", help="output pairs table (.tsv; .parquet with --features-only or for a parquet score table)")
    p_predict.add_argument("--model", help="model directory (model.json + model.joblib) or a model id installed under "
                           "the data dir [the model last installed by `setup --model`]")
    p_predict.add_argument("--data-dir", help="data directory of `slparalog setup` [$SLPARALOG_DATA or ~/.slparalog]")
    p_predict.add_argument("--allow-sklearn-mismatch", action="store_true",
                           help="load a model fitted with another scikit-learn minor version (unsafe)")
    p_predict.add_argument("--panel-config", default=None, help=argparse.SUPPRESS)  # panel.yaml (tests / dev checkout)
    p_predict.add_argument("-t", "--threads", type=int, default=1, help="CPU threads [1]")
    p_predict.add_argument("--features-only", action="store_true",
                           help="write caller + features only (no model, no score columns)")
    p_predict.add_argument("--workdir", help="working directory for MMseqs2 files (kept)")
    p_predict.add_argument("--keep-work", action="store_true",
                           help="keep the temporary working directory")
    p_predict.add_argument("--config", help=argparse.SUPPRESS)
    p_predict.add_argument("--tier-b", action="store_true",
                           help="--features-only: add Tier B reference-panel features (per-species panel search)")
    # hidden: training-time panel mask, removes panel species with d_s < c (algorithms_comparative change 7)
    p_predict.add_argument("--panel-dmin", type=float, default=None, help=argparse.SUPPRESS)
    p_predict.add_argument("--tier-c", action="store_true",
                           help="--features-only: add Tier C PLM embedding-pair features (embeds the proteome with ESM-2 in a child process)")
    p_predict.add_argument("--plm-cache", default=None,
                           help="embedding cache directory for Tier C (reused across runs) [temporary]")
    p_predict.add_argument("--plm-model", default=None, help=argparse.SUPPRESS)  # config/plm.yaml key (tests)
    p_predict.add_argument("-v", "--verbose", action="store_true")
    p_predict.set_defaults(func=_predict)

    p_setup = sub.add_parser("setup", help="fetch and build the species panel, ESM-2 weights and a model "
                             "(no flag: --panel --plm)")
    p_setup.add_argument("--panel", action="store_true",
                         help="download the panel proteomes (checksummed) and build the MMseqs2 databases")
    p_setup.add_argument("--plm", action="store_true", help="fetch the ESM-2 weights at the pinned revision")
    p_setup.add_argument("--model", metavar="DIR", help="install a model artifact from a local directory")
    p_setup.add_argument("--data-dir", help="data directory [$SLPARALOG_DATA or ~/.slparalog]")
    p_setup.add_argument("-t", "--threads", type=int, default=1, help="CPU threads for the database build [1]")
    p_setup.add_argument("--keep-downloads", action="store_true", help="keep the downloaded proteome archives")
    p_setup.add_argument("--allow-sklearn-mismatch", action="store_true", help=argparse.SUPPRESS)
    p_setup.add_argument("-v", "--verbose", action="store_true")
    p_setup.set_defaults(func=_setup)

    p_train = sub.add_parser("train", help="fit a model artifact from feature tables and development labels")
    p_train.add_argument("--features", action="append", metavar="SPECIES=PARQUET", required=True,
                         help="feature table written by `predict --features-only` for one species (repeat; "
                              "several tables of a species are joined by gene pair)")
    p_train.add_argument("--labels", action="append", metavar="SPECIES=TSV", required=True,
                         help="development label table of one species (repeat); held-out label paths are refused")
    p_train.add_argument("--feature-set", required=True,
                         help="registry tiers / Tier A groups joined by '+', e.g. A, A+B, A+B+C "
                              "(groups: config/models.yaml feature_groups)")
    p_train.add_argument("--model-type", choices=["gbdt", "lr"], required=True)
    p_train.add_argument("-o", "--output", required=True, help="model directory to create")
    p_train.add_argument("--tier", default="primary", help="label tier [primary]")
    p_train.add_argument("--param", action="append", metavar="KEY=JSON",
                         help="hyper-parameter value from the declared grid of config/models.yaml "
                              "[the configuration recorded there]")
    p_train.add_argument("--seed", type=int, default=None, help="a seed of config/seeds.yaml [the first]")
    p_train.add_argument("--overwrite", action="store_true")
    p_train.set_defaults(func=_train)
    __import__("slparalog.cli_embed", fromlist=["register"]).register(sub)  # `embed` subcommand (PLM, own process)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    code = args.func(args)
    sys.exit(code)


if __name__ == "__main__":
    main()
