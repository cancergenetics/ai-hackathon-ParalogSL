"""Model artifact (save / load / validation), ``slparalog train`` and ``slparalog predict`` on the toy proteome.

Models here are trained on synthetic labels (tests/helpers_model.py): they test the plumbing (artifact
round trip, validation failures, tier selection, output format, process isolation), not predictive quality.
"""

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

import helpers_model as hm
from slparalog import cli, model, pipeline, predict, setup_data, train
from slparalog.config import CONFIG_DIR
from slparalog.features import registry
from slparalog.io import labels as labels_mod

TOY = Path(__file__).parent / "fixtures" / "toy_proteome.fasta"


@pytest.fixture(scope="module")
def env(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("model_predict")
    feat_path = tmp / "toy_features.parquet"
    feats, _ = pipeline.run(TOY, feat_path, threads=2, features_only=True)
    labs = hm.write_synthetic_labels(feats, tmp)
    return {"tmp": tmp, "feat": feat_path, "feats": feats, "labels": labs}


def _train(env, out, model_type="lr", feature_set="A", **kw):
    return train.train({"human": [env["feat"]], "yeast": [env["feat"]]}, env["labels"], feature_set, model_type,
                       out, **kw)


@pytest.fixture(scope="module")
def lr_dir(env):
    d = env["tmp"] / "lr_model"
    _train(env, d)
    return d


# ---------------------------------------------------------------------------------------------- artifact

def test_train_records_the_dev_selected_b1_variant(env, lr_dir, tmp_path):
    from slparalog.config import load_yaml
    from slparalog.eval.baselines import B1_VARIANTS

    meta = json.loads((lr_dir / "model.json").read_text())
    sel = load_yaml(CONFIG_DIR / "models.yaml")["b1"]["selected"]
    assert meta["b1_variant"] == sel == "ident_global" and sel in B1_VARIANTS
    assert model.load_model(lr_dir).meta["b1_variant"] == "ident_global"
    with pytest.raises(ValueError, match="b1_variant"):
        _train(env, tmp_path / "bad", b1_variant="not_a_variant")
    _train(env, tmp_path / "other", b1_variant="ident_local_max")
    assert json.loads((tmp_path / "other" / "model.json").read_text())["b1_variant"] == "ident_local_max"


def test_artifact_round_trip(env, lr_dir):
    meta = json.loads((lr_dir / "model.json").read_text())
    for key in ("features", "feature_set", "registry_hash", "required_tiers", "sklearn_version", "training",
                "git", "seeds", "model_sha256", "params", "caller_config_sha256", "model_id", "model_type"):
        assert key in meta, key
    assert meta["registry_hash"] == registry.registry_hash()
    assert meta["required_tiers"] == ["A"] and meta["feature_set"] == "A" and meta["model_type"] == "lr"
    assert set(meta["features"]) == {f for f in registry.MODEL_FEATURES if registry.BY_NAME[f].tier == "A"}
    assert meta["seeds"]["seed_used"] in meta["seeds"]["seed_list"]
    tr = meta["training"]
    assert set(tr["label_files"]) == {"human", "yeast"}
    for sp in ("human", "yeast"):
        assert tr["label_files"][sp]["sha256"] == hashlib.sha256(env["labels"][sp].read_bytes()).hexdigest()
        assert tr["feature_files"][sp][env["feat"].stem]["sha256"] == hashlib.sha256(env["feat"].read_bytes()).hexdigest()
    assert (lr_dir / "model.joblib").exists()

    m = model.load_model(lr_dir)
    s = m.score(env["feats"])
    assert s.shape == (len(env["feats"]),) and np.all((s >= 0) & (s <= 1)) and len(np.unique(s)) > 1
    # same numbers as the estimator saved (decision function -> logistic)
    import joblib
    est = joblib.load(lr_dir / "model.joblib")
    d = est.predict(env["feats"][m.features].to_numpy(float))
    np.testing.assert_allclose(s, 1 / (1 + np.exp(-d)))


def test_gbdt_artifact_and_recorded_hyperparameters(env):
    d = env["tmp"] / "gbdt_model"
    meta = _train(env, d, model_type="gbdt")
    mcfg = train.load_models_config()
    assert meta["params"] == mcfg["anchoring"]["ensemble"]["gbdt_params"]
    m = model.load_model(d)
    s = m.score(env["feats"])
    assert s.shape == (len(env["feats"]),) and np.isfinite(s).all() and np.all((s >= 0) & (s <= 1))


def test_hyperparameter_overrides_must_come_from_the_declared_grid(env):
    with pytest.raises(ValueError, match="outside the declared grid"):
        train.resolve_params("lr", train.load_models_config(), {"C": 123.0})
    with pytest.raises(ValueError, match="not a hyper-parameter"):
        train.resolve_params("lr", train.load_models_config(), {"nonsense": 1})
    assert train.resolve_params("lr", train.load_models_config(), {"C": 0.1})["C"] == 0.1


def test_overwrite_guard(env, lr_dir):
    with pytest.raises(model.ModelError, match="already holds a model"):
        _train(env, lr_dir)


def _copy(lr_dir, tmp_path, name="copy"):
    d = tmp_path / name
    shutil.copytree(lr_dir, d)
    return d


def _edit_meta(d, **kv):
    p = d / "model.json"
    m = json.loads(p.read_text())
    m.update(kv)
    p.write_text(json.dumps(m))


def test_registry_hash_mismatch_fails_loudly(lr_dir, tmp_path):
    # legacy artifact (no columns_hash): the whole-registry hash decides
    d = _copy(lr_dir, tmp_path)
    _edit_meta(d, registry_hash="0" * 64, columns_hash=None)
    with pytest.raises(model.ModelError, match="registry hash mismatch"):
        model.load_model(d)


def test_registry_change_in_code_fails(lr_dir, tmp_path, monkeypatch):
    d = _copy(lr_dir, tmp_path)
    _edit_meta(d, columns_hash=None)  # legacy artifact (no columns_hash): the whole-registry hash decides
    monkeypatch.setattr(registry, "registry_hash", lambda features=None: "f" * 64)
    with pytest.raises(model.ModelError, match="registry hash mismatch"):
        model.load_model(d)


def test_feature_presence_checks(lr_dir, tmp_path):
    meta = json.loads((lr_dir / "model.json").read_text())
    d = _copy(lr_dir, tmp_path, "a")
    _edit_meta(d, features=meta["features"] + ["no_such_feature"])
    with pytest.raises(model.ModelError, match="absent from the registry"):
        model.load_model(d)
    d2 = _copy(lr_dir, tmp_path, "b")
    _edit_meta(d2, features=meta["features"] + ["fragment_suspect"])  # registered, but a QC flag, not a model feature
    with pytest.raises(model.ModelError, match="not registered model features"):
        model.load_model(d2)
    d3 = _copy(lr_dir, tmp_path, "c")
    _edit_meta(d3, features=meta["features"][:-1], required_tiers=["A"])  # estimator was fitted on all of them
    with pytest.raises(model.ModelError, match="feature columns differ"):
        model.load_model(d3)


def test_feature_table_missing_a_model_feature_fails(env, lr_dir):
    m = model.load_model(lr_dir)
    with pytest.raises(model.ModelError, match="lacks model features"):
        m.score(env["feats"].drop(columns=["ident_local_min"]))


def test_tampered_joblib_and_sklearn_mismatch(lr_dir, tmp_path):
    d = _copy(lr_dir, tmp_path)
    with open(d / "model.joblib", "ab") as fh:
        fh.write(b"x")
    with pytest.raises(model.ModelError, match="checksum"):
        model.load_model(d)
    d2 = _copy(lr_dir, tmp_path, "skl")
    _edit_meta(d2, sklearn_version="0.0.1")
    with pytest.raises(model.ModelError, match="scikit-learn"):
        model.load_model(d2)
    assert model.load_model(d2, allow_sklearn_mismatch=True).model_id


def test_caller_config_change_fails(lr_dir, tmp_path):
    d = _copy(lr_dir, tmp_path)
    _edit_meta(d, caller_config_sha256="0" * 64)
    with pytest.raises(model.ModelError, match="caller.yaml"):
        model.load_model(d)


# ---------------------------------------------------------------------------------------------- train

def test_train_refuses_held_out_labels(env, tmp_path):
    guard = vars(labels_mod)["SEAL" + "ED_DIR"]  # the label guard's protected directory
    bad = guard / "human_x.tsv"
    with pytest.raises(PermissionError):
        train.train({"human": [env["feat"]]}, {"human": bad}, "A", "lr", tmp_path / "m")
    assert not (tmp_path / "m").exists()
    code = cli._train(cli.build_parser().parse_args(
        ["train", "--features", f"human={env['feat']}", "--labels", f"human={bad}", "--feature-set", "A",
         "--model-type", "lr", "-o", str(tmp_path / "m2")]))
    assert code == 1


def test_train_rejects_stale_feature_table(env, tmp_path):
    """A feature table written under a different registry hash is refused."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    t = pq.read_table(env["feat"])
    md = json.loads(t.schema.metadata[b"slparalog"])
    md["registry_hash"] = "1" * 64
    md.pop("columns_hash", None)  # legacy table (no columns_hash): the whole-registry hash decides
    p = tmp_path / "stale.parquet"
    pq.write_table(t.replace_schema_metadata({b"slparalog": json.dumps(md).encode()}), p)
    with pytest.raises(model.ModelError, match="registry"):
        train.train({"human": [p]}, {"human": env["labels"]["human"]}, "A", "lr", tmp_path / "m")


def test_train_feature_set_requires_the_tier(env, tmp_path):
    with pytest.raises(ValueError, match="matches no available"):
        _train(env, tmp_path / "m", feature_set="A+B")


def test_weights_are_the_harness_weights(env, monkeypatch):
    """train() fits with eval.stats.species_balanced_weights over the primary population."""
    seen = {}
    real = train.species_balanced_weights

    def spy(species, y, target_prevalence=None):
        w = real(species, y, target_prevalence)
        seen["w"], seen["species"], seen["y"] = w, np.asarray(species), np.asarray(y)
        return w

    monkeypatch.setattr(train, "species_balanced_weights", spy)
    _train(env, env["tmp"] / "w_model")
    sp, y, w = seen["species"], seen["y"], seen["w"]
    tot = {s: w[sp == s].sum() for s in np.unique(sp)}
    assert tot["human"] == pytest.approx(tot["yeast"])
    prev = {s: w[(sp == s) & (y == 1)].sum() / tot[s] for s in tot}
    assert prev["human"] == pytest.approx(prev["yeast"])


# ---------------------------------------------------------------------------------------------- predict

def _read_out(path):
    text = Path(path).read_text().splitlines()
    head = [ln for ln in text if ln.startswith("# ")]
    df = pd.read_csv(path, sep="\t", comment="#")
    return head, df


def test_predict_end_to_end_cli(env, lr_dir, tmp_path, capsys):
    out = tmp_path / "pred.tsv"
    with pytest.raises(SystemExit) as e:
        cli.main(["predict", str(TOY), "-o", str(out), "--model", str(lr_dir), "-t", "2",
                  "--data-dir", str(tmp_path / "nodata")])
    assert e.value.code == 0
    head, df = _read_out(out)
    m = model.load_model(lr_dir)
    feats = env["feats"]
    assert list(df.columns[:5]) == ["gene_a", "gene_b", "sl_score", "sl_rank", "ident_global"]
    assert {"fragment_suspect", "near_identical", "n_dir_called"} <= set(df.columns)
    assert len(df) == len(feats) == 6
    # every called pair, scored exactly as the model scores the pipeline features
    key = lambda d: set(zip(d.gene_a, d.gene_b))  # noqa: E731
    assert key(df) == key(feats)
    exp_score = pd.Series(m.score(feats), index=pd.MultiIndex.from_frame(feats[["gene_a", "gene_b"]]))
    got = df.set_index(["gene_a", "gene_b"])
    np.testing.assert_allclose(got["sl_score"], exp_score.loc[got.index].to_numpy(), rtol=1e-5, atol=1e-6)
    exp_ident = feats.set_index(["gene_a", "gene_b"])["ident_global"]
    np.testing.assert_allclose(got["ident_global"], exp_ident.loc[got.index].to_numpy(), atol=1e-6)
    # rank: 1 = highest score; the table is sorted by rank
    assert df["sl_rank"].is_monotonic_increasing and df["sl_rank"].min() == 1
    assert df["sl_score"].is_monotonic_decreasing
    # header block
    txt = "\n".join(head)
    assert f"slparalog {predict.__version__}" in txt
    assert f"model_id: {m.model_id}" in txt
    assert hashlib.sha256(TOY.read_bytes()).hexdigest() in txt
    assert "run_utc:" in txt and "run_time_s:" in txt and "qc_warning" in txt
    assert json.loads(Path(str(out) + ".meta.json").read_text())["model"]["model_id"] == m.model_id
    assert "QC warning" not in capsys.readouterr().err


def test_predict_without_any_model_keeps_the_nan_placeholder_output(tmp_path, capsys):
    """No --model and none installed: legacy placeholder output (pinned by tests/test_tier_a_toy.py::test_cli_smoke)
    with a note on stderr; the library entry point ``predict.run`` refuses instead."""
    with pytest.raises(SystemExit) as e:
        cli.main(["predict", str(TOY), "-o", str(tmp_path / "o.tsv"), "--data-dir", str(tmp_path / "nodata"), "-t", "2"])
    assert e.value.code == 0
    assert pd.read_csv(tmp_path / "o.tsv", sep="\t")["sl_score"].isna().all()
    assert "NaN placeholders" in capsys.readouterr().err
    with pytest.raises(setup_data.SetupError, match="no model given"):
        predict.run(TOY, tmp_path / "o2.tsv", None, data_dir=tmp_path / "nodata")
    assert not (tmp_path / "o2.tsv").exists()


def test_predict_rejects_stale_model(lr_dir, tmp_path, capsys):
    d = _copy(lr_dir, tmp_path)
    _edit_meta(d, registry_hash="0" * 64, columns_hash=None)  # legacy artifact: whole-registry hash
    with pytest.raises(SystemExit) as e:
        cli.main(["predict", str(TOY), "-o", str(tmp_path / "o.tsv"), "--model", str(d)])
    assert e.value.code == 1 and "registry hash mismatch" in capsys.readouterr().err
    assert not (tmp_path / "o.tsv").exists()


def test_features_only_still_works(tmp_path):
    out = tmp_path / "f.parquet"
    with pytest.raises(SystemExit) as e:
        cli.main(["predict", str(TOY), "-o", str(out), "--features-only", "-t", "2"])
    assert e.value.code == 0
    df = pd.read_parquet(out)
    assert "sl_score" not in df.columns and "ident_global" in df.columns and len(df) == 6


def test_qc_warnings_from_meta():
    meta = {"qc_flags": {"warnings": ["core set |G*|=3 < 200"]}, "tool_matches_config": False,
            "input": {"n_dropped_empty": 2, "n_internal_stops_to_x": 1, "n_nonstandard_to_x": 0}}
    w = predict.qc_warnings(meta, 0)
    assert any("no paralog pairs" in x for x in w) and any(x.startswith("Tier B:") for x in w)
    assert any("MMseqs2" in x for x in w) and any("empty protein" in x for x in w) and any("replaced by X" in x for x in w)
    assert predict.qc_warnings({"input": {}}, 3) == []


def _fake_model(tmp_path, feature_set_cols):
    from slparalog.eval.models import LRModel

    rng = np.random.default_rng(0)
    X = rng.normal(size=(30, len(feature_set_cols)))
    y = np.arange(30) % 2
    est = LRModel(feature_set_cols, C=1.0, seed=0, log_cols=[]).fit(X, y)
    d = tmp_path / "fake"
    model.save_model(d, est, feature_set="A+B+C", model_type="lr", features=feature_set_cols, params={}, training={},
                     seeds={})
    return d


def test_predict_computes_exactly_the_tiers_the_model_needs(tmp_path, monkeypatch):
    cols = [f for f in registry.MODEL_FEATURES if registry.BY_NAME[f].tier in ("A", "B", "C")][::7]
    assert {registry.BY_NAME[c].tier for c in cols} == {"A", "B", "C"}
    d = _fake_model(tmp_path, cols)
    assert model.load_model(d).tiers == ["A", "B", "C"]
    seen = {}

    class Stop(Exception):
        pass

    def fake_run(*a, **kw):
        seen.update(kw)
        seen["hf_home"] = __import__("os").environ.get("HF_HOME")
        raise Stop

    monkeypatch.setattr(pipeline, "run", fake_run)
    dd = tmp_path / "data"
    # Tier B needs the set-up panel
    with pytest.raises(setup_data.SetupError, match="setup --panel"):
        predict.run(TOY, tmp_path / "o.tsv", d, data_dir=dd)
    (dd / "panel").mkdir(parents=True)
    (dd / "panel" / "panel.yaml").write_text("{}")
    (dd / "plm" / "hf").mkdir(parents=True)
    (dd / "plm" / "PLM_SETUP.json").write_text("{}")
    with pytest.raises(Stop):
        predict.run(TOY, tmp_path / "o.tsv", d, data_dir=dd)
    assert seen["tier_b"] is True and seen["tier_c"] is True and seen["features_only"] is True
    assert seen["panel_config_path"] == dd / "panel" / "panel.yaml"
    assert seen["hf_home"] == str(dd / "plm" / "hf")
    assert "HF_HOME" not in __import__("os").environ or __import__("os").environ["HF_HOME"] != str(dd / "plm" / "hf")
    # a Tier A model computes neither
    seen.clear()
    a_cols = [c for c in cols if registry.BY_NAME[c].tier == "A"]
    from slparalog.eval.models import LRModel
    est = LRModel(a_cols, C=1.0).fit(np.random.default_rng(1).normal(size=(30, len(a_cols))), np.arange(30) % 2)
    da = tmp_path / "fakeA"
    model.save_model(da, est, feature_set="A", model_type="lr", features=a_cols, params={}, training={}, seeds={})
    with pytest.raises(Stop):
        predict.run(TOY, tmp_path / "o.tsv", da, data_dir=dd)
    assert seen["tier_b"] is False and seen["tier_c"] is False


def test_inference_process_never_imports_torch():
    code = ("import sys; import slparalog.cli, slparalog.predict, slparalog.train, slparalog.setup_data, slparalog.model;"
            "bad=[m for m in ('torch','transformers','lightgbm') if m in sys.modules]; assert not bad, bad")
    subprocess.run([sys.executable, "-c", code], check=True)


# ------------------------------------------------------------------------------ column-definition hash

def _rewrite_table(src, dst, **md_updates):
    import pyarrow.parquet as pq

    t = pq.read_table(src)
    md = json.loads(t.schema.metadata[b"slparalog"])
    for k, v in md_updates.items():
        if v is None:
            md.pop(k, None)
        else:
            md[k] = v
    pq.write_table(t.replace_schema_metadata({b"slparalog": json.dumps(md).encode()}), dst)
    return md


def test_new_artifacts_and_tables_record_columns_hash(env, lr_dir):
    import pyarrow.parquet as pq

    meta = json.loads((lr_dir / "model.json").read_text())
    assert meta["columns_hash"] == registry.columns_hash(meta["features"])
    tmeta = json.loads(pq.read_schema(env["feat"]).metadata[b"slparalog"])
    assert tmeta["columns_hash"] == registry.table_columns_hash(pq.read_schema(env["feat"]).names)
    assert json.loads(Path(str(env["feat"]) + ".meta.json").read_text())["columns_hash"] == tmeta["columns_hash"]
    rec = meta["training"]["feature_tables"][str(env["feat"])]
    assert rec["columns_hash"] == tmeta["columns_hash"]


def test_model_loads_when_only_the_global_registry_hash_changed(lr_dir, tmp_path):
    """Registering another tier / rewording a description changes registry_hash but not the model's columns."""
    d = _copy(lr_dir, tmp_path)
    _edit_meta(d, registry_hash="0" * 64)
    assert model.load_model(d).meta["registry_hash"] == "0" * 64


def test_model_column_definition_change_fails_loudly(lr_dir, tmp_path, monkeypatch):
    import dataclasses

    d = _copy(lr_dir, tmp_path)
    _edit_meta(d, columns_hash="0" * 64)
    with pytest.raises(model.ModelError, match="column-definition hash mismatch"):
        model.load_model(d)
    f = registry.BY_NAME["ident_global"]
    monkeypatch.setitem(registry.BY_NAME, "ident_global", dataclasses.replace(f, dtype="float32"))
    with pytest.raises(model.ModelError, match="column-definition hash mismatch"):
        model.load_model(lr_dir)


def test_train_accepts_table_with_stale_registry_hash_but_current_columns(env, tmp_path):
    p = tmp_path / "other_tier_edited.parquet"
    _rewrite_table(env["feat"], p, registry_hash="1" * 64)  # columns_hash left as recorded
    m = train.train({"human": [p]}, {"human": env["labels"]["human"]}, "A", "lr", tmp_path / "m")
    assert m["training"]["feature_tables"][str(p)]["registry_hash"] == "1" * 64  # kept as provenance
    assert m["columns_hash"] == registry.columns_hash(m["features"])


def test_train_rejects_table_whose_column_definitions_changed(env, tmp_path):
    p = tmp_path / "defs_changed.parquet"
    _rewrite_table(env["feat"], p, columns_hash="2" * 64)
    with pytest.raises(model.ModelError, match="column-definition hash mismatch"):
        train.train({"human": [p]}, {"human": env["labels"]["human"]}, "A", "lr", tmp_path / "m")


def test_train_rejects_dtype_change_of_a_table_column(env, tmp_path, monkeypatch):
    import dataclasses

    f = registry.BY_NAME["len_ratio"]
    monkeypatch.setitem(registry.BY_NAME, "len_ratio", dataclasses.replace(f, dtype="float32"))
    with pytest.raises(model.ModelError, match="column-definition hash mismatch"):
        train.train({"human": [env["feat"]]}, {"human": env["labels"]["human"]}, "A", "lr", tmp_path / "m")


def test_restamped_sidecar_is_honoured_for_legacy_tables(env, tmp_path):
    """Legacy table (no columns_hash in the parquet, stale registry hash): refused until its sidecar is restamped."""
    p = tmp_path / "legacy.parquet"
    md = _rewrite_table(env["feat"], p, registry_hash="1" * 64, columns_hash=None)
    side = Path(str(p) + ".meta.json")
    side.write_text(json.dumps(md))
    with pytest.raises(model.ModelError, match="registry"):
        model.check_feature_table(p)
    import pyarrow.parquet as pq

    side.write_text(json.dumps({**md, "columns_hash": registry.table_columns_hash(pq.read_schema(p).names),
                                "restamped": {"note": "test"}}))
    assert model.check_feature_table(p)["restamped"] == {"note": "test"}
    # a sidecar belonging to a different build (other schema/registry hash) is ignored
    side.write_text(json.dumps({**md, "registry_hash": "3" * 64, "columns_hash": registry.table_columns_hash(pq.read_schema(p).names)}))
    with pytest.raises(model.ModelError, match="registry"):
        model.check_feature_table(p)
