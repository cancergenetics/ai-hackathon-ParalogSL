"""Guard tests for scripts/eval/final_test.py (gate G15, reviews/round2/agent_failure_modes.md).

Everything runs on fakes: a throw-away git repository built in a temp directory (its own code, configs, model,
and a *fake* set of test labels copied from tests/fixtures/final_test_toy), never the project's real test labels.
Checked here: each precondition failing makes the script refuse before any labels are read; the sealed reader
needs the access record and the manifest; a second access is refused; the toy dry-run reproduces hand-computed
metrics. Protected file (config/protected.sha256, by the PI); allowlisted in tests/test_data_access.py.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("final_test", ROOT / "scripts" / "eval" / "final_test.py")
ft = importlib.util.module_from_spec(spec)
sys.modules["final_test"] = ft
spec.loader.exec_module(ft)

TOY = ROOT / "tests" / "fixtures" / "final_test_toy"
SEALED_NAMES = ["human_test.tsv", "yeast_test.tsv", "unseen_arabidopsis.tsv", "unseen_worm.tsv"]
GIT_ENV = {**os.environ, "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t", "GIT_COMMITTER_NAME": "t",
           "GIT_COMMITTER_EMAIL": "t@t", "GIT_CONFIG_GLOBAL": "/dev/null", "GIT_CONFIG_SYSTEM": "/dev/null"}


def g(repo: Path, *args: str) -> str:
    return subprocess.run(["git", *args], cwd=repo, env=GIT_ENV, capture_output=True, text=True, check=True).stdout.strip()


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ fake repository

def _write_protected(repo: Path, files: list[str]) -> None:
    lines = [f"{sha(repo / f)}  {f}" for f in files if (repo / f).exists()]
    (repo / "config" / "protected.sha256").write_text("\n".join(lines) + "\n")


@pytest.fixture(scope="module")
def template(tmp_path_factory) -> Path:
    """A complete fake project: commit A (code + configs), commit B (fake labels, dev data), commit C (model),
    tag final-model-v1 at C, pytest marker for C. Features sit in the git-ignored data/interim as in the project."""
    from slparalog.train import train

    repo = tmp_path_factory.mktemp("fake") / "repo"
    (repo / "config").mkdir(parents=True)
    for rel in ("scripts/eval/final_test.py", "tests/test_final_test_guard.py", "src/slparalog/eval/metrics.py",
                "tests/test_data_access.py"):
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / rel, repo / rel)
    for n in ("success.yaml", "models.yaml", "labels_unseen.yaml", "seeds.yaml"):
        shutil.copy(ROOT / "config" / n, repo / "config" / n)
    (repo / ".gitignore").write_text("data/interim/\n__pycache__/\n")
    (repo / "RESEARCH_TREE.md").write_text("# tree\n")
    g(repo, "init", "-q", "-b", "main")

    sealed = repo / "data" / "labels" / "sealed"
    dev = repo / "data" / "labels" / "dev"
    sealed.mkdir(parents=True)
    dev.mkdir(parents=True)
    for n in SEALED_NAMES:
        shutil.copy(TOY / n, sealed / n)
    (sealed / "MANIFEST.sha256").write_text("".join(f"{sha(sealed / n)}  {n}\n" for n in SEALED_NAMES))
    for sp in ("human", "yeast"):
        shutil.copy(TOY / f"{sp}_test.tsv", dev / f"{sp}_dev.tsv")
    for sp in ("human", "yeast", "arabidopsis", "worm"):
        fd = repo / "data" / "interim" / "features" / sp
        fd.mkdir(parents=True)
        shutil.copy(TOY / "features" / sp / "tierA.parquet", fd / "tierA.parquet")
    prot = list(ft.REQUIRED_PROTECTED) + ["data/labels/sealed/MANIFEST.sha256"]
    _write_protected(repo, prot)

    # commit A: code + configs; commit B: labels
    g(repo, "add", "-A", "--", ".gitignore", "RESEARCH_TREE.md", "config", "scripts", "tests", "src")
    g(repo, "commit", "-q", "-m", "A: code and config")
    g(repo, "add", "-A")
    g(repo, "commit", "-q", "-m", "B: fake labels")
    commit_b = g(repo, "rev-parse", "HEAD")

    # model: a real LR artifact fitted on the toy tables, then re-stamped with the fake repository's commit
    mdir = repo / "results" / "frozen" / "model"
    feats = {sp: [repo / "data" / "interim" / "features" / sp / "tierA.parquet"] for sp in ("human", "yeast")}
    labs = {sp: dev / f"{sp}_dev.tsv" for sp in ("human", "yeast")}
    train(feats, labs, "A+B", "gbdt", mdir, models_config=repo / "config" / "models.yaml")
    meta = json.loads((mdir / "model.json").read_text())
    meta["git"] = {"commit": commit_b, "dirty": False}
    (mdir / "model.json").write_text(json.dumps(meta, indent=1))
    g(repo, "add", "-A")
    g(repo, "commit", "-q", "-m", "C: frozen model")
    g(repo, "tag", "final-model-v1")
    ft.marker_path(ft.Ctx(root=repo)).write_text(json.dumps({"head": g(repo, "rev-parse", "HEAD"), "passed": True}))
    return repo


@pytest.fixture()
def repo(template, tmp_path) -> Path:
    r = tmp_path / "repo"
    shutil.copytree(template, r)
    return r


@pytest.fixture()
def no_sealed_reads(monkeypatch):
    """Any parse of a sealed table is a test failure (used for every refusal); the hash-only check is allowed."""
    real = ft.read_sealed

    def boom(ctx, name, verify_only=False):
        if verify_only:
            return real(ctx, name, verify_only=True)
        raise AssertionError("the sealed reader parsed a table although a precondition failed")
    monkeypatch.setattr(ft, "read_sealed", boom)


def retag(repo: Path, message: str, *paths: str, marker: bool = True) -> None:
    """Commit the working-tree change, move the tag to it (fake repo only) and refresh the pytest marker."""
    g(repo, "add", "-A")
    g(repo, "commit", "-q", "-m", message)
    g(repo, "tag", "-f", "final-model-v1")
    if marker:
        ft.marker_path(ft.Ctx(root=repo)).write_text(json.dumps({"head": g(repo, "rev-parse", "HEAD"), "passed": True}))


def refused_with(repo: Path, needle: str, capsys=None) -> list[str]:
    ctx = ft.Ctx(root=repo)
    fails = ft.run_preconditions(ctx)
    assert any(needle in f for f in fails), f"{needle!r} not among {fails}"
    assert ft.main([], root=repo) == 3  # the script itself refuses...
    assert not ctx.log_path.exists()  # ...without writing an access record
    return fails


# ------------------------------------------------------------------------------------------------ preconditions

def test_happy_path_preconditions_pass(repo):
    assert ft.run_preconditions(ft.Ctx(root=repo)) == []


def test_refuses_dirty_tree(repo, no_sealed_reads):
    (repo / "scratch.txt").write_text("x")
    refused_with(repo, "not clean")


def test_refuses_modified_tracked_file(repo, no_sealed_reads):
    (repo / "RESEARCH_TREE.md").write_text("# changed\n")
    refused_with(repo, "not clean")


def test_refuses_missing_tag(repo, no_sealed_reads):
    g(repo, "tag", "-d", "final-model-v1")
    refused_with(repo, "does not exist")


def test_refuses_head_past_tag_with_code_change(repo, no_sealed_reads):
    (repo / "extra.py").write_text("x = 1\n")
    g(repo, "add", "-A")
    g(repo, "commit", "-q", "-m", "code after the tag")
    ft.marker_path(ft.Ctx(root=repo)).write_text(json.dumps({"head": g(repo, "rev-parse", "HEAD"), "passed": True}))
    refused_with(repo, "non-results paths")


def test_head_past_tag_with_results_only_is_allowed(repo):
    (repo / "results" / "note.txt").write_text("x")
    g(repo, "add", "-A")
    g(repo, "commit", "-q", "-m", "results only")
    ft.marker_path(ft.Ctx(root=repo)).write_text(json.dumps({"head": g(repo, "rev-parse", "HEAD"), "passed": True}))
    assert ft.run_preconditions(ft.Ctx(root=repo)) == []


def test_refuses_modified_model_after_tag(repo, no_sealed_reads):
    mj = repo / "results" / "frozen" / "model" / "model.joblib"
    mj.write_bytes(mj.read_bytes() + b"x")
    g(repo, "add", "-A")
    g(repo, "commit", "-q", "-m", "model touched after the tag")
    ft.marker_path(ft.Ctx(root=repo)).write_text(json.dumps({"head": g(repo, "rev-parse", "HEAD"), "passed": True}))
    fails = refused_with(repo, "sha256 differs")
    assert any("non-results paths" in f for f in fails)  # results/frozen is frozen at the tag too


def test_refuses_model_checksum_mismatch(repo, no_sealed_reads):
    p = repo / "results" / "frozen" / "model" / "model.json"
    meta = json.loads(p.read_text())
    meta["model_sha256"] = "0" * 64
    p.write_text(json.dumps(meta))
    retag(repo, "bad checksum")
    refused_with(repo, "sha256 differs")


def test_refuses_model_commit_not_ancestor(repo, no_sealed_reads):
    head = g(repo, "rev-parse", "HEAD")
    g(repo, "checkout", "-q", "-b", "side", "HEAD~2")
    (repo / "side.txt").write_text("s")
    g(repo, "add", "-A")
    g(repo, "commit", "-q", "-m", "side")
    side = g(repo, "rev-parse", "HEAD")
    g(repo, "checkout", "-q", "main")
    assert g(repo, "rev-parse", "HEAD") == head
    p = repo / "results" / "frozen" / "model" / "model.json"
    meta = json.loads(p.read_text())
    meta["git"]["commit"] = side
    p.write_text(json.dumps(meta))
    retag(repo, "model from another branch")
    refused_with(repo, "not an ancestor of HEAD")


def test_refuses_model_from_dirty_tree(repo, no_sealed_reads):
    p = repo / "results" / "frozen" / "model" / "model.json"
    meta = json.loads(p.read_text())
    meta["git"]["dirty"] = True
    p.write_text(json.dumps(meta))
    retag(repo, "dirty model")
    refused_with(repo, "dirty git tree")


def test_refuses_missing_b1_variant(repo, no_sealed_reads):
    p = repo / "results" / "frozen" / "model" / "model.json"
    meta = json.loads(p.read_text())
    del meta["b1_variant"]
    p.write_text(json.dumps(meta))
    retag(repo, "no B1 variant")
    refused_with(repo, "b1_variant")


def test_refuses_changed_training_data(repo, no_sealed_reads):
    f = repo / "data" / "interim" / "features" / "human" / "tierA.parquet"
    df = pd.read_parquet(f)
    df.loc[0, "bits_pair"] += 1.0
    df.to_parquet(f)  # ignored by git, so the tree stays clean: only the checksum check can notice
    refused_with(repo, "training feature table human/tierA")


def test_refuses_protected_hash_mismatch(repo, no_sealed_reads):
    (repo / "config" / "success.yaml").write_text((repo / "config" / "success.yaml").read_text() + "\n# edited\n")
    retag(repo, "success.yaml edited")
    refused_with(repo, "protected file changed without a D<n> note")


def test_protected_change_with_decision_note_verifies(repo):
    (repo / "config" / "success.yaml").write_text((repo / "config" / "success.yaml").read_text() + "\n# edited\n")
    (repo / "RESEARCH_TREE.md").write_text(f"- D99: success.yaml {sha(repo / 'config' / 'success.yaml')[:12]}\n")
    retag(repo, "success.yaml edited, noted")
    assert ft.check_protected(ft.Ctx(root=repo)) == []


def test_refuses_protected_list_without_the_script_itself(repo, no_sealed_reads):
    pf = repo / "config" / "protected.sha256"
    pf.write_text("".join(ln + "\n" for ln in pf.read_text().splitlines() if "final_test.py" not in ln))
    retag(repo, "script unprotected")
    refused_with(repo, "protected list lacks scripts/eval/final_test.py")


def test_refuses_without_pytest_marker(repo, no_sealed_reads):
    ft.marker_path(ft.Ctx(root=repo)).unlink()
    refused_with(repo, "no pytest marker")


def test_refuses_pytest_marker_for_other_commit(repo, no_sealed_reads):
    ft.marker_path(ft.Ctx(root=repo)).write_text(json.dumps({"head": "0" * 40, "passed": True}))
    refused_with(repo, "not a pass at HEAD")


def test_refuses_failed_pytest_marker(repo, no_sealed_reads):
    ft.marker_path(ft.Ctx(root=repo)).write_text(json.dumps({"head": g(repo, "rev-parse", "HEAD"), "passed": False}))
    refused_with(repo, "not a pass at HEAD")


def test_refuses_prior_started_access(repo, no_sealed_reads):
    ctx = ft.Ctx(root=repo)
    ft.log_append(ctx, "access_started")
    retag(repo, "log of an earlier access", marker=True)
    fails = ft.run_preconditions(ctx)
    assert any("already records an access" in f for f in fails)
    assert ft.main([], root=repo) == 3


def test_refuses_existing_results(repo, no_sealed_reads):
    (repo / "results" / "final").mkdir()
    (repo / "results" / "final" / "final_test.json").write_text("{}")
    retag(repo, "results exist")
    refused_with(repo, "final_test.json already exists")


def test_refuses_labels_committed_after_tag(repo, no_sealed_reads):
    first = g(repo, "rev-list", "--max-parents=0", "HEAD")
    g(repo, "tag", "-f", "final-model-v1", first)  # tag before the labels were committed
    ft.marker_path(ft.Ctx(root=repo)).write_text(json.dumps({"head": g(repo, "rev-parse", "HEAD"), "passed": True}))
    refused_with(repo, "was first committed after tag")


def test_refuses_changed_success_criterion_config(repo, no_sealed_reads):
    p = repo / "config" / "success.yaml"
    p.write_text(p.read_text().replace("percentile", "basic"))
    retag(repo, "ci type changed")
    refused_with(repo, "ci_type")


def test_refusals_are_recorded_outside_the_tree(repo, no_sealed_reads):
    (repo / "scratch.txt").write_text("x")
    ft.main([], root=repo)
    assert (repo / ".git" / "final_test_refusals.log").exists()
    assert not (repo / "results" / "final").exists()


def _mini_repo(tmp_path: Path, body: str) -> Path:
    r = tmp_path / "mini"
    (r / "tests").mkdir(parents=True)
    (r / "tests" / "test_x.py").write_text(body)
    (r / ".gitignore").write_text("__pycache__/\n")
    g(r, "init", "-q", "-b", "main")
    g(r, "add", "-A")
    g(r, "commit", "-q", "-m", "m")
    return r


def test_record_tests_writes_marker_only_on_a_clean_pass(tmp_path):
    r = _mini_repo(tmp_path, "def test_ok():\n    assert True\n")
    ctx = ft.Ctx(root=r)
    assert ft.check_tests_marker(ctx)
    assert ft.main(["--record-tests"], root=r) == 0
    assert ft.check_tests_marker(ctx) == []
    (r / "new.txt").write_text("x")  # dirty tree: no run, no marker refresh
    with pytest.raises(ft.Refusal, match="not clean"):
        ft.record_tests(ctx)
    g(r, "add", "-A")
    g(r, "commit", "-q", "-m", "next")
    assert ft.check_tests_marker(ctx)  # the marker names the previous commit


def test_record_tests_refuses_on_a_failing_suite(tmp_path):
    r = _mini_repo(tmp_path, "def test_bad():\n    assert False\n")
    assert ft.main(["--record-tests"], root=r) == 3
    assert not ft.marker_path(ft.Ctx(root=r)).exists()


def refused_in_preflight(repo: Path, needle: str, capsys) -> str:
    """Preconditions pass, the preflight refuses: exit 3, no access record, no sealed read."""
    ctx = ft.Ctx(root=repo)
    assert ft.run_preconditions(ctx) == []
    assert ft.main([], root=repo) == 3
    err = capsys.readouterr().err
    assert "preflight" in err and needle in err, err
    assert not ctx.log_path.exists()
    return err


def _edit_table(repo: Path, sp: str, fn) -> None:
    f = repo / "data" / "interim" / "features" / sp / "tierA.parquet"
    df = pd.read_parquet(f)
    fn(df)
    df.to_parquet(f)  # git-ignored: the tree stays clean, only the preflight can notice


def test_preflight_refuses_nan_b1_input_in_an_unseen_species(repo, no_sealed_reads, capsys):
    _edit_table(repo, "arabidopsis", lambda d: d.__setitem__("ident_global", d["ident_global"].where(d.index != 0)))
    refused_in_preflight(repo, "arabidopsis: non-finite values", capsys)


def test_preflight_refuses_nan_b2_input(repo, no_sealed_reads, capsys):
    _edit_table(repo, "worm", lambda d: d.__setitem__("rbh", d["rbh"].astype(float).where(d.index != 1)))
    assert "B2 input rbh" in refused_in_preflight(repo, "worm: non-finite values", capsys)


def test_preflight_refuses_missing_model_feature_column(repo, no_sealed_reads, capsys):
    _edit_table(repo, "worm", lambda d: d.__delitem__("len_min"))
    refused_in_preflight(repo, "worm: feature tables lack columns", capsys)


def test_preflight_refuses_non_finite_model_scores(repo, no_sealed_reads, capsys, monkeypatch):
    real = ft.RealBackend.score
    monkeypatch.setattr(ft.RealBackend, "score", lambda self, df: real(self, df) * np.nan)
    refused_in_preflight(repo, "'model'", capsys)


def test_preflight_refuses_non_finite_yeast_only_scores(repo, no_sealed_reads, capsys, monkeypatch):
    monkeypatch.setattr(ft.RealBackend, "fit_species_only", lambda self, sp: (lambda d: np.full(len(d), np.inf)))
    refused_in_preflight(repo, "yeast_only_model", capsys)


def test_preflight_covers_all_four_species(toy):
    assert set(toy[0]["preflight"]) == {"human", "yeast", "arabidopsis", "worm"}
    assert toy[0]["preflight"]["human"]["called_pairs_scored"] == 8
    with pytest.raises(ft.Refusal, match="preflight"):
        class Bad(ft.ToyBackend):
            def score(self, df):
                return df["bits_pair"].to_numpy(float) * np.nan
        ft.preflight(Bad())


def test_cli_options_are_checked():
    assert ft.main(["--dry-run-toy", "--record-tests"]) == 2
    assert ft.main(["--out-dir", "x"]) == 2


# ------------------------------------------------------------------------------------------------ R1: frozen spec

def _edit_model_json(repo: Path, fn) -> None:
    p = repo / "results" / "frozen" / "model" / "model.json"
    meta = json.loads(p.read_text())
    fn(meta)
    p.write_text(json.dumps(meta))
    retag(repo, "model.json edited")


def test_frozen_spec_matches_config_in_the_template(repo):
    meta = json.loads((repo / "results/frozen/model/model.json").read_text())
    cfg = ft.load_yaml(repo / "config" / "models.yaml")
    assert meta["feature_set"] == "A+B" and meta["model_type"] == "gbdt" and meta["training"]["tier"] == "primary"
    assert meta["params"] == cfg["anchoring"]["ensemble"]["gbdt_params"]
    assert meta["params"] == {"l2_regularization": 1.0, "learning_rate": 0.03, "max_iter": 100, "max_leaf_nodes": 7,
                              "min_samples_leaf": 50}  # GBDT_AB_17
    assert meta["seeds"]["seed_used"] == 11 and meta["b1_variant"] == "ident_global"
    assert not [f for f in meta["features"] if f.startswith("tb_asym")]
    assert ft.frozen_spec_failures(ft.Ctx(root=repo), meta) == []


@pytest.mark.parametrize("name,edit,needle", [
    ("feature_set", lambda m: m.__setitem__("feature_set", "A"), "feature_set"),
    ("model_type", lambda m: m.__setitem__("model_type", "lr"), "model_type"),
    ("params", lambda m: m["params"].__setitem__("max_iter", 300), "params"),
    ("params_extra", lambda m: m["params"].__setitem__("max_depth", 3), "params"),
    ("features_drop", lambda m: m["features"].pop(), "features"),
    ("features_asym", lambda m: m["features"].append("tb_asym_abs"), "features"),
    ("tier", lambda m: m["training"].__setitem__("tier", "lenient"), "training.tier"),
    ("seed", lambda m: m["seeds"].__setitem__("seed_used", 23), "seeds"),
    ("b1_variant", lambda m: m.__setitem__("b1_variant", "ident_local_max"), "b1_variant"),
])
def test_refuses_model_that_is_not_the_frozen_spec(repo, no_sealed_reads, name, edit, needle):
    _edit_model_json(repo, edit)
    refused_with(repo, "the frozen spec" if needle not in ("training.tier", "seeds") else needle)


def test_frozen_spec_follows_config_not_literals(repo, no_sealed_reads):
    cfg = repo / "config" / "models.yaml"
    cfg.write_text(cfg.read_text().replace("excluded_features:", "excluded_features:\n  ident_global: \"x\"", 1))
    retag(repo, "config changed")
    fails = ft.frozen_spec_failures(ft.Ctx(root=repo), json.loads((repo / "results/frozen/model/model.json").read_text()))
    assert any("features" in f for f in fails)


# ------------------------------------------------------------------------------------------------ columns_hash tables

def _stamp(path: Path, meta: dict | None, sidecar: dict | None) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    t = pa.Table.from_pandas(pd.read_parquet(path))
    md = dict(t.schema.metadata or {})
    if meta is not None:
        md[b"slparalog"] = json.dumps(meta).encode()
    pq.write_table(t.replace_schema_metadata(md), path)
    side = Path(str(path) + ".meta.json")
    side.unlink(missing_ok=True)
    if sidecar is not None:
        side.write_text(json.dumps(sidecar))


def test_table_validation_uses_columns_hash_not_the_global_registry_hash(repo):
    from slparalog.features import registry

    f = repo / "data/interim/features/worm/tierA.parquet"
    good = registry.table_columns_hash(pd.read_parquet(f).columns)
    old = {"registry_hash": "d" * 64, "schema_hash": None}  # embedded record predates columns_hash and is stale
    backend = ft.RealBackend(ft.Ctx(root=repo))
    _stamp(f, old, {**old, "columns_hash": good})  # restamp sidecar only: accepted although the global hash differs
    assert backend.table_failures("worm") == []
    _stamp(f, {**old, "columns_hash": good}, None)  # embedded columns_hash: accepted
    assert backend.table_failures("worm") == []
    _stamp(f, old, {**old, "columns_hash": "0" * 64})  # wrong column definitions
    assert any("column-definition hash mismatch" in x for x in backend.table_failures("worm"))
    _stamp(f, old, None)  # no columns_hash anywhere: falls back to the whole-registry hash (stale -> refused)
    assert any("registry" in x for x in backend.table_failures("worm"))


def test_training_table_record_in_model_json_is_checked(repo):
    backend = ft.RealBackend(ft.Ctx(root=repo))
    key = next(k for k in backend.meta["training"]["feature_tables"] if "/human/" in k)
    backend.meta["training"]["feature_tables"][key] = {"columns_hash": "0" * 64}
    assert any("columns_hash differs from model.json" in x for x in backend.table_failures("human"))
    backend.meta["training"]["feature_tables"][key] = {"schema_hash": "x" * 64}  # older artifacts
    assert any("schema hash differs" in x for x in backend.table_failures("human"))


# ------------------------------------------------------------------------------------------------ R3

def test_sealed_hashes_are_verified_before_the_access_record(repo, capsys):
    p = repo / "data/labels/sealed/unseen_worm.tsv"
    p.write_text(p.read_text() + "x\ty\t0\n")
    retag(repo, "unseen file changed after the manifest")
    ctx = ft.Ctx(root=repo)
    assert ft.run_preconditions(ctx) == []  # a hash problem is not visible to the git-level checks
    assert ft.main([], root=repo) == 3
    assert "unseen_worm.tsv" in capsys.readouterr().err
    assert not ctx.log_path.exists()  # the access was not consumed


def test_hash_only_verification_does_not_parse_or_keep_anything(repo):
    ctx = ft.Ctx(root=repo)
    d = ft.read_sealed(ctx, "human_test.tsv", verify_only=True)
    assert d == sha(repo / "data/labels/sealed/human_test.tsv") and ctx.sealed_read == {}
    assert ft.verify_sealed(ctx) == []
    with pytest.raises(ft.Refusal, match="before the access_started"):
        ft.read_sealed(ctx, "human_test.tsv")


def test_primary_is_written_before_any_secondary_and_survives_secondary_errors(repo, monkeypatch):
    seen = {}
    out = repo / "results" / "final" / "final_test.json"

    def first_secondary(st):
        d = json.loads(out.read_text())  # already on disk, before this analysis runs
        seen.update(status=d["status"], T=d["primary"]["T"]["point"], keys=set(d))
        raise RuntimeError("secondary broke")

    monkeypatch.setattr(ft, "sec_delta_vs_B2", first_secondary)
    monkeypatch.setattr(ft, "sec_end_to_end", lambda st: 1 / 0)
    assert ft.main([], root=repo) == 0  # the access completes
    assert seen["status"] == "primary_only" and "secondary" in seen["keys"]
    res = json.loads(out.read_text())
    assert res["status"] == "complete_with_secondary_errors"
    assert set(res["secondary_errors"]) == {"delta_vs_B2", "end_to_end"}
    assert "secondary broke" in res["secondary_errors"]["delta_vs_B2"]["error"]
    assert res["primary"]["T"]["point"] == pytest.approx(seen["T"])  # primary unchanged
    ok = set(ft.SECONDARY_NAMES) - {"delta_vs_B2", "end_to_end"}
    assert ok <= set(res["secondary"])  # every other secondary still ran
    assert "Secondary analyses that failed" in (repo / "results/final/final_test.md").read_text()
    rows = ft.read_log(ft.Ctx(root=repo))
    assert [r["event"] for r in rows] == ["access_started", "access_completed"]
    assert "complete_with_secondary_errors" in rows[1]["note"]
    assert rows[0]["sealed_manifest_sha256"] == sha(repo / "data/labels/sealed/MANIFEST.sha256")


def test_every_secondary_failing_leaves_the_primary(tmp_path, monkeypatch):
    for n in ft.SECONDARY_NAMES:
        monkeypatch.setattr(ft, "sec_" + n, lambda st: (_ for _ in ()).throw(ValueError("x")))
    res = ft.run_toy(tmp_path)
    assert res["status"] == "complete_with_secondary_errors" and set(res["secondary_errors"]) == set(ft.SECONDARY_NAMES)
    on_disk = json.loads((tmp_path / "final_test.json").read_text())
    assert on_disk["primary"]["T"]["point"] == pytest.approx(11 / 36)
    assert "failed" in (tmp_path / "final_test.md").read_text()


def test_primary_failure_before_hook_writes_nothing(tmp_path, monkeypatch):
    monkeypatch.setattr(ft, "primary_block", lambda *a, **k: (_ for _ in ()).throw(ValueError("primary broke")))
    with pytest.raises(ValueError):
        ft.run_toy(tmp_path)
    assert not (tmp_path / "final_test.json").exists()


# ------------------------------------------------------------------------------------------------ sealed reader

def test_sealed_reader_needs_the_access_record(repo):
    ctx = ft.Ctx(root=repo)
    with pytest.raises(ft.Refusal, match="before the access_started"):
        ft.read_sealed(ctx, "human_test.tsv")


def test_sealed_reader_checks_manifest_and_names(repo):
    ctx = ft.Ctx(root=repo, access_started=True)
    assert set(ft.read_sealed(ctx, "MANIFEST.sha256")) == set(SEALED_NAMES)
    df = ft.read_sealed(ctx, "human_test.tsv")
    assert len(df) == 10 and ctx.sealed_read["human_test.tsv"] == sha(repo / "data/labels/sealed/human_test.tsv")
    with pytest.raises(ft.Refusal, match="not a file directly"):
        ft.read_sealed(ctx, "../dev/human_dev.tsv")
    with pytest.raises(ft.Refusal, match="not listed"):
        ft.read_sealed(ctx, "other.tsv")
    p = repo / "data/labels/sealed/yeast_test.tsv"
    p.write_text(p.read_text() + "x\tx\tpos\tgx\n")
    with pytest.raises(ft.Refusal, match="does not match the sealed manifest"):
        ft.read_sealed(ctx, "yeast_test.tsv")


def test_sealed_reader_checks_manifest_against_protected_hash(repo):
    m = repo / "data/labels/sealed/MANIFEST.sha256"
    m.write_text(m.read_text() + "\n")
    with pytest.raises(ft.Refusal, match="protected"):
        ft.read_sealed(ft.Ctx(root=repo, access_started=True), "human_test.tsv")


def test_only_one_function_opens_sealed_files():
    src = (ROOT / "scripts" / "eval" / "final_test.py").read_text()
    assert src.count("sealed_dir") >= 1
    # every use of the sealed directory outside read_sealed is a path computation for git/protected checks
    body = src.split("def read_sealed")[1].split("\ndef ")[0]
    assert "read_bytes" in body
    rest = src.replace(body, "")
    assert "sealed / " not in rest and "sealed_dir /" not in rest


# ------------------------------------------------------------------------------------------------ the access itself

def test_full_access_on_fake_labels_then_second_access_refused(repo, monkeypatch):
    ctx = ft.Ctx(root=repo)
    calls = []
    real = ft.read_sealed

    def spy(c, name, verify_only=False):
        if verify_only:  # hash-only check, allowed (and done) before the access record
            assert ft.read_log(c) == []
            return real(c, name, verify_only=True)
        rows = ft.read_log(c)
        assert [r["event"] for r in rows] == ["access_started"], "sealed file read before the access record"
        calls.append(name)
        return real(c, name)

    monkeypatch.setattr(ft, "read_sealed", spy)
    assert ft.main([], root=repo) == 0
    assert sorted(calls) == sorted(SEALED_NAMES)  # each sealed table read once, through the reader
    rows = ft.read_log(ctx)
    assert [r["event"] for r in rows] == ["access_started", "access_completed"]
    assert rows[0]["purpose"] == "primary" and rows[0]["model_sha256"] == json.loads(
        (repo / "results/frozen/model/model.json").read_text())["model_sha256"]
    res = json.loads((repo / "results" / "final" / "final_test.json").read_text())
    assert (repo / "results" / "final" / "final_test.md").exists()
    p = res["primary"]
    for sp, (n, npos) in {"human": (8, 3), "yeast": (8, 2)}.items():
        r = p["per_species"][sp]
        assert (r["n"], r["n_pos"], r["n_neg"]) == (n, npos, n - npos)
        assert r["ap_B1"] == pytest.approx({"human": 13 / 18, "yeast": 0.5}[sp])  # B1 is scored by the same code
    assert p["T"]["point"] == pytest.approx(np.mean([p["per_species"][s]["delta_vs_B1"]["point"] for s in ("human", "yeast")]))
    assert set(res["secondary"]["unseen_species"]) == {"arabidopsis", "worm"}
    assert res["secondary"]["unseen_species"]["arabidopsis"]["strict"]["counts"]["called_pos"] == 2
    assert res["secondary"]["unseen_species"]["arabidopsis"]["strict"]["counts"]["not_called_neg"] == 1
    yo = res["secondary"]["yeast_only_model_on_human"]
    assert yo["n"] == 8 and set(yo["ap"]) == {"model_yeast_only", "model", "B1"}
    assert set(res["model"]["sealed_files_sha256"]) == set(SEALED_NAMES)
    assert {k: v["called_pairs_scored"] for k, v in res["preflight"].items()} == {
        "human": 8, "yeast": 8, "arabidopsis": 6, "worm": 5}
    # the access log and results are uncommitted: even a pytest marker cannot make the next run pass
    ft.marker_path(ctx).write_text(json.dumps({"head": g(repo, "rev-parse", "HEAD"), "passed": True}))
    monkeypatch.setattr(ft, "read_sealed", lambda *a, **k: (_ for _ in ()).throw(AssertionError("second access read labels")))
    assert ft.main([], root=repo) == 3
    assert [r["event"] for r in ft.read_log(ctx)] == ["access_started", "access_completed"]
    # committed, the second access is refused because of the log itself
    g(repo, "add", "-A")
    g(repo, "commit", "-q", "-m", "first access")
    ft.marker_path(ctx).write_text(json.dumps({"head": g(repo, "rev-parse", "HEAD"), "passed": True}))
    fails = ft.run_preconditions(ctx)
    assert [f for f in fails if "already records an access" in f] and [f for f in fails if "final_test.json already exists" in f]
    assert len(fails) == 2
    assert ft.main([], root=repo) == 3
    assert [r["event"] for r in ft.read_log(ctx)] == ["access_started", "access_completed"]


def test_crashed_access_still_counts_as_an_access(repo, monkeypatch):
    ctx = ft.Ctx(root=repo)

    def explode(*a, **k):
        raise RuntimeError("boom")

    monkeypatch.setattr(ft, "compute", explode)
    with pytest.raises(RuntimeError, match="boom"):
        ft.main([], root=repo)
    assert [r["event"] for r in ft.read_log(ctx)] == ["access_started", "access_failed"]
    g(repo, "add", "-A")
    g(repo, "commit", "-q", "-m", "failed access logged")
    ft.marker_path(ctx).write_text(json.dumps({"head": g(repo, "rev-parse", "HEAD"), "passed": True}))
    assert any("already records an access" in f for f in ft.run_preconditions(ctx))
    assert ft.main([], root=repo) == 3


# ------------------------------------------------------------------------------------------------ toy dry run

@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    out = tmp_path_factory.mktemp("toy")
    return ft.run_toy(out), out


def test_toy_known_answers(toy):
    res, out = toy
    p = res["primary"]
    h, y = p["per_species"]["human"], p["per_species"]["yeast"]
    # human: positives P0,P2,P5; model ranks them first (AP 1); B1 orders labels 1,0,1,0,0,1,0,0 (AP 13/18);
    # B2 orders 0,1,0,1,1,0,0,0 (AP 8/15)
    assert (h["n"], h["n_pos"], h["n_neg"]) == (8, 3, 5) and h["prevalence"] == pytest.approx(3 / 8)
    assert h["ap_model"] == pytest.approx(1.0, abs=1e-12)
    assert h["ap_B1"] == pytest.approx(13 / 18, abs=1e-12)
    assert h["ap_B2"] == pytest.approx(8 / 15, abs=1e-12)
    assert h["delta_vs_B1"]["point"] == pytest.approx(5 / 18, abs=1e-12)
    assert h["relative_gain_vs_B1"] == pytest.approx(5 / 13, abs=1e-12)
    # yeast: positives Q1,Q3; model orders 1,0,1,0.. (AP 5/6); B1 0,1,0,1.. (AP 1/2); B2 ranks both first (AP 1)
    assert (y["n"], y["n_pos"], y["n_neg"]) == (8, 2, 6) and y["prevalence"] == pytest.approx(0.25)
    assert y["ap_model"] == pytest.approx(5 / 6, abs=1e-12)
    assert y["ap_B1"] == pytest.approx(0.5, abs=1e-12)
    assert y["ap_B2"] == pytest.approx(1.0, abs=1e-12)
    assert y["delta_vs_B1"]["point"] == pytest.approx(1 / 3, abs=1e-12)
    assert y["delta_vs_B2"]["point"] == pytest.approx(-1 / 6, abs=1e-12)
    # T, relative gain, criteria
    assert p["T"]["point"] == pytest.approx(11 / 36, abs=1e-12)
    assert p["relative_gain_mean"] == pytest.approx((5 / 13 + 2 / 3) / 2, abs=1e-12)
    assert p["relative_gain_threshold"] == 0.20
    assert p["T_vs_B2"]["point"] == pytest.approx((7 / 15 - 1 / 6) / 2, abs=1e-12)
    assert p["T"]["resamples"] == 10000 and p["T"]["seed"] == 11 and p["T"]["ci_type"] == "percentile"
    assert p["T"]["strata"] == ["species_composition", "group_has_positive"]
    assert p["T"]["lo"] <= p["T"]["point"] <= p["T"]["hi"] and p["T"]["lo"] > 0
    assert p["criteria"] == {"ci_of_T_excludes_0": True, "delta_ap_positive_in_each_species": True,
                             "relative_gain_ge_threshold": True}
    assert p["pass"] and p["clearly_beats_identity"] and p["verdict"].startswith("CLEARLY BEATS")
    # population counts: one labelled-but-uncalled pair per species, one non-primary human pair dropped
    assert res["population"]["human"]["labelled"] == 9 and res["population"]["human"]["not_called_pos"] == 1
    assert res["population"]["yeast"]["not_called"] == 1


def test_toy_secondary_known_answers(toy):
    sec = toy[0]["secondary"]
    # >=3 studies: P0..P4 (positives P0,P2): model AP 1, B1 AP 5/6
    g3 = sec["human_ge3_studies"]
    assert (g3["n"], g3["n_pos"], g3["n_neg"]) == (5, 2, 3)
    assert g3["ap"]["model"] == pytest.approx(1.0) and g3["ap"]["B1"] == pytest.approx(5 / 6)
    assert g3["delta_vs_B1"]["model"]["point"] == pytest.approx(1 / 6)
    # exact study-set cells (sorted study strings): {P1,P2,P3} and {P0,P4,P6,P7}; P5 sits alone (no negative: excluded)
    c = sec["human_exact_study_set_cells"]
    assert (c["n_cells_used"], c["n_cells_total"], c["n_pairs"], c["n_pos"], c["n_pos_in_excluded_cells"]) == (2, 3, 7, 2, 1)
    assert c["ap"]["model"] == pytest.approx(1.0) and c["ap"]["B1"] == pytest.approx(0.75)
    assert c["delta_vs_B1"]["point"] == pytest.approx(0.25)
    # yeast-only stand-in ranks P2,P0,P1,P3,P5.. -> AP 13/15 on the human pairs
    yo = sec["yeast_only_model_on_human"]
    assert yo["ap"]["model_yeast_only"] == pytest.approx(13 / 15) and yo["ap"]["model"] == pytest.approx(1.0)
    assert yo["delta_vs_B1"]["model_yeast_only"]["point"] == pytest.approx(13 / 15 - 13 / 18)
    assert yo["n"] == 8 and yo["n_pos"] == 3
    # delta vs B2
    assert sec["delta_vs_B2"]["T"]["point"] == pytest.approx((7 / 15 - 1 / 6) / 2)
    # unseen: groups are gene-graph components (A0 and A1 share a gene), NA labels dropped, uncalled pairs counted
    ara = sec["unseen_species"]["arabidopsis"]
    s = ara["strict"]
    assert (s["n"], s["n_pos"], s["n_neg"], s["n_groups"]) == (5, 2, 3, 4)
    assert s["ap"]["model"] == pytest.approx(1.0) and s["ap"]["B1"] == pytest.approx(0.5)
    assert s["delta_vs_B1"]["model"]["point"] == pytest.approx(0.5)
    assert s["counts"]["not_called_neg"] == 1 and s["n_labelled"] == 6
    b = ara["broad"]
    assert (b["n"], b["n_pos"], b["n_groups"]) == (6, 3, 5)
    assert b["ap"]["model"] == pytest.approx(13 / 15) and b["ap"]["B1"] == pytest.approx(8 / 15)
    w = sec["unseen_species"]["worm"]["primary"]
    assert (w["n"], w["n_pos"], w["n_neg"]) == (5, 2, 3) and w["counts"]["not_called_neg"] == 1
    assert w["ap"]["model"] == pytest.approx(1.0) and w["ap"]["B1"] == pytest.approx(0.4)  # constant B1 -> prevalence
    assert "warning" in w  # fewer than 10 positives


def test_toy_end_to_end_tiers_and_loso_known_answers(toy):
    sec = toy[0]["secondary"]
    # end-to-end: labelled primary pairs P0..P8 (P8 = caller-missed positive, scored at the minimum, last)
    e = sec["end_to_end"]["human"]
    assert (e["n"], e["n_pos"], e["n_neg"], e["n_not_called"], e["n_pos_not_called"]) == (9, 4, 5, 1, 1)
    assert e["ap"]["model"] == pytest.approx((3 + 4 / 9) / 4)  # positives at ranks 1,2,3,9
    assert e["ap"]["B1"] == pytest.approx(47 / 72)  # positives at ranks 1,3,6,9
    assert sec["end_to_end"]["yeast"]["ap"]["model"] == pytest.approx(5 / 6)  # missed pair is a negative: no change
    assert sec["end_to_end"]["T_vs_B1"] == pytest.approx(np.mean([e["ap"]["model"] - e["ap"]["B1"], 1 / 3]))
    # the primary population still excludes the uncalled pairs (DESIGN 6.2)
    assert toy[0]["primary"]["per_species"]["human"]["n"] == 8
    # tiers (own bootstrap each; a tier whose label column is absent records its error and the others still run)
    t = sec["tiers"]
    assert t["lenient"]["per_species"]["human"]["ap_model"] == pytest.approx(11 / 12)
    assert t["lenient"]["per_species"]["human"]["ap_B1"] == pytest.approx(11 / 15)
    assert t["lenient"]["per_species"]["yeast"]["ap_model"] == pytest.approx(13 / 18)
    assert t["lenient"]["per_species"]["yeast"]["ap_B1"] == pytest.approx(0.5)
    assert t["two_studies"]["per_species"]["human"]["ap_B1"] == pytest.approx(5 / 6)
    assert t["strict"]["per_species"]["yeast"]["ap_model"] == pytest.approx(1.0)
    assert t["strict"]["per_species"]["yeast"]["ap_B1"] == pytest.approx(0.5)
    assert "error" in t["biogrid"]
    # LOSO: yeast-only on human (above) and human-only on yeast: stand-in ranks Q7..Q0, positives at ranks 5 and 7
    h2y = sec["human_only_model_on_yeast"]
    assert h2y["ap"]["model_human_only"] == pytest.approx((1 / 5 + 2 / 7) / 2)
    assert h2y["delta_vs_B1"]["model_human_only"]["point"] == pytest.approx((1 / 5 + 2 / 7) / 2 - 0.5)
    assert (h2y["n"], h2y["n_pos"]) == (8, 2)


def test_toy_d8_d10_extras_known_answers(toy):
    from sklearn.metrics import average_precision_score as aps

    sec = toy[0]["secondary"]
    d8 = sec["human_strata_D8"]
    # B4 (n_testable, diagnostic) reported next to every human AP
    h = pd.read_csv(TOY / "human_test.tsv", sep="\t")
    hc = h.iloc[:8]
    assert d8["by_n_studies"][">=3"]["ap"]["B4"] == pytest.approx(aps(hc.label_primary[:5], hc.n_testable[:5]))
    assert d8["stratum_averaged_ap"]["model"]["mean"] == pytest.approx(1.0)  # strata 2 and >=3 (1 has no positives)
    assert np.isnan(d8["stratum_averaged_ap"]["model"]["1"])
    assert d8["stratum_averaged_ap"]["B1"]["mean"] == pytest.approx((1.0 + 5 / 6) / 2)
    # without Harle/Chan-dependent labels: P0 (a positive) removed; model AP 1, B1 order 0,1,0,0,1,0,0 -> 0.45
    w = d8["without_harle_chan_dependent"]
    assert (w["n"], w["n_pos"]) == (7, 2) and w["ap"]["model"] == pytest.approx(1.0) and w["ap"]["B1"] == pytest.approx(0.45)
    assert d8["harle_chan_dependent_counts"] == {"n": 1, "n_pos": 1}
    # D10: yeast without orientation-conflict positives (Q1 removed): model order Q0,Q3,.. AP 1/2; B1 Q0,Q2,Q3 -> 1/3
    d10 = sec["yeast_and_ribosomal_D10"]
    o = d10["without_orientation_conflict_positives"]
    assert (o["n"], o["n_pos"]) == (7, 1) and o["ap"]["model"] == pytest.approx(0.5) and o["ap"]["B1"] == pytest.approx(1 / 3)
    assert d10["orientation_conflict_positives_removed"] == 1
    assert d10["ribosomal_KW-0689"]["human"]["status"].startswith("unavailable")
    # D4 (unseen): Arabidopsis without the positive whose pair has a lethal single mutant (A0)
    s = sec["unseen_species"]["arabidopsis"]["strict"]["without_lethal_single_positives"]
    assert (s["n"], s["n_pos"], s["n_positives_removed"]) == (4, 1, 1) and s["ap"]["B1"] == pytest.approx(1 / 3)


def test_toy_leakage_tertiles_match_an_independent_computation(toy):
    from sklearn.metrics import average_precision_score as aps

    lk = toy[0]["secondary"]["leakage_tertiles"]
    h = pd.read_csv(TOY / "human_test.tsv", sep="\t")
    y = pd.read_csv(TOY / "yeast_test.tsv", sep="\t")
    allv = np.concatenate([h.loc[h.label_primary.isin([0, 1]), "leakage_bsr"], y.loc[y.primary.isin(["pos", "neg"]), "leakage_bsr"]])
    cuts = np.quantile(allv[allv > 0], [1 / 3, 2 / 3])
    assert lk["cut_points"] == pytest.approx(list(cuts))
    fh = pd.read_parquet(TOY / "features/human/tierA.parquet")
    fy = pd.read_parquet(TOY / "features/yeast/tierA.parquet")
    parts = {}
    for sp, lab, feat, lcol in (("human", h.iloc[:8], fh, "label_primary"), ("yeast", y.iloc[:8], fy, "primary")):
        d = lab.reset_index(drop=True).assign(model=feat["bits_pair"].to_numpy(), B1=feat["ident_local_max"].to_numpy())
        d["y"] = (d[lcol] == 1) if sp == "human" else (d[lcol] == "pos")
        parts[sp] = d
    def tert(v):
        return "no_hit" if v <= 0 else "T1" if v <= cuts[0] else "T2" if v <= cuts[1] else "T3"
    T = {}
    for k in ("T1", "T2", "T3"):
        deltas, n, npos = [], 0, 0
        for sp, d in parts.items():
            m = d["leakage_bsr"].map(tert) == k
            n, npos = n + int(m.sum()), npos + int(d.y[m].sum())
            if d.y[m].sum() > 0:
                deltas.append(aps(d.y[m], d.model[m]) - aps(d.y[m], d.B1[m]))
        assert (lk["per_tertile"][k]["n"], lk["per_tertile"][k]["n_pos"]) == (n, npos)
        T[k] = np.mean(deltas)
        assert lk["per_tertile"][k]["T_species_mean"]["point"] == pytest.approx(T[k])
    assert lk["trend_slope_per_tertile"]["point"] == pytest.approx((T["T3"] - T["T1"]) / 2)
    assert lk["per_tertile"]["no_hit"]["n"] == sum(int((d["leakage_bsr"] <= 0).sum()) for d in parts.values())


def test_every_prereg_secondary_is_present_in_the_toy_run(toy):
    res = toy[0]
    assert res["status"] == "complete" and res["secondary_errors"] == {}
    assert set(res["secondary"]) == set(ft.SECONDARY_NAMES)


def test_shared_dev_functions_known_answers():
    from slparalog.eval import controls

    y = np.array([1, 0, 1, 0, 1, 0])
    s = np.array([.9, .8, .7, .1, .6, .5])
    nst = np.array([1, 1, 2, 2, 3, 3])
    out = controls.strat_avg_ap(y, s, nst, np.ones(6, bool))
    assert out["1"] == pytest.approx(1.0) and out["2"] == pytest.approx(1.0) and out[">=3"] == pytest.approx(1.0)
    lab = pd.DataFrame({"species": "human", "pair_id": list("abcd"), "label": [1, 0, 1, 0], "called": [True, True, True, False]})
    sc = pd.DataFrame({"m": [0.9, 0.1, 0.5]}, index=["a", "b", "c"])
    e = controls.end_to_end_ap(lab, sc, ["m"], ["human"])
    assert e["human"]["ap"]["m"] == pytest.approx((1 + 2 / 2) / 2)  # d (uncalled negative) sits last
    assert e["human"]["n_not_called"] == 1


def test_toy_is_deterministic_and_writes_only_its_out_dir(toy, tmp_path):
    res, out = toy
    res2 = ft.run_toy(tmp_path / "again")
    strip = lambda r: {k: v for k, v in r.items() if k != "timestamp_utc"}  # noqa: E731
    assert json.dumps(strip(res), sort_keys=True, default=str) == json.dumps(strip(res2), sort_keys=True, default=str)
    assert {p.name for p in out.iterdir()} == {"final_test.json", "final_test.md"}
    assert "CLEARLY BEATS" in (out / "final_test.md").read_text()
    assert json.loads((out / "final_test.json").read_text())["purpose"] == "dry-run-toy"


def test_toy_never_touches_sealed_paths_or_the_log(tmp_path, monkeypatch):
    def boom(*a, **k):
        raise AssertionError("toy mode called the sealed reader")
    monkeypatch.setattr(ft, "read_sealed", boom)
    seen = []
    real_read_csv = pd.read_csv

    def spy(path, *a, **k):
        seen.append(str(path))
        return real_read_csv(path, *a, **k)

    monkeypatch.setattr(pd, "read_csv", spy)
    log_before = (ROOT / "results" / "final").exists()
    ft.run_toy(tmp_path / "t")
    assert seen and all("sealed" not in s.lower() and "final_test_toy" in s for s in seen)
    assert (ROOT / "results" / "final").exists() == log_before
    with pytest.raises(ft.Refusal, match="sealed"):
        ft.ToyBackend(ROOT / "data" / "labels" / "sealed")
    with pytest.raises(ft.Refusal, match="results/final"):
        ft.run_toy(ROOT / "results" / "final")
    with pytest.raises(ft.Refusal, match="sealed"):
        ft.read_fixture(ROOT / "data" / "labels" / "sealed" / "x.tsv")


def test_toy_cli_runs_as_a_script(tmp_path):
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "eval" / "final_test.py"), "--dry-run-toy",
                        "--out-dir", str(tmp_path / "o")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "CLEARLY BEATS" in r.stdout and (tmp_path / "o" / "final_test.json").exists()


def test_relative_gain_threshold_comes_from_config(toy):
    succ = dict(ft.load_yaml(ROOT / "config" / "success.yaml"), relative_gain=0.9)
    p = ft.compute(ft.ToyBackend(), succ)["primary"]
    assert p["pass"] and not p["clearly_beats_identity"]
    assert p["verdict"].startswith("BEATS identity (B1) but not")


def test_identical_scores_give_zero_delta_and_not_demonstrated():
    from slparalog.eval.bootstrap import make_bootstrap

    rng = np.random.default_rng(0)
    rows = []
    for sp in ("human", "yeast"):
        for i in range(40):
            rows.append({"species": sp, "group": f"{sp}{i}", "label": int(i % 5 == 0), "pair_id": f"{sp}:{i}",
                         "model": rng.random()})
    S = pd.DataFrame(rows)
    S["B1"] = S["model"]
    S["B2"] = S["model"]
    succ = ft.load_yaml(ROOT / "config" / "success.yaml")
    bs = make_bootstrap(S, 500, succ["bootstrap_seed"])
    p = ft.primary_block(S, bs, succ, ["human", "yeast"], {})
    assert p["T"]["point"] == 0.0 and (p["T"]["lo"], p["T"]["hi"]) == (0.0, 0.0)
    assert not p["pass"] and not p["clearly_beats_identity"] and p["verdict"].startswith("NOT DEMONSTRATED")
    assert p["relative_gain_mean"] == 0.0


def test_gene_components():
    a = pd.Series(["g1", "g2", "g4", "g6"])
    b = pd.Series(["g2", "g3", "g5", "g1"])
    grp = ft.gene_components(a, b)
    assert grp[0] == grp[1] == grp[3] and grp[2] != grp[0] and len(set(grp)) == 2
