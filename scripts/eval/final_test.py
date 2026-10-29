#!/usr/bin/env python
"""The single test-set access (DESIGN v3 section 6.2-6.3; RESEARCH_TREE D8, D12, D17; gate G15 in
reviews/round2/agent_failure_modes.md). Run once, by the PI, after the model is frozen and tagged.

    .venv/bin/python scripts/eval/final_test.py --record-tests     # full pytest at HEAD, writes the marker
    .venv/bin/python scripts/eval/final_test.py                    # the access
    .venv/bin/python scripts/eval/final_test.py --dry-run-toy      # whole pipeline on tests/fixtures; no sealed path

Refuses to run unless ALL of these hold (all checked before any sealed file is touched):
  1. clean git tree (untracked files included);
  2. HEAD is the commit of tag ``final-model-v1``, or a descendant whose diff to the tag touches only
     ``results/`` (and not ``results/frozen/``);
  3. the model artifact (``results/frozen/model``) is tracked, its model.joblib sha256 equals model.json, the
     training label / feature checksums in model.json equal the files on disk, model.json's git commit is an
     ancestor of HEAD and was taken from a clean tree, and model.json IS the frozen spec (R1): GBDT, feature
     set A+B, the hyper-parameters of config/models.yaml anchoring.ensemble.gbdt_params (GBDT_AB_17), the
     registry's A+B model features minus models.yaml excluded_features, tier primary, config/seeds.yaml's first
     seed, and the dev-selected B1 variant (models.yaml b1.selected);
  4. ``config/protected.sha256`` verifies (hash equals, or a D<n> note in RESEARCH_TREE.md cites the new hash)
     and itself lists this script and its guard test;
  5. the full pytest suite passed at HEAD: marker ``<git-dir>/final_test_pytest_ok.json`` (written by
     ``--record-tests`` or ``--run-tests``) names HEAD;
  6. ``results/final/access_log.tsv`` holds no ``access_started`` record (single access) and
     ``results/final/final_test.json`` does not exist;
  7. every sealed file (names from config) is tracked and was first committed at or before the tag, and the
     success-criterion config is the pre-registered one;
  7b. (R3) every sealed file is HASH-verified against the manifest (bytes hashed, nothing parsed or kept) before
     the access record, so a corrupted file refuses without consuming the access.
  8. PREFLIGHT (still before any sealed file or access record): every called pair of every species evaluated
     (human, yeast, arabidopsis, worm; rows of the feature tables, no labels) is scored by the frozen model and by
     the yeast-only refit; B1 and B2 are computed too. Refuses, without consuming the access, unless the feature
     tables carry all model features and the registry/schema hashes recorded in model.json, and all scores are
     finite. NaN policy: model features may be NaN only as the estimators handle them (as in training); B1 inputs
     and the B2 inputs (family size, closest-pair flag) must be finite, there is no imputation for them.
Then it writes an ``access_started`` record (fsync'd) and only afterwards reads sealed files, through the one
function ``read_sealed`` (checks the access record, the manifest and each file's sha256). It then appends
``access_completed`` (or ``access_failed``) and writes results/final/final_test.json and final_test.md. The log
and the results are left uncommitted: the next access of any kind is refused until they are committed (clean-tree
check), and re-runs need a reviewed bug (DESIGN 6.3) and a PI decision, so only purpose=primary is supported.

Computations (metrics, bootstrap, baselines are reused from slparalog.eval; nothing is re-implemented):
  Primary  T = mean over {human, yeast} of AP(model) - AP(B1) on the frozen test pairs (labelled, primary tier,
           called by the caller; B1 and the model on exactly the same pairs); percentile 95% CI from the
           stratified group bootstrap (resamples, seed, strata from config/success.yaml); per-species dAP > 0;
           relative gain (mean over species of AP_model/AP_B1 - 1) on the point estimate vs relative_gain.
           pass = CI of T excludes 0 and dAP_s > 0 in both species; "clearly beats identity" = pass and
           relative gain >= threshold.
  Secondary (D17c, DESIGN 6.2): dAP vs B2; human dAP vs B1 within the >=3-studies stratum and within exact
           study-set cells; the yeast-only-trained model (refit on dev yeast, same hyper-parameters) on human
           test pairs; unseen species (Arabidopsis, worm): AP and dAP vs B1 with bootstrap CI, features from
           data/interim/features/<species>/ (descriptive, pre-declared underpowered, no pass/fail).
  Also (R2, evaluated on the same frozen model, dev-harness functions): end-to-end (caller-missed labelled pairs at
  the minimum score; the primary population still excludes them), lenient / strict / >=2-studies / BioGRID tiers,
  leakage tertiles with a trend CI, LOSO in both directions, D8 (B4, strata, stratum-averaged AP, without
  Harle/Chan-dependent), D10 (yeast without orientation-conflict positives, ribosomal), D4 (Arabidopsis without
  lethal-single positives), identity tertiles. The primary block is written to results/final first (fsync); each
  secondary runs in its own try/except and a failure is recorded under ``secondary_errors`` (R3).
  The nested-CV fallback of the round-2 review (B5) is NOT implemented: DESIGN v3 6.1 rules it out.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import getpass
import hashlib
import io
import json
import os
import re
import socket
import subprocess
import sys
import tempfile
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from slparalog.config import load_yaml
from slparalog.eval import baselines, controls
from slparalog.eval.bootstrap import GroupBootstrap, boot_metric, ci, make_bootstrap, species_boot, species_mean_delta
from slparalog.eval.controls import subset_compare
from slparalog.eval.dataset import (TIERS, Assembled, file_sha256, join_called, load_dev_labels, load_feature_table,
                                    pair_key)
from slparalog.eval.metrics import ap, ap_pair
from slparalog.eval.stats import species_balanced_weights, summary
from slparalog.io.labels import is_sealed

ROOT = Path(__file__).resolve().parents[2]
TAG = "final-model-v1"
SPECIES = ("human", "yeast")
UNSEEN = ("arabidopsis", "worm")
# unseen-species tiers: (name -> label column); the first is the headline endpoint (config/labels_unseen.yaml)
UNSEEN_TIERS = {"arabidopsis": {"strict": "label_strict", "broad": "label_broad"}, "worm": {"primary": "label"}}
MANIFEST = "MANIFEST.sha256"
LOG_NAME = "access_log.tsv"
RESULT_NAME = "final_test.json"
MARKER_NAME = "final_test_pytest_ok.json"
LOG_COLUMNS = ["timestamp_utc", "event", "purpose", "git_head", "tag_commit", "model_sha256", "config_sha256",
               "sealed_manifest_sha256", "user", "host", "note"]
# The frozen model (RESEARCH_TREE D17a: incumbent GBDT A+B = GBDT_AB_17). Hyper-parameters, tier, seed, excluded
# features and the feature universe are read from config/ and the registry, not written here.
FROZEN_FEATURE_SET = "A+B"
FROZEN_MODEL_TYPE = "gbdt"
FROZEN_TIER = "primary"
REQUIRED_PROTECTED = ("scripts/eval/final_test.py", "tests/test_final_test_guard.py", "config/success.yaml",
                      "src/slparalog/eval/metrics.py", "tests/test_data_access.py")
EXPECTED_STRATA = ["species_composition", "group_has_positive"]
SMALL_N_POS = 10  # CIs with fewer positives are flagged, not suppressed (pre-declared secondaries)


class Refusal(RuntimeError):
    """A precondition failed (or the sealed reader was called out of order): nothing was read."""

    def __init__(self, failures: list[str] | str):
        self.failures = [failures] if isinstance(failures, str) else list(failures)
        super().__init__("; ".join(self.failures))


# ====================================================================================== context + git helpers

@dataclass
class Ctx:
    root: Path
    model_rel: str = "results/frozen/model"
    tag: str = TAG
    access_started: bool = False
    manifest: dict | None = None
    sealed_read: dict = field(default_factory=dict)  # name -> sha256 of every sealed file read

    @property
    def config_dir(self) -> Path: return self.root / "config"
    @property
    def out_dir(self) -> Path: return self.root / "results" / "final"
    @property
    def log_path(self) -> Path: return self.out_dir / LOG_NAME
    @property
    def sealed_dir(self) -> Path: return self.root / "data" / "labels" / "sealed"
    @property
    def dev_dir(self) -> Path: return self.root / "data" / "labels" / "dev"
    @property
    def feature_dir(self) -> Path: return self.root / "data" / "interim" / "features"
    @property
    def model_dir(self) -> Path: return self.root / self.model_rel


def git(ctx: Ctx, *args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=ctx.root, capture_output=True, text=True, check=check)


def git_dir(ctx: Ctx) -> Path:
    p = Path(git(ctx, "rev-parse", "--git-dir").stdout.strip())
    return p if p.is_absolute() else ctx.root / p


def head_sha(ctx: Ctx) -> str:
    return git(ctx, "rev-parse", "HEAD").stdout.strip()


def tag_commit(ctx: Ctx) -> str | None:
    r = git(ctx, "rev-parse", "--verify", "--quiet", f"refs/tags/{ctx.tag}^{{commit}}", check=False)
    return r.stdout.strip() if r.returncode == 0 and r.stdout.strip() else None


def is_ancestor(ctx: Ctx, a: str, b: str) -> bool:
    return git(ctx, "merge-base", "--is-ancestor", a, b, check=False).returncode == 0


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def config_sha256(ctx: Ctx) -> str:
    h = hashlib.sha256()
    for n in ("success.yaml", "models.yaml"):
        h.update((ctx.config_dir / n).read_bytes())
    return h.hexdigest()


# ====================================================================================== preconditions
# Each check returns a list of failure strings (empty = ok). None of them opens a sealed file.

def check_clean_tree(ctx: Ctx) -> list[str]:
    st = git(ctx, "status", "--porcelain").stdout
    if st.strip():
        return [f"git tree is not clean ({len(st.splitlines())} changed/untracked paths); commit everything first"]
    return []


def check_tag(ctx: Ctx) -> list[str]:
    tc = tag_commit(ctx)
    if tc is None:
        return [f"tag {ctx.tag} does not exist"]
    head = head_sha(ctx)
    if head == tc:
        return []
    if not is_ancestor(ctx, tc, head):
        return [f"HEAD {head[:12]} is not tag {ctx.tag} ({tc[:12]}) nor a descendant of it"]
    changed = git(ctx, "diff", "--name-only", tc, head).stdout.split()
    bad = [p for p in changed if not p.startswith("results/") or p.startswith("results/frozen/")]
    if bad:
        return [f"HEAD is a descendant of {ctx.tag} but its diff touches non-results paths: {bad[:5]}"]
    return []


def frozen_spec_failures(ctx: Ctx, meta: dict) -> list[str]:
    """R1: model.json must BE the frozen model: GBDT on feature set A+B with the hyper-parameters of
    config/models.yaml anchoring.ensemble.gbdt_params (= GBDT_AB_17), the registry's A+B model features minus
    models.yaml excluded_features, tier primary, the first seed of config/seeds.yaml, and the dev-selected B1."""
    from slparalog.features import registry

    out: list[str] = []
    mcfg = load_yaml(ctx.config_dir / "models.yaml")
    seeds = load_yaml(ctx.config_dir / "seeds.yaml")["seeds"]
    excluded = set(mcfg.get("excluded_features") or {})
    want_feats = [f for f in registry.MODEL_FEATURES if registry.BY_NAME[f].tier in ("A", "B") and f not in excluded]
    want = {"feature_set": FROZEN_FEATURE_SET, "model_type": FROZEN_MODEL_TYPE,
            "params": dict(mcfg["anchoring"]["ensemble"]["gbdt_params"]),
            "features": want_feats, "b1_variant": mcfg["b1"]["selected"]}
    for k, w in want.items():
        g_ = meta.get(k)
        if g_ != w:
            out.append(f"model.json {k} {g_!r:.200} != the frozen spec {w!r:.200}")
    if (meta.get("training") or {}).get("tier") != FROZEN_TIER:
        out.append(f"model.json training.tier {(meta.get('training') or {}).get('tier')!r} != {FROZEN_TIER!r}")
    sd = meta.get("seeds") or {}
    if sd.get("seed_used") != seeds[0] or sd.get("seed_list") != seeds:
        out.append(f"model.json seeds {sd} != config/seeds.yaml (seed_used {seeds[0]}, list {seeds})")
    return out


def check_model(ctx: Ctx) -> list[str]:
    d = ctx.model_dir
    mj, mp = d / "model.json", d / "model.joblib"
    if not mj.exists() or not mp.exists():
        return [f"model artifact incomplete: {d} needs model.json and model.joblib"]
    out: list[str] = []
    try:
        meta = json.loads(mj.read_text())
    except ValueError as e:
        return [f"{mj}: invalid JSON ({e})"]
    for f in (mj, mp):
        rel = f.relative_to(ctx.root).as_posix()
        if git(ctx, "ls-files", "--error-unmatch", "--", rel, check=False).returncode != 0:
            out.append(f"{rel} is not tracked by git")
    if sha256_bytes(mp.read_bytes()) != meta.get("model_sha256"):
        out.append("model.joblib sha256 differs from model.json model_sha256")
    g = meta.get("git") or {}
    commit = g.get("commit")
    if not commit or commit == "unknown":
        out.append("model.json has no git commit")
    elif git(ctx, "cat-file", "-e", f"{commit}^{{commit}}", check=False).returncode != 0:
        out.append(f"model.json git commit {commit[:12]} does not exist in this repository")
    elif not is_ancestor(ctx, commit, head_sha(ctx)):
        out.append(f"model.json git commit {commit[:12]} is not an ancestor of HEAD")
    if g.get("dirty"):
        out.append("model.json records that the model was trained from a dirty git tree")
    v = meta.get("b1_variant")
    if v not in baselines.B1_VARIANTS:  # (the value itself is pinned in frozen_spec_failures)
        out.append(f"model.json lacks a valid frozen b1_variant (got {v!r}; one of {sorted(baselines.B1_VARIANTS)})")
    for k in ("features", "feature_set", "model_type", "params", "seeds", "training"):
        if k not in meta:
            out.append(f"model.json lacks {k!r}")
    out += frozen_spec_failures(ctx, meta)
    tr = meta.get("training") or {}
    if set(tr.get("label_files", {})) != set(SPECIES):
        out.append("model.json training.label_files must cover human and yeast")
    for sp, rec in (tr.get("label_files") or {}).items():
        f = ctx.dev_dir / Path(rec["path"]).name
        if not f.exists() or file_sha256(f) != rec["sha256"]:
            out.append(f"training label file for {sp} differs from model.json ({f.name})")
    for sp, tables in (tr.get("feature_files") or {}).items():
        for name, rec in tables.items():
            f = ctx.feature_dir / sp / Path(rec["path"]).name
            if not f.exists() or file_sha256(f) != rec["sha256"]:
                out.append(f"training feature table {sp}/{name} differs from model.json ({f.name})")
    return out


def _protected_entries(ctx: Ctx) -> dict[str, str]:
    out = {}
    for line in (ctx.config_dir / "protected.sha256").read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            digest, rel = line.split(None, 1)
            out[rel.strip()] = digest
    return out


def check_protected(ctx: Ctx) -> list[str]:
    pf = ctx.config_dir / "protected.sha256"
    if not pf.exists():
        return ["config/protected.sha256 is missing"]
    entries = _protected_entries(ctx)
    out = [f"protected list lacks {r}" for r in REQUIRED_PROTECTED if r not in entries]
    tree = ctx.root / "RESEARCH_TREE.md"
    notes = ([ln for ln in tree.read_text().splitlines() if re.search(r"\bD\d+\b", ln)] if tree.exists() else [])
    sealed_rel = ctx.sealed_dir.relative_to(ctx.root).as_posix() + "/"
    for rel, recorded in entries.items():
        if rel.startswith(sealed_rel):  # hashed inside read_sealed, after the access record
            continue
        f = ctx.root / rel
        if not f.is_file():
            out.append(f"protected file missing: {rel}")
            continue
        cur = file_sha256(f)
        if cur != recorded and not any(cur[:12] in n for n in notes):
            out.append(f"protected file changed without a D<n> note: {rel}")
    return out


def marker_path(ctx: Ctx) -> Path:
    return git_dir(ctx) / MARKER_NAME


def check_tests_marker(ctx: Ctx) -> list[str]:
    mp = marker_path(ctx)
    if not mp.exists():
        return [f"no pytest marker ({mp}); run with --record-tests (or --run-tests) first"]
    try:
        m = json.loads(mp.read_text())
    except ValueError:
        return [f"{mp}: invalid JSON"]
    head = head_sha(ctx)
    if not m.get("passed") or m.get("head") != head:
        return [f"pytest marker is not a pass at HEAD {head[:12]} (marker head {str(m.get('head'))[:12]})"]
    return []


def read_log(ctx: Ctx) -> list[dict]:
    if not ctx.log_path.exists():
        return []
    with open(ctx.log_path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check_no_prior_access(ctx: Ctx) -> list[str]:
    out = []
    started = [r for r in read_log(ctx) if r.get("event") in ("access_started", "access_completed")]
    if started:
        out.append(f"{ctx.log_path.relative_to(ctx.root)} already records an access ({started[0].get('event')} "
                   f"{started[0].get('timestamp_utc')}): single access, refused")
    if (ctx.out_dir / RESULT_NAME).exists():
        out.append(f"results/final/{RESULT_NAME} already exists")
    return out


def sealed_names(ctx: Ctx) -> list[str]:
    files = load_yaml(ctx.config_dir / "labels_unseen.yaml")["outputs"]["files"]
    return [MANIFEST] + [f"{sp}_test.tsv" for sp in SPECIES] + [files[sp] for sp in UNSEEN]


def check_sealed_committed(ctx: Ctx) -> list[str]:
    """Each sealed file is tracked and first committed at or before the tag (git metadata only)."""
    tc = tag_commit(ctx)
    if tc is None:
        return []  # reported by check_tag
    out = []
    for n in sealed_names(ctx):
        rel = ctx.sealed_dir.relative_to(ctx.root).as_posix() + "/" + n
        r = git(ctx, "log", "--diff-filter=A", "--format=%H", "--", rel).stdout.split()
        if not r:
            out.append(f"{rel} is not committed")
        elif not is_ancestor(ctx, r[-1], tc):
            out.append(f"{rel} was first committed after tag {ctx.tag}")
    return out


def check_config(ctx: Ctx) -> list[str]:
    s = load_yaml(ctx.config_dir / "success.yaml")
    out = []
    if s.get("bootstrap_strata") != EXPECTED_STRATA:
        out.append(f"success.yaml bootstrap_strata {s.get('bootstrap_strata')} != {EXPECTED_STRATA}")
    if s.get("ci_type") != "percentile":
        out.append("success.yaml ci_type is not 'percentile'")
    if s.get("relative_gain_on") != "point_estimate":
        out.append("success.yaml relative_gain_on is not 'point_estimate'")
    for k in ("relative_gain", "bootstrap_resamples", "bootstrap_seed", "ci_level"):
        if k not in s:
            out.append(f"success.yaml lacks {k}")
    return out


CHECKS = (check_clean_tree, check_tag, check_model, check_protected, check_tests_marker, check_no_prior_access,
          check_sealed_committed, check_config)


def run_preconditions(ctx: Ctx) -> list[str]:
    failures: list[str] = []
    for chk in CHECKS:
        try:
            failures += chk(ctx)
        except Exception as e:  # a check that cannot run is a failed check
            failures.append(f"{chk.__name__} could not run: {type(e).__name__}: {e}")
    return failures


def record_tests(ctx: Ctx) -> None:
    """Run the FULL pytest suite at HEAD from a clean tree; on success write the marker (inside the git dir)."""
    pre = check_clean_tree(ctx)
    if pre:
        raise Refusal(pre)
    head = head_sha(ctx)
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"], cwd=ctx.root)
    if r.returncode != 0:
        raise Refusal(f"pytest failed at {head[:12]} (exit {r.returncode}); no marker written")
    if head_sha(ctx) != head or check_clean_tree(ctx):
        raise Refusal("the tree changed while pytest ran; no marker written")
    marker_path(ctx).write_text(json.dumps({"head": head, "passed": True, "command": "python -m pytest -q",
                                            "timestamp_utc": _now()}) + "\n")


# ====================================================================================== access log + sealed reader

def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log_append(ctx: Ctx, event: str, note: str = "", model_sha: str = "") -> None:
    ctx.out_dir.mkdir(parents=True, exist_ok=True)
    new = not ctx.log_path.exists()
    row = {"timestamp_utc": _now(), "event": event, "purpose": "primary", "git_head": head_sha(ctx),
           "tag_commit": tag_commit(ctx) or "", "model_sha256": model_sha, "config_sha256": config_sha256(ctx),
           "sealed_manifest_sha256": _protected_entries(ctx).get(
               ctx.sealed_dir.relative_to(ctx.root).as_posix() + "/" + MANIFEST, ""),
           "user": getpass.getuser(), "host": socket.gethostname(), "note": note.replace("\t", " ").replace("\n", " ")}
    with open(ctx.log_path, "a", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        if new:
            w.writerow(LOG_COLUMNS)
        w.writerow([row[c] for c in LOG_COLUMNS])
        fh.flush()
        os.fsync(fh.fileno())


def log_refusal(ctx: Ctx, failures: list[str]) -> None:
    """Refusals go to an untracked file inside the git dir, so they never dirty the tree."""
    try:
        with open(git_dir(ctx) / "final_test_refusals.log", "a") as fh:
            fh.write(f"{_now()}\t" + " | ".join(failures) + "\n")
    except OSError:
        pass


def read_sealed(ctx: Ctx, name: str, verify_only: bool = False) -> pd.DataFrame | dict | str:
    """THE sealed reader: the only function in this script that opens a file under data/labels/sealed/.

    Parsing refuses unless the ``access_started`` record has been written. Verifies the manifest against the
    protected hash, and every table against the manifest, before parsing it. ``MANIFEST.sha256`` returns
    {name: sha256}. ``verify_only=True`` (allowed before the record, so that a corrupted file refuses without
    consuming the access) only hashes the bytes and returns the digest: nothing is parsed or kept.
    """
    if not ctx.access_started and not verify_only:
        raise Refusal("sealed reader called before the access_started record was written")
    sealed = ctx.sealed_dir.resolve()
    path = (sealed / name).resolve()
    if path.parent != sealed:
        raise Refusal(f"{name}: not a file directly inside the sealed directory")
    if ctx.manifest is None:
        raw = (sealed / MANIFEST).read_bytes()
        rec = _protected_entries(ctx).get(ctx.sealed_dir.relative_to(ctx.root).as_posix() + "/" + MANIFEST)
        if rec is None or sha256_bytes(raw) != rec:
            raise Refusal("sealed MANIFEST.sha256 does not match config/protected.sha256")
        ctx.manifest = {}
        for line in raw.decode().splitlines():
            if line.strip():
                digest, n = line.split(None, 1)
                ctx.manifest[n.strip().lstrip("*")] = digest
    if name == MANIFEST:
        return dict(ctx.manifest)
    if name not in ctx.manifest:
        raise Refusal(f"{name} is not listed in the sealed manifest")
    data = path.read_bytes()
    digest = sha256_bytes(data)
    if digest != ctx.manifest[name]:
        raise Refusal(f"{name} sha256 does not match the sealed manifest")
    if verify_only:
        return digest
    ctx.sealed_read[name] = digest
    return pd.read_csv(io.BytesIO(data), sep="\t", low_memory=False)


def verify_sealed(ctx: Ctx) -> list[str]:
    """R3: hash-verify every sealed file (names from config) against the manifest, before the access record."""
    out = []
    for n in sealed_names(ctx):
        try:
            read_sealed(ctx, n, verify_only=True)
        except (Refusal, OSError) as e:
            out.append(f"sealed file {n}: {e}")
    return out


def read_fixture(path: Path) -> pd.DataFrame:
    """Toy-mode label reader: fixture files only; refuses any sealed path."""
    if is_sealed(path, ROOT / "data" / "labels" / "sealed"):
        raise Refusal(f"toy mode must never touch sealed paths: {path}")
    return pd.read_csv(path, sep="\t", low_memory=False)


# ====================================================================================== population builders

def gene_components(a: pd.Series, b: pd.Series) -> np.ndarray:
    """Group id per pair = connected component of the gene graph (pairs sharing a gene share a group)."""
    parent: dict[str, str] = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for x, y in zip(a, b):
        parent[find(x)] = find(y)
    roots = [find(x) for x in a]
    ids = {r: f"c{i}" for i, r in enumerate(sorted(set(roots)))}
    return np.array([ids[r] for r in roots])


def unseen_labels(raw: pd.DataFrame, species: str, col: str) -> pd.DataFrame:
    """Unseen-species labels resolved to {0,1}; NA / other values dropped; groups = gene-graph components."""
    lab = pd.to_numeric(raw[col], errors="coerce")
    keep = lab.isin([0, 1])
    d = raw.loc[keep]
    lo, hi = pair_key(d["gene_a"], d["gene_b"])
    out = pd.DataFrame({"species": species, "gene_a": lo.to_numpy(), "gene_b": hi.to_numpy()})
    out["pair_id"] = species + ":" + out["gene_a"] + "|" + out["gene_b"]
    out["label"] = lab[keep].astype(int).to_numpy()
    if "single_lethal_n" in raw.columns:
        out["n_lethal_single"] = pd.to_numeric(raw.loc[keep, "single_lethal_n"], errors="coerce").to_numpy()
    out["group"] = gene_components(out["gene_a"], out["gene_b"])
    if out["pair_id"].duplicated().any():
        raise ValueError(f"{species}: duplicate unordered pairs")
    return out.reset_index(drop=True)


def b1_scores(df: pd.DataFrame, variant: str) -> np.ndarray:
    """Scores of the frozen B1 variant only (the registered definition in slparalog.eval.baselines)."""
    return np.asarray(baselines.B1_VARIANTS[variant][1](df), dtype=float)


def feature_paths(feature_dir: Path, species: str, names: list[str]) -> dict[str, Path]:
    out = {}
    for n in names:
        p = feature_dir / species / f"{n}.parquet"
        if not p.exists():
            raise FileNotFoundError(f"feature table missing: {p}")
        out[n] = p
    return out


def population(lab: pd.DataFrame, species: str, tier: str, paths: dict[str, Path]) -> dict:
    """labelled -> called join; counts via the harness' Assembled.counts."""
    lab, prim, cols, prov = join_called(lab, paths)
    A = Assembled(tier=tier, labelled=lab, primary=prim, feature_columns=cols)
    return {"labelled": lab, "primary": prim.reset_index(drop=True), "counts": A.counts().get(species, {}),
            "feature_files": prov}


# ====================================================================================== backends

class RealBackend:
    """The real model, dev data and sealed files."""

    toy = False

    def __init__(self, ctx: Ctx):
        from slparalog.model import load_model

        self.ctx = ctx
        self.model = load_model(ctx.model_dir)
        self.meta = self.model.meta
        self.b1_variant = self.meta["b1_variant"]
        tabs = {sp: list(self.meta["training"]["feature_files"][sp]) for sp in SPECIES}
        if tabs["human"] != tabs["yeast"]:
            raise ValueError(f"human and yeast training tables differ: {tabs}")
        self.table_names = tabs["human"]
        self.mcfg = load_yaml(ctx.config_dir / "models.yaml")
        self.features = list(self.meta["features"])
        b2c = self.mcfg["b2"]
        self.b2_cols = [b2c["family_size"], b2c["closest_pair"]]
        self._b2 = None
        self._only: dict = {}

    def b2_scorer(self):
        if self._b2 is None:
            self._b2 = self.fit_b2()
        return self._b2

    def species_scorer(self, sp: str):
        """Frozen model type / feature set / hyper-parameters refit on development pairs of one species only."""
        if sp not in self._only:
            self._only[sp] = self.fit_species_only(sp)
        return self._only[sp]

    def yeast_scorer(self):
        return self.species_scorer("yeast")

    def human_scorer(self):
        return self.species_scorer("human")

    def ribosomal_genes(self, sp: str):
        from slparalog.eval import annotations

        return annotations.keyword_genes(sp, self.mcfg["controls"]["ribosomal_keyword"])

    def table_failures(self, sp: str) -> list[str]:
        """Feature tables are validated by the column definitions of their own columns (``columns_hash``, read
        with slparalog.model.table_meta, i.e. embedded metadata plus the restamp sidecar, and compared with
        registry.table_columns_hash). The global registry hash is provenance only. For the training species the
        record in model.json (``columns_hash`` or, for older artifacts, ``schema_hash``) must agree too."""
        import pyarrow.parquet as pq

        from slparalog.features import registry
        from slparalog.model import ModelError, check_feature_table, table_meta

        out = []
        rec = (self.meta["training"].get("feature_tables") or {})
        for name, p in self.paths(sp).items():
            try:
                check_feature_table(p)  # raises ModelError on a changed column definition
            except ModelError as e:
                out.append(f"{sp}/{name}: {e}")
                continue
            m = table_meta(p) or {}
            cur = registry.table_columns_hash(pq.read_schema(p).names)
            for k, v in rec.items():
                if sp not in SPECIES or Path(k).name != p.name or Path(k).parent.name != sp:
                    continue
                if v.get("columns_hash") is not None:
                    if v["columns_hash"] != cur:
                        out.append(f"{sp}/{name}: columns_hash differs from model.json")
                elif v.get("schema_hash") and m.get("schema_hash") != v["schema_hash"]:
                    out.append(f"{sp}/{name}: schema hash differs from model.json")
        return out

    def read_test_labels(self, sp: str) -> pd.DataFrame:
        return read_sealed(self.ctx, f"{sp}_test.tsv")

    def read_unseen_labels(self, sp: str) -> pd.DataFrame:
        files = load_yaml(self.ctx.config_dir / "labels_unseen.yaml")["outputs"]["files"]
        return read_sealed(self.ctx, files[sp])

    def paths(self, sp: str) -> dict[str, Path]:
        return feature_paths(self.ctx.feature_dir, sp, self.table_names)

    def dev_assembled(self):
        from slparalog.eval.dataset import assemble

        return assemble({sp: self.paths(sp) for sp in SPECIES}, tier="primary",
                        label_paths={sp: self.ctx.dev_dir / f"{sp}_dev.tsv" for sp in SPECIES})

    def fit_b2(self):
        """B2 = identity (frozen B1 variant) + family size + closest-pair flag; standardised L2 logistic
        regression with species-balanced weights, fitted on all development pairs (config/models.yaml b2)."""
        from slparalog.eval.models import LRModel

        b2c = self.mcfg["b2"]
        seed = self.meta["seeds"]["seed_used"]
        dev = self.dev_assembled().primary
        cols = ["b1_selected", b2c["family_size"], b2c["closest_pair"]]

        def design(d: pd.DataFrame) -> np.ndarray:
            b1 = b1_scores(d, self.b1_variant)
            return np.column_stack([b1, d[b2c["family_size"]].to_numpy(float), d[b2c["closest_pair"]].to_numpy(float)])

        y = dev["label"].to_numpy(int)
        m = LRModel(cols, C=b2c["C"], seed=seed, log_cols=[b2c["family_size"]])
        m.fit(design(dev), y, sample_weight=species_balanced_weights(dev["species"].to_numpy(), y))
        return lambda d: np.asarray(m.predict(design(d)), dtype=float)

    def fit_species_only(self, sp: str):
        from slparalog.model import load_model
        from slparalog.train import train

        grid = self.mcfg[self.meta["model_type"]]["grid"]
        params = {k: v for k, v in self.meta["params"].items() if k in grid}
        with tempfile.TemporaryDirectory(dir=self.ctx.root / "results") as tmp:
            meta = train({sp: list(self.paths(sp).values())}, {sp: self.ctx.dev_dir / f"{sp}_dev.tsv"},
                         self.meta["feature_set"], self.meta["model_type"], Path(tmp) / "m", tier="primary",
                         params=params, seed=self.meta["seeds"]["seed_used"], models_config=self.ctx.config_dir / "models.yaml")
            if meta["params"] != self.meta["params"] or meta["features"] != self.meta["features"]:
                raise ValueError(f"{sp}-only refit does not reproduce the frozen hyper-parameters / feature list")
            m = load_model(Path(tmp) / "m")
        return m.score

    def provenance(self) -> dict:
        m = self.meta
        return {"model_id": m["model_id"], "model_sha256": m["model_sha256"], "feature_set": m["feature_set"],
                "model_type": m["model_type"], "params": m["params"], "n_features": len(m["features"]),
                "model_git_commit": (m.get("git") or {}).get("commit"), "seeds": m["seeds"],
                "b1_variant": self.b1_variant, "tables": self.table_names,
                "sealed_files_sha256": dict(self.ctx.sealed_read)}

    def score(self, df: pd.DataFrame) -> np.ndarray:
        return self.model.score(df)

    def mde(self) -> dict:
        p = self.ctx.root / "results" / "mde.json"
        return json.loads(p.read_text()) if p.exists() else {}


# ====================================================================================== statistics

def cell_weighted(bs: GroupBootstrap, sub: pd.DataFrame, cols: list[str], cell_col: str) -> tuple[dict, dict, dict]:
    """Positive-weighted mean over cells of the within-cell AP (cells with >=1 positive and >=1 negative),
    point (sklearn AP via slparalog.eval.metrics.ap) and per-resample (weighted AP from the same group
    multiplicities; a cell contributes in a resample only when it has positive and negative weight)."""
    cells = sub.groupby(cell_col)["label"].agg(["sum", "size"])
    ok = list(cells[(cells["sum"] > 0) & (cells["sum"] < cells["size"])].index)
    num_p = {c: 0.0 for c in cols}
    den_p = 0.0
    num_b = {c: np.zeros(bs.B) for c in cols}
    den_b = np.zeros(bs.B)
    for cell in ok:
        s = sub[sub[cell_col] == cell]
        y = s["label"].to_numpy(int)
        npos = float(y.sum())
        for c in cols:
            num_p[c] += ap(y, s[c].to_numpy()) * npos
        den_p += npos
        W = bs.counts[:, bs.pair_index(s["group"].to_numpy())].astype(float)
        P, N = W @ y, W @ (1 - y)
        valid = (P > 0) & (N > 0)
        boot = boot_metric(bs, s["group"].to_numpy(), y, {c: s[c].to_numpy() for c in cols}, check=False)
        for c in cols:
            num_b[c] += np.where(valid, np.nan_to_num(boot[c]) * P, 0.0)
        den_b += np.where(valid, P, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        point = {c: (num_p[c] / den_p if den_p else float("nan")) for c in cols}
        bootv = {c: np.where(den_b > 0, num_b[c] / np.where(den_b > 0, den_b, 1.0), np.nan) for c in cols}
    n_pairs = int(cells.loc[ok, "size"].sum()) if ok else 0
    n_pos = int(cells.loc[ok, "sum"].sum()) if ok else 0
    info = {"n_cells_used": len(ok), "n_cells_total": int(len(cells)), "n_pairs": n_pairs, "n_pos": n_pos,
            "n_neg": n_pairs - n_pos, "prevalence": (n_pos / n_pairs) if n_pairs else float("nan"),
            "n_pos_in_excluded_cells": int(cells["sum"].sum() - n_pos)}
    return point, bootv, info


def delta_block(df: pd.DataFrame, S: pd.DataFrame, mask: np.ndarray, model: str, base: str, bs: GroupBootstrap,
                level: float, extra: list[str] | None = None, min_pos: int = 1) -> dict:
    """AP of the methods on a subset and dAP vs ``base`` with bootstrap CI (harness' subset_compare).
    ``min_pos``: fewer positives -> counts and AP only, no CI (dev: controls.min_pos_for_ci; the pre-declared
    D17 analyses use 1 and carry a warning below SMALL_N_POS positives)."""
    meths = list(dict.fromkeys([model, base] + (extra or [])))
    r = subset_compare(df, S, mask, meths, base, bs, level, min_pos_for_ci=min_pos, with_ci=[m for m in meths if m != base])
    if r["n_pos"] < SMALL_N_POS and min_pos < SMALL_N_POS:
        r["warning"] = f"fewer than {SMALL_N_POS} positives: the CI is unreliable"
    return r


def primary_block(S: pd.DataFrame, bs: GroupBootstrap, succ: dict, species: list[str], mde: dict) -> dict:
    level = float(succ["ci_level"])
    thr = float(succ["relative_gain"])
    boot = species_boot(bs, S, ["model", "B1", "B2"], "ap", species=species)
    per = {}
    for sp in species:
        d = S[S["species"] == sp]
        mf = pd.DataFrame({"pair_id": d["pair_id"], "label": d["label"], "score": d["model"]})
        r1 = ap_pair(mf, pd.DataFrame({"pair_id": d["pair_id"], "label": d["label"], "score": d["B1"]}))
        r2 = ap_pair(mf, pd.DataFrame({"pair_id": d["pair_id"], "label": d["label"], "score": d["B2"]}))
        y = d["label"].to_numpy(int)
        per[sp] = {
            "n": r1["n"], "n_pos": r1["n_pos"], "n_neg": r1["n_neg"], "prevalence": r1["prevalence"],
            "ap_model": r1["ap_model"], "ap_B1": r1["ap_baseline"], "ap_B2": r2["ap_baseline"],
            "delta_vs_B1": {"point": r1["delta"], **ci(boot[sp]["model"] - boot[sp]["B1"], level)},
            "delta_vs_B2": {"point": r2["delta"], **ci(boot[sp]["model"] - boot[sp]["B2"], level)},
            "relative_gain_vs_B1": (r1["ap_model"] / r1["ap_baseline"] - 1) if r1["ap_baseline"] > 0 else float("nan"),
            "model_summary": summary(y, d["model"].to_numpy()), "B1_summary": summary(y, d["B1"].to_numpy()),
            "B2_summary": summary(y, d["B2"].to_numpy()),
        }
    T_pt = float(np.mean([per[sp]["delta_vs_B1"]["point"] for sp in species]))
    T_ci = ci(species_mean_delta(boot, "model", "B1"), level)
    T2_pt = float(np.mean([per[sp]["delta_vs_B2"]["point"] for sp in species]))
    T2_ci = ci(species_mean_delta(boot, "model", "B2"), level)
    rel = float(np.mean([per[sp]["relative_gain_vs_B1"] for sp in species]))
    excl0 = bool(np.isfinite(T_ci["lo"]) and T_ci["lo"] > 0)
    each_pos = bool(all(per[sp]["delta_vs_B1"]["point"] > 0 for sp in species))
    rel_ok = bool(np.isfinite(rel) and rel >= thr)
    passed = excl0 and each_pos
    clearly = passed and rel_ok
    verdict = ("CLEARLY BEATS identity (B1)" if clearly else
               "BEATS identity (B1) but not 'clearly' (relative gain below threshold)" if passed else
               "NOT DEMONSTRATED (the primary criterion is not met; this is not evidence of no effect)")
    return {
        "T": {"point": T_pt, **T_ci, "ci_level": level, "ci_type": succ["ci_type"],
              "resamples": bs.B, "seed": bs.seed, "strata": succ["bootstrap_strata"],
              "n_strata_groups": bs.strata["stratum"].value_counts().sort_index().to_dict()},
        "per_species": per,
        "relative_gain_mean": rel, "relative_gain_threshold": thr, "relative_gain_judged_on": succ["relative_gain_on"],
        "criteria": {"ci_of_T_excludes_0": excl0, "delta_ap_positive_in_each_species": each_pos,
                     "relative_gain_ge_threshold": rel_ok},
        "pass": passed, "clearly_beats_identity": clearly, "verdict": verdict,
        "predeclared_underpowered": mde.get("underpowered"),
        "T_vs_B2": {"point": T2_pt, **T2_ci},
    }


def raw_map(raw: pd.DataFrame, sp: str, col: str) -> dict:
    """pair_id -> value of a raw label-table column (for evaluation-only attributes such as studies, leakage)."""
    lo, hi = pair_key(raw["gene_a"], raw["gene_b"])
    return dict(zip(sp + ":" + lo + "|" + hi, raw[col]))


def score_frame(backend, df: pd.DataFrame) -> pd.DataFrame:
    """Model, B1 and B2 scores of the called pairs of ``df`` (all must be finite)."""
    S = df[["species", "group", "label", "pair_id"]].copy()
    S["model"] = backend.score(df)
    S["B1"] = b1_scores(df, backend.b1_variant)
    S["B2"] = np.asarray(backend.b2_scorer()(df), dtype=float)
    for c in ("model", "B1", "B2"):
        if not np.isfinite(S[c].to_numpy()).all():
            raise ValueError(f"non-finite {c} scores (G5: same pairs required)")
    return S


def unseen_block(backend, sp: str, succ: dict, mcfg: dict) -> dict:
    out: dict = {"note": "pre-declared underpowered and descriptive; no pass/fail (DESIGN 6.2)"}
    level = float(succ["ci_level"])
    B = int(succ["bootstrap_resamples"])
    raw = backend.read_unseen_labels(sp)  # one read per file
    for tier, col in UNSEEN_TIERS[sp].items():
        try:
            lab = unseen_labels(raw, sp, col)
            pop = population(lab, sp, tier, backend.paths(sp))
            df = pop["primary"]
            entry: dict = {"label_column": col, "counts": pop["counts"], "n_labelled": int(len(lab))}
            if len(df) == 0 or df["label"].sum() == 0:
                entry["status"] = "unavailable: no called positive pairs"
                out[tier] = entry
                continue
            S = score_frame(backend, df)
            bs = make_bootstrap(S, B, int(succ["bootstrap_seed"]))
            mask = np.ones(len(S), bool)
            entry.update(delta_block(S, S, mask, "model", "B1", bs, level, extra=["B2"]))
            entry["delta_vs_B2"] = delta_block(S, S, mask, "model", "B2", bs, level)["delta_vs_B2"]
            entry["n_groups"] = int(len(bs.groups))
            entry["group_definition"] = "connected components of the gene graph of the labelled pairs"
            if "n_lethal_single" in df.columns:  # D4: without positives whose pair has a lethal single mutant
                nls = pd.to_numeric(df["n_lethal_single"], errors="coerce").fillna(0).to_numpy()
                entry["without_lethal_single_positives"] = {
                    "n_positives_removed": int(((df["label"].to_numpy() == 1) & (nls > 0)).sum()),
                    **delta_block(S, S, ~((df["label"].to_numpy() == 1) & (nls > 0)), "model", "B1", bs, level,
                                  extra=["B2"])}
            out[tier] = entry
        except Exception as e:  # one tier failing must not hide the other
            out[tier] = {"error": f"{type(e).__name__}: {e}"}
    return out


def preflight(backend) -> dict:
    """Score every called pair of every evaluated species from the feature tables alone (no labels, no sealed
    file, no access record) and refuse on any problem, so that a feature / schema / NaN bug cannot burn the
    single access. Covers the frozen model, the single-species refits (yeast-only, human-only), B1 and B2."""
    fails: list[str] = []
    report: dict = {}
    for sp in SPECIES + UNSEEN:
        try:
            paths = backend.paths(sp)
            fails += backend.table_failures(sp)
            first = load_feature_table(Path(next(iter(paths.values()))))
            _, df, _, _ = join_called(first[["gene_a", "gene_b"]].copy(), paths)
            need = list(dict.fromkeys(backend.features + backend.b2_cols))
            missing = [c for c in need if c not in df.columns]
            if missing:
                fails.append(f"{sp}: feature tables lack columns {missing}")
                continue
            checks = {"B1": b1_scores(df, backend.b1_variant),
                      "model": np.asarray(backend.score(df), dtype=float),
                      "B2": np.asarray(backend.b2_scorer()(df), dtype=float),
                      "yeast_only_model": np.asarray(backend.yeast_scorer()(df), dtype=float),
                      "human_only_model": np.asarray(backend.human_scorer()(df), dtype=float)}
            for c in backend.b2_cols:
                checks[f"B2 input {c}"] = df[c].to_numpy(float)
            bad = {k: int((~np.isfinite(v)).sum()) for k, v in checks.items() if not np.isfinite(v).all()}
            if bad:
                fails.append(f"{sp}: non-finite values among {len(df)} called pairs: {bad}")
            report[sp] = {"called_pairs_scored": int(len(df)), "n_model_features": len(backend.features)}
        except Exception as e:  # any failure here is a refusal, not a crash after the record
            fails.append(f"{sp}: preflight error {type(e).__name__}: {e}")
    if fails:
        raise Refusal(["preflight: " + f for f in fails])
    return report


# ====================================================================================== secondary analyses
# Each takes the shared state ``st`` and returns a JSON-able dict. They run after the primary result is on disk,
# each inside its own try/except (compute): an error is recorded, the others continue.

def sec_delta_vs_B2(st) -> dict:
    p = st.primary
    return {"T": p["T_vs_B2"], "per_species": {sp: p["per_species"][sp]["delta_vs_B2"] for sp in st.species}}


def sec_human_ge3_studies(st) -> dict:
    return delta_block(st.S, st.S, st.hm & (st.nst >= 3), "model", "B1", st.bs, st.level)


def sec_human_exact_study_set_cells(st) -> dict:
    hs = st.S[st.hm].reset_index(drop=True)
    if not (hs["cell"] != "").any():
        return {"status": "unavailable: studies column missing"}
    pt, bt, info = cell_weighted(st.bs, hs, ["model", "B1"], "cell")
    return {**info, "statistic": "positive-weighted mean over exact study-set cells of the within-cell AP "
                                 "(CI from the resamples of the full test population)",
            "ap": pt, "delta_vs_B1": {"point": pt["model"] - pt["B1"], **ci(bt["model"] - bt["B1"], st.level)}}


def _only_model(st, trained_on: str, tested_on: str) -> dict:
    S = st.S.copy()
    m = (S["species"] == tested_on).to_numpy()
    scorer = st.backend.species_scorer(trained_on)
    S["m_only"] = np.nan
    S.loc[m, "m_only"] = np.asarray(scorer(st.df.loc[m].reset_index(drop=True)), dtype=float)
    if not np.isfinite(S.loc[m, "m_only"].to_numpy()).all():
        raise ValueError(f"non-finite {trained_on}-only scores")
    r = delta_block(S, S, m, "m_only", "B1", st.bs, st.level, extra=["model"])
    r["note"] = (f"m_only: same type, features and hyper-parameters as the frozen model, fitted on development "
                 f"{trained_on} pairs only; 'model' = the frozen model on the same pairs")
    return r


def sec_yeast_only_model_on_human(st) -> dict:
    r = _only_model(st, "yeast", "human")
    r["ap"]["model_yeast_only"] = r["ap"].pop("m_only")
    if "m_only" in r.get("delta_vs_B1", {}):
        r["delta_vs_B1"]["model_yeast_only"] = r["delta_vs_B1"].pop("m_only")
    return r


def sec_human_only_model_on_yeast(st) -> dict:
    r = _only_model(st, "human", "yeast")
    r["ap"]["model_human_only"] = r["ap"].pop("m_only")
    if "m_only" in r.get("delta_vs_B1", {}):
        r["delta_vs_B1"]["model_human_only"] = r["delta_vs_B1"].pop("m_only")
    return r


def _meths(st) -> list[str]:
    return ["model", "B1", "B2"] + (["B4"] if "B4" in st.S else [])


def sec_human_strata_D8(st) -> dict:
    """D8(a) B4 next to every human AP, D8(b) strata 1/2/>=3 and the stratum-averaged AP, D8(c) without
    Harle/Chan-dependent labels, statistics R2 item 13 n-lines strata (dev definitions)."""
    df, S, hm, nst, mp = st.df, st.S, st.hm, st.nst, st.min_pos
    meths = _meths(st)
    ci_m = ["model", "B2"]
    out = {"by_n_studies": {lab: subset_compare(df, S, hm & m, meths, "B1", st.bs, st.level, mp, ci_m)
                            for lab, m in (("1", nst == 1), ("2", nst == 2), (">=3", nst >= 3))}}
    y = S["label"].to_numpy()
    out["stratum_averaged_ap"] = {m: controls.strat_avg_ap(y, S[m].to_numpy(float), nst, hm) for m in meths}
    if "meta_harle_chan_dependent" in df.columns:
        hc = df["meta_harle_chan_dependent"].astype(str).str.lower().isin(["true", "1"]).to_numpy()
        out["without_harle_chan_dependent"] = subset_compare(df, S, hm & ~hc, meths, "B1", st.bs, st.level, mp, ci_m)
        out["harle_chan_dependent_counts"] = {"n": int((hm & hc).sum()), "n_pos": int(y[hm & hc].sum())}
    if "meta_n_testable" in df.columns:
        nl = pd.to_numeric(df["meta_n_testable"], errors="coerce").to_numpy()
        r13 = {}
        for lo, hi in st.mcfg["controls"]["n_lines_strata"]:
            r13[f"n_lines_{lo}-{hi}"] = subset_compare(df, S, hm & (nl >= lo) & (nl <= hi), meths, "B1", st.bs,
                                                       st.level, mp, ci_m)
        out["n_lines_strata"] = r13
    return out


def sec_yeast_and_ribosomal_D10(st) -> dict:
    """D10: yeast without orientation-conflict positives; ribosomal (UniProt keyword, evaluation only) pairs."""
    from slparalog.eval import annotations

    df, S, mp = st.df, st.S, st.min_pos
    meths = ["model", "B1", "B2"]
    ci_m = ["model", "B2"]
    ym = (S["species"] == "yeast").to_numpy()
    out: dict = {"all_yeast": subset_compare(df, S, ym, meths, "B1", st.bs, st.level, mp, ci_m)}
    if "meta_orientation_conflict" in df.columns:
        oc = pd.to_numeric(df["meta_orientation_conflict"], errors="coerce").fillna(0).to_numpy() == 1
        pos = S["label"].to_numpy() == 1
        out["without_orientation_conflict_positives"] = subset_compare(df, S, ym & ~(oc & pos), meths, "B1", st.bs,
                                                                       st.level, mp, ci_m)
        out["orientation_conflict_positives_removed"] = int((ym & oc & pos).sum())
    kw = st.mcfg["controls"]["ribosomal_keyword"]
    ribo = {}
    for sp in st.species:
        genes = st.backend.ribosomal_genes(sp)
        if genes is None:
            ribo[sp] = {"status": "unavailable (UniProt keyword file missing)"}
            continue
        m_sp = (S["species"] == sp).to_numpy()
        flag = annotations.pair_flag(df, genes, "both").to_numpy()
        ribo[sp] = {"status": "ok", "n_genes_with_keyword": len(genes),
                    "ribosomal_pairs": subset_compare(df, S, m_sp & flag, meths, "B1", st.bs, st.level, mp, ci_m),
                    "without_ribosomal_pairs": subset_compare(df, S, m_sp & ~flag, meths, "B1", st.bs, st.level, mp, ci_m)}
    out["ribosomal_" + kw] = ribo
    return out


def sec_identity_tertiles(st) -> dict:
    q = st.mcfg["controls"]["identity_strata"]
    out = {}
    for sp in st.species:
        m_sp = (st.S["species"] == sp).to_numpy()
        qs = controls.quantile_strata(st.S.loc[m_sp, "B1"], q)
        full_q = pd.Series(-1, index=st.S.index)
        full_q[m_sp] = qs.to_numpy()
        for k in range(q):
            m = (full_q == k).to_numpy()
            r = subset_compare(st.df, st.S, m, ["model", "B1", "B2"], "B1", st.bs, st.level, st.min_pos, ["model", "B2"])
            r["b1_range"] = [float(st.S.loc[m, "B1"].min()), float(st.S.loc[m, "B1"].max())]
            out[f"{sp}_q{k + 1}"] = r
    return out


def sec_end_to_end(st) -> dict:
    """DESIGN 6.4: all labelled primary-tier test pairs, caller-missed pairs at the minimum score (the dev
    harness' controls.end_to_end_ap). Point estimates, as on dev."""
    lab_all = pd.concat([st.pops[sp]["labelled"] for sp in st.species], ignore_index=True)
    sc = st.S[["pair_id", "model", "B1", "B2"]].set_index("pair_id")
    out = controls.end_to_end_ap(lab_all, sc, ["model", "B1", "B2"], st.species)
    for sp in st.species:
        o = out[sp]
        o["n_neg"] = o["n"] - o["n_pos"]
        o["prevalence"] = o["n_pos"] / o["n"]
        o["delta_vs_B1"] = o["ap"]["model"] - o["ap"]["B1"]
    out["T_vs_B1"] = float(np.mean([out[sp]["delta_vs_B1"] for sp in st.species]))
    out["note"] = "caller-missed pairs at the minimum score for every method, B1 included; point estimates"
    return out


def sec_tiers(st) -> dict:
    """Lenient / strict / >=2-studies / BioGRID tiers (dataset.TIERS): same frozen model, own bootstrap per tier."""
    out = {}
    for tier in ("lenient", "two_studies", "strict", "biogrid"):
        try:
            pops = {}
            for sp in st.species:
                if sp in TIERS[tier]:
                    pops[sp] = population(load_dev_labels(sp, tier, raw=st.raw[sp]), sp, tier, st.backend.paths(sp))
            df_t = pd.concat([p["primary"] for p in pops.values()], ignore_index=True)
            S_t = score_frame(st.backend, df_t)
            sps = [sp for sp in pops if S_t.loc[S_t["species"] == sp, "label"].sum() > 0]
            if not sps:
                out[tier] = {"status": "unavailable: no called positives"}
                continue
            bs_t = make_bootstrap(S_t, int(st.succ["bootstrap_resamples"]), int(st.succ["bootstrap_seed"]))
            pb = primary_block(S_t, bs_t, st.succ, sps, st.backend.mde())
            out[tier] = {"species": sps, "T": pb["T"], "T_vs_B2": pb["T_vs_B2"], "per_species": pb["per_species"],
                         "relative_gain_mean": pb["relative_gain_mean"],
                         "counts": {sp: pops[sp]["counts"] for sp in pops}, "tier_columns": {sp: TIERS[tier][sp][0] for sp in pops},
                         "note": "descriptive; T is the mean over the species listed"}
        except Exception as e:
            out[tier] = {"error": f"{type(e).__name__}: {e}"}
    return out


def sec_leakage_tertiles(st) -> dict:
    """DESIGN 6.1: dAP by homology-leakage tertile of the test pair (leakage_bsr, written by make_splits).
    Cut points = terciles of leakage_bsr>0 over all primary-labelled test pairs of both species; no-hit pairs
    (leakage 0) are reported separately. Trend = slope of the species-mean dAP over tertiles 1..3."""
    vals = []
    for sp in st.species:
        col, pos, neg = TIERS["primary"][sp]
        raw = st.raw[sp]
        keep = (raw[col] == pos) | (raw[col] == neg)
        vals.append(pd.to_numeric(raw.loc[keep, "leakage_bsr"], errors="coerce").to_numpy())
    allv = np.concatenate(vals)
    hit = allv[np.isfinite(allv) & (allv > 0)]
    cuts = [float(x) for x in np.quantile(hit, [1 / 3, 2 / 3])] if len(hit) else []
    lk = np.full(len(st.df), np.nan)
    for sp in st.species:
        mp_ = raw_map(st.raw[sp], sp, "leakage_bsr")
        idx = np.flatnonzero((st.S["species"] == sp).to_numpy())
        lk[idx] = pd.to_numeric(pd.Series([mp_.get(p) for p in st.S["pair_id"].to_numpy()[idx]]), errors="coerce").to_numpy()
    if not np.isfinite(lk).all():
        raise ValueError("leakage_bsr missing for some test pairs")
    tert = np.where(lk <= 0, "no_hit", np.where(lk <= cuts[0], "T1", np.where(lk <= cuts[1], "T2", "T3"))) if cuts \
        else np.where(lk <= 0, "no_hit", "T1")
    y = st.S["label"].to_numpy()
    out: dict = {"cut_points": cuts, "definition": "max over the pair's genes of bits/selfbits to any development gene",
                 "per_tertile": {}}
    Tpt, Tbt = {}, {}
    for k in ("no_hit", "T1", "T2", "T3"):
        m = tert == k
        blk = {sp: subset_compare(st.df, st.S, m & (st.S["species"] == sp).to_numpy(), ["model", "B1", "B2"], "B1",
                                  st.bs, st.level, st.min_pos, ["model"]) for sp in st.species}
        out["per_tertile"][k] = {"n": int(m.sum()), "n_pos": int(y[m].sum()), "per_species": blk}
        if k == "no_hit":
            continue
        pts = []
        for sp in st.species:
            mm = m & (st.S["species"] == sp).to_numpy()
            if y[mm].sum() > 0:
                pts.append(ap(y[mm], st.S.loc[mm, "model"]) - ap(y[mm], st.S.loc[mm, "B1"]))
        Tpt[k] = float(np.mean(pts)) if pts else float("nan")
        boot = species_boot(st.bs, st.S, ["model", "B1"], "ap", mask=m, species=st.species)
        # mean over the species that have a positive in the resample (the point estimate does the same)
        d_ = np.vstack([boot[sp]["model"] - boot[sp]["B1"] for sp in boot]) if boot else np.full((1, st.bs.B), np.nan)
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            Tbt[k] = np.nanmean(d_, axis=0)
        out["per_tertile"][k]["T_species_mean"] = {"point": Tpt[k], **ci(Tbt[k], st.level)}
    out["trend_slope_per_tertile"] = {"point": (Tpt["T3"] - Tpt["T1"]) / 2, **ci((Tbt["T3"] - Tbt["T1"]) / 2, st.level),
                                      "definition": "(T3 - T1) / 2 of the mean dAP over the species with a positive in the tertile; no-hit pairs excluded"}
    return out


def sec_unseen_species(st) -> dict:
    return {sp: unseen_block(st.backend, sp, st.succ, st.mcfg) for sp in UNSEEN}


SECONDARY_NAMES = ["delta_vs_B2", "human_ge3_studies", "human_exact_study_set_cells", "yeast_only_model_on_human",
                   "human_only_model_on_yeast", "human_strata_D8", "yeast_and_ribosomal_D10", "identity_tertiles",
                   "end_to_end", "tiers", "leakage_tertiles", "unseen_species"]


# ====================================================================================== the pipeline

def compute(backend, succ: dict, header: dict | None = None, on_primary=None) -> dict:
    """Everything the single access computes. ``backend`` supplies labels, features, the model and the refits.
    The primary block is computed first and handed to ``on_primary`` (which writes it to disk) before any
    secondary analysis runs; every secondary runs inside its own try/except and records its error."""
    import traceback
    from types import SimpleNamespace

    level = float(succ["ci_level"])
    B, seed = int(succ["bootstrap_resamples"]), int(succ["bootstrap_seed"])
    mcfg = backend.mcfg
    # ---- test population: labelled primary-tier pairs called by the caller (DESIGN 6.2)
    pops, raws = {}, {}
    for sp in SPECIES:
        raws[sp] = backend.read_test_labels(sp)
        pops[sp] = population(load_dev_labels(sp, "primary", raw=raws[sp]), sp, "primary", backend.paths(sp))
    df = pd.concat([pops[sp]["primary"] for sp in SPECIES], ignore_index=True)
    if df["pair_id"].duplicated().any():
        raise ValueError("duplicate pair IDs in the test population")
    S = score_frame(backend, df)
    hm = (S["species"] == "human").to_numpy()
    species = [sp for sp in SPECIES if (S["species"] == sp).any()]
    bs = make_bootstrap(S, B, seed)
    primary = primary_block(S, bs, succ, species, backend.mde())
    population_info = {sp: pops[sp]["counts"] for sp in SPECIES}
    res = {**(header or {}), "status": "primary_only", "primary": primary, "population": population_info,
           "labelled_pairs_not_called_excluded_from_primary": {
               sp: {k: population_info[sp].get(k) for k in ("not_called", "not_called_pos", "not_called_neg")}
               for sp in SPECIES},
           "model": backend.provenance(), "config": {k: succ[k] for k in succ},
           "secondary": {}, "secondary_errors": {}}
    if on_primary is not None:
        on_primary(res)  # the primary result is on disk before any secondary can fail

    # ---- shared state for the secondaries (evaluation-only attributes come from the raw test tables)
    if "meta_n_studies_testable" in df:
        nst = pd.to_numeric(df["meta_n_studies_testable"], errors="coerce").to_numpy()
    else:
        nst = np.full(len(df), np.nan)
    studies = raw_map(raws["human"], "human", "studies") if "studies" in raws["human"].columns else {}
    S["cell"] = [(",".join(sorted(str(studies[p]).split(","))) if isinstance(studies.get(p), str) else "")
                 if sp == "human" else "" for p, sp in zip(S["pair_id"], S["species"])]
    if "meta_n_testable" in df:  # B4: ascertainment baseline, human only, diagnostic (D8(a))
        S["B4"] = np.where(hm, pd.to_numeric(df["meta_n_testable"], errors="coerce"), np.nan)
    st = SimpleNamespace(backend=backend, succ=succ, mcfg=mcfg, level=level, df=df, S=S, bs=bs, hm=hm, nst=nst,
                         pops=pops, raw=raws, species=species, primary=primary,
                         min_pos=int(mcfg["controls"]["min_pos_for_ci"]))
    for name in SECONDARY_NAMES:
        try:
            res["secondary"][name] = globals()["sec_" + name](st)
        except Exception as e:
            res["secondary_errors"][name] = {"error": f"{type(e).__name__}: {e}",
                                             "traceback": traceback.format_exc().splitlines()[-6:]}
    res["status"] = "complete" if not res["secondary_errors"] else "complete_with_secondary_errors"
    res["model"] = backend.provenance()
    return res


# ====================================================================================== report

def _f(x, nd=3) -> str:
    return "nan" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{nd}f}"


def _d(d: dict, nd=3) -> str:
    return f"{_f(d['point'], nd)} [{_f(d['lo'], nd)}, {_f(d['hi'], nd)}]"


def render_md(res: dict) -> str:
    p, sec = res["primary"], res.get("secondary", {})
    L = [f"# Final test ({res['purpose']}; {res.get('status', '')})", "",
         f"**{p['verdict']}**", "",
         f"T = mean over species of dAP (model - B1) = {_d(p['T'])} (95% percentile CI, {p['T']['resamples']} stratified "
         f"group-bootstrap resamples, seed {p['T']['seed']}). Relative gain (point estimate, mean over species) "
         f"{_f(p['relative_gain_mean'])} vs threshold {_f(p['relative_gain_threshold'], 2)}. "
         f"Criteria: {json.dumps(p['criteria'])}. Pre-declared underpowered: {p['predeclared_underpowered']}.", "",
         "| species | n | n_pos | n_neg | prevalence | AP model | AP B1 | AP B2 | dAP vs B1 [95% CI] | dAP vs B2 [95% CI] |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for sp, r in p["per_species"].items():
        L.append(f"| {sp} | {r['n']} | {r['n_pos']} | {r['n_neg']} | {_f(r['prevalence'])} | {_f(r['ap_model'])} | "
                 f"{_f(r['ap_B1'])} | {_f(r['ap_B2'])} | {_d(r['delta_vs_B1'])} | {_d(r['delta_vs_B2'])} |")
    L += ["", f"T vs B2: {_d(p['T_vs_B2'])}", "", "## Secondary (pre-declared, D17 and DESIGN 6.1-6.4)", "",
          "| analysis | n | n_pos | n_neg | prevalence | AP model | AP B1 | dAP vs B1 [95% CI] |", "|---|---|---|---|---|---|---|---|"]

    def row(name: str, b: dict | None, m: str = "model") -> None:
        if not isinstance(b, dict):
            return
        if "ap" not in b or m not in b.get("ap", {}):
            L.append(f"| {name} | {b.get('n', '')} | {b.get('n_pos', '')} | {b.get('n_neg', '')} | | | | "
                     f"{b.get('status', b.get('error', 'no positives'))} |")
            return
        d = (b.get("delta_vs_B1") or {}).get(m)
        L.append(f"| {name} | {b['n']} | {b['n_pos']} | {b['n_neg']} | {_f(b['prevalence'])} | {_f(b['ap'][m])} | "
                 f"{_f(b['ap']['B1'])} | {_d(d) if d else ''} |")

    row("human, >=3 testing studies", sec.get("human_ge3_studies"))
    c = sec.get("human_exact_study_set_cells") or {}
    if "ap" in c:
        L.append(f"| human, exact study-set cells ({c['n_cells_used']} cells) | {c['n_pairs']} | {c['n_pos']} | {c['n_neg']} | "
                 f"{_f(c['prevalence'])} | {_f(c['ap']['model'])} | {_f(c['ap']['B1'])} | {_d(c['delta_vs_B1'])} |")
    y = sec.get("yeast_only_model_on_human") or {}
    row("human, yeast-only-trained model", y, "model_yeast_only")
    row("human, frozen model (same pairs)", y, "model")
    h = sec.get("human_only_model_on_yeast") or {}
    row("yeast, human-only-trained model", h, "model_human_only")
    row("yeast, frozen model (same pairs)", h, "model")
    d8 = sec.get("human_strata_D8") or {}
    for lab, b in (d8.get("by_n_studies") or {}).items():
        row(f"human, {lab} testing studies", b)
    row("human, without Harle/Chan-dependent labels", d8.get("without_harle_chan_dependent"))
    d10 = sec.get("yeast_and_ribosomal_D10") or {}
    row("yeast, without orientation-conflict positives", d10.get("without_orientation_conflict_positives"))
    for sp, tiers in (sec.get("unseen_species") or {}).items():
        for t, b in tiers.items():
            if isinstance(b, dict) and "counts" in b:
                head = ", headline" if (sp, t) == ("arabidopsis", "strict") else ""
                row(f"{sp} ({t}{head}, unseen)", b)
                if "without_lethal_single_positives" in b:
                    row(f"{sp} ({t}, without lethal-single positives)", b["without_lethal_single_positives"])
    L += ["", "| tier / analysis | species | n_pos | AP model | AP B1 | dAP vs B1 [95% CI] |", "|---|---|---|---|---|---|"]
    for tier, b in (sec.get("tiers") or {}).items():
        for sp, r in (b.get("per_species") or {}).items():
            L.append(f"| tier {tier} | {sp} | {r['n_pos']} | {_f(r['ap_model'])} | {_f(r['ap_B1'])} | {_d(r['delta_vs_B1'])} |")
        if "error" in b or "status" in b:
            L.append(f"| tier {tier} | | | | | {b.get('error', b.get('status'))} |")
    for sp, r in ((sec.get("end_to_end") or {}).items()):
        if isinstance(r, dict) and "ap" in r:
            L.append(f"| end-to-end ({r['n_not_called']} uncalled at min score) | {sp} | {r['n_pos']} | "
                     f"{_f(r['ap']['model'])} | {_f(r['ap']['B1'])} | {_f(r['delta_vs_B1'])} |")
    lk = sec.get("leakage_tertiles") or {}
    for k, r in (lk.get("per_tertile") or {}).items():
        if "T_species_mean" in r:
            L.append(f"| leakage {k} (species-mean) | both | {r['n_pos']} | | | {_d(r['T_species_mean'])} |")
    if "trend_slope_per_tertile" in lk:
        L.append(f"| leakage trend (slope per tertile) | both | | | | {_d(lk['trend_slope_per_tertile'])} |")
    if res.get("secondary_errors"):
        L += ["", "## Secondary analyses that failed (recorded, primary unaffected)", ""]
        L += [f"- {k}: {v['error']}" for k, v in res["secondary_errors"].items()]
    m = res.get("model", {})
    L += ["", f"B1 variant: {m.get('b1_variant')}; model {m.get('model_id')} "
              f"(sha256 {str(m.get('model_sha256'))[:12]}); git HEAD {str(res.get('git', {}).get('head'))[:12]}.", ""]
    return "\n".join(L)


# ====================================================================================== entry points

def _write_durable(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + ".part")
    with open(tmp, "w") as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def write_outputs(out_dir: Path, res: dict) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    jp, mp = out_dir / RESULT_NAME, out_dir / "final_test.md"
    _write_durable(jp, json.dumps(res, indent=2, sort_keys=False, default=str) + "\n")
    _write_durable(mp, render_md(res))
    return {"json_sha256": file_sha256(jp), "md_sha256": file_sha256(mp)}


def run_access(ctx: Ctx, run_tests: bool = False) -> dict:
    if run_tests and check_tests_marker(ctx):
        record_tests(ctx)
    failures = run_preconditions(ctx)
    if failures:
        log_refusal(ctx, failures)
        raise Refusal(failures)
    meta = json.loads((ctx.model_dir / "model.json").read_text())
    succ = load_yaml(ctx.config_dir / "success.yaml")
    bad = verify_sealed(ctx)  # hashes only; refuses without consuming the access
    if bad:
        log_refusal(ctx, bad)
        raise Refusal(bad)
    try:
        backend = RealBackend(ctx)
        pre = preflight(backend)
    except Refusal as e:
        log_refusal(ctx, e.failures)
        raise
    except Exception as e:
        fails = [f"preflight: could not load / score: {type(e).__name__}: {e}"]
        log_refusal(ctx, fails)
        raise Refusal(fails) from e
    log_append(ctx, "access_started", note="preconditions, sealed hashes and preflight passed", model_sha=meta["model_sha256"])
    ctx.access_started = True  # from here on read_sealed may parse
    header = {"purpose": "primary", "timestamp_utc": _now(), "script": "scripts/eval/final_test.py",
              "git": {"head": head_sha(ctx), "tag": ctx.tag, "tag_commit": tag_commit(ctx)},
              "config_sha256": config_sha256(ctx), "preflight": pre}
    try:
        res = compute(backend, succ, header=header, on_primary=lambda r: write_outputs(ctx.out_dir, r))
        res["model"]["sealed_files_sha256"] = dict(ctx.sealed_read)
        hashes = write_outputs(ctx.out_dir, res)
    except BaseException as e:
        log_append(ctx, "access_failed", note=f"{type(e).__name__}: {e}", model_sha=meta["model_sha256"])
        raise
    log_append(ctx, "access_completed", note=json.dumps({**hashes, "status": res["status"],
                                                         "secondary_errors": sorted(res["secondary_errors"])}),
               model_sha=meta["model_sha256"])
    return res


TOY_DIR = ROOT / "tests" / "fixtures" / "final_test_toy"


class ToyBackend:
    """Synthetic fixture labels and feature tables with hand-computable answers (tests/fixtures/final_test_toy).
    The 'models' are stand-ins that rank by one registered column each; everything downstream (population join,
    B1, bootstrap, strata, cells, unseen species, report) is the real code."""

    toy = True
    b1_variant = "ident_local_max"
    table_names = ["tierA"]
    features = ["bits_pair"]
    b2_cols = ["cov_min"]

    def table_failures(self, sp):
        return []

    def b2_scorer(self):
        return self.fit_b2()

    def species_scorer(self, sp):
        return lambda d: d["len_min"].to_numpy(float)  # stand-in for both single-species refits

    yeast_scorer = human_scorer = lambda self: self.species_scorer("x")  # noqa: E731

    def ribosomal_genes(self, sp):
        return None  # evaluation-only keyword files are project data: unavailable in toy mode

    @property
    def mcfg(self):
        return load_yaml(ROOT / "config" / "models.yaml")

    def __init__(self, fixture_dir: Path = TOY_DIR):
        self.dir = Path(fixture_dir)
        if is_sealed(self.dir, ROOT / "data" / "labels" / "sealed"):
            raise Refusal(f"toy mode must never touch sealed paths: {self.dir}")

    def read_test_labels(self, sp):
        return read_fixture(self.dir / f"{sp}_test.tsv")

    def read_unseen_labels(self, sp):
        return read_fixture(self.dir / f"unseen_{sp}.tsv")

    def paths(self, sp):
        return feature_paths(self.dir / "features", sp, self.table_names)

    def score(self, df):
        return df["bits_pair"].to_numpy(float)

    def fit_b2(self):
        return lambda d: d["cov_min"].to_numpy(float)

    def provenance(self):
        return {"model_id": "toy-standin", "model_sha256": "0" * 64, "b1_variant": self.b1_variant,
                "stand_ins": {"model": "bits_pair", "B2": "cov_min", "single_species_refits": "len_min"}}

    def mde(self):
        return {"underpowered": None}


def run_toy(out_dir: Path, fixture_dir: Path = TOY_DIR) -> dict:
    """Whole pipeline on the toy fixtures. Never touches sealed paths, the access log or results/final."""
    out_dir = Path(out_dir)
    if out_dir.resolve() == (ROOT / "results" / "final").resolve():
        raise Refusal("toy output must not go to results/final")
    succ = load_yaml(ROOT / "config" / "success.yaml")
    backend = ToyBackend(fixture_dir)
    pre = preflight(backend)
    header = {"purpose": "dry-run-toy", "timestamp_utc": _now(), "script": "scripts/eval/final_test.py",
              "git": {"head": "toy", "tag": None, "tag_commit": None}, "preflight": pre}
    res = compute(backend, succ, header=header, on_primary=lambda r: write_outputs(out_dir, r))
    write_outputs(out_dir, res)
    return res


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dry-run-toy", action="store_true", help="run the pipeline on tests/fixtures; never touches sealed paths")
    p.add_argument("--out-dir", default=None, help="toy mode only: output directory (default logs/final_toy (git-ignored))")
    p.add_argument("--record-tests", action="store_true", help="run the full pytest suite at HEAD, write the marker, exit")
    p.add_argument("--run-tests", action="store_true", help="if no valid pytest marker, run the full suite first")
    p.add_argument("--model-dir", default="results/frozen/model", help="frozen model artifact, relative to the project root")
    return p.parse_args(argv)


def main(argv=None, root: Path | None = None) -> int:
    args = parse_args(argv)
    if args.dry_run_toy:
        if args.record_tests or args.run_tests:
            print("--dry-run-toy cannot be combined with the real-access options", file=sys.stderr)
            return 2
        res = run_toy(Path(args.out_dir) if args.out_dir else ROOT / "logs" / "final_toy")
        print(render_md(res))
        return 0
    if args.out_dir:
        print("--out-dir is only valid with --dry-run-toy", file=sys.stderr)
        return 2
    ctx = Ctx(root=Path(root or ROOT), model_rel=args.model_dir)
    try:
        if args.record_tests:
            record_tests(ctx)
            print(f"pytest passed at HEAD; marker written to {marker_path(ctx)}")
            return 0
        res = run_access(ctx, run_tests=args.run_tests)
    except Refusal as e:
        print("REFUSED:\n  - " + "\n  - ".join(e.failures), file=sys.stderr)
        return 3
    print(render_md(res))
    print(f"\nwritten: {ctx.out_dir / RESULT_NAME} and final_test.md; the access log and results are uncommitted: commit them now.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
