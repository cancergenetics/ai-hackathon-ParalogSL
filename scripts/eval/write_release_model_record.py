#!/usr/bin/env python
"""Write results/release/release_model.json: counts, hashes and provenance of the release model.

Loads results/release/model with slparalog.model.load_model, compares its spec to the frozen model
(model.json), and records training counts from the model artifact. Run after
`slparalog train ... -o results/release/model` (command recorded below).
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

from slparalog import model as M

ROOT = Path(__file__).resolve().parents[2]
REL, FRZ = ROOT / "results/release/model", ROOT / "results/frozen/model"
CMD = ("slparalog train --features human=data/interim/features/human/tierA.parquet "
       "--features human=data/interim/features/human/tierB.parquet "
       "--features yeast=data/interim/features/yeast/tierA.parquet "
       "--features yeast=data/interim/features/yeast/tierB.parquet "
       "--labels human=data/labels/release/human_release.tsv --labels yeast=data/labels/release/yeast_release.tsv "
       "--feature-set A+B --model-type gbdt --tier primary --seed 11 "
       "--param l2_regularization=1.0 --param learning_rate=0.03 --param max_iter=100 "
       "--param max_leaf_nodes=7 --param min_samples_leaf=50 -o results/release/model")


def sha(p: Path) -> str:
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()


def main() -> None:
    meta = json.loads((REL / "model.json").read_text())
    fro = json.loads((FRZ / "model.json").read_text())
    loaded = M.load_model(REL)  # must not raise
    same = {k: meta[k] == fro[k] for k in ("model_type", "feature_set", "features", "params", "registry_hash",
                                           "columns_hash", "caller_config_sha256", "required_tiers")}
    same["seed_used"] = meta["seeds"]["seed_used"] == fro["seeds"]["seed_used"]
    tr = meta["training"]
    per_sp = tr.get("per_species") or {k: v for k, v in tr.items() if k in ("human", "yeast")}
    out = {
        "model_id": meta["model_id"], "model_dir": "results/release/model", "load_model_ok": loaded is not None,
        "model_sha256_recorded": meta["model_sha256"], "model_joblib_sha256": sha(REL / "model.joblib"),
        "model_json_sha256": sha(REL / "model.json"),
        "frozen_model_id": fro["model_id"], "frozen_model_sha256": fro["model_sha256"],
        "spec_identical_to_frozen": same, "n_features": len(meta["features"]),
        "params": meta["params"], "seed": meta["seeds"]["seed_used"], "tier": tr.get("tier"),
        "training": {k: tr[k] for k in tr if k not in ("label_files",)},
        "label_files": {sp: {"path": v["path"], "sha256": v["sha256"]} for sp, v in tr["label_files"].items()},
        "train_command": CMD, "git_head": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                                                         text=True).stdout.strip(),
        "statement": ("RELEASE MODEL: refit of the frozen specification on development + held-out test primary-tier "
                      "pairs of human and yeast. No held-out evaluation of this model exists or is possible: every "
                      "labelled human and yeast pair was used for training. Its performance is represented only by "
                      "the frozen model's final test (results/final/final_test.md). The unseen Arabidopsis and worm "
                      "labels and data/labels/sealed/ were not used."),
    }
    (ROOT / "results/release/release_model.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps({k: out[k] for k in ("model_id", "load_model_ok", "spec_identical_to_frozen", "training")}, indent=1, default=str))


if __name__ == "__main__":
    main()
