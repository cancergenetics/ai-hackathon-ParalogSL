"""``slparalog setup`` with mocked downloads: no network is used.

The planted toy panel of tests/fixtures/make_toy_panel.py is served as gzipped per-species FASTA files from an
in-memory "server"; checksums, manifest entries and the species table are generated for it. Covers checksum
verification, idempotency, failure on a corrupted source, panel assembly + MMseqs2 database build, ESM-2 weight
verification, model installation, and an end-to-end ``predict`` of a Tier A + B model through a set-up panel.
"""

import gzip
import hashlib
import json
import shutil
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

import helpers_model as hm
from slparalog import cli, model, panel, pipeline, predict, setup_data, train
from slparalog.config import CONFIG_DIR

FIX = Path(__file__).parent / "fixtures"
sys.path.insert(0, str(FIX))
import make_toy_panel as mk  # noqa: E402

sys.path.pop(0)
TOY = FIX / "toy_proteome.fasta"
BASE_URL = "https://qfo.example.invalid/Eukaryota/"
CORE, GAP = ["S0", "S1"], ["S2", "S3"]


def _gz_fasta(code, seqs) -> bytes:
    text = "".join(f">sp|A{code}{i:04d}|P{i}_{code} toy protein\n" + "\n".join(textwrap.wrap(s, 70)) + "\n"
                   for i, s in enumerate(seqs))
    return gzip.compress(text.encode(), mtime=0)


class Server:
    """URL -> bytes; counts requests; ``tamper`` replaces the body of one URL."""

    def __init__(self):
        self.files, self.calls, self.tamper = {}, [], {}

    def __call__(self, url, dest: Path):
        self.calls.append(url)
        dest.write_bytes(self.tamper.get(url, self.files[url]))


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("setup_toy")
    query, pan, _ = mk.build()
    qf, _, _ = mk.write(tmp, query, pan)
    srv = Server()
    rows, entries, core_files = [], [], {}
    expected = []
    for k, code in enumerate(mk.SPECIES):
        upid, taxid = f"UP9{k:08d}", str(100 + k)
        blob = _gz_fasta(code, pan[code])
        fn = f"{upid}_{taxid}.fasta.gz"
        is_core = code in CORE
        url = BASE_URL + fn if is_core else f"https://uniprot.example.invalid/{upid}/{fn}"
        srv.files[url] = blob
        sha = hashlib.sha256(blob).hexdigest()
        if is_core:
            core_files[fn] = {"sha256": sha, "bytes": len(blob)}
        else:
            entries.append({"id": f"panel_extra_{upid}", "filename": fn, "url": url, "sha256": sha})
        rows.append({"code": code, "taxid": taxid, "name": f"toy {code}", "source": "toy", "n_proteins": len(pan[code]),
                     "clade": "toy", "upid": upid, "selection_group": "core" if is_core else "gap"})
        for i, s in enumerate(pan[code]):  # independent expectation of the assembled panel
            expected.append(f">{code}|{i}\n" + "\n".join(textwrap.wrap(s, 60)) + "\n")
    st = tmp / "species.tsv"
    pd.DataFrame(rows).to_csv(st, sep="\t", index=False)
    manifest = tmp / "manifest.yaml"
    manifest.write_text(yaml.safe_dump({"entries": [{"id": "qfo_2026_02", "url": BASE_URL.rstrip("/")}] + entries}))
    exp_text = "".join(expected)
    sources = tmp / "sources.json"
    sources.write_text(json.dumps({"qfo_manifest_id": "qfo_2026_02", "files": core_files,
                                   "panel_fasta_sha256": hashlib.sha256(exp_text.encode()).hexdigest()}))
    pc = panel.load_panel_config()
    pc["tier_b"]["core"]["min_size"] = 10   # the toy has 30 core genes (real panel: 200)
    base = tmp / "base_panel.yaml"
    base.write_text(yaml.safe_dump(pc))
    return {"tmp": tmp, "srv": srv, "st": st, "manifest": manifest, "sources": sources, "base": base,
            "query": qf, "expected_fasta": exp_text}


def _panel(toy, dd, **kw):
    return setup_data.setup_panel(dd, threads=2, downloader=toy["srv"], species_table=toy["st"],
                                  manifest_path=toy["manifest"], sources_path=toy["sources"],
                                  base_config=toy["base"], **kw)


@pytest.fixture(scope="module")
def panel_dd(toy, tmp_path_factory):
    dd = tmp_path_factory.mktemp("data_dir")
    toy["srv"].calls.clear()
    r = _panel(toy, dd)
    return dd, r


# ---------------------------------------------------------------------------------------------- panel

def test_panel_setup_downloads_verifies_assembles_and_builds_databases(toy, panel_dd):
    dd, r = panel_dd
    assert r["panel_fasta"] == "assembled" and sorted(r["downloaded"]) == sorted(mk.SPECIES)
    assert len(toy["srv"].calls) == 4
    assert (dd / "panel" / "panel.fasta").read_text() == toy["expected_fasta"]
    pcfg = panel.load_panel_config(dd / "panel" / "panel.yaml")
    assert Path(pcfg["tier_b"]["db_dir"]) == dd / "panel" / "db" and Path(pcfg["panel"]["fasta"]).parent == dd / "panel"
    stamp = json.loads((dd / "panel" / "db" / "stamp.json").read_text())
    assert stamp["codes"] == mk.SPECIES and stamp["panel_fasta_sha256"] == json.loads(toy["sources"].read_text())["panel_fasta_sha256"]
    for c in mk.SPECIES:
        assert (dd / "panel" / "db" / c / "db.dbtype").exists()
    assert not (dd / "panel" / "downloads").exists()          # archives removed after a verified assembly
    assert not list((dd / "panel").rglob("*.part"))


def test_panel_setup_is_idempotent(toy, panel_dd):
    dd, _ = panel_dd
    toy["srv"].calls.clear()
    stamp_before = (dd / "panel" / "db" / "stamp.json").read_text()
    r = _panel(toy, dd)
    assert toy["srv"].calls == [] and r["panel_fasta"] == "present"
    assert (dd / "panel" / "db" / "stamp.json").read_text() == stamp_before   # databases not rebuilt


def test_panel_setup_repairs_a_damaged_panel(toy, tmp_path):
    toy["srv"].calls.clear()
    _panel(toy, tmp_path)
    f = tmp_path / "panel" / "panel.fasta"
    f.write_text(f.read_text() + ">x|0\nAAA\n")
    toy["srv"].calls.clear()
    r = _panel(toy, tmp_path)
    assert r["panel_fasta"] == "assembled" and f.read_text() == toy["expected_fasta"]
    assert len(toy["srv"].calls) == 4


def test_corrupted_download_is_rejected_and_not_kept(toy, tmp_path):
    srv = toy["srv"]
    url = next(u for u in srv.files if "UP900000001" in u)
    srv.tamper[url] = srv.files[url][:-5] + b"xxxxx"
    try:
        with pytest.raises(setup_data.SetupError, match="checksum mismatch"):
            _panel(toy, tmp_path)
    finally:
        srv.tamper.clear()
    dl = tmp_path / "panel" / "downloads"
    assert not (dl / url.rsplit("/", 1)[1]).exists() and not list(dl.glob("*.part"))
    assert not (tmp_path / "panel" / "panel.fasta").exists()


def test_wrong_assembled_checksum_is_rejected(toy, tmp_path):
    bad = tmp_path / "sources_bad.json"
    s = json.loads(toy["sources"].read_text())
    s["panel_fasta_sha256"] = "0" * 64
    bad.write_text(json.dumps(s))
    with pytest.raises(setup_data.SetupError, match="assembled panel.fasta checksum"):
        setup_data.setup_panel(tmp_path, threads=1, downloader=toy["srv"], species_table=toy["st"],
                               manifest_path=toy["manifest"], sources_path=bad, base_config=toy["base"])
    assert not (tmp_path / "panel" / "panel.fasta").exists()


def test_panel_sources_come_from_species_table_manifest_and_checksums(toy):
    sources, _ = setup_data.panel_sources(toy["st"], toy["manifest"], toy["sources"])
    assert [s["code"] for s in sources] == mk.SPECIES
    assert sources[0]["url"].startswith(BASE_URL) and sources[2]["url"].startswith("https://uniprot.example.invalid")
    assert all(len(s["sha256"]) == 64 for s in sources)


def test_real_config_is_consistent():
    """Every species of the real panel has a URL and a checksum in the committed config (no network)."""
    sources, src = setup_data.panel_sources()
    tab = pd.read_csv(CONFIG_DIR / "panel_species.tsv", sep="\t", dtype=str)
    assert [s["code"] for s in sources] == tab["code"].tolist() and len(sources) == 80
    assert len({s["filename"] for s in sources}) == 80 and len(src["panel_fasta_sha256"]) == 64
    assert all(s["url"].startswith("https://") and len(s["sha256"]) == 64 for s in sources)


def test_http_download_rejects_html(monkeypatch, tmp_path):
    import io

    class R(io.BytesIO):
        headers = {}

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    monkeypatch.setattr(setup_data.urllib.request, "urlopen", lambda *a, **k: R(b"<!DOCTYPE html><title>Just a moment"))
    monkeypatch.setattr(setup_data.time, "sleep", lambda s: None)
    with pytest.raises(setup_data.SetupError, match="HTML"):
        setup_data.http_download("https://x.invalid/a.gz", tmp_path / "a.gz", retries=2)
    assert not (tmp_path / "a.gz").exists()


# ---------------------------------------------------------------------------------------------- PLM

def _plm_fixture(tmp_path):
    files = {"config.json": b'{"a": 1}', "model.safetensors": b"weights" * 10, "vocab.txt": b"A\nC\n"}
    rev = "r" * 40
    wj = tmp_path / "plm_weights.json"
    wj.write_text(json.dumps({"toy": {"hf_id": "org/toy", "revision": rev,
                                      "files": {k: hashlib.sha256(v).hexdigest() for k, v in files.items()}}}))
    cfg = {"default_model": "toy", "models": {"toy": {"hf_id": "org/toy", "revision": rev}}}
    calls = []

    def dl(repo_id, revision, cache_dir, names, bodies=files):
        calls.append((repo_id, revision, sorted(names)))
        snap = Path(cache_dir) / "models--org--toy" / "snapshots" / revision
        snap.mkdir(parents=True, exist_ok=True)
        for n in names:
            (snap / n).write_bytes(bodies[n])

    return wj, cfg, calls, dl, files


def test_plm_setup_verifies_and_is_idempotent(tmp_path):
    wj, cfg, calls, dl, files = _plm_fixture(tmp_path)
    r = setup_data.setup_plm(tmp_path, snapshot_download=dl, weights_path=wj, plm_config=cfg)
    assert r["status"] == "downloaded" and len(calls) == 1 and calls[0][2] == sorted(files)
    assert setup_data.plm_home(tmp_path) == tmp_path / "plm" / "hf"
    assert setup_data.setup_plm(tmp_path, snapshot_download=dl, weights_path=wj, plm_config=cfg)["status"] == "present"
    assert len(calls) == 1
    # a damaged file is fetched again
    snap = tmp_path / "plm" / "hf" / "hub" / "models--org--toy" / "snapshots" / ("r" * 40)
    (snap / "vocab.txt").write_bytes(b"corrupt")
    assert setup_data.setup_plm(tmp_path, snapshot_download=dl, weights_path=wj, plm_config=cfg)["status"] == "downloaded"
    assert (snap / "vocab.txt").read_bytes() == files["vocab.txt"]


def test_plm_setup_rejects_wrong_files_and_stale_checksums(tmp_path):
    wj, cfg, calls, dl, files = _plm_fixture(tmp_path)
    bad = lambda *a: dl(*a, bodies={**files, "vocab.txt": b"tampered"})  # noqa: E731
    with pytest.raises(setup_data.SetupError, match="failed verification"):
        setup_data.setup_plm(tmp_path / "d1", snapshot_download=bad, weights_path=wj, plm_config=cfg)
    assert setup_data.plm_home(tmp_path / "d1") is None            # no record of an unverified download
    cfg2 = {"default_model": "toy", "models": {"toy": {"hf_id": "org/toy", "revision": "q" * 40}}}
    with pytest.raises(setup_data.SetupError, match="pins"):
        setup_data.setup_plm(tmp_path / "d2", snapshot_download=dl, weights_path=wj, plm_config=cfg2)


def test_real_plm_checksums_match_pinned_revisions():
    from slparalog import plm

    cfg = plm.load_plm_config()
    wj = json.loads((CONFIG_DIR / "plm_weights.json").read_text())
    assert wj["esm2_150m"]["revision"] == cfg["models"]["esm2_150m"]["revision"]
    assert "model.safetensors" in wj["esm2_150m"]["files"]


# ---------------------------------------------------------------------------------------------- model

@pytest.fixture(scope="module")
def trained(toy, panel_dd, tmp_path_factory):
    """Tier A + B model trained on synthetic labels over the toy query, features from the set-up panel."""
    dd, _ = panel_dd
    tmp = tmp_path_factory.mktemp("trained")
    fp = tmp / "toy_ab.parquet"
    feats, meta = pipeline.run(toy["query"], fp, threads=2, features_only=True, tier_b=True,
                               panel_config_path=dd / "panel" / "panel.yaml", workdir=tmp / "wd")
    labs = hm.write_synthetic_labels(feats, tmp)
    out = tmp / "model_ab"
    m = train.train({"human": [fp], "yeast": [fp]}, labs, "A+B", "lr", out)
    return {"dir": out, "meta": m, "feats": feats}


def test_model_setup_validates_copies_and_is_idempotent(trained, tmp_path):
    r = setup_data.setup_model(tmp_path, trained["dir"])
    assert r["status"] == "installed" and (tmp_path / "models" / r["model_id"] / "model.json").exists()
    assert setup_data.setup_model(tmp_path, trained["dir"])["status"] == "present"
    assert setup_data.resolve_model(None, tmp_path) == tmp_path / "models" / r["model_id"]
    assert setup_data.resolve_model(r["model_id"], tmp_path) == tmp_path / "models" / r["model_id"]
    assert model.load_model(setup_data.resolve_model(None, tmp_path)).tiers == ["A", "B"]


def test_model_setup_refuses_a_stale_model(trained, tmp_path):
    bad = tmp_path / "stale"
    shutil.copytree(trained["dir"], bad)
    m = json.loads((bad / "model.json").read_text())
    m["registry_hash"] = "0" * 64
    m.pop("columns_hash", None)  # legacy artifact: the whole-registry hash decides
    (bad / "model.json").write_text(json.dumps(m))
    with pytest.raises(setup_data.SetupError, match="registry hash mismatch"):
        setup_data.setup_model(tmp_path / "dd", bad)
    assert not (tmp_path / "dd" / "models").exists()


def test_setup_cli_model_and_missing_model_message(trained, tmp_path, capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["setup", "--model", str(trained["dir"]), "--data-dir", str(tmp_path)])
    assert e.value.code == 0 and "model:" in capsys.readouterr().out
    with pytest.raises(SystemExit) as e:
        cli.main(["setup", "--model", str(tmp_path / "nope"), "--data-dir", str(tmp_path)])
    assert e.value.code == 1


def test_data_dir_resolution(monkeypatch, tmp_path):
    monkeypatch.delenv("SLPARALOG_DATA", raising=False)
    assert setup_data.data_dir() == Path("~/.slparalog").expanduser()
    monkeypatch.setenv("SLPARALOG_DATA", str(tmp_path / "env"))
    assert setup_data.data_dir() == tmp_path / "env"
    assert setup_data.data_dir(tmp_path / "arg") == tmp_path / "arg"


# ---------------------------------------------------------------------------------------------- end to end

def test_predict_through_a_set_up_panel_with_qc_flags(toy, panel_dd, trained, tmp_path, capsys):
    """setup (mocked downloads) -> install model -> predict --model (Tier A + B computed from the data dir)."""
    dd, _ = panel_dd
    setup_data.setup_model(dd, trained["dir"])
    out = tmp_path / "pred.tsv"
    with pytest.raises(SystemExit) as e:
        cli.main(["predict", str(toy["query"]), "-o", str(out), "--data-dir", str(dd), "-t", "2"])
    assert e.value.code == 0
    df = pd.read_csv(out, sep="\t", comment="#")
    head = [ln for ln in out.read_text().splitlines() if ln.startswith("# ")]
    assert len(df) == len(trained["feats"]) and {"tb_fragment_suspect", "tb_n_outgroups"} <= set(df.columns)
    assert "tiers_computed: A,B" in "\n".join(head)
    # same scores as the model on the features computed in training
    m = model.load_model(trained["dir"])
    exp = pd.Series(m.score(trained["feats"]), index=pd.MultiIndex.from_frame(trained["feats"][["gene_a", "gene_b"]]))
    got = df.set_index(["gene_a", "gene_b"])["sl_score"]
    np.testing.assert_allclose(got.to_numpy(), exp.loc[got.index].to_numpy(), rtol=1e-4, atol=1e-5)
    assert df["sl_rank"].is_monotonic_increasing


def test_predict_without_panel_fails_with_setup_hint(trained, tmp_path, capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["predict", str(TOY), "-o", str(tmp_path / "o.tsv"), "--model", str(trained["dir"]),
                  "--data-dir", str(tmp_path / "empty")])
    assert e.value.code == 1 and "setup --panel" in capsys.readouterr().err


def test_panel_qc_warning_reaches_header_and_stderr(toy, panel_dd, trained, tmp_path, capsys):
    """A core set below tier_b.core.min_size is a proteome-level QC warning: header comment block + stderr."""
    dd, _ = panel_dd
    pc = yaml.safe_load((dd / "panel" / "panel.yaml").read_text())
    pc["tier_b"]["core"]["min_size"] = 1000      # more than the toy core set
    alt = tmp_path / "panel_alt.yaml"
    alt.write_text(yaml.safe_dump(pc))
    out = tmp_path / "pred.tsv"
    with pytest.raises(SystemExit) as e:
        cli.main(["predict", str(toy["query"]), "-o", str(out), "--model", str(trained["dir"]),
                  "--data-dir", str(dd), "--panel-config", str(alt), "-t", "2"])
    assert e.value.code == 0
    head = [ln for ln in out.read_text().splitlines() if ln.startswith("# ")]
    assert any(ln.startswith("# qc_warning: Tier B: core set") for ln in head)
    assert "QC warning: Tier B: core set" in capsys.readouterr().err


def test_schema_hash_does_not_depend_on_data_dir_paths(toy, panel_dd, tmp_path):
    """Two data dirs with identical contents (different paths) give the same schema_hash; paths are provenance."""
    dd1, _ = panel_dd
    dd2 = tmp_path / "other_data_dir"
    toy["srv"].calls.clear()
    _panel(toy, dd2)
    assert dd1 != dd2
    metas = []
    for dd in (dd1, dd2):
        _, meta = pipeline.run(toy["query"], None, threads=2, features_only=True, tier_b=True,
                               panel_config_path=dd / "panel" / "panel.yaml", workdir=tmp_path / f"wd_{dd.name}")
        metas.append(meta)
    assert metas[0]["schema_hash"] == metas[1]["schema_hash"]
    assert metas[0]["tier_b"]["settings_hash"] == metas[1]["tier_b"]["settings_hash"]
    assert metas[0]["tier_b"]["paths"]["db_dir"] != metas[1]["tier_b"]["paths"]["db_dir"]
    # a content change (panel setting) does change it
    pc = panel.load_panel_config(dd1 / "panel" / "panel.yaml")
    ccfg = pipeline.load_yaml(pipeline.DEFAULT_CALLER_CONFIG)
    h0 = panel.settings_hash(pc, ccfg, "x")
    pc["tier_b"]["presence_min_bits"] = 99
    assert panel.settings_hash(pc, ccfg, "x") != h0
