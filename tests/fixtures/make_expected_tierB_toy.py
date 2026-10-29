"""Regenerate tests/fixtures/expected_tierB_toy.parquet (stored Tier B table of the toy proteome).

Only run this when a Tier B feature definition is changed ON PURPOSE (and say so in the log); the
regression test exists to catch unintended value changes.

    .venv/bin/python tests/fixtures/make_expected_tierB_toy.py
"""

import copy
import sys
import tempfile
from pathlib import Path

import yaml

FIX = Path(__file__).resolve().parent
sys.path.insert(0, str(FIX))
import make_toy_panel as mk  # noqa: E402

from slparalog import panel, pipeline  # noqa: E402
from slparalog.features import registry  # noqa: E402

with tempfile.TemporaryDirectory() as t:
    tmp = Path(t)
    query, pan, _ = mk.build()
    qf, pf, st = mk.write(tmp, query, pan)
    pc = copy.deepcopy(panel.load_panel_config())
    pc["panel"].update({"fasta": str(pf), "species_table": str(st)})
    pc["tier_b"].update({"db_dir": str(tmp / "db"), "cache_dir": str(tmp / "hits")})
    pc["tier_b"]["core"]["min_size"] = 10
    (tmp / "panel.yaml").write_text(yaml.safe_dump(pc))
    df, _ = pipeline.run(qf, None, threads=2, features_only=True, tier_b=True,
                         panel_config_path=tmp / "panel.yaml", workdir=tmp / "work")
    out = df[["gene_a", "gene_b"] + registry.TIER_B_COLUMNS].reset_index(drop=True)
    out.to_parquet(FIX / "expected_tierB_toy.parquet", index=False)
    print(out.shape)
