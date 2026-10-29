"""G8 (reviews/round2/agent_failure_modes.md): feature registry checks."""

import dataclasses

import pytest

from slparalog.features import registry
from slparalog.features.registry import Feature
from slparalog.guards import FORBIDDEN_NAME_RE


def test_registry_valid():
    assert registry.validate() == []


def test_sources_allowed_and_names_clean():
    for f in registry.FEATURES:
        assert f.source in {"self_search", "panel", "plm", "structure", "sequence"}, f.name
        assert not FORBIDDEN_NAME_RE.search(f.name), f.name


def test_forbidden_regex_catches_examples():
    for bad in ("gtex_tau", "ppi_degree", "is_essential", "lowcomplexity", "ceres_min",
                "depmap", "go_sim", "pfam_shared", "interpro_n", "expr_corr", "ygob_flag",
                "ohnolog", "label", "sl"):
        assert FORBIDDEN_NAME_RE.search(bad), bad
        assert registry.validate([Feature(bad, "sequence", "x")]), bad
    assert registry.validate([Feature("ok_name", "annotation", "x")])  # bad source


def test_percentile_versions_registered():
    names = {f.name for f in registry.FEATURES}
    for f in registry.FEATURES:
        if f.pct:
            assert f"{f.name}_pct" in names


def test_registry_hash_changes_with_definition():
    h = registry.registry_hash()
    alt = list(registry.FEATURES[:-1])
    assert registry.registry_hash(alt) != h
    assert len(h) == 64


def _replace(name, **kw):
    return [dataclasses.replace(f, **kw) if f.name == name else f for f in registry.FEATURES]


def test_columns_hash_excludes_descriptions():
    cols = registry.TIER_A_COLUMNS
    h = registry.columns_hash(cols)
    assert len(h) == 64 and h == registry.columns_hash(reversed(cols))  # order-independent
    reworded = [dataclasses.replace(f, description=f.description + " (reworded)") for f in registry.FEATURES]
    assert registry.registry_hash(reworded) != registry.registry_hash()  # the global hash sees descriptions...
    assert registry.columns_hash(cols, reworded) == h  # ...the column hash does not


def test_columns_hash_detects_definition_changes():
    cols = registry.TIER_A_COLUMNS
    h = registry.columns_hash(cols)
    for field, value in (("dtype", "float32"), ("source", "panel"), ("kind", "qc"), ("tier", "B"), ("pct", True)):
        name = "len_ratio"
        assert getattr(registry.BY_NAME[name], field) != value
        assert registry.columns_hash(cols, _replace(name, **{field: value})) != h, field
    with pytest.raises(ValueError):
        registry.columns_hash(["no_such_column"])


def test_columns_hash_depends_only_on_its_own_columns():
    a, b = registry.TIER_A_COLUMNS, registry.TIER_B_COLUMNS
    ha = registry.columns_hash(a)
    # redefining a Tier B column, or dropping tiers B and C, leaves the Tier A hash unchanged
    assert registry.columns_hash(a, _replace(b[0], dtype="float32")) == ha
    assert registry.columns_hash(a, [f for f in registry.FEATURES if f.tier == "A"]) == ha
    assert registry.columns_hash(b, _replace(b[0], dtype="float32")) != registry.columns_hash(b)
    assert registry.table_columns_hash(["gene_a", "gene_b", *a]) == ha  # identifiers are skipped


def test_pipeline_columns_subset_of_registry():
    # every non-identifier column the feature stage can emit is registered
    assert set(registry.TIER_A_COLUMNS) <= set(registry.BY_NAME)
    # (was "<= TIER_A_COLUMNS" while Tier A was the only tier; Tier B columns are registered too)
    assert set(registry.MODEL_FEATURES) <= set(registry.BY_NAME)
    assert {registry.BY_NAME[c].tier for c in registry.TIER_A_COLUMNS} == {"A"}
