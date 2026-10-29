"""Known-answer tests for the unseen-species label parsers (counts only; nothing is printed per pair).

Tolerances were fixed before the first run of this test file:
  Cusack strict (Dataset_1): 51 positive / 110 negative  -- exact (tolerance 0), before the negative re-check
  Cusack broad  (Dataset_2): 190 positive / 110 negative -- exact (tolerance 0), before the negative re-check
  Tischler tested pairs: 143 -- exact (tolerance 0)
  Tischler synthetic positives: 16 (127 negatives) -- exact; the positives are not contained in the Additional
  files on disk (2 and 4), Table 4 is now parsed from the Europe PMC full text (scripts/data/extract_tischler_table4.py); the check counts
  pairs before identifier-mapping loss (D6).
"""
import importlib.util
from pathlib import Path

import pytest
import yaml

pytestmark = pytest.mark.data

ROOT = Path(__file__).resolve().parents[1]
UNIPROT = ROOT / "data" / "raw" / "uniprot_2026_03"
CFG = yaml.safe_load((ROOT / "config" / "labels_unseen.yaml").read_text())


def _load_builder():
    spec = importlib.util.spec_from_file_location("build_unseen_labels", ROOT / "scripts" / "labels" / "build_unseen_labels.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B = _load_builder()
needs_cusack = pytest.mark.skipif(not (ROOT / "data/raw/cusack2021/Dataset_1.txt").exists(), reason="raw Cusack files absent")
needs_tischler = pytest.mark.skipif(not (ROOT / "data/raw/tischler2006/gb-2006-7-8-r69-S2.doc").exists(), reason="raw Tischler files absent")


@needs_cusack
def test_cusack_strict_and_broad_known_answers():
    rows, cnt = B.build_arabidopsis(CFG, UNIPROT)
    ka = CFG["arabidopsis"]["known_answer"]
    assert cnt["raw_dataset1"] == {"pos": ka["strict_positive"], "neg": ka["strict_negative"]} == {"pos": 51, "neg": 110}
    assert cnt["raw_dataset2"] == {"pos": ka["broad_positive"], "neg": ka["broad_negative"]} == {"pos": 190, "neg": 110}
    # no mapping loss expected; mapped tiers before the re-check equal the published counts
    assert cnt["tiers_before_recheck_mapped"]["strict"] == {"pos": 51, "neg": 110}
    assert cnt["tiers_before_recheck_mapped"]["broad"] == {"pos": 190, "neg": 110}


@needs_cusack
def test_cusack_pairs_unordered_unique_and_recheck_only_removes_negatives():
    rows, cnt = B.build_arabidopsis(CFG, UNIPROT)
    keys = [(r[0], r[1]) for r in rows]
    assert all(a < b for a, b in keys) and len(set(keys)) == len(keys)
    assert cnt["tiers_after_recheck"]["strict"]["pos"] == 51
    assert cnt["tiers_after_recheck"]["strict"]["neg"] == 110 - cnt["negatives_excluded_lethal_single"]


@needs_tischler
def test_tischler_tested_pairs_143():
    _, cnt = B.build_worm(CFG, UNIPROT)
    assert cnt["tested_pairs"] == 143
    assert cnt["tested_in_s4"] == 143


@needs_tischler
def test_tischler_16_positives_127_negatives():
    # D6: the published known answer (16 synthetic positives among 143 tested pairs) is checked on the
    # parsed source, i.e. before identifier-mapping losses, like the Cusack checks before the re-check.
    _, cnt = B.build_worm(CFG, UNIPROT)
    assert cnt.get("positives") == 16
    assert cnt.get("negatives") + cnt.get("mapping_loss_pairs") == 127
