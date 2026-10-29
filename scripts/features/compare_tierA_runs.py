"""Compare two Tier A outputs of the same proteome (e.g. before/after a caller-rule change).

Usage: .venv/bin/python scripts/features/compare_tierA_runs.py <species> <old.parquet> <new.parquet> <tag>
Writes results/features/tierA_compare_<tag>_<species>.json: pair-set relation (identical / superset /
subset / neither) with counts, and per-feature max |difference| over common pairs.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main(species, old_p, new_p, tag):
    old = pd.read_parquet(old_p).set_index(["gene_a", "gene_b"])
    new = pd.read_parquet(new_p).set_index(["gene_a", "gene_b"])
    so, sn = set(old.index), set(new.index)
    common = sorted(so & sn)
    rel = ("identical" if so == sn else "new_superset" if sn > so else
           "new_subset" if sn < so else "neither")
    diffs = {}
    o, n = old.loc[common], new.loc[common]
    for c in new.columns.intersection(old.columns):
        d = np.abs(o[c].to_numpy(np.float64) - n[c].to_numpy(np.float64))
        diffs[c] = {"max_abs_diff": float(np.nanmax(d)) if len(d) else 0.0,
                    "n_pairs_differing": int((d > 1e-12).sum())}
    only_new = new.loc[sorted(sn - so)]
    only_old = old.loc[sorted(so - sn)]
    out = {"species": species, "old": str(old_p), "new": str(new_p), "relation": rel,
           "n_old": len(so), "n_new": len(sn), "n_common": len(common),
           "n_only_old": len(so - sn), "n_only_new": len(sn - so),
           "frac_new_not_in_old": len(sn - so) / max(len(sn), 1),
           "frac_old_not_in_new": len(so - sn) / max(len(so), 1),
           "only_new_bits_pair_median": float(only_new["bits_pair"].median()) if len(only_new) else None,
           "only_old_bits_pair_median": float(only_old["bits_pair"].median()) if len(only_old) else None,
           "feature_diffs_on_common_pairs": diffs}
    dst = ROOT / f"results/features/tierA_compare_{tag}_{species}.json"
    dst.write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "feature_diffs_on_common_pairs"}, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:5])
