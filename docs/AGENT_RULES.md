# AGENT RULES — preamble for every delegated task

Project: .. Read PROMPT.md (fixed goal and rules) and the parts of DESIGN.md relevant to your task.

1. **Environment:** use only `.venv/bin/python`; install only into that venv; binaries under `bin/`; model weights under `models_cache/` (HF_HOME/TORCH_HOME). Never use system/conda Python, sudo or brew. Check `df -h` before downloads >1 GB. Keep everything under the project directory.
2. **Forbidden inputs to features/models:** expression, PPI, essentiality, genetic interactions, GO/Pfam/any annotation, labels of orthologous pairs, any hand-entered knowledge of specific genes or pairs. Never hard-code gene names/IDs outside test fixtures and label-source configs. Never read `data/raw/_quarantine/` except via the documented column whitelist in label code.
3. **Sealed data:** never read `data/labels/sealed/` (test and unseen-species labels). Only `scripts/eval/final_test.py` may, and only the PI runs it.
4. **Protected files:** do not change `config/labels.yaml`, `config/seeds.yaml`, metric code, split files or the leakage audit without an explicit PI instruction. Seeds come from `config/seeds.yaml`.
5. **Metrics:** AUPRC = `sklearn.metrics.average_precision_score`. Always report prevalence, n positives, n negatives.
6. **Reproducibility:** scripts, not one-offs; results written as JSON/TSV under `results/` by the script that computed them; no numbers typed by hand into reports. Log each development experiment in `results/experiments.tsv` with its adoption rule declared before running.
7. **Verification:** test new code on toy data with known answers; add pytest tests. State clearly what you verified and what you did not. Never fabricate; mark uncertainty.
8. **Process isolation:** do not import torch and fit LightGBM in the same process (libomp segfault).
9. **Git:** do not commit unless the task says so; never commit raw/restricted data.
10. **Integrity:** never edit a failing test, its tolerances, pytest configuration or skip/xfail markers to make it pass; never rewrite git history; never edit RESEARCH_TREE.md or DESIGN.md (PI only); never read raw unseen-species sources (Cusack, Tischler, pombe, FlyBase) unless your task is building those sealed labels; every number you report cites the result file (and commit, if committed) it came from.
11. **Report back** concisely (≤300 words unless asked otherwise): what you did, key numbers with file paths, what is unverified, open problems.
