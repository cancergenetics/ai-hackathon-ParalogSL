"""Name guard for model features (gate G8, reviews/round2/agent_failure_modes.md).

Kept outside ``slparalog/features`` on purpose: the G1 guard (tests/test_data_access.py) forbids
feature code from mentioning evaluation-only annotation sources, and this regex has to name them
in order to reject them.
"""

import re

FORBIDDEN_NAME_RE = re.compile(
    r"(?i)gtex|ppi|essential|complex|ceres|depmap|go_|pfam|interpro|expr|ygob|ohnolog|label|^sl$"
)
