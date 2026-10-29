#!/usr/bin/env python
"""Extract the list of delegated agents (description, model, type) from the PI session transcript.

Usage: extract_agent_log.py TRANSCRIPT.jsonl -o results/process/agents.tsv
Only the PI's own Agent tool calls are read (sidechain/subagent messages are skipped). Agents spawned
without an explicit model inherit the PI model (Opus); they are recorded as `opus (inherited)`.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("transcript")
    ap.add_argument("-o", "--output", required=True)
    a = ap.parse_args()
    seen, rows = set(), []
    for line in open(a.transcript, encoding="utf-8"):
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        m = d.get("message") or {}
        if m.get("role") != "assistant" or d.get("isSidechain") or not isinstance(m.get("content"), list):
            continue
        for c in m["content"]:
            if c.get("type") == "tool_use" and c.get("name") == "Agent" and c["id"] not in seen:
                seen.add(c["id"])
                i = c["input"]
                model = {"opus": "opus", "sonnet": "sonnet"}.get(i.get("model"), "opus (inherited)")
                rows.append({"order": len(rows) + 1, "timestamp_utc": d.get("timestamp", ""),
                             "description": i.get("description", ""), "model": model})
    out = Path(a.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["order", "timestamp_utc", "description", "model"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} agents -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
