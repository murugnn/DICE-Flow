"""Render probe results as markdown + JSON."""
from __future__ import annotations
from datetime import datetime, timezone
import json
from typing import List

from .probe_base import ProbeResult


def _fmt_value(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    if isinstance(v, list):
        return f"[{', '.join(_fmt_value(x) for x in v[:10])}{'...' if len(v) > 10 else ''}]"
    return str(v)


def _fmt_dict(d, indent=0):
    lines = []
    pre = "  " * indent
    for k, v in d.items():
        if isinstance(v, dict):
            lines.append(f"{pre}- **{k}**:")
            lines.extend(_fmt_dict(v, indent + 1))
        else:
            lines.append(f"{pre}- **{k}**: {_fmt_value(v)}")
    return lines


def render_markdown(results: List[ProbeResult], commit: str) -> str:
    out = [
        "# Diagnostic Dossier",
        "",
        f"- **Commit**: `{commit}`",
        f"- **Generated**: {datetime.now(timezone.utc).isoformat()}",
        f"- **N probes**: {len(results)}",
        "",
        "---",
    ]
    by_stage = sorted(results, key=lambda r: r.stage)
    for r in by_stage:
        out.append(f"\n## Stage {r.stage} — `{r.name}`\n")
        if r.error:
            out.append(f"**ERROR**: `{r.error}`\n")
            continue
        out.append(f"- _duration_: {r.duration_s:.3f}s")
        out.extend(_fmt_dict(r.measurements))
    out.append("\n---\n")
    out.append("## Bug Registry\n")
    out.append("_Populate from per-stage findings; sort by impact × confidence × blast_radius._\n")
    out.append("| ID | Stage | Symptom | Evidence | Confidence | Blast radius | Hypothesized root cause |")
    out.append("|----|-------|---------|----------|------------|--------------|-------------------------|")
    return "\n".join(out) + "\n"


def render_json(results: List[ProbeResult], commit: str) -> str:
    return json.dumps({
        "commit": commit,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "probes": [
            {
                "name": r.name,
                "stage": r.stage,
                "duration_s": r.duration_s,
                "error": r.error,
                "measurements": r.measurements,
            }
            for r in results
        ],
    }, indent=2, default=str)
