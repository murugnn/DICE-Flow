"""CLI: run all probes, render markdown + JSON. Pure observation."""
from __future__ import annotations
import argparse
import logging
import subprocess
from pathlib import Path
from typing import Optional

from .probe_base import Probe, ProbeResult
from .probes.raw_inputs import RawInputsProbe
from .probes.parsers import ParsersProbe
from .probes.merge import MergeProbe
from .probes.filters import FiltersProbe
from .probes.weighting import WeightingProbe
from .probes.centrality import CentralityProbe
from .probes.differential import DifferentialProbe
from .probes.validation import ValidationProbe
from .probes.dnf import DNFProbe
from .render import render_markdown, render_json

log = logging.getLogger(__name__)


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def main(argv: Optional[list] = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output", required=True)
    p.add_argument("--expression-normal")
    p.add_argument("--expression-tumor")
    p.add_argument("--kegg-dir")
    p.add_argument("--biogrid-file")
    p.add_argument("--kegg-interactions")
    p.add_argument("--biogrid-interactions")
    p.add_argument("--merged-network")
    p.add_argument("--dea-table"); p.add_argument("--dea-genes")
    p.add_argument("--ig-genes"); p.add_argument("--candidates")
    p.add_argument("--weighted-network"); p.add_argument("--expression")
    p.add_argument("--cpp-diagnostics-json")
    p.add_argument("--centrality-normal"); p.add_argument("--centrality-tumor")
    p.add_argument("--cgc-list")
    p.add_argument("--dice-genes"); p.add_argument("--deg-genes")
    p.add_argument("--depmap-table")
    p.add_argument("--dnf-pairwise"); p.add_argument("--synlethdb")
    args = p.parse_args(argv)

    ctx = {k.replace("-", "_"): getattr(args, k.replace("-", "_"))
           for k in vars(args) if k != "output"}

    probes: list[Probe] = [
        RawInputsProbe(), ParsersProbe(), MergeProbe(), FiltersProbe(),
        WeightingProbe(), CentralityProbe(), DifferentialProbe(),
        ValidationProbe(), DNFProbe(),
    ]
    results: list[ProbeResult] = []
    for probe in probes:
        # Skip probes whose required ctx keys are absent
        try:
            results.append(probe.run(ctx))
        except Exception as e:  # noqa: BLE001
            log.exception(f"probe {probe.name} failed unexpectedly")

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    commit = _git_commit()
    (out / "diagnostic_dossier.md").write_text(render_markdown(results, commit))
    (out / "dossier.json").write_text(render_json(results, commit))
    print(f"wrote {out}/diagnostic_dossier.md and dossier.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
