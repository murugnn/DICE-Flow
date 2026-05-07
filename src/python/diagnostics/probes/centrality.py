"""Probe stage 6: consume C++ --diagnostics-out JSON; compute Spearman."""
from __future__ import annotations
import json
from typing import Any, Dict
from scipy import stats

from ..probe_base import Probe


class CentralityProbe(Probe):
    name = "centrality"
    stage = "6-centrality"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        with open(ctx["cpp_diagnostics_json"]) as f:
            d = json.load(f)
        full = d["sampled_full_betweenness"]
        smp = d["sampled_betweenness"]
        rho, p = stats.spearmanr(full, smp)
        return {
            "omp_determinism_max_relative_drift": d["omp_determinism_max_relative_drift"],
            "sampled_vs_full_subgraph_n": d["sampled_vs_full_subgraph_n"],
            "sampled_vs_full_spearman_rho": float(rho),
            "sampled_vs_full_spearman_p": float(p),
        }
