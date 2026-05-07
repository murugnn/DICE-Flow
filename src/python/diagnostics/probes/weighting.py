"""Probe stage 5: weight_network — explain every edge that disappears."""
from __future__ import annotations
from typing import Any, Dict
import numpy as np
import pandas as pd

from ..probe_base import Probe


class WeightingProbe(Probe):
    name = "weighting"
    stage = "5-weighting"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        merged = pd.read_csv(ctx["merged_network"], sep="\t")
        weighted = pd.read_csv(ctx["weighted_network"], sep="\t")
        expr = pd.read_csv(ctx["expression"], sep="\t", index_col=0)
        eps = float(ctx.get("epsilon", 1e-6))
        scale = int(ctx.get("scale_factor", 1000))

        merged_pairs = set(zip(merged["source"], merged["target"]))
        weighted_pairs = set(zip(weighted["source"], weighted["target"]))
        dropped_pairs = merged_pairs - weighted_pairs

        # Decompose drops by reason
        expr_genes = set(expr.index)
        unresolved = sum(1 for s, t in dropped_pairs if s not in expr_genes or t not in expr_genes)

        # Resolution rate: union of all genes in merged net vs expression
        all_merged_genes = set(merged["source"]).union(merged["target"])
        resolution_rate = len(all_merged_genes & expr_genes) / max(1, len(all_merged_genes))

        # Epsilon clamp hits: edges with weight near epsilon (within 1.5x)
        eps_clamp = int((weighted["weight"] <= eps * 1.5).sum())

        # Integer-truncation collisions
        int_w = (weighted["weight"] * scale).astype(int)
        unique_weights = weighted["weight"].nunique()
        unique_int = int_w.nunique()
        truncation_collisions = int(unique_weights - unique_int)

        # Bottom-clamp from weight_int==0 -> 1 in current code
        bottom_clamp_to_1 = int((weighted["weight_int"] == 1).sum())

        return {
            "n_input_edges": len(merged),
            "n_output_edges": len(weighted),
            "n_dropped": len(dropped_pairs),
            "edges_unresolved_in_expression": unresolved,
            "gene_symbol_resolution_rate": float(resolution_rate),
            "epsilon_clamped_edges": eps_clamp,
            "weight_truncation_collisions": truncation_collisions,
            "weight_int_eq_1_count": bottom_clamp_to_1,
            "weight_distribution": {
                "min": float(weighted["weight"].min()),
                "max": float(weighted["weight"].max()),
                "p10": float(np.quantile(weighted["weight"], 0.1)),
                "p50": float(np.quantile(weighted["weight"], 0.5)),
                "p90": float(np.quantile(weighted["weight"], 0.9)),
            },
        }
