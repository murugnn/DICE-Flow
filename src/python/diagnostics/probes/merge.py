"""Probe stage 3: merged network — directionality + structural stats."""
from __future__ import annotations
from typing import Any, Dict
import pandas as pd

from ..probe_base import Probe


class MergeProbe(Probe):
    name = "merge_networks"
    stage = "3-merge"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        path = ctx["merged_network"]
        df = pd.read_csv(path, sep="\t")
        n_edges = len(df)
        unique_genes = len(set(df["source"]).union(df["target"]))
        if "directed" in df.columns:
            n_dir = int(df["directed"].astype(str).str.lower().eq("true").sum())
            n_und = n_edges - n_dir
        else:
            n_dir, n_und = 0, n_edges

        # bidirectional pairs (X->Y and Y->X both present in directed subset)
        if "directed" in df.columns:
            d = df[df["directed"].astype(str).str.lower() == "true"]
        else:
            d = df.iloc[0:0]
        pairs = set(zip(d["source"], d["target"]))
        bid = sum(1 for (s, t) in pairs if (t, s) in pairs and s < t)

        return {
            "path": path,
            "n_edges": n_edges,
            "n_directed": n_dir,
            "n_undirected": n_und,
            "unique_genes": unique_genes,
            "bidirectional_pairs": bid,
            "self_loops": int((df["source"] == df["target"]).sum()),
        }
