"""Probe stage 2: parsed KEGG / BioGRID outputs."""
from __future__ import annotations
from typing import Any, Dict
import pandas as pd

from ..probe_base import Probe


def _edge_stats(path: str) -> Dict[str, Any]:
    df = pd.read_csv(path, sep="\t")
    src_col = next((c for c in df.columns if c.lower() in ("source", "src")), df.columns[0])
    tgt_col = next((c for c in df.columns if c.lower() in ("target", "tgt", "dst")), df.columns[1])
    pairs = df[[src_col, tgt_col]].astype(str)
    duplicates = pairs.duplicated().sum()
    unique_genes = len(set(pairs[src_col]).union(pairs[tgt_col]))
    self_loops = int((pairs[src_col] == pairs[tgt_col]).sum())
    return {
        "path": path,
        "n_edges": len(pairs),
        "duplicate_edges": int(duplicates),
        "unique_genes": int(unique_genes),
        "self_loops": self_loops,
    }


class ParsersProbe(Probe):
    name = "parsers"
    stage = "2-parsers"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        m: Dict[str, Any] = {}
        if ctx.get("kegg_interactions"):
            m["kegg"] = _edge_stats(ctx["kegg_interactions"])
        if ctx.get("biogrid_interactions"):
            m["biogrid"] = _edge_stats(ctx["biogrid_interactions"])
        return m
