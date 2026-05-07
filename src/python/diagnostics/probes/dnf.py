"""Probe stage 10: DNF — pair count, DSS distribution, SynLethDB overlap."""
from __future__ import annotations
from typing import Any, Dict
import numpy as np
import pandas as pd

from ..probe_base import Probe


class DNFProbe(Probe):
    name = "dnf"
    stage = "10-dnf"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        m: Dict[str, Any] = {}
        if ctx.get("dnf_pairwise"):
            df = pd.read_csv(ctx["dnf_pairwise"], sep="\t")
            score_col = next((c for c in df.columns if c.upper() == "DSS"), None)
            m["n_pairs"] = int(len(df))
            if score_col is not None and len(df):
                m["dss_min"] = float(df[score_col].min())
                m["dss_max"] = float(df[score_col].max())
                m["dss_p50"] = float(np.quantile(df[score_col], 0.5))
                m["dss_positive_pct"] = float((df[score_col] > 0).mean())
        if ctx.get("synlethdb") and ctx.get("dnf_pairwise"):
            sl = pd.read_csv(ctx["synlethdb"], sep="\t")
            sl_pairs = set(map(tuple, sl[["gene_a", "gene_b"]].astype(str).values))
            sl_pairs |= {(b, a) for (a, b) in sl_pairs}  # undirected
            df = pd.read_csv(ctx["dnf_pairwise"], sep="\t")
            dnf_pairs = set(map(tuple, df[["gene_a", "gene_b"]].astype(str).values))
            m["synlethdb_pairs"] = len(sl_pairs) // 2
            m["synlethdb_overlap"] = sum(1 for p in dnf_pairs if p in sl_pairs)
        return m
