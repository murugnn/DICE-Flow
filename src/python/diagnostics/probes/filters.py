"""Probe stage 4: DEA + IG + intersection."""
from __future__ import annotations
from typing import Any, Dict
import numpy as np
import pandas as pd

from ..probe_base import Probe


def _read_gene_list(path: str) -> set:
    with open(path) as f:
        return {ln.strip() for ln in f if ln.strip()}


class FiltersProbe(Probe):
    name = "filters"
    stage = "4-filters"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        m: Dict[str, Any] = {}
        if ctx.get("dea_table"):
            df = pd.read_csv(ctx["dea_table"], sep="\t")
            sig = df[(df["FDR"] < 0.05) & (df["log2FC"].abs() > 1.0)]
            m["dea"] = {
                "n_genes_total": len(df),
                "n_significant": len(sig),
                "fdr_p10": float(np.quantile(df["FDR"], 0.1)),
                "fdr_p50": float(np.quantile(df["FDR"], 0.5)),
                "fdr_p90": float(np.quantile(df["FDR"], 0.9)),
                "fdr_le_0_05_pct": float((df["FDR"] < 0.05).mean()),
                "abs_lfc_gt_1_pct": float((df["log2FC"].abs() > 1.0).mean()),
            }
        dea_set = _read_gene_list(ctx["dea_genes"]) if ctx.get("dea_genes") else set()
        ig_set = _read_gene_list(ctx["ig_genes"]) if ctx.get("ig_genes") else set()
        inter_set = _read_gene_list(ctx["candidates"]) if ctx.get("candidates") else set()
        if dea_set:
            m.setdefault("dea", {})["n_genes_listed"] = len(dea_set)
        if ig_set:
            m["ig"] = {"n_genes": len(ig_set)}
        if dea_set and ig_set:
            jac = len(dea_set & ig_set) / max(1, len(dea_set | ig_set))
            m["intersection"] = {
                "n_genes": len(inter_set),
                "n_dea_only": len(dea_set - ig_set),
                "n_ig_only": len(ig_set - dea_set),
                "jaccard_dea_ig": float(jac),
            }
        return m
