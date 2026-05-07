"""Probe stage 1: raw input files (KEGG, BioGRID, expression)."""
from __future__ import annotations
import re
from pathlib import Path
from typing import Any, Dict

import pandas as pd

from ..probe_base import Probe

ENTREZ_RE = re.compile(r"^\d+$")


def _detect_id_type(idx: pd.Index) -> str:
    sample = [str(x) for x in idx[:200]]
    if not sample:
        return "empty"
    entrez_hits = sum(1 for s in sample if ENTREZ_RE.match(s))
    if entrez_hits / len(sample) > 0.9:
        return "entrez"
    if entrez_hits / len(sample) > 0.1:
        return "mixed"
    return "symbol"


def _expr_stats(path: str) -> Dict[str, Any]:
    df = pd.read_csv(path, sep="\t", index_col=0)
    return {
        "path": path,
        "n_genes": len(df),
        "n_samples": df.shape[1],
        "nan_cells": int(df.isna().sum().sum()),
        "duplicate_genes": int(df.index.duplicated().sum()),
        "id_type": _detect_id_type(df.index),
    }


class RawInputsProbe(Probe):
    name = "raw_inputs"
    stage = "1-raw-inputs"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        m: Dict[str, Any] = {}
        if ctx.get("expression_normal"):
            m["expression_normal"] = _expr_stats(ctx["expression_normal"])
        if ctx.get("expression_tumor"):
            m["expression_tumor"] = _expr_stats(ctx["expression_tumor"])

        if ctx.get("kegg_dir"):
            kegg = Path(ctx["kegg_dir"])
            xmls = list(kegg.glob("*.xml")) if kegg.exists() else []
            m["kegg"] = {"path": str(kegg), "n_xml_files": len(xmls)}

        if ctx.get("biogrid_file"):
            bg = Path(ctx["biogrid_file"])
            size = bg.stat().st_size if bg.exists() else 0
            m["biogrid"] = {"path": str(bg), "exists": bg.exists(), "size_bytes": size}

        return m
