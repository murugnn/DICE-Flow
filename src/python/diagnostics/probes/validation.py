"""Probe stage 8: validation pipelines — DEG overlap, DepMap join, survival cohort."""
from __future__ import annotations
from typing import Any, Dict
import pandas as pd
from scipy import stats

from ..probe_base import Probe


def _read_set(path):
    with open(path) as f:
        return {ln.strip() for ln in f if ln.strip() and ln.strip() != "gene"}


class ValidationProbe(Probe):
    name = "validation"
    stage = "8-validation"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        m: Dict[str, Any] = {}
        if ctx.get("dice_genes") and ctx.get("deg_genes") and ctx.get("universe_size"):
            dice = _read_set(ctx["dice_genes"])
            deg  = _read_set(ctx["deg_genes"])
            obs = len(dice & deg)
            N = int(ctx["universe_size"])
            # Hypergeom: prob >= obs given draws=|dice|, hits=|deg|, total=N
            rv = stats.hypergeom(N, len(deg), len(dice))
            p = float(rv.sf(obs - 1))
            m["dice_deg_overlap"] = {
                "observed": int(obs),
                "expected": float(len(dice) * len(deg) / max(1, N)),
                "hypergeom_p": p,
            }
        elif ctx.get("dice_genes") and ctx.get("deg_genes"):
            dice = _read_set(ctx["dice_genes"])
            deg  = _read_set(ctx["deg_genes"])
            obs = len(dice & deg)
            N = max(1, len(dice | deg)) * 10  # crude default universe
            rv = stats.hypergeom(N, len(deg), len(dice))
            m["dice_deg_overlap"] = {
                "observed": int(obs),
                "expected": float(len(dice) * len(deg) / N),
                "hypergeom_p": float(rv.sf(obs - 1)),
                "universe_assumed": N,
            }

        if ctx.get("depmap_table"):
            dm = pd.read_csv(ctx["depmap_table"], sep="\t")
            dice = _read_set(ctx["dice_genes"])
            matched = dm[dm["gene"].isin(dice)]
            m["depmap"] = {
                "depmap_genes": int(len(dm)),
                "join_match_count": int(len(matched)),
                "join_match_rate": float(len(matched) / max(1, len(dice))),
            }

        if ctx.get("survival_cohorts"):
            m["survival"] = {
                "cohort_sizes": ctx["survival_cohorts"],   # passed in by orchestrator
            }

        return m
