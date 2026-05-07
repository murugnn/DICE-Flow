"""Probe stage 7: differential_centrality — rank-product double-count + hub kill list."""
from __future__ import annotations
from typing import Any, Dict
import math

import pandas as pd
from scipy import stats

from ..probe_base import Probe


def _read_set(path):
    with open(path) as f:
        return {ln.strip() for ln in f if ln.strip()}


class DifferentialProbe(Probe):
    name = "differential"
    stage = "7-differential"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        cn = pd.read_csv(ctx["centrality_normal"], sep="\t")
        ct = pd.read_csv(ctx["centrality_tumor"], sep="\t")
        merged = cn.merge(ct, on="gene", suffixes=("_n", "_t"))

        diff_bet = (merged["betweenness_t"] - merged["betweenness_n"]).abs()
        diff_eig = (merged["eigenvector_t"] - merged["eigenvector_n"]).abs()
        # spearmanr returns NaN if either input is constant. Treat that as
        # "independence" (rho=0) so the dossier always carries a number.
        rho, _ = stats.spearmanr(diff_bet, diff_eig)
        if rho is None or (isinstance(rho, float) and math.isnan(rho)):
            rho = 0.0

        # Top-12% eigenvector hub list (would-be-dropped under DiCE noise filter)
        eig_thresh_n = merged["eigenvector_n"].quantile(0.88)
        eig_thresh_t = merged["eigenvector_t"].quantile(0.88)
        kill_mask = (merged["eigenvector_n"] > eig_thresh_n) | (merged["eigenvector_t"] > eig_thresh_t)
        kill_genes = merged.loc[kill_mask, "gene"].tolist()

        cgc = _read_set(ctx["cgc_list"]) if ctx.get("cgc_list") else set()
        in_cgc = sum(1 for g in kill_genes if g in cgc)

        return {
            "n_genes": len(merged),
            "bet_eig_rank_correlation": float(rho),
            "top12_eig_dropped_n": len(kill_genes),
            "top12_eig_dropped_genes": kill_genes[:50],   # cap for dossier
            "top12_eig_dropped_in_cgc": in_cgc,
            "cgc_size": len(cgc),
        }
