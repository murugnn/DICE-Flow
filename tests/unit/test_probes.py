import pandas as pd
import pytest
from src.python.diagnostics.probes.raw_inputs import RawInputsProbe
from src.python.diagnostics.probes.parsers import ParsersProbe
from src.python.diagnostics.probes.filters import FiltersProbe


def test_raw_inputs_probe_counts_rows_and_detects_id_type(
    tmp_path, tiny_expr_normal, tiny_expr_tumor
):
    norm = tmp_path / "norm.tsv"
    tum = tmp_path / "tum.tsv"
    tiny_expr_normal.to_csv(norm, sep="\t")
    tiny_expr_tumor.to_csv(tum, sep="\t")

    r = RawInputsProbe().run(ctx={
        "expression_normal": str(norm),
        "expression_tumor": str(tum),
        "kegg_dir": None,           # accepted but skipped if None
        "biogrid_file": None,
    })
    assert r.error is None
    m = r.measurements
    assert m["expression_normal"]["n_genes"] == 10
    assert m["expression_normal"]["n_samples"] == 5
    assert m["expression_normal"]["nan_cells"] == 0
    assert m["expression_normal"]["id_type"] == "symbol"
    assert m["expression_tumor"]["n_genes"] == 10


def test_parsers_probe_counts_edges_and_unique_genes(tmp_path):
    kegg = tmp_path / "kegg_interactions.txt"
    kegg.write_text("source\ttarget\tdirected\nA\tB\ttrue\nB\tC\ttrue\n")
    bgr = tmp_path / "biogrid_interactions.txt"
    bgr.write_text("source\ttarget\nA\tD\nC\tD\nA\tD\n")  # one duplicate

    r = ParsersProbe().run(ctx={
        "kegg_interactions": str(kegg),
        "biogrid_interactions": str(bgr),
    })
    m = r.measurements
    assert m["kegg"]["n_edges"] == 2
    assert m["kegg"]["unique_genes"] == 3
    assert m["biogrid"]["n_edges"] == 3
    assert m["biogrid"]["duplicate_edges"] == 1
    assert m["biogrid"]["unique_genes"] == 3


from src.python.diagnostics.probes.merge import MergeProbe


def test_merge_probe_counts_directed_undirected(tmp_path):
    merged = tmp_path / "merged.txt"
    merged.write_text(
        "source\ttarget\tdirected\n"
        "A\tB\ttrue\n"
        "B\tA\ttrue\n"        # forms a 2-cycle but each entry is directed
        "C\tD\tfalse\n"
        "D\tE\tfalse\n"
    )
    r = MergeProbe().run(ctx={"merged_network": str(merged)})
    m = r.measurements
    assert m["n_edges"] == 4
    assert m["n_directed"] == 2
    assert m["n_undirected"] == 2
    assert m["unique_genes"] == 5
    assert m["bidirectional_pairs"] == 1   # A-B has both directions


def test_filters_probe_reports_dea_ig_intersection(tmp_path):
    dea_table = tmp_path / "dea.tsv"
    dea_table.write_text(
        "gene\tlog2FC\tpvalue\tFDR\n"
        "G1\t2.0\t0.001\t0.01\n"
        "G2\t-2.0\t0.001\t0.01\n"
        "G3\t0.5\t0.5\t0.5\n"
        "G4\t3.0\t0.001\t0.01\n"
    )
    dea_genes = tmp_path / "dea_genes.txt"
    dea_genes.write_text("G1\nG2\nG4\n")
    ig_genes = tmp_path / "ig_genes.txt"
    ig_genes.write_text("G2\nG4\nG5\n")
    inter = tmp_path / "candidates.txt"
    inter.write_text("G2\nG4\n")

    r = FiltersProbe().run(ctx={
        "dea_table": str(dea_table),
        "dea_genes": str(dea_genes),
        "ig_genes": str(ig_genes),
        "candidates": str(inter),
    })
    m = r.measurements
    assert m["dea"]["n_genes_total"] == 4
    assert m["dea"]["n_significant"] == 3
    assert m["dea"]["fdr_p10"] <= m["dea"]["fdr_p90"]
    assert m["ig"]["n_genes"] == 3
    assert m["intersection"]["n_genes"] == 2
    assert m["intersection"]["jaccard_dea_ig"] == pytest.approx(0.5)


from src.python.diagnostics.probes.weighting import WeightingProbe


def test_weighting_probe_explains_edge_loss(tmp_path):
    merged = tmp_path / "merged.txt"
    merged.write_text(
        "source\ttarget\n"
        "A\tB\n"
        "A\tZZUNKNOWN\n"   # dropped: target not in expression matrix
        "C\tD\n"
        "E\tF\n"
    )
    weighted = tmp_path / "weighted.txt"
    weighted.write_text(
        "source\ttarget\tweight\tweight_int\n"
        "A\tB\t0.05\t50\n"
        "C\tD\t0.000001\t1\n"     # epsilon-clamp survivor
        "E\tF\t0.5\t500\n"
    )
    expr = tmp_path / "expr.tsv"
    expr.write_text("gene\tS1\tS2\nA\t1\t2\nB\t2\t1\nC\t1\t1\nD\t1\t1\nE\t1\t2\nF\t2\t1\n")

    r = WeightingProbe().run(ctx={
        "merged_network": str(merged),
        "weighted_network": str(weighted),
        "expression": str(expr),
        "epsilon": 1e-6,
        "scale_factor": 1000,
    })
    m = r.measurements
    assert m["n_input_edges"] == 4
    assert m["n_output_edges"] == 3
    assert m["n_dropped"] == 1
    assert m["edges_unresolved_in_expression"] == 1
    assert m["epsilon_clamped_edges"] >= 1
    # merged net has 7 unique genes (A,B,ZZUNKNOWN,C,D,E,F); expr has 6 -> 6/7
    assert m["gene_symbol_resolution_rate"] == pytest.approx(6/7)


import json
from src.python.diagnostics.probes.centrality import CentralityProbe


def test_centrality_probe_computes_spearman(tmp_path):
    diag = tmp_path / "cpp_diagnostics.json"
    diag.write_text(json.dumps({
        "omp_determinism_max_relative_drift": 0.0,
        "sampled_vs_full_subgraph_n": 5,
        "sampled_full_betweenness": [0.0, 1.0, 2.0, 3.0, 4.0],
        "sampled_betweenness":      [0.1, 1.1, 2.1, 3.1, 4.1],
    }))
    r = CentralityProbe().run(ctx={"cpp_diagnostics_json": str(diag)})
    m = r.measurements
    assert m["omp_determinism_max_relative_drift"] == 0.0
    assert m["sampled_vs_full_spearman_rho"] == pytest.approx(1.0)


from pathlib import Path
from src.python.diagnostics.probes.differential import DifferentialProbe

CGC = Path(__file__).parents[1] / "fixtures" / "cgc_mini.txt"


def test_differential_probe_finds_hub_kills(tmp_path):
    centrality_normal = tmp_path / "cn.txt"
    centrality_normal.write_text(
        "gene\tbetweenness\teigenvector\n"
        + "\n".join(f"G{i}\t{i}\t{0.01*i}" for i in range(100)) + "\n"
        + "TP53\t10\t0.99\nEP300\t12\t0.95\n"
    )
    centrality_tumor = tmp_path / "ct.txt"
    centrality_tumor.write_text(
        "gene\tbetweenness\teigenvector\n"
        + "\n".join(f"G{i}\t{i+1}\t{0.01*i}" for i in range(100)) + "\n"
        + "TP53\t9\t0.97\nEP300\t11\t0.93\n"
    )
    r = DifferentialProbe().run(ctx={
        "centrality_normal": str(centrality_normal),
        "centrality_tumor": str(centrality_tumor),
        "cgc_list": str(CGC),
    })
    m = r.measurements
    # bet/eig rank correlation defined and bounded
    assert -1.0 <= m["bet_eig_rank_correlation"] <= 1.0
    # TP53 + EP300 are in the top-12% of eigenvector and in CGC
    assert "TP53" in m["top12_eig_dropped_genes"]
    assert "EP300" in m["top12_eig_dropped_genes"]
    assert m["top12_eig_dropped_in_cgc"] >= 2


from src.python.diagnostics.probes.validation import ValidationProbe


def test_validation_probe_reports_overlap_and_join_rate(tmp_path):
    dice = tmp_path / "dice.tsv"
    dice.write_text("gene\nG1\nG2\nG3\nG4\nG5\n")
    deg  = tmp_path / "deg.tsv"
    deg.write_text("gene\nG3\nG4\nG6\n")
    depmap = tmp_path / "depmap.tsv"
    depmap.write_text("gene\tdependency_score\nG1\t-0.5\nG3\t-0.7\nG7\t0.0\n")

    r = ValidationProbe().run(ctx={
        "dice_genes": str(dice),
        "deg_genes": str(deg),
        "depmap_table": str(depmap),
    })
    m = r.measurements
    assert m["dice_deg_overlap"]["observed"] == 2
    assert 0.0 < m["dice_deg_overlap"]["hypergeom_p"] <= 1.0
    assert m["depmap"]["join_match_count"] == 2
    assert m["depmap"]["join_match_rate"] == pytest.approx(0.4)
