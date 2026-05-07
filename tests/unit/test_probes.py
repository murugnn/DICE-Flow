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
