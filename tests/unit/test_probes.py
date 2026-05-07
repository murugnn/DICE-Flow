import pandas as pd
from src.python.diagnostics.probes.raw_inputs import RawInputsProbe
from src.python.diagnostics.probes.parsers import ParsersProbe


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
