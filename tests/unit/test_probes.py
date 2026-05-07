import pandas as pd
from src.python.diagnostics.probes.raw_inputs import RawInputsProbe


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
