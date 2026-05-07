"""Bug 2.4: DiCE 1-|r| weighting drops the sign. We continue using |r|-distance
for SSSP, but additionally emit (source, target, signed_correlation) so
downstream analysis can distinguish coactivation from antagonism."""
import pathlib
import subprocess
import sys

import pandas as pd


def test_signed_correlation_file_emitted(tmp_path):
    net = tmp_path / "net.txt"
    net.write_text("source\ttarget\nA\tB\nA\tC\n")
    expr = tmp_path / "expr.tsv"
    expr.write_text(
        "gene\tS1\tS2\tS3\tS4\tS5\n"
        "A\t1\t2\t3\t4\t5\n"
        "B\t5\t4\t3\t2\t1\n"
        "C\t1\t2\t3\t4\t5\n"
    )
    out = tmp_path / "weighted.txt"
    r = subprocess.run([
        sys.executable, "src/python/weight_network.py",
        "--network", str(net), "--expression", str(expr),
        "--output", str(out),
    ], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr

    signed = pathlib.Path(str(out).replace(".txt", "_signed.tsv"))
    assert signed.exists(), f"expected signed-correlation file at {signed}"
    df = pd.read_csv(signed, sep="\t")
    # A-B is anti-correlated (~ -1), A-C is positively correlated (~ +1)
    ab = df[(df["source"] == "A") & (df["target"] == "B")]["signed_correlation"].iloc[0]
    ac = df[(df["source"] == "A") & (df["target"] == "C")]["signed_correlation"].iloc[0]
    assert ab < -0.9, f"expected A-B anti-correlated, got {ab}"
    assert ac > 0.9, f"expected A-C positively correlated, got {ac}"
