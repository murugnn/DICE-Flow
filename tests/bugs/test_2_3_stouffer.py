"""Bug 2.3: rank-product double-counts because betweenness and eigenvector
are correlated. Add Stouffer combined score (rank-percentiles -> z-scores -> sum/sqrt(2))
as a decorrelated alternative."""
import subprocess
import sys

import pandas as pd


def test_stouffer_score_present(tmp_path):
    cn = tmp_path / "cn.txt"
    ct = tmp_path / "ct.txt"
    rows_n = ["gene\tbetweenness\teigenvector"]
    rows_t = ["gene\tbetweenness\teigenvector"]
    for i in range(50):
        # Vary the deltas so rank-percentiles span (0, 1) and Stouffer z-scores
        # cover both signs — otherwise constant input collapses the Stouffer
        # column to a single value and the sign assertions can't be met.
        rows_n.append(f"G{i}\t{i}\t{0.01*i}")
        rows_t.append(f"G{i}\t{i*2}\t{0.005*i}")
    cn.write_text("\n".join(rows_n) + "\n")
    ct.write_text("\n".join(rows_t) + "\n")
    out = tmp_path / "dice.txt"
    r = subprocess.run([
        sys.executable, "src/python/differential_centrality.py",
        "--normal", str(cn), "--tumor", str(ct), "--output", str(out),
        "--top-n", "100",
    ], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    df = pd.read_csv(out, sep="\t")
    assert "ensemble_score_stouffer" in df.columns
    # Stouffer is a z-score; should span both signs across the full set
    assert df["ensemble_score_stouffer"].max() > 0
    assert df["ensemble_score_stouffer"].min() < 0
