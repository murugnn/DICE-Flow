"""Bug 2.2: top-12% hub removal silently kills TP53/EP300.
Default must keep them (correct-by-default); --dice-faithful drops them."""
import subprocess
import sys

import pandas as pd


def _make_fake_centralities(tmp_path):
    cn = tmp_path / "cn.txt"
    ct = tmp_path / "ct.txt"
    rows_n = ["gene\tbetweenness\teigenvector"]
    rows_t = ["gene\tbetweenness\teigenvector"]
    for i in range(100):
        rows_n.append(f"G{i}\t{i}\t{0.005*i}")
        rows_t.append(f"G{i}\t{i+1}\t{0.005*i}")
    rows_n.append("TP53\t10\t0.99")
    rows_t.append("TP53\t9\t0.97")
    cn.write_text("\n".join(rows_n) + "\n")
    ct.write_text("\n".join(rows_t) + "\n")
    return cn, ct


def test_default_keeps_tp53(tmp_path):
    cn, ct = _make_fake_centralities(tmp_path)
    out = tmp_path / "dice.txt"
    r = subprocess.run([
        sys.executable, "src/python/differential_centrality.py",
        "--normal", str(cn), "--tumor", str(ct), "--output", str(out),
        "--top-n", "200",
    ], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    df = pd.read_csv(out, sep="\t")
    assert "TP53" in set(df["gene"])
    assert "ensemble_score_raw" in df.columns
    assert "ensemble_score_hub_corrected" in df.columns


def test_dice_faithful_drops_tp53_when_top12_eig(tmp_path):
    cn, ct = _make_fake_centralities(tmp_path)
    out = tmp_path / "dice.txt"
    r = subprocess.run([
        sys.executable, "src/python/differential_centrality.py",
        "--normal", str(cn), "--tumor", str(ct), "--output", str(out),
        "--top-n", "200", "--dice-faithful",
    ], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    df = pd.read_csv(out, sep="\t")
    assert "TP53" not in set(df["gene"])
    assert "ensemble_score" in df.columns
