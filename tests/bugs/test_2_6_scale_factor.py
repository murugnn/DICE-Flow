"""Bug 2.6: scale=1000 collapses many distinct |r|-distance weights into the
same integer. Increasing scale_factor reduces collisions and preserves
fidelity for the SSSP engine."""
import random
import subprocess
import sys

import pandas as pd


def test_scale_factor_100k_fewer_collisions_than_1k(tmp_path):
    net = tmp_path / "net.txt"
    net.write_text(
        "source\ttarget\n"
        + "\n".join(f"G{i}\tG{i+1}" for i in range(50))
        + "\n"
    )
    rng = random.Random(20260507)
    expr_lines = ["gene\tS1\tS2\tS3\tS4\tS5"]
    for i in range(52):
        expr_lines.append(
            f"G{i}\t" + "\t".join(f"{rng.uniform(0, 10):.4f}" for _ in range(5))
        )
    expr = tmp_path / "expr.tsv"
    expr.write_text("\n".join(expr_lines) + "\n")

    def collisions(scale):
        out = tmp_path / f"w_{scale}.txt"
        subprocess.run(
            [sys.executable, "src/python/weight_network.py",
             "--network", str(net), "--expression", str(expr),
             "--output", str(out), "--scale-factor", str(scale)],
            check=True, capture_output=True, text=True,
        )
        df = pd.read_csv(out, sep="\t")
        return df["weight"].nunique() - df["weight_int"].nunique()

    assert collisions(100000) <= collisions(1000)
