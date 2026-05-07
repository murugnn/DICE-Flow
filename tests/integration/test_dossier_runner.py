import json
import subprocess
import sys
from pathlib import Path


def test_dossier_runner_writes_files(tmp_path, tiny_expr_normal, tiny_expr_tumor):
    norm = tmp_path / "norm.tsv"; tiny_expr_normal.to_csv(norm, sep="\t")
    tum  = tmp_path / "tum.tsv";  tiny_expr_tumor.to_csv(tum, sep="\t")
    out  = tmp_path / "out"

    r = subprocess.run([
        sys.executable, "-m", "src.python.diagnostics.dossier_runner",
        "--expression-normal", str(norm),
        "--expression-tumor", str(tum),
        "--output", str(out),
    ], capture_output=True, text=True, cwd=Path(__file__).parents[2])
    assert r.returncode == 0, r.stderr

    md = (out / "diagnostic_dossier.md").read_text()
    js = json.loads((out / "dossier.json").read_text())
    assert "Diagnostic Dossier" in md
    assert any(p["name"] == "raw_inputs" for p in js["probes"])
