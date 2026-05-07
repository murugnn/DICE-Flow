from src.python.diagnostics.probe_base import ProbeResult
from src.python.diagnostics.render import render_markdown, render_json


def _make_results():
    return [
        ProbeResult(name="raw_inputs", stage="1-raw-inputs",
                    measurements={"a": 1, "b": {"c": 2}}, duration_s=0.01),
        ProbeResult(name="weighting", stage="5-weighting",
                    measurements={"n_dropped": 100}, duration_s=0.02),
    ]


def test_render_markdown_emits_per_stage_sections():
    md = render_markdown(_make_results(), commit="abc123")
    assert "# Diagnostic Dossier" in md
    assert "## Stage 1-raw-inputs" in md
    assert "## Stage 5-weighting" in md
    assert "n_dropped" in md
    assert "abc123" in md


def test_render_json_round_trip():
    import json
    s = render_json(_make_results(), commit="abc123")
    d = json.loads(s)
    assert d["commit"] == "abc123"
    assert len(d["probes"]) == 2
    assert d["probes"][0]["name"] == "raw_inputs"
