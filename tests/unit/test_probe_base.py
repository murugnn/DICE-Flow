from src.python.diagnostics.probe_base import Probe, ProbeResult


class _TinyProbe(Probe):
    name = "tiny"
    stage = "stage-0"

    def measure(self, ctx):
        return {"answer": 42}


def test_probe_yields_probe_result_with_metadata():
    p = _TinyProbe()
    r = p.run(ctx={"foo": "bar"})
    assert isinstance(r, ProbeResult)
    assert r.name == "tiny"
    assert r.stage == "stage-0"
    assert r.measurements == {"answer": 42}
    assert r.error is None


def test_probe_captures_exceptions_into_result():
    class _Boom(Probe):
        name = "boom"
        stage = "explosive"

        def measure(self, ctx):
            raise RuntimeError("nope")

    r = _Boom().run(ctx={})
    assert r.error == "RuntimeError: nope"
    assert r.measurements == {}
