"""Base class for pipeline boundary probes. Pure observation; never mutates."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, Optional
import logging
import time

log = logging.getLogger(__name__)


@dataclass
class ProbeResult:
    name: str
    stage: str
    measurements: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None
    duration_s: float = 0.0


class Probe:
    name: str = "unnamed"
    stage: str = "unknown"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        raise NotImplementedError

    def run(self, ctx: Dict[str, Any]) -> ProbeResult:
        t0 = time.perf_counter()
        try:
            m = self.measure(ctx)
            return ProbeResult(name=self.name, stage=self.stage, measurements=m,
                               duration_s=time.perf_counter() - t0)
        except Exception as e:  # noqa: BLE001 — probes must never crash the runner
            log.exception(f"probe {self.name} crashed")
            return ProbeResult(
                name=self.name,
                stage=self.stage,
                error=f"{type(e).__name__}: {e}",
                duration_s=time.perf_counter() - t0,
            )
