"""Logged row-filter wrapper. No silent drops."""
from __future__ import annotations
from dataclasses import dataclass
import logging

import pandas as pd

log = logging.getLogger(__name__)


@dataclass
class FilterResult:
    df: pd.DataFrame
    in_count: int
    out_count: int
    dropped_count: int
    reason: str
    stage: str


def safe_filter(
    df: pd.DataFrame,
    mask: pd.Series,
    reason: str,
    stage: str,
    strict: bool = False,
) -> FilterResult:
    """Apply mask, log counts, return FilterResult. No silent drops."""
    in_count = len(df)
    filtered = df[mask].copy()
    out_count = len(filtered)
    dropped = in_count - out_count

    msg = f"safe_filter[{stage}] reason='{reason}' in={in_count} out={out_count} dropped={dropped}"
    if dropped == 0:
        log.info(msg)
    elif out_count == 0:
        if strict:
            raise ValueError(f"safe_filter[{stage}] dropped {dropped}/{in_count} rows ({reason})")
        log.warning(f"{msg} (dropped {dropped}/{in_count})")
    else:
        log.info(msg)

    return FilterResult(
        df=filtered,
        in_count=in_count,
        out_count=out_count,
        dropped_count=dropped,
        reason=reason,
        stage=stage,
    )
