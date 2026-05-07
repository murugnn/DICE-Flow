import logging
import pandas as pd
import pytest
from src.python.utils.safe_filter import safe_filter, FilterResult


def test_safe_filter_records_in_out_counts(caplog):
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5]})
    mask = df["x"] > 2

    with caplog.at_level(logging.INFO):
        result = safe_filter(df, mask, reason="x > 2", stage="unit-test")

    assert isinstance(result, FilterResult)
    assert len(result.df) == 3
    assert result.in_count == 5
    assert result.out_count == 3
    assert result.dropped_count == 2
    assert result.reason == "x > 2"
    assert any("safe_filter[unit-test]" in r.message for r in caplog.records)


def test_safe_filter_raises_on_all_dropped_when_strict():
    df = pd.DataFrame({"x": [1, 2, 3]})
    mask = df["x"] > 99
    with pytest.raises(ValueError):
        safe_filter(df, mask, reason="impossible", stage="unit-test", strict=True)


def test_safe_filter_warns_on_all_dropped_by_default(caplog):
    df = pd.DataFrame({"x": [1, 2, 3]})
    mask = df["x"] > 99
    with caplog.at_level(logging.WARNING):
        safe_filter(df, mask, reason="impossible", stage="unit-test")
    assert any("dropped 3/3" in r.message for r in caplog.records)
