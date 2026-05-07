from pathlib import Path
import pytest
from src.python.utils.id_resolver import IDResolver, ResolveResult

FIX = Path(__file__).parents[1] / "fixtures" / "hgnc_mini.tsv"


@pytest.fixture(scope="module")
def resolver():
    return IDResolver.from_hgnc(FIX)


def test_resolves_canonical_symbol(resolver):
    r = resolver.resolve("TP53")
    assert r.symbol == "TP53"
    assert r.entrez_id == "7157"
    assert r.match_type == "symbol"


def test_resolves_prev_symbol(resolver):
    r = resolver.resolve("P53")
    assert r.symbol == "TP53"
    assert r.match_type == "prev_symbol"


def test_resolves_alias(resolver):
    r = resolver.resolve("c-Myc")
    assert r.symbol == "MYC"
    assert r.match_type == "alias"


def test_resolves_entrez(resolver):
    r = resolver.resolve("7157")
    assert r.symbol == "TP53"
    assert r.match_type == "entrez"


def test_unresolved_returns_none(resolver):
    r = resolver.resolve("ZZZNOTAGENE")
    assert r is None


def test_batch_resolution_reports_miss_rate(resolver):
    syms = ["TP53", "P53", "ZZZNOTAGENE", "MYC"]
    rs = resolver.resolve_batch(syms)
    assert sum(1 for r in rs if r is not None) == 3
    assert resolver.last_miss_rate == pytest.approx(0.25)
