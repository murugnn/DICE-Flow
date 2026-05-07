"""HGNC-based gene symbol resolver. Single canonical authority for ID joins."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Iterable, List
import logging

import pandas as pd

log = logging.getLogger(__name__)


@dataclass
class ResolveResult:
    symbol: str
    entrez_id: Optional[str]
    match_type: str  # "symbol" | "prev_symbol" | "alias" | "entrez"


class IDResolver:
    def __init__(
        self,
        symbol_to_canon: dict,
        prev_to_canon: dict,
        alias_to_canon: dict,
        entrez_to_canon: dict,
        canon_to_entrez: dict,
    ):
        self._sym = symbol_to_canon
        self._prev = prev_to_canon
        self._alias = alias_to_canon
        self._entrez = entrez_to_canon
        self._canon_entrez = canon_to_entrez
        self.last_miss_rate: float = 0.0

    @classmethod
    def from_hgnc(cls, path: Path) -> "IDResolver":
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        sym = {row["symbol"]: row["symbol"] for _, row in df.iterrows()}
        canon_entrez = {row["symbol"]: row["entrez_id"] or None for _, row in df.iterrows()}
        prev = {}
        alias = {}
        entrez = {}
        for _, row in df.iterrows():
            for p in filter(None, row["prev_symbol"].split("|")):
                prev.setdefault(p, row["symbol"])
            for a in filter(None, row["alias_symbol"].split("|")):
                alias.setdefault(a, row["symbol"])
            if row["entrez_id"]:
                entrez[row["entrez_id"]] = row["symbol"]
        return cls(sym, prev, alias, entrez, canon_entrez)

    def resolve(self, key: str) -> Optional[ResolveResult]:
        if key in self._sym:
            return ResolveResult(self._sym[key], self._canon_entrez.get(key), "symbol")
        if key in self._prev:
            canon = self._prev[key]
            return ResolveResult(canon, self._canon_entrez.get(canon), "prev_symbol")
        if key in self._alias:
            canon = self._alias[key]
            return ResolveResult(canon, self._canon_entrez.get(canon), "alias")
        if key in self._entrez:
            canon = self._entrez[key]
            return ResolveResult(canon, key, "entrez")
        return None

    def resolve_batch(self, keys: Iterable[str]) -> List[Optional[ResolveResult]]:
        keys = list(keys)
        results = [self.resolve(k) for k in keys]
        if keys:
            misses = sum(1 for r in results if r is None)
            self.last_miss_rate = misses / len(keys)
            log.info(
                f"id_resolver: resolved {len(keys) - misses}/{len(keys)} "
                f"(miss_rate={self.last_miss_rate:.3f})"
            )
        return results
