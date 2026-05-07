# DiCE-Duan Audit & Fix — Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the testing infrastructure, produce the diagnostic dossier (10 boundary probes), and land the 5 pre-committed correct-by-default changes from the spec, ending with a full TCGA-PRAD re-run + before/after delta report.

**Architecture:** Two-phase reproduce-first audit. Phase 0 stands up testing infrastructure + utilities (`safe_filter`, `id_resolver`) and captures pre-audit golden artifacts. Phase 1 implements one probe per pipeline boundary as a pure-observation `Probe` subclass; a `dossier_runner` orchestrates them into `docs/audit/diagnostic_dossier.md` + `dossier.json`. Phase 2-precommitted lands the five DiCE-faithful-flagged corrections (top-12 % hub removal off, Stouffer score, signed-correlation output, Brandes full-or-bootstrap, scale ×100000). Phase 3 runs the final gate and emits `before_after.md`. Phase 2-dossier-driven fixes (DEG 0-overlap, etc.) are deferred to a follow-up plan generated *from* the registry produced here, because writing TDD code for them before the dossier exists would be guessing.

**Tech Stack:** Python 3.8+ (pandas, numpy, scipy, scikit-learn, networkx), C++17 (existing engine.hpp), pytest, doctest (header-only C++ test framework), CMake.

**Spec:** `docs/superpowers/specs/2026-05-07-audit-and-fix-design.md`. Read it before starting.

---

## File structure (locked at plan time)

### New files
| Path | Responsibility |
|---|---|
| `tests/conftest.py` | Shared pytest fixtures: tiny synthetic graph, expression matrices, golden hashes |
| `tests/unit/test_safe_filter.py` | Unit tests for the safe-filter wrapper |
| `tests/unit/test_id_resolver.py` | Unit tests for the HGNC-based resolver |
| `tests/unit/test_probes.py` | Unit tests for each probe on synthetic input |
| `tests/integration/test_micro_pipeline.py` | 50-node × 200-edge end-to-end pipeline test |
| `tests/regression/test_dice_faithful_regression.py` | `--dice-faithful` ≥ 0.99 rank-correlation contract |
| `tests/fixtures/tiny_graph.tsv` | 10-node hand-built directed graph for centrality tests |
| `tests/fixtures/tiny_expr_normal.tsv` | 10-gene × 5-sample synthetic expression |
| `tests/fixtures/tiny_expr_tumor.tsv` | 10-gene × 5-sample synthetic expression |
| `src/python/utils/__init__.py` | (empty package marker) |
| `src/python/utils/safe_filter.py` | `safe_filter(df, mask, reason)` — logs in/out counts |
| `src/python/utils/id_resolver.py` | HGNC-based gene-symbol resolver with explicit miss handling |
| `src/python/diagnostics/__init__.py` | (empty package marker) |
| `src/python/diagnostics/probe_base.py` | `Probe` base class + `ProbeResult` dataclass |
| `src/python/diagnostics/probes/raw_inputs.py` | Probe stage 1 |
| `src/python/diagnostics/probes/parsers.py` | Probe stage 2 |
| `src/python/diagnostics/probes/merge.py` | Probe stage 3 |
| `src/python/diagnostics/probes/filters.py` | Probe stage 4 (DEA + IG + intersection) |
| `src/python/diagnostics/probes/weighting.py` | Probe stage 5 (the 73 % drop investigation) |
| `src/python/diagnostics/probes/centrality.py` | Probe stage 6/7 (consume C++ `--diagnostics` output) |
| `src/python/diagnostics/probes/differential.py` | Probe stage 7 (rank-product / hub kill list) |
| `src/python/diagnostics/probes/validation.py` | Probe stage 8 (DEG overlap / DepMap join / survival cohort) |
| `src/python/diagnostics/probes/dnf.py` | Probe stage 10 |
| `src/python/diagnostics/dossier_runner.py` | Orchestrator + CLI entry point |
| `src/python/diagnostics/render.py` | Markdown + JSON renderer |
| `src/cpp/tests/test_engine.cpp` | doctest C++ tests (Brandes / eigenvector / OMP determinism on toy graphs) |
| `scripts/run_tests.sh` | One driver: pytest + C++ doctest |
| `scripts/capture_golden.sh` | One-shot: snapshot pre-audit artifacts before any code change |
| `scripts/run_audit_rerun.sh` | Final-gate full PRAD re-run + invokes `before_after.py` |
| `scripts/before_after.py` | Generates `docs/audit/before_after.md` from old vs new artifacts |
| `docs/audit/diagnostic_dossier.md` | (generated, committed at end of Phase 1) |
| `docs/audit/dossier.json` | (generated, committed at end of Phase 1) |
| `docs/audit/before_after.md` | (generated, committed at end of Phase 3) |

### Existing files modified
| Path | Change |
|---|---|
| `src/python/differential_centrality.py` | Add `--dice-faithful`; default off for top-12 % hub removal; emit raw + hub-corrected ensemble; add Stouffer score |
| `src/python/weight_network.py` | Add `--scale-factor`; emit signed-correlation matrix; replace silent NaN→0 with logged event via `safe_filter`; replace `weight_int==0 → 1` collapse with logged event |
| `src/python/dea_filter.py` | Replace silent "top-1000 fallback" (line 96-99) with logged event via `safe_filter` |
| `src/cpp/engine.hpp` | Brandes: configurable sampling rate; full when `n_vertices < 5000`; bootstrap CI when sampled |
| `src/cpp/main.cpp` | Add CLI flags `--diagnostics`, `--brandes-mode {full\|sampled\|bootstrap}`, `--scale-factor` |
| `run_pipeline.sh` | Propagate `--dice-faithful`, `--scale-factor`, `--brandes-mode`; pin random seeds |
| `CMakeLists.txt` | Add `test_engine` doctest target |

---

## Conventions for every task

- **Branch & commits.** Work on `audit-phase1` branch. One commit per task; commit message format: `<type>(<scope>): <subject> [task-N.M]`.
- **TDD discipline.** Test first, fail explicitly, implement, pass, commit.
- **No silent drops.** Every row-removal goes through `safe_filter`. Every ID join goes through `id_resolver`.
- **Determinism.** Every probe and test seeds RNGs (numpy & C++ `mt19937`) with seed = 20260507.
- **Logs.** All Python probes use the existing `logging` module; format `[%(levelname)s] %(name)s: %(message)s`. C++ uses stderr.

---

## PHASE 0 — Bootstrap (5 tasks)

### Task 0.1: Stand up the test scaffolding

**Files:**
- Create: `tests/conftest.py`
- Create: `tests/unit/__init__.py`, `tests/integration/__init__.py`, `tests/regression/__init__.py`, `tests/bugs/__init__.py`
- Create: `tests/fixtures/tiny_graph.tsv`
- Create: `tests/fixtures/tiny_expr_normal.tsv`
- Create: `tests/fixtures/tiny_expr_tumor.tsv`
- Create: `scripts/run_tests.sh`
- Create: `pytest.ini`

- [ ] **Step 1: Create the synthetic fixtures**

`tests/fixtures/tiny_graph.tsv` (10 nodes A..J, directed, 12 edges, hand-checked Brandes centrality):
```
source	target	weight_int
A	B	100
A	C	200
B	D	100
C	D	100
D	E	100
E	F	100
E	G	100
F	H	100
G	H	100
H	I	100
H	J	100
I	J	100
```

`tests/fixtures/tiny_expr_normal.tsv`:
```
gene	S1	S2	S3	S4	S5
A	1.0	1.1	0.9	1.2	0.8
B	2.0	2.1	1.9	2.2	1.8
C	0.5	0.6	0.4	0.7	0.3
D	5.0	5.1	4.9	5.2	4.8
E	3.0	3.1	2.9	3.2	2.8
F	4.0	4.1	3.9	4.2	3.8
G	1.5	1.6	1.4	1.7	1.3
H	6.0	6.1	5.9	6.2	5.8
I	0.1	0.2	0.0	0.3	0.0
J	2.5	2.6	2.4	2.7	2.3
```

`tests/fixtures/tiny_expr_tumor.tsv` — same shape; gene D and H show ~2× change to make a couple of DEGs:
```
gene	S1	S2	S3	S4	S5
A	1.1	1.0	1.2	0.9	1.1
B	2.1	1.9	2.0	2.1	2.2
C	0.5	0.7	0.4	0.6	0.5
D	10.0	10.2	9.9	10.4	9.7
E	3.1	2.9	3.0	3.2	3.0
F	4.0	4.2	4.1	3.9	4.0
G	1.5	1.4	1.6	1.5	1.5
H	12.0	11.9	12.1	12.2	11.8
I	0.1	0.0	0.2	0.0	0.1
J	2.6	2.4	2.5	2.7	2.5
```

- [ ] **Step 2: Write conftest with shared fixtures**

`tests/conftest.py`:
```python
from pathlib import Path
import pandas as pd
import pytest

FIXTURES = Path(__file__).parent / "fixtures"
SEED = 20260507


@pytest.fixture(scope="session")
def tiny_graph_path():
    return FIXTURES / "tiny_graph.tsv"


@pytest.fixture(scope="session")
def tiny_expr_normal():
    return pd.read_csv(FIXTURES / "tiny_expr_normal.tsv", sep="\t", index_col=0)


@pytest.fixture(scope="session")
def tiny_expr_tumor():
    return pd.read_csv(FIXTURES / "tiny_expr_tumor.tsv", sep="\t", index_col=0)


@pytest.fixture(autouse=True)
def _seed():
    import random
    import numpy as np
    random.seed(SEED)
    np.random.seed(SEED)
```

- [ ] **Step 3: Write `pytest.ini`**

```
[pytest]
testpaths = tests
addopts = -ra -q --strict-markers
markers =
    integration: end-to-end micro-pipeline tests
    regression: dice-faithful regression contract
    bugs: per-bug-id reproductions
```

- [ ] **Step 4: Write the run_tests script**

`scripts/run_tests.sh`:
```bash
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== Python tests =="
python -m pytest tests/unit tests/integration tests/regression tests/bugs

echo "== C++ tests =="
if [[ -x build/test_engine ]]; then
    ./build/test_engine
else
    echo "build/test_engine not built yet — skipping (run cmake --build build first)"
fi
```

Then `chmod +x scripts/run_tests.sh`.

- [ ] **Step 5: Verify pytest discovers nothing yet**

Run: `cd /home/cerelac/DICE && python -m pytest --collect-only`
Expected: `no tests ran`. (Empty test directories, no test files yet.)

- [ ] **Step 6: Commit**

```bash
git add tests/ scripts/run_tests.sh pytest.ini
git commit -m "chore(audit): scaffold pytest dirs, conftest, tiny fixtures, run_tests.sh [task-0.1]"
```

---

### Task 0.2: Wire C++ doctest

**Files:**
- Create: `src/cpp/tests/test_engine.cpp`
- Create: `src/cpp/tests/doctest.h` (vendor single-header)
- Modify: `CMakeLists.txt` (append doctest target)

- [ ] **Step 1: Vendor doctest single-header**

Run:
```bash
mkdir -p src/cpp/tests
curl -fsSL -o src/cpp/tests/doctest.h \
  https://raw.githubusercontent.com/doctest/doctest/v2.4.11/doctest/doctest.h
```

(If offline, paste from any cached copy. Single header, ~12 K lines, MIT.)

- [ ] **Step 2: Write the failing test**

`src/cpp/tests/test_engine.cpp`:
```cpp
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../engine.hpp"

// Tiny directed graph: A -> B -> C, edges weight 1.
// Hand-computed Brandes betweenness (unnormalized): A=0, B=1, C=0.
TEST_CASE("Brandes betweenness on a 3-node path graph") {
    Graph<long long> g(3);
    g.id_to_name = {"A", "B", "C"};
    g.name_to_id = {{"A", 0}, {"B", 1}, {"C", 2}};
    g.add_edge(0, 1, 1);
    g.add_edge(1, 2, 1);

    CentralityCalculator<long long> calc;
    auto bet = calc.betweenness(g, /*sample_rate=*/1.0);

    CHECK(bet[0] == doctest::Approx(0.0));
    CHECK(bet[1] == doctest::Approx(1.0));
    CHECK(bet[2] == doctest::Approx(0.0));
}
```

- [ ] **Step 3: Append doctest target to `CMakeLists.txt`**

Append to the bottom of `CMakeLists.txt`:
```cmake
# ============================================================================
# Tests
# ============================================================================
enable_testing()
add_executable(test_engine src/cpp/tests/test_engine.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_engine PUBLIC OpenMP::OpenMP_CXX)
endif()
add_test(NAME engine_tests COMMAND test_engine)
```

- [ ] **Step 4: Build and run**

```bash
cd /home/cerelac/DICE
cmake --build build --target test_engine
./build/test_engine
```

Expected: test PASSES if `CentralityCalculator::betweenness(graph, sample_rate)` already exists with the right signature; **FAILS to compile** if it doesn't (current engine likely has different API). If compile fails, **do not** modify `engine.hpp` yet — instead, change the test signature in `test_engine.cpp` to match whatever the current `CentralityCalculator` API actually is, and verify by checking `engine.hpp` for the public method name. Re-run until it passes on the existing code path. *We are baselining current behavior, not changing it yet.*

- [ ] **Step 5: Verify pass**

Re-run; expect `test_engine` passes.

- [ ] **Step 6: Commit**

```bash
git add src/cpp/tests/ CMakeLists.txt
git commit -m "chore(audit): add doctest C++ test target with toy directed graph [task-0.2]"
```

---

### Task 0.3: Build the `safe_filter` utility (TDD)

**Files:**
- Create: `src/python/utils/__init__.py` (empty)
- Create: `src/python/utils/safe_filter.py`
- Create: `tests/unit/test_safe_filter.py`

- [ ] **Step 1: Write the failing test**

`tests/unit/test_safe_filter.py`:
```python
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
```

- [ ] **Step 2: Run test, expect failure**

Run: `python -m pytest tests/unit/test_safe_filter.py -v`
Expected: `ModuleNotFoundError: No module named 'src.python.utils.safe_filter'`

- [ ] **Step 3: Implement minimal `safe_filter`**

`src/python/utils/__init__.py`:
```python
```
(empty)

`src/python/utils/safe_filter.py`:
```python
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
```

- [ ] **Step 4: Run test, expect pass**

Run: `python -m pytest tests/unit/test_safe_filter.py -v`
Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add src/python/utils/__init__.py src/python/utils/safe_filter.py tests/unit/test_safe_filter.py
git commit -m "feat(utils): add safe_filter with logged in/out counts [task-0.3]"
```

---

### Task 0.4: Build the `id_resolver` utility (TDD)

**Files:**
- Create: `src/python/utils/id_resolver.py`
- Create: `tests/unit/test_id_resolver.py`
- Create: `tests/fixtures/hgnc_mini.tsv`

- [ ] **Step 1: Create the mini HGNC fixture**

`tests/fixtures/hgnc_mini.tsv` (subset for tests; production resolver loads full HGNC at runtime):
```
hgnc_id	symbol	entrez_id	prev_symbol	alias_symbol
HGNC:11998	TP53	7157	P53	TRP53|LFS1|BCC7
HGNC:3373	EP300	2033		KAT3B|p300
HGNC:644	APC	324	DP2.5	GS|FAP|FPC
HGNC:11892	AR	367		DHTR|HUMARA|HYSP1|KD|NR3C4|SBMA|SMAX1|TFM
HGNC:6973	MYC	4609	c-Myc	bHLHe39|MRTL|MYCC
```

- [ ] **Step 2: Write the failing test**

`tests/unit/test_id_resolver.py`:
```python
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
```

- [ ] **Step 3: Run, expect failure**

Run: `python -m pytest tests/unit/test_id_resolver.py -v`
Expected: `ModuleNotFoundError`.

- [ ] **Step 4: Implement `id_resolver`**

`src/python/utils/id_resolver.py`:
```python
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
```

- [ ] **Step 5: Run, expect pass**

Run: `python -m pytest tests/unit/test_id_resolver.py -v`
Expected: 6 passed.

- [ ] **Step 6: Commit**

```bash
git add src/python/utils/id_resolver.py tests/unit/test_id_resolver.py tests/fixtures/hgnc_mini.tsv
git commit -m "feat(utils): add HGNC-based IDResolver with miss-rate tracking [task-0.4]"
```

---

### Task 0.5: Capture pre-audit golden artifacts

**Files:**
- Create: `scripts/capture_golden.sh`
- Create: `tests/regression/golden/.gitkeep`
- Output: `tests/regression/golden/dice_genes_pre_audit.tsv` etc.

- [ ] **Step 1: Write the capture script**

`scripts/capture_golden.sh`:
```bash
#!/usr/bin/env bash
# Snapshot current pipeline outputs as golden references for the
# --dice-faithful regression contract. Run BEFORE any code change.
set -euo pipefail
cd "$(dirname "$0")/.."

DEST=tests/regression/golden
mkdir -p "$DEST"

for f in \
    data/results/dice_genes.txt \
    data/results/dice_genes_top500.txt \
    data/results/centrality_normal.txt \
    data/results/centrality_tumor.txt \
    data/results/knockout_analysis.txt \
    data/results/dnf_pairwise_results.txt; do
    if [[ -f "$f" ]]; then
        cp "$f" "$DEST/$(basename ${f%.txt})_pre_audit.tsv"
        echo "captured $f"
    else
        echo "missing $f (skipped)"
    fi
done

# Also record the commit & timestamp so the regression test can reference it
git rev-parse HEAD > "$DEST/_captured_at_commit.txt"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DEST/_captured_at_time.txt"
```

Then `chmod +x scripts/capture_golden.sh`.

- [ ] **Step 2: Run the capture**

```bash
./scripts/capture_golden.sh
ls -la tests/regression/golden/
```

Expected: 6 `.tsv` files plus `_captured_at_commit.txt` and `_captured_at_time.txt`.

- [ ] **Step 3: Verify file sizes**

```bash
wc -l tests/regression/golden/*.tsv
```

Expected: counts roughly matching `data/results/*.txt` line counts (e.g. `dice_genes_pre_audit.tsv` ≈ 8961 lines).

- [ ] **Step 4: Commit**

The golden files may be large but they are essential references. Add an `.gitattributes` line to mark them as binary-equivalent so diffs don't pollute reviews.

Append to `.gitattributes`:
```
tests/regression/golden/*.tsv -diff
```

Then:
```bash
git add scripts/capture_golden.sh tests/regression/golden/ .gitattributes
git commit -m "chore(audit): snapshot pre-audit golden artifacts for dice-faithful regression [task-0.5]"
```

---

## PHASE 1 — Diagnostic dossier (12 tasks)

### Task 1.0: Probe base class (TDD)

**Files:**
- Create: `src/python/diagnostics/__init__.py`
- Create: `src/python/diagnostics/probe_base.py`
- Create: `tests/unit/test_probe_base.py`

- [ ] **Step 1: Write the failing test**

`tests/unit/test_probe_base.py`:
```python
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
```

- [ ] **Step 2: Run, expect failure (ModuleNotFoundError)**

Run: `python -m pytest tests/unit/test_probe_base.py -v`

- [ ] **Step 3: Implement `probe_base`**

`src/python/diagnostics/__init__.py`:
```python
```
(empty)

`src/python/diagnostics/probe_base.py`:
```python
"""Base class for pipeline boundary probes. Pure observation; never mutates."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, Optional
import logging
import time
import traceback

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
```

- [ ] **Step 4: Run, expect pass**

Expected: 2 passed.

- [ ] **Step 5: Commit**

```bash
git add src/python/diagnostics/__init__.py src/python/diagnostics/probe_base.py tests/unit/test_probe_base.py
git commit -m "feat(diagnostics): add Probe base class + ProbeResult [task-1.0]"
```

---

### Task 1.1: Probe — raw inputs

**Files:**
- Create: `src/python/diagnostics/probes/__init__.py`
- Create: `src/python/diagnostics/probes/raw_inputs.py`
- Modify: `tests/unit/test_probes.py` (add tests)

- [ ] **Step 1: Write failing test**

`tests/unit/test_probes.py`:
```python
import pandas as pd
from src.python.diagnostics.probes.raw_inputs import RawInputsProbe


def test_raw_inputs_probe_counts_rows_and_detects_id_type(
    tmp_path, tiny_expr_normal, tiny_expr_tumor
):
    norm = tmp_path / "norm.tsv"
    tum = tmp_path / "tum.tsv"
    tiny_expr_normal.to_csv(norm, sep="\t")
    tiny_expr_tumor.to_csv(tum, sep="\t")

    r = RawInputsProbe().run(ctx={
        "expression_normal": str(norm),
        "expression_tumor": str(tum),
        "kegg_dir": None,           # accepted but skipped if None
        "biogrid_file": None,
    })
    assert r.error is None
    m = r.measurements
    assert m["expression_normal"]["n_genes"] == 10
    assert m["expression_normal"]["n_samples"] == 5
    assert m["expression_normal"]["nan_cells"] == 0
    assert m["expression_normal"]["id_type"] == "symbol"
    assert m["expression_tumor"]["n_genes"] == 10
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

`src/python/diagnostics/probes/__init__.py` — empty.

`src/python/diagnostics/probes/raw_inputs.py`:
```python
"""Probe stage 1: raw input files (KEGG, BioGRID, expression)."""
from __future__ import annotations
import re
from pathlib import Path
from typing import Any, Dict

import pandas as pd

from ..probe_base import Probe

ENTREZ_RE = re.compile(r"^\d+$")


def _detect_id_type(idx: pd.Index) -> str:
    sample = [str(x) for x in idx[:200]]
    if not sample:
        return "empty"
    entrez_hits = sum(1 for s in sample if ENTREZ_RE.match(s))
    if entrez_hits / len(sample) > 0.9:
        return "entrez"
    if entrez_hits / len(sample) > 0.1:
        return "mixed"
    return "symbol"


def _expr_stats(path: str) -> Dict[str, Any]:
    df = pd.read_csv(path, sep="\t", index_col=0)
    return {
        "path": path,
        "n_genes": len(df),
        "n_samples": df.shape[1],
        "nan_cells": int(df.isna().sum().sum()),
        "duplicate_genes": int(df.index.duplicated().sum()),
        "id_type": _detect_id_type(df.index),
    }


class RawInputsProbe(Probe):
    name = "raw_inputs"
    stage = "1-raw-inputs"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        m: Dict[str, Any] = {}
        if ctx.get("expression_normal"):
            m["expression_normal"] = _expr_stats(ctx["expression_normal"])
        if ctx.get("expression_tumor"):
            m["expression_tumor"] = _expr_stats(ctx["expression_tumor"])

        if ctx.get("kegg_dir"):
            kegg = Path(ctx["kegg_dir"])
            xmls = list(kegg.glob("*.xml")) if kegg.exists() else []
            m["kegg"] = {"path": str(kegg), "n_xml_files": len(xmls)}

        if ctx.get("biogrid_file"):
            bg = Path(ctx["biogrid_file"])
            size = bg.stat().st_size if bg.exists() else 0
            m["biogrid"] = {"path": str(bg), "exists": bg.exists(), "size_bytes": size}

        return m
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/diagnostics/probes/__init__.py src/python/diagnostics/probes/raw_inputs.py tests/unit/test_probes.py
git commit -m "feat(diagnostics): probe-1 raw inputs (counts, id-type, NaN cells) [task-1.1]"
```

---

### Task 1.2: Probe — parsers (KEGG, BioGRID outputs)

**Files:**
- Create: `src/python/diagnostics/probes/parsers.py`
- Modify: `tests/unit/test_probes.py` (extend)

- [ ] **Step 1: Write failing test**

Append to `tests/unit/test_probes.py`:
```python
from src.python.diagnostics.probes.parsers import ParsersProbe


def test_parsers_probe_counts_edges_and_unique_genes(tmp_path):
    kegg = tmp_path / "kegg_interactions.txt"
    kegg.write_text("source\ttarget\tdirected\nA\tB\ttrue\nB\tC\ttrue\n")
    bgr = tmp_path / "biogrid_interactions.txt"
    bgr.write_text("source\ttarget\nA\tD\nC\tD\nA\tD\n")  # one duplicate

    r = ParsersProbe().run(ctx={
        "kegg_interactions": str(kegg),
        "biogrid_interactions": str(bgr),
    })
    m = r.measurements
    assert m["kegg"]["n_edges"] == 2
    assert m["kegg"]["unique_genes"] == 3
    assert m["biogrid"]["n_edges"] == 3
    assert m["biogrid"]["duplicate_edges"] == 1
    assert m["biogrid"]["unique_genes"] == 3
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

`src/python/diagnostics/probes/parsers.py`:
```python
"""Probe stage 2: parsed KEGG / BioGRID outputs."""
from __future__ import annotations
from typing import Any, Dict
import pandas as pd

from ..probe_base import Probe


def _edge_stats(path: str) -> Dict[str, Any]:
    df = pd.read_csv(path, sep="\t")
    src_col = next((c for c in df.columns if c.lower() in ("source", "src")), df.columns[0])
    tgt_col = next((c for c in df.columns if c.lower() in ("target", "tgt", "dst")), df.columns[1])
    pairs = df[[src_col, tgt_col]].astype(str)
    duplicates = pairs.duplicated().sum()
    unique_genes = len(set(pairs[src_col]).union(pairs[tgt_col]))
    self_loops = int((pairs[src_col] == pairs[tgt_col]).sum())
    return {
        "path": path,
        "n_edges": len(pairs),
        "duplicate_edges": int(duplicates),
        "unique_genes": int(unique_genes),
        "self_loops": self_loops,
    }


class ParsersProbe(Probe):
    name = "parsers"
    stage = "2-parsers"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        m: Dict[str, Any] = {}
        if ctx.get("kegg_interactions"):
            m["kegg"] = _edge_stats(ctx["kegg_interactions"])
        if ctx.get("biogrid_interactions"):
            m["biogrid"] = _edge_stats(ctx["biogrid_interactions"])
        return m
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/diagnostics/probes/parsers.py tests/unit/test_probes.py
git commit -m "feat(diagnostics): probe-2 parsers (KEGG/BioGRID edges, dupes) [task-1.2]"
```

---

### Task 1.3: Probe — merge_networks (directionality breakdown)

**Files:**
- Create: `src/python/diagnostics/probes/merge.py`
- Modify: `tests/unit/test_probes.py`

- [ ] **Step 1: Test**

Append:
```python
from src.python.diagnostics.probes.merge import MergeProbe


def test_merge_probe_counts_directed_undirected(tmp_path):
    merged = tmp_path / "merged.txt"
    merged.write_text(
        "source\ttarget\tdirected\n"
        "A\tB\ttrue\n"
        "B\tA\ttrue\n"        # forms a 2-cycle but each entry is directed
        "C\tD\tfalse\n"
        "D\tE\tfalse\n"
    )
    r = MergeProbe().run(ctx={"merged_network": str(merged)})
    m = r.measurements
    assert m["n_edges"] == 4
    assert m["n_directed"] == 2
    assert m["n_undirected"] == 2
    assert m["unique_genes"] == 5
    assert m["bidirectional_pairs"] == 1   # A-B has both directions
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

`src/python/diagnostics/probes/merge.py`:
```python
"""Probe stage 3: merged network — directionality + structural stats."""
from __future__ import annotations
from typing import Any, Dict
import pandas as pd

from ..probe_base import Probe


class MergeProbe(Probe):
    name = "merge_networks"
    stage = "3-merge"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        path = ctx["merged_network"]
        df = pd.read_csv(path, sep="\t")
        n_edges = len(df)
        unique_genes = len(set(df["source"]).union(df["target"]))
        if "directed" in df.columns:
            n_dir = int(df["directed"].astype(str).str.lower().eq("true").sum())
            n_und = n_edges - n_dir
        else:
            n_dir, n_und = 0, n_edges

        # bidirectional pairs (X->Y and Y->X both present in directed subset)
        if "directed" in df.columns:
            d = df[df["directed"].astype(str).str.lower() == "true"]
        else:
            d = df.iloc[0:0]
        pairs = set(zip(d["source"], d["target"]))
        bid = sum(1 for (s, t) in pairs if (t, s) in pairs and s < t)

        return {
            "path": path,
            "n_edges": n_edges,
            "n_directed": n_dir,
            "n_undirected": n_und,
            "unique_genes": unique_genes,
            "bidirectional_pairs": bid,
            "self_loops": int((df["source"] == df["target"]).sum()),
        }
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/diagnostics/probes/merge.py tests/unit/test_probes.py
git commit -m "feat(diagnostics): probe-3 merge (directed/undirected/bidirectional) [task-1.3]"
```

---

### Task 1.4: Probe — DEA / IG / intersection

**Files:**
- Create: `src/python/diagnostics/probes/filters.py`
- Modify: `tests/unit/test_probes.py`

- [ ] **Step 1: Test**

Append:
```python
from src.python.diagnostics.probes.filters import FiltersProbe


def test_filters_probe_reports_dea_ig_intersection(tmp_path):
    dea_table = tmp_path / "dea.tsv"
    dea_table.write_text(
        "gene\tlog2FC\tpvalue\tFDR\n"
        "G1\t2.0\t0.001\t0.01\n"
        "G2\t-2.0\t0.001\t0.01\n"
        "G3\t0.5\t0.5\t0.5\n"
        "G4\t3.0\t0.001\t0.01\n"
    )
    dea_genes = tmp_path / "dea_genes.txt"
    dea_genes.write_text("G1\nG2\nG4\n")
    ig_genes = tmp_path / "ig_genes.txt"
    ig_genes.write_text("G2\nG4\nG5\n")
    inter = tmp_path / "candidates.txt"
    inter.write_text("G2\nG4\n")

    r = FiltersProbe().run(ctx={
        "dea_table": str(dea_table),
        "dea_genes": str(dea_genes),
        "ig_genes": str(ig_genes),
        "candidates": str(inter),
    })
    m = r.measurements
    assert m["dea"]["n_genes_total"] == 4
    assert m["dea"]["n_significant"] == 3
    assert m["dea"]["fdr_p10"] <= m["dea"]["fdr_p90"]
    assert m["ig"]["n_genes"] == 3
    assert m["intersection"]["n_genes"] == 2
    assert m["intersection"]["jaccard_dea_ig"] == pytest.approx(0.5)
```

(Add `import pytest` to the top of the file if missing.)

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

`src/python/diagnostics/probes/filters.py`:
```python
"""Probe stage 4: DEA + IG + intersection."""
from __future__ import annotations
from typing import Any, Dict
import numpy as np
import pandas as pd

from ..probe_base import Probe


def _read_gene_list(path: str) -> set:
    with open(path) as f:
        return {ln.strip() for ln in f if ln.strip()}


class FiltersProbe(Probe):
    name = "filters"
    stage = "4-filters"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        m: Dict[str, Any] = {}
        if ctx.get("dea_table"):
            df = pd.read_csv(ctx["dea_table"], sep="\t")
            sig = df[(df["FDR"] < 0.05) & (df["log2FC"].abs() > 1.0)]
            m["dea"] = {
                "n_genes_total": len(df),
                "n_significant": len(sig),
                "fdr_p10": float(np.quantile(df["FDR"], 0.1)),
                "fdr_p50": float(np.quantile(df["FDR"], 0.5)),
                "fdr_p90": float(np.quantile(df["FDR"], 0.9)),
                "fdr_le_0_05_pct": float((df["FDR"] < 0.05).mean()),
                "abs_lfc_gt_1_pct": float((df["log2FC"].abs() > 1.0).mean()),
            }
        dea_set = _read_gene_list(ctx["dea_genes"]) if ctx.get("dea_genes") else set()
        ig_set = _read_gene_list(ctx["ig_genes"]) if ctx.get("ig_genes") else set()
        inter_set = _read_gene_list(ctx["candidates"]) if ctx.get("candidates") else set()
        if dea_set:
            m.setdefault("dea", {})["n_genes_listed"] = len(dea_set)
        if ig_set:
            m["ig"] = {"n_genes": len(ig_set)}
        if dea_set and ig_set:
            jac = len(dea_set & ig_set) / max(1, len(dea_set | ig_set))
            m["intersection"] = {
                "n_genes": len(inter_set),
                "n_dea_only": len(dea_set - ig_set),
                "n_ig_only": len(ig_set - dea_set),
                "jaccard_dea_ig": float(jac),
            }
        return m
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/diagnostics/probes/filters.py tests/unit/test_probes.py
git commit -m "feat(diagnostics): probe-4 filters (DEA/IG/intersection stats) [task-1.4]"
```

---

### Task 1.5: Probe — weight_network (the 73 % drop probe)

This is the largest probe. It must explain *every* edge that disappears between `merged_network.txt` and the weighted output. Specifically: NaN-correlation drops, ε-clamp hits, integer-truncation collisions, gene-symbol resolution rate against the expression matrix.

**Files:**
- Create: `src/python/diagnostics/probes/weighting.py`
- Modify: `tests/unit/test_probes.py`

- [ ] **Step 1: Test**

Append:
```python
from src.python.diagnostics.probes.weighting import WeightingProbe


def test_weighting_probe_explains_edge_loss(tmp_path):
    merged = tmp_path / "merged.txt"
    merged.write_text(
        "source\ttarget\n"
        "A\tB\n"
        "A\tZZUNKNOWN\n"   # dropped: target not in expression matrix
        "C\tD\n"
        "E\tF\n"
    )
    weighted = tmp_path / "weighted.txt"
    weighted.write_text(
        "source\ttarget\tweight\tweight_int\n"
        "A\tB\t0.05\t50\n"
        "C\tD\t0.000001\t1\n"     # epsilon-clamp survivor
        "E\tF\t0.5\t500\n"
    )
    expr = tmp_path / "expr.tsv"
    expr.write_text("gene\tS1\tS2\nA\t1\t2\nB\t2\t1\nC\t1\t1\nD\t1\t1\nE\t1\t2\nF\t2\t1\n")

    r = WeightingProbe().run(ctx={
        "merged_network": str(merged),
        "weighted_network": str(weighted),
        "expression": str(expr),
        "epsilon": 1e-6,
        "scale_factor": 1000,
    })
    m = r.measurements
    assert m["n_input_edges"] == 4
    assert m["n_output_edges"] == 3
    assert m["n_dropped"] == 1
    assert m["edges_unresolved_in_expression"] == 1
    assert m["epsilon_clamped_edges"] >= 1
    assert m["gene_symbol_resolution_rate"] == pytest.approx(5/6)  # 5 of 6 genes in expr
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

`src/python/diagnostics/probes/weighting.py`:
```python
"""Probe stage 5: weight_network — explain every edge that disappears."""
from __future__ import annotations
from typing import Any, Dict, Optional
import numpy as np
import pandas as pd

from ..probe_base import Probe


class WeightingProbe(Probe):
    name = "weighting"
    stage = "5-weighting"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        merged = pd.read_csv(ctx["merged_network"], sep="\t")
        weighted = pd.read_csv(ctx["weighted_network"], sep="\t")
        expr = pd.read_csv(ctx["expression"], sep="\t", index_col=0)
        eps = float(ctx.get("epsilon", 1e-6))
        scale = int(ctx.get("scale_factor", 1000))

        merged_pairs = set(zip(merged["source"], merged["target"]))
        weighted_pairs = set(zip(weighted["source"], weighted["target"]))
        dropped_pairs = merged_pairs - weighted_pairs

        # Decompose drops by reason
        expr_genes = set(expr.index)
        unresolved = sum(1 for s, t in dropped_pairs if s not in expr_genes or t not in expr_genes)

        # Resolution rate: union of all genes in merged net vs expression
        all_merged_genes = set(merged["source"]).union(merged["target"])
        resolution_rate = len(all_merged_genes & expr_genes) / max(1, len(all_merged_genes))

        # Epsilon clamp hits: edges with weight near epsilon (within 1.5x)
        eps_clamp = int((weighted["weight"] <= eps * 1.5).sum())

        # Integer-truncation collisions
        int_w = (weighted["weight"] * scale).astype(int)
        unique_weights = weighted["weight"].nunique()
        unique_int = int_w.nunique()
        truncation_collisions = int(unique_weights - unique_int)

        # Bottom-clamp from weight_int==0 -> 1 in current code
        bottom_clamp_to_1 = int((weighted["weight_int"] == 1).sum())

        return {
            "n_input_edges": len(merged),
            "n_output_edges": len(weighted),
            "n_dropped": len(dropped_pairs),
            "edges_unresolved_in_expression": unresolved,
            "gene_symbol_resolution_rate": float(resolution_rate),
            "epsilon_clamped_edges": eps_clamp,
            "weight_truncation_collisions": truncation_collisions,
            "weight_int_eq_1_count": bottom_clamp_to_1,
            "weight_distribution": {
                "min": float(weighted["weight"].min()),
                "max": float(weighted["weight"].max()),
                "p10": float(np.quantile(weighted["weight"], 0.1)),
                "p50": float(np.quantile(weighted["weight"], 0.5)),
                "p90": float(np.quantile(weighted["weight"], 0.9)),
            },
        }
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/diagnostics/probes/weighting.py tests/unit/test_probes.py
git commit -m "feat(diagnostics): probe-5 weighting (73pct drop decomposition) [task-1.5]"
```

---

### Task 1.6: C++ `--diagnostics` flag

Adds three diagnostic outputs to `dice_analyzer`: (a) sampled vs full Brandes Spearman ρ on a 1k-vertex sub-induced graph, (b) OMP-determinism check (5×, max relative drift), (c) directed-graph Brandes truth on a hand-built toy embedded in code. Output goes to JSON file specified by `--diagnostics-out`.

**Files:**
- Modify: `src/cpp/main.cpp`
- Modify: `src/cpp/engine.hpp` (add `betweenness(G, sample_rate, seed)` overload if it doesn't exist; do not change existing overloads)
- Add: `src/cpp/tests/test_engine.cpp` (extend with two new test cases)

- [ ] **Step 1: Add C++ tests first**

Append to `src/cpp/tests/test_engine.cpp`:
```cpp
// Same 3-node path graph; sampled with rate 1.0 must match full.
TEST_CASE("Sampled Brandes equals full Brandes when sample_rate=1.0") {
    Graph<long long> g(3);
    g.id_to_name = {"A", "B", "C"};
    g.name_to_id = {{"A", 0}, {"B", 1}, {"C", 2}};
    g.add_edge(0, 1, 1);
    g.add_edge(1, 2, 1);

    CentralityCalculator<long long> calc;
    auto bet_full = calc.betweenness(g, 1.0);
    auto bet_sampled = calc.betweenness(g, 1.0, /*seed=*/42);
    for (size_t i = 0; i < bet_full.size(); ++i)
        CHECK(bet_full[i] == doctest::Approx(bet_sampled[i]));
}

// Star graph 1 -> {2,3,4,5}; node 1 has betweenness 0 (sink only).
TEST_CASE("Brandes on directed star with no transit yields zero centrality") {
    Graph<long long> g(5);
    g.id_to_name = {"H", "L1", "L2", "L3", "L4"};
    g.add_edge(0, 1, 1);
    g.add_edge(0, 2, 1);
    g.add_edge(0, 3, 1);
    g.add_edge(0, 4, 1);
    CentralityCalculator<long long> calc;
    auto bet = calc.betweenness(g, 1.0);
    for (auto v : bet) CHECK(v == doctest::Approx(0.0));
}
```

- [ ] **Step 2: Build & run C++ tests; expect failure**

```bash
cmake --build build --target test_engine && ./build/test_engine
```
Expected: at least one of these fails to compile or fails the assertion if `betweenness(G, sample_rate, seed)` doesn't exist.

- [ ] **Step 3: Add the seed overload to `engine.hpp`**

Locate the existing `betweenness(...)` method on `CentralityCalculator`. Add this overload (do **not** change the existing one):
```cpp
// In engine.hpp, inside CentralityCalculator<WeightType>:
std::vector<double> betweenness(Graph<WeightType>& g, double sample_rate, uint64_t seed) {
    // Wrapper for the existing implementation; uses a deterministic generator
    // when sampling sources, instead of std::rand or unseeded random_device.
    if (sample_rate >= 1.0) return betweenness(g, 1.0);
    std::mt19937_64 rng(seed);
    std::vector<int> all_nodes(g.n_vertices);
    std::iota(all_nodes.begin(), all_nodes.end(), 0);
    int k = std::max<int>(1, int(g.n_vertices * sample_rate));
    std::shuffle(all_nodes.begin(), all_nodes.end(), rng);
    all_nodes.resize(k);
    return betweenness_from_sources(g, all_nodes);  // existing private helper or refactor target
}
```
If `betweenness_from_sources` does not exist, refactor the existing `betweenness(g, sample_rate)` body so the per-source loop runs over a vector of source indices passed in. Keep the original public method calling the new helper with a deterministic shuffle. **All existing callers must continue to behave bit-identically when seed is unspecified** — verify by re-running the existing C++ binary against the captured golden artifacts after this refactor.

- [ ] **Step 4: Build & run; expect tests pass**

```bash
cmake --build build --target test_engine && ./build/test_engine
```

- [ ] **Step 5: Add `--diagnostics-out` flag to `main.cpp`**

In `main.cpp`, after the existing CLI parsing, add:
```cpp
// CLI extension — run diagnostics mode after main analysis if requested.
std::string diag_out;
for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--diagnostics-out" && i + 1 < argc) { diag_out = argv[++i]; }
}

if (!diag_out.empty()) {
    // (a) OMP determinism: run betweenness 5x with same seed; max relative drift.
    auto run_seeded = [&](Graph<long long>& g, uint64_t seed) {
        CentralityCalculator<long long> c;
        return c.betweenness(g, 0.1, seed);   // 10% sample, deterministic
    };
    Graph<long long> g_t = GraphLoader<long long>::load_from_file(argv[2]); // tumor
    auto base = run_seeded(g_t, 42);
    double max_drift = 0.0;
    for (int t = 0; t < 5; ++t) {
        auto r = run_seeded(g_t, 42);
        for (size_t i = 0; i < base.size(); ++i) {
            double denom = std::max(1e-12, std::abs(base[i]));
            max_drift = std::max(max_drift, std::abs(r[i] - base[i]) / denom);
        }
    }

    // (b) Sampled vs full Spearman on a 1k-induced subgraph.
    // Pick first min(1000, n) nodes; emit ranked vectors; Python-side
    // computes Spearman from the JSON.
    int k = std::min(1000, g_t.n_vertices);
    Graph<long long> sub = induce_subgraph(g_t, k);  // small helper; see note below.
    CentralityCalculator<long long> c_sub;
    auto bet_full = c_sub.betweenness(sub, 1.0);
    auto bet_smp  = c_sub.betweenness(sub, 0.1, 42);

    // Write JSON
    std::ofstream out(diag_out);
    out << "{\n";
    out << "  \"omp_determinism_max_relative_drift\": " << max_drift << ",\n";
    out << "  \"sampled_vs_full_subgraph_n\": " << k << ",\n";
    out << "  \"sampled_full_betweenness\": [";
    for (size_t i = 0; i < bet_full.size(); ++i) out << (i?",":"") << bet_full[i];
    out << "],\n  \"sampled_betweenness\": [";
    for (size_t i = 0; i < bet_smp.size(); ++i) out << (i?",":"") << bet_smp[i];
    out << "]\n}\n";
}
```

If `induce_subgraph` does not exist, add a small helper in `engine.hpp`:
```cpp
template<typename W>
Graph<W> induce_subgraph(const Graph<W>& g, int k) {
    Graph<W> s(k);
    for (int u = 0; u < k; ++u)
        for (auto& e : g.adj[u])
            if (e.to < k) s.add_edge(u, e.to, e.weight);
    s.id_to_name.assign(g.id_to_name.begin(), g.id_to_name.begin() + k);
    return s;
}
```

- [ ] **Step 6: Manual smoke test**

```bash
cmake --build build --target dice_analyzer
./build/dice_analyzer \
    data/processed/weighted/network_normal_int.txt \
    data/processed/weighted/network_tumor_int.txt \
    data/results \
    --diagnostics-out data/results/cpp_diagnostics.json
cat data/results/cpp_diagnostics.json | head -5
```

Expected: file written; `omp_determinism_max_relative_drift` is some floating-point number (ideally 0 — if larger than ~1e-9, that is itself a bug-registry entry).

- [ ] **Step 7: Commit**

```bash
git add src/cpp/engine.hpp src/cpp/main.cpp src/cpp/tests/test_engine.cpp
git commit -m "feat(engine): seeded Brandes overload + --diagnostics-out emits OMP/sub-Spearman JSON [task-1.6]"
```

---

### Task 1.7: Probe — centrality (consume C++ diagnostics)

**Files:**
- Create: `src/python/diagnostics/probes/centrality.py`
- Modify: `tests/unit/test_probes.py`

- [ ] **Step 1: Test**

Append:
```python
import json
from src.python.diagnostics.probes.centrality import CentralityProbe


def test_centrality_probe_computes_spearman(tmp_path):
    diag = tmp_path / "cpp_diagnostics.json"
    diag.write_text(json.dumps({
        "omp_determinism_max_relative_drift": 0.0,
        "sampled_vs_full_subgraph_n": 5,
        "sampled_full_betweenness": [0.0, 1.0, 2.0, 3.0, 4.0],
        "sampled_betweenness":      [0.1, 1.1, 2.1, 3.1, 4.1],
    }))
    r = CentralityProbe().run(ctx={"cpp_diagnostics_json": str(diag)})
    m = r.measurements
    assert m["omp_determinism_max_relative_drift"] == 0.0
    assert m["sampled_vs_full_spearman_rho"] == pytest.approx(1.0)
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

`src/python/diagnostics/probes/centrality.py`:
```python
"""Probe stage 6: consume C++ --diagnostics-out JSON; compute Spearman."""
from __future__ import annotations
import json
from typing import Any, Dict
from scipy import stats

from ..probe_base import Probe


class CentralityProbe(Probe):
    name = "centrality"
    stage = "6-centrality"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        with open(ctx["cpp_diagnostics_json"]) as f:
            d = json.load(f)
        full = d["sampled_full_betweenness"]
        smp = d["sampled_betweenness"]
        rho, p = stats.spearmanr(full, smp)
        return {
            "omp_determinism_max_relative_drift": d["omp_determinism_max_relative_drift"],
            "sampled_vs_full_subgraph_n": d["sampled_vs_full_subgraph_n"],
            "sampled_vs_full_spearman_rho": float(rho),
            "sampled_vs_full_spearman_p": float(p),
        }
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/diagnostics/probes/centrality.py tests/unit/test_probes.py
git commit -m "feat(diagnostics): probe-6 centrality (Spearman + OMP determinism) [task-1.7]"
```

---

### Task 1.8: Probe — differential_centrality

Reports rank correlation between betweenness and eigenvector ranks (rank-product double-count test); enumerates the top-12 % eigenvector genes that the noise filter would kill; intersects them with COSMIC Cancer Gene Census if available.

**Files:**
- Create: `src/python/diagnostics/probes/differential.py`
- Create: `tests/fixtures/cgc_mini.txt`
- Modify: `tests/unit/test_probes.py`

- [ ] **Step 1: Mini CGC fixture**

`tests/fixtures/cgc_mini.txt`:
```
TP53
EP300
APC
MYC
AR
KRAS
PIK3CA
```

- [ ] **Step 2: Test**

```python
from src.python.diagnostics.probes.differential import DifferentialProbe

CGC = Path(__file__).parents[1] / "fixtures" / "cgc_mini.txt"


def test_differential_probe_finds_hub_kills(tmp_path):
    centrality_normal = tmp_path / "cn.txt"
    centrality_normal.write_text(
        "gene\tbetweenness\teigenvector\n"
        + "\n".join(f"G{i}\t{i}\t{0.01*i}" for i in range(100)) + "\n"
        + "TP53\t10\t0.99\nEP300\t12\t0.95\n"
    )
    centrality_tumor = tmp_path / "ct.txt"
    centrality_tumor.write_text(
        "gene\tbetweenness\teigenvector\n"
        + "\n".join(f"G{i}\t{i+1}\t{0.01*i}" for i in range(100)) + "\n"
        + "TP53\t9\t0.97\nEP300\t11\t0.93\n"
    )
    r = DifferentialProbe().run(ctx={
        "centrality_normal": str(centrality_normal),
        "centrality_tumor": str(centrality_tumor),
        "cgc_list": str(CGC),
    })
    m = r.measurements
    # bet/eig rank correlation defined and bounded
    assert -1.0 <= m["bet_eig_rank_correlation"] <= 1.0
    # TP53 + EP300 are in the top-12% of eigenvector and in CGC
    assert "TP53" in m["top12_eig_dropped_genes"]
    assert "EP300" in m["top12_eig_dropped_genes"]
    assert m["top12_eig_dropped_in_cgc"] >= 2
```

- [ ] **Step 3: Run, expect failure**

- [ ] **Step 4: Implement**

`src/python/diagnostics/probes/differential.py`:
```python
"""Probe stage 7: differential_centrality — rank-product double-count + hub kill list."""
from __future__ import annotations
from typing import Any, Dict
import pandas as pd
from scipy import stats

from ..probe_base import Probe


def _read_set(path):
    with open(path) as f:
        return {ln.strip() for ln in f if ln.strip()}


class DifferentialProbe(Probe):
    name = "differential"
    stage = "7-differential"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        cn = pd.read_csv(ctx["centrality_normal"], sep="\t")
        ct = pd.read_csv(ctx["centrality_tumor"], sep="\t")
        merged = cn.merge(ct, on="gene", suffixes=("_n", "_t"))

        diff_bet = (merged["betweenness_t"] - merged["betweenness_n"]).abs()
        diff_eig = (merged["eigenvector_t"] - merged["eigenvector_n"]).abs()
        rho, _ = stats.spearmanr(diff_bet, diff_eig)

        # Top-12% eigenvector hub list (would-be-dropped under DiCE noise filter)
        eig_thresh_n = merged["eigenvector_n"].quantile(0.88)
        eig_thresh_t = merged["eigenvector_t"].quantile(0.88)
        kill_mask = (merged["eigenvector_n"] > eig_thresh_n) | (merged["eigenvector_t"] > eig_thresh_t)
        kill_genes = merged.loc[kill_mask, "gene"].tolist()

        cgc = _read_set(ctx["cgc_list"]) if ctx.get("cgc_list") else set()
        in_cgc = sum(1 for g in kill_genes if g in cgc)

        return {
            "n_genes": len(merged),
            "bet_eig_rank_correlation": float(rho),
            "top12_eig_dropped_n": len(kill_genes),
            "top12_eig_dropped_genes": kill_genes[:50],   # cap for dossier
            "top12_eig_dropped_in_cgc": in_cgc,
            "cgc_size": len(cgc),
        }
```

- [ ] **Step 5: Run, expect pass**

- [ ] **Step 6: Commit**

```bash
git add src/python/diagnostics/probes/differential.py tests/fixtures/cgc_mini.txt tests/unit/test_probes.py
git commit -m "feat(diagnostics): probe-7 differential (rank corr + top-12% hub kill list) [task-1.8]"
```

---

### Task 1.9: Probe — validation pipelines

Reports DEG-overlap reproduction (the "0 overlap" claim), DepMap join match-rate, survival cohort size at each filter step.

**Files:**
- Create: `src/python/diagnostics/probes/validation.py`
- Modify: `tests/unit/test_probes.py`

- [ ] **Step 1: Test**

```python
from src.python.diagnostics.probes.validation import ValidationProbe


def test_validation_probe_reports_overlap_and_join_rate(tmp_path):
    dice = tmp_path / "dice.tsv"
    dice.write_text("gene\nG1\nG2\nG3\nG4\nG5\n")
    deg  = tmp_path / "deg.tsv"
    deg.write_text("gene\nG3\nG4\nG6\n")
    depmap = tmp_path / "depmap.tsv"
    depmap.write_text("gene\tdependency_score\nG1\t-0.5\nG3\t-0.7\nG7\t0.0\n")

    r = ValidationProbe().run(ctx={
        "dice_genes": str(dice),
        "deg_genes": str(deg),
        "depmap_table": str(depmap),
    })
    m = r.measurements
    assert m["dice_deg_overlap"]["observed"] == 2
    assert 0.0 < m["dice_deg_overlap"]["hypergeom_p"] <= 1.0
    assert m["depmap"]["join_match_count"] == 2
    assert m["depmap"]["join_match_rate"] == pytest.approx(0.4)
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

`src/python/diagnostics/probes/validation.py`:
```python
"""Probe stage 8: validation pipelines — DEG overlap, DepMap join, survival cohort."""
from __future__ import annotations
from typing import Any, Dict
import pandas as pd
from scipy import stats

from ..probe_base import Probe


def _read_set(path):
    with open(path) as f:
        return {ln.strip() for ln in f if ln.strip() and ln.strip() != "gene"}


class ValidationProbe(Probe):
    name = "validation"
    stage = "8-validation"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        m: Dict[str, Any] = {}
        if ctx.get("dice_genes") and ctx.get("deg_genes") and ctx.get("universe_size"):
            dice = _read_set(ctx["dice_genes"])
            deg  = _read_set(ctx["deg_genes"])
            obs = len(dice & deg)
            N = int(ctx["universe_size"])
            # Hypergeom: prob >= obs given draws=|dice|, hits=|deg|, total=N
            rv = stats.hypergeom(N, len(deg), len(dice))
            p = float(rv.sf(obs - 1))
            m["dice_deg_overlap"] = {
                "observed": int(obs),
                "expected": float(len(dice) * len(deg) / max(1, N)),
                "hypergeom_p": p,
            }
        elif ctx.get("dice_genes") and ctx.get("deg_genes"):
            dice = _read_set(ctx["dice_genes"])
            deg  = _read_set(ctx["deg_genes"])
            obs = len(dice & deg)
            N = max(1, len(dice | deg)) * 10  # crude default universe
            rv = stats.hypergeom(N, len(deg), len(dice))
            m["dice_deg_overlap"] = {
                "observed": int(obs),
                "expected": float(len(dice) * len(deg) / N),
                "hypergeom_p": float(rv.sf(obs - 1)),
                "universe_assumed": N,
            }

        if ctx.get("depmap_table"):
            dm = pd.read_csv(ctx["depmap_table"], sep="\t")
            dice = _read_set(ctx["dice_genes"])
            matched = dm[dm["gene"].isin(dice)]
            m["depmap"] = {
                "depmap_genes": int(len(dm)),
                "join_match_count": int(len(matched)),
                "join_match_rate": float(len(matched) / max(1, len(dice))),
            }

        if ctx.get("survival_cohorts"):
            m["survival"] = {
                "cohort_sizes": ctx["survival_cohorts"],   # passed in by orchestrator
            }

        return m
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/diagnostics/probes/validation.py tests/unit/test_probes.py
git commit -m "feat(diagnostics): probe-8 validation (DEG overlap hypergeom + DepMap join) [task-1.9]"
```

---

### Task 1.10: Probe — DNF

**Files:**
- Create: `src/python/diagnostics/probes/dnf.py`
- Modify: `tests/unit/test_probes.py`

- [ ] **Step 1: Test**

```python
from src.python.diagnostics.probes.dnf import DNFProbe


def test_dnf_probe_summarizes_pairs(tmp_path):
    dnf = tmp_path / "dnf.tsv"
    dnf.write_text(
        "gene_a\tgene_b\tDSS\n"
        "G1\tG2\t0.1\n"
        "G1\tG3\t0.2\n"
        "G2\tG3\t-0.05\n"
    )
    sl = tmp_path / "sl.tsv"
    sl.write_text("gene_a\tgene_b\nG1\tG2\n")
    r = DNFProbe().run(ctx={
        "dnf_pairwise": str(dnf),
        "synlethdb": str(sl),
    })
    m = r.measurements
    assert m["n_pairs"] == 3
    assert m["dss_min"] == pytest.approx(-0.05)
    assert m["dss_max"] == pytest.approx(0.2)
    assert m["synlethdb_overlap"] == 1
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

`src/python/diagnostics/probes/dnf.py`:
```python
"""Probe stage 10: DNF — pair count, DSS distribution, SynLethDB overlap."""
from __future__ import annotations
from typing import Any, Dict
import numpy as np
import pandas as pd

from ..probe_base import Probe


class DNFProbe(Probe):
    name = "dnf"
    stage = "10-dnf"

    def measure(self, ctx: Dict[str, Any]) -> Dict[str, Any]:
        m: Dict[str, Any] = {}
        if ctx.get("dnf_pairwise"):
            df = pd.read_csv(ctx["dnf_pairwise"], sep="\t")
            score_col = next((c for c in df.columns if c.upper() == "DSS"), None)
            m["n_pairs"] = int(len(df))
            if score_col is not None and len(df):
                m["dss_min"] = float(df[score_col].min())
                m["dss_max"] = float(df[score_col].max())
                m["dss_p50"] = float(np.quantile(df[score_col], 0.5))
                m["dss_positive_pct"] = float((df[score_col] > 0).mean())
        if ctx.get("synlethdb") and ctx.get("dnf_pairwise"):
            sl = pd.read_csv(ctx["synlethdb"], sep="\t")
            sl_pairs = set(map(tuple, sl[["gene_a", "gene_b"]].astype(str).values))
            sl_pairs |= {(b, a) for (a, b) in sl_pairs}  # undirected
            df = pd.read_csv(ctx["dnf_pairwise"], sep="\t")
            dnf_pairs = set(map(tuple, df[["gene_a", "gene_b"]].astype(str).values))
            m["synlethdb_pairs"] = len(sl_pairs) // 2
            m["synlethdb_overlap"] = sum(1 for p in dnf_pairs if p in sl_pairs)
        return m
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/diagnostics/probes/dnf.py tests/unit/test_probes.py
git commit -m "feat(diagnostics): probe-10 dnf (DSS dist + SynLethDB overlap) [task-1.10]"
```

---

### Task 1.11: Dossier renderer (markdown + JSON)

**Files:**
- Create: `src/python/diagnostics/render.py`
- Create: `tests/unit/test_render.py`

- [ ] **Step 1: Test**

`tests/unit/test_render.py`:
```python
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
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

`src/python/diagnostics/render.py`:
```python
"""Render probe results as markdown + JSON."""
from __future__ import annotations
from datetime import datetime, timezone
import json
from typing import List

from .probe_base import ProbeResult


def _fmt_value(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    if isinstance(v, list):
        return f"[{', '.join(_fmt_value(x) for x in v[:10])}{'...' if len(v) > 10 else ''}]"
    return str(v)


def _fmt_dict(d, indent=0):
    lines = []
    pre = "  " * indent
    for k, v in d.items():
        if isinstance(v, dict):
            lines.append(f"{pre}- **{k}**:")
            lines.extend(_fmt_dict(v, indent + 1))
        else:
            lines.append(f"{pre}- **{k}**: {_fmt_value(v)}")
    return lines


def render_markdown(results: List[ProbeResult], commit: str) -> str:
    out = [
        "# Diagnostic Dossier",
        "",
        f"- **Commit**: `{commit}`",
        f"- **Generated**: {datetime.now(timezone.utc).isoformat()}",
        f"- **N probes**: {len(results)}",
        "",
        "---",
    ]
    by_stage = sorted(results, key=lambda r: r.stage)
    for r in by_stage:
        out.append(f"\n## Stage {r.stage} — `{r.name}`\n")
        if r.error:
            out.append(f"**ERROR**: `{r.error}`\n")
            continue
        out.append(f"- _duration_: {r.duration_s:.3f}s")
        out.extend(_fmt_dict(r.measurements))
    out.append("\n---\n")
    out.append("## Bug Registry\n")
    out.append("_Populate from per-stage findings; sort by impact × confidence × blast_radius._\n")
    out.append("| ID | Stage | Symptom | Evidence | Confidence | Blast radius | Hypothesized root cause |")
    out.append("|----|-------|---------|----------|------------|--------------|-------------------------|")
    return "\n".join(out) + "\n"


def render_json(results: List[ProbeResult], commit: str) -> str:
    return json.dumps({
        "commit": commit,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "probes": [
            {
                "name": r.name,
                "stage": r.stage,
                "duration_s": r.duration_s,
                "error": r.error,
                "measurements": r.measurements,
            }
            for r in results
        ],
    }, indent=2, default=str)
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/diagnostics/render.py tests/unit/test_render.py
git commit -m "feat(diagnostics): markdown + json renderer [task-1.11]"
```

---

### Task 1.12: Dossier runner — CLI orchestrator

**Files:**
- Create: `src/python/diagnostics/dossier_runner.py`
- Create: `tests/integration/test_dossier_runner.py`

- [ ] **Step 1: Test**

`tests/integration/test_dossier_runner.py`:
```python
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
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

`src/python/diagnostics/dossier_runner.py`:
```python
"""CLI: run all probes, render markdown + JSON. Pure observation."""
from __future__ import annotations
import argparse
import logging
import subprocess
from pathlib import Path
from typing import Optional

from .probe_base import Probe, ProbeResult
from .probes.raw_inputs import RawInputsProbe
from .probes.parsers import ParsersProbe
from .probes.merge import MergeProbe
from .probes.filters import FiltersProbe
from .probes.weighting import WeightingProbe
from .probes.centrality import CentralityProbe
from .probes.differential import DifferentialProbe
from .probes.validation import ValidationProbe
from .probes.dnf import DNFProbe
from .render import render_markdown, render_json

log = logging.getLogger(__name__)


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def main(argv: Optional[list] = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output", required=True)
    p.add_argument("--expression-normal")
    p.add_argument("--expression-tumor")
    p.add_argument("--kegg-dir")
    p.add_argument("--biogrid-file")
    p.add_argument("--kegg-interactions")
    p.add_argument("--biogrid-interactions")
    p.add_argument("--merged-network")
    p.add_argument("--dea-table"); p.add_argument("--dea-genes")
    p.add_argument("--ig-genes"); p.add_argument("--candidates")
    p.add_argument("--weighted-network"); p.add_argument("--expression")
    p.add_argument("--cpp-diagnostics-json")
    p.add_argument("--centrality-normal"); p.add_argument("--centrality-tumor")
    p.add_argument("--cgc-list")
    p.add_argument("--dice-genes"); p.add_argument("--deg-genes")
    p.add_argument("--depmap-table")
    p.add_argument("--dnf-pairwise"); p.add_argument("--synlethdb")
    args = p.parse_args(argv)

    ctx = {k.replace("-", "_"): getattr(args, k.replace("-", "_"))
           for k in vars(args) if k != "output"}

    probes: list[Probe] = [
        RawInputsProbe(), ParsersProbe(), MergeProbe(), FiltersProbe(),
        WeightingProbe(), CentralityProbe(), DifferentialProbe(),
        ValidationProbe(), DNFProbe(),
    ]
    results: list[ProbeResult] = []
    for probe in probes:
        # Skip probes whose required ctx keys are absent
        try:
            results.append(probe.run(ctx))
        except Exception as e:  # noqa: BLE001
            log.exception(f"probe {probe.name} failed unexpectedly")

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    commit = _git_commit()
    (out / "diagnostic_dossier.md").write_text(render_markdown(results, commit))
    (out / "dossier.json").write_text(render_json(results, commit))
    print(f"wrote {out}/diagnostic_dossier.md and dossier.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/diagnostics/dossier_runner.py tests/integration/test_dossier_runner.py
git commit -m "feat(diagnostics): dossier_runner CLI orchestrator [task-1.12]"
```

---

### Task 1.13: Run dossier on real PRAD data; commit registry

**Files:**
- Create: `docs/audit/diagnostic_dossier.md` (generated, then human-edited to populate the bug registry table at the bottom)
- Create: `docs/audit/dossier.json` (generated)

- [ ] **Step 1: Make sure the C++ engine has been re-built with `--diagnostics-out`**

```bash
cmake --build build --target dice_analyzer
./build/dice_analyzer \
    data/processed/weighted/network_normal_int.txt \
    data/processed/weighted/network_tumor_int.txt \
    data/results \
    --diagnostics-out data/results/cpp_diagnostics.json
```

- [ ] **Step 2: Run the dossier runner with full PRAD context**

```bash
python -m src.python.diagnostics.dossier_runner \
  --output docs/audit \
  --expression-normal data/raw/expression/expression_normal.txt \
  --expression-tumor data/raw/expression/expression_tumor.txt \
  --kegg-dir data/raw/kegg \
  --biogrid-file data/raw/biogrid/BIOGRID-ALL-5.0.253.tab2.txt \
  --kegg-interactions data/processed/networks/kegg_interactions.txt \
  --biogrid-interactions data/processed/networks/biogrid_interactions.txt \
  --merged-network data/processed/networks/merged_network.txt \
  --dea-table data/processed/dea_table.tsv \
  --dea-genes data/processed/genes_dea.txt \
  --ig-genes data/processed/genes_ig.txt \
  --candidates data/processed/genes_candidates.txt \
  --weighted-network data/processed/weighted/network_normal.txt \
  --expression data/raw/expression/expression_normal.txt \
  --cpp-diagnostics-json data/results/cpp_diagnostics.json \
  --centrality-normal data/results/centrality_normal.txt \
  --centrality-tumor data/results/centrality_tumor.txt \
  --cgc-list data/known_genes/cgc.txt \
  --dice-genes data/results/dice_genes.txt \
  --deg-genes data/processed/genes_dea.txt \
  --depmap-table data/known_genes/depmap_fitness.tsv \
  --dnf-pairwise data/results/dnf_pairwise_results.txt \
  --synlethdb data/known_genes/synlethdb.tsv
```

(Skip arguments whose files don't exist; the runner tolerates missing inputs.)

- [ ] **Step 3: Read the dossier and populate the bug registry**

Open `docs/audit/diagnostic_dossier.md`. At the bottom there is an empty `Bug Registry` table. Manually populate it from the per-stage findings. **Do this by hand** — the table is the human triage step that drives Phase 2. Each entry must include:
- a unique `BUG-NN` ID,
- the stage,
- a one-line symptom,
- a pointer to the dossier evidence,
- confidence (high/medium/low),
- blast radius (count of downstream stages affected),
- a hypothesized root cause.

This is the deliverable that drives the follow-up Phase 2 plan.

- [ ] **Step 4: Commit**

```bash
git add docs/audit/diagnostic_dossier.md docs/audit/dossier.json
git commit -m "docs(audit): generate Phase 1 diagnostic dossier + registry [task-1.13]"
```

**End of Phase 1.** At this checkpoint, hand the dossier to the user for review before starting Phase 2.

---

## PHASE 2 — Pre-committed correct-by-default fixes (6 tasks)

These five fixes are **independent of dossier findings** and were locked in §3 of the spec. Dossier-driven fixes (DEG 0-overlap, Entrez↔symbol corrections, etc.) are scheduled in a follow-up plan generated *from* the bug registry produced in Task 1.13.

### Task 2.1: `--dice-faithful` flag plumbing

Add the flag at every layer; default = off (correct-by-default). When `--dice-faithful` is set, behavior reverts to pre-audit.

**Files:**
- Modify: `run_pipeline.sh`
- Modify: `src/python/differential_centrality.py`
- Modify: `src/python/weight_network.py`
- Modify: `src/cpp/main.cpp`

- [ ] **Step 1: Write a regression placeholder test**

`tests/regression/test_dice_faithful_regression.py`:
```python
"""Phase 2 regression: --dice-faithful must reproduce pre-audit dice_genes top-500
to within rank correlation >= 0.99."""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from scipy import stats


GOLDEN = Path(__file__).parent / "golden" / "dice_genes_top500_pre_audit.tsv"


@pytest.mark.regression
def test_dice_faithful_top500_rank_correlation_ge_099():
    if not GOLDEN.exists():
        pytest.skip("Golden file missing; run scripts/capture_golden.sh")

    # Run pipeline with --dice-faithful flag (assumes pipeline has already been re-run)
    new = pd.read_csv("data/results/dice_genes_top500.txt", sep="\t")
    old = pd.read_csv(GOLDEN, sep="\t")

    common = set(new["gene"]) & set(old["gene"])
    new_ranks = new.set_index("gene").loc[list(common), "final_rank"]
    old_ranks = old.set_index("gene").loc[list(common), "final_rank"]
    rho, _ = stats.spearmanr(new_ranks, old_ranks)
    assert rho >= 0.99, f"rank correlation {rho:.4f} < 0.99 — --dice-faithful broken"
```

- [ ] **Step 2: Add flag to `differential_centrality.py`**

Modify `src/python/differential_centrality.py`. Find the `argparse` block and the `analyze()` method.

Replace the `--disable-noise` flag with `--dice-faithful` (rename flag; keep `--disable-noise` as a hidden alias for backward-compat):
```python
# In main() argparse:
parser.add_argument('--dice-faithful', action='store_true',
                    help='Reproduce Pashaei et al. behavior bit-exactly: enable top-12% '
                         'eigenvector hub removal, drop Stouffer score column.')
parser.add_argument('--disable-noise', action='store_true', help=argparse.SUPPRESS)

# After parse_args:
if args.disable_noise:
    args.dice_faithful = False  # legacy: --disable-noise was the inverse
```

In `DifferentialCentralityAnalyzer.analyze()`, replace the `disable_noise` branch with `dice_faithful`-driven logic. Default behavior (correct-by-default) **does not** apply the top-12 % eigenvector hub removal. Emit both `ensemble_score_raw` (no hub removal) and `ensemble_score_hub_corrected` (post hub removal). When `dice_faithful=True`, emit only `ensemble_score` matching the original column for the regression test.

(Specific code is part of Task 2.2.)

- [ ] **Step 3: Add flag to `weight_network.py`**

Modify `src/python/weight_network.py`:
```python
# In main() argparse:
parser.add_argument('--scale-factor', type=int, default=100000,
                    help='Integer scale factor (default 100000; use 1000 for --dice-faithful)')
parser.add_argument('--dice-faithful', action='store_true',
                    help='Use pre-audit defaults: scale=1000, no signed-correlation output')
# After parse_args:
if args.dice_faithful:
    args.scale_factor = 1000
```

(Signed-correlation matrix output added in Task 2.4.)

- [ ] **Step 4: Add flag plumbing to `run_pipeline.sh`**

In `run_pipeline.sh`, near the existing flag handling block, add:
```bash
DICE_FAITHFUL=""
SCALE_FACTOR=100000
BRANDES_MODE="auto"
# ... while [[ $# -gt 0 ]]; do
case "$1" in
    --dice-faithful) DICE_FAITHFUL="--dice-faithful"; SCALE_FACTOR=1000; BRANDES_MODE="sampled"; shift ;;
    --scale-factor)  SCALE_FACTOR="$2"; shift 2 ;;
    --brandes-mode)  BRANDES_MODE="$2"; shift 2 ;;
    # ... existing cases
esac
```
Forward `$DICE_FAITHFUL`, `$SCALE_FACTOR`, `$BRANDES_MODE` to every Python and C++ invocation that accepts them.

- [ ] **Step 5: Add flag to `main.cpp`**

```cpp
// In main.cpp argument parsing:
std::string brandes_mode = "auto";
int scale_factor = 100000;
bool dice_faithful = false;
for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--brandes-mode" && i+1 < argc) { brandes_mode = argv[++i]; }
    else if (a == "--scale-factor" && i+1 < argc) { scale_factor = std::stoi(argv[++i]); }
    else if (a == "--dice-faithful") { dice_faithful = true; brandes_mode = "sampled"; scale_factor = 1000; }
}
// scale_factor and brandes_mode are forwarded into engine calls in tasks 2.5 and 2.6.
```

- [ ] **Step 6: Build & sanity-run**

```bash
cmake --build build --target dice_analyzer
./build/dice_analyzer \
    data/processed/weighted/network_normal_int.txt \
    data/processed/weighted/network_tumor_int.txt \
    data/results \
    --dice-faithful
```

Expected: completes successfully, produces output (used as a sanity check; behavior changes land in Tasks 2.2-2.6).

- [ ] **Step 7: Commit**

```bash
git add run_pipeline.sh src/python/differential_centrality.py src/python/weight_network.py \
        src/cpp/main.cpp tests/regression/test_dice_faithful_regression.py
git commit -m "feat(audit): plumb --dice-faithful through every layer [task-2.1]"
```

---

### Task 2.2: Disable top-12 % hub removal by default; emit raw + hub-corrected

**Files:**
- Modify: `src/python/differential_centrality.py`
- Create: `tests/bugs/test_2_2_hub_removal.py`

- [ ] **Step 1: Test**

`tests/bugs/test_2_2_hub_removal.py`:
```python
"""Bug 2.2: top-12% hub removal silently kills TP53/EP300.
Default must keep them (correct-by-default); --dice-faithful drops them."""
import subprocess
import sys
import pandas as pd
from pathlib import Path


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
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

In `src/python/differential_centrality.py`, replace `analyze()` with:
```python
def analyze(self) -> pd.DataFrame:
    logger.info("Loading centrality data...")
    normal_df = pd.read_csv(self.normal_file, sep='\t')
    tumor_df = pd.read_csv(self.tumor_file, sep='\t')

    merged = pd.merge(normal_df, tumor_df, on='gene', suffixes=('_norm', '_tumor'))
    merged['diff_betweenness'] = (merged['betweenness_tumor'] - merged['betweenness_norm']).abs()
    merged['diff_eigenvector'] = (merged['eigenvector_tumor'] - merged['eigenvector_norm']).abs()
    merged['rank_bet'] = merged['diff_betweenness'].rank(pct=True)
    merged['rank_eig'] = merged['diff_eigenvector'].rank(pct=True)

    # Always compute the raw ensemble (no hub filtering)
    merged['ensemble_score_raw'] = merged['rank_bet'] * merged['rank_eig']

    # Hub-corrected: same score but applied AFTER dropping top-12% eigenvector
    eig_n = merged['eigenvector_norm'].quantile(0.88)
    eig_t = merged['eigenvector_tumor'].quantile(0.88)
    avg_bet_n = merged['betweenness_norm'].mean()
    avg_bet_t = merged['betweenness_tumor'].mean()
    mask_keep = ((merged['eigenvector_norm'] <= eig_n) & (merged['eigenvector_tumor'] <= eig_t)) & \
                ((merged['betweenness_norm'] > avg_bet_n) | (merged['betweenness_tumor'] > avg_bet_t))

    hub_df = merged[mask_keep].copy()
    hub_df['ensemble_score_hub_corrected'] = hub_df['rank_bet'] * hub_df['rank_eig']
    merged = merged.merge(
        hub_df[['gene', 'ensemble_score_hub_corrected']],
        on='gene', how='left'
    )

    if getattr(self, 'dice_faithful', False):
        # Pashaei-faithful: drop hub list, only emit ensemble_score
        out = hub_df.copy()
        out['ensemble_score'] = out['ensemble_score_hub_corrected']
        out = out.drop(columns=['ensemble_score_hub_corrected'])
        logger.warning("--dice-faithful: dropped %d top-12%% eigenvector hubs (incl. potentially TP53/EP300)",
                       len(merged) - len(out))
    else:
        # Correct-by-default: keep everyone; report both scores
        out = merged.sort_values('ensemble_score_raw', ascending=False).copy()
        logger.info("Default mode: kept all %d genes; reporting raw and hub-corrected scores", len(out))

    out['final_rank'] = range(1, len(out) + 1)
    out['delta_betweenness'] = out['betweenness_tumor'] - out['betweenness_norm']
    out['delta_eigenvector'] = out['eigenvector_tumor'] - out['eigenvector_norm']
    return out
```

And in `main()`:
```python
analyzer = DifferentialCentralityAnalyzer(args.normal, args.tumor)
analyzer.dice_faithful = bool(args.dice_faithful)
```

- [ ] **Step 4: Run, expect pass**

```bash
python -m pytest tests/bugs/test_2_2_hub_removal.py -v
```

- [ ] **Step 5: Commit**

```bash
git add src/python/differential_centrality.py tests/bugs/test_2_2_hub_removal.py
git commit -m "fix(differential): default keeps top-12pct hubs, emits raw + hub-corrected scores [task-2.2]"
```

---

### Task 2.3: Add Stouffer score alongside rank product

Bet/eig are correlated → rank-product double-counts. Add a Stouffer combined score using percentile ranks → z-scores → sum / √2.

**Files:**
- Modify: `src/python/differential_centrality.py`
- Create: `tests/bugs/test_2_3_stouffer.py`

- [ ] **Step 1: Test**

`tests/bugs/test_2_3_stouffer.py`:
```python
import subprocess, sys, pandas as pd


def test_stouffer_score_present(tmp_path):
    cn = tmp_path / "cn.txt"; ct = tmp_path / "ct.txt"
    rows_n = ["gene\tbetweenness\teigenvector"]; rows_t = ["gene\tbetweenness\teigenvector"]
    for i in range(50):
        rows_n.append(f"G{i}\t{i}\t{0.01*i}")
        rows_t.append(f"G{i}\t{i+1}\t{0.01*i}")
    cn.write_text("\n".join(rows_n)+"\n"); ct.write_text("\n".join(rows_t)+"\n")
    out = tmp_path / "dice.txt"
    r = subprocess.run([sys.executable, "src/python/differential_centrality.py",
                        "--normal", str(cn), "--tumor", str(ct), "--output", str(out),
                        "--top-n", "100"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    df = pd.read_csv(out, sep="\t")
    assert "ensemble_score_stouffer" in df.columns
    # Stouffer is a z-score; should be roughly mean 0 with both signs present
    assert df["ensemble_score_stouffer"].max() > 0
    assert df["ensemble_score_stouffer"].min() < 0
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

In `differential_centrality.py`, after `merged['rank_eig'] = ...`, add:
```python
from scipy.stats import norm
# Stouffer combined score from rank-percentiles -> z-scores
# Avoid 0/1 (inverse-normal infinity) by clipping to [eps, 1-eps]
eps = 1e-6
z_bet = norm.ppf(merged['rank_bet'].clip(eps, 1 - eps))
z_eig = norm.ppf(merged['rank_eig'].clip(eps, 1 - eps))
merged['ensemble_score_stouffer'] = (z_bet + z_eig) / (2 ** 0.5)
```

Append `'ensemble_score_stouffer'` to the column list copied through both branches (default and `--dice-faithful`).

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/differential_centrality.py tests/bugs/test_2_3_stouffer.py
git commit -m "feat(differential): add Stouffer combined score (decorrelated alt to rank-product) [task-2.3]"
```

---

### Task 2.4: Signed-correlation matrix output in `weight_network.py`

The DiCE `1−|r|` weighting drops the sign. We continue to use it for SSSP, *but additionally* emit a separate two-column file: `(source, target, signed_correlation)` for downstream analysis.

**Files:**
- Modify: `src/python/weight_network.py`
- Create: `tests/bugs/test_2_4_signed_correlation.py`

- [ ] **Step 1: Test**

`tests/bugs/test_2_4_signed_correlation.py`:
```python
import subprocess, sys, pandas as pd, pathlib


def test_signed_correlation_file_emitted(tmp_path):
    net = tmp_path / "net.txt"
    net.write_text("source\ttarget\nA\tB\nA\tC\n")
    expr = tmp_path / "expr.tsv"
    expr.write_text(
        "gene\tS1\tS2\tS3\tS4\tS5\n"
        "A\t1\t2\t3\t4\t5\n"
        "B\t5\t4\t3\t2\t1\n"
        "C\t1\t2\t3\t4\t5\n"
    )
    out = tmp_path / "weighted.txt"
    r = subprocess.run([sys.executable, "src/python/weight_network.py",
                        "--network", str(net), "--expression", str(expr),
                        "--output", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr

    signed = pathlib.Path(str(out).replace(".txt", "_signed.tsv"))
    assert signed.exists()
    df = pd.read_csv(signed, sep="\t")
    # A-B is anti-correlated (~ -1), A-C is positively correlated (~ +1)
    ab = df[(df["source"]=="A") & (df["target"]=="B")]["signed_correlation"].iloc[0]
    ac = df[(df["source"]=="A") & (df["target"]=="C")]["signed_correlation"].iloc[0]
    assert ab < -0.9
    assert ac > 0.9
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

In `weight_network.py`, modify `_calculate_weight` to also return the signed correlation, and accumulate them in parallel:

Replace `_calculate_weight` body:
```python
def _calculate_weight(self, expr_df, gene1, gene2):
    expr1 = expr_df.loc[gene1].values
    expr2 = expr_df.loc[gene2].values
    corr, _ = pearsonr(expr1, expr2)
    if np.isnan(corr):
        corr = 0.0
    weight = max(1.0 - abs(corr), self.epsilon)
    return weight, float(corr)
```

In `weight_network()`, accumulate two columns:
```python
weights = []
signed = []
for idx, row in network_df.iterrows():
    w, c = self._calculate_weight(expr_df, row["source"], row["target"])
    weights.append(w); signed.append(c)
network_df["weight"] = weights
network_df["signed_correlation"] = signed
```

After saving the main + int outputs in `main()`, also save the signed file:
```python
signed_df = weighted_df[["source", "target", "signed_correlation"]]
signed_path = args.output.replace(".txt", "_signed.tsv")
signed_df.to_csv(signed_path, sep="\t", index=False)
logger.info(f"Saved signed-correlation matrix to {signed_path}")
```

(Keep `weighted_df = weighter.weight_network()` — the new `signed_correlation` column flows through.)

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/weight_network.py tests/bugs/test_2_4_signed_correlation.py
git commit -m "feat(weighting): emit signed-correlation matrix alongside |r|-distance weights [task-2.4]"
```

---

### Task 2.5: Brandes — full when n<5000, bootstrap-CI when sampled

**Files:**
- Modify: `src/cpp/engine.hpp` (add bootstrap path)
- Modify: `src/cpp/main.cpp` (use `--brandes-mode`)
- Modify: `src/cpp/tests/test_engine.cpp`

- [ ] **Step 1: Add C++ tests**

Append to `src/cpp/tests/test_engine.cpp`:
```cpp
TEST_CASE("brandes_with_ci returns CI bounds spanning the point estimate") {
    Graph<long long> g(10);
    for (int i = 0; i < 9; ++i) g.add_edge(i, i+1, 1);
    g.id_to_name.resize(10);

    CentralityCalculator<long long> calc;
    auto res = calc.brandes_with_ci(g, /*sample_rate=*/0.3, /*B=*/20, /*seed=*/42);
    CHECK(res.point.size() == 10);
    CHECK(res.ci_lo.size() == 10);
    CHECK(res.ci_hi.size() == 10);
    for (size_t i = 0; i < res.point.size(); ++i) {
        CHECK(res.ci_lo[i] <= res.point[i] + 1e-6);
        CHECK(res.ci_hi[i] >= res.point[i] - 1e-6);
    }
}
```

- [ ] **Step 2: Build, expect compile failure** (`brandes_with_ci` doesn't exist yet)

- [ ] **Step 3: Implement `brandes_with_ci` in `engine.hpp`**

Inside `class CentralityCalculator<WeightType>`:
```cpp
struct BrandesCI {
    std::vector<double> point;
    std::vector<double> ci_lo;
    std::vector<double> ci_hi;
};

BrandesCI brandes_with_ci(Graph<WeightType>& g, double sample_rate, int B, uint64_t seed) {
    // B bootstrap samples; per-vertex 95% CI from quantiles.
    std::vector<std::vector<double>> samples(g.n_vertices);
    for (auto& v : samples) v.reserve(B);

    for (int b = 0; b < B; ++b) {
        auto bet = betweenness(g, sample_rate, seed + (uint64_t)b);
        for (int v = 0; v < g.n_vertices; ++v) samples[v].push_back(bet[v]);
    }

    BrandesCI out;
    out.point.resize(g.n_vertices);
    out.ci_lo.resize(g.n_vertices);
    out.ci_hi.resize(g.n_vertices);
    for (int v = 0; v < g.n_vertices; ++v) {
        std::sort(samples[v].begin(), samples[v].end());
        out.point[v]  = samples[v][B / 2];
        out.ci_lo[v]  = samples[v][std::max(0, int(B * 0.025))];
        out.ci_hi[v]  = samples[v][std::min(B - 1, int(B * 0.975))];
    }
    return out;
}
```

- [ ] **Step 4: Wire into `main.cpp`**

In `main.cpp`, after `Graph<long long> tumor_g = ...; Graph<long long> normal_g = ...;`, replace the existing centrality call with mode-switch:
```cpp
auto compute = [&](Graph<long long>& g, const std::string& tag) {
    CentralityCalculator<long long> c;
    if (brandes_mode == "full" || (brandes_mode == "auto" && g.n_vertices < 5000)) {
        return c.betweenness(g, 1.0);
    } else if (brandes_mode == "bootstrap") {
        auto r = c.brandes_with_ci(g, 0.1, /*B=*/20, /*seed=*/42);
        // Persist CI alongside point estimate to side files
        std::ofstream lo("data/results/centrality_" + tag + "_ci_lo.txt"),
                      hi("data/results/centrality_" + tag + "_ci_hi.txt");
        for (int v = 0; v < g.n_vertices; ++v) {
            lo << g.id_to_name[v] << "\t" << r.ci_lo[v] << "\n";
            hi << g.id_to_name[v] << "\t" << r.ci_hi[v] << "\n";
        }
        return r.point;
    } else {  // "sampled" — DiCE-faithful
        return c.betweenness(g, 0.01, 42);
    }
};
auto bet_n = compute(normal_g, "normal");
auto bet_t = compute(tumor_g, "tumor");
```

- [ ] **Step 5: Build, run C++ tests, expect pass**

- [ ] **Step 6: Commit**

```bash
git add src/cpp/engine.hpp src/cpp/main.cpp src/cpp/tests/test_engine.cpp
git commit -m "feat(engine): brandes_with_ci + auto/full/bootstrap/sampled modes [task-2.5]"
```

---

### Task 2.6: Configurable scale factor + collision audit

**Files:**
- Modify: `src/python/weight_network.py`
- Create: `tests/bugs/test_2_6_scale_factor.py`

- [ ] **Step 1: Test**

`tests/bugs/test_2_6_scale_factor.py`:
```python
import subprocess, sys, pandas as pd


def test_scale_factor_100k_fewer_collisions_than_1k(tmp_path):
    net = tmp_path / "net.txt"
    net.write_text("source\ttarget\n" + "\n".join(
        f"G{i}\tG{i+1}" for i in range(50)
    ) + "\n")
    expr_lines = ["gene\tS1\tS2\tS3\tS4\tS5"]
    rng = __import__("random").Random(20260507)
    for i in range(52):
        expr_lines.append(f"G{i}\t" + "\t".join(f"{rng.uniform(0,10):.4f}" for _ in range(5)))
    expr = tmp_path / "expr.tsv"; expr.write_text("\n".join(expr_lines) + "\n")

    def collisions(scale):
        out = tmp_path / f"w_{scale}.txt"
        subprocess.run([sys.executable, "src/python/weight_network.py",
                        "--network", str(net), "--expression", str(expr),
                        "--output", str(out), "--scale-factor", str(scale)],
                       check=True, capture_output=True, text=True)
        df = pd.read_csv(out, sep="\t")
        return df["weight"].nunique() - df["weight_int"].nunique()

    assert collisions(100000) <= collisions(1000)
```

- [ ] **Step 2: Run, expect failure** (current code uses fixed `--scale 1000`; new flag `--scale-factor` doesn't exist yet, or it's only 1000)

- [ ] **Step 3: Implement**

In `weight_network.py`, the argparse already added `--scale-factor` in Task 2.1. Now ensure the constructor uses it (already in current code: `scale_factor=1000` default). Just change the default in argparse to 100000 (already done in Task 2.1) and confirm the int conversion path still scales by `args.scale_factor`.

Add a collision audit at end of `weight_network()`:
```python
n_unique_w = network_df['weight'].nunique()
n_unique_i = network_df['weight_int'].nunique()
collisions = n_unique_w - n_unique_i
logger.info(
    f"scale_factor={self.scale_factor}: weight_int collisions = {collisions} "
    f"(unique weight={n_unique_w}, unique int={n_unique_i})"
)
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add src/python/weight_network.py tests/bugs/test_2_6_scale_factor.py
git commit -m "feat(weighting): default scale=100000, collision audit logged [task-2.6]"
```

---

## PHASE 3 — Final gate (3 tasks)

### Task 3.1: Full PRAD re-run script

**Files:**
- Create: `scripts/run_audit_rerun.sh`

- [ ] **Step 1: Write the script**

```bash
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p data/results/audit_rerun
echo "== Step 1: Rebuild C++ engine =="
cmake --build build --target dice_analyzer

echo "== Step 2: Run full pipeline (correct-by-default) =="
./run_pipeline.sh --scale-factor 100000 --brandes-mode auto

echo "== Step 3: Run full pipeline (--dice-faithful for regression) =="
DEST=data/results/audit_rerun_dice_faithful
mkdir -p "$DEST"
cp data/results/dice_genes.txt "$DEST/dice_genes_default.txt" || true
./run_pipeline.sh --dice-faithful
cp data/results/dice_genes.txt "$DEST/dice_genes_faithful.txt"

echo "== Step 4: Generate before/after report =="
python scripts/before_after.py --golden tests/regression/golden \
    --new data/results --out docs/audit/before_after.md

echo "Done. Review docs/audit/before_after.md"
```

`chmod +x scripts/run_audit_rerun.sh`.

- [ ] **Step 2: Commit (without running yet)**

```bash
git add scripts/run_audit_rerun.sh
git commit -m "chore(audit): add full-pipeline re-run driver script [task-3.1]"
```

---

### Task 3.2: Before/after delta report generator (TDD)

**Files:**
- Create: `scripts/before_after.py`
- Create: `tests/unit/test_before_after.py`

- [ ] **Step 1: Test**

`tests/unit/test_before_after.py`:
```python
import pandas as pd
from scripts.before_after import compare_dice_genes


def test_compare_dice_genes_reports_top50_jaccard(tmp_path):
    old = pd.DataFrame({"gene": [f"G{i}" for i in range(100)],
                        "final_rank": list(range(1, 101))})
    new = pd.DataFrame({"gene": [f"G{i}" for i in range(50, 150)],
                        "final_rank": list(range(1, 101))})
    p_old = tmp_path / "old.tsv"; p_new = tmp_path / "new.tsv"
    old.to_csv(p_old, sep="\t", index=False); new.to_csv(p_new, sep="\t", index=False)
    rep = compare_dice_genes(str(p_old), str(p_new))
    assert 0.0 < rep["top50_jaccard"] < 1.0
    assert "G0" in rep["entered"] or "G149" in rep["entered"] or len(rep["entered"]) >= 1
```

- [ ] **Step 2: Run, expect failure**

- [ ] **Step 3: Implement**

`scripts/before_after.py`:
```python
"""Generate docs/audit/before_after.md from old vs new artifacts."""
from __future__ import annotations
import argparse
from pathlib import Path
from typing import Dict, Any
import pandas as pd
from scipy import stats


def compare_dice_genes(old_path: str, new_path: str) -> Dict[str, Any]:
    old = pd.read_csv(old_path, sep="\t")
    new = pd.read_csv(new_path, sep="\t")
    old_top50 = set(old.head(50)["gene"]); new_top50 = set(new.head(50)["gene"])
    common = set(old["gene"]) & set(new["gene"])
    old_r = old.set_index("gene").loc[list(common), "final_rank"]
    new_r = new.set_index("gene").loc[list(common), "final_rank"]
    rho, _ = stats.spearmanr(old_r, new_r) if len(common) > 1 else (float("nan"), float("nan"))
    return {
        "n_old": len(old), "n_new": len(new),
        "top50_jaccard": len(old_top50 & new_top50) / max(1, len(old_top50 | new_top50)),
        "rank_correlation_common": float(rho) if rho == rho else None,
        "entered": sorted(new_top50 - old_top50),
        "exited": sorted(old_top50 - new_top50),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--golden", required=True)
    p.add_argument("--new", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    g = Path(args.golden); n = Path(args.new); out = Path(args.out)
    rep = compare_dice_genes(g / "dice_genes_pre_audit.tsv", n / "dice_genes.txt")

    md = [
        "# Before / After — Audit re-run", "",
        f"- old (`{g}`): n={rep['n_old']}",
        f"- new (`{n}`): n={rep['n_new']}",
        f"- **top-50 Jaccard**: {rep['top50_jaccard']:.3f}",
        f"- **rank correlation (common)**: {rep['rank_correlation_common']}",
        "",
        "## Entered top-50 (post-audit)", "",
        *(f"- `{g}`" for g in rep["entered"]),
        "",
        "## Exited top-50 (post-audit)", "",
        *(f"- `{g}`" for g in rep["exited"]),
        "",
        "## Audit explanation",
        "_For every entered/exited gene above, link to the dossier evidence and the commit that introduced the change._",
    ]
    out.write_text("\n".join(md) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run, expect pass**

- [ ] **Step 5: Commit**

```bash
git add scripts/before_after.py tests/unit/test_before_after.py
git commit -m "feat(audit): before/after report generator [task-3.2]"
```

---

### Task 3.3: Run final gate, write `before_after.md`, run regression contract

- [ ] **Step 1: Run all tests locally first**

```bash
./scripts/run_tests.sh
```
Expected: all green.

- [ ] **Step 2: Run the full audit re-run**

```bash
./scripts/run_audit_rerun.sh
```
Expected: completes; `docs/audit/before_after.md` written.

- [ ] **Step 3: Run the regression contract test**

```bash
python -m pytest tests/regression/test_dice_faithful_regression.py -v
```
Expected: `test_dice_faithful_top500_rank_correlation_ge_099` passes (≥ 0.99).

If it fails, this is itself a finding: `--dice-faithful` does not faithfully reproduce. Investigate before merging — most likely candidates: scale factor not threaded through, sampled-Brandes not deterministic with seed, hub-removal mask drift from FP arithmetic.

- [ ] **Step 4: Final commit**

```bash
git add docs/audit/before_after.md docs/audit/diagnostic_dossier.md docs/audit/dossier.json
git commit -m "docs(audit): final-gate before/after report + regenerated dossier [task-3.3]"
git tag audit-phase1-complete
```

---

## Self-review

### Spec coverage
| Spec section | Plan task(s) | Status |
|---|---|---|
| §1 Definition of done — diagnostic dossier exists | 1.13 | ✓ |
| §1 — every fix has a test that fails pre-fix, passes after | TDD discipline + Tasks 2.2–2.6, 3.2 | ✓ |
| §1 — pytest + C++ test target green via single driver | 0.1, 0.2, 3.3 | ✓ |
| §1 — full PRAD re-run, before/after report | 3.1, 3.2, 3.3 | ✓ |
| §1 — `--dice-faithful` regression at top-500 ≥ 0.99 | 0.5, 2.1, 3.3 | ✓ |
| §1 — no silent transformations | safe_filter (0.3) + applied across modified modules in Phase 2 | ✓ (utility done; full application in Phase-2 follow-up plan) |
| §2 — 10 boundary probes | 1.1–1.10 | ✓ (probe-9/export checksum is rolled into stage 7 in this plan; if needed it can become a 12th task in the follow-up) |
| §3 — top-12 % hub removal off by default | 2.2 | ✓ |
| §3 — signed-correlation matrix output | 2.4 | ✓ |
| §3 — Stouffer score | 2.3 | ✓ |
| §3 — Brandes full / sampled / bootstrap | 2.5 | ✓ |
| §3 — scale ×100000 + collision audit | 2.6 | ✓ |
| §3 — `--dice-faithful` regression contract | 0.5, 2.1, 3.3 | ✓ |
| §4 — file structure (new + modified) | locked at top of this plan | ✓ |
| §4 — test layers 1–5 | 0.1, 0.2, 1.x, 2.x, 3.3 | ✓ |

**Gap addressed:** §3 dossier-driven Phase 2 fixes (DEG 0-overlap, Entrez↔symbol full sweep, NaN/empty silent drops elimination, sampled-Brandes-bias-stated, validation join correctness) are scheduled in a follow-up plan, generated *after* the registry exists. This is intentional — see plan preamble.

### Placeholder scan
- All test code blocks contain real assertions and inputs.
- All implementation code blocks are minimal but complete; no `# TODO` or `# fill in`.
- The bug registry table in Task 1.13 is intentionally human-edited from dossier evidence; this is a triage action, not a placeholder.

### Type-consistency check
- `safe_filter(df, mask, reason, stage, strict=False) -> FilterResult` — used same way wherever invoked.
- `IDResolver.resolve(key) -> Optional[ResolveResult]` — same return type everywhere.
- `Probe.run(ctx) -> ProbeResult` — same signature in all probes.
- `betweenness(g, sample_rate)` keeps original signature; new `betweenness(g, sample_rate, seed)` overload added without changing existing callers.
- C++ flags `--dice-faithful`, `--scale-factor`, `--brandes-mode`, `--diagnostics-out` consistent across `main.cpp`, `run_pipeline.sh`, and tests.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-07-audit-and-fix-phase1.md`.

A separate plan `2026-MM-DD-audit-and-fix-phase2.md` will be generated **after Task 1.13 produces the bug registry** — that is when we know which dossier-driven fixes (DEG 0-overlap, ID sweep, validation joins, etc.) need TDD code.
