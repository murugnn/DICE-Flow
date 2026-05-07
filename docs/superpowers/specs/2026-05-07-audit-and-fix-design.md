# Audit & Fix — Design Spec

**Date:** 2026-05-07
**Sub-project:** Audit & Fix (1 of N for the DiCE-Duan publication track)
**Approach:** Reproduce-first / two-phase, correct-by-default
**Time budget:** 2–3 weeks
**Status:** Approved by user; awaiting spec review.

---

## 1. Goal & scope

### Goal
Produce a TCGA-PRAD pipeline whose every numerical output is defensible at peer review. No fabricated, hidden, or silently-altered values. Every contested DiCE design choice is documented and either corrected or flagged behind `--dice-faithful`.

### In scope
- All four bug categories: ID/data, statistical/ranking, numerical/perf, biology/validation.
- Two-phase process: (1) instrument + diagnostic dossier, (2) TDD fixes.
- Correct-by-default for DiCE choices that are statistically problematic; DiCE-faithful behavior available behind explicit flags.
- One full PRAD re-run after Phase 2; before/after delta tables for every primary artifact.

### Out of scope (deferred to other sub-projects)
- DNF formalization beyond keeping its current behavior intact.
- Multi-cohort runs (BRCA / LUAD / COAD / GBM) — single-cohort PRAD only.
- Pearson → MI / dCor methodology replacement.
- Docker / CI / SSSP-vs-Boost benchmarks.
- Paper text. The dossier and the before/after report are paper-quality supplementary, but the manuscript is its own sub-project.

### Definition of done
1. `docs/audit/diagnostic_dossier.md` exists, with one finding per bug, evidence, and a reproduction command.
2. Every fix has at least one test that fails on pre-fix code and passes after.
3. `pytest` and the C++ test target both green via `scripts/run_tests.sh`.
4. Full PRAD re-run completes; `docs/audit/before_after.md` reports deltas for every primary artifact, with explanations for non-trivial changes.
5. `--dice-faithful` reproduces pre-audit behavior to a documented tolerance (top-500 rank correlation ≥ 0.99); this is a regression test.
6. No silent transformations: every filter logs counts in/out with a reason string.

### Open input
`finalreview.pdf` cannot be rendered on this host (no `poppler-utils`). User to either install poppler or paste bug-relevant sections so its findings can be folded into the dossier. Spec does not block on this; registry is append-only.

---

## 2. Phase 1 — Hunt / Diagnostic Dossier

### Output artifacts
- `docs/audit/diagnostic_dossier.md` — human-readable per-stage report.
- `docs/audit/dossier.json` — machine-readable measurements; future runs diff against this.

### Mechanism
A new `src/python/diagnostics/` package with one probe per pipeline boundary plus `dossier_runner.py` that runs them all from a single command and produces both files. Each probe is pure-observation: it only records, never mutates pipeline outputs. C++ side: a `--diagnostics` flag in `main.cpp` emits boundary stats (Brandes-source distribution, sampled-vs-full correlation on a 1k subgraph, OMP determinism check) without changing pipeline numerics.

### Boundaries probed (10 stages)

1. **Raw inputs** — KEGG XML, BioGRID TAB2, expression normal/tumor: row counts, ID-type detection (Entrez vs symbol vs mixed), NaN cells, sample counts.
2. **`parse_kegg`, `parse_biogrid`** — edge counts, unique genes, directionality preservation.
3. **`merge_networks`** — directed-vs-undirected breakdown (origin of directed-centrality validity).
4. **`dea_filter`, `pre_filter`, intersection** — DEG count, FDR-distribution sanity, IG distribution, intersection size; Welch-vs-Mann-Whitney sensitivity on a sample.
5. **`weight_network`** — *the 73 % edge drop must be explained here* — counts of: edges with NaN correlation, edges hitting the ε=1e-6 clamp, integer-truncation collisions, gene-symbol resolution rate against the expression matrix.
6. **C++ engine** — Brandes sampled-vs-full Spearman ρ on a 1k-node subgraph; OMP determinism (same input × 5 runs, max relative drift); directed-graph handling verified against a toy graph with known analytical centrality.
7. **`differential_centrality`** — bet/eig rank correlation (rank-product double-count test); list of the top-12 % eigenvector genes the noise filter kills, intersected with COSMIC Cancer Gene Census (this is the TP53 / EP300 evidence).
8. **Validation pipelines** — DEG-overlap reproduction (the "0 overlap" claim), DepMap join match-rate, survival cohort size at each filter step, hypergeometric-test denominator sanity.
9. **Export pipeline** — figures regenerated checksum-stable.
10. **DNF** — pair count, DSS distribution, SynLethDB ID-resolution rate.

### Phase 1 deliverable: ranked bug registry
At the bottom of the dossier:

| ID | Stage | Symptom | Evidence | Confidence | Blast radius | Hypothesized root cause |

Sort key: `impact × confidence × blast_radius`. Drives Phase 2 triage. **Nothing in Phase 1 modifies pipeline outputs.**

### Reproducibility
Single command: `python -m diagnostics.dossier_runner --output docs/audit/`. Deterministic: every probe seeds RNGs and pins thread counts.

---

## 3. Phase 2 — TDD Fix Process + Correct-by-Default

### TDD discipline
Every fix:
1. From dossier evidence, write a failing test (unit or integration) reproducing the bug deterministically.
2. Fix the smallest amount of code that turns the test green.
3. Run affected stage's full suite + a focused re-run of the downstream stage.
4. Commit atomically. One commit = one bug-id. Message: `fix(<stage>): <one-line> [bug-#NN]`.

### Triage
Bug registry sorted by `impact × confidence × blast_radius`. Highest-impact first. Items in same category batched only if their fixes don't interact.

### Phase gates (4 + 1)
After each cluster of related fixes, run a partial re-run from the earliest affected stage through `dice_genes.txt`. Compare to previous gate's output: rank correlation, gene-set Jaccard, top-50 stability. Any regression on a previously-passing test = stop and investigate.

- **Gate 1**: ID/data fixes → re-run from `weight_network` onward.
- **Gate 2**: statistical/ranking → re-run `differential_centrality` + validation.
- **Gate 3**: numerical/perf → re-run C++ engine.
- **Gate 4**: biology/validation → re-run validation suite.
- **Final gate**: full PRAD pipeline + before/after delta report.

### Correct-by-default policy

For each contested DiCE choice, the audit installs a statistically-justified default and a `--dice-faithful` flag that restores Pashaei et al. behavior bit-exactly (or to ≥ 0.99 rank correlation if numerical changes propagate):

| DiCE behavior | New default | Justification |
|---|---|---|
| Top-12 % eigenvector hub removal | **Off**; report both raw and hub-corrected ensemble | Drops TP53 / EP300 silently; documented in dossier with CGC overlap |
| Pearson `1−|r|` (sign-blind) | Keep `1−|r|` for SSSP weights, **and** additionally output a signed-correlation matrix for downstream | Anti-correlated regulators are biologically meaningful; sign loss is the gap |
| Rank-product over correlated bet/eig | Report rank-product **and** a decorrelated Stouffer score | Bet/eig are positively correlated; rank-product double-counts |
| Sampled Brandes at `n/100` | Bootstrap-CI sampled estimate, with full Brandes when n < 5000 | Sampling bias currently un-quantified |
| Integer scale ×1000 | ×100000, with collision audit | ×1000 collapses distinct edges (verified in dossier) |

### Regression contract
`--dice-faithful` produces output whose top-500 rank correlation with pre-audit `dice_genes.txt` is ≥ 0.99. This contract is itself a test in `tests/regression/test_dice_faithful_regression.py`.

### Final PRAD re-run
One full pipeline execution at the end of Phase 2. Deltas reported in `docs/audit/before_after.md` for:
- `dice_genes.txt` — rank correlation, top-50 Jaccard, gene-set entry/exit list with rationale.
- `centrality_*` — distribution shifts.
- `knockout_analysis` — top-vitality changes.
- All `validation/*` tables.
- `dnf_pairwise_results` — DSS-rank stability.

---

## 4. Components, file layout, tests, risks

### New modules
- `src/python/diagnostics/` — `dossier_runner.py` + one probe per boundary.
- `src/python/utils/id_resolver.py` — single canonical HGNC-based resolver, used by every join.
- `src/python/utils/safe_filter.py` — wrapper that turns every "drop rows" call into a logged event with in/out counts + reason. No silent drops anywhere.
- `tests/` — `unit/`, `integration/`, `bugs/` (one file per registry entry), `regression/`.
- `src/cpp/tests/` — doctest harness wired into CMake; verifies Brandes / eigenvector against analytical truth on toy directed graphs.

### Existing files modified (touch list)
**Python:** `dea_filter`, `pre_filter`, `parse_kegg`, `parse_biogrid`, `merge_networks`, `weight_network`, `differential_centrality`, `validate_candidates`, `validation_pipeline`, `survival_validation`, `fitness_analysis`, `dnf_validation*`. Each gets logging via `safe_filter` + ID resolution via `id_resolver`. `differential_centrality` gains `--dice-faithful` and the Stouffer-score branch.

**C++:** `engine.hpp` (configurable scale, Brandes sampling rate + bootstrap CI, OMP correctness audit), `main.cpp` (`--diagnostics`, `--brandes-mode`, `--scale-factor`).

**Shell / build:** `run_pipeline.sh` propagates new flags, pins seeds. `CMakeLists.txt` adds doctest target.

### Test layers
1. **Unit** — utility functions + every probe on synthetic input.
2. **Integration** — 50-node × 200-edge micro-pipeline, < 30 s, asserts no NaN survives any filter, ID-resolution ≥ 0.95, rank correlations within expected ranges.
3. **Regression** — `test_dice_faithful_regression.py`: with `--dice-faithful`, top-500 rank correlation ≥ 0.99 against a committed golden file from current outputs.
4. **Bug-specific** — every registry entry gets a test that fails pre-fix, passes post-fix.
5. **C++** — Brandes / eigenvector on hand-computed directed toy graphs; OMP determinism (5×, max relative drift < 1e-10 with seeded RNG).

Driver: `scripts/run_tests.sh`. CI itself is a later sub-project.

### Risks & mitigations
- **Full PRAD re-run > 4 h** → profile after Gate 3, drop to stratified subset for final gate, document.
- **Top-50 changes substantially** → expected; before/after report explains every entry/exit.
- **`--dice-faithful` not bit-exact** under scale change → add `--legacy-scale` sub-flag if needed.
- **DEG 0-overlap may not be a bug** → dossier reports hypergeometric expected vs observed; only flag as bug if observed < expected.
- **No OMP race actually present** → no atomic added; we still prove it via the determinism test.
- **`finalreview.pdf` adds late items** → registry is append-only; Phase 2 absorbs them.

### Explicitly out of scope
DNF math formalization, multi-cohort runs, Pearson → MI replacement, Docker / CI, paper text.

---

## 5. Spec self-review

- **Placeholders:** none. Every numeric (12 %, 1e-6, ×1000, n/100, 0.99, 0.95) is grounded in observed code or stated as a target.
- **Internal consistency:** correct-by-default table (§3) matches the probes that produce its evidence (§2). Touch-list (§4) covers every module referenced by probes (§2) and fixes (§3).
- **Scope check:** single sub-project, one cohort, one approach. Decomposition agreed at parent level. Other sub-projects explicitly listed as out-of-scope.
- **Ambiguity check:**
  - "Bit-exact reproducibility" softened to "top-500 rank correlation ≥ 0.99" with `--legacy-scale` fallback to avoid an unachievable contract under scale changes.
  - "No silent drops" operationalized as `safe_filter` wrapper, not a vague intent.
  - "Phase gate" defined as a concrete partial re-run + comparison, not a meeting checkpoint.

---

## 6. Hand-off

Once user approves this spec, hand off to the **writing-plans** skill to produce the implementation plan. Plan should sequence Phase 1 probes, Phase 2 fixes (ordered by registry triage), and the four phase-gate re-runs.
