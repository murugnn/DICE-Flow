"""Phase-3 regression: `./run_pipeline.sh --global` must reproduce the
captured reference dice_genes top-500 to within rank correlation >= 0.99.

The reference was recaptured at the close of Phase 3 (commit recorded in
`tests/regression/golden/_captured_at_commit.txt`); see
`docs/audit/before_after.md` for why this contract replaces the original
`--dice-faithful` framing (which compared apples to oranges: pre-audit was
captured under `--global --disable-noise`, but `--dice-faithful` does
neither)."""
from pathlib import Path

import pandas as pd
import pytest
from scipy import stats


REFERENCE = Path(__file__).parent / "golden" / "dice_genes_top500_reference.tsv"


@pytest.mark.regression
def test_global_default_top500_rank_correlation_ge_099():
    if not REFERENCE.exists():
        pytest.skip("Reference file missing; run scripts/capture_golden.sh after a vetted ./run_pipeline.sh --global run.")

    current = Path("data/results/dice_genes_top500.txt")
    if not current.exists():
        pytest.skip("data/results/dice_genes_top500.txt missing; run ./run_pipeline.sh --global")

    new = pd.read_csv(current, sep="\t")
    ref = pd.read_csv(REFERENCE, sep="\t")

    common = set(new["gene"]) & set(ref["gene"])
    if not common:
        pytest.fail("No genes in common between current and reference top-500")
    new_ranks = new.set_index("gene").loc[list(common), "final_rank"]
    ref_ranks = ref.set_index("gene").loc[list(common), "final_rank"]
    rho, _ = stats.spearmanr(new_ranks, ref_ranks)
    assert rho >= 0.99, (
        f"rank correlation {rho:.4f} < 0.99 — pipeline outputs drifted vs "
        f"reference. If this drift is intentional, recapture via "
        f"scripts/capture_golden.sh after a vetted --global run."
    )
