"""Phase 2 regression: --dice-faithful must reproduce pre-audit dice_genes top-500
to within rank correlation >= 0.99."""
from pathlib import Path

import pandas as pd
import pytest
from scipy import stats


GOLDEN = Path(__file__).parent / "golden" / "dice_genes_top500_pre_audit.tsv"


@pytest.mark.regression
def test_dice_faithful_top500_rank_correlation_ge_099():
    if not GOLDEN.exists():
        pytest.skip("Golden file missing; run scripts/capture_golden.sh")

    current = Path("data/results/dice_genes_top500.txt")
    if not current.exists():
        pytest.skip("data/results/dice_genes_top500.txt missing; run pipeline with --dice-faithful")

    new = pd.read_csv(current, sep="\t")
    old = pd.read_csv(GOLDEN, sep="\t")

    common = set(new["gene"]) & set(old["gene"])
    if not common:
        pytest.fail("No genes in common between current and golden top-500")
    new_ranks = new.set_index("gene").loc[list(common), "final_rank"]
    old_ranks = old.set_index("gene").loc[list(common), "final_rank"]
    rho, _ = stats.spearmanr(new_ranks, old_ranks)
    assert rho >= 0.99, f"rank correlation {rho:.4f} < 0.99 — --dice-faithful broken"
