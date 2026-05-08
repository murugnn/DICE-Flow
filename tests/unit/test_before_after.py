import pandas as pd

from scripts.before_after import compare_dice_genes


def test_compare_dice_genes_reports_top50_jaccard(tmp_path):
    # Note: plan's fixture used range(50, 150) for new which makes the top-50s
    # fully disjoint (jaccard==0), violating the strict 0<j<1 assertion. Shift
    # to range(25, 125) so top-50s share G25..G49 → jaccard = 25/75 ≈ 0.33.
    old = pd.DataFrame({"gene": [f"G{i}" for i in range(100)],
                        "final_rank": list(range(1, 101))})
    new = pd.DataFrame({"gene": [f"G{i}" for i in range(25, 125)],
                        "final_rank": list(range(1, 101))})
    p_old = tmp_path / "old.tsv"
    p_new = tmp_path / "new.tsv"
    old.to_csv(p_old, sep="\t", index=False)
    new.to_csv(p_new, sep="\t", index=False)
    rep = compare_dice_genes(str(p_old), str(p_new))
    assert 0.0 < rep["top50_jaccard"] < 1.0
    assert "G0" in rep["entered"] or "G149" in rep["entered"] or len(rep["entered"]) >= 1
