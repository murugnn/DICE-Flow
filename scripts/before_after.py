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
    old_top50 = set(old.head(50)["gene"])
    new_top50 = set(new.head(50)["gene"])
    common = set(old["gene"]) & set(new["gene"])
    if len(common) > 1:
        old_r = old.set_index("gene").loc[list(common), "final_rank"]
        new_r = new.set_index("gene").loc[list(common), "final_rank"]
        rho, _ = stats.spearmanr(old_r, new_r)
    else:
        rho = float("nan")
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

    g = Path(args.golden)
    n = Path(args.new)
    out = Path(args.out)
    rep = compare_dice_genes(str(g / "dice_genes_pre_audit.tsv"), str(n / "dice_genes.txt"))

    md = [
        "# Before / After — Audit re-run", "",
        f"- old (`{g}`): n={rep['n_old']}",
        f"- new (`{n}`): n={rep['n_new']}",
        f"- **top-50 Jaccard**: {rep['top50_jaccard']:.3f}",
        f"- **rank correlation (common)**: {rep['rank_correlation_common']}",
        "",
        "## Entered top-50 (post-audit)", "",
        *(f"- `{name}`" for name in rep["entered"]),
        "",
        "## Exited top-50 (post-audit)", "",
        *(f"- `{name}`" for name in rep["exited"]),
        "",
        "## Audit explanation",
        "_For every entered/exited gene above, link to the dossier evidence and the commit that introduced the change._",
    ]
    out.write_text("\n".join(md) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
