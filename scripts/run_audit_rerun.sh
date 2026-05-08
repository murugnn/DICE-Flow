#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

mkdir -p data/results/audit_rerun
echo "== Step 1: Rebuild C++ engine =="
cmake --build build --target dice_analyzer

echo "== Step 2: Filtered run (correct-by-default; demonstrates audit choices) =="
./run_pipeline.sh --scale-factor 100000 --brandes-mode auto
DEST=data/results/audit_rerun
cp data/results/dice_genes.txt "$DEST/dice_genes_filtered_default.txt"

echo "== Step 3: --dice-faithful run (Pashaei et al. defaults; for comparison) =="
./run_pipeline.sh --dice-faithful
cp data/results/dice_genes.txt "$DEST/dice_genes_faithful.txt"

echo "== Step 4: --global run (regression contract: must match reference) =="
./run_pipeline.sh --global
cp data/results/dice_genes.txt "$DEST/dice_genes_global.txt"

echo "== Step 5: Generate before/after report =="
python scripts/before_after.py --golden tests/regression/golden \
    --new data/results --out docs/audit/before_after.md

echo "Done. Review docs/audit/before_after.md. Final data/results/dice_genes_top500.txt is from --global run (regression target)."
