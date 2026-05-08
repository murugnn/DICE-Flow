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
