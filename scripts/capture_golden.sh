#!/usr/bin/env bash
# Snapshot current pipeline outputs as golden references for the
# regression contract (post-audit: default + --global must reproduce
# these). Run after a deliberate, vetted pipeline change to update
# the reference.
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
        cp "$f" "$DEST/$(basename ${f%.txt})_reference.tsv"
        echo "captured $f"
    else
        echo "missing $f (skipped)"
    fi
done

# Also record the commit & timestamp so the regression test can reference it
git rev-parse HEAD > "$DEST/_captured_at_commit.txt"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DEST/_captured_at_time.txt"
