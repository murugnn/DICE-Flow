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
