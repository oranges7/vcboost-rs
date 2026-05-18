#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DIST_DIR="${PROJECT_DIR}/dist"

echo "=== Quick Build: x86_64 only ==="
mkdir -p "${DIST_DIR}"

if command -v maturin &> /dev/null; then
    echo "Using local maturin..."
    maturin build --release --out "${DIST_DIR}"
else
    echo "Using Docker + maturin..."
    docker run --rm -v "${PROJECT_DIR}":/io --workdir /io ghcr.io/pyo3/maturin:latest build --release --out /io/dist
fi

echo "=== Build Complete ==="
ls -la "${DIST_DIR}/"*.whl 2>/dev/null
echo "Install with: pip install ${DIST_DIR}/vcboostrs-*.whl"
