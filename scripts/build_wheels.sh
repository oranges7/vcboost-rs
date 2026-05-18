#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DIST_DIR="${PROJECT_DIR}/dist"

echo "=== vcboostrs Wheel Builder ==="
mkdir -p "${DIST_DIR}"

build_wheel() {
    local python_version="$1"
    local target="$2"
    echo ">>> Building wheel for Python ${python_version} (${target})..."
    docker run --rm -v "${PROJECT_DIR}":/io --workdir /io ghcr.io/pyo3/maturin:latest build --release --interpreter python${python_version} --target ${target} --out /io/dist
    echo ">>> Done: Python ${python_version} (${target})"
}

echo "[1/2] Building x86_64 wheels..."
for py in 3.9 3.10 3.11 3.12; do build_wheel "${py}" "x86_64"; done

echo "[2/2] Building aarch64 wheels (cross-compile)..."
for py in 3.9 3.10 3.11 3.12; do build_wheel "${py}" "aarch64"; done

echo "=== Build Complete ==="
ls -la "${DIST_DIR}/"*.whl 2>/dev/null || echo "No wheels found"
