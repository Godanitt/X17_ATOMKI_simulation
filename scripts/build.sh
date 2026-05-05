#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${PROJECT_DIR}/build"
JOBS=4
BUILD_TYPE="Release"

cd "${PROJECT_DIR}"

command -v cmake >/dev/null 2>&1 || { echo "ERROR: cmake not found. Source your Geant4/CMake environment first."; exit 1; }

rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"

echo "==> Project: ${PROJECT_DIR}"
echo "==> Cleaning old build/"
echo "==> Configuring Geant4 project"
cmake -S "${PROJECT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"

echo "==> Building with ${JOBS} jobs"
cmake --build "${BUILD_DIR}" --parallel "${JOBS}"

[[ -x "${BUILD_DIR}/x17sim" ]] || { echo "ERROR: executable not found: ${BUILD_DIR}/x17sim"; exit 1; }

echo "==> Build finished: ${BUILD_DIR}/x17sim"
