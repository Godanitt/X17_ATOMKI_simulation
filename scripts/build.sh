#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${PROJECT_DIR}/build"
INPUT="${PROJECT_DIR}/data/data_pair_creation.txt"

cd "${PROJECT_DIR}"
mkdir -p "${BUILD_DIR}"

if [[ ! -f "${INPUT}" ]]; then
  echo "ERROR: input data table not found: ${INPUT}" >&2
  exit 2
fi

cmake -S "${PROJECT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_DIR}" --parallel "$(nproc)"

# Launch detector visualization automatically for screenshots.
"${BUILD_DIR}/x17sim" --vis -i "${INPUT}" -o "${PROJECT_DIR}/sampled.root"
