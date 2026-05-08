#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${PROJECT_DIR}/build"
LOG_DIR="${PROJECT_DIR}/logs"

NEVENTS="${1:-10000}"
OUTPUT="${2:-${PROJECT_DIR}/background_sampled.root}"
MACRO="${LOG_DIR}/run_background.mac"

cd "${PROJECT_DIR}"
mkdir -p "${BUILD_DIR}" "${LOG_DIR}"

cmake -S "${PROJECT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_DIR}" --parallel "$(nproc)"

cat > "${MACRO}" <<MACRO_EOF
/control/verbose 0
/run/verbose 1
/event/verbose 0
/tracking/verbose 0
/run/numberOfThreads 1
/run/initialize
/run/beamOn ${NEVENTS}
MACRO_EOF

rm -f "${OUTPUT}"
echo "Background mode: smooth IPC-like toy template"
echo "Output ROOT file: ${OUTPUT}"
"${BUILD_DIR}/x17sim" -m "${MACRO}" -o "${OUTPUT}" --mode background | tee "${LOG_DIR}/run_background.log"

echo "Wrote ${OUTPUT} with exactly two trees: generated, hits"
