#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${PROJECT_DIR}/build"

NEVENTS="${1:-10000}"
OUTPUT="${2:-${PROJECT_DIR}/background_sampled.root}"
GEOMETRY="${3:-current}"
LOG_DIR="${4:-${PROJECT_DIR}/logs/${GEOMETRY}}"
MACRO="${LOG_DIR}/run_background.mac"

cd "${PROJECT_DIR}"

if [[ -f "${BUILD_DIR}/CMakeCache.txt" ]] && ! grep -q "CMAKE_HOME_DIRECTORY:INTERNAL=${PROJECT_DIR}" "${BUILD_DIR}/CMakeCache.txt"; then
  echo "Removing stale build directory from a different source path: ${BUILD_DIR}"
  rm -rf "${BUILD_DIR}"
fi

mkdir -p "${BUILD_DIR}" "${LOG_DIR}" "$(dirname "${OUTPUT}")"

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
echo "Geometry: ${GEOMETRY}"
echo "Output ROOT file: ${OUTPUT}"
"${BUILD_DIR}/x17sim" -m "${MACRO}" -o "${OUTPUT}" --mode background --geometry "${GEOMETRY}" | tee "${LOG_DIR}/run_background.log"

echo "Wrote ${OUTPUT} with exactly two trees: generated, hits"
