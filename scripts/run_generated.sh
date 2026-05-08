#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${PROJECT_DIR}/build"
LOG_DIR="${PROJECT_DIR}/logs"

NEVENTS="${1:-10000}"
OUTPUT="${2:-${PROJECT_DIR}/sampled.root}"
INPUT="${3:-${PROJECT_DIR}/data/data_pair_creation.txt}"
MACRO="${LOG_DIR}/run_sampled.mac"

cd "${PROJECT_DIR}"
mkdir -p "${BUILD_DIR}" "${LOG_DIR}"

if [[ ! -f "${INPUT}" ]]; then
  echo "ERROR: input data table not found: ${INPUT}" >&2
  echo "Pass it as the 3rd argument, e.g.: bash scripts/run_generated.sh 10000 sampled.root data/data_pair_creation.txt" >&2
  exit 2
fi

# Always run CMake/build. This avoids accidentally using a stale executable
# from an older version that wrote -999 in the generated tree.
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
echo "Using input data table: ${INPUT}"
echo "First data row [thetaEE thetaEe energyEe thetaEp energyEp]: $(grep -v '^#' "${INPUT}" | head -n 1)"
echo "Output ROOT file: ${OUTPUT}"
"${BUILD_DIR}/x17sim" -m "${MACRO}" -i "${INPUT}" -o "${OUTPUT}" --mode signal | tee "${LOG_DIR}/run_generated.log"

echo "Wrote ${OUTPUT} with exactly two trees: generated, hits"
