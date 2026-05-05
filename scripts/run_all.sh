#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${PROJECT_DIR}/build"
LOG_DIR="${PROJECT_DIR}/logs"
PLOTS_DIR="${PROJECT_DIR}/plots"

# Safe defaults for a laptop. Increase later only after checking temperatures.
JOBS=4
THREADS=2
EVENTS=10000

cd "${PROJECT_DIR}"
mkdir -p "${LOG_DIR}" "${PLOTS_DIR}"

rm -f x17_output.root x17_analysis.root
rm -f plots/*.pdf

echo "==> Project: ${PROJECT_DIR}"
echo "==> Cleaning old build/"
rm -rf "${BUILD_DIR}"

command -v cmake >/dev/null 2>&1 || { echo "ERROR: cmake not found. Source your Geant4/CMake environment first."; exit 1; }

echo "==> Configuring Geant4 project"
cmake -S "${PROJECT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release

echo "==> Building with ${JOBS} jobs"
cmake --build "${BUILD_DIR}" --parallel "${JOBS}"

[[ -x "build/x17sim" ]] || { echo "ERROR: build/x17sim not found"; exit 1; }
[[ -f "data/data_pair_creation.txt" ]] || { echo "ERROR: data/data_pair_creation.txt not found"; exit 1; }
[[ -f "scripts/run_analysis.sh" ]] || { echo "ERROR: scripts/run_analysis.sh not found"; exit 1; }

echo "==> Creating batch macro with ${EVENTS} events and ${THREADS} threads"
cat > "logs/run_default.mac" <<MACRO
/run/numberOfThreads ${THREADS}
/run/initialize
/run/beamOn ${EVENTS}
MACRO

echo "==> Running Geant4 simulation"
./build/x17sim logs/run_default.mac 2>&1 | tee logs/geant4_run.log

[[ -f "x17_output.root" ]] || { echo "ERROR: x17_output.root was not created"; exit 1; }

echo "==> Running analysis"
bash scripts/run_analysis.sh

if [[ -f "scripts/make_root_plots.sh" ]]; then
    echo "==> Making input-table ROOT plots"
    bash scripts/make_root_plots.sh || true
fi

echo "==> Finished full chain"
echo

echo "ROOT files:"
ls -lh x17_output.root x17_analysis.root 2>/dev/null || true
echo

echo "PDF plots:"
ls -lh plots/*.pdf 2>/dev/null || true
