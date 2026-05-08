#!/usr/bin/env bash
# fit_ipc_excess.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

SIG_RECO="${1:-analysis_hits_detector_effects.root}"
BKG_RECO="${2:-background_analysis_hits_detector_effects.root}"
OUT_ROOT="${3:-ipc_excess_fit.root}"
OUT_DIR="${4:-plots_ipc_excess_fit}"

N_SIG="${5:--1}"
N_BKG="${6:--1}"
SEED="${7:-12345}"

THETA_MU="${8:-140}"
THETA_SIGMA="${9:-12}"

MASS_MU="${10:-17}"
MASS_SIGMA="${11:-1}"

ENERGY_MU="${12:-17.5}"
ENERGY_SIGMA="${13:-0.8}"

MACRO="analysis/fit_ipc_excess.C"

if ! command -v root >/dev/null 2>&1; then
    echo "[ERROR] ROOT is not available in PATH."
    exit 1
fi

if [[ ! -f "${MACRO}" ]]; then
    echo "[ERROR] Missing macro: ${MACRO}"
    exit 1
fi

if [[ ! -f "${SIG_RECO}" ]]; then
    echo "[ERROR] Missing signal file: ${SIG_RECO}"
    exit 1
fi

if [[ ! -f "${BKG_RECO}" ]]; then
    echo "[ERROR] Missing background file: ${BKG_RECO}"
    exit 1
fi

mkdir -p "${OUT_DIR}"

echo "=== Running IPC-like background + X17 excess fit ==="
echo "Signal reco file       : ${SIG_RECO}"
echo "Background reco file   : ${BKG_RECO}"
echo "Output ROOT file       : ${OUT_ROOT}"
echo "Output plot directory  : ${OUT_DIR}"
echo "Injected signal events : ${N_SIG} (-1 = all)"
echo "Injected bkg events    : ${N_BKG} (-1 = all)"
echo "Seed                   : ${SEED}"
echo "Theta init             : mu=${THETA_MU}, sigma=${THETA_SIGMA}"
echo "Mass init              : mu=${MASS_MU}, sigma=${MASS_SIGMA}"
echo "Energy init            : mu=${ENERGY_MU}, sigma=${ENERGY_SIGMA}"
echo

root -l -q "${MACRO}(\"${SIG_RECO}\",\"${BKG_RECO}\",\"${OUT_ROOT}\",\"${OUT_DIR}\",${N_SIG},${N_BKG},${SEED},${THETA_MU},${THETA_SIGMA},${MASS_MU},${MASS_SIGMA},${ENERGY_MU},${ENERGY_SIGMA})"

echo
echo "=== Done ==="
echo "Output ROOT file : ${OUT_ROOT}"
echo "Plots directory  : ${OUT_DIR}"
