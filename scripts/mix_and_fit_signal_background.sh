#!/usr/bin/env bash
# mix_and_fit_signal_background.sh
#
# Combine signal and background detector-effects samples into one pseudo-data sample
# and estimate signal/background yields with a template fit.
#
# Default usage from project root:
#
#   bash scripts/mix_and_fit_signal_background.sh
#
# Full usage:
#
#   bash scripts/mix_and_fit_signal_background.sh \
#     analysis_hits_detector_effects.root \
#     background_analysis_hits_detector_effects.root \
#     signal_background_fit.root \
#     plots_signal_background_fit \
#     -1 \
#     -1 \
#     120 \
#     180 \
#     12345
#
# Arguments:
#   1 signal detector-effects ROOT file
#   2 background detector-effects ROOT file
#   3 output ROOT file
#   4 output plot directory
#   5 number of signal events to inject; -1 = all signal entries
#   6 number of background events to inject; -1 = all background entries
#   7 theta signal-region lower edge in deg
#   8 theta signal-region upper edge in deg
#   9 random seed

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

SIG_RECO="${1:-analysis_hits_detector_effects.root}"
BKG_RECO="${2:-background_analysis_hits_detector_effects.root}"
OUT_ROOT="${3:-signal_background_fit.root}"
OUT_DIR="${4:-plots_signal_background_fit}"
N_SIG="${5:--1}"
N_BKG="${6:--1}"
THETA_MIN="${7:-120}"
THETA_MAX="${8:-180}"
SEED="${9:-12345}"

MACRO="analysis/mix_and_fit_signal_background_fixed_v2.C"

if ! command -v root >/dev/null 2>&1; then
    echo "[ERROR] ROOT is not available in PATH. Source your ROOT environment first."
    exit 1
fi

if [[ ! -f "${MACRO}" ]]; then
    echo "[ERROR] Missing macro: ${MACRO}"
    exit 1
fi

if [[ ! -f "${SIG_RECO}" ]]; then
    echo "[ERROR] Missing signal input file: ${SIG_RECO}"
    exit 1
fi

if [[ ! -f "${BKG_RECO}" ]]; then
    echo "[ERROR] Missing background input file: ${BKG_RECO}"
    exit 1
fi

mkdir -p "${OUT_DIR}"

echo "=== Mixing and fitting signal/background ==="
echo "Signal detector-effects file     : ${SIG_RECO}"
echo "Background detector-effects file : ${BKG_RECO}"
echo "Output ROOT file                 : ${OUT_ROOT}"
echo "Output plot directory            : ${OUT_DIR}"
echo "Injected signal events           : ${N_SIG} (-1 = all)"
echo "Injected background events       : ${N_BKG} (-1 = all)"
echo "Theta signal region              : [${THETA_MIN}, ${THETA_MAX}] deg"
echo "Seed                             : ${SEED}"
echo

root -l -q "${MACRO}(\"${SIG_RECO}\",\"${BKG_RECO}\",\"${OUT_ROOT}\",\"${OUT_DIR}\",${N_SIG},${N_BKG},${THETA_MIN},${THETA_MAX},${SEED})"

echo
echo "=== Done ==="
echo "Output ROOT file : ${OUT_ROOT}"
echo "Plots directory  : ${OUT_DIR}"
