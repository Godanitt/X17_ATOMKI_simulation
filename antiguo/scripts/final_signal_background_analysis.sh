#!/usr/bin/env bash
# final_signal_background_analysis_fixed.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

SIG_RECO="${1:-analysis_hits_detector_effects.root}"
BKG_RECO="${2:-background_analysis_hits_detector_effects.root}"
OUT_ROOT="${3:-final_signal_background.root}"
SIG_SCALE="${4:-1.0}"
BKG_SCALE="${5:-1.0}"

MACRO="analysis/final_signal_background_analysis_fixed.C"

if ! command -v root >/dev/null 2>&1; then
    echo "[ERROR] ROOT is not available in PATH."
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

root -l -q "${MACRO}(\"${SIG_RECO}\",\"${BKG_RECO}\",\"${OUT_ROOT}\",${SIG_SCALE},${BKG_SCALE})"
