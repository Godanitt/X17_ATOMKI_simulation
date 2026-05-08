#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
INPUT="${1:-${PROJECT_DIR}/background_analysis_hits.root}"
OUTPUT="${2:-${PROJECT_DIR}/background_analysis_hits_detector_effects.root}"
EFFICIENCY="${3:-0.90}"
SIGMA_THETA_DEG="${4:-2.0}"
SIGMA_PHI_DEG="${5:-2.0}"
RELATIVE_ENERGY_RESOLUTION="${6:-0.05}"
ENERGY_THRESHOLD_MEV="${7:-1.0}"
SEED="${8:-23456}"
exec "${SCRIPT_DIR}/apply_detector_effects.sh" "${INPUT}" "${OUTPUT}" "${EFFICIENCY}" "${SIGMA_THETA_DEG}" "${SIGMA_PHI_DEG}" "${RELATIVE_ENERGY_RESOLUTION}" "${ENERGY_THRESHOLD_MEV}" "${SEED}"
