#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUT="${1:-${PROJECT_DIR}/analysis_hits.root}"
OUTPUT="${2:-${PROJECT_DIR}/analysis_hits_detector_effects.root}"
EFFICIENCY="${3:-0.90}"
SIGMA_THETA_DEG="${4:-2.0}"
SIGMA_PHI_DEG="${5:-2.0}"
RELATIVE_ENERGY_RESOLUTION="${6:-0.05}"
ENERGY_THRESHOLD_MEV="${7:-1.0}"
SEED="${8:-12345}"

if ! command -v root >/dev/null 2>&1; then
  echo "ERROR: ROOT executable not found in PATH." >&2
  exit 2
fi

if [[ ! -f "${INPUT}" ]]; then
  echo "ERROR: input analysis ROOT file not found: ${INPUT}" >&2
  echo "Run scripts/analyze_hits.sh first." >&2
  exit 2
fi

cd "${PROJECT_DIR}"
root -l -b -q "analysis/apply_detector_effects.C(\"${INPUT}\",\"${OUTPUT}\",${EFFICIENCY},${SIGMA_THETA_DEG},${SIGMA_PHI_DEG},${RELATIVE_ENERGY_RESOLUTION},${ENERGY_THRESHOLD_MEV},${SEED})"
