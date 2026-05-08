#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUT="${1:-${PROJECT_DIR}/sampled.root}"
OUTPUT="${2:-${PROJECT_DIR}/analysis_hits.root}"
VOLUME_ID="${3:-0}"

if ! command -v root >/dev/null 2>&1; then
  echo "ERROR: ROOT executable not found in PATH." >&2
  exit 2
fi

if [[ ! -f "${INPUT}" ]]; then
  echo "ERROR: input ROOT file not found: ${INPUT}" >&2
  exit 2
fi

cd "${PROJECT_DIR}"
root -l -b -q "analysis/analyze_hits.C(\"${INPUT}\",\"${OUTPUT}\",${VOLUME_ID})"
