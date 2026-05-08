#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
INPUT="${1:-${PROJECT_DIR}/background_sampled.root}"
OUTPUT="${2:-${PROJECT_DIR}/background_analysis_hits.root}"
VOLUME_ID="${3:-0}"
exec "${SCRIPT_DIR}/analyze_hits.sh" "${INPUT}" "${OUTPUT}" "${VOLUME_ID}"
