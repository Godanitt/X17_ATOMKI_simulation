#!/usr/bin/env bash
set -Eeuo pipefail

# Detector-effects ROOT study. No options.
# Usage from anywhere:
#   bash scripts/run_detector_effects.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${PROJECT_DIR}"

echo "==> Running detector-effects study"

command -v root >/dev/null 2>&1 || { echo "ERROR: ROOT not found. Source ROOT first."; exit 1; }
[[ -f "x17_output.root" ]] || { echo "ERROR: x17_output.root not found. Run: bash scripts/run_all.sh"; exit 1; }
[[ -f "scripts/study_detector_effects.C" ]] || { echo "ERROR: scripts/study_detector_effects.C not found"; exit 1; }
[[ -f "scripts/X17Style.C" ]] || { echo "ERROR: scripts/X17Style.C not found"; exit 1; }

mkdir -p plots

root -l -b -q -I scripts <<'ROOTEOF'
.L scripts/study_detector_effects.C
study_detector_effects("x17_output.root", "plots");
.q
ROOTEOF

[[ -f "x17_detector_effects.root" ]] || { echo "ERROR: x17_detector_effects.root was not created"; exit 1; }
[[ -f "plots/detector_acceptance_vs_theta.pdf" ]] || { echo "ERROR: plots/detector_acceptance_vs_theta.pdf was not created"; exit 1; }
[[ -f "plots/detector_thetaee_distortion.pdf" ]] || { echo "ERROR: plots/detector_thetaee_distortion.pdf was not created"; exit 1; }
[[ -f "plots/detector_thetaee_resolution.pdf" ]] || { echo "ERROR: plots/detector_thetaee_resolution.pdf was not created"; exit 1; }
[[ -f "plots/detector_threshold_scan.pdf" ]] || { echo "ERROR: plots/detector_threshold_scan.pdf was not created"; exit 1; }

echo "==> Detector-effects study done"
echo "  x17_detector_effects.root"
echo "  plots/detector_acceptance_vs_theta.pdf"
echo "  plots/detector_thetaee_distortion.pdf"
echo "  plots/detector_thetaee_resolution.pdf"
echo "  plots/detector_threshold_scan.pdf"
