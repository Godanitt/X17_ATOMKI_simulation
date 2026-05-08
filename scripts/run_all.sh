#!/usr/bin/env bash
# run_all.sh
#
# Pipeline completo:
#   1. señal
#   2. background
#   3. análisis de hits
#   4. efectos de detector
#   5. análisis signal/background
#   6. fit por templates
#   7. fit fenomenológico IPC-like + exceso
#
# Uso:
#   bash scripts/run_all.sh
#
# Opcional:
#   bash scripts/run_all.sh 100000 0
#
# Argumentos:
#   $1 = número de eventos por muestra
#   $2 = volumeID para análisis de hits
#
# volumeID:
#   0  = SiliconStripLV
#   1  = ScintillatorLV
#  -1  = cualquier volumen

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

N_EVENTS="${1:-10000}"
VOLUME_ID="${2:-0}"

# Detector effects
EFFICIENCY="${EFFICIENCY:-0.90}"
SIGMA_THETA="${SIGMA_THETA:-2.0}"
SIGMA_PHI="${SIGMA_PHI:-2.0}"
ENERGY_RESOLUTION="${ENERGY_RESOLUTION:-0.05}"
ENERGY_THRESHOLD="${ENERGY_THRESHOLD:-1.0}"
DETECTOR_SEED="${DETECTOR_SEED:-12345}"

# Signal/background scales
SIGNAL_SCALE="${SIGNAL_SCALE:-1.0}"
BACKGROUND_SCALE="${BACKGROUND_SCALE:-1.0}"

# Template fit settings
N_SIGNAL_MIX="${N_SIGNAL_MIX:--1}"
N_BACKGROUND_MIX="${N_BACKGROUND_MIX:--1}"
THETA_REGION_MIN="${THETA_REGION_MIN:-120}"
THETA_REGION_MAX="${THETA_REGION_MAX:-180}"
FIT_SEED="${FIT_SEED:-12345}"

# IPC-like fit initial guesses
THETA_MU="${THETA_MU:-140}"
THETA_SIGMA="${THETA_SIGMA:-12}"
MASS_MU="${MASS_MU:-17}"
MASS_SIGMA="${MASS_SIGMA:-1}"
ENERGY_MU="${ENERGY_MU:-17.5}"
ENERGY_SIGMA="${ENERGY_SIGMA:-0.8}"

echo
echo "============================================================"
echo " X17 full pipeline"
echo "============================================================"
echo "Project dir              : ${PROJECT_DIR}"
echo "Events per sample        : ${N_EVENTS}"
echo "Hit volumeID             : ${VOLUME_ID}"
echo "Detector efficiency      : ${EFFICIENCY}"
echo "Angular smearing theta   : ${SIGMA_THETA} deg"
echo "Angular smearing phi     : ${SIGMA_PHI} deg"
echo "Energy resolution        : ${ENERGY_RESOLUTION}"
echo "Energy threshold         : ${ENERGY_THRESHOLD} MeV"
echo "============================================================"
echo

run_step() {
    local title="$1"
    shift

    echo
    echo "------------------------------------------------------------"
    echo "${title}"
    echo "------------------------------------------------------------"
    echo "+ $*"
    "$@"
}

# Optional cleanup
if [[ "${CLEAN:-0}" == "1" ]]; then
    echo "Cleaning previous outputs..."
    rm -f sampled.root
    rm -f background_sampled.root
    rm -f analysis_hits.root
    rm -f background_analysis_hits.root
    rm -f analysis_hits_detector_effects.root
    rm -f background_analysis_hits_detector_effects.root
    rm -f final_signal_background.root
    rm -f signal_background_fit.root
    rm -f ipc_excess_fit.root
    rm -rf plots
    rm -rf plots_final_signal_background
    rm -rf plots_signal_background_fit
    rm -rf plots_ipc_excess_fit
fi

# Optional detector visualization.
# Not run by default because build.sh may open the Geant4 visualizer and block the pipeline.
if [[ "${RUN_VIS:-0}" == "1" ]]; then
    run_step "Building and launching detector visualization" \
        bash scripts/build.sh
fi

run_step "1/9 - Generate signal sample" \
    bash scripts/run_generated.sh "${N_EVENTS}" sampled.root

run_step "2/9 - Generate background sample" \
    bash scripts/run_background.sh "${N_EVENTS}" background_sampled.root

run_step "3/9 - Analyze signal hits" \
    bash scripts/analyze_hits.sh sampled.root analysis_hits.root "${VOLUME_ID}"

run_step "4/9 - Analyze background hits" \
    bash scripts/analyze_background_hits.sh background_sampled.root background_analysis_hits.root "${VOLUME_ID}"

run_step "5/9 - Apply detector effects to signal" \
    bash scripts/apply_detector_effects.sh \
        analysis_hits.root \
        analysis_hits_detector_effects.root \
        "${EFFICIENCY}" \
        "${SIGMA_THETA}" \
        "${SIGMA_PHI}" \
        "${ENERGY_RESOLUTION}" \
        "${ENERGY_THRESHOLD}" \
        "${DETECTOR_SEED}"

run_step "6/9 - Apply detector effects to background" \
    bash scripts/apply_background_detector_effects.sh \
        background_analysis_hits.root \
        background_analysis_hits_detector_effects.root \
        "${EFFICIENCY}" \
        "${SIGMA_THETA}" \
        "${SIGMA_PHI}" \
        "${ENERGY_RESOLUTION}" \
        "${ENERGY_THRESHOLD}" \
        "${DETECTOR_SEED}"

run_step "7/9 - Final signal/background comparison" \
    bash scripts/final_signal_background_analysis.sh \
        analysis_hits_detector_effects.root \
        background_analysis_hits_detector_effects.root \
        final_signal_background.root \
        "${SIGNAL_SCALE}" \
        "${BACKGROUND_SCALE}"

run_step "8/9 - Template fit signal/background" \
    bash scripts/mix_and_fit_signal_background.sh \
        analysis_hits_detector_effects.root \
        background_analysis_hits_detector_effects.root \
        signal_background_fit.root \
        plots_signal_background_fit \
        "${N_SIGNAL_MIX}" \
        "${N_BACKGROUND_MIX}" \
        "${THETA_REGION_MIN}" \
        "${THETA_REGION_MAX}" \
        "${FIT_SEED}"

run_step "9/9 - IPC-like smooth background + X17 excess fit" \
    bash scripts/fit_ipc_excess.sh \
        analysis_hits_detector_effects.root \
        background_analysis_hits_detector_effects.root \
        ipc_excess_fit.root \
        plots_ipc_excess_fit \
        "${N_SIGNAL_MIX}" \
        "${N_BACKGROUND_MIX}" \
        "${FIT_SEED}" \
        "${THETA_MU}" \
        "${THETA_SIGMA}" \
        "${MASS_MU}" \
        "${MASS_SIGMA}" \
        "${ENERGY_MU}" \
        "${ENERGY_SIGMA}"

echo
echo "============================================================"
echo " Pipeline finished"
echo "============================================================"
echo
echo "Main outputs:"
echo "  sampled.root"
echo "  background_sampled.root"
echo "  analysis_hits.root"
echo "  background_analysis_hits.root"
echo "  analysis_hits_detector_effects.root"
echo "  background_analysis_hits_detector_effects.root"
echo "  final_signal_background.root"
echo "  signal_background_fit.root"
echo "  ipc_excess_fit.root"
echo
echo "Plot folders:"
echo "  plots_final_signal_background/"
echo "  plots_signal_background_fit/"
echo "  plots_ipc_excess_fit/"
echo