#!/usr/bin/env bash
# run_all.sh
#
# Full X17 pipeline with predefined geometries and per-geometry output folders.
#
# Usage:
#   bash scripts/run_all.sh all 10000 0
#   bash scripts/run_all.sh current 10000 0
#   bash scripts/run_all.sh 2pi 10000 0
#   bash scripts/run_all.sh 4pi 10000 0
#   bash scripts/run_all.sh padplane 10000 0
#
# Backwards-compatible usage:
#   bash scripts/run_all.sh 10000 0       # runs all geometries
#
# Arguments:
#   $1 = geometry name or "all". If numeric, interpreted as N_EVENTS.
#   $2 = number of events per sample
#   $3 = volumeID for hit analysis
#
# volumeID:
#   0  = SiliconStripLV
#   1  = ScintillatorLV
#  -1  = any detector volume

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

GEOMETRIES=("current" "2pi" "4pi" "padplane")

is_integer() {
    [[ "${1:-}" =~ ^-?[0-9]+$ ]]
}

normalize_geometry() {
    local g="${1,,}"
    case "${g}" in
        current|baseline|atomki) echo "current" ;;
        2pi|coverage2pi|two_pi) echo "2pi" ;;
        4pi|coverage4pi|four_pi) echo "4pi" ;;
        padplane|pad_plane|pads|frontpads) echo "padplane" ;;
        all|ALL) echo "all" ;;
        *) echo "" ;;
    esac
}

if [[ $# -gt 0 ]] && is_integer "$1"; then
    GEOMETRY_REQUEST="all"
    N_EVENTS="$1"
    VOLUME_ID="${2:-0}"
else
    GEOMETRY_REQUEST="${1:-all}"
    N_EVENTS="${2:-10000}"
    VOLUME_ID="${3:-0}"
fi

GEOMETRY_REQUEST="$(normalize_geometry "${GEOMETRY_REQUEST}")"
if [[ -z "${GEOMETRY_REQUEST}" ]]; then
    echo "ERROR: unknown geometry. Use one of: all, ${GEOMETRIES[*]}" >&2
    exit 2
fi

if [[ "${GEOMETRY_REQUEST}" == "all" ]]; then
    echo "=== Running ALL predefined geometries ==="
    printf 'Geometries:'
    for GEO in "${GEOMETRIES[@]}"; do
        printf ' %s' "${GEO}"
    done
    printf '\n'
    echo "Events per sample: ${N_EVENTS}"
    echo "Hit volumeID     : ${VOLUME_ID}"
    echo "Results root     : ${PROJECT_DIR}/results/<geometry>/"

    FAILED_GEOMETRIES=()
    for GEO in "${GEOMETRIES[@]}"; do
        echo
        echo "############################################################"
        echo "# Geometry: ${GEO}"
        echo "############################################################"
        if bash "${SCRIPT_DIR}/run_all.sh" "${GEO}" "${N_EVENTS}" "${VOLUME_ID}"; then
            echo "# Geometry finished successfully: ${GEO}"
        else
            status=$?
            echo "# ERROR: geometry failed: ${GEO} (exit status ${status})" >&2
            FAILED_GEOMETRIES+=("${GEO}:${status}")
        fi
    done

    echo
    echo "============================================================"
    echo " Geometry scan summary"
    echo "============================================================"
    echo "Requested geometries: ${GEOMETRIES[*]}"
    echo "Results directory   : ${PROJECT_DIR}/results/<geometry>/"
    if (( ${#FAILED_GEOMETRIES[@]} > 0 )); then
        echo "Failed geometries   : ${FAILED_GEOMETRIES[*]}" >&2
        echo "Completed geometries are still available under results/<geometry>/" >&2
        exit 1
    fi
    echo "All geometry pipelines finished successfully."
    exit 0
fi

GEOMETRY="${GEOMETRY_REQUEST}"
OUT_DIR="${PROJECT_DIR}/results/${GEOMETRY}"
LOG_DIR="${OUT_DIR}/logs"
INPUT_TABLE="${PROJECT_DIR}/data/data_pair_creation.txt"

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

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

if [[ "${CLEAN:-0}" == "1" ]]; then
    echo "Cleaning previous outputs for geometry ${GEOMETRY}..."
    rm -rf "${OUT_DIR}"/*
    mkdir -p "${OUT_DIR}" "${LOG_DIR}"
fi

cat > "${OUT_DIR}/config.txt" <<CONFIG_EOF
geometry = ${GEOMETRY}
n_events = ${N_EVENTS}
volume_id = ${VOLUME_ID}
efficiency = ${EFFICIENCY}
sigma_theta_deg = ${SIGMA_THETA}
sigma_phi_deg = ${SIGMA_PHI}
energy_resolution = ${ENERGY_RESOLUTION}
energy_threshold_MeV = ${ENERGY_THRESHOLD}
detector_seed = ${DETECTOR_SEED}
signal_scale = ${SIGNAL_SCALE}
background_scale = ${BACKGROUND_SCALE}
n_signal_mix = ${N_SIGNAL_MIX}
n_background_mix = ${N_BACKGROUND_MIX}
theta_region_min_deg = ${THETA_REGION_MIN}
theta_region_max_deg = ${THETA_REGION_MAX}
fit_seed = ${FIT_SEED}
CONFIG_EOF

echo
echo "============================================================"
echo " X17 full pipeline"
echo "============================================================"
echo "Project dir              : ${PROJECT_DIR}"
echo "Geometry                 : ${GEOMETRY}"
echo "Output dir               : ${OUT_DIR}"
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

# Build once before running this geometry. The executable itself supports all geometries.
run_step "0/9 - Build executable" \
    bash scripts/build.sh "${GEOMETRY}" --build-only

run_step "1/9 - Generate signal sample" \
    bash scripts/run_generated.sh \
        "${N_EVENTS}" \
        "${OUT_DIR}/sampled.root" \
        "${INPUT_TABLE}" \
        "${GEOMETRY}" \
        "${LOG_DIR}"

run_step "2/9 - Generate background sample" \
    bash scripts/run_background.sh \
        "${N_EVENTS}" \
        "${OUT_DIR}/background_sampled.root" \
        "${GEOMETRY}" \
        "${LOG_DIR}"

run_step "3/9 - Analyze signal hits" \
    bash scripts/analyze_hits.sh \
        "${OUT_DIR}/sampled.root" \
        "${OUT_DIR}/analysis_hits.root" \
        "${VOLUME_ID}" \
        "${GEOMETRY}"

run_step "4/9 - Analyze background hits" \
    bash scripts/analyze_background_hits.sh \
        "${OUT_DIR}/background_sampled.root" \
        "${OUT_DIR}/background_analysis_hits.root" \
        "${VOLUME_ID}" \
        "${GEOMETRY}"

run_step "5/9 - Apply detector effects to signal" \
    bash scripts/apply_detector_effects.sh \
        "${OUT_DIR}/analysis_hits.root" \
        "${OUT_DIR}/analysis_hits_detector_effects.root" \
        "${EFFICIENCY}" \
        "${SIGMA_THETA}" \
        "${SIGMA_PHI}" \
        "${ENERGY_RESOLUTION}" \
        "${ENERGY_THRESHOLD}" \
        "${DETECTOR_SEED}"

run_step "6/9 - Apply detector effects to background" \
    bash scripts/apply_background_detector_effects.sh \
        "${OUT_DIR}/background_analysis_hits.root" \
        "${OUT_DIR}/background_analysis_hits_detector_effects.root" \
        "${EFFICIENCY}" \
        "${SIGMA_THETA}" \
        "${SIGMA_PHI}" \
        "${ENERGY_RESOLUTION}" \
        "${ENERGY_THRESHOLD}" \
        "${DETECTOR_SEED}"

run_step "7/9 - Final signal/background comparison" \
    bash scripts/final_signal_background_analysis.sh \
        "${OUT_DIR}/analysis_hits_detector_effects.root" \
        "${OUT_DIR}/background_analysis_hits_detector_effects.root" \
        "${OUT_DIR}/final_signal_background.root" \
        "${SIGNAL_SCALE}" \
        "${BACKGROUND_SCALE}"

run_step "8/9 - Template fit signal/background" \
    bash scripts/mix_and_fit_signal_background.sh \
        "${OUT_DIR}/analysis_hits_detector_effects.root" \
        "${OUT_DIR}/background_analysis_hits_detector_effects.root" \
        "${OUT_DIR}/signal_background_fit.root" \
        "${OUT_DIR}/plots_signal_background_fit" \
        "${N_SIGNAL_MIX}" \
        "${N_BACKGROUND_MIX}" \
        "${THETA_REGION_MIN}" \
        "${THETA_REGION_MAX}" \
        "${FIT_SEED}"

run_step "9/9 - IPC-like smooth background + X17 excess fit" \
    bash scripts/fit_ipc_excess.sh \
        "${OUT_DIR}/analysis_hits_detector_effects.root" \
        "${OUT_DIR}/background_analysis_hits_detector_effects.root" \
        "${OUT_DIR}/ipc_excess_fit.root" \
        "${OUT_DIR}/plots_ipc_excess_fit" \
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
echo " Pipeline finished for geometry: ${GEOMETRY}"
echo "============================================================"
echo "Main outputs:"
echo "  ${OUT_DIR}/sampled.root"
echo "  ${OUT_DIR}/background_sampled.root"
echo "  ${OUT_DIR}/analysis_hits.root"
echo "  ${OUT_DIR}/background_analysis_hits.root"
echo "  ${OUT_DIR}/analysis_hits_detector_effects.root"
echo "  ${OUT_DIR}/background_analysis_hits_detector_effects.root"
echo "  ${OUT_DIR}/final_signal_background.root"
echo "  ${OUT_DIR}/signal_background_fit.root"
echo "  ${OUT_DIR}/ipc_excess_fit.root"
echo
echo "Plot folders:"
echo "  ${OUT_DIR}/plots_signal_background_fit/"
echo "  ${OUT_DIR}/plots_ipc_excess_fit/"
echo
