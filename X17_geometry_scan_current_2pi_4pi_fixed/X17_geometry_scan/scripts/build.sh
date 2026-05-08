#!/usr/bin/env bash
set -Eeuo pipefail

# Build and, by default, open the Geant4 visualisation for one predefined geometry.
#
# Usage:
#   ./build.sh current
#   ./build.sh 2pi
#   ./build.sh 4pi
#   ./build.sh padplane
#
# Options:
#   --build-only / --no-vis   Build only, do not launch the viewer.
#   --vis                     Explicitly launch the viewer after building.
#   --clean                   Remove build/ before configuring.
#   --geometry GEOM           Alternative way to pass the geometry.
#
# Available geometries: current, 2pi, 4pi, padplane

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${PROJECT_DIR}/build"
INPUT="${PROJECT_DIR}/data/data_pair_creation.txt"

GEOMETRY="current"
RUN_VIS=1
CLEAN=0

normalize_geometry() {
  local g="${1,,}"
  case "${g}" in
    current|baseline|atomki) echo "current" ;;
    2pi|coverage2pi|two_pi) echo "2pi" ;;
    4pi|coverage4pi|four_pi) echo "4pi" ;;
    padplane|pad_plane|pads|frontpads) echo "padplane" ;;
    *) echo "" ;;
  esac
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --geometry|--geom)
      GEOMETRY="${2:-}"
      shift 2
      ;;
    --vis)
      RUN_VIS=1
      shift
      ;;
    --build-only|--no-vis)
      RUN_VIS=0
      shift
      ;;
    --clean)
      CLEAN=1
      shift
      ;;
    current|baseline|atomki|2pi|coverage2pi|two_pi|4pi|coverage4pi|four_pi|padplane|pad_plane|pads|frontpads)
      GEOMETRY="$1"
      shift
      ;;
    -h|--help)
      sed -n '1,28p' "$0"
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      echo "Use: ./build.sh [current|2pi|4pi|padplane] [--build-only] [--clean]" >&2
      exit 2
      ;;
  esac
done

GEOMETRY="$(normalize_geometry "${GEOMETRY}")"
if [[ -z "${GEOMETRY}" ]]; then
  echo "ERROR: unknown geometry. Use one of: current, 2pi, 4pi, padplane" >&2
  exit 2
fi

cd "${PROJECT_DIR}"

if [[ "${CLEAN}" == "1" ]]; then
  echo "Cleaning build directory: ${BUILD_DIR}"
  rm -rf "${BUILD_DIR}"
fi

if [[ -f "${BUILD_DIR}/CMakeCache.txt" ]] && ! grep -q "CMAKE_HOME_DIRECTORY:INTERNAL=${PROJECT_DIR}" "${BUILD_DIR}/CMakeCache.txt"; then
  echo "Removing stale build directory from a different source path: ${BUILD_DIR}"
  rm -rf "${BUILD_DIR}"
fi

mkdir -p "${BUILD_DIR}" "${PROJECT_DIR}/results/${GEOMETRY}" "${PROJECT_DIR}/results/${GEOMETRY}/logs"

if [[ ! -f "${INPUT}" ]]; then
  echo "ERROR: input data table not found: ${INPUT}" >&2
  exit 2
fi

echo "=== Building X17 simulation ==="
echo "Project dir : ${PROJECT_DIR}"
echo "Build dir   : ${BUILD_DIR}"
echo "Geometry    : ${GEOMETRY}"
echo "Viewer      : $([[ "${RUN_VIS}" == "1" ]] && echo yes || echo no)"
echo

cmake -S "${PROJECT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_DIR}" --parallel "$(nproc)"

echo
echo "Build complete: ${BUILD_DIR}/x17sim"

if [[ "${RUN_VIS}" == "1" ]]; then
  echo
  echo "=== Opening Geant4 viewer for geometry: ${GEOMETRY} ==="
  echo "Close the viewer/session to return to the terminal."
  echo
  "${BUILD_DIR}/x17sim" --vis --geometry "${GEOMETRY}" -i "${INPUT}" -o "${PROJECT_DIR}/results/${GEOMETRY}/preview.root"
else
  echo "Open the viewer with: ./build.sh ${GEOMETRY} --vis"
fi
