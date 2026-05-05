#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOG_DIR="${PROJECT_DIR}/logs"

cd "${PROJECT_DIR}"
mkdir -p plots "${LOG_DIR}"

command -v root >/dev/null 2>&1 || { echo "ERROR: ROOT not found. Source ROOT first."; exit 1; }
[[ -f "x17_output.root" ]] || { echo "ERROR: x17_output.root not found. Run: bash scripts/run_all.sh"; exit 1; }
[[ -f "scripts/analyze_X17.C" ]] || { echo "ERROR: scripts/analyze_X17.C not found"; exit 1; }
[[ -f "scripts/X17Style.C" ]] || { echo "ERROR: scripts/X17Style.C not found"; exit 1; }

rm -f x17_analysis.root
rm -f plots/*.pdf

cat > "${LOG_DIR}/run_analysis_driver.C" <<'ROOTEOF'
{
    gSystem->AddIncludePath("-Iscripts");

    Int_t error = 0;
    gROOT->ProcessLine(".L scripts/analyze_X17.C", &error);

    if (error != 0) {
        std::cerr << "ERROR: could not load scripts/analyze_X17.C" << std::endl;
        gSystem->Exit(10);
    }

    gROOT->ProcessLine("analyze_x17(\"x17_output.root\", \"plots\");", &error);

    if (error != 0) {
        std::cerr << "ERROR: analyze_x17 failed" << std::endl;
        gSystem->Exit(11);
    }
}
ROOTEOF

echo "==> Running offline analysis"
root -l -b -q "${LOG_DIR}/run_analysis_driver.C" | tee "${LOG_DIR}/analysis.log"

[[ -f "x17_analysis.root" ]] || { echo "ERROR: x17_analysis.root was not created"; exit 1; }

required=(
  plots/thetaee_components_generated.pdf
  plots/thetaee_components_detector.pdf
  plots/thetaee_generated_ideal_geometry_detector.pdf
  plots/thetaee_ideal_resolution.pdf
  plots/thetaee_detector_resolution.pdf
  plots/thetaee_resolution_comparison.pdf
  plots/acceptance_ideal_geometry_detector_vs_thetaee.pdf
  plots/detector_efficiency_given_geometry_vs_thetaee.pdf
  plots/hit_detector_occupancy_raw.pdf
  plots/hit_detector_occupancy_eventlevel.pdf
  plots/hit_detector_occupancy_reconstruction_used.pdf
  plots/scint_energy_em_ep.pdf
  plots/scint_energy_sum.pdf
  plots/thetaee_vs_Esum_reco.pdf
  plots/thetaee_vs_Y_reco.pdf
  plots/cutflow.pdf
)

for f in "${required[@]}"; do
    [[ -f "$f" ]] || { echo "ERROR: missing expected plot: $f"; exit 1; }
done

echo "==> Analysis finished"
echo "==> ROOT output: x17_analysis.root"
echo "==> PDF plots: plots/"
