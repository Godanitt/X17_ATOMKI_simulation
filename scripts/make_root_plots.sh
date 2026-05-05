#!/usr/bin/env bash
set -Eeuo pipefail

# ROOT plots for input angular table. No options.
# Usage from anywhere:
#   bash scripts/make_root_plots.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${PROJECT_DIR}"

echo "==> Making ROOT angular input plots"

command -v root >/dev/null 2>&1 || { echo "ERROR: ROOT not found. Source ROOT first."; exit 1; }
[[ -f "data/data_pair_creation.txt" ]] || { echo "ERROR: data/data_pair_creation.txt not found"; exit 1; }
[[ -f "scripts/plot_angles_ROOT.C" ]] || { echo "ERROR: scripts/plot_angles_ROOT.C not found"; exit 1; }
[[ -f "scripts/X17Style.C" ]] || { echo "ERROR: scripts/X17Style.C not found"; exit 1; }

mkdir -p plots

root -l -b -q -I scripts <<'ROOTEOF'
.L scripts/plot_angles_ROOT.C
plot_angles_ROOT("data/data_pair_creation.txt", "plots");
.q
ROOTEOF

[[ -f "plots/thetaee_distribution_ROOT.pdf" ]] || { echo "ERROR: plots/thetaee_distribution_ROOT.pdf was not created"; exit 1; }
[[ -f "plots/angular_distributions_ROOT.pdf" ]] || { echo "ERROR: plots/angular_distributions_ROOT.pdf was not created"; exit 1; }

echo "==> ROOT angular plots done"
echo "  plots/thetaee_distribution_ROOT.pdf"
echo "  plots/angular_distributions_ROOT.pdf"
