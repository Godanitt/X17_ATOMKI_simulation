#!/usr/bin/env bash
set -Eeuo pipefail

# Same transport-only run as run_generated.sh.
# It writes sampled.root with exactly two trees: generated and hits.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${SCRIPT_DIR}/run_generated.sh" "$@"
