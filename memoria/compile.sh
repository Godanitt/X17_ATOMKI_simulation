#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
typst compile paper.typ paper.pdf
