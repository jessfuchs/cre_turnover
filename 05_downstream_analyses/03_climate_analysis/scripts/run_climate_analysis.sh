#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo
echo "============================================================"
echo "Climate and focal Tier-1 analysis"
echo "============================================================"
echo

echo "Step 1/4: Turnover rate by climate"
python "$SCRIPT_DIR/01_plot_turnover_by_climate.py"

echo
echo "Step 2/4: Phylogenetically controlled climate test"
Rscript "$SCRIPT_DIR/02_test_turnover_by_climate.R"

echo
echo "Step 3/4: Focal Tier-1 enrichment"
python "$SCRIPT_DIR/03_test_focal_tier1_enrichment.py"

echo
echo "Step 4/4: Focal Tier-1 enrichment plot"
python "$SCRIPT_DIR/04_plot_focal_tier1_enrichment.py"

echo
echo "============================================================"
echo "ALL ANALYSES COMPLETE"
echo "============================================================"
