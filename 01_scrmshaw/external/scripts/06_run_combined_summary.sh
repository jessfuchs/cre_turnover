#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/external_config.sh"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$POST_ENV"
python "$OG_PIPELINE_ROOT/scripts/make_summary.py" \
  --manifest "$COMBINED_MANIFEST" \
  --results "$COMBINED_RESULTS_DIR"
