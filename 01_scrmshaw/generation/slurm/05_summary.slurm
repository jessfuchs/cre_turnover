#!/usr/bin/env bash
#SBATCH --job-name=scrm_summary
#SBATCH --nodelist=abacus-2
#SBATCH --mem=5G
#SBATCH --output=logs/summary_%j.out
#SBATCH --error=logs/summary_%j.err
set -euo pipefail

source config/config.sh
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$POST_ENV"

python scripts/05_make_summary.py --manifest "$MANIFEST" --results "$RESULTS_DIR"
