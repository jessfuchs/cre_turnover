#!/usr/bin/env bash
#SBATCH --job-name=scrm_offsets
#SBATCH --nodelist=abacus-2
#SBATCH --mem=5G
#SBATCH --output=logs/offset_%A_%a.out
#SBATCH --error=logs/offset_%A_%a.err
set -euo pipefail
scripts/run_one_offset.sh
