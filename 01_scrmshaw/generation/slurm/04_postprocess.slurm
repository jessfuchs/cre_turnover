#!/usr/bin/env bash
#SBATCH --job-name=scrm_post
#SBATCH --nodelist=abacus-2
#SBATCH --mem=5G
#SBATCH --output=logs/post_%A_%a.out
#SBATCH --error=logs/post_%A_%a.err
set -euo pipefail
scripts/04_postprocess_one_species.sh
