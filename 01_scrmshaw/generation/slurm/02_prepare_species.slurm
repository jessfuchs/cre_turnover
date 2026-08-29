#!/usr/bin/env bash
#SBATCH --job-name=scrm_prep
#SBATCH --nodelist=abacus-2
#SBATCH --mem=5G
#SBATCH --output=logs/prep_%A_%a.out
#SBATCH --error=logs/prep_%A_%a.err
set -euo pipefail
scripts/prepare_one_species.sh
