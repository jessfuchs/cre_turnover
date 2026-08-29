#!/usr/bin/env bash
#SBATCH --job-name=post_external
#SBATCH --nodelist=abacus-2
#SBATCH --mem=5G
#SBATCH --output=logs/post_external_%A_%a.out
#SBATCH --error=logs/post_external_%A_%a.err

set -euo pipefail

scripts/postprocess_external_species.sh
