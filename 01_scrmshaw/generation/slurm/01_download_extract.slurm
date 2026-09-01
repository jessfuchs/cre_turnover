#!/usr/bin/env bash
#SBATCH --job-name=scrm_download
#SBATCH --nodelist=abacus-2
#SBATCH --mem=5G
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err
set -euo pipefail

source config/config.sh
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$SCRM_ENV"

python scripts/01_download_extract_zenodo.py \
    --species "$SPECIES_FILE" \
    --data-dir "$DATA_DIR" \
    --record-id "$ZENODO_RECORD_ID" \
    --genomes-md5 "$GENOMES_MD5" \
    --annotations-md5 "$ANNOTATIONS_MD5"
