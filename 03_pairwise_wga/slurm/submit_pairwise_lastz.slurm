#!/usr/bin/env bash
#SBATCH --job-name=dmel_lastz
#SBATCH --nodelist=abacus-2
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-39%8
#SBATCH --output=/home/jessica/cre_turnover/03_pairwise_wga/logs/lastz_%A_%a.out
#SBATCH --error=/home/jessica/cre_turnover/03_pairwise_wga/logs/lastz_%A_%a.err

set -euo pipefail


# ============================================================
# 02 - Submit pairwise whole-genome alignments
#
# Purpose:
#   Run the pairwise D. melanogaster-to-target LASTZ alignment
#   for one species selected by the SLURM array index.
#
# Input:
#   - target_species.txt
#
# Output:
#   - alignments_dmel/<species>/
#   - logs/lastz_<job>_<task>.out/.err
#
# Configuration:
#   config/wga_config.sh
# ============================================================


# ============================================================
# Configuration
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

source "$PIPELINE_ROOT/config/wga_config.sh"

export PATH="$BIN_DIR:$PATH"


# ============================================================
# Select target species
# ============================================================

SP="$(
    sed -n "${SLURM_ARRAY_TASK_ID}p" \
        "$TARGET_SPECIES_FILE"
)"

[[ -n "$SP" ]] || {
    echo "ERROR: no species found for array task ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
}


# ============================================================
# Run pairwise alignment
# ============================================================

echo "Species: $SP"
echo "Node   : $(hostname)"
echo "Start  : $(date)"

"$PIPELINE_ROOT/scripts/02_run_pairwise_lastz.sh" "$SP"

echo "End    : $(date)"
