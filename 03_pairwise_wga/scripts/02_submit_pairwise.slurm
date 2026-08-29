#!/usr/bin/env bash
#SBATCH --job-name=dmel_lastz
#SBATCH --nodelist=abacus-2
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-39%8
#SBATCH --output=/home/jessica/cre_turnover/project/pairwise_wga/logs/lastz_%A_%a.out
#SBATCH --error=/home/jessica/cre_turnover/project/pairwise_wga/logs/lastz_%A_%a.err

set -euo pipefail

export PATH="$HOME/cre_turnover/project/pairwise_wga/bin:$PATH"

ROOT="$HOME/cre_turnover/project/pairwise_wga"

SP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ROOT/target_species.txt")

echo "Species: $SP"
echo "Node: $(hostname)"
echo "Start: $(date)"

"$ROOT/scripts/02_run_pairwise_lastz.sh" "$SP"

echo "End: $(date)"
