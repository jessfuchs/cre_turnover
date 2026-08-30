#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Pairwise whole-genome alignment pipeline
#
# Purpose:
#   Prepare genome inputs, submit pairwise D. melanogaster
#   alignments for all target species, and run CRE liftOver
#   after all alignment jobs have completed successfully.
#
# Configuration:
#   config/wga_config.sh
# ============================================================

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "$ROOT/config/wga_config.sh"

export PATH="$BIN_DIR:$PATH"


# ============================================================
# 1. Prepare genomes
# ============================================================

echo "[1/3] Preparing genome inputs..."

bash "$ROOT/scripts/01_prepare_wga_genomes.sh"

echo "[1/3] PASSED"
echo


# ============================================================
# 2. Submit pairwise alignments
# ============================================================

echo "[2/3] Submitting pairwise LASTZ alignments..."

N_TARGETS="$(
    awk '
        /^[[:space:]]*#/ {next}
        NF {n++}
        END {print n+0}
    ' "$TARGET_SPECIES_FILE"
)"

if [[ "$N_TARGETS" -eq 0 ]]; then
    echo "ERROR: no target species found." >&2
    exit 1
fi

ALIGNMENT_JOB="$(
    sbatch \
        --parsable \
        --array="1-${N_TARGETS}%${MAX_WGA_JOBS}" \
        --nodelist="$SLURM_NODE" \
        --mem="$WGA_MEM" \
        --cpus-per-task="$WGA_CPUS" \
        --output="$LOG_DIR/lastz_%A_%a.out" \
        --error="$LOG_DIR/lastz_%A_%a.err" \
        "$ROOT/slurm/02_submit_pairwise.slurm"
)"

echo "Submitted alignment job: $ALIGNMENT_JOB"
echo


# ============================================================
# 3. Submit liftOver after successful alignments
# ============================================================

echo "[3/3] Submitting liftOver dependency..."

LIFTOVER_JOB="$(
    sbatch \
        --parsable \
        --dependency="afterok:${ALIGNMENT_JOB}" \
        --job-name=dmel_liftover \
        --nodelist="$SLURM_NODE" \
        --mem="$WGA_MEM" \
        --cpus-per-task=1 \
        --output="$LOG_DIR/liftover_%j.out" \
        --error="$LOG_DIR/liftover_%j.err" \
        --wrap="bash '$ROOT/scripts/03_run_liftover.sh'"
)"

echo "Submitted liftOver job: $LIFTOVER_JOB"
echo

echo "============================================================"
echo "Pairwise WGA pipeline submitted"
echo "============================================================"
echo "Alignment job : $ALIGNMENT_JOB"
echo "liftOver job  : $LIFTOVER_JOB"
echo
echo "liftOver will start only if all alignment tasks finish"
echo "successfully."
echo "============================================================"
