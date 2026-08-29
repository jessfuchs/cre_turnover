#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Paths
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_DIR="${HOME}/cre_turnover/project"

CANDIDATE_DIR="${PROJECT_DIR}/downstream_analyses/candidate_analysis"

RESULTS_DIR="${CANDIDATE_DIR}/results"

LOG_DIR="${RESULTS_DIR}/logs"

mkdir -p "${LOG_DIR}"


# ============================================================
# Log file
# ============================================================

TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"

LOG_FILE="${LOG_DIR}/candidate_analysis_pipeline_${TIMESTAMP}.log"


# ============================================================
# Pipeline scripts
# ============================================================

SCRIPTS=(
    "01_analyze_focal_clades.py"
    "02_analyze_secondary_tier1.py"
    "03_qc_tier1_candidates.py"
    "04_prioritize_tier1_candidates.py"
    "05_annotate_candidate_gene_distances.py"
    "06_assign_candidate_genes.py"
    "07_explode_candidate_fbgns.py"
    "08_build_gene_level_summary.py"
)


# ============================================================
# Helper
# ============================================================

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}


log() {
    echo "[$(timestamp)] $*" | tee -a "${LOG_FILE}"
}


# ============================================================
# Start
# ============================================================

log "============================================================"
log "Candidate-analysis pipeline started"
log "============================================================"

log "Script directory: ${SCRIPT_DIR}"
log "Project directory: ${PROJECT_DIR}"
log "Log file: ${LOG_FILE}"

log ""

log "Python executable: $(command -v python)"
log "Python version: $(python --version 2>&1)"

if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
    log "Conda environment: ${CONDA_DEFAULT_ENV}"
else
    log "Conda environment: not detected"
fi

log ""


# ============================================================
# Check all scripts before starting
# ============================================================

log "Checking pipeline scripts..."

for script in "${SCRIPTS[@]}"; do

    script_path="${SCRIPT_DIR}/${script}"

    if [[ ! -f "${script_path}" ]]; then

        log "ERROR: required script not found:"
        log "${script_path}"

        exit 1
    fi

    log "FOUND: ${script}"

done

log ""
log "All required scripts found."
log ""


# ============================================================
# Run pipeline
# ============================================================

N_STEPS="${#SCRIPTS[@]}"
STEP=0


for script in "${SCRIPTS[@]}"; do

    STEP=$((STEP + 1))

    script_path="${SCRIPT_DIR}/${script}"

    log "============================================================"
    log "Step ${STEP}/${N_STEPS}: ${script}"
    log "============================================================"

    START_SECONDS="$(date +%s)"

    # --------------------------------------------------------
    # Run Python script
    #
    # stderr and stdout both go to:
    #   - terminal
    #   - pipeline log
    # --------------------------------------------------------

    if python "${script_path}" 2>&1 | tee -a "${LOG_FILE}"; then

        END_SECONDS="$(date +%s)"

        RUNTIME=$(
            awk \
                -v start="${START_SECONDS}" \
                -v end="${END_SECONDS}" \
                'BEGIN { printf "%.1f", end-start }'
        )

        log ""
        log "PASS: ${script}"
        log "Runtime: ${RUNTIME} s"

    else

        STATUS="${PIPESTATUS[0]}"

        END_SECONDS="$(date +%s)"

        RUNTIME=$(
            awk \
                -v start="${START_SECONDS}" \
                -v end="${END_SECONDS}" \
                'BEGIN { printf "%.1f", end-start }'
        )

        log ""
        log "ERROR: ${script} failed"
        log "Exit status: ${STATUS}"
        log "Runtime before failure: ${RUNTIME} s"

        log ""
        log "Pipeline stopped."

        exit "${STATUS}"
    fi

    log ""

done


# ============================================================
# Final output checks
# ============================================================

log "============================================================"
log "Checking final outputs"
log "============================================================"


EXPECTED_OUTPUTS=(
    "${RESULTS_DIR}/focal_clade_summary.tsv"
    "${RESULTS_DIR}/secondary_clades/secondary_tier1_all_clades.tsv"
    "${RESULTS_DIR}/qc/tier1_candidate_qc_summary.tsv"
    "${RESULTS_DIR}/tables/tier1_candidates_prioritized.tsv"
    "${RESULTS_DIR}/tables/tier1_candidate_gene_distances.tsv"
    "${RESULTS_DIR}/tables/tier1_candidate_gene_assignments.tsv"
    "${RESULTS_DIR}/tables/tier1_candidate_fbgn_exploded.tsv"
    "${RESULTS_DIR}/tables/tier1_gene_level_summary.tsv"
)


N_MISSING=0


for outfile in "${EXPECTED_OUTPUTS[@]}"; do

    if [[ -s "${outfile}" ]]; then

        log "FOUND: ${outfile}"

    else

        log "MISSING/EMPTY: ${outfile}"

        N_MISSING=$((N_MISSING + 1))
    fi

done


log ""


if [[ "${N_MISSING}" -gt 0 ]]; then

    log "WARNING: pipeline scripts finished, but ${N_MISSING} expected output(s) are missing or empty."
    exit 1

fi


# ============================================================
# Finish
# ============================================================

log "============================================================"
log "Candidate-analysis pipeline completed successfully"
log "============================================================"

log ""
log "Final gene-level table:"
log "${RESULTS_DIR}/tables/tier1_gene_level_summary.tsv"

log ""
log "Full pipeline log:"
log "${LOG_FILE}"
