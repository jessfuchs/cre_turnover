#!/usr/bin/env bash

set -euo pipefail


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_DIR="${HOME}/cre_turnover/project"

ANALYSIS_DIR="${PROJECT_DIR}/downstream_analyses/sensitivity_analysis"

RESULTS_DIR="${ANALYSIS_DIR}/results"

LOG_DIR="${RESULTS_DIR}/logs"

mkdir -p \
    "${RESULTS_DIR}" \
    "${LOG_DIR}"


TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"

LOG_FILE="${LOG_DIR}/sensitivity_analysis_${TIMESTAMP}.log"


SCRIPTS=(
    "01_generate_sensitivity_scenarios.py"
    "02_check_primary_regression.py"
    "03_summarize_global_sensitivity.py"
    "04_summarize_candidate_sensitivity.py"
)


timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}


log() {
    echo "[$(timestamp)] $*" \
        | tee -a "${LOG_FILE}"
}


run_step() {

    local step="$1"
    local total="$2"
    local script="$3"

    local path="${SCRIPT_DIR}/${script}"


    log "============================================================"
    log "Step ${step}/${total}: ${script}"
    log "============================================================"


    local start

    start="$(date +%s)"


    set +e

    python "${path}" \
        2>&1 \
        | tee -a "${LOG_FILE}"

    local status="${PIPESTATUS[0]}"

    set -e


    local end

    end="$(date +%s)"

    local runtime=$((end - start))


    if [[ "${status}" -ne 0 ]]; then

        log ""
        log "ERROR: ${script} failed"
        log "Exit status: ${status}"
        log "Runtime before failure: ${runtime} s"
        log ""
        log "Pipeline stopped."

        exit "${status}"

    fi


    log ""
    log "PASS: ${script}"
    log "Runtime: ${runtime} s"
    log ""
}


# ============================================================
# Start
# ============================================================

log "============================================================"
log "Sensitivity-analysis pipeline started"
log "============================================================"

log ""

log "Script directory: ${SCRIPT_DIR}"
log "Results directory: ${RESULTS_DIR}"
log "Log file: ${LOG_FILE}"

log ""

log "Python executable: $(command -v python)"
log "Python version: $(python --version 2>&1)"

if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
    log "Conda environment: ${CONDA_DEFAULT_ENV}"
fi

log ""


# ============================================================
# Check scripts
# ============================================================

for script in "${SCRIPTS[@]}"; do

    if [[ ! -f "${SCRIPT_DIR}/${script}" ]]; then

        log "ERROR: required script missing:"
        log "${SCRIPT_DIR}/${script}"

        exit 1
    fi

done


# ============================================================
# Run
# ============================================================

TOTAL="${#SCRIPTS[@]}"

STEP=0


for script in "${SCRIPTS[@]}"; do

    STEP=$((STEP + 1))

    run_step \
        "${STEP}" \
        "${TOTAL}" \
        "${script}"

done


# ============================================================
# Expected final outputs
# ============================================================

EXPECTED_OUTPUTS=(
    "${RESULTS_DIR}/primary_regression_check.tsv"
    "${RESULTS_DIR}/primary_regression_state_changes.tsv"

    "${RESULTS_DIR}/sensitivity_global_summary.tsv"
    "${RESULTS_DIR}/sensitivity_species_stability.tsv"
    "${RESULTS_DIR}/sensitivity_state_transitions.tsv"
    "${RESULTS_DIR}/sensitivity_cre_stability.tsv"

    "${RESULTS_DIR}/candidate_sensitivity_long.tsv"
    "${RESULTS_DIR}/candidate_sensitivity_summary.tsv"
    "${RESULTS_DIR}/candidate_sensitivity_priority_matrix.tsv"
    "${RESULTS_DIR}/candidate_sensitivity_retention_matrix.tsv"
)


log "============================================================"
log "Checking final sensitivity outputs"
log "============================================================"


MISSING=0


for file in "${EXPECTED_OUTPUTS[@]}"; do

    if [[ -s "${file}" ]]; then

        log "FOUND: ${file}"

    else

        log "MISSING/EMPTY: ${file}"

        MISSING=$((MISSING + 1))

    fi

done


if [[ "${MISSING}" -gt 0 ]]; then

    log ""
    log "ERROR: ${MISSING} expected output(s) missing."

    exit 1

fi


log ""

log "============================================================"
log "Sensitivity-analysis pipeline completed successfully"
log "============================================================"

log ""

log "Global sensitivity:"
log "${RESULTS_DIR}/sensitivity_global_summary.tsv"

log ""

log "Candidate sensitivity:"
log "${RESULTS_DIR}/candidate_sensitivity_summary.tsv"

log ""

log "Full log:"
log "${LOG_FILE}"
