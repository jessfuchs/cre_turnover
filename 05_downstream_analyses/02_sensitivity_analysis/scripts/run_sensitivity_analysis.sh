#!/usr/bin/env bash

# ============================================================
# Run sensitivity-analysis pipeline
#
# Purpose:
#   Execute the complete CRE-classification sensitivity
#   workflow and validate the expected final outputs.
#
# Pipeline steps:
#   1. Generate sensitivity classifications
#   2. Verify primary-scenario regression
#   3. Summarize global sensitivity
#   4. Summarize candidate sensitivity
#   5. Add sensitivity evidence to prioritized candidates
#
# Paths and analysis parameters are defined in:
#   config/sensitivity_config.sh
# ============================================================

set -euo pipefail


# ============================================================
# Pipeline root and configuration
# ============================================================

# run_sensitivity_analysis.sh is located directly in the
# sensitivity-analysis root directory.
PIPELINE_ROOT="$(
    cd "$(dirname "${BASH_SOURCE[0]}")"
    pwd
)"

CONFIG_FILE="${PIPELINE_ROOT}/config/sensitivity_config.sh"


if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: sensitivity configuration not found:"
    echo "${CONFIG_FILE}"
    exit 1
fi


# shellcheck source=/dev/null
source "${CONFIG_FILE}"


# ============================================================
# Logging
# ============================================================

LOG_DIR="${RESULTS_DIR}/logs"

mkdir -p \
    "${RESULTS_DIR}" \
    "${LOG_DIR}"


TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"

LOG_FILE="${LOG_DIR}/sensitivity_analysis_${TIMESTAMP}.log"


timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}


log() {
    echo "[$(timestamp)] $*" \
        | tee -a "${LOG_FILE}"
}


# ============================================================
# Run one pipeline step
# ============================================================

run_step() {

    local step="$1"
    local total="$2"
    local label="$3"

    shift 3

    log "============================================================"
    log "Step ${step}/${total}: ${label}"
    log "============================================================"

    local start
    start="$(date +%s)"

    set +e

    "$@" 2>&1 \
        | tee -a "${LOG_FILE}"

    local statuses=("${PIPESTATUS[@]}")
    local command_status="${statuses[0]:-1}"
    local tee_status="${statuses[1]:-1}"

    set -e

    local end
    end="$(date +%s)"

    local runtime=$((end - start))


    if [[ "${command_status}" -ne 0 ]]; then

        log ""
        log "ERROR: ${label} failed"
        log "Exit status: ${command_status}"
        log "Runtime before failure: ${runtime} s"
        log "Pipeline stopped."

        exit "${command_status}"
    fi


    if [[ "${tee_status}" -ne 0 ]]; then

        log ""
        log "ERROR: logging failed during ${label}"
        log "tee exit status: ${tee_status}"

        exit "${tee_status}"
    fi


    log ""
    log "PASS: ${label}"
    log "Runtime: ${runtime} s"
    log ""
}


# ============================================================
# Pipeline scripts
# ============================================================

SCRIPT_01="${PIPELINE_ROOT}/01_generate_sensitivity_scenarios.py"
SCRIPT_02="${PIPELINE_ROOT}/02_check_primary_regression.py"
SCRIPT_03="${PIPELINE_ROOT}/03_summarize_global_sensitivity.py"
SCRIPT_04="${PIPELINE_ROOT}/04_summarize_candidate_sensitivity.py"
SCRIPT_05="${PIPELINE_ROOT}/05_add_sensitivity_to_prioritized_candidates.py"


for script in \
    "${SCRIPT_01}" \
    "${SCRIPT_02}" \
    "${SCRIPT_03}" \
    "${SCRIPT_04}" \
    "${SCRIPT_05}"
do

    if [[ ! -f "${script}" ]]; then
        echo "ERROR: required pipeline script not found:"
        echo "${script}"
        exit 1
    fi

done


# ============================================================
# Start
# ============================================================

log "============================================================"
log "Sensitivity-analysis pipeline started"
log "============================================================"

log ""
log "Pipeline root: ${PIPELINE_ROOT}"
log "Config:        ${CONFIG_FILE}"
log "Results:       ${RESULTS_DIR}"
log "Log file:      ${LOG_FILE}"
log ""

log "Python executable: $(command -v python3)"
log "Python version: $(python3 --version 2>&1)"

if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
    log "Conda environment: ${CONDA_DEFAULT_ENV}"
fi

log ""


# ============================================================
# Step 1: Generate sensitivity scenarios
# ============================================================

run_step \
    1 \
    5 \
    "Generate sensitivity scenarios" \
    python3 "${SCRIPT_01}" \
        --classifier "${CLASSIFIER_SCRIPT}" \
        --targets "${TARGET_SPECIES_FILE}" \
        --reference "${REFERENCE_CRES_TSV}" \
        --predictions "${SO_ALL_SPECIES_FBGN}" \
        --lifted-dir "${LIFTED_CRES_DIR}" \
        --out-dir "${SENSITIVITY_SCENARIOS_DIR}" \
        --scenario-manifest "${SENSITIVITY_SCENARIO_MANIFEST}" \
        --metadata-out "${SENSITIVITY_GENERATION_METADATA}" \
        --overlaps "${SENSITIVITY_OVERLAPS[@]}" \
        --distances "${SENSITIVITY_DISTANCES[@]}" \
        --expected-reference-cres "${EXPECTED_REFERENCE_CRES}"


# ============================================================
# Step 2: Primary regression check
# ============================================================

run_step \
    2 \
    5 \
    "Check primary regression" \
    python3 "${SCRIPT_02}" \
        --baseline-dir "${BASELINE_TURNOVER_DIR}" \
        --sensitivity-dir "${SENSITIVITY_SCENARIOS_DIR}" \
        --targets "${TARGET_SPECIES_FILE}" \
        --scenario-manifest "${SENSITIVITY_SCENARIO_MANIFEST}" \
        --primary-scenario "${PRIMARY_SENSITIVITY_SCENARIO}" \
        --out-summary "${PRIMARY_REGRESSION_SUMMARY}" \
        --out-changes "${PRIMARY_REGRESSION_CHANGES}" \
        --metadata-out "${PRIMARY_REGRESSION_METADATA}" \
        --expected-reference-cres "${EXPECTED_REFERENCE_CRES}"


# ============================================================
# Step 3: Global sensitivity summary
# ============================================================

run_step \
    3 \
    5 \
    "Summarize global sensitivity" \
    python3 "${SCRIPT_03}" \
        --sensitivity-dir "${SENSITIVITY_SCENARIOS_DIR}" \
        --targets "${TARGET_SPECIES_FILE}" \
        --scenario-manifest "${SENSITIVITY_SCENARIO_MANIFEST}" \
        --primary-scenario "${PRIMARY_SENSITIVITY_SCENARIO}" \
        --out-global "${GLOBAL_SENSITIVITY_SUMMARY}" \
        --out-species "${SENSITIVITY_SPECIES_STABILITY}" \
        --out-transitions "${SENSITIVITY_STATE_TRANSITIONS}" \
        --out-cre "${SENSITIVITY_CRE_STABILITY}" \
        --metadata-out "${GLOBAL_SENSITIVITY_METADATA}" \
        --expected-reference-cres "${EXPECTED_REFERENCE_CRES}"


# ============================================================
# Step 4: Candidate sensitivity
# ============================================================

run_step \
    4 \
    5 \
    "Summarize candidate sensitivity" \
    python3 "${SCRIPT_04}" \
        --sensitivity-dir "${SENSITIVITY_SCENARIOS_DIR}" \
        --scenario-manifest "${SENSITIVITY_SCENARIO_MANIFEST}" \
        --groups "${FOCAL_CLADES_FILE}" \
        --candidates "${CANDIDATE_TIER1_TABLE}" \
        --primary-scenario "${PRIMARY_SENSITIVITY_SCENARIO}" \
        --out-long "${CANDIDATE_SENSITIVITY_LONG}" \
        --out-summary "${CANDIDATE_SENSITIVITY_SUMMARY}" \
        --out-priority-matrix "${CANDIDATE_SENSITIVITY_PRIORITY_MATRIX}" \
        --out-retention-matrix "${CANDIDATE_SENSITIVITY_RETENTION_MATRIX}" \
        --metadata-out "${CANDIDATE_SENSITIVITY_METADATA}" \
        --focal-recurrent-min-clades "${FOCAL_RECURRENCE_MIN_CLADES}" \
        --secondary-recurrent-min-clades "${SECONDARY_RECURRENCE_MIN_CLADES}" \
        --expected-reference-cres "${EXPECTED_REFERENCE_CRES}"


# ============================================================
# Step 5: Add sensitivity to prioritized candidates
# ============================================================

run_step \
    5 \
    5 \
    "Add sensitivity to prioritized candidates" \
    python3 "${SCRIPT_05}" \
        --candidates "${CANDIDATE_TIER1_TABLE}" \
        --sensitivity "${CANDIDATE_SENSITIVITY_SUMMARY}" \
        --out "${CANDIDATES_WITH_SENSITIVITY}" \
        --metadata-out "${CANDIDATE_SENSITIVITY_MERGE_METADATA}" \
        --moderate-robustness-min "${CANDIDATE_MODERATE_ROBUSTNESS_MIN}"


# ============================================================
# Expected final outputs
# ============================================================

EXPECTED_OUTPUTS=(

    # Scenario generation
    "${SENSITIVITY_SCENARIO_MANIFEST}"
    "${SENSITIVITY_GENERATION_METADATA}"

    # Primary regression
    "${PRIMARY_REGRESSION_SUMMARY}"
    "${PRIMARY_REGRESSION_CHANGES}"
    "${PRIMARY_REGRESSION_METADATA}"

    # Global sensitivity
    "${GLOBAL_SENSITIVITY_SUMMARY}"
    "${SENSITIVITY_SPECIES_STABILITY}"
    "${SENSITIVITY_STATE_TRANSITIONS}"
    "${SENSITIVITY_CRE_STABILITY}"
    "${GLOBAL_SENSITIVITY_METADATA}"

    # Candidate sensitivity
    "${CANDIDATE_SENSITIVITY_LONG}"
    "${CANDIDATE_SENSITIVITY_SUMMARY}"
    "${CANDIDATE_SENSITIVITY_PRIORITY_MATRIX}"
    "${CANDIDATE_SENSITIVITY_RETENTION_MATRIX}"
    "${CANDIDATE_SENSITIVITY_METADATA}"

    # Sensitivity-annotated primary candidates
    "${CANDIDATES_WITH_SENSITIVITY}"
    "${CANDIDATE_SENSITIVITY_MERGE_METADATA}"
)


log "============================================================"
log "Checking final sensitivity outputs"
log "============================================================"


N_MISSING=0


for file in "${EXPECTED_OUTPUTS[@]}"; do

    if [[ -s "${file}" ]]; then
        log "FOUND: ${file}"
    else
        log "MISSING/EMPTY: ${file}"
        N_MISSING=$((N_MISSING + 1))
    fi

done


if [[ "${N_MISSING}" -gt 0 ]]; then

    log ""
    log "ERROR: ${N_MISSING} expected output(s) missing."

    exit 1
fi


# ============================================================
# Finish
# ============================================================

log ""
log "============================================================"
log "Sensitivity-analysis pipeline completed successfully"
log "============================================================"

log ""
log "Global sensitivity:"
log "${GLOBAL_SENSITIVITY_SUMMARY}"

log ""
log "Candidate sensitivity:"
log "${CANDIDATE_SENSITIVITY_SUMMARY}"

log ""
log "Sensitivity-annotated candidates:"
log "${CANDIDATES_WITH_SENSITIVITY}"

log ""
log "Full log:"
log "${LOG_FILE}"
