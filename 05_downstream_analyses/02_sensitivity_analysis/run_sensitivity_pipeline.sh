#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# CRE-classification sensitivity-analysis pipeline
#
# Steps:
#   1. Generate sensitivity scenarios
#   2. Check primary-scenario regression
#   3. Summarize global sensitivity
#   4. Summarize candidate sensitivity
#   5. Add sensitivity to prioritized candidates
#   6. Summarize parameter-specific sensitivity
#   7. Plot candidate sensitivity matrix
#   8. Plot ranked Tier-1 robustness
#   9. Plot recurrence versus robustness
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${PIPELINE_ROOT}/scripts"
CONFIG_FILE="${PIPELINE_ROOT}/config/sensitivity_config.sh"

[[ -s "$CONFIG_FILE" ]] || {
    echo "ERROR: configuration file missing or empty:" >&2
    echo "  $CONFIG_FILE" >&2
    exit 1
}

# shellcheck source=/dev/null
source "$CONFIG_FILE"

LOG_DIR="${RESULTS_DIR}/logs"

mkdir -p \
    "$RESULTS_DIR" \
    "$SENSITIVITY_SCENARIOS_DIR" \
    "$FIGURES_DIR" \
    "$LOG_DIR"

LOG_FILE="${LOG_DIR}/sensitivity_analysis_$(date '+%Y%m%d_%H%M%S').log"


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

run_step() {
    local step="$1"
    local label="$2"
    shift 2

    log "============================================================"
    log "[${step}/9] ${label}"
    log "============================================================"

    local start status tee_status
    start="$(date +%s)"

    set +e
    "$@" 2>&1 | tee -a "$LOG_FILE"
    status="${PIPESTATUS[0]}"
    tee_status="${PIPESTATUS[1]}"
    set -e

    if [[ "$status" -ne 0 ]]; then
        log "ERROR: ${label} failed (exit ${status})."
        exit "$status"
    fi

    if [[ "$tee_status" -ne 0 ]]; then
        log "ERROR: logging failed during ${label}."
        exit "$tee_status"
    fi

    log "PASS: ${label} ($(( $(date +%s) - start )) s)"
    log ""
}


# ------------------------------------------------------------
# Pipeline scripts
# ------------------------------------------------------------

SCRIPT_01="${SCRIPTS_DIR}/01_generate_sensitivity_scenarios.py"
SCRIPT_02="${SCRIPTS_DIR}/02_check_primary_regression.py"
SCRIPT_03="${SCRIPTS_DIR}/03_summarize_global_sensitivity.py"
SCRIPT_04="${SCRIPTS_DIR}/04_summarize_candidate_sensitivity.py"
SCRIPT_05="${SCRIPTS_DIR}/05_add_sensitivity_to_prioritized_candidates.py"
SCRIPT_06="${SCRIPTS_DIR}/06_summarize_parameter_sensitivity.py"
SCRIPT_07="${SCRIPTS_DIR}/07_plot_sensitivity_matrix.py"
SCRIPT_08="${SCRIPTS_DIR}/08_plot_ranked_tier1.py"
SCRIPT_09="${SCRIPTS_DIR}/09_plot_recurrence_robustness.py"


# ------------------------------------------------------------
# Preconditions
# ------------------------------------------------------------

command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 not found." >&2
    exit 1
}

REQUIRED_FILES=(
    "$SCRIPT_01"
    "$SCRIPT_02"
    "$SCRIPT_03"
    "$SCRIPT_04"
    "$SCRIPT_05"
    "$SCRIPT_06"
    "$SCRIPT_07"
    "$SCRIPT_08"
    "$SCRIPT_09"

    "$CLASSIFIER_SCRIPT"
    "$TARGET_SPECIES_FILE"
    "$REFERENCE_CRES_TSV"
    "$SO_ALL_SPECIES_FBGN"
    "$CANDIDATE_TIER1_TABLE"
    "$CANDIDATE_TIER1_METADATA"
    "$FOCAL_CLADES_FILE"
)

for file in "${REQUIRED_FILES[@]}"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required file missing or empty:" >&2
        echo "  $file" >&2
        exit 1
    }
done

[[ -d "$LIFTED_CRES_DIR" ]] || {
    echo "ERROR: lifted-CRE directory not found:" >&2
    echo "  $LIFTED_CRES_DIR" >&2
    exit 1
}

[[ -d "$BASELINE_TURNOVER_DIR" ]] || {
    echo "ERROR: baseline turnover directory not found:" >&2
    echo "  $BASELINE_TURNOVER_DIR" >&2
    exit 1
}


# ============================================================
# Start
# ============================================================

log "============================================================"
log "Sensitivity-analysis pipeline"
log "============================================================"
log "Started: $(date)"
log "Pipeline root: $PIPELINE_ROOT"
log "Primary scenario: $PRIMARY_SENSITIVITY_SCENARIO"
log "Overlap grid: ${SENSITIVITY_OVERLAPS[*]}"
log "Distance grid: ${SENSITIVITY_DISTANCES[*]}"
log ""


# ============================================================
# 1. Generate sensitivity scenarios
# ============================================================

run_step 1 "Generate sensitivity scenarios" \
    python3 "$SCRIPT_01" \
        --classifier "$CLASSIFIER_SCRIPT" \
        --targets "$TARGET_SPECIES_FILE" \
        --reference "$REFERENCE_CRES_TSV" \
        --predictions "$SO_ALL_SPECIES_FBGN" \
        --lifted-dir "$LIFTED_CRES_DIR" \
        --out-dir "$SENSITIVITY_SCENARIOS_DIR" \
        --scenario-manifest "$SENSITIVITY_SCENARIO_MANIFEST" \
        --metadata-out "$SENSITIVITY_GENERATION_METADATA" \
        --overlaps "${SENSITIVITY_OVERLAPS[@]}" \
        --distances "${SENSITIVITY_DISTANCES[@]}" \
        --expected-reference-cres "$EXPECTED_REFERENCE_CRES"


# ============================================================
# 2. Primary regression
# ============================================================

run_step 2 "Check primary regression" \
    python3 "$SCRIPT_02" \
        --baseline-dir "$BASELINE_TURNOVER_DIR" \
        --sensitivity-dir "$SENSITIVITY_SCENARIOS_DIR" \
        --targets "$TARGET_SPECIES_FILE" \
        --scenario-manifest "$SENSITIVITY_SCENARIO_MANIFEST" \
        --primary-scenario "$PRIMARY_SENSITIVITY_SCENARIO" \
        --out-summary "$PRIMARY_REGRESSION_SUMMARY" \
        --out-changes "$PRIMARY_REGRESSION_CHANGES" \
        --metadata-out "$PRIMARY_REGRESSION_METADATA" \
        --expected-reference-cres "$EXPECTED_REFERENCE_CRES"


# ============================================================
# 3. Global sensitivity
# ============================================================

run_step 3 "Summarize global sensitivity" \
    python3 "$SCRIPT_03" \
        --sensitivity-dir "$SENSITIVITY_SCENARIOS_DIR" \
        --targets "$TARGET_SPECIES_FILE" \
        --scenario-manifest "$SENSITIVITY_SCENARIO_MANIFEST" \
        --primary-scenario "$PRIMARY_SENSITIVITY_SCENARIO" \
        --out-global "$GLOBAL_SENSITIVITY_SUMMARY" \
        --out-species "$SENSITIVITY_SPECIES_STABILITY" \
        --out-transitions "$SENSITIVITY_STATE_TRANSITIONS" \
        --out-cre "$SENSITIVITY_CRE_STABILITY" \
        --metadata-out "$GLOBAL_SENSITIVITY_METADATA" \
        --expected-reference-cres "$EXPECTED_REFERENCE_CRES"


# ============================================================
# 4. Candidate sensitivity
# ============================================================

run_step 4 "Summarize candidate sensitivity" \
    python3 "$SCRIPT_04" \
        --sensitivity-dir "$SENSITIVITY_SCENARIOS_DIR" \
        --scenario-manifest "$SENSITIVITY_SCENARIO_MANIFEST" \
        --groups "$FOCAL_CLADES_FILE" \
        --candidates "$CANDIDATE_TIER1_TABLE" \
        --candidate-metadata "$CANDIDATE_TIER1_METADATA" \
        --primary-scenario "$PRIMARY_SENSITIVITY_SCENARIO" \
        --out-long "$CANDIDATE_SENSITIVITY_LONG" \
        --out-summary "$CANDIDATE_SENSITIVITY_SUMMARY" \
        --out-priority-matrix "$CANDIDATE_SENSITIVITY_PRIORITY_MATRIX" \
        --out-retention-matrix "$CANDIDATE_SENSITIVITY_RETENTION_MATRIX" \
        --metadata-out "$CANDIDATE_SENSITIVITY_METADATA" \
        --focal-recurrent-min-clades "$FOCAL_RECURRENCE_MIN_CLADES" \
        --secondary-recurrent-min-clades "$SECONDARY_RECURRENCE_MIN_CLADES" \
        --expected-reference-cres "$EXPECTED_REFERENCE_CRES"


# ============================================================
# 5. Add sensitivity to candidate table
# ============================================================

run_step 5 "Add sensitivity to prioritized candidates" \
    python3 "$SCRIPT_05" \
        --candidates "$CANDIDATE_TIER1_TABLE" \
        --sensitivity "$CANDIDATE_SENSITIVITY_SUMMARY" \
        --out "$CANDIDATES_WITH_SENSITIVITY" \
        --metadata-out "$CANDIDATE_SENSITIVITY_MERGE_METADATA" \
        --moderate-robustness-min "$CANDIDATE_MODERATE_ROBUSTNESS_MIN"


# ============================================================
# 6. Parameter-specific sensitivity
# ============================================================

run_step 6 "Summarize parameter-specific sensitivity" \
    python3 "$SCRIPT_06" \
        --retention-matrix "$CANDIDATE_SENSITIVITY_RETENTION_MATRIX" \
        --global-summary "$GLOBAL_SENSITIVITY_SUMMARY" \
        --scenario-manifest "$SENSITIVITY_SCENARIO_MANIFEST" \
        --primary-scenario "$PRIMARY_SENSITIVITY_SCENARIO" \
        --out-candidates "$PARAMETER_SENSITIVITY" \
        --out-summary "$PARAMETER_SENSITIVITY_SUMMARY" \
        --out-classes "$PARAMETER_SENSITIVITY_CLASS_SUMMARY" \
        --metadata-out "$PARAMETER_SENSITIVITY_METADATA"


# ============================================================
# 7. Candidate sensitivity matrix
# ============================================================

run_step 7 "Plot candidate sensitivity matrix" \
    python3 "$SCRIPT_07" \
        --matrix "$CANDIDATE_SENSITIVITY_RETENTION_MATRIX" \
        --summary "$CANDIDATE_SENSITIVITY_SUMMARY" \
        --scenario-manifest "$SENSITIVITY_SCENARIO_MANIFEST" \
        --primary-scenario "$PRIMARY_SENSITIVITY_SCENARIO" \
        --out-png "$TIER1_SENSITIVITY_MATRIX_PNG" \
        --out-pdf "$TIER1_SENSITIVITY_MATRIX_PDF"


# ============================================================
# 8. Ranked Tier-1 robustness
# ============================================================

run_step 8 "Plot ranked Tier-1 robustness" \
    python3 "$SCRIPT_08" \
        --input "$CANDIDATES_WITH_SENSITIVITY" \
        --out-png "$TIER1_ROBUSTNESS_RANKED_PNG" \
        --out-pdf "$TIER1_ROBUSTNESS_RANKED_PDF" \
        --moderate-robustness-min "$CANDIDATE_MODERATE_ROBUSTNESS_MIN"


# ============================================================
# 9. Recurrence versus robustness
# ============================================================

run_step 9 "Plot recurrence versus robustness" \
    python3 "$SCRIPT_09" \
        --input "$CANDIDATES_WITH_SENSITIVITY" \
        --out-png "$RECURRENCE_ROBUSTNESS_PNG" \
        --out-pdf "$RECURRENCE_ROBUSTNESS_PDF" \
        --moderate-robustness-min "$CANDIDATE_MODERATE_ROBUSTNESS_MIN"


# ============================================================
# Final output QC
# ============================================================

EXPECTED_OUTPUTS=(
    "$SENSITIVITY_SCENARIO_MANIFEST"
    "$SENSITIVITY_GENERATION_METADATA"

    "$PRIMARY_REGRESSION_SUMMARY"
    "$PRIMARY_REGRESSION_CHANGES"
    "$PRIMARY_REGRESSION_METADATA"

    "$GLOBAL_SENSITIVITY_SUMMARY"
    "$SENSITIVITY_SPECIES_STABILITY"
    "$SENSITIVITY_STATE_TRANSITIONS"
    "$SENSITIVITY_CRE_STABILITY"
    "$GLOBAL_SENSITIVITY_METADATA"

    "$CANDIDATE_SENSITIVITY_LONG"
    "$CANDIDATE_SENSITIVITY_SUMMARY"
    "$CANDIDATE_SENSITIVITY_PRIORITY_MATRIX"
    "$CANDIDATE_SENSITIVITY_RETENTION_MATRIX"
    "$CANDIDATE_SENSITIVITY_METADATA"

    "$CANDIDATES_WITH_SENSITIVITY"
    "$CANDIDATE_SENSITIVITY_MERGE_METADATA"

    "$PARAMETER_SENSITIVITY"
    "$PARAMETER_SENSITIVITY_SUMMARY"
    "$PARAMETER_SENSITIVITY_CLASS_SUMMARY"
    "$PARAMETER_SENSITIVITY_METADATA"

    "$TIER1_SENSITIVITY_MATRIX_PNG"
    "$TIER1_SENSITIVITY_MATRIX_PDF"
    "$TIER1_ROBUSTNESS_RANKED_PNG"
    "$TIER1_ROBUSTNESS_RANKED_PDF"
    "$RECURRENCE_ROBUSTNESS_PNG"
    "$RECURRENCE_ROBUSTNESS_PDF"
)

for file in "${EXPECTED_OUTPUTS[@]}"; do
    [[ -s "$file" ]] || {
        log "ERROR: expected output missing or empty:"
        log "  $file"
        exit 1
    }
done


# ============================================================
# Finish
# ============================================================

log "============================================================"
log "SENSITIVITY ANALYSIS COMPLETE"
log "============================================================"
log "Finished: $(date)"
log ""
log "Main outputs:"
log "  $GLOBAL_SENSITIVITY_SUMMARY"
log "  $CANDIDATE_SENSITIVITY_SUMMARY"
log "  $CANDIDATES_WITH_SENSITIVITY"
log "  $PARAMETER_SENSITIVITY_SUMMARY"
log ""
log "Main figures:"
log "  $TIER1_SENSITIVITY_MATRIX_PNG"
log "  $TIER1_ROBUSTNESS_RANKED_PNG"
log "  $RECURRENCE_ROBUSTNESS_PNG"
log ""
log "Log:"
log "  $LOG_FILE"
