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

#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Run the complete sensitivity-analysis workflow
# ============================================================

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${PIPELINE_ROOT}/config/sensitivity_config.sh"
[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: missing config: $CONFIG_FILE" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG_FILE"

LOG_DIR="${RESULTS_DIR}/logs"
mkdir -p "$RESULTS_DIR" "$FIGURES_DIR" "$LOG_DIR"
LOG_FILE="${LOG_DIR}/sensitivity_analysis_$(date '+%Y%m%d_%H%M%S').log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

run_step() {
    local step="$1" total="$2" label="$3"; shift 3
    log "============================================================"
    log "Step ${step}/${total}: ${label}"
    log "============================================================"
    local start status tee_status
    start="$(date +%s)"
    set +e
    "$@" 2>&1 | tee -a "$LOG_FILE"
    status="${PIPESTATUS[0]}"; tee_status="${PIPESTATUS[1]}"
    set -e
    [[ "$status" -eq 0 ]] || { log "ERROR: ${label} failed (exit ${status})."; exit "$status"; }
    [[ "$tee_status" -eq 0 ]] || { log "ERROR: logging failed during ${label}."; exit "$tee_status"; }
    log "PASS: ${label} ($(( $(date +%s) - start )) s)"
    log ""
}

SCRIPT_01="${PIPELINE_ROOT}/01_generate_sensitivity_scenarios.py"
SCRIPT_02="${PIPELINE_ROOT}/02_check_primary_regression.py"
SCRIPT_03="${PIPELINE_ROOT}/03_summarize_global_sensitivity.py"
SCRIPT_04="${PIPELINE_ROOT}/04_summarize_candidate_sensitivity.py"
SCRIPT_05="${PIPELINE_ROOT}/05_add_sensitivity_to_prioritized_candidates.py"
SCRIPT_06="${PIPELINE_ROOT}/06_summarize_parameter_sensitivity.py"
PLOT_MATRIX="${PIPELINE_ROOT}/07_plot_sensitivity_matrix.py"
PLOT_RANKED="${PIPELINE_ROOT}/08_plot_ranked_tier1.py"
PLOT_RECURRENCE="${PIPELINE_ROOT}/09_plot_recurrence_robustness.py"

required=(
    "$SCRIPT_01" "$SCRIPT_02" "$SCRIPT_03" "$SCRIPT_04" "$SCRIPT_05" "$SCRIPT_06"
    "$PLOT_MATRIX" "$PLOT_RANKED" "$PLOT_RECURRENCE"
    "$CLASSIFIER_SCRIPT" "$TARGET_SPECIES_FILE" "$REFERENCE_CRES_TSV"
    "$SO_ALL_SPECIES_FBGN" "$CANDIDATE_TIER1_TABLE" "$CANDIDATE_TIER1_METADATA" "$FOCAL_CLADES_FILE"
)
for path in "${required[@]}"; do
    [[ -f "$path" ]] || { echo "ERROR: required file not found: $path" >&2; exit 1; }
done
[[ -d "$LIFTED_CRES_DIR" ]] || { echo "ERROR: lifted-CRE directory not found: $LIFTED_CRES_DIR" >&2; exit 1; }
[[ -d "$BASELINE_TURNOVER_DIR" ]] || { echo "ERROR: baseline turnover directory not found: $BASELINE_TURNOVER_DIR" >&2; exit 1; }

log "Sensitivity-analysis pipeline started"
log "Pipeline root: $PIPELINE_ROOT"
log "Primary scenario: $PRIMARY_SENSITIVITY_SCENARIO"
log "Overlap grid: ${SENSITIVITY_OVERLAPS[*]}"
log "Distance grid: ${SENSITIVITY_DISTANCES[*]}"
log ""

run_step 1 9 "Generate sensitivity scenarios" \
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

run_step 2 9 "Check primary regression" \
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

run_step 3 9 "Summarize global sensitivity" \
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

run_step 4 9 "Summarize candidate sensitivity" \
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

run_step 5 9 "Add sensitivity to prioritized candidates" \
    python3 "$SCRIPT_05" \
        --candidates "$CANDIDATE_TIER1_TABLE" \
        --sensitivity "$CANDIDATE_SENSITIVITY_SUMMARY" \
        --out "$CANDIDATES_WITH_SENSITIVITY" \
        --metadata-out "$CANDIDATE_SENSITIVITY_MERGE_METADATA" \
        --moderate-robustness-min "$CANDIDATE_MODERATE_ROBUSTNESS_MIN"

run_step 6 9 "Summarize parameter-specific sensitivity" \
    python3 "$SCRIPT_06" \
        --retention-matrix "$CANDIDATE_SENSITIVITY_RETENTION_MATRIX" \
        --global-summary "$GLOBAL_SENSITIVITY_SUMMARY" \
        --scenario-manifest "$SENSITIVITY_SCENARIO_MANIFEST" \
        --primary-scenario "$PRIMARY_SENSITIVITY_SCENARIO" \
        --out-candidates "$PARAMETER_SENSITIVITY" \
        --out-summary "$PARAMETER_SENSITIVITY_SUMMARY" \
        --out-classes "$PARAMETER_SENSITIVITY_CLASS_SUMMARY" \
        --metadata-out "$PARAMETER_SENSITIVITY_METADATA"

run_step 7 9 "Plot candidate sensitivity matrix" \
    python3 "$PLOT_MATRIX" \
        --matrix "$CANDIDATE_SENSITIVITY_RETENTION_MATRIX" \
        --summary "$CANDIDATE_SENSITIVITY_SUMMARY" \
        --scenario-manifest "$SENSITIVITY_SCENARIO_MANIFEST" \
        --primary-scenario "$PRIMARY_SENSITIVITY_SCENARIO" \
        --out-png "$TIER1_SENSITIVITY_MATRIX_PNG" \
        --out-pdf "$TIER1_SENSITIVITY_MATRIX_PDF"

run_step 8 9 "Plot ranked Tier-1 robustness" \
    python3 "$PLOT_RANKED" \
        --input "$CANDIDATES_WITH_SENSITIVITY" \
        --out-png "$TIER1_ROBUSTNESS_RANKED_PNG" \
        --out-pdf "$TIER1_ROBUSTNESS_RANKED_PDF" \
        --moderate-robustness-min "$CANDIDATE_MODERATE_ROBUSTNESS_MIN"

run_step 9 9 "Plot recurrence versus robustness" \
    python3 "$PLOT_RECURRENCE" \
        --input "$CANDIDATES_WITH_SENSITIVITY" \
        --out-png "$RECURRENCE_ROBUSTNESS_PNG" \
        --out-pdf "$RECURRENCE_ROBUSTNESS_PDF" \
        --moderate-robustness-min "$CANDIDATE_MODERATE_ROBUSTNESS_MIN"

expected=(
    "$SENSITIVITY_SCENARIO_MANIFEST" "$SENSITIVITY_GENERATION_METADATA"
    "$PRIMARY_REGRESSION_SUMMARY" "$PRIMARY_REGRESSION_CHANGES" "$PRIMARY_REGRESSION_METADATA"
    "$GLOBAL_SENSITIVITY_SUMMARY" "$SENSITIVITY_SPECIES_STABILITY" "$SENSITIVITY_STATE_TRANSITIONS"
    "$SENSITIVITY_CRE_STABILITY" "$GLOBAL_SENSITIVITY_METADATA"
    "$CANDIDATE_SENSITIVITY_LONG" "$CANDIDATE_SENSITIVITY_SUMMARY"
    "$CANDIDATE_SENSITIVITY_PRIORITY_MATRIX" "$CANDIDATE_SENSITIVITY_RETENTION_MATRIX"
    "$CANDIDATE_SENSITIVITY_METADATA" "$CANDIDATES_WITH_SENSITIVITY" "$CANDIDATE_SENSITIVITY_MERGE_METADATA"
    "$PARAMETER_SENSITIVITY" "$PARAMETER_SENSITIVITY_SUMMARY" "$PARAMETER_SENSITIVITY_CLASS_SUMMARY"
    "$PARAMETER_SENSITIVITY_METADATA" "$TIER1_SENSITIVITY_MATRIX_PNG" "$TIER1_SENSITIVITY_MATRIX_PDF"
    "$TIER1_ROBUSTNESS_RANKED_PNG" "$TIER1_ROBUSTNESS_RANKED_PDF"
    "$RECURRENCE_ROBUSTNESS_PNG" "$RECURRENCE_ROBUSTNESS_PDF"
)

missing=0
for path in "${expected[@]}"; do
    if [[ ! -s "$path" ]]; then
        log "MISSING/EMPTY: $path"
        missing=$((missing + 1))
    fi
done
[[ "$missing" -eq 0 ]] || { log "ERROR: ${missing} expected output(s) missing."; exit 1; }

log "Sensitivity-analysis pipeline completed successfully"
log "Global summary: $GLOBAL_SENSITIVITY_SUMMARY"
log "Candidate summary: $CANDIDATE_SENSITIVITY_SUMMARY"
log "Annotated candidates: $CANDIDATES_WITH_SENSITIVITY"
log "Log: $LOG_FILE"
