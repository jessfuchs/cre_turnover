#!/usr/bin/env bash

set -euo pipefail


# ============================================================
# Candidate-analysis pipeline
#
# Purpose:
#   Run the complete Tier-1 candidate-analysis workflow from
#   focal-clade detection to the final gene-level summary.
#
# Pipeline:
#   01  Focal-clade analysis
#   02  Secondary singleton Tier-1 analysis
#   03  Tier-1 candidate QC
#   04  CRE-level candidate prioritization
#   05  Candidate gene-distance annotation
#   06  Candidate gene assignment
#   07  CRE x FBgn expansion
#   08  Gene-level candidate summary
#
# Configuration:
#   All central input/output paths and analysis parameters are
#   supplied through config/candidate_config.sh.
# ============================================================


# ============================================================
# Pipeline root and configuration
# ============================================================

PIPELINE_ROOT="$(
    cd "$(dirname "${BASH_SOURCE[0]}")"
    pwd
)"

CONFIG_FILE="${PIPELINE_ROOT}/config/candidate_config.sh"


if [[ ! -f "${CONFIG_FILE}" ]]; then

    echo "ERROR: candidate-analysis config not found:"
    echo "${CONFIG_FILE}"

    exit 1
fi


# shellcheck source=/dev/null
source "${CONFIG_FILE}"


# ============================================================
# Runtime settings
# ============================================================

PYTHON_BIN="${PYTHON_BIN:-python3}"

LOG_DIR="${RESULTS_DIR}/logs"

mkdir -p \
    "${RESULTS_DIR}" \
    "${FOCAL_CLADES_RESULTS_DIR}" \
    "${SECONDARY_RESULTS_DIR}" \
    "${TIER1_QC_DIR}" \
    "${TIER1_TABLES_DIR}" \
    "${LOG_DIR}"


TIMESTAMP="$(
    date '+%Y%m%d_%H%M%S'
)"

LOG_FILE="${LOG_DIR}/candidate_analysis_pipeline_${TIMESTAMP}.log"


# ============================================================
# Pipeline scripts
# ============================================================

SCRIPT_01="${PIPELINE_ROOT}/01_analyze_focal_clades.py"
SCRIPT_02="${PIPELINE_ROOT}/02_analyze_secondary_tier1.py"
SCRIPT_03="${PIPELINE_ROOT}/03_qc_tier1_candidates.py"
SCRIPT_04="${PIPELINE_ROOT}/04_prioritize_tier1_candidates.py"
SCRIPT_05="${PIPELINE_ROOT}/05_annotate_candidate_gene_distances.py"
SCRIPT_06="${PIPELINE_ROOT}/06_assign_candidate_genes.py"
SCRIPT_07="${PIPELINE_ROOT}/07_explode_candidate_fbgns.py"
SCRIPT_08="${PIPELINE_ROOT}/08_build_gene_level_summary.py"


SCRIPTS=(
    "${SCRIPT_01}"
    "${SCRIPT_02}"
    "${SCRIPT_03}"
    "${SCRIPT_04}"
    "${SCRIPT_05}"
    "${SCRIPT_06}"
    "${SCRIPT_07}"
    "${SCRIPT_08}"
)


# ============================================================
# Helpers
# ============================================================

timestamp() {

    date '+%Y-%m-%d %H:%M:%S'
}


log() {

    echo "[$(timestamp)] $*" \
        | tee -a "${LOG_FILE}"
}


require_file() {

    local file="$1"
    local label="$2"

    if [[ ! -f "${file}" ]]; then

        log "ERROR: ${label} not found:"
        log "${file}"

        exit 1
    fi
}


run_step() {

    local step_number="$1"
    local step_name="$2"

    shift 2

    log "============================================================"
    log "Step ${step_number}/8: ${step_name}"
    log "============================================================"

    local start_seconds
    local end_seconds
    local runtime

    start_seconds="$(
        date +%s
    )"


    if "$@" 2>&1 | tee -a "${LOG_FILE}"; then

        end_seconds="$(
            date +%s
        )"

        runtime="$(
            (
                end_seconds
                - start_seconds
            )
        )"

        log ""
        log "PASS: ${step_name}"
        log "Runtime: ${runtime} s"
        log ""

    else

        local pipeline_status=(
            "${PIPESTATUS[@]}"
        )

        local command_status="${
            pipeline_status[0]
        }"

        local tee_status="${
            pipeline_status[1]:-0
        }"

        local status="${command_status}"

        if (
            ( status == 0 )
            &&
            ( tee_status != 0 )
        ); then

            status="${tee_status}"
        fi


        end_seconds="$(
            date +%s
        )"

        runtime="$(
            (
                end_seconds
                - start_seconds
            )
        )"


        log ""
        log "ERROR: ${step_name} failed"
        log "Exit status: ${status}"
        log "Runtime before failure: ${runtime} s"
        log ""
        log "Pipeline stopped."

        exit "${status}"
    fi
}


# ============================================================
# Start
# ============================================================

log "============================================================"
log "Candidate-analysis pipeline started"
log "============================================================"

log "Pipeline root: ${PIPELINE_ROOT}"
log "Configuration: ${CONFIG_FILE}"
log "Results directory: ${RESULTS_DIR}"
log "Log file: ${LOG_FILE}"

log ""


# ============================================================
# Python environment
# ============================================================

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then

    log "ERROR: Python executable not found:"
    log "${PYTHON_BIN}"

    exit 1
fi


log "Python executable: $(command -v "${PYTHON_BIN}")"
log "Python version: $("${PYTHON_BIN}" --version 2>&1)"


if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then

    log "Conda environment: ${CONDA_DEFAULT_ENV}"

else

    log "Conda environment: not detected"
fi


log ""


# ============================================================
# Validate pipeline scripts
# ============================================================

log "Checking pipeline scripts..."


for script in "${SCRIPTS[@]}"; do

    require_file \
        "${script}" \
        "pipeline script"

    log "FOUND: $(basename "${script}")"

done


log ""
log "All required pipeline scripts found."
log ""


# ============================================================
# Validate primary pipeline inputs
# ============================================================

log "Checking primary inputs..."


require_file \
    "${CRE_TURNOVER_MATRIX}" \
    "CRE turnover matrix"

require_file \
    "${FOCAL_CLADES_FILE}" \
    "focal-clade definition"

require_file \
    "${COMBINED_MANIFEST}" \
    "combined species manifest"

require_file \
    "${SPECIES_TRAITS}" \
    "species trait table"

require_file \
    "${REFERENCE_CRES_TSV}" \
    "D. melanogaster reference CRE table"

require_file \
    "${SO_ALL_SPECIES_FBGN}" \
    "combined SO/FBgn table"


log ""
log "All primary inputs found."
log ""


# ============================================================
# Step 01
# Focal-clade analysis
# ============================================================

run_step \
    "1" \
    "Focal-clade analysis" \
    "${PYTHON_BIN}" \
    "${SCRIPT_01}" \
    --matrix "${CRE_TURNOVER_MATRIX}" \
    --groups "${FOCAL_CLADES_FILE}" \
    --manifest "${COMBINED_MANIFEST}" \
    --traits "${SPECIES_TRAITS}" \
    --out-dir "${FOCAL_CLADES_RESULTS_DIR}" \
    --summary-out "${FOCAL_CLADE_SUMMARY}" \
    --reference-species "${REFERENCE_SPECIES}" \
    --expected-reference-cres "${EXPECTED_REFERENCE_CRES}"


# ============================================================
# Step 02
# Secondary singleton Tier-1 analysis
# ============================================================

run_step \
    "2" \
    "Secondary singleton Tier-1 analysis" \
    "${PYTHON_BIN}" \
    "${SCRIPT_02}" \
    --matrix "${CRE_TURNOVER_MATRIX}" \
    --groups "${FOCAL_CLADES_FILE}" \
    --outdir "${SECONDARY_RESULTS_DIR}" \
    --out-long "${SECONDARY_TIER1_ALL}" \
    --out-summary "${SECONDARY_TIER1_RECURRENT}" \
    --out-group-summary "${SECONDARY_CLADE_SUMMARY}" \
    --metadata-out "${SECONDARY_RUN_METADATA}" \
    --expected-reference-cres "${EXPECTED_REFERENCE_CRES}" \
    --recurrence-min-clades "${SECONDARY_RECURRENCE_MIN_CLADES}"


# ============================================================
# Step 03
# Tier-1 candidate QC
# ============================================================

run_step \
    "3" \
    "Tier-1 candidate QC" \
    "${PYTHON_BIN}" \
    "${SCRIPT_03}" \
    --focal-dir "${FOCAL_CLADES_RESULTS_DIR}" \
    --secondary "${SECONDARY_TIER1_ALL}" \
    --reference "${REFERENCE_CRES_TSV}" \
    --turnover-dir "${TURNOVER_BY_SPECIES_DIR}" \
    --out-summary "${TIER1_QC_SUMMARY}" \
    --out-detail "${TIER1_QC_DETAILS}" \
    --metadata-out "${TIER1_QC_METADATA}" \
    --min-reciprocal-overlap "${TIER1_QC_MIN_RECIPROCAL_OVERLAP}" \
    --expected-reference-cres "${EXPECTED_REFERENCE_CRES}" \
    --recurrence-min-clades "${SECONDARY_RECURRENCE_MIN_CLADES}"


# ============================================================
# Step 04
# CRE-level Tier-1 prioritization
# ============================================================

run_step \
    "4" \
    "Tier-1 CRE prioritization" \
    "${PYTHON_BIN}" \
    "${SCRIPT_04}" \
    --qc "${TIER1_QC_SUMMARY}" \
    --reference "${REFERENCE_CRES_TSV}" \
    --out-evidence "${TIER1_CLADE_EVIDENCE}" \
    --out-prioritized "${TIER1_PRIORITIZED}" \
    --metadata-out "${TIER1_PRIORITIZATION_METADATA}" \
    --focal-recurrent-min-clades "${FOCAL_RECURRENCE_MIN_CLADES}" \
    --secondary-recurrent-min-clades "${SECONDARY_RECURRENCE_MIN_CLADES}" \
    --expected-reference-cres "${EXPECTED_REFERENCE_CRES}"


# ============================================================
# Step 05
# Candidate gene-distance annotation
# ============================================================

run_step \
    "5" \
    "Candidate gene-distance annotation" \
    "${PYTHON_BIN}" \
    "${SCRIPT_05}" \
    --candidates "${TIER1_PRIORITIZED}" \
    --so-table "${SO_ALL_SPECIES_FBGN}" \
    --out "${TIER1_GENE_DISTANCES}" \
    --metadata-out "${TIER1_GENE_DISTANCE_METADATA}"


# ============================================================
# Step 06
# Candidate gene assignment
# ============================================================

run_step \
    "6" \
    "Candidate gene assignment" \
    "${PYTHON_BIN}" \
    "${SCRIPT_06}" \
    --input "${TIER1_GENE_DISTANCES}" \
    --out "${TIER1_GENE_ASSIGNMENTS}" \
    --metadata-out "${TIER1_GENE_ASSIGNMENT_METADATA}"


# ============================================================
# Step 07
# CRE x FBgn expansion
# ============================================================

run_step \
    "7" \
    "CRE x FBgn expansion" \
    "${PYTHON_BIN}" \
    "${SCRIPT_07}" \
    --input "${TIER1_GENE_ASSIGNMENTS}" \
    --out "${TIER1_FBGN_EXPLODED}" \
    --metadata-out "${TIER1_FBGN_EXPLODED_METADATA}"


# ============================================================
# Step 08
# Gene-level candidate summary
# ============================================================

run_step \
    "8" \
    "Gene-level candidate summary" \
    "${PYTHON_BIN}" \
    "${SCRIPT_08}" \
    --input "${TIER1_FBGN_EXPLODED}" \
    --out "${TIER1_GENE_LEVEL_SUMMARY}" \
    --metadata-out "${TIER1_GENE_LEVEL_METADATA}"


# ============================================================
# Final output checks
# ============================================================

log "============================================================"
log "Checking final outputs"
log "============================================================"


EXPECTED_OUTPUTS=(
    "${FOCAL_CLADE_SUMMARY}"

    "${SECONDARY_TIER1_ALL}"
    "${SECONDARY_TIER1_RECURRENT}"
    "${SECONDARY_CLADE_SUMMARY}"
    "${SECONDARY_RUN_METADATA}"

    "${TIER1_QC_SUMMARY}"
    "${TIER1_QC_DETAILS}"
    "${TIER1_QC_METADATA}"

    "${TIER1_CLADE_EVIDENCE}"
    "${TIER1_PRIORITIZED}"
    "${TIER1_PRIORITIZATION_METADATA}"

    "${TIER1_GENE_DISTANCES}"
    "${TIER1_GENE_DISTANCE_METADATA}"

    "${TIER1_GENE_ASSIGNMENTS}"
    "${TIER1_GENE_ASSIGNMENT_METADATA}"

    "${TIER1_FBGN_EXPLODED}"
    "${TIER1_FBGN_EXPLODED_METADATA}"

    "${TIER1_GENE_LEVEL_SUMMARY}"
    "${TIER1_GENE_LEVEL_METADATA}"
)


N_MISSING=0


for outfile in "${EXPECTED_OUTPUTS[@]}"; do

    if [[ -s "${outfile}" ]]; then

        log "FOUND: ${outfile}"

    else

        log "MISSING/EMPTY: ${outfile}"

        N_MISSING=$(
            (
                N_MISSING
                + 1
            )
        )
    fi

done


log ""


if [[ "${N_MISSING}" -gt 0 ]]; then

    log "ERROR: pipeline steps finished, but ${N_MISSING} expected output(s) are missing or empty."

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
log "${TIER1_GENE_LEVEL_SUMMARY}"

log ""
log "Prioritized CRE-level candidate table:"
log "${TIER1_PRIORITIZED}"

log ""
log "Full pipeline log:"
log "${LOG_FILE}"
