#!/usr/bin/env bash

# ======================================================================
# Sensitivity-analysis configuration
# ======================================================================

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOWNSTREAM_ROOT="$(cd "${PIPELINE_ROOT}/.." && pwd)"
PROJECT_ROOT="$(cd "${DOWNSTREAM_ROOT}/.." && pwd)"

CLASSIFICATION_ROOT="${PROJECT_ROOT}/04_cre_classification"
MAPPING_ROOT="${PROJECT_ROOT}/02_mapping_orthologs"
WGA_ROOT="${PROJECT_ROOT}/03_pairwise_wga"

CANDIDATE_ANALYSIS_ROOT="${DOWNSTREAM_ROOT}/01_candidate_analysis"


# ======================================================================
# Inputs
# ======================================================================

CLASSIFIER_SCRIPT="${CLASSIFICATION_ROOT}/scripts/01_classify_target_cres.py"
TARGET_SPECIES_FILE="${WGA_ROOT}/target_species.txt"
REFERENCE_CRES_TSV="${MAPPING_ROOT}/reference_cres/dmel_reference_cres.tsv"
SO_ALL_SPECIES_FBGN="${MAPPING_ROOT}/results/SO_all_species_fbgn.tsv"
LIFTED_CRES_DIR="${WGA_ROOT}/lifted_cres_dmel"
BASELINE_TURNOVER_DIR="${CLASSIFICATION_ROOT}/results/turnover_by_species"

CANDIDATE_TIER1_TABLE="${CANDIDATE_ANALYSIS_ROOT}/results/tables/tier1_candidates_prioritized.tsv"
FOCAL_CLADES_FILE="${CANDIDATE_ANALYSIS_ROOT}/config/focal_clades.tsv"

# ======================================================================
# Sensitivity parameters
# ======================================================================

SENSITIVITY_OVERLAPS=(
    0.25
    0.50
    0.75
)

SENSITIVITY_DISTANCES=(
    12000
    24000
    48000
)

PRIMARY_SENSITIVITY_SCENARIO="${PRIMARY_SENSITIVITY_SCENARIO:-ov050_dist24000}"
EXPECTED_REFERENCE_CRES="${EXPECTED_REFERENCE_CRES:-337}"

# Candidate-priority recurrence thresholds.
# These must match the thresholds used in candidate_analysis.
FOCAL_RECURRENCE_MIN_CLADES="${FOCAL_RECURRENCE_MIN_CLADES:-2}"
SECONDARY_RECURRENCE_MIN_CLADES="${SECONDARY_RECURRENCE_MIN_CLADES:-2}"

# Minimum retention percentage classified as moderate robustness.
CANDIDATE_MODERATE_ROBUSTNESS_MIN="${CANDIDATE_MODERATE_ROBUSTNESS_MIN:-66.6}"

# ======================================================================
# Results
# ======================================================================

RESULTS_DIR="${PIPELINE_ROOT}/results"
SENSITIVITY_SCENARIOS_DIR="${RESULTS_DIR}/scenarios"
SENSITIVITY_SCENARIO_MANIFEST="${RESULTS_DIR}/sensitivity_scenarios.tsv"
SENSITIVITY_GENERATION_METADATA="${RESULTS_DIR}/sensitivity_generation_metadata.tsv"

# Primary regression

PRIMARY_REGRESSION_SUMMARY="${RESULTS_DIR}/primary_regression_check.tsv"
PRIMARY_REGRESSION_CHANGES="${RESULTS_DIR}/primary_regression_state_changes.tsv"
PRIMARY_REGRESSION_METADATA="${RESULTS_DIR}/primary_regression_metadata.tsv"

# Global sensitivity

GLOBAL_SENSITIVITY_SUMMARY="${RESULTS_DIR}/sensitivity_global_summary.tsv"
SENSITIVITY_SPECIES_STABILITY="${RESULTS_DIR}/sensitivity_species_stability.tsv"
SENSITIVITY_STATE_TRANSITIONS="${RESULTS_DIR}/sensitivity_state_transitions.tsv"
SENSITIVITY_CRE_STABILITY="${RESULTS_DIR}/sensitivity_cre_stability.tsv"
GLOBAL_SENSITIVITY_METADATA="${RESULTS_DIR}/sensitivity_global_summary_metadata.tsv"

# Candidate sensitivity

CANDIDATE_SENSITIVITY_LONG="${RESULTS_DIR}/candidate_sensitivity_long.tsv"
CANDIDATE_SENSITIVITY_SUMMARY="${RESULTS_DIR}/candidate_sensitivity_summary.tsv"
CANDIDATE_SENSITIVITY_PRIORITY_MATRIX="${RESULTS_DIR}/candidate_sensitivity_priority_matrix.tsv"
CANDIDATE_SENSITIVITY_RETENTION_MATRIX="${RESULTS_DIR}/candidate_sensitivity_retention_matrix.tsv"
CANDIDATE_SENSITIVITY_METADATA="${RESULTS_DIR}/candidate_sensitivity_metadata.tsv"
CANDIDATE_SENSITIVITY_MERGE_METADATA="${RESULTS_DIR}/candidate_sensitivity_merge_metadata.tsv"

CANDIDATES_WITH_SENSITIVITY="${CANDIDATE_ANALYSIS_ROOT}/results/tables/tier1_candidates_prioritized_with_sensitivity.tsv"

# ======================================================================
# Figures
# ======================================================================

FIGURES_DIR="${RESULTS_DIR}/figures"

TIER1_ROBUSTNESS_RANKED_PNG="${FIGURES_DIR}/tier1_candidate_robustness_ranked.png"
TIER1_ROBUSTNESS_RANKED_PDF="${FIGURES_DIR}/tier1_candidate_robustness_ranked.pdf"

