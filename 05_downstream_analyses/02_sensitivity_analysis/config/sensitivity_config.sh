#!/usr/bin/env bash

# ======================================================================
# Sensitivity-analysis configuration
# ======================================================================

PIPELINE_ROOT="$(
    cd "$(dirname "${BASH_SOURCE[0]}")/.."
    pwd
)"

DOWNSTREAM_ROOT="$(
    cd "${PIPELINE_ROOT}/.."
    pwd
)"

PROJECT_ROOT="$(
    cd "${DOWNSTREAM_ROOT}/.."
    pwd
)"


# ======================================================================
# Upstream pipeline roots
# ======================================================================

CLASSIFICATION_ROOT="${PROJECT_ROOT}/04_cre_classification"
MAPPING_ROOT="${PROJECT_ROOT}/02_mapping_orthologs"
WGA_ROOT="${PROJECT_ROOT}/03_pairwise_wga"


# ======================================================================
# Input
# ======================================================================

CLASSIFIER_SCRIPT="${CLASSIFICATION_ROOT}/scripts/01_classify_target_cres.py"
TARGET_SPECIES_FILE="${WGA_ROOT}/target_species.txt"
REFERENCE_CRES_TSV="${MAPPING_ROOT}/reference_cres/dmel_reference_cres.tsv"
SO_ALL_SPECIES_FBGN="${MAPPING_ROOT}/results/SO_all_species_fbgn.tsv"
LIFTED_CRES_DIR="${WGA_ROOT}/lifted_cres_dmel"


# ======================================================================
# Output
# ======================================================================

RESULTS_DIR="${PIPELINE_ROOT}/results"
SENSITIVITY_SCENARIOS_DIR="${RESULTS_DIR}/scenarios"
SENSITIVITY_SCENARIO_MANIFEST="${RESULTS_DIR}/sensitivity_scenarios.tsv"
SENSITIVITY_GENERATION_METADATA="${RESULTS_DIR}/sensitivity_generation_metadata.tsv"


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

EXPECTED_REFERENCE_CRES="${EXPECTED_REFERENCE_CRES:-337}"


# ======================================================================
# Primary-scenario regression check
# ======================================================================

# Baseline CRE classifications from the primary analysis.
BASELINE_TURNOVER_DIR="${CLASSIFICATION_ROOT}/results/turnover_by_species"

# Sensitivity scenario expected to reproduce the baseline analysis.
PRIMARY_SENSITIVITY_SCENARIO="${PRIMARY_SENSITIVITY_SCENARIO:-ov050_dist24000}"

# Regression-check outputs.
PRIMARY_REGRESSION_SUMMARY="${RESULTS_DIR}/primary_regression_check.tsv"

PRIMARY_REGRESSION_CHANGES="${RESULTS_DIR}/primary_regression_state_changes.tsv"

PRIMARY_REGRESSION_METADATA="${RESULTS_DIR}/primary_regression_metadata.tsv"


# ======================================================================
# Export
# ======================================================================

export PIPELINE_ROOT
export PROJECT_ROOT

export CLASSIFICATION_ROOT
export MAPPING_ROOT
export WGA_ROOT

export CLASSIFIER_SCRIPT
export TARGET_SPECIES_FILE
export REFERENCE_CRES_TSV
export SO_ALL_SPECIES_FBGN
export LIFTED_CRES_DIR

export RESULTS_DIR
export SENSITIVITY_SCENARIOS_DIR
export SENSITIVITY_SCENARIO_MANIFEST
export SENSITIVITY_GENERATION_METADATA

export EXPECTED_REFERENCE_CRES

export BASELINE_TURNOVER_DIR
export PRIMARY_SENSITIVITY_SCENARIO

export PRIMARY_REGRESSION_SUMMARY
export PRIMARY_REGRESSION_CHANGES
export PRIMARY_REGRESSION_METADATA
