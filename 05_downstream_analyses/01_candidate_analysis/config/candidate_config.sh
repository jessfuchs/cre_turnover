#!/usr/bin/env bash

# ======================================================================
# Candidate-analysis configuration
# ======================================================================

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOWNSTREAM_ROOT="$(cd "$PIPELINE_ROOT/.." && pwd)"
PROJECT_ROOT="$(cd "$DOWNSTREAM_ROOT/.." && pwd)"

CLASSIFICATION_ROOT="${PROJECT_ROOT}/04_cre_classification"
SCRM_EXTERNAL_ROOT="${PROJECT_ROOT}/01_scrmshaw/external"

CONFIG_DIR="${PIPELINE_ROOT}/config"
RESULTS_DIR="${PIPELINE_ROOT}/results"


# ======================================================================
# Input
# ======================================================================

# Final CRE-state matrix.
CRE_TURNOVER_MATRIX="${CLASSIFICATION_ROOT}/results/cre_turnover_matrix.tsv"

# Combined species manifest.
COMBINED_MANIFEST="${SCRM_EXTERNAL_ROOT}/combined_manifest.tsv"

# Species-level trait annotations, including climatic zone.
SPECIES_TRAITS="${CLASSIFICATION_ROOT}/phylogeny/data/species_traits.tsv"

# Focal-clade definitions.
FOCAL_CLADES_FILE="${CONFIG_DIR}/focal_clades.tsv"


# ======================================================================
# Output
# ======================================================================

FOCAL_CLADES_RESULTS_DIR="${RESULTS_DIR}/focal_clades"

FOCAL_CLADE_SUMMARY="${RESULTS_DIR}/focal_clade_summary.tsv"


# ======================================================================
# Analysis parameters
# ======================================================================

REFERENCE_SPECIES="${REFERENCE_SPECIES:-d_melanogaster}"

EXPECTED_REFERENCE_CRES="${EXPECTED_REFERENCE_CRES:-337}"


# ======================================================================
# Export
# ======================================================================

export PIPELINE_ROOT
export PROJECT_ROOT

export CRE_TURNOVER_MATRIX
export COMBINED_MANIFEST
export SPECIES_TRAITS
export FOCAL_CLADES_FILE

export FOCAL_CLADES_RESULTS_DIR
export FOCAL_CLADE_SUMMARY

export REFERENCE_SPECIES
export EXPECTED_REFERENCE_CRES
