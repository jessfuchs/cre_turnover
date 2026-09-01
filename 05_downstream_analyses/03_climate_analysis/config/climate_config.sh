#!/usr/bin/env bash

# ======================================================================
# Climate-analysis configuration
# ======================================================================

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOWNSTREAM_ROOT="$(cd "${PIPELINE_ROOT}/.." && pwd)"
PROJECT_ROOT="$(cd "${DOWNSTREAM_ROOT}/.." && pwd)"

CLASSIFICATION_ROOT="${PROJECT_ROOT}/04_cre_classification"
SCRM_EXTERNAL_ROOT="${PROJECT_ROOT}/01_scrmshaw/external"

RESULTS_DIR="${PIPELINE_ROOT}/results"
FIGURES_DIR="${RESULTS_DIR}/figures"


# ======================================================================
# Inputs
# ======================================================================

SPECIES_SUMMARY="${CLASSIFICATION_ROOT}/results/species_summary.tsv"
COMBINED_MANIFEST="${SCRM_EXTERNAL_ROOT}/combined_manifest.tsv"
SPECIES_TRAITS="${CLASSIFICATION_ROOT}/phylogeny/data/species_traits.tsv"


# ======================================================================
# Species-level turnover by climatic zone
# ======================================================================

CLIMATE_TURNOVER_SUMMARY="${RESULTS_DIR}/climate_turnover_summary.tsv"
CLIMATE_TURNOVER_DESCRIPTIVE="${RESULTS_DIR}/climate_turnover_descriptive_statistics.tsv"
CLIMATE_TURNOVER_KRUSKAL="${RESULTS_DIR}/climate_turnover_kruskal.tsv"

TURNOVER_BY_CLIMATE_PNG="${FIGURES_DIR}/turnover_rate_by_climate.png"
TURNOVER_BY_CLIMATE_PDF="${FIGURES_DIR}/turnover_rate_by_climate.pdf"
