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

# ======================================================================
# Phylogenetically controlled climate analysis
# ======================================================================

PHYLOGENY_TREE="${CLASSIFICATION_ROOT}/phylogeny/results/301Fly_HOG_UCLDtree_40species.nw"
PGLS_CLIMATE_MODEL_COMPARISON="${RESULTS_DIR}/pgls_climate_model_comparison.tsv"
PGLS_CLIMATE_COEFFICIENTS="${RESULTS_DIR}/pgls_climate_coefficients.tsv"
PGLS_CLIMATE_SUMMARY="${RESULTS_DIR}/pgls_climate_summary.txt"

# ======================================================================
# Focal Tier-1 enrichment
# ======================================================================

CANDIDATE_ANALYSIS_ROOT="${DOWNSTREAM_ROOT}/01_candidate_analysis"
CRE_TURNOVER_MATRIX="${CLASSIFICATION_ROOT}/results/cre_turnover_matrix.tsv"
FOCAL_CLADES_FILE="${CANDIDATE_ANALYSIS_ROOT}/config/focal_clades.tsv"
FOCAL_TIER1_CANDIDATES="${RESULTS_DIR}/focal_tier1_candidates.tsv"
FOCAL_TIER1_ENRICHMENT="${RESULTS_DIR}/focal_tier1_enrichment.tsv"

# ======================================================================
# Azteca focal candidates
# ======================================================================

CANDIDATE_ANALYSIS_ROOT="${DOWNSTREAM_ROOT}/01_candidate_analysis"
CANDIDATES_WITH_SENSITIVITY="${CANDIDATE_ANALYSIS_ROOT}/results/tables/tier1_candidates_prioritized_with_sensitivity.tsv"
AZTECA_FOCAL_CLADE="${AZTECA_FOCAL_CLADE:-azteca_affinis_miranda_group}"
AZTECA_FOCAL_CANDIDATES="${RESULTS_DIR}/azteca_focal_tier1_candidates.tsv"

# ======================================================================
# Figures
# ======================================================================

FOCAL_TIER1_ENRICHMENT_PNG="${FIGURES_DIR}/focal_tier1_enrichment.png"
FOCAL_TIER1_ENRICHMENT_PDF="${FIGURES_DIR}/focal_tier1_enrichment.pdf"
