#!/usr/bin/env bash
set -euo pipefail


# ======================================================================
# Pipeline roots
# ======================================================================

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_ROOT="$(cd "$PIPELINE_ROOT/.." && pwd)"

MAPPING_ROOT="${PROJECT_ROOT}/02_mapping_orthologs"
WGA_ROOT="${PROJECT_ROOT}/03_pairwise_wga"


# ======================================================================
# Input data
# ======================================================================

REFERENCE_CRES_TSV="${MAPPING_ROOT}/reference_cres/dmel_reference_cres.tsv"

SO_ALL_SPECIES_FBGN="${MAPPING_ROOT}/ortholog_results/SO_all_species_fbgn.tsv"

LIFTED_CRES_DIR="${WGA_ROOT}/lifted_cres_dmel"

TARGET_SPECIES_FILE="${WGA_ROOT}/target_species.txt"


# ======================================================================
# CRE classification
# ======================================================================

RESULTS_DIR="${PIPELINE_ROOT}/results"
SENSITIVITY_DIR="${PIPELINE_ROOT}/sensitivity"

TURNOVER_BY_SPECIES_DIR="${RESULTS_DIR}/turnover_by_species"

CRE_TURNOVER_ALL_SPECIES="${RESULTS_DIR}/cre_turnover_all_species.tsv"
CRE_TURNOVER_MATRIX="${RESULTS_DIR}/cre_turnover_matrix.tsv"


# ======================================================================
# Classification parameters
# ======================================================================

RECIPROCAL_OVERLAP="${RECIPROCAL_OVERLAP:-0.50}"

LOCAL_GENE_DISTANCE="${LOCAL_GENE_DISTANCE:-24000}"


# ======================================================================
# Phylogeny
# ======================================================================

PHYLOGENY_ROOT="${PIPELINE_ROOT}/phylogeny"

PHYLOGENY_DATA_DIR="${PHYLOGENY_ROOT}/data"
PHYLOGENY_RESULTS_DIR="${PHYLOGENY_ROOT}/results"
PHYLOGENY_SCRIPTS_DIR="${PHYLOGENY_ROOT}/scripts"

# Adjust filenames to the files actually present in your project.
TREE_FILE="${PHYLOGENY_RESULTS_DIR}/301Fly_HOG_UCLDtree_40species.nw"

SPECIES_ORDER_FILE="${PHYLOGENY_RESULTS_DIR}/species_order_40_tree_names.txt"

SPECIES_METADATA="${PHYLOGENY_DATA_DIR}/species_metadata.tsv"


# ======================================================================
# Export
# ======================================================================

export PIPELINE_ROOT
export PROJECT_ROOT
export MAPPING_ROOT
export WGA_ROOT

export REFERENCE_CRES_TSV
export SO_ALL_SPECIES_FBGN
export LIFTED_CRES_DIR
export TARGET_SPECIES_FILE

export RESULTS_DIR
export SENSITIVITY_DIR
export TURNOVER_BY_SPECIES_DIR
export CRE_TURNOVER_ALL_SPECIES
export CRE_TURNOVER_MATRIX

export RECIPROCAL_OVERLAP
export LOCAL_GENE_DISTANCE

export PHYLOGENY_ROOT
export PHYLOGENY_DATA_DIR
export PHYLOGENY_RESULTS_DIR
export PHYLOGENY_SCRIPTS_DIR
export TREE_FILE
export SPECIES_ORDER_FILE
export SPECIES_METADATA
