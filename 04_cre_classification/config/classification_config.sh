#!/usr/bin/env bash
set -euo pipefail

# ======================================================================
# Pipeline roots
# ======================================================================

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_ROOT="$(cd "$PIPELINE_ROOT/.." && pwd)"

SCRIPTS_DIR="${PIPELINE_ROOT}/scripts"

MAPPING_ROOT="${PROJECT_ROOT}/02_mapping_orthologs"
WGA_ROOT="${PROJECT_ROOT}/03_pairwise_wga"

# ======================================================================
# Input data
# ======================================================================

REFERENCE_CRES_TSV="${MAPPING_ROOT}/reference_cres/dmel_reference_cres.tsv"

SO_ALL_SPECIES_FBGN="${MAPPING_ROOT}/results/SO_all_species_fbgn.tsv"

LIFTED_CRES_DIR="${WGA_ROOT}/lifted_cres_dmel"

TARGET_SPECIES_FILE="${WGA_ROOT}/target_species.txt"

# ======================================================================
# CRE classification
# ======================================================================

RESULTS_DIR="${PIPELINE_ROOT}/results"
SENSITIVITY_DIR="${PIPELINE_ROOT}/sensitivity"

# ======================================================================
# Classification parameters
# ======================================================================

RECIPROCAL_OVERLAP="${RECIPROCAL_OVERLAP:-0.50}"

LOCAL_GENE_DISTANCE="${LOCAL_GENE_DISTANCE:-24000}"

# ======================================================================
# Classification outputs
# ======================================================================

TURNOVER_BY_SPECIES_DIR="${RESULTS_DIR}/turnover_by_species"

CRE_TURNOVER_ALL_SPECIES="${RESULTS_DIR}/cre_turnover_all_species.tsv"
CRE_TURNOVER_MATRIX="${RESULTS_DIR}/cre_turnover_matrix.tsv"
SPECIES_SUMMARY="${RESULTS_DIR}/species_summary.tsv"

# ======================================================================
# Expected dataset dimensions
# ======================================================================

EXPECTED_REFERENCE_CRES="${EXPECTED_REFERENCE_CRES:-337}"

# ======================================================================
# Phylogeny
# ======================================================================

PHYLOGENY_ROOT="${PIPELINE_ROOT}/phylogeny"

PHYLOGENY_DATA_DIR="${PHYLOGENY_ROOT}/data"
PHYLOGENY_RESULTS_DIR="${PHYLOGENY_ROOT}/results"
PHYLOGENY_SCRIPTS_DIR="${PHYLOGENY_ROOT}/scripts"

# Adjust filenames to the files actually positional_match in your project.
TREE_FILE="${PHYLOGENY_RESULTS_DIR}/301Fly_HOG_UCLDtree_40species.nw"

# Species-level trait annotations, e.g. climate zone.
SPECIES_TRAITS="${PHYLOGENY_ROOT}/data/species_traits.tsv"

SPECIES_ORDER_FILE="${PHYLOGENY_RESULTS_DIR}/species_order_40_tree_names.txt"

# ======================================================================
# Species-level QC parameters
# ======================================================================

QC_MIN_MAPPING_RATE="${QC_MIN_MAPPING_RATE:-0.50}"
QC_MIN_SCRMSHAW_PEAKS="${QC_MIN_SCRMSHAW_PEAKS:-50}"
QC_MIN_SEQID_OVERLAP="${QC_MIN_SEQID_OVERLAP:-0.25}"
QC_MIN_MAPPED_FOR_OVERLAP_CHECK="${QC_MIN_MAPPED_FOR_OVERLAP_CHECK:-100}"
QC_HIGH_MAPPING_RATE="${QC_HIGH_MAPPING_RATE:-0.90}"

SPECIES_QC_SUMMARY="${RESULTS_DIR}/species_qc_summary.tsv"

# ======================================================================
# Analysis outputs
# ======================================================================

FIGURES_DIR="${RESULTS_DIR}/figures"
FOCAL_CLADE_CRE_TABLE="${RESULTS_DIR}/focal_clade_heatmap_CREs.tsv"

CRE_CONSERVATION_SUMMARY="${RESULTS_DIR}/CRE_conservation_summary.tsv"
SPECIES_PLOT_QC_SUMMARY="${RESULTS_DIR}/species_plot_qc_summary.tsv"

# D. melanogaster CRE-to-gene distance summary.
DMEL_GENE_DISTANCE_SUMMARY="${RESULTS_DIR}/dmel_cre_gene_distance_summary.tsv"

