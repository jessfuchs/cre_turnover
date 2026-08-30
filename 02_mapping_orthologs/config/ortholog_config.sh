#!/usr/bin/env bash
set -euo pipefail

# ======================================================================
# Pipeline roots
# ======================================================================

# Root directory of the ortholog-mapping pipeline: ~/cre_turnover/02_mapping_orthologs
PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Root directory of the complete CRE turnover project: ~/cre_turnover
PROJECT_ROOT="$(cd "$PIPELINE_ROOT/.." && pwd)"

# ======================================================================
# Pipeline directories
# ======================================================================

SCRIPTS_DIR="${PIPELINE_ROOT}/scripts"

ORTHOLOG_RESULTS_DIR="${PIPELINE_ROOT}/ortholog_results"
REFERENCE_CRES_DIR="${PIPELINE_ROOT}/reference_cres"

# ======================================================================
# Upstream SCRMshaw data
# ======================================================================

SCRM_GENERATED_ROOT="${PROJECT_ROOT}/01_scrmshaw/generation"
SCRM_EXTERNAL_ROOT="${PROJECT_ROOT}/01_scrmshaw/external"

COMBINED_MANIFEST="${SCRM_EXTERNAL_ROOT}/combined_manifest.tsv"

# Species-specific GFF3 files used for ortholog-annotation QC.
GENERATED_GFF_ROOT="${SCRM_GENERATED_ROOT}/data/selected"
EXTERNAL_GFF_ROOT="${SCRM_EXTERNAL_ROOT}/external_data/selected"

# D. melanogaster annotation used later for FBgn normalization.
DMEL_GFF="${EXTERNAL_GFF_ROOT}/d_melanogaster/annotation.gff3"

# ======================================================================
# Main ortholog-mapping outputs
# ======================================================================

SO_ALL_SPECIES="${ORTHOLOG_RESULTS_DIR}/SO_all_species.tsv"
SO_ALL_SPECIES_FBGN="${ORTHOLOG_RESULTS_DIR}/SO_all_species_fbgn.tsv"

UNRESOLVED_DMEL="${ORTHOLOG_RESULTS_DIR}/unresolved_dmel_identifiers.tsv"
UNRESOLVED_QC_SUMMARY="${ORTHOLOG_RESULTS_DIR}/unresolved_dmel_qc_summary.tsv"

DMEL_ORTHOLOG_QC="${ORTHOLOG_RESULTS_DIR}/dmel_ortholog_annotation_qc.tsv"

# ======================================================================
# Reference CRE outputs
# ======================================================================

REFERENCE_CRES_TSV="${REFERENCE_CRES_DIR}/dmel_reference_cres.tsv"
REFERENCE_CRES_BED="${REFERENCE_CRES_DIR}/dmel_reference_cres.bed"

# ======================================================================
# Expected dataset dimensions and QC parameters
# ======================================================================

EXPECTED_SPECIES="${EXPECTED_SPECIES:-40}"
EXPECTED_REFERENCE_CRES="${EXPECTED_REFERENCE_CRES:-337}"

# Conservative technical QC threshold (not a biological filtering criterion)
MIN_ORTHOLOG_FRACTION="${MIN_ORTHOLOG_FRACTION:-0.10}"

# ======================================================================
# Export variables for child Bash/Python processes
# ======================================================================

export PIPELINE_ROOT
export PROJECT_ROOT

export SCRIPTS_DIR
export ORTHOLOG_RESULTS_DIR
export REFERENCE_CRES_DIR

export SCRM_GENERATED_ROOT
export SCRM_EXTERNAL_ROOT
export GENERATED_GFF_ROOT
export EXTERNAL_GFF_ROOT

export COMBINED_MANIFEST
export DMEL_GFF

export SO_ALL_SPECIES
export SO_ALL_SPECIES_FBGN
export UNRESOLVED_DMEL
export UNRESOLVED_QC_SUMMARY
export DMEL_ORTHOLOG_QC

export REFERENCE_CRES_TSV
export REFERENCE_CRES_BED

export EXPECTED_SPECIES
export EXPECTED_REFERENCE_CRES

export MIN_ORTHOLOG_FRACTION
