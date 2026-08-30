#!/usr/bin/env bash

# ======================================================================
# Analysis configuration
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
# Secondary Tier-1 analysis
# ======================================================================

SECONDARY_RESULTS_DIR="${RESULTS_DIR}/secondary_clades"
SECONDARY_TIER1_ALL="${SECONDARY_RESULTS_DIR}/secondary_tier1_all_clades.tsv"
SECONDARY_TIER1_RECURRENT="${SECONDARY_RESULTS_DIR}/secondary_tier1_recurrent_cres.tsv"
SECONDARY_CLADE_SUMMARY="${SECONDARY_RESULTS_DIR}/secondary_clade_summary.tsv"
SECONDARY_RUN_METADATA="${SECONDARY_RESULTS_DIR}/secondary_tier1_run_metadata.tsv"

# Minimum number of clades required to classify a CRE as recurrent.
SECONDARY_RECURRENCE_MIN_CLADES="${SECONDARY_RECURRENCE_MIN_CLADES:-2}"


# ======================================================================
# Tier-1 candidate QC
# ======================================================================

MAPPING_ROOT="${PROJECT_ROOT}/02_mapping_orthologs"

# D. melanogaster reference CRE annotation.
REFERENCE_CRES_TSV="${MAPPING_ROOT}/reference_cres/dmel_reference_cres.tsv"

# Detailed CRE classifications for individual target species.
TURNOVER_BY_SPECIES_DIR="${CLASSIFICATION_ROOT}/results/turnover_by_species"

# Tier-1 QC outputs.
TIER1_QC_DIR="${RESULTS_DIR}/qc"
TIER1_QC_SUMMARY="${TIER1_QC_DIR}/tier1_candidate_qc_summary.tsv"
TIER1_QC_DETAILS="${TIER1_QC_DIR}/tier1_candidate_qc_details.tsv"
TIER1_QC_METADATA="${TIER1_QC_DIR}/tier1_candidate_qc_metadata.tsv"

# Minimum reciprocal overlap required for a positional CRE match.
# This should match the threshold used during CRE classification.
TIER1_QC_MIN_RECIPROCAL_OVERLAP="${TIER1_QC_MIN_RECIPROCAL_OVERLAP:-0.50}"


# ======================================================================
# Tier-1 candidate prioritization
# ======================================================================

TIER1_TABLES_DIR="${RESULTS_DIR}/tables"
TIER1_CLADE_EVIDENCE="${TIER1_TABLES_DIR}/tier1_candidate_clade_evidence.tsv"
TIER1_PRIORITIZED="${TIER1_TABLES_DIR}/tier1_candidates_prioritized.tsv"
TIER1_PRIORITIZATION_METADATA="${TIER1_TABLES_DIR}/tier1_candidate_prioritization_metadata.tsv"

# Minimum number of focal Tier-1 clades required to classify
# a candidate as recurrent focal support.
FOCAL_RECURRENCE_MIN_CLADES="${FOCAL_RECURRENCE_MIN_CLADES:-2}"


# ======================================================================
# Candidate gene-distance annotation
# ======================================================================

# Combined SCRMshaw orthology table with FBgn annotations.
SO_ALL_SPECIES_FBGN="${MAPPING_ROOT}/results/SO_all_species_fbgn.tsv"

# Gene-distance annotation for prioritized Tier-1 candidates.
TIER1_GENE_DISTANCES="${TIER1_TABLES_DIR}/tier1_candidate_gene_distances.tsv"

TIER1_GENE_DISTANCE_METADATA="${TIER1_TABLES_DIR}/tier1_candidate_gene_distances_metadata.tsv"


# ======================================================================
# Candidate gene assignment
# ======================================================================

TIER1_GENE_ASSIGNMENTS="${TIER1_TABLES_DIR}/tier1_candidate_gene_assignments.tsv"
TIER1_GENE_ASSIGNMENT_METADATA="${TIER1_TABLES_DIR}/tier1_candidate_gene_assignments_metadata.tsv"


# ======================================================================
# Candidate FBgn expansion
# ======================================================================

TIER1_FBGN_EXPLODED="${TIER1_TABLES_DIR}/tier1_candidate_fbgn_exploded.tsv"
TIER1_FBGN_EXPLODED_METADATA="${TIER1_TABLES_DIR}/tier1_candidate_fbgn_exploded_metadata.tsv"


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

export SECONDARY_RESULTS_DIR
export SECONDARY_TIER1_ALL
export SECONDARY_TIER1_RECURRENT
export SECONDARY_CLADE_SUMMARY
export SECONDARY_RUN_METADATA
export SECONDARY_RECURRENCE_MIN_CLADES

export REFERENCE_CRES_TSV
export TURNOVER_BY_SPECIES_DIR

export TIER1_QC_DIR
export TIER1_QC_SUMMARY
export TIER1_QC_DETAILS
export TIER1_QC_METADATA

export TIER1_QC_MIN_RECIPROCAL_OVERLAP

export TIER1_TABLES_DIR
export TIER1_CLADE_EVIDENCE
export TIER1_PRIORITIZED
export TIER1_PRIORITIZATION_METADATA

export FOCAL_RECURRENCE_MIN_CLADES

export SO_ALL_SPECIES_FBGN
export TIER1_GENE_DISTANCES
export TIER1_GENE_DISTANCE_METADATA

export TIER1_GENE_ASSIGNMENTS
export TIER1_GENE_ASSIGNMENT_METADATA

export TIER1_FBGN_EXPLODED
export TIER1_FBGN_EXPLODED_METADATA

