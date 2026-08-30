#!/usr/bin/env bash
set -euo pipefail

# ======================================================================
# Pipeline roots
# ======================================================================

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_ROOT="$(cd "$PIPELINE_ROOT/.." && pwd)"

SCRM_GENERATED_ROOT="${PROJECT_ROOT}/01_scrmshaw/generation"
SCRM_EXTERNAL_ROOT="${PROJECT_ROOT}/01_scrmshaw/external"
MAPPING_ROOT="${PROJECT_ROOT}/02_mapping_orthologs"

# ======================================================================
# Upstream input
# ======================================================================

COMBINED_MANIFEST="${SCRM_EXTERNAL_ROOT}/combined_manifest.tsv"

GENERATED_GENOME_ROOT="${SCRM_GENERATED_ROOT}/data/selected"
EXTERNAL_GENOME_ROOT="${SCRM_EXTERNAL_ROOT}/external_data/selected"

REFERENCE_CRES_BED="${MAPPING_ROOT}/reference_cres/dmel_reference_cres.bed"

# ======================================================================
# Pairwise WGA directories
# ======================================================================

BIN_DIR="${PIPELINE_ROOT}/bin"

TWOBIT_DIR="${PIPELINE_ROOT}/twobit"
CHROM_SIZES_DIR="${PIPELINE_ROOT}/chrom_sizes"

ALIGNMENT_DIR="${PIPELINE_ROOT}/alignments_dmel"
LIFTED_CRES_DIR="${PIPELINE_ROOT}/lifted_cres_dmel"

LOG_DIR="${PIPELINE_ROOT}/logs"

TARGET_SPECIES_FILE="${PIPELINE_ROOT}/target_species.txt"

# ======================================================================
# Reference species
# ======================================================================

REFERENCE_SPECIES="${REFERENCE_SPECIES:-d_melanogaster}"

# ======================================================================
# LASTZ parameters
# ======================================================================

LASTZ_FORMAT="${LASTZ_FORMAT:-axt}"
LASTZ_AMBIGUOUS="${LASTZ_AMBIGUOUS:-iupac}"
LASTZ_NOTRANSITION="${LASTZ_NOTRANSITION:-1}"

LASTZ_STEP="${LASTZ_STEP:-20}"
LASTZ_SEED="${LASTZ_SEED:-12of19}"
LASTZ_HSPTHRESH="${LASTZ_HSPTHRESH:-2200}"
LASTZ_GAPPEDTHRESH="${LASTZ_GAPPEDTHRESH:-4000}"
LASTZ_YDROP="${LASTZ_YDROP:-3400}"
LASTZ_INNER="${LASTZ_INNER:-2000}"

# ======================================================================
# UCSC chain/net parameters
# ======================================================================

AXTCHAIN_LINEAR_GAP="${AXTCHAIN_LINEAR_GAP:-medium}"
AXTCHAIN_MIN_SCORE="${AXTCHAIN_MIN_SCORE:-3000}"

# ======================================================================
# liftOver parameters
# ======================================================================

LIFTOVER_MIN_MATCH="${LIFTOVER_MIN_MATCH:-0.50}"

# ======================================================================
# SLURM resources
# ======================================================================

SLURM_NODE="${SLURM_NODE:-abacus-2}"

WGA_MEM="${WGA_MEM:-5G}"
WGA_CPUS="${WGA_CPUS:-1}"
MAX_WGA_JOBS="${MAX_WGA_JOBS:-10}"

# ======================================================================
# Export
# ======================================================================

export PIPELINE_ROOT
export PROJECT_ROOT

export SCRM_GENERATED_ROOT
export SCRM_EXTERNAL_ROOT
export MAPPING_ROOT

export COMBINED_MANIFEST
export GENERATED_GENOME_ROOT
export EXTERNAL_GENOME_ROOT
export REFERENCE_CRES_BED

export BIN_DIR
export TWOBIT_DIR
export CHROM_SIZES_DIR
export ALIGNMENT_DIR
export LIFTED_CRES_DIR
export LOG_DIR
export TARGET_SPECIES_FILE

export REFERENCE_SPECIES

export LASTZ_FORMAT
export LASTZ_AMBIGUOUS
export LASTZ_NOTRANSITION
export LASTZ_STEP
export LASTZ_SEED
export LASTZ_HSPTHRESH
export LASTZ_GAPPEDTHRESH
export LASTZ_YDROP
export LASTZ_INNER

export AXTCHAIN_LINEAR_GAP
export AXTCHAIN_MIN_SCORE

export LIFTOVER_MIN_MATCH

export SLURM_NODE
export WGA_MEM
export WGA_CPUS
export MAX_WGA_JOBS
