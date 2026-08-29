```bash
#!/usr/bin/env bash

# Central configuration for the multi-species SCRMshaw pipeline.
# This file is sourced by all pipeline scripts and SLURM jobs.

# ======================================================================
# Pipeline root
# ======================================================================

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ======================================================================
# Source dataset
# ======================================================================

# Genome and annotation dataset from Dhakad et al. (2026), Zenodo record.
ZENODO_RECORD_ID="18453526"
ZENODO_BASE_URL="https://zenodo.org/records/${ZENODO_RECORD_ID}/files"

GENOMES_ARCHIVE="genomes.tar.gz"
ANNOTATIONS_ARCHIVE="annotations.tar.gz"
SPECIES_SUMMARY="Species_summary_301Fly.xlsx"

# MD5 checksums for Zenodo record 18453526.
GENOMES_MD5="bca079304da4dbe8e0c9998fc049eb03"
ANNOTATIONS_MD5="d7cd2d6d0b98b4d51036b05c619c590b"

# ======================================================================
# Pipeline input and output directories
# ======================================================================

SPECIES_FILE="${PIPELINE_ROOT}/species.txt"

DATA_DIR="${PIPELINE_ROOT}/data"
DOWNLOAD_DIR="${DATA_DIR}/downloads"
EXTRACT_DIR="${DATA_DIR}/selected"

# Species manifest generated during genome preparation
MANIFEST="${DATA_DIR}/manifest.tsv"

# Temporary SCRMshaw working directories
RUNS_DIR="${PIPELINE_ROOT}/runs"

RESULTS_DIR="${PIPELINE_ROOT}/results"

# Third-party software required by the pipeline
SOFTWARE_DIR="${PIPELINE_ROOT}/software"

LOG_DIR="${PIPELINE_ROOT}/logs"

# ======================================================================
# Conda environments
# ======================================================================

CONDA_BASE="${CONDA_BASE:-${HOME}/miniforge3}"

SCRM_ENV="${SCRM_ENV:-scrmshaw}"
POST_ENV="${POST_ENV:-scrm_postproc}"

# ======================================================================
# SCRMshaw software and training data
# ======================================================================

SCRM_ROOT="${SOFTWARE_DIR}/SCRMshaw_HD"
UTILITY_ROOT="${SOFTWARE_DIR}/UtilityPrograms"
POSTPROC_ROOT="${SOFTWARE_DIR}/post_processing_SCRMshaw_pipeline"

TRAINING_ROOT="${SOFTWARE_DIR}/dmel_training_sets"

# SCRMshaw training set used for all generated species
TRAINING_SET="${TRAINING_ROOT}/final_combined_48Tsets/adult_muscle"

# Training-set configuration supplied to SCRMshaw
TRAINING_LIST="${PIPELINE_ROOT}/config/trainingSet_adult_muscle.lst"

# Tandem Repeat Finder executable
TRF_BIN="${SOFTWARE_DIR}/trf409.linux64"

# ======================================================================
# SCRMshaw analysis parameters
# ======================================================================

# SCORING_FLAGS="--imm --hexmcd --pac"
SCORING_FLAGS="--imm"

# Number of top SCRMshaw hits retained per offset before post-processing
THITW="${THITW:-10000}"

# Tandem Repeat Finder masking
RUN_TRF="${RUN_TRF:-1}"

# ======================================================================
# SLURM configuration
# ======================================================================

# Leave partition, QOS, and account empty to use cluster defaults.
SLURM_PARTITION="${SLURM_PARTITION:-}"
SLURM_QOS="${SLURM_QOS:-}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"

# ----------------------------------------------------------------------
# Download and extraction
# ----------------------------------------------------------------------

DOWNLOAD_TIME="${DOWNLOAD_TIME:-12:00:00}"
DOWNLOAD_MEM="${DOWNLOAD_MEM:-5G}"
DOWNLOAD_CPUS="${DOWNLOAD_CPUS:-1}"

# ----------------------------------------------------------------------
# Genome preparation
# ----------------------------------------------------------------------

PREP_TIME="${PREP_TIME:-12:00:00}"
PREP_MEM="${PREP_MEM:-5G}"
PREP_CPUS="${PREP_CPUS:-1}"
MAX_PREP_JOBS="${MAX_PREP_JOBS:-10}"

# ----------------------------------------------------------------------
# SCRMshaw offset jobs
# ----------------------------------------------------------------------

OFFSET_TIME="${OFFSET_TIME:-72:00:00}"
OFFSET_MEM="${OFFSET_MEM:-5G}"
OFFSET_CPUS="${OFFSET_CPUS:-1}"
MAX_OFFSET_JOBS="${MAX_OFFSET_JOBS:-10}"

# ----------------------------------------------------------------------
# SCRMshaw post-processing
# ----------------------------------------------------------------------

POST_TIME="${POST_TIME:-12:00:00}"
POST_MEM="${POST_MEM:-5G}"
POST_CPUS="${POST_CPUS:-1}"
MAX_POST_JOBS="${MAX_POST_JOBS:-8}"

# ----------------------------------------------------------------------
# Summary generation
# ----------------------------------------------------------------------

SUMMARY_TIME="${SUMMARY_TIME:-01:00:00}"
SUMMARY_MEM="${SUMMARY_MEM:-4G}"

# ======================================================================
# Resume behavior
# ======================================================================

# Reuse existing downloads and completed outputs when set to 1.
RESUME="${RESUME:-1}"
```
