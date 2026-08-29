#!/usr/bin/env bash
# Zentrale Konfiguration für die Multi-Spezies-SCRMshaw-Pipeline.
# Dieses File wird von allen Slurm-Jobs eingelesen.

# Projektwurzel: standardmäßig der Ordner, in dem dieses Paket liegt.
PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Datenquelle des Papers Dhakad et al. (2026), aktuelle Zenodo-Version.
ZENODO_RECORD_ID="18453526"
ZENODO_BASE_URL="https://zenodo.org/records/${ZENODO_RECORD_ID}/files"
GENOMES_ARCHIVE="genomes.tar.gz"
ANNOTATIONS_ARCHIVE="annotations.tar.gz"
SPECIES_SUMMARY="Species_summary_301Fly.xlsx"

# Checksummen der Zenodo-Version 18453526.
GENOMES_MD5="bca079304da4dbe8e0c9998fc049eb03"
ANNOTATIONS_MD5="d7cd2d6d0b98b4d51036b05c619c590b"

# Eingaben/Ausgaben
SPECIES_FILE="${PIPELINE_ROOT}/species.txt"
DATA_DIR="${PIPELINE_ROOT}/data"
DOWNLOAD_DIR="${DATA_DIR}/downloads"
EXTRACT_DIR="${DATA_DIR}/selected"
MANIFEST="${DATA_DIR}/manifest.tsv"
RUNS_DIR="${PIPELINE_ROOT}/runs"
RESULTS_DIR="${PIPELINE_ROOT}/results"
SOFTWARE_DIR="${PIPELINE_ROOT}/software"
LOG_DIR="${PIPELINE_ROOT}/logs"

# Conda/Miniforge
CONDA_BASE="${CONDA_BASE:-${HOME}/miniforge3}"
SCRM_ENV="${SCRM_ENV:-scrmshaw}"
POST_ENV="${POST_ENV:-scrm_postproc}"

# SCRMshaw
SCRM_ROOT="${SOFTWARE_DIR}/SCRMshaw_HD"
UTILITY_ROOT="${SOFTWARE_DIR}/UtilityPrograms"
POSTPROC_ROOT="${SOFTWARE_DIR}/post_processing_SCRMshaw_pipeline"
TRAINING_ROOT="${SOFTWARE_DIR}/dmel_training_sets"
TRAINING_SET="${TRAINING_ROOT}/final_combined_48Tsets/adult_muscle"
TRAINING_LIST="${PIPELINE_ROOT}/config/trainingSet_adult_muscle.lst"
TRF_BIN="${SOFTWARE_DIR}/trf409.linux64"

# Standardmäßig nur IMM: deutlich schneller für 22 Spezies.
# Für alle Methoden:
# SCORING_FLAGS="--imm --hexmcd --pac"
SCORING_FLAGS="--imm"
THITW="${THITW:-10000}"

# Tandem-repeat masking zusätzlich zu den im Paper bereits soft-maskierten Genomen.
RUN_TRF="${RUN_TRF:-1}"

# Slurm-Ressourcen. Partition/QOS leer lassen, wenn die Cluster-Defaults gelten.
SLURM_PARTITION="${SLURM_PARTITION:-}"
SLURM_QOS="${SLURM_QOS:-}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"

DOWNLOAD_TIME="${DOWNLOAD_TIME:-12:00:00}"
DOWNLOAD_MEM="${DOWNLOAD_MEM:-5G}"
DOWNLOAD_CPUS="${DOWNLOAD_CPUS:-1}"

PREP_TIME="${PREP_TIME:-12:00:00}"
PREP_MEM="${PREP_MEM:-5G}"
PREP_CPUS="${PREP_CPUS:-1}"
MAX_PREP_JOBS="${MAX_PREP_JOBS:-10}"

OFFSET_TIME="${OFFSET_TIME:-72:00:00}"
OFFSET_MEM="${OFFSET_MEM:-5G}"
OFFSET_CPUS="${OFFSET_CPUS:-1}"
MAX_OFFSET_JOBS="${MAX_OFFSET_JOBS:-10}"

POST_TIME="${POST_TIME:-12:00:00}"
POST_MEM="${POST_MEM:-5G}"
POST_CPUS="${POST_CPUS:-1}"
MAX_POST_JOBS="${MAX_POST_JOBS:-8}"

SUMMARY_TIME="${SUMMARY_TIME:-01:00:00}"
SUMMARY_MEM="${SUMMARY_MEM:-4G}"

# Bei 1 werden bestehende Downloads/Outputs wiederverwendet.
RESUME="${RESUME:-1}"
