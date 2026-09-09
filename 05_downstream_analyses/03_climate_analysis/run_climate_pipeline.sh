#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Climate-analysis pipeline
#
# Steps:
#   1. Species-level turnover by climatic zone
#   2. Phylogenetically controlled climate test
#   3. PGLS climate figure
#   4. Focal Tier-1 enrichment
#   5. Focal Tier-1 enrichment figure
#   6. Focal Tier-1 composition figure
# ============================================================

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${PIPELINE_ROOT}/config/climate_config.sh"

[[ -s "${CONFIG_FILE}" ]] || {
    echo "ERROR: configuration file missing or empty:" >&2
    echo "  ${CONFIG_FILE}" >&2
    exit 1
}

# shellcheck source=/dev/null
source "${CONFIG_FILE}"

mkdir -p "${RESULTS_DIR}" "${FIGURES_DIR}"


# ============================================================
# Required files
# ============================================================

REQUIRED_FILES=(
    "${SPECIES_SUMMARY}"
    "${COMBINED_MANIFEST}"
    "${SPECIES_TRAITS}"
    "${PHYLOGENY_TREE}"
    "${CRE_TURNOVER_MATRIX}"
    "${FOCAL_CLADES_FILE}"

    "${PIPELINE_ROOT}/01_plot_turnover_by_climate.py"
    "${PIPELINE_ROOT}/02_test_turnover_by_climate.R"
    "${PIPELINE_ROOT}/03_plot_climate_pgls.R"
    "${PIPELINE_ROOT}/04_test_focal_tier1_enrichment.py"
    "${PIPELINE_ROOT}/05_plot_focal_tier1_enrichment.py"
    "${PIPELINE_ROOT}/06_plot_focal_summary.py"
)

for file in "${REQUIRED_FILES[@]}"; do
    [[ -s "${file}" ]] || {
        echo "ERROR: required file missing or empty:" >&2
        echo "  ${file}" >&2
        exit 1
    }
done


# ============================================================
# Start
# ============================================================

echo
echo "============================================================"
echo "Climate and focal Tier-1 analysis"
echo "============================================================"
echo "Started: $(date)"
echo


# ============================================================
# 1. Species-level turnover by climatic zone
# ============================================================

echo "[1/6] Species-level turnover by climatic zone"

python3 "${PIPELINE_ROOT}/01_plot_turnover_by_climate.py" \
    --species-summary "${SPECIES_SUMMARY}" \
    --manifest "${COMBINED_MANIFEST}" \
    --traits "${SPECIES_TRAITS}" \
    --out-summary "${CLIMATE_TURNOVER_SUMMARY}" \
    --out-descriptive "${CLIMATE_TURNOVER_DESCRIPTIVE}" \
    --out-kruskal "${CLIMATE_TURNOVER_KRUSKAL}" \
    --out-png "${TURNOVER_BY_CLIMATE_PNG}" \
    --out-pdf "${TURNOVER_BY_CLIMATE_PDF}"

echo


# ============================================================
# 2. Phylogenetically controlled climate analysis
# ============================================================

echo "[2/6] Phylogenetically controlled climate analysis"

Rscript "${PIPELINE_ROOT}/02_test_turnover_by_climate.R" \
    --input "${CLIMATE_TURNOVER_SUMMARY}" \
    --tree "${PHYLOGENY_TREE}" \
    --out-model "${PGLS_CLIMATE_MODEL_COMPARISON}" \
    --out-coefficients "${PGLS_CLIMATE_COEFFICIENTS}" \
    --out-summary "${PGLS_CLIMATE_SUMMARY}" \
    --out-predictions "${PGLS_CLIMATE_PREDICTIONS}" \
    --out-stats "${PGLS_CLIMATE_STATS}"

echo


# ============================================================
# 3. PGLS climate figure
# ============================================================

echo "[3/6] PGLS climate figure"

Rscript "${PIPELINE_ROOT}/03_plot_climate_pgls.R" \
    --input "${CLIMATE_TURNOVER_SUMMARY}" \
    --predictions "${PGLS_CLIMATE_PREDICTIONS}" \
    --stats "${PGLS_CLIMATE_STATS}" \
    --out-png "${CLIMATE_PGLS_PNG}" \
    --out-pdf "${CLIMATE_PGLS_PDF}"

echo


# ============================================================
# 4. Focal Tier-1 enrichment
# ============================================================

echo "[4/6] Focal Tier-1 enrichment"

python3 "${PIPELINE_ROOT}/04_test_focal_tier1_enrichment.py" \
    --matrix "${CRE_TURNOVER_MATRIX}" \
    --groups "${FOCAL_CLADES_FILE}" \
    --out-candidates "${FOCAL_TIER1_CANDIDATES}" \
    --out-enrichment "${FOCAL_TIER1_ENRICHMENT}"

echo


# ============================================================
# 5. Focal Tier-1 enrichment figure
# ============================================================

echo "[5/6] Focal Tier-1 enrichment figure"

python3 "${PIPELINE_ROOT}/05_plot_focal_tier1_enrichment.py" \
    --input "${FOCAL_TIER1_ENRICHMENT}" \
    --out-png "${FOCAL_TIER1_ENRICHMENT_PNG}" \
    --out-pdf "${FOCAL_TIER1_ENRICHMENT_PDF}"

echo


# ============================================================
# 6. Focal Tier-1 composition figure
# ============================================================

echo "[6/6] Focal Tier-1 composition figure"

python3 "${PIPELINE_ROOT}/06_plot_focal_summary.py" \
    --input "${FOCAL_TIER1_ENRICHMENT}" \
    --out-png "${FOCAL_TIER1_SUMMARY_PNG}" \
    --out-pdf "${FOCAL_TIER1_SUMMARY_PDF}"


# ============================================================
# Final output QC
# ============================================================

EXPECTED_OUTPUTS=(
    "${CLIMATE_TURNOVER_SUMMARY}"
    "${CLIMATE_TURNOVER_DESCRIPTIVE}"
    "${CLIMATE_TURNOVER_KRUSKAL}"

    "${PGLS_CLIMATE_MODEL_COMPARISON}"
    "${PGLS_CLIMATE_COEFFICIENTS}"
    "${PGLS_CLIMATE_SUMMARY}"
    "${PGLS_CLIMATE_PREDICTIONS}"
    "${PGLS_CLIMATE_STATS}"

    "${FOCAL_TIER1_CANDIDATES}"
    "${FOCAL_TIER1_ENRICHMENT}"

    "${TURNOVER_BY_CLIMATE_PNG}"
    "${TURNOVER_BY_CLIMATE_PDF}"
    "${CLIMATE_PGLS_PNG}"
    "${CLIMATE_PGLS_PDF}"
    "${FOCAL_TIER1_ENRICHMENT_PNG}"
    "${FOCAL_TIER1_ENRICHMENT_PDF}"
    "${FOCAL_TIER1_SUMMARY_PNG}"
    "${FOCAL_TIER1_SUMMARY_PDF}"
)

for file in "${EXPECTED_OUTPUTS[@]}"; do
    [[ -s "${file}" ]] || {
        echo "ERROR: expected output missing or empty:" >&2
        echo "  ${file}" >&2
        exit 1
    }
done


# ============================================================
# Finish
# ============================================================

echo
echo "============================================================"
echo "CLIMATE ANALYSIS COMPLETE"
echo "============================================================"
echo "Finished: $(date)"
echo
echo "Main results:"
echo "  ${CLIMATE_TURNOVER_SUMMARY}"
echo "  ${PGLS_CLIMATE_MODEL_COMPARISON}"
echo "  ${FOCAL_TIER1_ENRICHMENT}"
echo
echo "Main figures:"
echo "  ${CLIMATE_PGLS_PNG}"
echo "  ${FOCAL_TIER1_ENRICHMENT_PNG}"
echo
