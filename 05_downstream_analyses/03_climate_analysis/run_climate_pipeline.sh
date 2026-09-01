#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Climate-analysis pipeline
# ============================================================

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${PIPELINE_ROOT}/config/climate_config.sh"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: config file not found:"
    echo "${CONFIG_FILE}"
    exit 1
fi

source "${CONFIG_FILE}"

mkdir -p "${RESULTS_DIR}" "${FIGURES_DIR}"


echo
echo "============================================================"
echo "Climate and focal Tier-1 analysis"
echo "============================================================"
echo


# ============================================================
# Step 1: Species-level turnover by climatic zone
# ============================================================

echo "Step 1/7: Species-level turnover by climatic zone"

python3 "${PIPELINE_ROOT}/01_plot_turnover_by_climate.py" \
    --species-summary "${SPECIES_SUMMARY}" \
    --manifest "${COMBINED_MANIFEST}" \
    --traits "${SPECIES_TRAITS}" \
    --out-summary "${CLIMATE_TURNOVER_SUMMARY}" \
    --out-descriptive "${CLIMATE_TURNOVER_DESCRIPTIVE}" \
    --out-kruskal "${CLIMATE_TURNOVER_KRUSKAL}" \
    --out-png "${TURNOVER_BY_CLIMATE_PNG}" \
    --out-pdf "${TURNOVER_BY_CLIMATE_PDF}"


# ============================================================
# Step 2: Phylogenetically controlled climate test
# ============================================================

echo
echo "Step 2/7: Phylogenetically controlled climate test"

Rscript "${PIPELINE_ROOT}/02_test_turnover_by_climate.R" \
    --input "${CLIMATE_TURNOVER_SUMMARY}" \
    --tree "${PHYLOGENY_TREE}" \
    --out-model "${PGLS_CLIMATE_MODEL_COMPARISON}" \
    --out-coefficients "${PGLS_CLIMATE_COEFFICIENTS}" \
    --out-summary "${PGLS_CLIMATE_SUMMARY}" \
    --out-predictions "${PGLS_CLIMATE_PREDICTIONS}" \
    --out-stats "${PGLS_CLIMATE_STATS}"


# ============================================================
# Step 3: PGLS climate figure
# ============================================================

echo
echo "Step 3/7: PGLS climate figure"

Rscript "${PIPELINE_ROOT}/03_plot_climate_pgls.R" \
    --input "${CLIMATE_TURNOVER_SUMMARY}" \
    --predictions "${PGLS_CLIMATE_PREDICTIONS}" \
    --stats "${PGLS_CLIMATE_STATS}" \
    --out-png "${CLIMATE_PGLS_PNG}" \
    --out-pdf "${CLIMATE_PGLS_PDF}"


# ============================================================
# Step 4: Focal Tier-1 enrichment
# ============================================================

echo
echo "Step 4/7: Focal Tier-1 enrichment"

python3 "${PIPELINE_ROOT}/04_test_focal_tier1_enrichment.py" \
    --matrix "${CRE_TURNOVER_MATRIX}" \
    --groups "${FOCAL_CLADES_FILE}" \
    --out-candidates "${FOCAL_TIER1_CANDIDATES}" \
    --out-enrichment "${FOCAL_TIER1_ENRICHMENT}"


# ============================================================
# Step 5: Focal Tier-1 enrichment figure
# ============================================================

echo
echo "Step 5/7: Focal Tier-1 enrichment figure"

python3 "${PIPELINE_ROOT}/05_plot_focal_tier1_enrichment.py" \
    --input "${FOCAL_TIER1_ENRICHMENT}" \
    --out-png "${FOCAL_TIER1_ENRICHMENT_PNG}" \
    --out-pdf "${FOCAL_TIER1_ENRICHMENT_PDF}"


# ============================================================
# Step 6: Focal Tier-1 composition figure
# ============================================================

echo
echo "Step 6/7: Focal Tier-1 summary figure"

python3 "${PIPELINE_ROOT}/06_plot_focal_summary.py" \
    --input "${FOCAL_TIER1_ENRICHMENT}" \
    --out-png "${FOCAL_TIER1_SUMMARY_PNG}" \
    --out-pdf "${FOCAL_TIER1_SUMMARY_PDF}"


# ============================================================
# Step 7: Export Azteca focal candidates
# ============================================================

echo
echo "Step 7/7: Export Azteca focal Tier-1 candidates"

python3 "${PIPELINE_ROOT}/07_export_azteca_focal_candidates.py" \
    --candidates "${CANDIDATES_WITH_SENSITIVITY}" \
    --focal-events "${FOCAL_TIER1_CANDIDATES}" \
    --clade "${AZTECA_FOCAL_CLADE}" \
    --out "${AZTECA_FOCAL_CANDIDATES}"


echo
echo "============================================================"
echo "CLIMATE ANALYSIS COMPLETE"
echo "============================================================"
echo

echo "Main figures:"
echo "  ${CLIMATE_PGLS_PNG}"
echo "  ${FOCAL_TIER1_ENRICHMENT_PNG}"

echo
echo "Additional figures:"
echo "  ${TURNOVER_BY_CLIMATE_PNG}"
echo "  ${FOCAL_TIER1_SUMMARY_PNG}"

echo
echo "Azteca candidate table:"
echo "  ${AZTECA_FOCAL_CANDIDATES}"
echo
