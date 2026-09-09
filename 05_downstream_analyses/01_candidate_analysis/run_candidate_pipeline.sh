#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Tier-1 CRE candidate-analysis pipeline
#
# Steps:
#   1. Analyze focal-clade CRE-state contrasts
#   2. Analyze secondary singleton contrasts
#   3. QC Tier-1 candidate evidence
#   4. Prioritize Tier-1 CRE candidates
#   5. Annotate Dmel candidate-gene distances
#   6. Assign primary and secondary candidate genes
#   7. Expand CRE x FBgn associations
#   8. Build gene-level candidate summary
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${PIPELINE_ROOT}/scripts"
CONFIG_FILE="${PIPELINE_ROOT}/config/candidate_config.sh"

[[ -s "$CONFIG_FILE" ]] || {
    echo "ERROR: configuration file missing or empty:" >&2
    echo "  $CONFIG_FILE" >&2
    exit 1
}

# shellcheck source=/dev/null
source "$CONFIG_FILE"

mkdir -p \
    "$FOCAL_CLADES_RESULTS_DIR" \
    "$SECONDARY_RESULTS_DIR" \
    "$TIER1_QC_DIR" \
    "$TIER1_TABLES_DIR"


# ------------------------------------------------------------
# Preconditions
# ------------------------------------------------------------

command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 not found." >&2
    exit 1
}

REQUIRED_FILES=(
    "$CRE_TURNOVER_MATRIX"
    "$COMBINED_MANIFEST"
    "$SPECIES_TRAITS"
    "$FOCAL_CLADES_FILE"
    "$REFERENCE_CRES_TSV"
    "$SO_ALL_SPECIES_FBGN"

    "$SCRIPTS_DIR/01_analyze_focal_clades.py"
    "$SCRIPTS_DIR/02_analyze_secondary_tier1.py"
    "$SCRIPTS_DIR/03_qc_tier1_candidates.py"
    "$SCRIPTS_DIR/04_prioritize_tier1_candidates.py"
    "$SCRIPTS_DIR/05_annotate_candidate_gene_distances.py"
    "$SCRIPTS_DIR/06_assign_candidate_genes.py"
    "$SCRIPTS_DIR/07_explode_candidate_fbgns.py"
    "$SCRIPTS_DIR/08_build_gene_level_summary.py"
)

for file in "${REQUIRED_FILES[@]}"; do
    [[ -s "$file" ]] || {
        echo "ERROR: required file missing or empty:" >&2
        echo "  $file" >&2
        exit 1
    }
done

[[ -d "$TURNOVER_BY_SPECIES_DIR" ]] || {
    echo "ERROR: turnover-by-species directory not found:" >&2
    echo "  $TURNOVER_BY_SPECIES_DIR" >&2
    exit 1
}


# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------

run_step() {
    local step="$1"
    local label="$2"
    shift 2

    echo
    echo "============================================================"
    echo "[${step}/8] ${label}"
    echo "============================================================"

    "$@"
}


# ============================================================
# Start
# ============================================================

echo
echo "============================================================"
echo "Tier-1 CRE candidate analysis"
echo "============================================================"
echo "Started: $(date)"
echo


# ============================================================
# 1. Focal-clade analysis
# ============================================================

run_step 1 "Focal-clade CRE-state contrasts" \
    python3 "$SCRIPTS_DIR/01_analyze_focal_clades.py" \
        --matrix "$CRE_TURNOVER_MATRIX" \
        --groups "$FOCAL_CLADES_FILE" \
        --manifest "$COMBINED_MANIFEST" \
        --traits "$SPECIES_TRAITS" \
        --out-dir "$FOCAL_CLADES_RESULTS_DIR" \
        --summary-out "$FOCAL_CLADE_SUMMARY" \
        --reference-species "$REFERENCE_SPECIES" \
        --expected-reference-cres "$EXPECTED_REFERENCE_CRES"


# ============================================================
# 2. Secondary singleton analysis
# ============================================================

run_step 2 "Secondary singleton CRE-state contrasts" \
    python3 "$SCRIPTS_DIR/02_analyze_secondary_tier1.py" \
        --matrix "$CRE_TURNOVER_MATRIX" \
        --groups "$FOCAL_CLADES_FILE" \
        --outdir "$SECONDARY_RESULTS_DIR" \
        --out-long "$SECONDARY_TIER1_ALL" \
        --out-summary "$SECONDARY_TIER1_RECURRENT" \
        --out-group-summary "$SECONDARY_CLADE_SUMMARY" \
        --metadata-out "$SECONDARY_RUN_METADATA" \
        --expected-reference-cres "$EXPECTED_REFERENCE_CRES" \
        --recurrence-min-clades "$SECONDARY_RECURRENCE_MIN_CLADES"


# ============================================================
# 3. Tier-1 candidate QC
# ============================================================

run_step 3 "Tier-1 candidate QC" \
    python3 "$SCRIPTS_DIR/03_qc_tier1_candidates.py" \
        --focal-dir "$FOCAL_CLADES_RESULTS_DIR" \
        --secondary "$SECONDARY_TIER1_ALL" \
        --reference "$REFERENCE_CRES_TSV" \
        --turnover-dir "$TURNOVER_BY_SPECIES_DIR" \
        --out-summary "$TIER1_QC_SUMMARY" \
        --out-detail "$TIER1_QC_DETAILS" \
        --metadata-out "$TIER1_QC_METADATA" \
        --min-reciprocal-overlap "$TIER1_QC_MIN_RECIPROCAL_OVERLAP" \
        --expected-reference-cres "$EXPECTED_REFERENCE_CRES" \
        --recurrence-min-clades "$SECONDARY_RECURRENCE_MIN_CLADES"


# ============================================================
# 4. Candidate prioritization
# ============================================================

run_step 4 "Prioritize Tier-1 CRE candidates" \
    python3 "$SCRIPTS_DIR/04_prioritize_tier1_candidates.py" \
        --qc "$TIER1_QC_SUMMARY" \
        --reference "$REFERENCE_CRES_TSV" \
        --out-evidence "$TIER1_CLADE_EVIDENCE" \
        --out-prioritized "$TIER1_PRIORITIZED" \
        --metadata-out "$TIER1_PRIORITIZATION_METADATA" \
        --focal-recurrent-min-clades "$FOCAL_RECURRENCE_MIN_CLADES" \
        --secondary-recurrent-min-clades "$SECONDARY_RECURRENCE_MIN_CLADES" \
        --expected-reference-cres "$EXPECTED_REFERENCE_CRES"


# ============================================================
# 5. Candidate-gene distance annotation
# ============================================================

run_step 5 "Annotate candidate-gene distances" \
    python3 "$SCRIPTS_DIR/05_annotate_candidate_gene_distances.py" \
        --candidates "$TIER1_PRIORITIZED" \
        --so-table "$SO_ALL_SPECIES_FBGN" \
        --out "$TIER1_GENE_DISTANCES" \
        --metadata-out "$TIER1_GENE_DISTANCE_METADATA"


# ============================================================
# 6. Candidate-gene assignment
# ============================================================

run_step 6 "Assign candidate genes" \
    python3 "$SCRIPTS_DIR/06_assign_candidate_genes.py" \
        --input "$TIER1_GENE_DISTANCES" \
        --out "$TIER1_GENE_ASSIGNMENTS" \
        --metadata-out "$TIER1_GENE_ASSIGNMENT_METADATA"


# ============================================================
# 7. CRE x FBgn expansion
# ============================================================

run_step 7 "Expand CRE x FBgn associations" \
    python3 "$SCRIPTS_DIR/07_explode_candidate_fbgns.py" \
        --input "$TIER1_GENE_ASSIGNMENTS" \
        --out "$TIER1_FBGN_EXPLODED" \
        --metadata-out "$TIER1_FBGN_EXPLODED_METADATA"


# ============================================================
# 8. Gene-level summary
# ============================================================

run_step 8 "Build gene-level candidate summary" \
    python3 "$SCRIPTS_DIR/08_build_gene_level_summary.py" \
        --input "$TIER1_FBGN_EXPLODED" \
        --out "$TIER1_GENE_LEVEL_SUMMARY" \
        --metadata-out "$TIER1_GENE_LEVEL_METADATA"


# ============================================================
# Final output QC
# ============================================================

EXPECTED_OUTPUTS=(
    "$FOCAL_CLADE_SUMMARY"

    "$SECONDARY_TIER1_ALL"
    "$SECONDARY_TIER1_RECURRENT"
    "$SECONDARY_CLADE_SUMMARY"
    "$SECONDARY_RUN_METADATA"

    "$TIER1_QC_SUMMARY"
    "$TIER1_QC_DETAILS"
    "$TIER1_QC_METADATA"

    "$TIER1_CLADE_EVIDENCE"
    "$TIER1_PRIORITIZED"
    "$TIER1_PRIORITIZATION_METADATA"

    "$TIER1_GENE_DISTANCES"
    "$TIER1_GENE_DISTANCE_METADATA"
    "$TIER1_GENE_ASSIGNMENTS"
    "$TIER1_GENE_ASSIGNMENT_METADATA"
    "$TIER1_FBGN_EXPLODED"
    "$TIER1_FBGN_EXPLODED_METADATA"
    "$TIER1_GENE_LEVEL_SUMMARY"
    "$TIER1_GENE_LEVEL_METADATA"
)

for file in "${EXPECTED_OUTPUTS[@]}"; do
    [[ -s "$file" ]] || {
        echo
        echo "ERROR: expected output missing or empty:" >&2
        echo "  $file" >&2
        exit 1
    }
done


# ============================================================
# Finish
# ============================================================

echo
echo "============================================================"
echo "CANDIDATE ANALYSIS COMPLETE"
echo "============================================================"
echo "Finished: $(date)"
echo
echo "Main outputs:"
echo "  $TIER1_QC_SUMMARY"
echo "  $TIER1_PRIORITIZED"
echo "  $TIER1_GENE_ASSIGNMENTS"
echo "  $TIER1_GENE_LEVEL_SUMMARY"
echo
