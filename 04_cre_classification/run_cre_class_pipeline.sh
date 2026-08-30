#!/usr/bin/env bash
set -euo pipefail


# ============================================================
# CRE-state classification pipeline
#
# Purpose:
#   Classify D. melanogaster reference CREs across all target
#   species, combine species-specific classifications, and
#   generate species-level QC summaries.
#
# Output:
#   - results/turnover_by_species/*.tsv
#   - results/cre_turnover_all_species.tsv
#   - results/cre_turnover_matrix.tsv
#   - results/species_summary.tsv
#   - results/species_qc_summary.tsv
#
# Configuration:
#   config/classification_config.sh
# ============================================================


# ============================================================
# Configuration
# ============================================================

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "$ROOT/config/classification_config.sh"


# ============================================================
# Prepare output directories
# ============================================================

mkdir -p \
    "$RESULTS_DIR" \
    "$TURNOVER_BY_SPECIES_DIR"


# ============================================================
# Header
# ============================================================

echo "============================================================"
echo "CRE-state classification pipeline"
echo "============================================================"
echo "Started                : $(date)"
echo "Reciprocal overlap     : $RECIPROCAL_OVERLAP"
echo "Local gene distance    : $LOCAL_GENE_DISTANCE bp"
echo "Expected reference CREs: $EXPECTED_REFERENCE_CRES"
echo


# ============================================================
# 1. Check required inputs
# ============================================================

echo "[1/4] Checking classification inputs..."

required_files=(

    "$TARGET_SPECIES_FILE"
    "$REFERENCE_CRES_TSV"
    "$SO_ALL_SPECIES_FBGN"

    "$SCRIPTS_DIR/01_classify_target_cres.py"
    "$SCRIPTS_DIR/02_combine_turnover_results.py"
    "$SCRIPTS_DIR/03_build_species_qc_summary.py"
)

for file in "${required_files[@]}"; do

    [[ -s "$file" ]] || {
        echo "ERROR: required file missing or empty:" >&2
        echo "  $file" >&2
        exit 1
    }

done


N_REF_CRES="$(
    awk '
        NR > 1 && NF > 0 {n++}
        END {print n+0}
    ' "$REFERENCE_CRES_TSV"
)"

echo "Reference CREs: $N_REF_CRES"

if [[ "$N_REF_CRES" -ne "$EXPECTED_REFERENCE_CRES" ]]; then

    echo "ERROR: expected $EXPECTED_REFERENCE_CRES reference CREs," >&2
    echo "but found $N_REF_CRES." >&2
    exit 1

fi

echo
echo "[1/4] PASSED"
echo


# ============================================================
# 2. Classify all target species
# ============================================================

echo "[2/4] Classifying target species..."
echo

n_done=0

while IFS= read -r sp; do

    [[ -n "$sp" ]] || continue

    lifted="$LIFTED_CRES_DIR/$sp/dmel_reference_cres.${sp}.bed"

    if [[ ! -f "$lifted" ]]; then
        echo "ERROR: lifted BED missing for $sp:" >&2
        echo "  $lifted" >&2
        exit 1
    fi

    out="$TURNOVER_BY_SPECIES_DIR/dmel_to_${sp}_cre_turnover.tsv"

    echo "============================================================"
    echo "$sp"
    echo "============================================================"

    python3 "$SCRIPTS_DIR/01_classify_target_cres.py" \
        --species "$sp" \
        --reference "$REFERENCE_CRES_TSV" \
        --lifted "$lifted" \
        --predictions "$SO_ALL_SPECIES_FBGN" \
        --min-lifted-overlap "$RECIPROCAL_OVERLAP" \
        --max-local-gene-distance "$LOCAL_GENE_DISTANCE" \
        --out "$out"


    # --------------------------------------------------------
    # Per-species output QC
    # --------------------------------------------------------

    n_rows="$(
        awk '
            NR > 1 && NF > 0 {n++}
            END {print n+0}
        ' "$out"
    )"

    if [[ "$n_rows" -ne "$N_REF_CRES" ]]; then

        echo "ERROR: $sp produced $n_rows rows;" >&2
        echo "expected $N_REF_CRES." >&2
        exit 1

    fi

    ((n_done += 1))

    echo

done < "$TARGET_SPECIES_FILE"


echo "Species classified: $n_done"
echo
echo "[2/4] PASSED"
echo


# ============================================================
# 3. Combine species-specific classifications
# ============================================================

echo "[3/4] Combining classification results..."
echo

python3 "$SCRIPTS_DIR/02_combine_turnover_results.py" \
    "$TURNOVER_BY_SPECIES_DIR" \
    "$TARGET_SPECIES_FILE" \
    "$CRE_TURNOVER_ALL_SPECIES" \
    "$CRE_TURNOVER_MATRIX" \
    "$SPECIES_SUMMARY" \
    "$EXPECTED_REFERENCE_CRES"


for file in \
    "$CRE_TURNOVER_ALL_SPECIES" \
    "$CRE_TURNOVER_MATRIX" \
    "$SPECIES_SUMMARY"
do

    [[ -s "$file" ]] || {
        echo "ERROR: expected combined result missing:" >&2
        echo "  $file" >&2
        exit 1
    }

done

echo
echo "[3/4] PASSED"
echo


# ============================================================
# 4. Build species-level QC summary
# ============================================================

echo "[4/4] Building species QC summary..."
echo

python3 "$SCRIPTS_DIR/03_build_species_qc_summary.py" \
    --predictions "$SO_ALL_SPECIES_FBGN" \
    --targets "$TARGET_SPECIES_FILE" \
    --turnover-dir "$TURNOVER_BY_SPECIES_DIR" \
    --lifted-dir "$LIFTED_CRES_DIR" \
    --out "$SPECIES_QC_SUMMARY" \
    --overlap-threshold "$RECIPROCAL_OVERLAP" \
    --min-mapping-rate "$QC_MIN_MAPPING_RATE" \
    --min-scrmshaw-peaks "$QC_MIN_SCRMSHAW_PEAKS" \
    --min-seqid-overlap "$QC_MIN_SEQID_OVERLAP" \
    --min-mapped-for-overlap-check "$QC_MIN_MAPPED_FOR_OVERLAP_CHECK" \
    --high-mapping-rate "$QC_HIGH_MAPPING_RATE"


[[ -s "$SPECIES_QC_SUMMARY" ]] || {

    echo "ERROR: species_qc_summary.tsv was not created." >&2
    exit 1

}

echo
echo "[4/4] PASSED"
echo


# ============================================================
# Final summary
# ============================================================

echo "============================================================"
echo "CRE classification complete"
echo "============================================================"
echo "Species classified: $n_done"
echo "Reference CREs     : $N_REF_CRES"
echo
echo "Main outputs:"
echo "  $CRE_TURNOVER_ALL_SPECIES"
echo "  $CRE_TURNOVER_MATRIX"
echo "  $SPECIES_SUMMARY"
echo "  $SPECIES_QC_SUMMARY"
echo
echo "Finished: $(date)"
echo "============================================================"
