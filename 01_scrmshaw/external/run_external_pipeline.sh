#!/usr/bin/env bash
set -euo pipefail


# ======================================================================
# Paths
# ======================================================================

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

source "$ROOT/config/external_config.sh"

SPECIES_FILE="$ROOT/species_external.txt"


# ======================================================================
# Preconditions
# ======================================================================

echo "============================================================"
echo "External SCRMshaw import pipeline"
echo "============================================================"
echo "Root           : $ROOT"
echo "Species file   : $SPECIES_FILE"
echo "External BEDs  : $EXTERNAL_BED_DIR"
echo "External GFFs  : $EXTERNAL_GFF_DIR"
echo "Training set   : $TRAINING_SET"
echo "Method         : $METHOD"
echo


[[ -s "$SPECIES_FILE" ]] || {
    echo "ERROR: Species file is missing or empty:" >&2
    echo "  $SPECIES_FILE" >&2
    exit 1
}

[[ -d "$EXTERNAL_BED_DIR" ]] || {
    echo "ERROR: External BED directory is missing:" >&2
    echo "  $EXTERNAL_BED_DIR" >&2
    exit 1
}

[[ -d "$SOURCE_GFF_ROOT" ]] || {
    echo "ERROR: selected directory is missing:" >&2
    echo "  $SOURCE_GFF_ROOT" >&2
    exit 1
}

mkdir -p \
    "$ROOT/logs" \
    "$EXTERNAL_GFF_DIR" \
    "$FILTERED_BED_DIR" \
    "$EXTERNAL_RUNS_DIR" \
    "$EXTERNAL_RESULTS_DIR"


# ======================================================================
# Count species
# ======================================================================

N_SPECIES="$(
    awk '
        /^[[:space:]]*#/ {next}
        NF == 0 {next}
        {n++}
        END {print n+0}
    ' "$SPECIES_FILE"
)"


if [[ "$N_SPECIES" -lt 1 ]]; then
    echo "ERROR: No species found in $SPECIES_FILE." >&2
    exit 1
fi


echo "Species listed : $N_SPECIES"
echo


# ======================================================================
# Step 1: Extract GFF annotations
# ======================================================================

echo "============================================================"
echo "[1/5] Extracting external GFF files"
echo "============================================================"

bash "$ROOT/scripts/01_extract_external_gffs.sh" \
    "$SPECIES_FILE"

echo


# ======================================================================
# Step 2: Filter external BED files
# ======================================================================

echo "============================================================"
echo "[2/5] Filtering external SCRMshaw BEDs"
echo "      training = $TRAINING_SET"
echo "      method   = $METHOD"
echo "============================================================"

bash "$ROOT/scripts/02_filter_external_beds.sh"

echo


# ======================================================================
# Step 3: Validate BED/GFF/rank structure and create manifest
# ======================================================================

echo "============================================================"
echo "[3/5] Validating external inputs"
echo "============================================================"

python "$ROOT/scripts/03_validate_bed_gff_and_manifest.py" \
    --bed-dir "$FILTERED_BED_DIR" \
    --gff-dir "$EXTERNAL_GFF_DIR" \
    --species-file "$SPECIES_FILE" \
    --manifest "$EXTERNAL_MANIFEST" \
    --training "$TRAINING_SET" \
    --method "$METHOD" \
    --expected-offsets 25 \
    --max-rank 5000


[[ -s "$EXTERNAL_MANIFEST" ]] || {
    echo "ERROR: Manifest was not created:" >&2
    echo "  $EXTERNAL_MANIFEST" >&2
    exit 1
}


N_MANIFEST="$(
    awk '
        NR > 1 && NF > 0 {n++}
        END {print n+0}
    ' "$EXTERNAL_MANIFEST"
)"


if [[ "$N_MANIFEST" -ne "$N_SPECIES" ]]; then
    echo "ERROR: Species count does not match." >&2
    echo "  species_external.txt : $N_SPECIES" >&2
    echo "  external manifest    : $N_MANIFEST" >&2
    exit 1
fi


echo
echo "Manifest successfully created:"
echo "  $EXTERNAL_MANIFEST"
echo "Species: $N_MANIFEST"
echo


# ======================================================================
# Step 4: Run external SCRMshaw postprocessing directly
# ======================================================================

echo "============================================================"
echo "[4/5] Running external postprocessing"
echo "============================================================"
echo

N_DONE=0

while IFS=$'\t' read -r index slug species bed gff; do

    [[ "$index" == "index" ]] && continue
    [[ -n "$slug" ]] || continue

    echo "------------------------------------------------------------"
    echo "Postprocessing: $slug"
    echo "Species       : $species"
    echo "------------------------------------------------------------"

    bash "$ROOT/scripts/04_postprocess_external_species.sh" "$slug"

    N_DONE=$((N_DONE + 1))

    echo

done < "$EXTERNAL_MANIFEST"


if [[ "$N_DONE" -ne "$N_MANIFEST" ]]; then
    echo "ERROR: Not all external species were processed." >&2
    echo "  expected : $N_MANIFEST" >&2
    echo "  processed: $N_DONE" >&2
    exit 1
fi


echo "External postprocessing completed."
echo "Species processed: $N_DONE"
echo


# ======================================================================
# Step 5: Build combined results
# ======================================================================

echo "============================================================"
echo "[5/5] Building combined results"
echo "============================================================"
echo


python "$ROOT/scripts/05_build_combined_results.py" \
    --generated-manifest "$GENERATED_MANIFEST" \
    --external-manifest "$EXTERNAL_MANIFEST" \
    --generated-results "$GENERATED_RESULTS_DIR" \
    --external-results "$EXTERNAL_RESULTS_DIR" \
    --combined-manifest "$COMBINED_MANIFEST" \
    --combined-results "$COMBINED_RESULTS_DIR"


echo
echo "Running combined summary..."
echo

bash "$ROOT/scripts/06_run_combined_summary.sh"


# ======================================================================
# Final summary
# ======================================================================

echo
echo "============================================================"
echo "External SCRMshaw import pipeline COMPLETE"
echo "============================================================"
echo
echo "External species : $N_MANIFEST"
echo "Training set     : $TRAINING_SET"
echo "Method           : $METHOD"
echo
echo "External results : $EXTERNAL_RESULTS_DIR"
echo "Combined manifest:"
echo "  $COMBINED_MANIFEST"
echo "Combined results:"
echo "  $COMBINED_RESULTS_DIR"
echo
echo "Finished: $(date)"
echo "============================================================"
