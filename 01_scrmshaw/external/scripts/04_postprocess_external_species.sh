#!/usr/bin/env bash
set -euo pipefail


# ============================================================
# 04 - Postprocess one external SCRMshaw species
#
# Usage:
#   bash scripts/04_postprocess_external_species.sh <species_slug>
#
# Example:
#   bash scripts/04_postprocess_external_species.sh d_affinis
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "$ROOT/config/external_config.sh"

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$POST_ENV"

# MACS2 compatibility wrapper
export PATH="$OG_PIPELINE_ROOT/bin:$PATH"


# ------------------------------------------------------------
# Species argument
# ------------------------------------------------------------

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 SPECIES_SLUG" >&2
    echo "Example:" >&2
    echo "  $0 d_affinis" >&2
    exit 1
fi

slug="$1"


# ------------------------------------------------------------
# Find species in external manifest
# ------------------------------------------------------------

row="$(
    awk -F'\t' -v slug="$slug" '
        NR > 1 && $2 == slug {
            print
            exit
        }
    ' "$EXTERNAL_MANIFEST"
)"

if [[ -z "$row" ]]; then
    echo "FEHLER: Species '$slug' nicht im Manifest gefunden:" >&2
    echo "  $EXTERNAL_MANIFEST" >&2
    exit 1
fi


# Manifest columns:
# index | slug | species | bed | gff

IFS=$'\t' read -r index manifest_slug species bed gff <<< "$row"


# ------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------

if [[ "$manifest_slug" != "$slug" ]]; then
    echo "FEHLER: Manifest-Slug stimmt nicht überein." >&2
    echo "  requested: $slug" >&2
    echo "  manifest : $manifest_slug" >&2
    exit 1
fi

[[ -s "$bed" ]] || {
    echo "FEHLER: BED fehlt oder ist leer:" >&2
    echo "  $bed" >&2
    exit 1
}

[[ -s "$gff" ]] || {
    echo "FEHLER: GFF fehlt oder ist leer:" >&2
    echo "  $gff" >&2
    exit 1
}


# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------

run="$EXTERNAL_RUNS_DIR/$slug"
result="$EXTERNAL_RESULTS_DIR/$slug"

mkdir -p "$run" "$result"

cd "$run"


# ------------------------------------------------------------
# Clean previous postprocessing outputs for this species
# ------------------------------------------------------------

rm -f \
    scrmshawOutput_offset_0to240.bed \
    scrmshawOutput_peaksCalled_* \
    peaks_AllSets.bed \
    sumScoredsorted_*

rm -rf tmp


# ------------------------------------------------------------
# Link filtered external SCRMshaw input
# ------------------------------------------------------------

ln -s "$bed" scrmshawOutput_offset_0to240.bed


# ------------------------------------------------------------
# Run SCRMshaw postprocessing
# ------------------------------------------------------------

echo "============================================================"
echo "External SCRMshaw postprocessing"
echo "============================================================"
echo "Slug    : $slug"
echo "Species : $species"
echo "BED     : $bed"
echo "GFF     : $gff"
echo "Run dir : $run"
echo


python "$POSTPROC_ROOT/postProcessingScrmshawPipeline.py" \
    -num 5000 \
    -topN Median \
    -so scrmshawOutput_offset_0to240.bed \
    -gff "$gff"


# ------------------------------------------------------------
# Collect peak files
# ------------------------------------------------------------

shopt -s nullglob

peaks=(scrmshawOutput_peaksCalled_*)

if (( ${#peaks[@]} == 0 )); then
    echo "FEHLER: Keine Peaks für $slug erzeugt." >&2
    exit 1
fi


cat "${peaks[@]}" > peaks_AllSets.bed


[[ -s peaks_AllSets.bed ]] || {
    echo "FEHLER: peaks_AllSets.bed ist leer für $slug." >&2
    exit 1
}


# ------------------------------------------------------------
# Copy final results
# ------------------------------------------------------------

cp peaks_AllSets.bed "$result/"
cp "${peaks[@]}" "$result/"


# ------------------------------------------------------------
# QC / provenance
# ------------------------------------------------------------

wc -l peaks_AllSets.bed > "$result/peak_count.txt"

printf 'finished\n' > POSTPROCESS_FINISHED.ok


n_peaks=$(wc -l < peaks_AllSets.bed)


echo
echo "Externes Postprocessing abgeschlossen: $slug"
echo "Final peaks: $n_peaks"
echo "Result dir:"
echo "  $result"
