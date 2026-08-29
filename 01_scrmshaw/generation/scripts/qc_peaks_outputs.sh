#!/usr/bin/env bash
set -euo pipefail

# Default:
#   original pipeline results
#
# Optional:
#   bash qc_peaks_outputs.sh /pfad/zu/results

RESULTS_DIR="${1:-$HOME/cre_turnover/project/scrmshaw_pipeline/results}"

[[ -d "$RESULTS_DIR" ]] || {
    echo "FEHLER: Results-Verzeichnis nicht gefunden: $RESULTS_DIR" >&2
    exit 1
}

printf "species\tpeaks_allsets_exists\tpeaks_allsets_lines\tn_peak_files\tn_empty_peak_files\tpeak_lines_sum\tstatus\n"

for dir in "$RESULTS_DIR"/*; do
    [[ -d "$dir" ]] || continue

    species="$(basename "$dir")"

    peaks_allsets="$dir/peaks_AllSets.bed"

    allsets_exists="NO"
    allsets_lines=0

    if [[ -s "$peaks_allsets" ]]; then
        allsets_exists="YES"
        allsets_lines="$(wc -l < "$peaks_allsets")"
    fi

    shopt -s nullglob
    peak_files=("$dir"/scrmshawOutput_peaksCalled_*)
    shopt -u nullglob

    n_peak_files=${#peak_files[@]}
    n_empty_peak_files=0
    peak_sum=0

    if (( n_peak_files > 0 )); then
        for f in "${peak_files[@]}"; do
            if [[ ! -s "$f" ]]; then
                n_empty_peak_files=$((n_empty_peak_files + 1))
                continue
            fi

            n="$(wc -l < "$f")"
            peak_sum=$((peak_sum + n))
        done
    fi

    status="OK"

    if [[ "$allsets_exists" != "YES" ]]; then
        status="NO_PEAKS_ALLSETS"
    elif (( n_peak_files == 0 )); then
        status="NO_PEAK_FILES"
    elif (( n_empty_peak_files > 0 )); then
        status="EMPTY_PEAK_FILE"
    elif (( peak_sum != allsets_lines )); then
        status="LINECOUNT_MISMATCH"
    fi

    printf "%s\t%s\t%d\t%d\t%d\t%d\t%s\n" \
        "$species" \
        "$allsets_exists" \
        "$allsets_lines" \
        "$n_peak_files" \
        "$n_empty_peak_files" \
        "$peak_sum" \
        "$status"
done
