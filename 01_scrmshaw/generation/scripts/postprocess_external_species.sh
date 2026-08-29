#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../config/config.sh"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$POST_ENV"

export PATH="$PIPELINE_ROOT/bin:$PATH"

manifest="$PIPELINE_ROOT/data/external_manifest.tsv"
species_idx="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID fehlt}"

row="$(
    awk -F'\t' -v idx="$species_idx" '
        NR > 1 && $1 == idx {
            print
            exit
        }
    ' "$manifest"
)"

[[ -n "$row" ]] || {
    echo "Kein Eintrag für Index $species_idx in $manifest" >&2
    exit 1
}

IFS=$'\t' read -r idx slug species bed gff <<< "$row"

[[ -s "$bed" ]] || {
    echo "BED fehlt oder ist leer: $bed" >&2
    exit 1
}

[[ -s "$gff" ]] || {
    echo "GFF fehlt oder ist leer: $gff" >&2
    exit 1
}

run_dir="$PIPELINE_ROOT/external_runs/$slug"
result_dir="$PIPELINE_ROOT/external_results/$slug"

mkdir -p "$run_dir" "$result_dir"
cd "$run_dir"

rm -f scrmshawOutput_offset_0to240.bed
rm -f scrmshawOutput_peaksCalled_*
rm -f peaks_AllSets.bed
rm -f sumScoredsorted_*
rm -rf tmp

ln -s "$bed" scrmshawOutput_offset_0to240.bed

python "$POSTPROC_ROOT/postProcessingScrmshawPipeline.py" \
    -num 5000 \
    -topN Median \
    -so scrmshawOutput_offset_0to240.bed \
    -gff "$gff"

shopt -s nullglob
peak_files=(scrmshawOutput_peaksCalled_*)

if (( ${#peak_files[@]} == 0 )); then
    echo "Keine Peak-Dateien erzeugt: $slug" >&2
    exit 1
fi

cat "${peak_files[@]}" > peaks_AllSets.bed

[[ -s peaks_AllSets.bed ]] || {
    echo "peaks_AllSets.bed ist leer: $slug" >&2
    exit 1
}

cp peaks_AllSets.bed "$result_dir/"
cp "${peak_files[@]}" "$result_dir/"

wc -l peaks_AllSets.bed > "$result_dir/peak_count.txt"
printf 'finished\n' > POSTPROCESS_FINISHED.ok

echo "Externes Postprocessing abgeschlossen: $slug"
