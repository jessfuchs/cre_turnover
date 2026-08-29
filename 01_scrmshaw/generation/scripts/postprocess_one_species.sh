#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$POST_ENV"
export PATH="$PIPELINE_ROOT/bin:$PATH"

species_idx="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID fehlt}"
row="$(awk -F'\t' -v i="$species_idx" 'NR>1 && $1==i {print; exit}' "$MANIFEST")"
[[ -n "$row" ]] || { echo "Kein Manifest-Eintrag für Index $species_idx" >&2; exit 1; }
IFS=$'\t' read -r idx slug species genome annotation annotation_format <<< "$row"

run_dir="$RUNS_DIR/$slug"
gff="$(cat "$run_dir/gff.path")"
cd "$run_dir"

for offset_idx in $(seq 0 24); do
    lb=$((offset_idx * 10))
    task_num=$((offset_idx + 1))
    test -f "task_offset_${lb}_${task_num}/SCRM_FINISHED.ok" || {
        echo "Fehlender Offset: task_offset_${lb}_${task_num}" >&2
        exit 1
    }
done

# Reruns sicher machen; keine task_offset-Ordner löschen.
rm -f scrmshawOutput_offset_0to240.bed peaks_AllSets.bed
rm -f scrmshawOutput_peaksCalled_* sumScoredsorted_*
rm -rf tmp scrmsIndividualHits_0to240offset
mkdir -p scrmsIndividualHits_0to240offset

individual=()
for offset_idx in $(seq 0 24); do
    lb=$((offset_idx * 10))
    task_num=$((offset_idx + 1))
    outfile="scrmshawOutput_offset_${lb}.5000scrms"
    perl "$UTILITY_ROOT/Generate_top_N_SCRMhits.pl" \
        -d "task_offset_${lb}_${task_num}" \
        -n 5000 \
        -o "$outfile"
    test -s "$outfile"
    individual+=("$outfile")
done

cat "${individual[@]}" > scrmshawOutput_offset_0to240.bed
mv "${individual[@]}" scrmsIndividualHits_0to240offset/

python "$POSTPROC_ROOT/postProcessingScrmshawPipeline.py" \
    -num 5000 \
    -topN Median \
    -so scrmshawOutput_offset_0to240.bed \
    -gff "$gff"

shopt -s nullglob
peak_files=(scrmshawOutput_peaksCalled_*)
(( ${#peak_files[@]} > 0 )) || {
    echo "Postprocessing erzeugte keine Peak-Datei." >&2
    exit 1
}
cat "${peak_files[@]}" > peaks_AllSets.bed
test -s peaks_AllSets.bed
rm -rf tmp

mkdir -p "$RESULTS_DIR/$slug"
cp peaks_AllSets.bed "$RESULTS_DIR/$slug/"
cp "${peak_files[@]}" "$RESULTS_DIR/$slug/"
wc -l peaks_AllSets.bed > "$RESULTS_DIR/$slug/peak_count.txt"
touch POSTPROCESS_FINISHED.ok

echo "Postprocessing abgeschlossen: $slug"
