#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$SCRM_ENV"

task_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID fehlt}"
row="$(awk -F'\t' -v i="$task_id" 'NR>1 && $1==i {print; exit}' "$MANIFEST")"
[[ -n "$row" ]] || { echo "Kein Manifest-Eintrag für Index $task_id" >&2; exit 1; }

IFS=$'\t' read -r idx slug species genome annotation annotation_format <<< "$row"
run_dir="$RUNS_DIR/$slug"
input_dir="$run_dir/input"
mkdir -p "$input_dir" "$run_dir"

ln -sfn "$genome" "$input_dir/source_genome.fa"

gff="$input_dir/annotation.gff3"
if [[ "$annotation_format" == "gtf" ]]; then
    echo "Konvertiere GTF nach GFF3: $annotation"
    rm -f "$gff"
    agat_convert_sp_gxf2gxf.pl --gff "$annotation" -o "$gff"
else
    ln -sfn "$annotation" "$gff"
fi

echo "Preflight für $species"
perl "$UTILITY_ROOT/preflight-scrmshaw.pl" \
    --gff "$gff" \
    --fasta "$genome" \
    --suppress_log \
    > "$run_dir/preflight.log" 2>&1

grep -q "GFF3: number of genes:" "$run_dir/preflight.log"
grep -q "GFF3: number of exons:" "$run_dir/preflight.log"
grep -q "FASTA: all sequences have proper characters" "$run_dir/preflight.log"
grep -q "All seqids in GFF are also in FASTA" "$run_dir/preflight.log"

mapped="$input_dir/mappedOnly_${slug}.fa"
python "$PIPELINE_ROOT/scripts/filter_fasta_by_gff.py" \
    --gff "$gff" \
    --fasta "$genome" \
    --output "$mapped" \
    > "$run_dir/filter_fasta.log" 2>&1

if [[ "$RUN_TRF" == "1" ]]; then
    echo "TRF-Masking für $species"

    set +e
    (
        cd "$input_dir"
        "$TRF_BIN" "$(basename "$mapped")" \
            2 7 7 80 10 50 500 -m -h
    )
    trf_rc=$?
    set -e

    masked="$(
        find "$input_dir" \
            -maxdepth 1 \
            -type f \
            -name 'mappedOnly_*.mask' \
            -print -quit
    )"

    if [[ -z "$masked" || ! -s "$masked" ]]; then
        echo "FEHLER: TRF erzeugte keine verwendbare .mask-Datei." >&2
        echo "TRF Exit-Code: $trf_rc" >&2
        exit 1
    fi

    if [[ "$trf_rc" -ne 0 ]]; then
        echo "WARNUNG: TRF Exit-Code $trf_rc, aber Mask-Datei existiert."
    fi
else
    masked="$mapped"
fi

printf '%s\n' "$(realpath "$masked")" > "$run_dir/masked_genome.path"
printf '%s\n' "$(realpath "$gff")" > "$run_dir/gff.path"
touch "$run_dir/PREPARED.ok"

echo "Vorbereitung abgeschlossen: $slug"
