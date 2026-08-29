#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$SCRM_ENV"

global_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID fehlt}"
species_idx=$(( global_id / 25 ))
offset_idx=$(( global_id % 25 ))
lb=$(( offset_idx * 10 ))
task_num=$(( offset_idx + 1 ))

row="$(awk -F'\t' -v i="$species_idx" 'NR>1 && $1==i {print; exit}' "$MANIFEST")"
[[ -n "$row" ]] || { echo "Kein Manifest-Eintrag für Speziesindex $species_idx" >&2; exit 1; }
IFS=$'\t' read -r idx slug species genome annotation annotation_format <<< "$row"

run_dir="$RUNS_DIR/$slug"
[[ -f "$run_dir/PREPARED.ok" ]] || { echo "$slug ist nicht vorbereitet." >&2; exit 1; }
masked="$(cat "$run_dir/masked_genome.path")"
gff="$(cat "$run_dir/gff.path")"
outdir="$run_dir/task_offset_${lb}_${task_num}"
done_file="$outdir/SCRM_FINISHED.ok"

if [[ "$RESUME" == "1" && -s "$done_file" ]]; then
    echo "Bereits fertig: $slug offset=$lb"
    exit 0
fi

rm -rf "$outdir"
mkdir -p "$outdir"

read -r -a scoring <<< "$SCORING_FLAGS"

echo "Spezies: $species"
echo "Offset: $lb"
echo "Output: $outdir"

perl "$SCRM_ROOT/code/scrm.pl" \
    --thitw "$THITW" \
    --gff "$gff" \
    --genome "$masked" \
    --traindirlst "$TRAINING_LIST" \
    "${scoring[@]}" \
    --lb "$lb" \
    --step 123 \
    --outdir "$outdir"

# Mindestens IMM wird in der Standardkonfiguration erwartet.
if [[ " $SCORING_FLAGS " == *" --imm "* ]]; then
    test -s "$outdir/hits/imm/adult_muscle/adult_muscle.hits.ranked"
fi
if [[ " $SCORING_FLAGS " == *" --hexmcd "* ]]; then
    test -s "$outdir/hits/hexmcd/adult_muscle/adult_muscle.hits.ranked"
fi
if [[ " $SCORING_FLAGS " == *" --pac "* ]]; then
    test -s "$outdir/hits/pac/adult_muscle/adult_muscle.hits.ranked"
fi

touch "$done_file"
echo "Fertig: $slug offset=$lb"
