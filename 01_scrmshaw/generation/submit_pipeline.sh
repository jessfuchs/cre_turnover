#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/config/config.sh"
cd "$PIPELINE_ROOT"

mkdir -p "$LOG_DIR" "$DATA_DIR" "$RUNS_DIR" "$RESULTS_DIR"

n_species="$(grep -Ev '^[[:space:]]*(#|$)' "$SPECIES_FILE" | wc -l)"
if (( n_species == 0 )); then
    echo "FEHLER: $SPECIES_FILE enthält keine Spezies." >&2
    exit 1
fi

for required in \
    "$SCRM_ROOT/code/scrm.pl" \
    "$UTILITY_ROOT/preflight-scrmshaw.pl" \
    "$UTILITY_ROOT/Generate_top_N_SCRMhits.pl" \
    "$POSTPROC_ROOT/postProcessingScrmshawPipeline.py" \
    "$TRAINING_SET/crms.fasta" \
    "$TRAINING_SET/neg.fasta" \
    "$TRF_BIN"; do
    [[ -e "$required" ]] || {
        echo "FEHLER: $required fehlt. Zuerst scripts/setup_software.sh ausführen." >&2
        exit 1
    }
done

printf '%s\n' "$(realpath "$TRAINING_SET")" > "$TRAINING_LIST"

common=(--nodelist=abacus-2)
[[ -n "$SLURM_PARTITION" ]] && common+=(--partition="$SLURM_PARTITION")
[[ -n "$SLURM_QOS" ]] && common+=(--qos="$SLURM_QOS")
[[ -n "$SLURM_ACCOUNT" ]] && common+=(--account="$SLURM_ACCOUNT")

jid_download="$(
    sbatch --parsable "${common[@]}" \
        --time="$DOWNLOAD_TIME" --mem="$DOWNLOAD_MEM" --cpus-per-task="$DOWNLOAD_CPUS" \
        slurm/01_download_extract.slurm
)"

jid_prep="$(
    sbatch --parsable "${common[@]}" \
        --dependency="afterok:${jid_download}" \
        --array="0-$((n_species-1))%${MAX_PREP_JOBS}" \
        --time="$PREP_TIME" --mem="$PREP_MEM" --cpus-per-task="$PREP_CPUS" \
        slurm/02_prepare_species.slurm
)"

n_offsets=$((n_species * 25))
jid_offsets="$(
    sbatch --parsable "${common[@]}" \
        --dependency="afterok:${jid_prep}" \
        --array="0-$((n_offsets-1))%${MAX_OFFSET_JOBS}" \
        --time="$OFFSET_TIME" --mem="$OFFSET_MEM" --cpus-per-task="$OFFSET_CPUS" \
        slurm/03_scrmshaw_offsets.slurm
)"

jid_post="$(
    sbatch --parsable "${common[@]}" \
        --dependency="afterok:${jid_offsets}" \
        --array="0-$((n_species-1))%${MAX_POST_JOBS}" \
        --time="$POST_TIME" --mem="$POST_MEM" --cpus-per-task="$POST_CPUS" \
        slurm/04_postprocess.slurm
)"

jid_summary="$(
    sbatch --parsable "${common[@]}" \
        --dependency="afterok:${jid_post}" \
        --time="$SUMMARY_TIME" --mem="$SUMMARY_MEM" \
        slurm/05_summary.slurm
)"

cat <<EOF
Pipeline eingereicht:
  Download/Extraktion : $jid_download
  Preprocessing       : $jid_prep
  SCRMshaw-Array      : $jid_offsets  (${n_offsets} Tasks)
  Postprocessing      : $jid_post
  Zusammenfassung     : $jid_summary

Status:
  squeue -u "$USER"

Nach Abschluss:
  column -t -s \$'\t' results/summary.tsv
EOF
