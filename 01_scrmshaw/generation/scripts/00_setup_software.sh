#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

mkdir -p "$SOFTWARE_DIR" "$LOG_DIR" "$(dirname "$TRAINING_LIST")"

clone_or_update() {
    local url="$1"
    local dest="$2"
    if [[ -d "$dest/.git" ]]; then
        echo "Updating $dest"
        git -C "$dest" pull --ff-only
    elif [[ -e "$dest" ]]; then
        echo "ERROR: $dest exists but is not a Git repository." >&2
        exit 1
    else
        git clone "$url" "$dest"
    fi
}

clone_or_update https://github.com/HalfonLab/SCRMshaw_HD.git "$SCRM_ROOT"
clone_or_update https://github.com/HalfonLab/UtilityPrograms.git "$UTILITY_ROOT"
clone_or_update https://github.com/HalfonLab/post_processing_SCRMshaw_pipeline.git "$POSTPROC_ROOT"
clone_or_update https://github.com/HalfonLab/dmel_training_sets.git "$TRAINING_ROOT"

if [[ ! -x "$TRF_BIN" ]]; then
    curl -L --fail --retry 5 \
      -o "$TRF_BIN" \
      https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
    chmod +x "$TRF_BIN"
fi

for f in crms.fasta neg.fasta; do
    [[ -s "$TRAINING_SET/$f" ]] || {
        echo "ERROR: Training file is missing: $TRAINING_SET/$f" >&2
        exit 1
    }
done

printf '%s\n' "$(realpath "$TRAINING_SET")" > "$TRAINING_LIST"
chmod +x "$PIPELINE_ROOT/bin/macs2"
find "$PIPELINE_ROOT/scripts" "$PIPELINE_ROOT/slurm" -type f -name '*.sh' -o -name '*.slurm' | xargs chmod +x

echo
echo "Software and training set are prepared."
echo "Training list: $TRAINING_LIST"
