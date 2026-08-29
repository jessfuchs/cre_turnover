#!/usr/bin/env bash
set -euo pipefail

PROJECT="$HOME/cre_turnover/project"

WGA="$PROJECT/pairwise_wga"
MAP="$PROJECT/mapping_orthologs"
CLASS="$PROJECT/cre_classification"

TARGETS="$WGA/target_species.txt"

REF_TSV="$MAP/reference_cres/dmel_reference_cres.tsv"
PREDICTIONS="$MAP/ortholog_results/SO_all_species_fbgn.tsv"

OUTDIR="$CLASS/results/turnover_by_species"


echo "=================================================="
echo "CRE-state classification pipeline"
echo "Started: $(date)"
echo "=================================================="
echo


# --------------------------------------------------
# 1. Check required inputs
# --------------------------------------------------

echo "[1/3] Checking classification inputs..."

for f in \
    "$TARGETS" \
    "$REF_TSV" \
    "$PREDICTIONS" \
    "$CLASS/scripts/01_classify_target_cres.py" \
    "$CLASS/scripts/02_combine_turnover_results.py"
do
    [[ -s "$f" ]] || {
        echo "ERROR: missing required file:"
        echo "  $f"
        exit 1
    }
done

N_REF_CRES=$(
    awk 'NR > 1 && NF > 0 {n++} END {print n+0}' \
        "$REF_TSV"
)

echo "Reference CREs: $N_REF_CRES"
echo


# --------------------------------------------------
# 2. Check lifted CREs and classify
# --------------------------------------------------

echo "[2/3] Classifying all target species..."

mkdir -p "$OUTDIR"

n_done=0

while read -r sp
do
    [[ -z "$sp" ]] && continue

    lifted="$WGA/lifted_cres_dmel/$sp/dmel_reference_cres.${sp}.bed"

    if [[ ! -f "$lifted" ]]; then
        echo "ERROR: lifted BED missing for $sp:"
        echo "  $lifted"
        exit 1
    fi

    out="$OUTDIR/dmel_to_${sp}_cre_turnover.tsv"

    echo
    echo "===== $sp ====="

    python3 "$CLASS/scripts/01_classify_target_cres.py" \
        --species "$sp" \
        --reference "$REF_TSV" \
        --lifted "$lifted" \
        --predictions "$PREDICTIONS" \
        --min-lifted-overlap 0.50 \
        --max-local-gene-distance 24000 \
        --out "$out"

    n_rows=$(
        awk 'NR > 1 && NF > 0 {n++} END {print n+0}' \
            "$out"
    )

    if [[ "$n_rows" -ne "$N_REF_CRES" ]]; then
        echo "ERROR: $sp produced $n_rows rows; expected $N_REF_CRES."
        exit 1
    fi

    ((n_done+=1))

done < "$TARGETS"

echo
echo "Species classified: $n_done"
echo


# --------------------------------------------------
# 3. Combine results
# --------------------------------------------------

echo "[3/3] Combining classification results..."

python3 "$CLASS/scripts/02_combine_turnover_results.py"

echo
echo "Combined results:"
echo "  $CLASS/results/cre_turnover_all_species.tsv"
echo "  $CLASS/results/cre_turnover_matrix.tsv"
echo "  $CLASS/results/species_summary.tsv"
echo

echo "=================================================="
echo "CRE classification complete"
echo "Finished: $(date)"
echo "=================================================="


# --------------------------------------------------
# 4. Build species QC summary
# --------------------------------------------------

echo "[4/4] Building species QC summary..."
echo

python3 "$CLASS/scripts/03_build_species_qc_summary.py"

[[ -s "$CLASS/results/species_qc_summary.tsv" ]] || {
    echo "ERROR: species_qc_summary.tsv was not created."
    exit 1
}

echo
echo "Species QC summary:"
echo "  $CLASS/results/species_qc_summary.tsv"
echo
