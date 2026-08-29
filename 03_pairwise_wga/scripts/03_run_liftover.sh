#!/usr/bin/env bash
set -euo pipefail

WGA="$HOME/cre_turnover/project/pairwise_wga"
MAP="$HOME/cre_turnover/project/mapping_orthologs"

TARGETS="$WGA/target_species.txt"
REF_BED="$MAP/reference_cres/dmel_reference_cres.bed"

export PATH="$WGA/bin:$PATH"

echo "=================================================="
echo "Dmel pairwise WGA / liftOver pipeline"
echo "Started: $(date)"
echo "=================================================="
echo


# --------------------------------------------------
# 1. Check inputs
# --------------------------------------------------

echo "[1/3] Checking inputs..."

for f in \
    "$TARGETS" \
    "$REF_BED"
do
    [[ -s "$f" ]] || {
        echo "ERROR: missing required file:"
        echo "  $f"
        exit 1
    }
done

N_REF_CRES=$(grep -vcE '^(#|$)' "$REF_BED")

echo "Reference CREs: $N_REF_CRES"
echo


# --------------------------------------------------
# 2. Check all alignment chains
# --------------------------------------------------

echo "[2/3] Checking pairwise alignment chains..."

missing=0
valid=0
total=0

while read -r sp
do
    [[ -z "$sp" ]] && continue
    ((total+=1))

    chain="$WGA/alignments_dmel/$sp/d_melanogaster.${sp}.liftover.chain.gz"

    if [[ ! -s "$chain" ]]; then
        echo "MISSING  $sp"
        ((missing+=1))
        continue
    fi

    if ! gzip -t "$chain" 2>/dev/null; then
        echo "CORRUPT  $sp"
        ((missing+=1))
        continue
    fi

    echo "OK       $sp"
    ((valid+=1))

done < "$TARGETS"

echo
echo "Targets : $total"
echo "Valid   : $valid"
echo "Missing : $missing"
echo

if [[ "$missing" -ne 0 ]]; then
    echo "ERROR: Not all pairwise alignments are available."
    exit 1
fi


# --------------------------------------------------
# 3. liftOver + QC
# --------------------------------------------------

echo "[3/3] Lifting Dmel reference CREs..."

qc_fail=0

printf "%-20s %8s %10s %8s\n" \
    "species" "mapped" "unmapped" "total"

while read -r sp
do
    [[ -z "$sp" ]] && continue

    chain="$WGA/alignments_dmel/$sp/d_melanogaster.${sp}.liftover.chain.gz"
    outdir="$WGA/lifted_cres_dmel/$sp"

    mapped="$outdir/dmel_reference_cres.${sp}.bed"
    unmapped="$outdir/dmel_reference_cres.${sp}.unmapped.bed"

    mkdir -p "$outdir"

    liftOver \
        -minMatch=0.50 \
        "$REF_BED" \
        "$chain" \
        "$mapped" \
        "$unmapped"

    n_mapped=$(grep -vcE '^(#|$)' "$mapped" || true)
    n_unmapped=$(grep -vcE '^(#|$)' "$unmapped" || true)

    total_cre=$((n_mapped + n_unmapped))

    printf "%-20s %8d %10d %8d\n" \
        "$sp" "$n_mapped" "$n_unmapped" "$total_cre"

    if [[ "$total_cre" -ne "$N_REF_CRES" ]]; then
        echo "WARNING: $sp total=$total_cre expected=$N_REF_CRES"
        ((qc_fail+=1))
    fi

    duplicates=$(
        cut -f4 "$mapped" \
        | sort \
        | uniq -d \
        | wc -l
    )

    if [[ "$duplicates" -ne 0 ]]; then
        echo "WARNING: $sp has duplicated mapped CRE IDs"
        ((qc_fail+=1))
    fi

done < "$TARGETS"

echo

if [[ "$qc_fail" -ne 0 ]]; then
    echo "ERROR: liftOver QC failed."
    exit 1
fi

echo "=================================================="
echo "WGA / liftOver pipeline complete"
echo "Finished: $(date)"
echo "=================================================="
