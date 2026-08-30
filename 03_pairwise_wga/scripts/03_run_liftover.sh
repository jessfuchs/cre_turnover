#!/usr/bin/env bash
set -euo pipefail


# ============================================================
# 03 - Lift D. melanogaster reference CREs to target genomes
#
# Purpose:
#   Project D. melanogaster reference CRE coordinates to each
#   target genome using the synteny-filtered pairwise chain files.
#
# Input:
#   - target_species.txt
#   - dmel_reference_cres.bed
#   - alignments_dmel/<species>/*.liftover.chain.gz
#
# Output:
#   - lifted_cres_dmel/<species>/*.bed
#   - lifted_cres_dmel/<species>/*.unmapped.bed
#
# Configuration:
#   config/wga_config.sh
# ============================================================


# ============================================================
# Configuration
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

source "$PIPELINE_ROOT/config/wga_config.sh"

export PATH="$BIN_DIR:$PATH"


# ============================================================
# Header
# ============================================================

echo "============================================================"
echo "D. melanogaster reference CRE liftOver"
echo "============================================================"
echo "Started    : $(date)"
echo "Reference  : $REFERENCE_SPECIES"
echo "minMatch   : $LIFTOVER_MIN_MATCH"
echo


# ============================================================
# 1. Check inputs
# ============================================================

echo "[1/3] Checking inputs..."

for file in \
    "$TARGET_SPECIES_FILE" \
    "$REFERENCE_CRES_BED"
do
    [[ -s "$file" ]] || {
        echo "ERROR: missing required file:" >&2
        echo "  $file" >&2
        exit 1
    }
done

command -v liftOver >/dev/null 2>&1 || {
    echo "ERROR: liftOver not found." >&2
    exit 1
}

N_REF_CRES="$(
    grep -vcE '^(#|$)' "$REFERENCE_CRES_BED"
)"

echo "Reference CREs: $N_REF_CRES"
echo


# ============================================================
# 2. Check pairwise alignment chains
# ============================================================

echo "[2/3] Checking pairwise alignment chains..."

missing=0
valid=0
total=0

while IFS= read -r sp; do

    [[ -n "$sp" ]] || continue

    ((total += 1))

    chain="$ALIGNMENT_DIR/$sp/${REFERENCE_SPECIES}.${sp}.liftover.chain.gz"

    if [[ ! -s "$chain" ]]; then

        echo "MISSING  $sp"
        ((missing += 1))
        continue

    fi

    if ! gzip -t "$chain" 2>/dev/null; then

        echo "CORRUPT  $sp"
        ((missing += 1))
        continue

    fi

    echo "OK       $sp"
    ((valid += 1))

done < "$TARGET_SPECIES_FILE"


echo
echo "Targets : $total"
echo "Valid   : $valid"
echo "Missing : $missing"
echo

if [[ "$missing" -ne 0 ]]; then
    echo "ERROR: not all pairwise alignment chains are available." >&2
    exit 1
fi


# ============================================================
# 3. liftOver and QC
# ============================================================

echo "[3/3] Lifting D. melanogaster reference CREs..."
echo

qc_fail=0

printf "%-20s %8s %10s %8s\n" \
    "species" "mapped" "unmapped" "total"


while IFS= read -r sp; do

    [[ -n "$sp" ]] || continue

    chain="$ALIGNMENT_DIR/$sp/${REFERENCE_SPECIES}.${sp}.liftover.chain.gz"

    outdir="$LIFTED_CRES_DIR/$sp"

    mapped="$outdir/dmel_reference_cres.${sp}.bed"
    unmapped="$outdir/dmel_reference_cres.${sp}.unmapped.bed"

    mkdir -p "$outdir"


    # --------------------------------------------------------
    # Project reference CRE coordinates
    # --------------------------------------------------------

    liftOver \
        -minMatch="$LIFTOVER_MIN_MATCH" \
        "$REFERENCE_CRES_BED" \
        "$chain" \
        "$mapped" \
        "$unmapped"


    # --------------------------------------------------------
    # QC: mapped + unmapped must equal reference CRE count
    # --------------------------------------------------------

    n_mapped="$(
        grep -vcE '^(#|$)' "$mapped" || true
    )"

    n_unmapped="$(
        grep -vcE '^(#|$)' "$unmapped" || true
    )"

    total_cre=$((n_mapped + n_unmapped))

    printf "%-20s %8d %10d %8d\n" \
        "$sp" \
        "$n_mapped" \
        "$n_unmapped" \
        "$total_cre"

    if [[ "$total_cre" -ne "$N_REF_CRES" ]]; then

        echo "WARNING: $sp total=$total_cre expected=$N_REF_CRES" >&2
        ((qc_fail += 1))

    fi


    # --------------------------------------------------------
    # QC: mapped CRE identifiers must be unique
    # --------------------------------------------------------

    duplicates="$(
        cut -f4 "$mapped" \
            | sort \
            | uniq -d \
            | wc -l
    )"

    if [[ "$duplicates" -ne 0 ]]; then

        echo "WARNING: $sp has duplicated mapped CRE IDs." >&2
        ((qc_fail += 1))

    fi

done < "$TARGET_SPECIES_FILE"


echo

if [[ "$qc_fail" -ne 0 ]]; then
    echo "ERROR: liftOver QC failed." >&2
    exit 1
fi


# ============================================================
# Final summary
# ============================================================

echo "============================================================"
echo "Reference CRE liftOver completed"
echo "============================================================"
echo "Reference CREs : $N_REF_CRES"
echo "Target species : $total"
echo "minMatch       : $LIFTOVER_MIN_MATCH"
echo "Output         : $LIFTED_CRES_DIR"
echo "Finished       : $(date)"
echo "============================================================"
