#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Filter external SCRMshaw BED files
#
# Input:
#   species_external.txt
#
# Expected format:
#   slug    scientific_name
#
# Example:
#   d_affinis    Drosophila affinis
#
# Raw BED:
#   external_data/scrmshaw_output/<slug>.bed
#
# Output:
#   filtered/<slug>.<training>.<method>.bed
#
# Existing filtered BEDs are skipped only if they pass QC.
# ============================================================


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

source "$ROOT/config/external_config.sh"

SPECIES_FILE="${1:-$ROOT/species_external.txt}"

EXPECTED_OFFSETS="${EXPECTED_OFFSETS:-25}"
MAX_RANK="${MAX_RANK:-5000}"

THEORETICAL_MAX=$((EXPECTED_OFFSETS * MAX_RANK))


# ------------------------------------------------------------
# Preconditions
# ------------------------------------------------------------

[[ -s "$SPECIES_FILE" ]] || {
    echo "ERROR: Species file is missing or empty:" >&2
    echo "  $SPECIES_FILE" >&2
    exit 1
}

[[ -d "$EXTERNAL_BED_DIR" ]] || {
    echo "ERROR: External BED directory is missing:" >&2
    echo "  $EXTERNAL_BED_DIR" >&2
    exit 1
}

mkdir -p "$FILTERED_BED_DIR"


# ============================================================
# QC function
#
# Returns 0 only if an existing filtered BED is valid.
# ============================================================

validate_filtered_bed() {

    bed="$1"

    [[ -s "$bed" ]] || return 1

    awk \
        -F'\t' \
        -v training="$TRAINING_SET" \
        -v method="$METHOD" \
        -v expected_offsets="$EXPECTED_OFFSETS" \
        -v max_rank_allowed="$MAX_RANK" \
        -v theoretical_max="$THEORETICAL_MAX" '
        BEGIN {
            rows = 0
            rank1 = 0
            min_rank = -1
            max_rank = -1
        }

        /^#/ {
            next
        }

        NF != 17 {
            exit 10
        }

        $15 != training {
            exit 11
        }

        $16 != method {
            exit 12
        }

        $17 !~ /^[0-9]+$/ {
            exit 13
        }

        {
            rank = $17 + 0

            if (rank < 1) {
                exit 14
            }

            rows++

            if (rank == 1) {
                rank1++
            }

            if (min_rank == -1 || rank < min_rank) {
                min_rank = rank
            }

            if (rank > max_rank) {
                max_rank = rank
            }
        }

        END {
            if (rows == 0) {
                exit 20
            }

            if (min_rank != 1) {
                exit 21
            }

            if (rank1 != expected_offsets) {
                exit 22
            }

            if (max_rank > max_rank_allowed) {
                exit 23
            }

            if (rows > theoretical_max) {
                exit 24
            }
        }
    ' "$bed"
}


# ============================================================
# Main
# ============================================================

n_total=0
n_filtered=0
n_skipped=0
n_rebuilt=0


echo "============================================================"
echo "Filter external SCRMshaw BED files"
echo "============================================================"
echo "Species file      : $SPECIES_FILE"
echo "Raw BED dir       : $EXTERNAL_BED_DIR"
echo "Filtered BED dir  : $FILTERED_BED_DIR"
echo "Training set      : $TRAINING_SET"
echo "Method            : $METHOD"
echo "Expected offsets  : $EXPECTED_OFFSETS"
echo "Maximum rank      : $MAX_RANK"
echo "Theoretical max   : $THEORETICAL_MAX"
echo "============================================================"
echo


while IFS= read -r line || [[ -n "$line" ]]; do

    # Skip blank lines
    [[ -z "$(printf '%s' "$line" | tr -d '[:space:]')" ]] && continue

    # Skip comments
    case "$line" in
        \#*) continue ;;
    esac

    # Expected:
    # d_affinis   Drosophila affinis

    read -r slug genus epithet extra <<EOF
$line
EOF

    if [[ -z "${slug:-}" ||
          -z "${genus:-}" ||
          -z "${epithet:-}" ]]; then

        echo "ERROR: Invalid species line:" >&2
        echo "  $line" >&2
        echo "Expected:" >&2
        echo "  slug Genus species" >&2
        exit 1
    fi

    if [[ -n "${extra:-}" ]]; then
        echo "ERROR: Too many fields in species line:" >&2
        echo "  $line" >&2
        exit 1
    fi


    scientific_name="$genus $epithet"

    input_bed="$EXTERNAL_BED_DIR/${slug}.bed"

    output_bed="$FILTERED_BED_DIR/${slug}.${TRAINING_SET}.${METHOD}.bed"

    tmp_bed="${output_bed}.tmp.$$"

    n_total=$((n_total + 1))


    echo "------------------------------------------------------------"
    echo "[$slug] $scientific_name"
    echo "Input : $input_bed"
    echo "Output: $output_bed"


    # ========================================================
    # Robust skip
    # ========================================================

    if [[ -s "$output_bed" ]]; then

        echo "  existing filtered BED found"
        echo "  checking QC..."

        if validate_filtered_bed "$output_bed"; then

            rows="$(awk -F'\t' '!/^#/ && NF {n++} END {print n+0}' "$output_bed")"

            min_rank="$(awk -F'\t' '
                !/^#/ && NF {
                    r=$17+0
                    if (!seen || r<min) min=r
                    seen=1
                }
                END {print min+0}
            ' "$output_bed")"

            max_rank="$(awk -F'\t' '
                !/^#/ && NF {
                    r=$17+0
                    if (!seen || r>max) max=r
                    seen=1
                }
                END {print max+0}
            ' "$output_bed")"

            rank1="$(awk -F'\t' '
                !/^#/ && NF && $17==1 {n++}
                END {print n+0}
            ' "$output_bed")"

            echo "  OK: existing file is valid"
            echo "      rows     : $rows"
            echo "      rank1    : $rank1"
            echo "      min rank : $min_rank"
            echo "      max rank : $max_rank"
            echo "  -> skipping filtering"

            n_skipped=$((n_skipped + 1))

            echo
            continue

        else

            echo "  WARNING: existing file does not pass QC"
            echo "  -> will be rebuilt"

            n_rebuilt=$((n_rebuilt + 1))
        fi

    fi


    # ========================================================
    # Check raw BED
    # ========================================================

    if [[ ! -s "$input_bed" ]]; then
        echo "ERROR: Raw BED is missing or empty:" >&2
        echo "  $input_bed" >&2
        exit 1
    fi


    # ========================================================
    # Filter
    #
    # One pass through the raw BED.
    #
    # Expected columns:
    #   15 = training set
    #   16 = method
    #   17 = rank
    # ========================================================

    echo "  filtering ${TRAINING_SET} / ${METHOD} ..."

    rm -f "$tmp_bed"

    awk \
        -F'\t' \
        -v OFS='\t' \
        -v training="$TRAINING_SET" \
        -v method="$METHOD" '
        /^#/ {
            next
        }

        NF != 17 {
            printf(
                "ERROR: Line %d has %d fields instead of 17\n",
                NR,
                NF
            ) > "/dev/stderr"

            exit 10
        }

        $15 == training && $16 == method {
            print
        }
    ' "$input_bed" > "$tmp_bed"


    # ========================================================
    # Validate newly filtered BED
    # ========================================================

    if ! validate_filtered_bed "$tmp_bed"; then

        rm -f "$tmp_bed"

        echo "ERROR: Newly filtered BED does not pass QC:" >&2
        echo "  $input_bed" >&2
        echo "  Training: $TRAINING_SET" >&2
        echo "  Method  : $METHOD" >&2

        exit 1
    fi


    mv "$tmp_bed" "$output_bed"


    # ========================================================
    # Summary statistics
    # ========================================================

    read -r rows rank1 min_rank max_rank <<EOF
$(awk -F'\t' '
    !/^#/ && NF {
        rows++

        r=$17+0

        if (r==1)
            rank1++

        if (!seen || r<min)
            min=r

        if (!seen || r>max)
            max=r

        seen=1
    }

    END {
        print rows+0, rank1+0, min+0, max+0
    }
' "$output_bed")
EOF


    echo "  OK:"
    echo "      rows     : $rows"
    echo "      rank1    : $rank1"
    echo "      min rank : $min_rank"
    echo "      max rank : $max_rank"
    echo "      output   : $output_bed"

    n_filtered=$((n_filtered + 1))

    echo

done < "$SPECIES_FILE"


echo "============================================================"
echo "External BED filtering completed"
echo "============================================================"
echo "Species total       : $n_total"
echo "Newly filtered      : $n_filtered"
echo "Valid files skipped : $n_skipped"
echo "Invalid rebuilt     : $n_rebuilt"
echo "============================================================"
