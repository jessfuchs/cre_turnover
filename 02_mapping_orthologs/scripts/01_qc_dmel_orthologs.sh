#!/usr/bin/env bash
set -euo pipefail


# ============================================================
# 01 - QC of D. melanogaster ortholog annotations
#
# Purpose:
#   Validate that all species GFF3 files contain usable
#   dmel_orthologs annotations before ortholog mapping starts.
#
# Source of truth:
#   external_scrmshaw/combined_manifest.tsv
#
# Per species, this script checks:
#   - GFF exists and is non-empty
#   - mRNA features are present
#   - dmel_orthologs= attribute occurs
#   - at least one non-empty/non-NA Dmel ortholog exists
#   - fraction of mRNAs with a Dmel ortholog exceeds a
#     conservative technical minimum
#
# Output:
#   mapping_orthologs/ortholog_results/
#       dmel_ortholog_annotation_qc.tsv
#
# Exit status:
#   0 = QC passed
#   1 = QC failed; downstream mapping should not continue
# ============================================================


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

PROJECT="$HOME/cre_turnover/project"

GENERATED_ROOT="$PROJECT/scrmshaw_pipeline"
EXTERNAL_ROOT="$PROJECT/external_scrmshaw"
MAPPING_ROOT="$PROJECT/mapping_orthologs"

MANIFEST="$EXTERNAL_ROOT/combined_manifest.tsv"

ORTHO_RESULTS="$MAPPING_ROOT/ortholog_results"

OUTPUT="$ORTHO_RESULTS/dmel_ortholog_annotation_qc.tsv"


# ------------------------------------------------------------
# Expected final project state
# ------------------------------------------------------------

EXPECTED_SPECIES=40

# Conservative technical sanity threshold only.
# This is NOT intended as a biological filtering criterion.
MIN_FRACTION=0.10


# ============================================================
# Input checks
# ============================================================

[[ -s "$MANIFEST" ]] || {
    echo "ERROR: combined manifest missing or empty:" >&2
    echo "  $MANIFEST" >&2
    exit 1
}

mkdir -p "$ORTHO_RESULTS"


# ============================================================
# Output header
# ============================================================

printf '%s\n' \
"slug	species	source	total_mRNA	with_dmel_attr	with_dmel_ortholog	without_dmel_ortholog	fraction_with_dmel_ortholog	status" \
> "$OUTPUT"


# ============================================================
# Counters
# ============================================================

n_species=0
n_missing_gff=0
n_warn=0
n_failed=0


# ============================================================
# Header
# ============================================================

echo "Manifest      : $MANIFEST"
echo "Output        : $OUTPUT"
echo "Expected species : $EXPECTED_SPECIES"
echo "Min. fraction : $MIN_FRACTION"
echo


# ============================================================
# Read manifest robustly
#
# Python parses the TSV.
# A pipe character is used only for transfer into Bash to avoid
# accidental problems with Bash IFS/tab handling.
# ============================================================

while IFS='|' read -r slug species source; do

    [[ -n "$slug" ]] || continue

    n_species=$((n_species + 1))


    # --------------------------------------------------------
    # Locate species GFF
    # --------------------------------------------------------

    case "$source" in

        generated)
            gff="$GENERATED_ROOT/data/selected/$slug/annotation.gff3"
            ;;

        external)
            gff="$EXTERNAL_ROOT/external_data/selected/$slug/annotation.gff3"
            ;;

        *)
            echo "ERROR: unknown source='$source' for '$slug'" >&2

            printf "%s\t%s\t%s\t0\t0\t0\t0\t0.0000\tFAIL_UNKNOWN_SOURCE\n" \
                "$slug" \
                "$species" \
                "$source" \
                >> "$OUTPUT"

            n_failed=$((n_failed + 1))
            continue
            ;;
    esac


    # --------------------------------------------------------
    # Check GFF existence
    # --------------------------------------------------------

    if [[ ! -s "$gff" ]]; then

        echo "ERROR: GFF missing or empty for $slug:" >&2
        echo "  $gff" >&2

        printf "%s\t%s\t%s\t0\t0\t0\t0\t0.0000\tFAIL_MISSING_GFF\n" \
            "$slug" \
            "$species" \
            "$source" \
            >> "$OUTPUT"

        n_missing_gff=$((n_missing_gff + 1))
        n_failed=$((n_failed + 1))

        continue
    fi


    # --------------------------------------------------------
    # Count mRNA and Dmel ortholog annotations
    # --------------------------------------------------------

    result="$(
    awk -F'\t' '
        BEGIN {
            total  = 0
            attr   = 0
            mapped = 0
        }

        /^#/ {
            next
        }

        NF >= 9 && $3 == "mRNA" {

            total++

            if ($9 ~ /(^|;)dmel_orthologs=/) {

                attr++

                if ($9 !~ /(^|;)dmel_orthologs=NA([;]|$)/ && $9 !~ /(^|;)dmel_orthologs=([;]|$)/) {
                    mapped++
                }
            }
        }

        END {

            unmapped = total - mapped

            if (total > 0) {
                fraction = mapped / total
            } else {
                fraction = 0
            }

            printf "%d\t%d\t%d\t%d\t%.6f\n",
                total,
                attr,
                mapped,
                unmapped,
                fraction
        }
    ' "$gff"
)"

    IFS=$'\t' read -r \
        total \
        attr \
        mapped \
        unmapped \
        fraction \
        <<< "$result"


    # --------------------------------------------------------
    # Determine QC status
    # --------------------------------------------------------

    status="PASS"

    # D. melanogaster is the reference species.
    # It does not need dmel_orthologs annotations to map to itself.
    if [[ "$slug" == "d_melanogaster" ]]; then

        if [[ "$total" -eq 0 ]]; then
            status="FAIL_NO_MRNA"
            n_failed=$((n_failed + 1))
        else
            status="PASS_REFERENCE_SPECIES"
        fi

    elif [[ "$total" -eq 0 ]]; then

        status="FAIL_NO_MRNA"
        n_failed=$((n_failed + 1))

    elif [[ "$attr" -eq 0 ]]; then

        status="FAIL_NO_DMEL_ATTRIBUTE"
        n_failed=$((n_failed + 1))

    elif [[ "$mapped" -eq 0 ]]; then

        status="FAIL_NO_DMEL_ORTHOLOG"
        n_failed=$((n_failed + 1))

    elif ! awk \
        -v f="$fraction" \
        -v min="$MIN_FRACTION" \
        'BEGIN { exit !(f >= min) }'
    then

        status="FAIL_LOW_ORTHOLOG_FRACTION"
        n_failed=$((n_failed + 1))

    elif [[ "$attr" -lt "$total" ]]; then

        status="PASS"
        n_warn=$((n_warn + 1))

    fi


    # --------------------------------------------------------
    # Write QC row
    # --------------------------------------------------------

    printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%s\n" \
        "$slug" \
        "$species" \
        "$source" \
        "$total" \
        "$attr" \
        "$mapped" \
        "$unmapped" \
        "$fraction" \
        "$status" \
        >> "$OUTPUT"


done < <(
    python3 - "$MANIFEST" <<'PY'
import csv
import sys

manifest = sys.argv[1]

with open(
    manifest,
    encoding="utf-8-sig",
    newline=""
) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )

    required = {
        "slug",
        "species",
        "source",
    }

    missing = required - set(reader.fieldnames or [])

    if missing:
        raise SystemExit(
            "ERROR: missing manifest columns: "
            + ", ".join(sorted(missing))
        )

    for row in reader:

        slug = row["slug"].strip()
        species = row["species"].strip()
        source = row["source"].strip()

        if not slug:
            continue

        # Internal separator only.
        # Species names in this project do not contain "|".
        print(f"{slug}|{species}|{source}")
PY
)


# ============================================================
# Final project-level checks
# ============================================================

if [[ "$n_species" -ne "$EXPECTED_SPECIES" ]]; then

    echo "ERROR: expected $EXPECTED_SPECIES species, found $n_species." >&2
    n_failed=$((n_failed + 1))

fi


# ============================================================
# Summary
# ============================================================

echo
echo "============================================================"
echo "QC summary"
echo "============================================================"
echo "Species in manifest : $n_species"
echo "Missing GFFs        : $n_missing_gff"
echo "Warnings            : $n_warn"
echo "Failures            : $n_failed"
echo


# ------------------------------------------------------------
# Show result table
# ------------------------------------------------------------

if command -v column >/dev/null 2>&1; then
    column -t -s $'\t' "$OUTPUT"
else
    cat "$OUTPUT"
fi


# ============================================================
# Gate downstream pipeline
# ============================================================

echo

if [[ "$n_failed" -ne 0 ]]; then

    echo "============================================================"
    echo "QC FAILED"
    echo "============================================================"
    echo "Ortholog mapping will NOT continue."
    echo
    echo "Inspect:"
    echo "  $OUTPUT"
    echo "============================================================"

    exit 1
fi


echo "============================================================"
echo "QC PASSED"
echo "============================================================"
echo "Dmel ortholog annotations are technically sufficient"
echo "for downstream ortholog mapping."
echo
echo "Warnings are non-fatal and are retained in the QC report."
echo
echo "Report:"
echo "  $OUTPUT"
echo "============================================================"
