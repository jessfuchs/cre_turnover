#!/usr/bin/env bash
set -euo pipefail


# ============================================================
# 01 - Prepare genomes for pairwise whole-genome alignment
#
# Purpose:
#   Convert all species genome FASTA files to UCSC 2bit format,
#   generate chromosome-size files, and define the target
#   species used for pairwise alignments against D. melanogaster.
#
# Input:
#   - combined_manifest.tsv
#   - species-specific genome FASTA files
#   - dmel_reference_cres.bed
#
# Output:
#   - twobit/<species>.2bit
#   - chrom_sizes/<species>.sizes
#   - target_species.txt
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
# Input and tool checks
# ============================================================

[[ -s "$COMBINED_MANIFEST" ]] || {
    echo "ERROR: combined manifest missing or empty:" >&2
    echo "  $COMBINED_MANIFEST" >&2
    exit 1
}

[[ -s "$REFERENCE_CRES_BED" ]] || {
    echo "ERROR: reference CRE BED missing or empty:" >&2
    echo "  $REFERENCE_CRES_BED" >&2
    exit 1
}

command -v faToTwoBit >/dev/null 2>&1 || {
    echo "ERROR: faToTwoBit not found." >&2
    exit 1
}

command -v twoBitInfo >/dev/null 2>&1 || {
    echo "ERROR: twoBitInfo not found." >&2
    exit 1
}


# ============================================================
# Read combined species manifest
# ============================================================

MANIFEST_ROWS="$(
python3 - "$COMBINED_MANIFEST" <<'PY'
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
        "source",
    }

    missing = required - set(
        reader.fieldnames or []
    )

    if missing:
        raise SystemExit(
            "ERROR: combined manifest is missing required columns: "
            + ", ".join(sorted(missing))
        )

    for row in reader:

        slug = row["slug"].strip()
        source = row["source"].strip()

        if not slug:
            continue

        print(f"{slug}\t{source}")
PY
)"


# ============================================================
# Prepare output directories
# ============================================================

mkdir -p \
    "$TWOBIT_DIR" \
    "$CHROM_SIZES_DIR" \
    "$ALIGNMENT_DIR" \
    "$LIFTED_CRES_DIR" \
    "$LOG_DIR"

rm -f "$TARGET_SPECIES_FILE"


# ============================================================
# Header
# ============================================================

echo "============================================================"
echo "Prepare genomes for pairwise WGA"
echo "============================================================"
echo "Manifest  : $COMBINED_MANIFEST"
echo "2bit      : $TWOBIT_DIR"
echo "Sizes     : $CHROM_SIZES_DIR"
echo "Reference : $REFERENCE_SPECIES"
echo


# ============================================================
# Process all species
# ============================================================

while IFS=$'\t' read -r slug source; do

    [[ -n "$slug" ]] || continue

    case "$source" in

        generated)
            genome="$GENERATED_GENOME_ROOT/$slug/genome.fa"
            ;;

        external)
            genome="$EXTERNAL_GENOME_ROOT/$slug/genome.fa"
            ;;

        *)
            echo "ERROR: unknown source='$source' for '$slug'." >&2
            exit 1
            ;;
    esac


    # --------------------------------------------------------
    # Validate genome FASTA
    # --------------------------------------------------------

    if [[ ! -s "$genome" ]]; then
        echo "ERROR: genome FASTA missing for $slug:" >&2
        echo "  $genome" >&2
        exit 1
    fi


    twobit="$TWOBIT_DIR/${slug}.2bit"
    sizes="$CHROM_SIZES_DIR/${slug}.sizes"

    echo "[$slug]"
    echo "  source : $source"
    echo "  genome : $genome"


    # --------------------------------------------------------
    # FASTA -> 2bit
    # --------------------------------------------------------

    rm -f "$twobit" "$sizes"

    faToTwoBit \
        "$genome" \
        "$twobit"

    [[ -s "$twobit" ]] || {
        echo "ERROR: 2bit file was not created for $slug." >&2
        exit 1
    }


    # --------------------------------------------------------
    # 2bit -> chromosome sizes
    # --------------------------------------------------------

    twoBitInfo \
        "$twobit" \
        "$sizes"

    [[ -s "$sizes" ]] || {
        echo "ERROR: chromosome-size file was not created for $slug." >&2
        exit 1
    }


    # --------------------------------------------------------
    # Target-species list
    # --------------------------------------------------------

    if [[ "$slug" != "$REFERENCE_SPECIES" ]]; then
        printf '%s\n' "$slug" >> "$TARGET_SPECIES_FILE"
    fi

    echo "  2bit  : $twobit"
    echo "  sizes : $sizes"
    echo

done <<< "$MANIFEST_ROWS"


# ============================================================
# QC
# ============================================================

echo "============================================================"
echo "QC"
echo "============================================================"


# ------------------------------------------------------------
# Number of prepared species
# ------------------------------------------------------------

n_twobit="$(
    find "$TWOBIT_DIR" \
        -maxdepth 1 \
        -type f \
        -name '*.2bit' \
        | wc -l
)"

n_sizes="$(
    find "$CHROM_SIZES_DIR" \
        -maxdepth 1 \
        -type f \
        -name '*.sizes' \
        | wc -l
)"

n_targets="$(
    awk '
        /^[[:space:]]*#/ {next}
        NF == 0 {next}
        {n++}
        END {print n+0}
    ' "$TARGET_SPECIES_FILE"
)"

echo "2bit files : $n_twobit"
echo "Size files : $n_sizes"
echo "Targets    : $n_targets"


# ------------------------------------------------------------
# Check duplicate targets
# ------------------------------------------------------------

duplicates="$(
    sort "$TARGET_SPECIES_FILE" \
        | uniq -d
)"

if [[ -n "$duplicates" ]]; then

    echo "ERROR: duplicate target species detected:" >&2
    echo "$duplicates" >&2
    exit 1
fi


# ------------------------------------------------------------
# Check reference genome
# ------------------------------------------------------------

REF_2BIT="$TWOBIT_DIR/${REFERENCE_SPECIES}.2bit"
REF_SIZES="$CHROM_SIZES_DIR/${REFERENCE_SPECIES}.sizes"

[[ -s "$REF_2BIT" ]] || {
    echo "ERROR: D. melanogaster reference 2bit file missing." >&2
    exit 1
}

[[ -s "$REF_SIZES" ]] || {
    echo "ERROR: D. melanogaster chromosome-size file missing." >&2
    exit 1
}

echo
echo "D. melanogaster reference sequence IDs:"
head "$REF_SIZES"


# ------------------------------------------------------------
# Check reference CRE sequence IDs
# ------------------------------------------------------------

missing_ref_seqids="$(
    comm -23 \
        <(cut -f1 "$REFERENCE_CRES_BED" | sort -u) \
        <(cut -f1 "$REF_SIZES" | sort -u)
)"

if [[ -n "$missing_ref_seqids" ]]; then

    echo
    echo "ERROR: reference CRE BED contains sequence IDs that are" >&2
    echo "absent from the D. melanogaster reference genome:" >&2
    echo "$missing_ref_seqids" >&2
    exit 1
fi

echo
echo "Reference BED sequence IDs: OK"


# ------------------------------------------------------------
# Reference CRE count
# ------------------------------------------------------------

n_ref_cres="$(
    grep -vcE '^(#|$)' "$REFERENCE_CRES_BED"
)"

echo "Reference CREs: $n_ref_cres"


# ============================================================
# Final summary
# ============================================================

echo
echo "============================================================"
echo "WGA genome preparation completed"
echo "============================================================"
echo "Prepared genomes : $n_twobit"
echo "Reference        : $REFERENCE_SPECIES"
echo "Target species   : $n_targets"
echo "Target list      : $TARGET_SPECIES_FILE"
echo "============================================================"
