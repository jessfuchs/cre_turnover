#!/usr/bin/env bash
set -euo pipefail


# ============================================================
# Paths
# ============================================================

PROJECT="$HOME/cre_turnover/project"

WGA="$PROJECT/pairwise_wga"
GENERATED_ROOT="$PROJECT/scrmshaw_pipeline"
EXTERNAL_ROOT="$PROJECT/external_scrmshaw"

COMBINED_MANIFEST="$EXTERNAL_ROOT/combined_manifest.tsv"

TWOBIT_DIR="$WGA/twobit"
SIZES_DIR="$WGA/chrom_sizes"

TARGETS="$WGA/target_species.txt"

REF="d_melanogaster"

export PATH="$WGA/bin:$PATH"


# ============================================================
# Input / tool checks
# ============================================================

[[ -s "$COMBINED_MANIFEST" ]] || {
    echo "FEHLER: Combined manifest fehlt:"
    echo "  $COMBINED_MANIFEST"
    exit 1
}

command -v faToTwoBit >/dev/null 2>&1 || {
    echo "FEHLER: faToTwoBit nicht gefunden."
    exit 1
}

command -v twoBitInfo >/dev/null 2>&1 || {
    echo "FEHLER: twoBitInfo nicht gefunden."
    exit 1
}


# ============================================================
# Read combined manifest
# ============================================================

MANIFEST_ROWS="$(
python3 - "$COMBINED_MANIFEST" <<'PY'
import csv
import sys

manifest = sys.argv[1]

with open(manifest, encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    required = {"slug", "source"}

    missing = required - set(reader.fieldnames or [])

    if missing:
        raise SystemExit(
            "FEHLER: Fehlende Manifest-Spalten: "
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

mkdir -p "$TWOBIT_DIR"
mkdir -p "$SIZES_DIR"

rm -f "$TARGETS"


echo "============================================================"
echo "Prepare genomes for pairwise WGA"
echo "============================================================"
echo "Manifest : $COMBINED_MANIFEST"
echo "2bit     : $TWOBIT_DIR"
echo "sizes    : $SIZES_DIR"
echo "reference: $REF"
echo


# ============================================================
# Process all species
# ============================================================

while IFS=$'\t' read -r slug source; do

    [[ -n "$slug" ]] || continue

    case "$source" in

        generated)
            genome="$GENERATED_ROOT/data/selected/$slug/genome.fa"
            ;;

        external)
            genome="$EXTERNAL_ROOT/external_data/selected/$slug/genome.fa"
            ;;

        *)
            echo "FEHLER: Unbekannte source='$source' für '$slug'" >&2
            exit 1
            ;;
    esac


    # --------------------------------------------------------
    # Validate FASTA
    # --------------------------------------------------------

    if [[ ! -s "$genome" ]]; then
        echo "FEHLER: Genome fehlt für $slug:" >&2
        echo "  $genome" >&2
        exit 1
    fi


    twobit="$TWOBIT_DIR/${slug}.2bit"
    sizes="$SIZES_DIR/${slug}.sizes"

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
        echo "FEHLER: 2bit wurde nicht erzeugt für $slug" >&2
        exit 1
    }


    # --------------------------------------------------------
    # 2bit -> chromosome sizes
    # --------------------------------------------------------

    twoBitInfo \
        "$twobit" \
        "$sizes"

    [[ -s "$sizes" ]] || {
        echo "FEHLER: sizes-Datei wurde nicht erzeugt für $slug" >&2
        exit 1
    }


    # --------------------------------------------------------
    # Target list
    # --------------------------------------------------------

    if [[ "$slug" != "$REF" ]]; then
        printf '%s\n' "$slug" >> "$TARGETS"
    fi

    echo "  2bit   : $twobit"
    echo "  sizes  : $sizes"
    echo

done <<< "$MANIFEST_ROWS"


# ============================================================
# QC
# ============================================================

echo "============================================================"
echo "QC"
echo "============================================================"


# ------------------------------------------------------------
# Number of species
# ------------------------------------------------------------

n_twobit="$(
    find "$TWOBIT_DIR" \
        -maxdepth 1 \
        -type f \
        -name '*.2bit' \
        | wc -l
)"

n_sizes="$(
    find "$SIZES_DIR" \
        -maxdepth 1 \
        -type f \
        -name '*.sizes' \
        | wc -l
)"

n_targets="$(
    grep -vcE '^(#|$)' "$TARGETS"
)"


echo "2bit files : $n_twobit"
echo "size files : $n_sizes"
echo "targets    : $n_targets"

# ------------------------------------------------------------
# Duplicate targets
# ------------------------------------------------------------

duplicates="$(
    sort "$TARGETS" \
        | uniq -d
)"

if [[ -n "$duplicates" ]]; then
    echo "FEHLER: Doppelte Target-Arten:"
    echo "$duplicates"
    exit 1
fi


# ------------------------------------------------------------
# Check reference
# ------------------------------------------------------------

REF_2BIT="$TWOBIT_DIR/${REF}.2bit"
REF_SIZES="$SIZES_DIR/${REF}.sizes"

[[ -s "$REF_2BIT" ]] || {
    echo "FEHLER: Dmel reference 2bit fehlt."
    exit 1
}

[[ -s "$REF_SIZES" ]] || {
    echo "FEHLER: Dmel reference sizes fehlen."
    exit 1
}


echo
echo "Dmel reference sequence IDs:"
head "$REF_SIZES"


# ------------------------------------------------------------
# Check reference BED sequence IDs against Dmel 2bit
# ------------------------------------------------------------

REF_BED="$PROJECT/mapping_orthologs/reference_cres/dmel_reference_cres.bed"

[[ -s "$REF_BED" ]] || {
    echo "FEHLER: Reference CRE BED fehlt:"
    echo "  $REF_BED"
    exit 1
}

missing_ref_seqids="$(
    comm -23 \
        <(cut -f1 "$REF_BED" | sort -u) \
        <(cut -f1 "$REF_SIZES" | sort -u)
)"

if [[ -n "$missing_ref_seqids" ]]; then

    echo
    echo "FEHLER: CRE-BED enthält SeqIDs, die im Dmel-Genome fehlen:"
    echo "$missing_ref_seqids"
    exit 1
fi


echo
echo "Reference BED sequence IDs: OK"


# ------------------------------------------------------------
# Check reference CRE count
# ------------------------------------------------------------

n_ref_cres="$(
    grep -vcE '^(#|$)' "$REF_BED"
)"

echo "Reference CREs: $n_ref_cres"


# --------------------------------------------------
# Prepare downstream output directories
# --------------------------------------------------

mkdir -p "$ROOT/alignments_dmel"
mkdir -p "$ROOT/lifted_cres_dmel"
mkdir -p "$ROOT/logs"

echo
echo "============================================================"
echo "WGA genome preparation completed"
echo "============================================================"
echo "All genomes prepared"
echo "Reference: $REF"
echo "Targets: $TARGETS"
echo "============================================================"
