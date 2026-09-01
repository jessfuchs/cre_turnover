#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

source "$ROOT/config/external_config.sh"

SPECIES_LIST="${1:-$ROOT/species_external.txt}"

[[ -s "$SPECIES_LIST" ]] || {
    echo "ERROR: Species file is missing or empty: $SPECIES_LIST" >&2
    exit 1
}

[[ -d "$SOURCE_GFF_ROOT" ]] || {
    echo "ERROR: selected directory is missing: $SOURCE_GFF_ROOT" >&2
    exit 1
}

mkdir -p "$EXTERNAL_GFF_DIR"

n_total=0
n_existing=0
n_copied=0


validate_gff() {
    gff="$1"

    [[ -s "$gff" ]] || return 1

    awk -F'\t' '
        /^#/ {next}
        NF >= 9 {
            found=1
            exit
        }
        END {
            if (!found) exit 1
        }
    ' "$gff"
}


echo "============================================================"
echo "Prepare external GFF files"
echo "============================================================"
echo "Species file : $SPECIES_LIST"
echo "Source GFFs  : $SOURCE_GFF_ROOT"
echo "Target GFFs  : $EXTERNAL_GFF_DIR"
echo


while IFS= read -r line || [[ -n "$line" ]]; do

    [[ -z "$(printf '%s' "$line" | tr -d '[:space:]')" ]] && continue

    case "$line" in
        \#*) continue ;;
    esac

    read -r slug genus epithet extra <<EOF
$line
EOF

    if [[ -z "${slug:-}" ||
          -z "${genus:-}" ||
          -z "${epithet:-}" ]]; then

        echo "ERROR: Invalid species line:" >&2
        echo "  $line" >&2
        exit 1
    fi

    if [[ -n "${extra:-}" ]]; then
        echo "ERROR: Too many fields in species line:" >&2
        echo "  $line" >&2
        exit 1
    fi

    scientific_name="$genus $epithet"

    source_gff="$SOURCE_GFF_ROOT/${slug}/annotation.gff3"
    target_gff="$EXTERNAL_GFF_DIR/${slug}.gff3"

    n_total=$((n_total + 1))

    echo "[$slug] $scientific_name"

    # Already prepared?
    if validate_gff "$target_gff"; then
        echo "  OK: GFF already prepared:"
        echo "      $target_gff"

        n_existing=$((n_existing + 1))
        echo
        continue
    fi

    # Check source GFF
    if ! validate_gff "$source_gff"; then
        echo "ERROR: No valid source GFF found:" >&2
        echo "  Species: $scientific_name ($slug)" >&2
        echo "  Expected: $source_gff" >&2
        exit 1
    fi

    echo "  Source:"
    echo "      $source_gff"

    # Remove existing/stale target file or symlink
    if [[ -e "$target_gff" || -L "$target_gff" ]]; then
        rm -f "$target_gff"
    fi

    cp "$source_gff" "$target_gff"

    if ! validate_gff "$target_gff"; then
        rm -f "$target_gff"

        echo "ERROR: Copied GFF is invalid:" >&2
        echo "  $target_gff" >&2
        exit 1
    fi

    echo "  copied as:"
    echo "      $target_gff"

    n_copied=$((n_copied + 1))
    echo

done < "$SPECIES_LIST"


echo "============================================================"
echo "External GFF preparation completed"
echo "============================================================"
echo "Species total    : $n_total"
echo "Already prepared : $n_existing"
echo "Newly copied     : $n_copied"
echo "============================================================"
