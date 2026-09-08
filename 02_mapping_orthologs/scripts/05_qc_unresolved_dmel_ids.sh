#!/usr/bin/env bash
set -euo pipefail


# ============================================================
# 05 - QC of unresolved D. melanogaster identifiers
#
# Purpose:
#   Investigate Dmel identifiers that could not be converted to
#   FlyBase FBgn IDs and test whether alternative identifier
#   normalization resolves them against the Dmel GFF3.
#
# Input:
#   - unresolved_dmel_identifiers.tsv
#   - D. melanogaster GFF3 annotation
#
# Output:
#   - unresolved_dmel_qc_summary.tsv
#   - detailed QC files in results/unresolved_dmel_qc/
#
# Configuration:
#   config/ortholog_config.sh
# ============================================================


# ============================================================
# Configuration
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

source "$PIPELINE_ROOT/config/ortholog_config.sh"

UNRESOLVED="$UNRESOLVED_DMEL"
SUMMARY="$UNRESOLVED_QC_SUMMARY"
TMPDIR="$UNRESOLVED_QC_DIR"

mkdir -p "$TMPDIR"


# ============================================================
# Input checks
# ============================================================

[[ -s "$UNRESOLVED" ]] || {
    echo "ERROR: unresolved identifier file missing or empty:" >&2
    echo "  $UNRESOLVED" >&2
    exit 1
}

[[ -s "$DMEL_GFF" ]] || {
    echo "ERROR: D. melanogaster GFF3 missing or empty:" >&2
    echo "  $DMEL_GFF" >&2
    exit 1
}


# ============================================================
# Extract unique unresolved identifiers
# ============================================================

python3 - "$UNRESOLVED" "$TMPDIR/unresolved_ids.txt" <<'PY'
import csv
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

ids = set()

with open(infile, encoding="utf-8-sig", newline="") as handle:
    rows = list(csv.reader(handle, delimiter="\t"))

if not rows:
    raise SystemExit("ERROR: unresolved identifier table is empty")

header = rows[0]

candidate_names = {
    "identifier",
    "dmel_identifier",
    "unresolved_identifier",
    "dmel_id",
    "id",
}

id_col = None

for i, name in enumerate(header):
    if name.strip().lower() in candidate_names:
        id_col = i
        break

if id_col is None:
    id_col = 0

for row in rows[1:]:

    if len(row) <= id_col:
        continue

    value = row[id_col].strip()

    if not value:
        continue

    for identifier in value.split("|"):

        identifier = identifier.strip()

        if identifier and identifier not in {"NA", "."}:
            ids.add(identifier)

with open(outfile, "w", encoding="utf-8", newline="\n") as out:

    for identifier in sorted(ids):
        out.write(identifier + "\n")
PY


# ============================================================
# Classify unresolved identifiers
# ============================================================

grep '^Dmelanogast_M' \
    "$TMPDIR/unresolved_ids.txt" \
    > "$TMPDIR/dmelanogast_M_ids.txt" || true

grep -v '^Dmelanogast_M' \
    "$TMPDIR/unresolved_ids.txt" \
    > "$TMPDIR/other_unresolved_ids.txt" || true


# ============================================================
# Normalize Dmelanogast_M identifiers
#
# Example:
#   Dmelanogast_M00000001234
#   -> Dmelanogast_00000001234
# ============================================================

sed 's/^Dmelanogast_M/Dmelanogast_/' \
    "$TMPDIR/dmelanogast_M_ids.txt" \
    > "$TMPDIR/dmelanogast_normalized_ids.txt"


# ============================================================
# Compare unresolved and normalized IDs against the Dmel GFF3
# ============================================================

python3 - \
    "$DMEL_GFF" \
    "$TMPDIR/unresolved_ids.txt" \
    "$TMPDIR/dmelanogast_normalized_ids.txt" \
    "$TMPDIR/unresolved_literal_in_gff.txt" \
    "$TMPDIR/normalized_in_gff.txt" \
    "$TMPDIR/normalized_fbgn.txt" <<'PY'

import re
import sys

(
    gff_file,
    unresolved_file,
    normalized_file,
    literal_out,
    normalized_out,
    fbgn_out,
) = sys.argv[1:]


def read_ids(path):
    with open(path, encoding="utf-8") as handle:
        return {
            line.strip()
            for line in handle
            if line.strip()
        }


unresolved = read_ids(unresolved_file)
normalized = read_ids(normalized_file)

literal_found = set()
normalized_found = set()
normalized_fbgn = set()


with open(
    gff_file,
    encoding="utf-8",
    errors="replace",
) as handle:

    for line in handle:

        if line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")

        if len(parts) < 9:
            continue

        attributes = parts[8]

        for identifier in unresolved:
            if identifier in attributes:
                literal_found.add(identifier)

        matched_normalized = [
            identifier
            for identifier in normalized
            if identifier in attributes
        ]

        if matched_normalized:

            normalized_found.update(
                matched_normalized
            )

            fbgns = re.findall(
                r"FBgn\d+",
                attributes
            )

            normalized_fbgn.update(
                fbgns
            )


for path, values in [
    (literal_out, literal_found),
    (normalized_out, normalized_found),
    (fbgn_out, normalized_fbgn),
]:

    with open(
        path,
        "w",
        encoding="utf-8",
        newline="\n",
    ) as out:

        for value in sorted(values):
            out.write(value + "\n")
PY


# ============================================================
# Count QC categories
# ============================================================

count_lines() {

    local file="$1"

    if [[ -s "$file" ]]; then
        wc -l < "$file"
    else
        echo 0
    fi
}


n_total=$(count_lines "$TMPDIR/unresolved_ids.txt")
n_m=$(count_lines "$TMPDIR/dmelanogast_M_ids.txt")
n_other=$(count_lines "$TMPDIR/other_unresolved_ids.txt")
n_literal=$(count_lines "$TMPDIR/unresolved_literal_in_gff.txt")
n_normalized=$(count_lines "$TMPDIR/dmelanogast_normalized_ids.txt")
n_normalized_found=$(count_lines "$TMPDIR/normalized_in_gff.txt")
n_fbgn=$(count_lines "$TMPDIR/normalized_fbgn.txt")


# ============================================================
# Write QC summary
# ============================================================

printf "metric\tcount\n" > "$SUMMARY"

printf "unique_unresolved_identifiers\t%d\n" \
    "$n_total" >> "$SUMMARY"

printf "Dmelanogast_M_identifiers\t%d\n" \
    "$n_m" >> "$SUMMARY"

printf "other_unresolved_identifiers\t%d\n" \
    "$n_other" >> "$SUMMARY"

printf "unresolved_ids_found_literally_in_gff\t%d\n" \
    "$n_literal" >> "$SUMMARY"

printf "normalized_Dmelanogast_ids\t%d\n" \
    "$n_normalized" >> "$SUMMARY"

printf "normalized_Dmelanogast_ids_found_in_gff\t%d\n" \
    "$n_normalized_found" >> "$SUMMARY"

printf "unique_FBgn_for_normalized_models\t%d\n" \
    "$n_fbgn" >> "$SUMMARY"


# ============================================================
# Report
# ============================================================

echo
echo "============================================================"
echo "Unresolved Dmel identifier QC completed"
echo "============================================================"
echo

if command -v column >/dev/null 2>&1; then
    column -t -s $'\t' "$SUMMARY"
else
    cat "$SUMMARY"
fi

echo
echo "Detailed QC files:"
echo "  $TMPDIR"
echo
echo "Summary:"
echo "  $SUMMARY"
echo
