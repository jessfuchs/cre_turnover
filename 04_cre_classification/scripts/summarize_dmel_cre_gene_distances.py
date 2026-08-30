#!/usr/bin/env python3

# ============================================================
# Summarize D. melanogaster CRE-to-gene distances
#
# Purpose:
#   Summarize distances between D. melanogaster SCRMshaw CREs
#   and their associated target genes to derive an empirical
#   local regulatory-neighborhood threshold.
#
# Input:
#   - SO_all_species_fbgn.tsv
#
# Output:
#   - dmel_cre_gene_distance_summary.tsv
# ============================================================

import argparse
import csv
import math
from pathlib import Path


# ============================================================
# Arguments
# ============================================================

parser = argparse.ArgumentParser(
    description=(
        "Summarize D. melanogaster CRE-to-associated-gene "
        "distances."
    )
)

parser.add_argument(
    "--input",
    required=True,
)

parser.add_argument(
    "--out",
    required=True,
)

parser.add_argument(
    "--reference-species",
    default="d_melanogaster",
)

args = parser.parse_args()


# ============================================================
# Paths
# ============================================================

INFILE = Path(args.input)
OUTFILE = Path(args.out)

if not INFILE.is_file():
    raise SystemExit(
        f"ERROR: required input file does not exist: {INFILE}"
    )

OUTFILE.parent.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Percentiles to report
# ============================================================

PERCENTILES = [
    25,
    50,
    75,
    90,
    95,
    97.5,
    99,
    100,
]


# ============================================================
# Helper function
# ============================================================

def percentile(values, p):
    """
    Calculate a percentile using linear interpolation.
    """

    if not values:
        return float("nan")

    k = (len(values) - 1) * p / 100

    lo = math.floor(k)
    hi = math.ceil(k)

    if lo == hi:
        return values[lo]

    return (
        values[lo] * (hi - k)
        + values[hi] * (k - lo)
    )


# ============================================================
# Collect D. melanogaster CRE-to-gene distances
# ============================================================

distances = []

with INFILE.open(
    encoding="utf-8-sig",
    newline=""
) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t",
    )

    required_columns = {
        "species_key",
        "dmel_fbgn_flanking_gene",
        "distance_flanking_gene",
        "dmel_fbgn_next_flanking_gene",
        "distance_next_gene",
    }

    missing_columns = (
        required_columns
        - set(reader.fieldnames or [])
    )

    if missing_columns:
        raise SystemExit(
            "ERROR: input table is missing required columns:\n"
            + "\n".join(
                sorted(missing_columns)
            )
        )

    for row in reader:

        if row["species_key"] != args.reference_species:
            continue

        candidates = [
            (
                row.get(
                    "dmel_fbgn_flanking_gene",
                    "",
                ),
                row.get(
                    "distance_flanking_gene",
                    "",
                ),
            ),
            (
                row.get(
                    "dmel_fbgn_next_flanking_gene",
                    "",
                ),
                row.get(
                    "distance_next_gene",
                    "",
                ),
            ),
        ]

        for fbgn_field, distance_field in candidates:

            fbgns = [
                value.strip()
                for value in fbgn_field.split("|")
                if value.strip().startswith("FBgn")
            ]

            if not fbgns:
                continue

            try:
                distance = abs(
                    float(distance_field)
                )

            except (TypeError, ValueError):
                continue

            # One distance is retained for each valid FBgn
            # association represented in the input table.
            distances.extend(
                distance
                for _ in fbgns
            )


# ============================================================
# QC
# ============================================================

if not distances:
    raise SystemExit(
        "ERROR: no valid CRE-to-gene associations found for "
        f"{args.reference_species}."
    )

distances.sort()


# ============================================================
# Write percentile summary
# ============================================================

with OUTFILE.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow([
        "percentile",
        "distance_bp",
    ])

    for p in PERCENTILES:

        writer.writerow([
            p,
            f"{percentile(distances, p):.1f}",
        ])


# ============================================================
# Report
# ============================================================

print(
    f"Reference species: "
    f"{args.reference_species}"
)

print(
    f"N CRE-gene associations: "
    f"{len(distances)}"
)

for p in PERCENTILES:

    print(
        f"{p:5g}th percentile: "
        f"{percentile(distances, p):.1f} bp"
    )

print(
    f"Wrote: {OUTFILE}"
)
