#!/usr/bin/env python3

# ============================================================
# 06 - Build the D. melanogaster reference CRE set
#
# Purpose:
#   Extract D. melanogaster SCRMshaw CREs from the combined
#   FBgn-normalized table and generate stable reference CRE IDs
#   and associated target-gene annotations.
#
# Input:
#   - SO_all_species_fbgn.tsv
#
# Output:
#   - reference_cres/dmel_reference_cres.tsv
#   - reference_cres/dmel_reference_cres.bed
#   - QC summary written to stdout
#
# Configuration:
#   Paths are supplied by run_ortholog_pipeline.sh using
#   config/ortholog_config.sh.
# ============================================================


import csv
import sys
from pathlib import Path


# ============================================================
# Command-line arguments
# ============================================================

if len(sys.argv) != 5:
    sys.exit(
        "Usage: 06_make_dmel_reference_cres.py "
        "SO_all_species_fbgn.tsv "
        "OUTPUT.tsv OUTPUT.bed DMEL_SLUG"
    )

input_file = Path(sys.argv[1])
output_tsv = Path(sys.argv[2])
output_bed = Path(sys.argv[3])
dmel_slug = sys.argv[4]


# ============================================================
# Input validation
# ============================================================

if not input_file.is_file():
    sys.exit(
        f"ERROR: input TSV does not exist: {input_file}"
    )

output_tsv.parent.mkdir(
    parents=True,
    exist_ok=True
)

output_bed.parent.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helper functions
# ============================================================

def split_fbgn(value):
    """
    Return unique valid FBgn identifiers from a pipe-separated field.
    """

    if value is None:
        return []

    value = value.strip()

    if not value or value in {"NA", ".", "None"}:
        return []

    fbgns = []

    for item in value.split("|"):

        item = item.strip()

        if item.startswith("FBgn") and item not in fbgns:
            fbgns.append(item)

    return fbgns


# ============================================================
# Read D. melanogaster CREs
# ============================================================

with input_file.open(
    encoding="utf-8-sig",
    newline=""
) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )

    required = {
        "species_key",
        "chrom",
        "start0",
        "end0",
        "peak_amplitude",
        "scrmshaw_score",
        "training_set",
        "method",
        "rank",
        "dmel_fbgn_flanking_gene",
        "dmel_fbgn_next_flanking_gene",
    }

    missing = (
        required
        - set(reader.fieldnames or [])
    )

    if missing:
        sys.exit(
            "ERROR: missing required columns: "
            + ", ".join(sorted(missing))
        )

    rows = []

    for row in reader:

        if row["species_key"].strip() != dmel_slug:
            continue

        try:
            start = int(row["start0"])
            end = int(row["end0"])

        except ValueError:
            sys.exit(
                f"ERROR: invalid coordinates: "
                f"{row['chrom']}:{row['start0']}-{row['end0']}"
            )

        if end <= start:
            sys.exit(
                f"ERROR: invalid interval: "
                f"{row['chrom']}:{start}-{end}"
            )

        genes = []

        for field in (
            "dmel_fbgn_flanking_gene",
            "dmel_fbgn_next_flanking_gene",
        ):

            for fbgn in split_fbgn(row[field]):

                if fbgn not in genes:
                    genes.append(fbgn)

        rows.append(
            {
                "chrom": row["chrom"].strip(),
                "start0": start,
                "end0": end,
                "peak_amplitude": row["peak_amplitude"].strip(),
                "scrmshaw_score": row["scrmshaw_score"].strip(),
                "fbgn_target_genes": (
                    "|".join(genes)
                    if genes
                    else "NA"
                ),
                "n_fbgn_target_genes": len(genes),
                "training_set": row["training_set"].strip(),
                "method": row["method"].strip(),
                "rank": row["rank"].strip(),
            }
        )


# ============================================================
# Validate reference set
# ============================================================

if not rows:
    sys.exit(
        f"ERROR: no CREs found for species_key={dmel_slug!r}"
    )


# ============================================================
# Stable ordering
# ============================================================

rows.sort(
    key=lambda row: (
        row["chrom"],
        row["start0"],
        row["end0"],
        row["training_set"],
        row["method"],
        int(row["rank"]),
    )
)


# ============================================================
# Assign stable CRE identifiers
# ============================================================

for index, row in enumerate(
    rows,
    start=1
):
    row["dmel_cre_id"] = (
        f"DMEL_CRE_{index:05d}"
    )


# ============================================================
# Write reference CRE TSV
# ============================================================

fieldnames = [
    "dmel_cre_id",
    "chrom",
    "start0",
    "end0",
    "peak_amplitude",
    "scrmshaw_score",
    "fbgn_target_genes",
    "n_fbgn_target_genes",
    "training_set",
    "method",
    "rank",
]

with output_tsv.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=fieldnames,
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(rows)


# ============================================================
# Write BED6 reference file for liftOver
# ============================================================

with output_bed.open(
    "w",
    encoding="utf-8"
) as handle:

    for row in rows:

        handle.write(
            f"{row['chrom']}\t"
            f"{row['start0']}\t"
            f"{row['end0']}\t"
            f"{row['dmel_cre_id']}\t"
            f"0\t.\n"
        )


# ============================================================
# QC summary
# ============================================================

n_with_gene = sum(
    row["n_fbgn_target_genes"] > 0
    for row in rows
)

n_multi = sum(
    row["n_fbgn_target_genes"] > 1
    for row in rows
)

n_without_gene = (
    len(rows)
    - n_with_gene
)

print(f"Dmel CREs: {len(rows)}")
print(f"With >=1 FBgn target: {n_with_gene}")
print(f"With >1 FBgn target: {n_multi}")
print(f"Without FBgn target: {n_without_gene}")
print(f"Wrote: {output_tsv}")
print(f"Wrote: {output_bed}")
