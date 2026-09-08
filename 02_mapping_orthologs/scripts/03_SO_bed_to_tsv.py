#!/usr/bin/env python3

# ============================================================
# 03 - Combine species-specific ortholog-mapped SCRMshaw peaks
#
# Purpose:
#   Combine the species-specific SO_all_peaks.bed files into a
#   single tab-separated table and annotate each row with species
#   metadata from the combined manifest.
#
# Input:
#   - results/<species>/SO_all_peaks.bed
#   - combined_manifest.tsv
#
# Output:
#   - SO_all_species.tsv (written to stdout)
#   - warnings/QC messages (written to stderr)
#
# Configuration:
#   Input and output paths are supplied by run_ortholog_pipeline.sh
#   using config/ortholog_config.sh.
# ============================================================


from pathlib import Path
import csv
import sys


# ============================================================
# Output columns
# ============================================================

HEADER = [
    "species_key",
    "species",
    "source",
    "source_file",
    "chrom",
    "start0",
    "end0",
    "peak_amplitude",
    "scrmshaw_score",
    "flanking_gene",
    "dmel_ortholog_flanking_gene",
    "distance_flanking_gene",
    "location_flanking_gene",
    "local_rank_flanking_gene",
    "next_flanking_gene",
    "dmel_ortholog_next_flanking_gene",
    "distance_next_gene",
    "location_next_gene",
    "local_rank_next_gene",
    "training_set",
    "method",
    "rank",
]


# ============================================================
# Command-line arguments
# ============================================================

if len(sys.argv) != 3:
    sys.exit(
        "Usage: 03_SO_bed_to_tsv.py "
        "RESULTS_DIR COMBINED_MANIFEST"
    )

results_dir = Path(sys.argv[1])
manifest_file = Path(sys.argv[2])


# ============================================================
# Input checks
# ============================================================

if not results_dir.is_dir():
    sys.exit(
        f"ERROR: ortholog-results directory does not exist: "
        f"{results_dir}"
    )

if not manifest_file.is_file():
    sys.exit(
        f"ERROR: combined manifest does not exist: "
        f"{manifest_file}"
    )


# ============================================================
# Read species manifest
# ============================================================

with manifest_file.open(
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
        sys.exit(
            "ERROR: combined manifest is missing required columns: "
            + ", ".join(sorted(missing))
        )

    manifest = {
        row["slug"].strip(): row
        for row in reader
        if row.get("slug", "").strip()
    }


# ============================================================
# Write combined table
# ============================================================

print("\t".join(HEADER))

seen = set()

for path in sorted(
    results_dir.glob("*/SO_all_peaks.bed")
):

    slug = path.parent.name

    if slug not in manifest:
        print(
            f"WARNING: {slug} is not present in the combined manifest",
            file=sys.stderr,
        )
        continue

    info = manifest[slug]

    species = info["species"].strip()
    source = info.get("source", "").strip()

    seen.add(slug)

    # Store a portable path rather than an absolute server path.
    source_file = str(
        path.relative_to(results_dir)
    )

    with path.open(errors="replace") as handle:

        for line_number, line in enumerate(
            handle,
            start=1
        ):

            if not line.strip() or line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")

            if len(cols) != 18:
                print(
                    f"WARNING: {path}:{line_number} contains "
                    f"{len(cols)} columns instead of 18; row skipped",
                    file=sys.stderr,
                )
                continue

            row = [
                slug,
                species,
                source,
                source_file,
                *cols,
            ]

            print("\t".join(row))


# ============================================================
# Report species without an SO_all_peaks.bed file
# ============================================================

missing_species = sorted(
    set(manifest) - seen
)

for slug in missing_species:
    print(
        f"WARNING: no SO_all_peaks.bed found for {slug}",
        file=sys.stderr,
    )
