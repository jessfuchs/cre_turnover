#!/usr/bin/env python3

from pathlib import Path
import csv
import sys

header = [
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

if len(sys.argv) != 3:
    sys.exit(
        "Usage: SO_bed_to_tsv.py COMBINED_RESULTS COMBINED_MANIFEST"
    )

results_dir = Path(sys.argv[1])
manifest_file = Path(sys.argv[2])

with manifest_file.open() as handle:
    manifest = {
        row["slug"].strip(): row
        for row in csv.DictReader(handle, delimiter="\t")
    }

print("\t".join(header))

seen = set()

for path in sorted(results_dir.glob("*/SO_all_peaks.bed")):

    slug = path.parent.name

    if slug not in manifest:
        print(
            f"WARNUNG: {slug} fehlt im Manifest",
            file=sys.stderr,
        )
        continue

    info = manifest[slug]
    species = info["species"].strip()
    source = info.get("source", "").strip()

    seen.add(slug)

    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):

            if not line.strip() or line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")

            if len(cols) != 18:
                print(
                    f"WARNUNG: {path}:{line_number}: "
                    f"{len(cols)} Spalten statt 18",
                    file=sys.stderr,
                )
                continue

            row = [
                slug,
                species,
                source,
                str(path.absolute()),
                *cols,
            ]

            print("\t".join(row))

missing = sorted(set(manifest) - seen)

for slug in missing:
    print(
        f"WARNUNG: keine SO_all_peaks.bed für {slug}",
        file=sys.stderr,
    )
