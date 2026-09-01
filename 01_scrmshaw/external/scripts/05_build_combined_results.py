#!/usr/bin/env python3

# ============================================================
# Combine generated and external SCRMshaw results
#
# This script:
#   1. reads the manifests for generated and external species,
#   2. combines them into a single manifest,
#   3. verifies that species slugs are unique,
#   4. creates a combined results directory, and
#   5. links each species-specific peaks_AllSets.bed file
#      into the combined results structure.
#
# Existing combined results are removed and rebuilt.
# ============================================================

from pathlib import Path
import argparse
import csv
import shutil


# ============================================================
# Arguments
# ============================================================

ap = argparse.ArgumentParser()

for name in [
    "generated_manifest",
    "external_manifest",
    "generated_results",
    "external_results",
    "combined_manifest",
    "combined_results",
]:
    ap.add_argument(
        "--" + name.replace("_", "-"),
        dest=name,
        type=Path,
        required=True,
    )

args = ap.parse_args()


# ============================================================
# Read manifests
# ============================================================

rows = []

with args.generated_manifest.open() as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        rows.append(
            (
                row["slug"].strip(),
                row["species"].strip(),
                "generated",
            )
        )

with args.external_manifest.open() as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        rows.append(
            (
                row["slug"].strip(),
                row["species"].strip(),
                "external",
            )
        )


# ============================================================
# Check for duplicate species slugs
# ============================================================

slugs = [row[0] for row in rows]

if len(slugs) != len(set(slugs)):
    raise SystemExit("ERROR: Duplicate species slugs")


# ============================================================
# Write combined manifest
# ============================================================

args.combined_manifest.parent.mkdir(
    parents=True,
    exist_ok=True,
)

with args.combined_manifest.open("w") as handle:
    handle.write("index\tslug\tspecies\tsource\n")

    for index, (slug, species, source) in enumerate(rows):
        handle.write(
            f"{index}\t{slug}\t{species}\t{source}\n"
        )


# ============================================================
# Prepare combined results directory
# ============================================================

if args.combined_results.exists():
    shutil.rmtree(args.combined_results)

args.combined_results.mkdir(parents=True)


# ============================================================
# Link species-specific peak files
# ============================================================

missing = []

for slug, species, source_type in rows:

    source_root = (
        args.generated_results
        if source_type == "generated"
        else args.external_results
    )

    source = (
        source_root
        / slug
        / "peaks_AllSets.bed"
    )

    if not source.is_file() or source.stat().st_size == 0:
        missing.append(source)
        continue

    destination = args.combined_results / slug
    destination.mkdir()

    (
        destination
        / "peaks_AllSets.bed"
    ).symlink_to(source.resolve())


# ============================================================
# Summary
# ============================================================

print(
    f"Manifest species={len(rows)} "
    f"Symlinks={len(rows) - len(missing)}"
)

if missing:
    print("Missing:")

    for path in missing:
        print(path)

    raise SystemExit(1)
