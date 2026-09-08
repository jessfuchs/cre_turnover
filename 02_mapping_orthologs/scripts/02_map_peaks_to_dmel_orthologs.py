#!/usr/bin/env python3


# ============================================================
# 02 - Map SCRMshaw peak-associated genes to D. melanogaster
#      orthologs
#
# Purpose:
#   Map the flanking genes of each SCRMshaw CRE prediction to
#   their D. melanogaster orthologs using the `dmel_orthologs`
#   attribute from species-specific GFF3 annotations.
#
# Input:
#   - combined_manifest.tsv
#   - combined_results/<species>/peaks_AllSets.bed
#   - species-specific GFF3 annotations
#
# Output:
#   - results/<species>/SO_all_peaks.bed
#   - results/ortholog_mapping_qc.tsv
#
# Configuration:
#   02_mapping_orthologs/config/ortholog_config.sh
# ============================================================

import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


# ============================================================
# Configuration
# ============================================================

def get_env_path(name):
    """
    Read a required path from an exported environment variable.
    """

    value = os.environ.get(name, "").strip()

    if not value:
        sys.exit(
            f"ERROR: required environment variable {name} is not set. "
            "Run this script through run_ortholog_pipeline.sh or source "
            "config/ortholog_config.sh before execution."
        )

    return Path(value)


COMBINED_MANIFEST = get_env_path("COMBINED_MANIFEST")
COMBINED_RESULTS_DIR = get_env_path("COMBINED_RESULTS_DIR")

GENERATED_GFF_ROOT = get_env_path("GENERATED_GFF_ROOT")
EXTERNAL_GFF_ROOT = get_env_path("EXTERNAL_GFF_ROOT")

RESULTS_DIR = get_env_path("RESULTS_DIR")


# ============================================================
# Helper functions
# ============================================================

def parse_attributes(text):
    """
    Parse a GFF3 attribute field into a dictionary.
    """

    attrs = {}

    for item in text.strip().split(";"):

        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value

    return attrs


def build_gene_to_dmel(gff):
    """
    Build a mapping from target-species genes to D. melanogaster
    ortholog identifiers.

    Ortholog assignments are derived from mRNA entries in the
    301Fly GFF3 annotations, for example:

        Parent=gene-G...
        dmel_orthologs=rna-NM_...|rna-NM_...

    Multiple orthologs associated with the same target gene are
    retained as a pipe-separated string.

    Returns
    -------
    dict
        target_gene -> Dmel ortholog identifier(s)
    """

    mapping = defaultdict(set)

    with gff.open(errors="replace") as handle:

        for line in handle:

            if not line.strip() or line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")

            if len(cols) < 9:
                continue

            if cols[2] != "mRNA":
                continue

            attrs = parse_attributes(cols[8])

            parent = attrs.get("Parent")
            dmel = attrs.get("dmel_orthologs")

            if not parent:
                continue

            # A transcript may theoretically contain multiple
            # parent-gene identifiers.
            parents = [
                x.strip()
                for x in parent.split(",")
                if x.strip()
            ]

            if not dmel or dmel == "NA":
                continue

            # Multiple Dmel orthologs may be separated by
            # "|" or ",".
            orthologs = [
                x.strip()
                for x in re.split(r"[|,]", dmel)
                if x.strip() and x.strip() != "NA"
            ]

            for gene in parents:
                for ortholog in orthologs:
                    mapping[gene].add(ortholog)

    return {
        gene: "|".join(sorted(values))
        for gene, values in mapping.items()
    }


def read_manifest(path):
    """
    Read the combined tab-separated species manifest.

    Rows are returned as a dictionary keyed by species slug.
    """

    if not path.is_file():
        sys.exit(
            f"ERROR: combined manifest does not exist: {path}"
        )

    result = {}

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        if reader.fieldnames is None:
            sys.exit(
                f"ERROR: combined manifest has no header: {path}"
            )

        required = {
            "slug",
            "species",
            "source",
        }

        missing = required - set(reader.fieldnames)

        if missing:
            sys.exit(
                "ERROR: combined manifest is missing required columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:

            slug = row.get("slug", "").strip()

            if not slug:
                continue

            result[slug] = row

    return result


# ============================================================
# Read combined species manifest
# ============================================================

combined = read_manifest(COMBINED_MANIFEST)

if not combined:
    sys.exit(
        f"ERROR: combined manifest is empty: "
        f"{COMBINED_MANIFEST}"
    )


# ============================================================
# Prepare output directory
# ============================================================

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# QC header
#
# stderr is redirected by run_ortholog_pipeline.sh to:
#
#   results/ortholog_mapping_qc.tsv
#
# Keeping stderr restricted to the table ensures that the QC
# output remains machine-readable.
# ============================================================

print(
    "species_key\t"
    "rows\t"
    "mapped_flanking\t"
    "mapped_next\t"
    "unmapped_flanking\t"
    "unmapped_next\t"
    "flanking_mapped_pct\t"
    "next_mapped_pct",
    file=sys.stderr
)


# ============================================================
# Process all species
# ============================================================

for slug, info in combined.items():

    species = info["species"].strip()
    source = info.get("source", "").strip()


    # --------------------------------------------------------
    # Locate species-specific GFF3 annotation
    # --------------------------------------------------------

    if source in {
        "generated",
        "original",
        "original22",
    }:

        gff = (
            GENERATED_GFF_ROOT
            / slug
            / "annotation.gff3"
        )

    elif source in {
        "external",
        "external12",
    }:

        gff = (
            EXTERNAL_GFF_ROOT
            / slug
            / "annotation.gff3"
        )

    else:

        print(
            f"WARNING: unknown source={source!r} for {slug}; "
            "species skipped.",
            file=sys.stdout
        )

        continue


    # --------------------------------------------------------
    # Locate standardized SCRMshaw predictions
    #
    # All generated and external predictions have already been
    # harmonized by 01_scrmshaw/external and are therefore read
    # from the unified combined_results directory.
    # --------------------------------------------------------

    peaks = (
        COMBINED_RESULTS_DIR
        / slug
        / "peaks_AllSets.bed"
    )


    # --------------------------------------------------------
    # Validate input files
    # --------------------------------------------------------

    if not gff.is_file():

        print(
            f"WARNING: GFF3 missing for {slug}: {gff}",
            file=sys.stdout
        )

        continue


    if not peaks.is_file():

        print(
            f"WARNING: peaks_AllSets.bed missing for {slug}: "
            f"{peaks}",
            file=sys.stdout
        )

        continue


    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------

    output_dir = (
        RESULTS_DIR
        / slug
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    output = (
        output_dir
        / "SO_all_peaks.bed"
    )


    # --------------------------------------------------------
    # D. melanogaster special case
    #
    # D. melanogaster is the reference species and therefore
    # generally does not require a dmel_orthologs annotation.
    #
    # For the reference species:
    #
    #   output column 7  = input column 6
    #   output column 12 = input column 11
    # --------------------------------------------------------

    is_dmel = (
        species.lower()
        == "drosophila melanogaster"
    )

    if is_dmel:
        gene_to_dmel = {}
    else:
        gene_to_dmel = build_gene_to_dmel(gff)


    # --------------------------------------------------------
    # Counters
    # --------------------------------------------------------

    total = 0

    mapped_flanking = 0
    mapped_next = 0

    unmapped_flanking = 0
    unmapped_next = 0


    # --------------------------------------------------------
    # Process SCRMshaw peaks
    # --------------------------------------------------------

    with peaks.open(errors="replace") as inp, \
         output.open("w") as out:

        for line_number, line in enumerate(
            inp,
            start=1
        ):

            if not line.strip():
                continue

            if line.startswith("#"):
                out.write(line)
                continue

            cols = line.rstrip("\n").split("\t")


            # ------------------------------------------------
            # Validate SCRMshaw BED structure
            # ------------------------------------------------

            if len(cols) != 18:

                print(
                    f"WARNING: {peaks}:{line_number} contains "
                    f"{len(cols)} columns instead of 18; "
                    "row skipped.",
                    file=sys.stdout
                )

                continue


            total += 1


            # ------------------------------------------------
            # SCRMshaw gene fields
            #
            # Column 6  = flanking gene
            # Column 11 = next flanking gene
            # ------------------------------------------------

            flanking_gene = cols[5]
            next_gene = cols[10]


            # ------------------------------------------------
            # Flanking gene -> Dmel ortholog
            #
            # Output column 7
            # ------------------------------------------------

            if is_dmel:

                cols[6] = flanking_gene

                if flanking_gene not in {
                    "",
                    "NA",
                    ".",
                }:
                    mapped_flanking += 1
                else:
                    unmapped_flanking += 1

            else:

                ortholog = gene_to_dmel.get(
                    flanking_gene
                )

                if ortholog:

                    cols[6] = ortholog
                    mapped_flanking += 1

                else:

                    cols[6] = "NA"
                    unmapped_flanking += 1


            # ------------------------------------------------
            # Next flanking gene -> Dmel ortholog
            #
            # Output column 12
            # ------------------------------------------------

            if is_dmel:

                cols[11] = next_gene

                if next_gene not in {
                    "",
                    "NA",
                    ".",
                }:
                    mapped_next += 1
                else:
                    unmapped_next += 1

            else:

                ortholog = gene_to_dmel.get(
                    next_gene
                )

                if ortholog:

                    cols[11] = ortholog
                    mapped_next += 1

                else:

                    cols[11] = "NA"
                    unmapped_next += 1


            out.write(
                "\t".join(cols)
                + "\n"
            )


    # --------------------------------------------------------
    # Mapping percentages
    # --------------------------------------------------------

    if total > 0:

        flanking_pct = (
            100.0
            * mapped_flanking
            / total
        )

        next_pct = (
            100.0
            * mapped_next
            / total
        )

    else:

        flanking_pct = 0.0
        next_pct = 0.0


    # --------------------------------------------------------
    # QC output
    #
    # Keep this strictly tab-separated because stderr is
    # redirected to ortholog_mapping_qc.tsv.
    # --------------------------------------------------------

    print(
        f"{slug}\t"
        f"{total}\t"
        f"{mapped_flanking}\t"
        f"{mapped_next}\t"
        f"{unmapped_flanking}\t"
        f"{unmapped_next}\t"
        f"{flanking_pct:.1f}\t"
        f"{next_pct:.1f}",
        file=sys.stderr
    )


    # --------------------------------------------------------
    # Informational output
    # --------------------------------------------------------

    print(
        f"{slug}: {output}",
        file=sys.stdout
    )
