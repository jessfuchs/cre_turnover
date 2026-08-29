#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import sys
from collections import defaultdict


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path.home() / "cre_turnover" / "project"

GENERATED_ROOT = PROJECT_ROOT / "scrmshaw_pipeline"
EXTERNAL_ROOT = PROJECT_ROOT / "external_scrmshaw"
MAPPING_ROOT = PROJECT_ROOT / "mapping_orthologs"

GENERATED_MANIFEST = GENERATED_ROOT / "data" / "manifest.tsv"
EXTERNAL_MANIFEST = EXTERNAL_ROOT / "external_data" / "external_manifest.tsv"
COMBINED_MANIFEST = EXTERNAL_ROOT / "combined_manifest.tsv"

GENERATED_RESULTS = GENERATED_ROOT / "results"
EXTERNAL_RESULTS = EXTERNAL_ROOT / "external_results"

# Ortholog-annotierte BEDs werden separat gespeichert.
ORTHOLOG_RESULTS = MAPPING_ROOT / "ortholog_results"


# ============================================================
# Helper functions
# ============================================================

def parse_attributes(text):
    """
    Parse GFF3 attribute field into a dictionary.
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
    Build target-gene -> D. melanogaster ortholog mapping.

    The mapping is derived from mRNA entries in the 301Fly GFF:

        Parent=gene-G...
        dmel_orthologs=rna-NM_...|rna-NM_...

    Returns
    -------
    dict
        {
            target_gene:
                "dmel_ortholog1|dmel_ortholog2"
        }
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

            # Parent can theoretically contain multiple genes.
            parents = [
                x.strip()
                for x in parent.split(",")
                if x.strip()
            ]

            if not dmel or dmel == "NA":
                continue

            # Multiple orthologs can occur, for example:
            # rna-NM_x|rna-NM_y
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
    Read a tab-separated manifest and return rows keyed by slug.
    """

    if not path.is_file():
        return {}

    result = {}

    with path.open(encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            return {}

        for row in reader:

            slug = row.get("slug", "").strip()

            if not slug:
                continue

            result[slug] = row

    return result


# ============================================================
# Read manifests
# ============================================================

generated = read_manifest(GENERATED_MANIFEST)
external = read_manifest(EXTERNAL_MANIFEST)
combined = read_manifest(COMBINED_MANIFEST)

if not combined:
    sys.exit(
        f"FEHLER: Combined manifest fehlt oder ist leer: "
        f"{COMBINED_MANIFEST}"
    )


# ============================================================
# Prepare output directory
# ============================================================

ORTHOLOG_RESULTS.mkdir(parents=True, exist_ok=True)


# ============================================================
# QC header
#
# stderr can be redirected to:
# ortholog_mapping_qc.tsv
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
    # Determine annotation and peaks_AllSets.bed
    # --------------------------------------------------------

    if source in {"generated", "original", "original22"}:

        if slug not in generated:
            print(
                f"WARNUNG: {slug} fehlt im Generated-Manifest",
                file=sys.stderr
            )
            continue

        row = generated[slug]

        # Current generated manifest uses "annotation".
        # "gff" is retained as a fallback.
        gff_text = (
            row.get("annotation", "").strip()
            or row.get("gff", "").strip()
        )

        if not gff_text:
            print(
                f"WARNUNG: Keine Annotation für {slug}",
                file=sys.stderr
            )
            continue

        gff = Path(gff_text)

        peaks = (
            GENERATED_RESULTS
            / slug
            / "peaks_AllSets.bed"
        )

    elif source in {"external", "external12"}:

        if slug not in external:
            print(
                f"WARNUNG: {slug} fehlt im External-Manifest",
                file=sys.stderr
            )
            continue

        row = external[slug]

        # Current external manifest uses "gff".
        # "annotation" is retained as a fallback.
        gff_text = (
            row.get("gff", "").strip()
            or row.get("annotation", "").strip()
        )

        if not gff_text:
            print(
                f"WARNUNG: Keine Annotation für {slug}",
                file=sys.stderr
            )
            continue

        gff = Path(gff_text)

        peaks = (
            EXTERNAL_RESULTS
            / slug
            / "peaks_AllSets.bed"
        )

    else:
        print(
            f"WARNUNG: Unbekannte source={source!r} für {slug}",
            file=sys.stderr
        )
        continue


    # --------------------------------------------------------
    # Validate input files
    # --------------------------------------------------------

    if not gff.is_file():
        print(
            f"WARNUNG: GFF fehlt für {slug}: {gff}",
            file=sys.stderr
        )
        continue

    if not peaks.is_file():
        print(
            f"WARNUNG: peaks_AllSets.bed fehlt für {slug}: "
            f"{peaks}",
            file=sys.stderr
        )
        continue


    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------

    output_dir = ORTHOLOG_RESULTS / slug
    output_dir.mkdir(parents=True, exist_ok=True)

    output = output_dir / "SO_all_peaks.bed"


    # --------------------------------------------------------
    # D. melanogaster special case
    #
    # D. melanogaster itself usually does not require
    # dmel_orthologs= annotation.
    #
    # Therefore:
    #   column 7  = column 6
    #   column 12 = column 11
    # --------------------------------------------------------

    is_dmel = species.lower() == "drosophila melanogaster"

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
    # Process peaks
    # --------------------------------------------------------

    with peaks.open(errors="replace") as inp, \
         output.open("w") as out:

        for line_number, line in enumerate(inp, start=1):

            if not line.strip():
                continue

            if line.startswith("#"):
                out.write(line)
                continue

            cols = line.rstrip("\n").split("\t")

            if len(cols) != 18:
                print(
                    f"WARNUNG: {peaks}:{line_number} "
                    f"hat {len(cols)} statt 18 Spalten",
                    file=sys.stderr
                )
                continue

            total += 1

            # SCRMshaw:
            # column 6  = flanking gene
            # column 11 = next flanking gene
            flanking_gene = cols[5]
            next_gene = cols[10]


            # ------------------------------------------------
            # Flanking gene -> Dmel ortholog
            # Output column 7
            # ------------------------------------------------

            if is_dmel:

                cols[6] = flanking_gene

                if flanking_gene not in {"", "NA", "."}:
                    mapped_flanking += 1
                else:
                    unmapped_flanking += 1

            else:

                ortholog = gene_to_dmel.get(flanking_gene)

                if ortholog:
                    cols[6] = ortholog
                    mapped_flanking += 1

                else:
                    cols[6] = "NA"
                    unmapped_flanking += 1


            # ------------------------------------------------
            # Next flanking gene -> Dmel ortholog
            # Output column 12
            # ------------------------------------------------

            if is_dmel:

                cols[11] = next_gene

                if next_gene not in {"", "NA", "."}:
                    mapped_next += 1
                else:
                    unmapped_next += 1

            else:

                ortholog = gene_to_dmel.get(next_gene)

                if ortholog:
                    cols[11] = ortholog
                    mapped_next += 1

                else:
                    cols[11] = "NA"
                    unmapped_next += 1


            out.write("\t".join(cols) + "\n")


    # --------------------------------------------------------
    # Mapping percentages
    # --------------------------------------------------------

    if total > 0:

        flanking_pct = (
            100.0 * mapped_flanking / total
        )

        next_pct = (
            100.0 * mapped_next / total
        )

    else:

        flanking_pct = 0.0
        next_pct = 0.0


    # --------------------------------------------------------
    # QC output
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
    # Standard output
    # --------------------------------------------------------

    print(
        f"{slug}: {output}",
        file=sys.stdout
    )
