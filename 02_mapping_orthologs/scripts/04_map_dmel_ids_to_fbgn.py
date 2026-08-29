#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import sys
from collections import defaultdict


# ============================================================
# Arguments
# ============================================================

if len(sys.argv) != 4:
    sys.exit(
        "Usage: map_dmel_ids_to_fbgn.py "
        "SO_all_species.tsv drosophila_melanogaster.gff3 OUTPUT.tsv"
    )

input_tsv = Path(sys.argv[1])
gff = Path(sys.argv[2])
output_tsv = Path(sys.argv[3])

if not input_tsv.is_file():
    sys.exit(f"FEHLER: Input fehlt: {input_tsv}")

if not gff.is_file():
    sys.exit(f"FEHLER: D. melanogaster GFF fehlt: {gff}")

output_tsv.parent.mkdir(parents=True, exist_ok=True)

unresolved_file = (
    output_tsv.parent
    / "unresolved_dmel_identifiers.tsv"
)


# ============================================================
# Helper functions
# ============================================================

def parse_attrs(text):
    """
    Parse GFF3 attributes into a dictionary.
    """

    attrs = {}

    for item in text.strip().split(";"):

        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value

    return attrs


def get_fbgn(attrs):
    """
    Extract one or more FBgn identifiers from GFF attributes.

    The 301Fly D. melanogaster annotation stores FBgn identifiers
    primarily in dbxref.
    """

    dbxref = attrs.get("dbxref", "")

    return sorted(
        set(
            re.findall(
                r"FBgn\d+",
                dbxref
            )
        )
    )


# ============================================================
# Build Dmel identifier -> FBgn mapping
# ============================================================

id_to_fbgn = defaultdict(set)

with gff.open(
    encoding="utf-8",
    errors="replace"
) as handle:

    for line in handle:

        if not line.strip() or line.startswith("#"):
            continue

        cols = line.rstrip("\n").split("\t")

        if len(cols) < 9:
            continue

        attrs = parse_attrs(cols[8])

        fbgns = get_fbgn(attrs)

        if not fbgns:
            continue


        # ----------------------------------------------------
        # ID=
        #
        # Examples:
        #   ID=rna-NM_...
        #   ID=gene-Dmel_CG...
        # ----------------------------------------------------

        feature_id = attrs.get("ID")

        if feature_id:

            for fbgn in fbgns:
                id_to_fbgn[feature_id].add(fbgn)


        # ----------------------------------------------------
        # transcript_id=
        #
        # Store both:
        #   NM_...
        #   rna-NM_...
        # ----------------------------------------------------

        transcript_id = attrs.get("transcript_id")

        if transcript_id:

            for fbgn in fbgns:

                id_to_fbgn[transcript_id].add(fbgn)

                id_to_fbgn[
                    "rna-" + transcript_id
                ].add(fbgn)


        # ----------------------------------------------------
        # locus_tag=
        #
        # Store both:
        #   Dmel_CG...
        #   gene-Dmel_CG...
        # ----------------------------------------------------

        locus_tag = attrs.get("locus_tag")

        if locus_tag:

            for fbgn in fbgns:

                id_to_fbgn[locus_tag].add(fbgn)

                id_to_fbgn[
                    "gene-" + locus_tag
                ].add(fbgn)


# Convert sets to sorted lists
id_to_fbgn = {
    key: sorted(values)
    for key, values in id_to_fbgn.items()
}

print(
    f"Loaded Dmel identifier mappings: "
    f"{len(id_to_fbgn)}",
    file=sys.stderr
)


# ============================================================
# Map one field to FBgn
# ============================================================

def map_field(value):
    """
    Map one or multiple D. melanogaster identifiers to FBgn.

    Input examples:
        rna-NM_001...
        gene-Dmel_CG1234
        Dmel_CG1234
        id1|id2|id3

    Output:
        FBgn0000001
        FBgn0000001|FBgn0000002

    Returns
    -------
    mapped : str
        Pipe-separated FBgn identifiers or NA.

    identifiers : list
        Input identifiers that were evaluated.

    unresolved : list
        Identifiers for which no FBgn mapping was found.
    """

    if not value or value in {"NA", "."}:
        return "NA", [], []

    identifiers = [
        x.strip()
        for x in value.split("|")
        if x.strip()
        and x.strip() != "NA"
    ]

    fbgns = set()
    unresolved = []

    for ident in identifiers:

        hits = id_to_fbgn.get(ident)

        if hits:
            fbgns.update(hits)

        else:
            unresolved.append(ident)

    if fbgns:
        mapped = "|".join(
            sorted(fbgns)
        )

    else:
        mapped = "NA"

    return mapped, identifiers, unresolved


# ============================================================
# Read input and write FBgn-annotated TSV
# ============================================================

with input_tsv.open(
    encoding="utf-8-sig"
) as inp:

    reader = csv.DictReader(
        inp,
        delimiter="\t"
    )

    if reader.fieldnames is None:
        sys.exit(
            f"FEHLER: TSV hat keinen Header: "
            f"{input_tsv}"
        )

    required = {
        "dmel_ortholog_flanking_gene",
        "dmel_ortholog_next_flanking_gene",
    }

    missing = (
        required
        - set(reader.fieldnames)
    )

    if missing:
        sys.exit(
            "FEHLER: Fehlende Spalten: "
            + ", ".join(
                sorted(missing)
            )
        )


    # --------------------------------------------------------
    # Output columns
    # --------------------------------------------------------

    output_fields = list(
        reader.fieldnames
    )

    output_fields += [
        "dmel_fbgn_flanking_gene",
        "dmel_fbgn_next_flanking_gene",
    ]


    # --------------------------------------------------------
    # QC counters
    # --------------------------------------------------------

    total = 0

    flank_has_dmel_id = 0
    next_has_dmel_id = 0

    flank_fbgn_mapped = 0
    next_fbgn_mapped = 0

    unresolved_counts = defaultdict(int)


    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------

    with output_tsv.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as out:

        writer = csv.DictWriter(
            out,
            fieldnames=output_fields,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()


        for row in reader:

            total += 1


            # -----------------------------------------------
            # Flanking ortholog
            # -----------------------------------------------

            original_flank = (
                row[
                    "dmel_ortholog_flanking_gene"
                ].strip()
            )

            fbgn1, ids1, unresolved1 = (
                map_field(original_flank)
            )

            if ids1:
                flank_has_dmel_id += 1


            # -----------------------------------------------
            # Next ortholog
            # -----------------------------------------------

            original_next = (
                row[
                    "dmel_ortholog_next_flanking_gene"
                ].strip()
            )

            fbgn2, ids2, unresolved2 = (
                map_field(original_next)
            )

            if ids2:
                next_has_dmel_id += 1


            # -----------------------------------------------
            # Add FBgn columns
            # -----------------------------------------------

            row[
                "dmel_fbgn_flanking_gene"
            ] = fbgn1

            row[
                "dmel_fbgn_next_flanking_gene"
            ] = fbgn2


            if fbgn1 != "NA":
                flank_fbgn_mapped += 1

            if fbgn2 != "NA":
                next_fbgn_mapped += 1


            # -----------------------------------------------
            # Track unresolved identifiers
            # -----------------------------------------------

            for ident in (
                unresolved1
                + unresolved2
            ):
                unresolved_counts[
                    ident
                ] += 1


            writer.writerow(row)


# ============================================================
# Write unresolved identifiers
# ============================================================

with unresolved_file.open(
    "w",
    encoding="utf-8",
    newline=""
) as out:

    out.write(
        "identifier\toccurrences\n"
    )

    for ident, count in sorted(
        unresolved_counts.items(),
        key=lambda x: (
            -x[1],
            x[0]
        ),
    ):

        out.write(
            f"{ident}\t{count}\n"
        )


# ============================================================
# QC percentages
# ============================================================

if total > 0:

    flank_fbgn_pct = (
        100.0
        * flank_fbgn_mapped
        / total
    )

    next_fbgn_pct = (
        100.0
        * next_fbgn_mapped
        / total
    )

    flank_conversion_pct = (
        100.0
        * flank_fbgn_mapped
        / flank_has_dmel_id
        if flank_has_dmel_id > 0
        else 0.0
    )

    next_conversion_pct = (
        100.0
        * next_fbgn_mapped
        / next_has_dmel_id
        if next_has_dmel_id > 0
        else 0.0
    )

else:

    flank_fbgn_pct = 0.0
    next_fbgn_pct = 0.0
    flank_conversion_pct = 0.0
    next_conversion_pct = 0.0


# ============================================================
# QC report
# ============================================================

print(
    f"Rows: {total}",
    file=sys.stderr
)

print(
    f"Flanking with Dmel identifier: "
    f"{flank_has_dmel_id}",
    file=sys.stderr
)

print(
    f"Flanking with FBgn: "
    f"{flank_fbgn_mapped} "
    f"({flank_fbgn_pct:.1f}% of all rows; "
    f"{flank_conversion_pct:.1f}% of rows with Dmel ID)",
    file=sys.stderr
)

print(
    f"Next with Dmel identifier: "
    f"{next_has_dmel_id}",
    file=sys.stderr
)

print(
    f"Next with FBgn: "
    f"{next_fbgn_mapped} "
    f"({next_fbgn_pct:.1f}% of all rows; "
    f"{next_conversion_pct:.1f}% of rows with Dmel ID)",
    file=sys.stderr
)

print(
    f"Unique unresolved identifiers: "
    f"{len(unresolved_counts)}",
    file=sys.stderr
)

print(
    f"Unresolved identifier table: "
    f"{unresolved_file}",
    file=sys.stderr
)

print(
    f"Output: {output_tsv}",
    file=sys.stderr
)
