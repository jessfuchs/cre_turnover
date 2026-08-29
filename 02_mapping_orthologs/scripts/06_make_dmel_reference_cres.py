#!/usr/bin/env python3

import csv
from pathlib import Path


# ============================================================
# Paths
# ============================================================

ROOT = Path.home() / "cre_turnover" / "project" / "mapping_orthologs"

INFILE = ROOT / "ortholog_results" / "SO_all_species_fbgn.tsv"

OUTDIR = ROOT / "reference_cres"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = OUTDIR / "dmel_reference_cres.tsv"
BEDFILE = OUTDIR / "dmel_reference_cres.bed"

DMEL_SLUG = "d_melanogaster"


# ============================================================
# Helper
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

    out = []

    for x in value.split("|"):
        x = x.strip()

        if x.startswith("FBgn") and x not in out:
            out.append(x)

    return out


# ============================================================
# Input validation
# ============================================================

if not INFILE.is_file():
    raise SystemExit(
        f"FEHLER: Input fehlt: {INFILE}"
    )


# ============================================================
# Read D. melanogaster CREs
# ============================================================

with INFILE.open(
    encoding="utf-8-sig"
) as f:

    reader = csv.DictReader(
        f,
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
        raise SystemExit(
            "FEHLER: Missing required columns: "
            + ", ".join(sorted(missing))
        )

    rows = []

    for row in reader:

        if row["species_key"].strip() != DMEL_SLUG:
            continue

        try:
            start = int(row["start0"])
            end = int(row["end0"])

        except ValueError:
            raise SystemExit(
                f"FEHLER: Invalid coordinates: "
                f"{row['chrom']}:{row['start0']}-{row['end0']}"
            )

        if end <= start:
            raise SystemExit(
                f"FEHLER: Invalid interval: "
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
    raise SystemExit(
        f"FEHLER: Keine CREs für species_key={DMEL_SLUG!r} gefunden."
    )


# ============================================================
# Stable ordering
# ============================================================

rows.sort(
    key=lambda x: (
        x["chrom"],
        x["start0"],
        x["end0"],
        x["training_set"],
        x["method"],
        int(x["rank"]),
    )
)


# ============================================================
# Stable CRE IDs
# ============================================================

for i, row in enumerate(rows, start=1):
    row["dmel_cre_id"] = (
        f"DMEL_CRE_{i:05d}"
    )


# ============================================================
# Output TSV
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


with OUTFILE.open(
    "w",
    encoding="utf-8",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fieldnames,
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(rows)


# ============================================================
# BED6 for liftOver
# ============================================================

with BEDFILE.open(
    "w",
    encoding="utf-8"
) as f:

    for row in rows:

        f.write(
            f"{row['chrom']}\t"
            f"{row['start0']}\t"
            f"{row['end0']}\t"
            f"{row['dmel_cre_id']}\t"
            f"0\t.\n"
        )


# ============================================================
# QC
# ============================================================

n_with_gene = sum(
    1
    for r in rows
    if r["n_fbgn_target_genes"] > 0
)

n_multi = sum(
    1
    for r in rows
    if r["n_fbgn_target_genes"] > 1
)

n_without_gene = (
    len(rows)
    - n_with_gene
)


print(
    f"Dmel CREs: {len(rows)}"
)

print(
    f"With >=1 FBgn target: "
    f"{n_with_gene}"
)

print(
    f"With >1 FBgn target: "
    f"{n_multi}"
)

print(
    f"Without FBgn target: "
    f"{n_without_gene}"
)

print(
    f"Wrote: {OUTFILE}"
)

print(
    f"Wrote: {BEDFILE}"
)
