#!/usr/bin/env python3

import gzip
import re
import shutil
import tarfile
from pathlib import Path

import pandas as pd


PIPELINE_ROOT = Path(
    "/home/jessica/cre_turnover/project/scrmshaw_pipeline"
)

MANIFEST = (
    PIPELINE_ROOT /
    "data/external_manifest.tsv"
)

GENOME_ARCHIVE = Path(
    "/home/jessica/cre_turnover/project/"
    "scrmshaw_pipeline_6species/data/downloads/genomes.tar.gz"
)

OUTDIR = (
    PIPELINE_ROOT /
    "data/external_fastas"
)

OUT_MANIFEST = (
    OUTDIR /
    "external_fasta_manifest.tsv"
)


FASTA_EXTENSIONS = (
    ".fa",
    ".fasta",
    ".fna",
    ".fas",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
    ".fas.gz",
)


def normalize_taxon_name(species: str) -> str:
    """
    Drosophila affinis
    -> Drosophila_affinis
    """
    return re.sub(
        r"[^A-Za-z0-9]+",
        "_",
        species.strip()
    ).strip("_")


def genome_matches_species(
    basename: str,
    species: str,
) -> bool:

    taxon = normalize_taxon_name(
        species
    )

    # Expected archive style:
    # Drosophila_affinis.GCA_....rm.fna.gz
    #
    # Require exact taxon prefix followed by "."
    # to avoid accidentally matching subspecies.
    return bool(
        re.match(
            rf"^{re.escape(taxon)}\.",
            basename,
            flags=re.IGNORECASE,
        )
    )


# ============================================================
# Input checks
# ============================================================

if not MANIFEST.exists():
    raise SystemExit(
        f"ERROR: missing manifest:\n{MANIFEST}"
    )

if not GENOME_ARCHIVE.exists():
    raise SystemExit(
        f"ERROR: missing genome archive:\n{GENOME_ARCHIVE}"
    )


OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Read external manifest
# ============================================================

df = pd.read_csv(
    MANIFEST,
    sep="\t",
    dtype=str,
)


required = {
    "slug",
    "species",
}

missing = (
    required
    - set(df.columns)
)

if missing:
    raise SystemExit(
        "ERROR: external_manifest.tsv missing columns:\n"
        + "\n".join(sorted(missing))
    )


if df["slug"].duplicated().any():
    raise SystemExit(
        "ERROR: duplicated slugs in external_manifest.tsv"
    )


print(
    f"External species requested: {len(df)}"
)

for _, row in df.iterrows():

    print(
        f"  {row['slug']:16s}  "
        f"{row['species']}"
    )


# ============================================================
# Read genome archive
# ============================================================

print()
print(
    f"Opening archive:\n{GENOME_ARCHIVE}"
)


with tarfile.open(
    GENOME_ARCHIVE,
    "r:gz",
) as tar:

    genome_members = [
        member
        for member in tar.getmembers()
        if (
            member.isfile()
            and Path(
                member.name
            ).name.lower().endswith(
                FASTA_EXTENSIONS
            )
        )
    ]

    print(
        f"Genome-like files in archive: "
        f"{len(genome_members)}"
    )


    output_rows = []


    # ========================================================
    # Resolve and extract each species
    # ========================================================

    for _, row in df.iterrows():

        slug = row["slug"]
        species = row["species"]

        matches = []

        for member in genome_members:

            basename = Path(
                member.name
            ).name

            if genome_matches_species(
                basename,
                species,
            ):
                matches.append(
                    member
                )


        print()
        print(
            f"{slug}: {species}"
        )


        if len(matches) == 0:

            print(
                "  ERROR: no matching genome found"
            )

            output_rows.append({
                "slug": slug,
                "species": species,
                "status": "NO_MATCH",
                "archive_member": "",
                "fasta": "",
            })

            continue


        if len(matches) > 1:

            print(
                "  ERROR: multiple matching genomes:"
            )

            for m in matches:
                print(
                    f"    {m.name}"
                )

            output_rows.append({
                "slug": slug,
                "species": species,
                "status": "MULTIPLE_MATCHES",
                "archive_member": ";".join(
                    m.name
                    for m in matches
                ),
                "fasta": "",
            })

            continue


        member = matches[0]

        print(
            f"  archive member: "
            f"{member.name}"
        )


        source = tar.extractfile(
            member
        )

        if source is None:
            raise RuntimeError(
                f"Could not read archive member: "
                f"{member.name}"
            )


        output_path = (
            OUTDIR /
            f"{slug}.fa"
        )


        if member.name.lower().endswith(
            ".gz"
        ):

            with gzip.GzipFile(
                fileobj=source
            ) as src, output_path.open(
                "wb"
            ) as dst:

                shutil.copyfileobj(
                    src,
                    dst
                )

        else:

            with output_path.open(
                "wb"
            ) as dst:

                shutil.copyfileobj(
                    source,
                    dst
                )


        print(
            f"  wrote: {output_path}"
        )


        output_rows.append({
            "slug": slug,
            "species": species,
            "status": "OK",
            "archive_member": member.name,
            "fasta": str(output_path),
        })


# ============================================================
# Write extraction manifest
# ============================================================

out_df = pd.DataFrame(
    output_rows
)

out_df.to_csv(
    OUT_MANIFEST,
    sep="\t",
    index=False,
)


print()
print(
    "========================================"
)

print(
    f"Successful: "
    f"{(out_df['status'] == 'OK').sum()}"
    f"/{len(out_df)}"
)

print(
    f"Manifest:\n{OUT_MANIFEST}"
)

print(
    "========================================"
)


if not (
    out_df["status"] == "OK"
).all():

    raise SystemExit(
        "ERROR: not all species were uniquely resolved."
    )
