#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv
import sys
from dataclasses import dataclass


# ======================================================================
# Data structures
# ======================================================================

@dataclass(frozen=True)
class Species:
    slug: str
    scientific_name: str


@dataclass
class BedStats:
    n_rows: int
    n_rank1: int
    min_rank: int
    max_rank: int
    seqids: set[str]


# ======================================================================
# Helpers
# ======================================================================

def fail(message: str) -> None:
    raise SystemExit(f"FEHLER: {message}")


def warn(message: str) -> None:
    print(f"WARNUNG: {message}", file=sys.stderr)


# ======================================================================
# Species table
# ======================================================================

def read_species_file(path: Path) -> list[Species]:

    if not path.is_file() or path.stat().st_size == 0:
        fail(f"Species-Datei fehlt oder ist leer: {path}")

    species = []

    seen_slugs = set()
    seen_bed_ids = set()

    with path.open(errors="replace") as handle:

        for line_no, raw_line in enumerate(handle, start=1):

            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            # Split only once:
            # ["d_affinis", "Drosophila affinis"]
            fields = line.split(maxsplit=1)

            if len(fields) != 2:
                fail(
                    f"{path}:{line_no}: erwartet 2 Felder "
                    f"(slug, scientific_name), erhalten: {line}"
                )

            slug, scientific_name = fields

            scientific_parts = scientific_name.split()

            if len(scientific_parts) != 2:
                fail(
                    f"{path}:{line_no}: wissenschaftlicher Name muss "
                    f"'Genus species' sein: {scientific_name}"
                )

            if slug in seen_slugs:
                fail(
                    f"{path}:{line_no}: doppelter Slug: {slug}"
                )

            seen_slugs.add(slug)

            species.append(
                Species(
                    slug=slug,
                    scientific_name=scientific_name,
                )
            )

    if not species:
        fail(f"Keine Arten in {path} gefunden")

    return species


# ======================================================================
# BED validation
# ======================================================================

def inspect_bed(
    path: Path,
    expected_training: str,
    expected_method: str,
) -> BedStats:

    n_rows = 0
    n_rank1 = 0

    min_rank = None
    max_rank = None

    seqids = set()

    with path.open(errors="replace") as handle:

        for line_no, raw_line in enumerate(handle, start=1):

            if not raw_line.strip() or raw_line.startswith("#"):
                continue

            fields = raw_line.rstrip("\n").split("\t")

            if len(fields) != 17:
                fail(
                    f"{path}:{line_no}: erwartet 17 BED-Felder, "
                    f"gefunden: {len(fields)}"
                )

            training = fields[14]
            method = fields[15]

            if training != expected_training:
                fail(
                    f"{path}:{line_no}: Training set ist "
                    f"'{training}' statt '{expected_training}'"
                )

            if method != expected_method:
                fail(
                    f"{path}:{line_no}: Methode ist "
                    f"'{method}' statt '{expected_method}'"
                )

            try:
                rank = int(fields[16])
            except ValueError:
                fail(
                    f"{path}:{line_no}: Rank in Spalte 17 "
                    f"ist nicht ganzzahlig: {fields[16]}"
                )

            if rank < 1:
                fail(
                    f"{path}:{line_no}: ungültiger Rank: {rank}"
                )

            seqids.add(fields[0])

            n_rows += 1

            if rank == 1:
                n_rank1 += 1

            if min_rank is None or rank < min_rank:
                min_rank = rank

            if max_rank is None or rank > max_rank:
                max_rank = rank

    if n_rows == 0:
        fail(f"BED-Datei ist leer: {path}")

    return BedStats(
        n_rows=n_rows,
        n_rank1=n_rank1,
        min_rank=min_rank,
        max_rank=max_rank,
        seqids=seqids,
    )


# ======================================================================
# GFF validation
# ======================================================================

def seqids_gff(path: Path) -> set[str]:

    ids = set()

    with path.open(errors="replace") as handle:

        for raw_line in handle:

            if not raw_line.strip() or raw_line.startswith("#"):
                continue

            fields = raw_line.rstrip("\n").split("\t")

            if len(fields) >= 9:
                ids.add(fields[0])

    if not ids:
        fail(f"Keine SeqIDs im GFF gefunden: {path}")

    return ids


# ======================================================================
# Resolve files
# ======================================================================

def resolve_filtered_bed(
    bed_dir: Path,
    species: Species,
    training: str,
    method: str,
) -> Path:

    # Preferred new naming:
    #
    # daff.adult_muscle.imm.bed
    preferred = (
        bed_dir
        / f"{species.slug}.{training}.{method}.bed"
    )

    if preferred.is_file() and preferred.stat().st_size > 0:
        return preferred


def resolve_gff(
    gff_dir: Path,
    species: Species,
) -> Path:

    # New preferred naming generated by 01_extract_external_gffs.sh
    preferred = gff_dir / f"{species.slug}.gff3"

    if preferred.is_file() and preferred.stat().st_size > 0:
        return preferred

    # Compatibility with old extraction format:
    genus, epithet = species.scientific_name.split()

    legacy = (
        gff_dir
        / f"{genus.lower()}_{epithet.lower()}.gff3"
    )

    if legacy.is_file() and legacy.stat().st_size > 0:

        warn(
            f"{species.slug}: verwende Legacy-GFF-Namen "
            f"'{legacy.name}'. Besser wäre '{preferred.name}'."
        )

        return legacy

    fail(
        f"GFF fehlt für {species.scientific_name}.\n"
        f"  erwartet: {preferred}\n"
        f"  alternativ: {legacy}"
    )


# ======================================================================
# Main
# ======================================================================

def main():

    ap = argparse.ArgumentParser(
        description=(
            "Validate external SCRMshaw BEDs and matching GFF files "
            "and create a slug-based manifest."
        )
    )

    ap.add_argument(
        "--bed-dir",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--gff-dir",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--species-file",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--manifest",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--training",
        default="adult_muscle",
    )

    ap.add_argument(
        "--method",
        default="imm",
    )

    ap.add_argument(
        "--expected-offsets",
        type=int,
        default=25,
        help="Expected number of SCRMshaw offsets (default: 25)",
    )

    ap.add_argument(
        "--max-rank",
        type=int,
        default=5000,
        help="Maximum allowed rank per offset (default: 5000)",
    )

    args = ap.parse_args()

    species_list = read_species_file(args.species_file)

    print("=" * 72)
    print("External SCRMshaw validation")
    print("=" * 72)
    print(f"Species        : {len(species_list)}")
    print(f"Training set   : {args.training}")
    print(f"Method         : {args.method}")
    print(f"Expected offsets: {args.expected_offsets}")
    print(f"Maximum rank   : {args.max_rank}")
    print()

    manifest_rows = []

    for sp in species_list:

        print("-" * 72)
        print(f"{sp.slug}: {sp.scientific_name}")

        bed = resolve_filtered_bed(
            args.bed_dir,
            sp,
            args.training,
            args.method,
        )

        gff = resolve_gff(
            args.gff_dir,
            sp,
        )

        stats = inspect_bed(
            bed,
            args.training,
            args.method,
        )

        # --------------------------------------------------------------
        # Offset / rank structure
        # --------------------------------------------------------------

        if stats.min_rank != 1:
            fail(
                f"{sp.slug}: min rank = {stats.min_rank}; "
                f"erwartet wurde 1"
            )

        if stats.n_rank1 != args.expected_offsets:
            fail(
                f"{sp.slug}: Rank 1 kommt {stats.n_rank1}x vor; "
                f"erwartet wurden {args.expected_offsets} Offsets"
            )

        if stats.max_rank > args.max_rank:
            fail(
                f"{sp.slug}: max rank = {stats.max_rank}; "
                f"erlaubt sind maximal {args.max_rank}"
            )

        # We deliberately do NOT require exactly 125,000 rows.
        #
        # D_affinis showed:
        # 124,389 rows
        # rank1 = 25
        # max rank = 4,982
        #
        # Slightly fewer than 5,000 candidates per offset are therefore
        # valid and should not be discarded.
        theoretical_max = (
            args.expected_offsets
            * args.max_rank
        )

        if stats.n_rows > theoretical_max:
            fail(
                f"{sp.slug}: {stats.n_rows} Zeilen überschreiten "
                f"das theoretische Maximum von {theoretical_max}"
            )

        print(f"  BED            : {bed}")
        print(f"  GFF            : {gff}")
        print(f"  rows           : {stats.n_rows}")
        print(f"  rank1          : {stats.n_rank1}")
        print(f"  min rank       : {stats.min_rank}")
        print(f"  max rank       : {stats.max_rank}")
        print(
            f"  theoretical max: {theoretical_max}"
        )

        # --------------------------------------------------------------
        # BED / GFF sequence-ID compatibility
        # --------------------------------------------------------------

        gff_ids = seqids_gff(gff)

        missing = stats.seqids - gff_ids

        print(f"  BED SeqIDs     : {len(stats.seqids)}")
        print(f"  GFF SeqIDs     : {len(gff_ids)}")
        print(f"  missing SeqIDs : {len(missing)}")

        if missing:

            example = ", ".join(
                sorted(missing)[:10]
            )

            fail(
                f"{sp.slug}: BED/GFF-Mismatch. "
                f"{len(missing)} BED-SeqIDs fehlen im GFF. "
                f"Beispiele: {example}"
            )

        manifest_rows.append(
            (
                sp.slug,
                sp.scientific_name,
                str(bed.absolute()),
                str(gff.absolute()),
            )
        )

    # ------------------------------------------------------------------
    # Write manifest
    # ------------------------------------------------------------------

    args.manifest.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.manifest.open(
        "w",
        newline="",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "index",
                "slug",
                "species",
                "bed",
                "gff",
            ]
        )

        for index, row in enumerate(manifest_rows):

            writer.writerow(
                [index, *row]
            )

    print()
    print("=" * 72)
    print("Validation successful")
    print("=" * 72)
    print(f"Manifest : {args.manifest}")
    print(f"Species  : {len(manifest_rows)}")


if __name__ == "__main__":
    main()
