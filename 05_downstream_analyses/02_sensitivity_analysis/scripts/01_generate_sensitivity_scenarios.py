#!/usr/bin/env python3

from pathlib import Path
import argparse
import subprocess
import sys


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = (
    Path.home()
    / "cre_turnover"
    / "project"
)

WGA_DIR = (
    PROJECT_DIR
    / "pairwise_wga"
)

MAP_DIR = (
    PROJECT_DIR
    / "mapping_orthologs"
)

CLASS_DIR = (
    PROJECT_DIR
    / "cre_classification"
)

CLASS_SCRIPT = (
    CLASS_DIR
    / "scripts"
    / "01_classify_target_cres.py"
)

TARGETS = (
    WGA_DIR
    / "target_species.txt"
)

REFERENCE = (
    MAP_DIR
    / "reference_cres"
    / "dmel_reference_cres.tsv"
)

PREDICTIONS = (
    MAP_DIR
    / "ortholog_results"
    / "SO_all_species_fbgn.tsv"
)

LIFTED_BASE = (
    WGA_DIR
    / "lifted_cres_dmel"
)

DEFAULT_OUT_DIR = (
    CLASS_DIR
    / "sensitivity"
)


# ============================================================
# Sensitivity grid
# ============================================================

OVERLAPS = [
    0.25,
    0.50,
    0.75,
]

DISTANCES = [
    12000,
    24000,
    48000,
]


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Generate CRE-state classifications for all "
            "overlap x local-distance sensitivity scenarios."
        )
    )

    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Sensitivity output directory.",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Recompute files even if a non-empty output "
            "already exists."
        ),
    )

    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_file(path):

    if not path.is_file():

        raise SystemExit(
            "ERROR: required file not found:\n"
            f"{path}"
        )


def load_species(path):

    species = [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip()
        and not line.lstrip().startswith("#")
    ]

    if not species:

        raise SystemExit(
            "ERROR: target species file is empty:\n"
            f"{path}"
        )

    if len(species) != len(set(species)):

        raise SystemExit(
            "ERROR: duplicate species in target file."
        )

    return species


def count_data_rows(path):

    with path.open() as handle:

        return sum(
            1
            for i, line in enumerate(handle)
            if i > 0
            and line.strip()
        )


def scenario_name(
    overlap,
    distance,
):

    overlap_int = int(
        round(
            overlap * 100
        )
    )

    return (
        f"ov{overlap_int:03d}"
        f"_dist{distance}"
    )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    # --------------------------------------------------------
    # Check inputs
    # --------------------------------------------------------

    for path in [
        CLASS_SCRIPT,
        TARGETS,
        REFERENCE,
        PREDICTIONS,
    ]:
        require_file(path)

    species = load_species(
        TARGETS
    )

    n_reference_cres = count_data_rows(
        REFERENCE
    )

    if n_reference_cres == 0:

        raise SystemExit(
            "ERROR: reference CRE table contains no data rows."
        )

    args.out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    print()
    print("=" * 72)
    print("CRE CLASSIFICATION SENSITIVITY GRID")
    print("=" * 72)
    print()

    print(
        f"Target species:        "
        f"{len(species)}"
    )

    print(
        f"Reference CREs:        "
        f"{n_reference_cres}"
    )

    print(
        f"Overlap thresholds:    "
        f"{OVERLAPS}"
    )

    print(
        f"Distance thresholds:   "
        f"{DISTANCES}"
    )

    print(
        f"Scenarios:             "
        f"{len(OVERLAPS) * len(DISTANCES)}"
    )

    print(
        f"Output directory:      "
        f"{args.out_dir}"
    )

    print()


    n_generated = 0
    n_skipped = 0


    # ========================================================
    # Sensitivity scenarios
    # ========================================================

    for overlap in OVERLAPS:

        for distance in DISTANCES:

            scenario = scenario_name(
                overlap,
                distance,
            )

            scenario_dir = (
                args.out_dir
                / scenario
            )

            scenario_dir.mkdir(
                parents=True,
                exist_ok=True,
            )


            print()
            print("=" * 72)

            print(
                f"Scenario: {scenario}"
            )

            print(
                f"Reciprocal overlap: "
                f"{overlap:.2f}"
            )

            print(
                f"Local gene distance: "
                f"{distance} bp"
            )

            print("=" * 72)


            for index, sp in enumerate(
                species,
                start=1,
            ):

                lifted = (
                    LIFTED_BASE
                    / sp
                    / f"dmel_reference_cres.{sp}.bed"
                )

                require_file(
                    lifted
                )


                outfile = (
                    scenario_dir
                    / f"dmel_to_{sp}_cre_turnover.tsv"
                )


                # --------------------------------------------
                # Skip complete existing outputs
                # --------------------------------------------

                if (
                    outfile.is_file()
                    and outfile.stat().st_size > 0
                    and not args.force
                ):

                    n_rows = count_data_rows(
                        outfile
                    )

                    if n_rows == n_reference_cres:

                        print(
                            f"[{index:02d}/{len(species):02d}] "
                            f"{sp}: already complete -> skip"
                        )

                        n_skipped += 1

                        continue


                    print(
                        f"[{index:02d}/{len(species):02d}] "
                        f"{sp}: incomplete existing output "
                        f"({n_rows}/{n_reference_cres}) -> recompute"
                    )


                else:

                    print(
                        f"[{index:02d}/{len(species):02d}] "
                        f"{sp}"
                    )


                # --------------------------------------------
                # Classification command
                # --------------------------------------------

                command = [
                    sys.executable,
                    str(CLASS_SCRIPT),

                    "--species",
                    sp,

                    "--reference",
                    str(REFERENCE),

                    "--lifted",
                    str(lifted),

                    "--predictions",
                    str(PREDICTIONS),

                    "--min-lifted-overlap",
                    str(overlap),

                    "--max-local-gene-distance",
                    str(distance),

                    "--out",
                    str(outfile),
                ]


                result = subprocess.run(
                    command,
                    check=False,
                )


                if result.returncode != 0:

                    raise SystemExit(
                        "\nERROR: classification failed\n"
                        f"Scenario: {scenario}\n"
                        f"Species:  {sp}\n"
                        f"Exit code: {result.returncode}"
                    )


                # --------------------------------------------
                # Validate output
                # --------------------------------------------

                if not outfile.is_file():

                    raise SystemExit(
                        "ERROR: classifier finished but "
                        "output file was not created:\n"
                        f"{outfile}"
                    )


                n_rows = count_data_rows(
                    outfile
                )


                if n_rows != n_reference_cres:

                    raise SystemExit(
                        "\nERROR: unexpected number of rows\n"
                        f"Scenario: {scenario}\n"
                        f"Species:  {sp}\n"
                        f"Observed: {n_rows}\n"
                        f"Expected: {n_reference_cres}"
                    )


                n_generated += 1


    # ========================================================
    # Final completeness check
    # ========================================================

    print()
    print("=" * 72)
    print("FINAL COMPLETENESS CHECK")
    print("=" * 72)
    print()


    missing = []


    for overlap in OVERLAPS:

        for distance in DISTANCES:

            scenario = scenario_name(
                overlap,
                distance,
            )

            for sp in species:

                outfile = (
                    args.out_dir
                    / scenario
                    / f"dmel_to_{sp}_cre_turnover.tsv"
                )


                if not outfile.is_file():

                    missing.append(
                        f"{scenario}\t{sp}\tmissing"
                    )

                    continue


                n_rows = count_data_rows(
                    outfile
                )


                if n_rows != n_reference_cres:

                    missing.append(
                        f"{scenario}\t{sp}\t"
                        f"{n_rows}/{n_reference_cres}"
                    )


    expected_files = (
        len(OVERLAPS)
        * len(DISTANCES)
        * len(species)
    )


    if missing:

        print(
            "ERROR: sensitivity grid is incomplete:"
        )

        print()

        for item in missing:

            print(item)

        raise SystemExit(1)


    print(
        f"Expected scenario files: "
        f"{expected_files}"
    )

    print(
        f"Newly generated:        "
        f"{n_generated}"
    )

    print(
        f"Existing/skipped:       "
        f"{n_skipped}"
    )

    print()

    print(
        "PASS: all sensitivity classifications "
        "are complete."
    )

    print()

    print(
        f"Results:\n{args.out_dir}"
    )


if __name__ == "__main__":
    main()
