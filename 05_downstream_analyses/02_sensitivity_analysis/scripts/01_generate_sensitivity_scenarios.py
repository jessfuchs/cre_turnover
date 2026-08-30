#!/usr/bin/env python3

# ============================================================
# Generate CRE-classification sensitivity scenarios
#
# Purpose:
#   Re-run CRE-state classification across combinations of
#   reciprocal-overlap and local-gene-distance thresholds.
#
# Analysis:
#   Each target species is classified independently for every
#   sensitivity scenario using the standard CRE-classification
#   worker script.
#
# Validation:
#   - all required inputs must exist
#   - the reference CRE table must contain the expected number
#     of D. melanogaster reference CREs
#   - overlap and distance thresholds are validated
#   - each generated species table must contain one row per
#     reference CRE
#   - the complete scenario x species grid is checked after
#     classification
#
# Existing results:
#   Complete existing species outputs are reused unless
#   --force is supplied. Incomplete outputs are recomputed.
#
# Outputs:
#   - one classification directory per sensitivity scenario
#   - sensitivity scenario manifest
#   - run metadata
#
# Input/output paths and sensitivity parameters:
#   Supplied by the pipeline wrapper using
#   config/sensitivity_config.sh.
# ============================================================

import argparse
import csv
from datetime import datetime
import hashlib
from pathlib import Path
import platform
import subprocess
import sys


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Generate CRE-state classifications across "
            "reciprocal-overlap and local-gene-distance "
            "sensitivity scenarios."
        )
    )

    parser.add_argument(
        "--classifier",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--targets",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--reference",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--predictions",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--lifted-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--scenario-manifest",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--metadata-out",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--overlaps",
        type=float,
        nargs="+",
        required=True,
    )

    parser.add_argument(
        "--distances",
        type=int,
        nargs="+",
        required=True,
    )

    parser.add_argument(
        "--expected-reference-cres",
        type=int,
        required=True,
    )

    parser.add_argument(
        "--force",
        action="store_true",
    )

    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_file(
    path,
    label,
):
    """
    Abort if a required file is missing.
    """

    if not path.is_file():

        raise SystemExit(
            f"ERROR: {label} not found:\n"
            f"{path}"
        )


def require_directory(
    path,
    label,
):
    """
    Abort if a required directory is missing.
    """

    if not path.is_dir():

        raise SystemExit(
            f"ERROR: {label} not found:\n"
            f"{path}"
        )


def load_species(path):
    """
    Load unique target-species identifiers from a text file.
    """

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

        duplicates = sorted({
            species_name
            for species_name in species
            if species.count(
                species_name
            ) > 1
        })

        raise SystemExit(
            "ERROR: duplicate species in target file:\n"
            + "\n".join(
                duplicates
            )
        )

    return species


def count_data_rows(path):
    """
    Count non-empty data rows after the header line.
    """

    with path.open() as handle:

        return sum(
            1
            for index, line in enumerate(
                handle
            )
            if index > 0
            and line.strip()
        )


def file_sha256(path):
    """
    Calculate SHA256 checksum for reproducibility metadata.
    """

    sha = hashlib.sha256()

    with path.open(
        "rb"
    ) as handle:

        for block in iter(
            lambda:
                handle.read(
                    1024 * 1024
                ),
            b"",
        ):

            sha.update(
                block
            )

    return sha.hexdigest()


def scenario_name(
    overlap,
    distance,
):
    """
    Construct a deterministic sensitivity-scenario identifier.
    """

    overlap_int = int(
        round(
            overlap
            * 100
        )
    )

    return (
        f"ov{overlap_int:03d}"
        f"_dist{distance}"
    )


def build_scenarios(
    overlaps,
    distances,
):
    """
    Build the complete sensitivity grid.
    """

    scenarios = []

    for overlap in overlaps:

        for distance in distances:

            scenarios.append({
                "scenario":
                    scenario_name(
                        overlap,
                        distance,
                    ),

                "reciprocal_overlap":
                    overlap,

                "local_gene_distance_bp":
                    distance,
            })

    names = [
        scenario[
            "scenario"
        ]
        for scenario in scenarios
    ]

    if len(names) != len(set(names)):

        duplicates = sorted({
            name
            for name in names
            if names.count(
                name
            ) > 1
        })

        raise SystemExit(
            "ERROR: sensitivity parameters produce "
            "duplicate scenario names:\n"
            + "\n".join(
                duplicates
            )
        )

    return scenarios


def validate_parameters(args):
    """
    Validate sensitivity-analysis parameters.
    """

    if args.expected_reference_cres < 1:

        raise SystemExit(
            "ERROR: expected-reference-cres must be >= 1."
        )


    if len(args.overlaps) != len(set(args.overlaps)):

        raise SystemExit(
            "ERROR: duplicate overlap thresholds."
        )


    invalid_overlaps = [
        value
        for value in args.overlaps
        if value <= 0
        or value > 1
    ]

    if invalid_overlaps:

        raise SystemExit(
            "ERROR: overlap thresholds must be > 0 and <= 1:\n"
            + "\n".join(
                str(value)
                for value in invalid_overlaps
            )
        )


    if len(args.distances) != len(set(args.distances)):

        raise SystemExit(
            "ERROR: duplicate distance thresholds."
        )


    invalid_distances = [
        value
        for value in args.distances
        if value < 0
    ]

    if invalid_distances:

        raise SystemExit(
            "ERROR: distance thresholds must be >= 0:\n"
            + "\n".join(
                str(value)
                for value in invalid_distances
            )
        )


def write_scenario_manifest(
    scenarios,
    path,
    out_dir,
    n_species,
):
    """
    Write an explicit sensitivity-scenario definition table.
    """

    fieldnames = [
        "scenario",
        "reciprocal_overlap",
        "local_gene_distance_bp",
        "n_target_species",
        "expected_species_files",
        "scenario_directory",
    ]

    with path.open(
        "w",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        for scenario in scenarios:

            scenario_dir = (
                out_dir
                / scenario[
                    "scenario"
                ]
            )

            writer.writerow({
                "scenario":
                    scenario[
                        "scenario"
                    ],

                "reciprocal_overlap":
                    scenario[
                        "reciprocal_overlap"
                    ],

                "local_gene_distance_bp":
                    scenario[
                        "local_gene_distance_bp"
                    ],

                "n_target_species":
                    n_species,

                "expected_species_files":
                    n_species,

                "scenario_directory":
                    str(
                        scenario_dir.resolve()
                    ),
            })


def write_metadata(
    args,
    scenarios,
    species,
    n_reference_cres,
    expected_files,
    n_generated,
    n_skipped,
):
    """
    Write one-row run metadata for reproducibility.
    """

    metadata = {
        "script":
            Path(
                __file__
            ).name,

        "run_timestamp":
            datetime.now()
            .astimezone()
            .isoformat(),

        "python_version":
            sys.version.split()[0],

        "platform":
            platform.platform(),

        "classifier_script":
            str(
                args.classifier.resolve()
            ),

        "classifier_sha256":
            file_sha256(
                args.classifier
            ),

        "target_species_file":
            str(
                args.targets.resolve()
            ),

        "target_species_sha256":
            file_sha256(
                args.targets
            ),

        "reference_cre_file":
            str(
                args.reference.resolve()
            ),

        "reference_cre_sha256":
            file_sha256(
                args.reference
            ),

        "predictions_file":
            str(
                args.predictions.resolve()
            ),

        "predictions_sha256":
            file_sha256(
                args.predictions
            ),

        "lifted_cre_directory":
            str(
                args.lifted_dir.resolve()
            ),

        "output_directory":
            str(
                args.out_dir.resolve()
            ),

        "scenario_manifest":
            str(
                args.scenario_manifest.resolve()
            ),

        "scenario_manifest_sha256":
            file_sha256(
                args.scenario_manifest
            ),

        "expected_reference_cres":
            args.expected_reference_cres,

        "observed_reference_cres":
            n_reference_cres,

        "n_target_species":
            len(
                species
            ),

        "target_species":
            "|".join(
                species
            ),

        "n_overlap_thresholds":
            len(
                args.overlaps
            ),

        "overlap_thresholds":
            "|".join(
                f"{value:g}"
                for value in args.overlaps
            ),

        "n_distance_thresholds":
            len(
                args.distances
            ),

        "distance_thresholds_bp":
            "|".join(
                str(value)
                for value in args.distances
            ),

        "n_scenarios":
            len(
                scenarios
            ),

        "expected_scenario_files":
            expected_files,

        "n_generated_files":
            n_generated,

        "n_skipped_existing_files":
            n_skipped,

        "force_recompute":
            (
                "yes"
                if args.force
                else "no"
            ),

        "run_status":
            "complete",
    }


    with args.metadata_out.open(
        "w",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                metadata.keys()
            ),
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        writer.writerow(
            metadata
        )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    validate_parameters(
        args
    )


    # ========================================================
    # Validate primary inputs
    # ========================================================

    require_file(
        args.classifier,
        "CRE-classification script",
    )

    require_file(
        args.targets,
        "target-species file",
    )

    require_file(
        args.reference,
        "reference CRE table",
    )

    require_file(
        args.predictions,
        "combined SCRMshaw orthology table",
    )

    require_directory(
        args.lifted_dir,
        "lifted-CRE directory",
    )


    # ========================================================
    # Load target species
    # ========================================================

    species = load_species(
        args.targets
    )


    # ========================================================
    # Validate reference CRE count
    # ========================================================

    n_reference_cres = count_data_rows(
        args.reference
    )

    if (
        n_reference_cres
        != args.expected_reference_cres
    ):

        raise SystemExit(
            "ERROR: unexpected number of reference CREs.\n"
            f"Expected: {args.expected_reference_cres}\n"
            f"Observed: {n_reference_cres}"
        )


    # ========================================================
    # Validate all lifted-CRE inputs before starting
    # ========================================================

    lifted_files = {}

    for sp in species:

        lifted = (
            args.lifted_dir
            / sp
            / f"dmel_reference_cres.{sp}.bed"
        )

        require_file(
            lifted,
            f"lifted CRE file for {sp}",
        )

        lifted_files[
            sp
        ] = lifted


    # ========================================================
    # Build sensitivity grid
    # ========================================================

    scenarios = build_scenarios(
        args.overlaps,
        args.distances,
    )


    # ========================================================
    # Prepare output directories
    # ========================================================

    args.out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.scenario_manifest.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.metadata_out.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


    # ========================================================
    # Write scenario manifest
    # ========================================================

    write_scenario_manifest(
        scenarios,
        args.scenario_manifest,
        args.out_dir,
        len(
            species
        ),
    )


    # ========================================================
    # Run summary
    # ========================================================

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
        f"{args.overlaps}"
    )

    print(
        f"Distance thresholds:   "
        f"{args.distances}"
    )

    print(
        f"Scenarios:             "
        f"{len(scenarios)}"
    )

    print(
        f"Expected output files: "
        f"{len(scenarios) * len(species)}"
    )

    print(
        f"Force recomputation:   "
        f"{'yes' if args.force else 'no'}"
    )

    print(
        f"Output directory:      "
        f"{args.out_dir}"
    )

    print(
        f"Scenario manifest:     "
        f"{args.scenario_manifest}"
    )

    print()


    # ========================================================
    # Generate sensitivity classifications
    # ========================================================

    n_generated = 0
    n_skipped = 0


    for scenario_info in scenarios:

        scenario = scenario_info[
            "scenario"
        ]

        overlap = scenario_info[
            "reciprocal_overlap"
        ]

        distance = scenario_info[
            "local_gene_distance_bp"
        ]


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


        # ----------------------------------------------------
        # Run every target species
        # ----------------------------------------------------

        for index, sp in enumerate(
            species,
            start=1,
        ):

            lifted = lifted_files[
                sp
            ]

            outfile = (
                scenario_dir
                / f"dmel_to_{sp}_cre_turnover.tsv"
            )


            # ------------------------------------------------
            # Reuse complete existing output
            # ------------------------------------------------

            if (
                outfile.is_file()
                and not args.force
            ):

                n_rows = count_data_rows(
                    outfile
                )

                if (
                    outfile.stat().st_size > 0
                    and
                    n_rows == n_reference_cres
                ):

                    print(
                        f"[{index:02d}/{len(species):02d}] "
                        f"{sp}: already complete -> skip"
                    )

                    n_skipped += 1

                    continue


                print(
                    f"[{index:02d}/{len(species):02d}] "
                    f"{sp}: incomplete existing output "
                    f"({n_rows}/{n_reference_cres}) "
                    f"-> recompute"
                )


            elif (
                outfile.is_file()
                and args.force
            ):

                print(
                    f"[{index:02d}/{len(species):02d}] "
                    f"{sp}: force -> recompute"
                )


            else:

                print(
                    f"[{index:02d}/{len(species):02d}] "
                    f"{sp}: generate"
                )


            # ------------------------------------------------
            # Remove stale or incomplete output before rerun
            # ------------------------------------------------

            if outfile.is_file():

                outfile.unlink()


            # ------------------------------------------------
            # Classification command
            # ------------------------------------------------

            command = [
                sys.executable,
                str(
                    args.classifier
                ),

                "--species",
                sp,

                "--reference",
                str(
                    args.reference
                ),

                "--lifted",
                str(
                    lifted
                ),

                "--predictions",
                str(
                    args.predictions
                ),

                "--min-lifted-overlap",
                str(
                    overlap
                ),

                "--max-local-gene-distance",
                str(
                    distance
                ),

                "--out",
                str(
                    outfile
                ),
            ]


            result = subprocess.run(
                command,
                check=False,
            )


            if result.returncode != 0:

                raise SystemExit(
                    "\nERROR: CRE classification failed.\n"
                    f"Scenario: {scenario}\n"
                    f"Species:  {sp}\n"
                    f"Exit code: {result.returncode}"
                )


            # ------------------------------------------------
            # Validate generated species output
            # ------------------------------------------------

            if not outfile.is_file():

                raise SystemExit(
                    "ERROR: classifier finished successfully "
                    "but output file was not created:\n"
                    f"{outfile}"
                )


            n_rows = count_data_rows(
                outfile
            )


            if n_rows != n_reference_cres:

                raise SystemExit(
                    "\nERROR: unexpected number of "
                    "classification rows.\n"
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


    incomplete = []


    for scenario_info in scenarios:

        scenario = scenario_info[
            "scenario"
        ]

        for sp in species:

            outfile = (
                args.out_dir
                / scenario
                / f"dmel_to_{sp}_cre_turnover.tsv"
            )


            if not outfile.is_file():

                incomplete.append(
                    (
                        scenario,
                        sp,
                        "missing",
                    )
                )

                continue


            n_rows = count_data_rows(
                outfile
            )


            if n_rows != n_reference_cres:

                incomplete.append(
                    (
                        scenario,
                        sp,
                        (
                            f"{n_rows}/"
                            f"{n_reference_cres}"
                        ),
                    )
                )


    expected_files = (
        len(
            scenarios
        )
        * len(
            species
        )
    )


    if incomplete:

        print(
            "ERROR: sensitivity grid is incomplete:"
        )

        print()

        for (
            scenario,
            sp,
            status,
        ) in incomplete:

            print(
                f"{scenario}\t"
                f"{sp}\t"
                f"{status}"
            )

        raise SystemExit(1)


    # --------------------------------------------------------
    # Every expected file must have been either generated
    # during this run or reused as an existing complete file.
    # --------------------------------------------------------

    observed_files = (
        n_generated
        + n_skipped
    )

    if observed_files != expected_files:

        raise SystemExit(
            "ERROR: internal sensitivity-grid accounting "
            "is inconsistent.\n"
            f"Expected files:     {expected_files}\n"
            f"Generated + reused: {observed_files}"
        )


    # ========================================================
    # Write run metadata
    # ========================================================

    write_metadata(
        args,
        scenarios,
        species,
        n_reference_cres,
        expected_files,
        n_generated,
        n_skipped,
    )


    # ========================================================
    # Final report
    # ========================================================

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
        f"Scenario manifest:\n"
        f"{args.scenario_manifest}"
    )

    print()

    print(
        f"Run metadata:\n"
        f"{args.metadata_out}"
    )

    print()

    print(
        f"Results:\n"
        f"{args.out_dir}"
    )


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    main()
