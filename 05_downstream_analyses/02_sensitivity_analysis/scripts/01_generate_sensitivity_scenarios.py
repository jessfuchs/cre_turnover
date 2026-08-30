#!/usr/bin/env python3

# ============================================================
# Generate CRE-classification sensitivity scenarios
#
# Purpose:
#   Re-run CRE-state classification across combinations of
#   reciprocal-overlap and local-gene-distance thresholds.
#
# Existing complete outputs are reused unless --force is used.
#
# Input/output paths and sensitivity parameters are supplied
# by config/sensitivity_config.sh via the pipeline wrapper.
# ============================================================

from pathlib import Path
from datetime import datetime
import argparse
import csv
import subprocess
import sys


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Generate CRE-state classifications across "
            "overlap and local-distance sensitivity scenarios."
        )
    )

    parser.add_argument("--classifier", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--lifted-dir", type=Path, required=True)

    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--scenario-manifest", type=Path, required=True)
    parser.add_argument("--metadata-out", type=Path, required=True)

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

def require_file(path):

    if not path.is_file():
        raise SystemExit(
            f"ERROR: required file not found:\n{path}"
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
            f"ERROR: no target species found in:\n{path}"
        )

    if len(species) != len(set(species)):
        raise SystemExit(
            "ERROR: duplicate target species."
        )

    return species


def count_data_rows(path):

    with path.open() as handle:
        return sum(
            1
            for i, line in enumerate(handle)
            if i > 0 and line.strip()
        )


def scenario_name(overlap, distance):

    overlap_percent = int(
        round(overlap * 100)
    )

    return (
        f"ov{overlap_percent:03d}"
        f"_dist{distance}"
    )


def write_manifest(
    scenarios,
    out_path,
):

    with out_path.open(
        "w",
        newline="",
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "scenario",
            "reciprocal_overlap",
            "local_gene_distance_bp",
        ])

        for scenario, overlap, distance in scenarios:
            writer.writerow([
                scenario,
                overlap,
                distance,
            ])


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    # --------------------------------------------------------
    # Validate inputs
    # --------------------------------------------------------

    for path in [
        args.classifier,
        args.targets,
        args.reference,
        args.predictions,
    ]:
        require_file(path)

    if not args.lifted_dir.is_dir():
        raise SystemExit(
            f"ERROR: lifted-CRE directory not found:\n"
            f"{args.lifted_dir}"
        )

    if any(
        value <= 0 or value > 1
        for value in args.overlaps
    ):
        raise SystemExit(
            "ERROR: overlap thresholds must be > 0 and <= 1."
        )

    if any(
        value < 0
        for value in args.distances
    ):
        raise SystemExit(
            "ERROR: distance thresholds must be >= 0."
        )

    if len(args.overlaps) != len(set(args.overlaps)):
        raise SystemExit(
            "ERROR: duplicate overlap thresholds."
        )

    if len(args.distances) != len(set(args.distances)):
        raise SystemExit(
            "ERROR: duplicate distance thresholds."
        )

    species = load_species(
        args.targets
    )

    n_reference_cres = count_data_rows(
        args.reference
    )

    if n_reference_cres != args.expected_reference_cres:
        raise SystemExit(
            "ERROR: unexpected number of reference CREs.\n"
            f"Expected: {args.expected_reference_cres}\n"
            f"Observed: {n_reference_cres}"
        )

    # --------------------------------------------------------
    # Build sensitivity grid
    # --------------------------------------------------------

    scenarios = [
        (
            scenario_name(overlap, distance),
            overlap,
            distance,
        )
        for overlap in args.overlaps
        for distance in args.distances
    ]

    scenario_names = [
        scenario
        for scenario, _, _ in scenarios
    ]

    if len(scenario_names) != len(set(scenario_names)):
        raise SystemExit(
            "ERROR: duplicate scenario names generated."
        )

    # --------------------------------------------------------
    # Check lifted CRE files before starting
    # --------------------------------------------------------

    lifted_files = {}

    for sp in species:

        path = (
            args.lifted_dir
            / sp
            / f"dmel_reference_cres.{sp}.bed"
        )

        require_file(path)
        lifted_files[sp] = path

    # --------------------------------------------------------
    # Prepare outputs
    # --------------------------------------------------------

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

    write_manifest(
        scenarios,
        args.scenario_manifest,
    )

    # --------------------------------------------------------
    # Run sensitivity grid
    # --------------------------------------------------------

    print()
    print("=" * 72)
    print("CRE CLASSIFICATION SENSITIVITY GRID")
    print("=" * 72)
    print(f"Target species:      {len(species)}")
    print(f"Reference CREs:      {n_reference_cres}")
    print(f"Scenarios:           {len(scenarios)}")
    print(f"Overlap thresholds:  {args.overlaps}")
    print(f"Distance thresholds: {args.distances}")
    print()

    n_generated = 0
    n_skipped = 0

    for scenario, overlap, distance in scenarios:

        scenario_dir = (
            args.out_dir
            / scenario
        )

        scenario_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        print("=" * 72)
        print(
            f"{scenario}: overlap={overlap}, "
            f"distance={distance} bp"
        )
        print("=" * 72)

        for index, sp in enumerate(
            species,
            start=1,
        ):

            outfile = (
                scenario_dir
                / f"dmel_to_{sp}_cre_turnover.tsv"
            )

            # ------------------------------------------------
            # Reuse complete output
            # ------------------------------------------------

            if (
                outfile.is_file()
                and not args.force
                and count_data_rows(outfile)
                == n_reference_cres
            ):

                print(
                    f"[{index:02d}/{len(species):02d}] "
                    f"{sp}: complete -> skip"
                )

                n_skipped += 1
                continue

            print(
                f"[{index:02d}/{len(species):02d}] "
                f"{sp}: classify"
            )

            if outfile.is_file():
                outfile.unlink()

            # ------------------------------------------------
            # Run standard classifier
            # ------------------------------------------------

            command = [
                sys.executable,
                str(args.classifier),

                "--species",
                sp,

                "--reference",
                str(args.reference),

                "--lifted",
                str(lifted_files[sp]),

                "--predictions",
                str(args.predictions),

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
                    "ERROR: classification failed.\n"
                    f"Scenario: {scenario}\n"
                    f"Species:  {sp}\n"
                    f"Exit code: {result.returncode}"
                )

            if (
                not outfile.is_file()
                or count_data_rows(outfile)
                != n_reference_cres
            ):
                raise SystemExit(
                    "ERROR: incomplete sensitivity output.\n"
                    f"Scenario: {scenario}\n"
                    f"Species:  {sp}"
                )

            n_generated += 1

    # --------------------------------------------------------
    # Final completeness check
    # --------------------------------------------------------

    missing = []

    for scenario, _, _ in scenarios:
        for sp in species:

            outfile = (
                args.out_dir
                / scenario
                / f"dmel_to_{sp}_cre_turnover.tsv"
            )

            if (
                not outfile.is_file()
                or count_data_rows(outfile)
                != n_reference_cres
            ):
                missing.append(
                    f"{scenario}\t{sp}"
                )

    if missing:
        raise SystemExit(
            "ERROR: incomplete sensitivity grid:\n"
            + "\n".join(missing)
        )

    # --------------------------------------------------------
    # Metadata
    # --------------------------------------------------------

    metadata = {
        "script":
            Path(__file__).name,

        "run_timestamp":
            datetime.now().astimezone().isoformat(),

        "n_target_species":
            len(species),

        "n_reference_cres":
            n_reference_cres,

        "n_scenarios":
            len(scenarios),

        "n_expected_files":
            len(scenarios) * len(species),

        "n_generated_files":
            n_generated,

        "n_skipped_files":
            n_skipped,

        "force":
            args.force,
    }

    with args.metadata_out.open(
        "w",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=metadata.keys(),
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerow(metadata)

    # --------------------------------------------------------
    # Finish
    # --------------------------------------------------------

    print()
    print("=" * 72)
    print("Sensitivity grid complete")
    print("=" * 72)
    print(f"Generated: {n_generated}")
    print(f"Skipped:   {n_skipped}")
    print(f"Results:   {args.out_dir}")
    print()


if __name__ == "__main__":
    main()
