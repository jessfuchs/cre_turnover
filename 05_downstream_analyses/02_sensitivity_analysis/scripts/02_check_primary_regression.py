#!/usr/bin/env python3

# ============================================================
# Check primary sensitivity-scenario regression
#
# Purpose:
#   Verify that the sensitivity scenario corresponding to the
#   primary CRE-classification parameters exactly reproduces
#   the baseline CRE-state assignments.
#
# Comparison:
#   For every target species, baseline and sensitivity tables
#   are matched one-to-one by D. melanogaster reference CRE ID
#   and their CRE-state classifications are compared.
#
# Validation:
#   - all required input files must exist
#   - each species table must contain the expected number of
#     reference CREs
#   - CRE IDs must be unique within each table
#   - baseline and sensitivity tables must contain identical
#     CRE-ID sets
#   - CRE-state values must use the expected vocabulary
#   - the primary scenario must occur exactly once in the
#     sensitivity scenario manifest
#
# Outputs:
#   - per-species regression summary
#   - detailed table of changed CRE states
#   - run metadata
#
# The script exits with a non-zero status if any CRE state
# differs between the baseline and primary sensitivity run.
#
# Input/output paths:
#   Supplied by the pipeline wrapper using
#   config/sensitivity_config.sh.
# ============================================================

import argparse
from datetime import datetime
import hashlib
from pathlib import Path
import platform
import sys

import pandas as pd


# ============================================================
# CRE-state vocabulary
# ============================================================

VALID_STATES = {
    "present",
    "turnover_candidate",
    "no_detected_CRE",
    "uncertain",
}


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Verify that the primary sensitivity scenario "
            "exactly reproduces the baseline CRE-state "
            "classification."
        )
    )

    parser.add_argument(
        "--baseline-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--sensitivity-dir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--targets",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--scenario-manifest",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--primary-scenario",
        required=True,
    )

    parser.add_argument(
        "--out-summary",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--out-changes",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--metadata-out",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--expected-reference-cres",
        type=int,
        required=True,
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


def require_columns(
    df,
    required,
    label,
):
    """
    Verify that all required columns are present.
    """

    missing = sorted(
        set(required)
        - set(df.columns)
    )

    if missing:

        raise SystemExit(
            f"ERROR: required columns missing from {label}:\n"
            + "\n".join(
                missing
            )
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


def load_species(path):
    """
    Load unique target-species identifiers.
    """

    require_file(
        path,
        "target-species file",
    )

    species = [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip()
        and not line.lstrip().startswith("#")
    ]

    if not species:

        raise SystemExit(
            "ERROR: target-species file contains no species:\n"
            f"{path}"
        )


    if len(species) != len(set(species)):

        duplicates = sorted({
            sp
            for sp in species
            if species.count(
                sp
            ) > 1
        })

        raise SystemExit(
            "ERROR: duplicate species in target-species file:\n"
            + "\n".join(
                duplicates
            )
        )


    return species


def load_scenario_definition(
    path,
    primary_scenario,
):
    """
    Retrieve the primary scenario from the scenario manifest.
    """

    require_file(
        path,
        "sensitivity scenario manifest",
    )

    manifest = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )


    require_columns(
        manifest,
        {
            "scenario",
            "reciprocal_overlap",
            "local_gene_distance_bp",
        },
        "sensitivity scenario manifest",
    )


    if manifest[
        "scenario"
    ].duplicated().any():

        duplicated = (
            manifest.loc[
                manifest[
                    "scenario"
                ].duplicated(
                    keep=False
                ),
                "scenario",
            ]
            .drop_duplicates()
            .tolist()
        )

        raise SystemExit(
            "ERROR: duplicate scenario names in sensitivity "
            "manifest:\n"
            + "\n".join(
                sorted(
                    duplicated
                )
            )
        )


    matched = manifest.loc[
        manifest[
            "scenario"
        ]
        == primary_scenario
    ]


    if len(matched) != 1:

        raise SystemExit(
            "ERROR: primary scenario must occur exactly once "
            "in sensitivity manifest.\n"
            f"Scenario: {primary_scenario}\n"
            f"Matches:  {len(matched)}"
        )


    row = matched.iloc[
        0
    ]


    try:

        overlap = float(
            row[
                "reciprocal_overlap"
            ]
        )

        distance = int(
            row[
                "local_gene_distance_bp"
            ]
        )

    except ValueError as exc:

        raise SystemExit(
            "ERROR: invalid primary-scenario parameter values "
            "in sensitivity manifest."
        ) from exc


    return overlap, distance


def load_states(
    path,
    expected_reference_cres,
    label,
):
    """
    Load and validate one species-level CRE-state table.
    """

    require_file(
        path,
        label,
    )


    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )


    require_columns(
        df,
        {
            "dmel_cre_id",
            "class",
        },
        label,
    )


    # --------------------------------------------------------
    # Required identifiers and states
    # --------------------------------------------------------

    missing_cre_id = (
        df[
            "dmel_cre_id"
        ]
        .astype(str)
        .str.strip()
        == ""
    )

    if missing_cre_id.any():

        raise SystemExit(
            "ERROR: empty dmel_cre_id values in:\n"
            f"{path}"
        )


    missing_class = (
        df[
            "class"
        ]
        .astype(str)
        .str.strip()
        == ""
    )

    if missing_class.any():

        raise SystemExit(
            "ERROR: empty CRE-state values in:\n"
            f"{path}"
        )


    # --------------------------------------------------------
    # One row per reference CRE
    # --------------------------------------------------------

    if df[
        "dmel_cre_id"
    ].duplicated().any():

        duplicated = (
            df.loc[
                df[
                    "dmel_cre_id"
                ].duplicated(
                    keep=False
                ),
                "dmel_cre_id",
            ]
            .drop_duplicates()
            .tolist()
        )

        raise SystemExit(
            "ERROR: duplicate CRE IDs in:\n"
            f"{path}\n\n"
            + "\n".join(
                sorted(
                    duplicated
                )
            )
        )


    if len(df) != expected_reference_cres:

        raise SystemExit(
            "ERROR: unexpected number of CRE rows in:\n"
            f"{path}\n"
            f"Expected: {expected_reference_cres}\n"
            f"Observed: {len(df)}"
        )


    # --------------------------------------------------------
    # CRE-state vocabulary
    # --------------------------------------------------------

    observed_states = set(
        df[
            "class"
        ]
    )

    invalid_states = sorted(
        observed_states
        - VALID_STATES
    )

    if invalid_states:

        raise SystemExit(
            "ERROR: unexpected CRE-state values in:\n"
            f"{path}\n\n"
            + "\n".join(
                invalid_states
            )
        )


    return (
        df[
            [
                "dmel_cre_id",
                "class",
            ]
        ]
        .copy()
    )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()


    # ========================================================
    # Parameter validation
    # ========================================================

    if args.expected_reference_cres < 1:

        raise SystemExit(
            "ERROR: expected-reference-cres must be >= 1."
        )


    if not args.primary_scenario.strip():

        raise SystemExit(
            "ERROR: primary-scenario must not be empty."
        )


    # ========================================================
    # Input validation
    # ========================================================

    require_directory(
        args.baseline_dir,
        "baseline classification directory",
    )

    require_directory(
        args.sensitivity_dir,
        "sensitivity scenario directory",
    )

    require_file(
        args.targets,
        "target-species file",
    )

    require_file(
        args.scenario_manifest,
        "sensitivity scenario manifest",
    )


    species = load_species(
        args.targets
    )


    # ========================================================
    # Primary scenario definition
    # ========================================================

    (
        primary_overlap,
        primary_distance,
    ) = load_scenario_definition(
        args.scenario_manifest,
        args.primary_scenario,
    )


    scenario_dir = (
        args.sensitivity_dir
        / args.primary_scenario
    )


    require_directory(
        scenario_dir,
        "primary sensitivity-scenario directory",
    )


    # ========================================================
    # Prepare outputs
    # ========================================================

    args.out_summary.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.out_changes.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.metadata_out.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


    # ========================================================
    # Compare baseline and primary sensitivity scenario
    # ========================================================

    summary_rows = []
    changed_parts = []


    for index, sp in enumerate(
        species,
        start=1,
    ):

        baseline_file = (
            args.baseline_dir
            / f"dmel_to_{sp}_cre_turnover.tsv"
        )

        scenario_file = (
            scenario_dir
            / f"dmel_to_{sp}_cre_turnover.tsv"
        )


        baseline = load_states(
            baseline_file,
            args.expected_reference_cres,
            f"baseline classification for {sp}",
        )


        scenario = load_states(
            scenario_file,
            args.expected_reference_cres,
            f"primary sensitivity classification for {sp}",
        )


        # ----------------------------------------------------
        # One-to-one comparison by reference CRE
        # ----------------------------------------------------

        merged = baseline.merge(
            scenario,
            on="dmel_cre_id",
            how="outer",
            suffixes=(
                "_baseline",
                "_scenario",
            ),
            indicator=True,
            validate="one_to_one",
        )


        if not merged[
            "_merge"
        ].eq(
            "both"
        ).all():

            missing_baseline = sorted(
                merged.loc[
                    merged[
                        "_merge"
                    ]
                    == "right_only",
                    "dmel_cre_id",
                ]
                .astype(str)
                .tolist()
            )

            missing_scenario = sorted(
                merged.loc[
                    merged[
                        "_merge"
                    ]
                    == "left_only",
                    "dmel_cre_id",
                ]
                .astype(str)
                .tolist()
            )


            message = [
                f"ERROR: CRE-ID mismatch for species {sp}.",
            ]


            if missing_baseline:

                message.extend([
                    "",
                    "Missing from baseline:",
                    *missing_baseline,
                ])


            if missing_scenario:

                message.extend([
                    "",
                    "Missing from sensitivity scenario:",
                    *missing_scenario,
                ])


            raise SystemExit(
                "\n".join(
                    message
                )
            )


        merged = merged.drop(
            columns=[
                "_merge",
            ]
        )


        # ----------------------------------------------------
        # State differences
        # ----------------------------------------------------

        changed = (
            merged.loc[
                merged[
                    "class_baseline"
                ]
                != merged[
                    "class_scenario"
                ]
            ]
            .copy()
        )


        n_cre = len(
            merged
        )

        n_changed = len(
            changed
        )

        n_identical = (
            n_cre
            - n_changed
        )

        percent_identical = (
            100.0
            * n_identical
            / n_cre
        )


        summary_rows.append({
            "species":
                sp,

            "n_cre":
                n_cre,

            "n_identical":
                n_identical,

            "n_changed":
                n_changed,

            "percent_identical":
                percent_identical,
        })


        # ----------------------------------------------------
        # Detailed changed-state records
        # ----------------------------------------------------

        if not changed.empty:

            changed.insert(
                1,
                "species",
                sp,
            )


            changed[
                "state_transition"
            ] = (
                changed[
                    "class_baseline"
                ]
                + "->"
                + changed[
                    "class_scenario"
                ]
            )


            changed_parts.append(
                changed[
                    [
                        "dmel_cre_id",
                        "species",
                        "class_baseline",
                        "class_scenario",
                        "state_transition",
                    ]
                ]
            )


        print(
            f"[{index:02d}/{len(species):02d}] "
            f"{sp}: "
            f"{n_identical}/{n_cre} identical"
        )


    # ========================================================
    # Build regression summary
    # ========================================================

    summary = pd.DataFrame(
        summary_rows
    )


    if summary.empty:

        raise SystemExit(
            "ERROR: no species-level regression results generated."
        )


    # ========================================================
    # Build changed-state table
    # ========================================================

    if changed_parts:

        changes = pd.concat(
            changed_parts,
            ignore_index=True,
        )

    else:

        changes = pd.DataFrame(
            columns=[
                "dmel_cre_id",
                "species",
                "class_baseline",
                "class_scenario",
                "state_transition",
            ]
        )


    # ========================================================
    # Final regression statistics
    # ========================================================

    total_comparisons = int(
        summary[
            "n_cre"
        ].sum()
    )


    total_identical = int(
        summary[
            "n_identical"
        ].sum()
    )


    total_changed = int(
        summary[
            "n_changed"
        ].sum()
    )


    overall_percent_identical = (
        100.0
        * total_identical
        / total_comparisons
    )


    expected_comparisons = (
        len(
            species
        )
        * args.expected_reference_cres
    )


    if (
        total_comparisons
        != expected_comparisons
    ):

        raise SystemExit(
            "ERROR: unexpected number of CRE x species "
            "comparisons.\n"
            f"Expected: {expected_comparisons}\n"
            f"Observed: {total_comparisons}"
        )


    # ========================================================
    # Stable output ordering
    # ========================================================

    summary = (
        summary
        .reset_index(
            drop=True
        )
    )


    if not changes.empty:

        species_order = {
            sp:
                index
            for index, sp in enumerate(
                species
            )
        }


        changes[
            "_species_rank"
        ] = (
            changes[
                "species"
            ]
            .map(
                species_order
            )
        )


        changes = (
            changes
            .sort_values(
                [
                    "_species_rank",
                    "dmel_cre_id",
                    "class_baseline",
                    "class_scenario",
                ],
                kind="mergesort",
            )
            .drop(
                columns=[
                    "_species_rank",
                ]
            )
            .reset_index(
                drop=True
            )
        )


    # ========================================================
    # Write outputs
    # ========================================================

    summary.to_csv(
        args.out_summary,
        sep="\t",
        index=False,
    )


    changes.to_csv(
        args.out_changes,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Run metadata
    # ========================================================

    regression_status = (
        "PASS"
        if total_changed == 0
        else "FAIL"
    )


    metadata = pd.DataFrame([
        {
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

            "pandas_version":
                pd.__version__,

            "platform":
                platform.platform(),

            "baseline_directory":
                str(
                    args.baseline_dir.resolve()
                ),

            "sensitivity_directory":
                str(
                    args.sensitivity_dir.resolve()
                ),

            "primary_scenario":
                args.primary_scenario,

            "primary_reciprocal_overlap":
                primary_overlap,

            "primary_local_gene_distance_bp":
                primary_distance,

            "target_species_file":
                str(
                    args.targets.resolve()
                ),

            "target_species_sha256":
                file_sha256(
                    args.targets
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

            "n_target_species":
                len(
                    species
                ),

            "expected_cre_species_comparisons":
                expected_comparisons,

            "observed_cre_species_comparisons":
                total_comparisons,

            "n_identical_states":
                total_identical,

            "n_changed_states":
                total_changed,

            "percent_identical":
                overall_percent_identical,

            "regression_status":
                regression_status,

            "summary_output":
                str(
                    args.out_summary.resolve()
                ),

            "changes_output":
                str(
                    args.out_changes.resolve()
                ),
        }
    ])


    metadata.to_csv(
        args.metadata_out,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Console summary
    # ========================================================

    print()
    print("=" * 72)
    print("PRIMARY REGRESSION CHECK")
    print("=" * 72)
    print()


    print(
        summary.to_string(
            index=False,
        )
    )


    print()

    print(
        f"Primary scenario:      "
        f"{args.primary_scenario}"
    )

    print(
        f"Reciprocal overlap:    "
        f"{primary_overlap:g}"
    )

    print(
        f"Local gene distance:   "
        f"{primary_distance} bp"
    )

    print(
        f"Target species:        "
        f"{len(species)}"
    )

    print(
        f"Reference CREs:        "
        f"{args.expected_reference_cres}"
    )

    print(
        f"Total CRE x species:   "
        f"{total_comparisons}"
    )

    print(
        f"Identical states:      "
        f"{total_identical}"
    )

    print(
        f"Changed states:        "
        f"{total_changed}"
    )

    print(
        f"Percent identical:     "
        f"{overall_percent_identical:.4f}%"
    )


    print()
    print(
        f"Wrote regression summary:\n"
        f"{args.out_summary}"
    )

    print()

    print(
        f"Wrote changed states:\n"
        f"{args.out_changes}"
    )

    print()

    print(
        f"Wrote metadata:\n"
        f"{args.metadata_out}"
    )


    # ========================================================
    # Regression result
    # ========================================================

    if total_changed:

        print()
        print(
            "FAIL: primary sensitivity scenario does not "
            "reproduce the baseline classification."
        )

        raise SystemExit(1)


    print()
    print(
        f"PASS: {args.primary_scenario} reproduces "
        "all baseline CRE states."
    )


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    main()
