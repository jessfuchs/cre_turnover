#!/usr/bin/env python3

# ============================================================
# Summarize global CRE-classification sensitivity
#
# Purpose:
#   Quantify how CRE-state classifications change across the
#   complete sensitivity grid relative to the primary analysis
#   scenario.
#
# Analysis:
#   For each sensitivity scenario, CRE states are compared with
#   the primary scenario across all target species.
#
# Outputs:
#   - global scenario-level stability summary
#   - species-level stability summary
#   - CRE-level stability summary
#   - state-transition summary
#   - run metadata
#
# Validation:
#   - scenarios are read from the scenario manifest generated
#     in Step 01
#   - all required scenario x species files must exist
#   - each file must contain exactly the expected number of
#     reference CREs
#   - CRE IDs must be unique within species tables
#   - all scenario tables must contain identical CRE x species
#     combinations
#   - only expected CRE-state values are accepted
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
            "Summarize CRE-state stability across the complete "
            "sensitivity grid relative to the primary scenario."
        )
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
        "--out-global",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--out-species",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--out-transitions",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--out-cre",
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
            "ERROR: no target species found in:\n"
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


def load_scenario_manifest(
    path,
    primary_scenario,
):
    """
    Load and validate sensitivity-scenario definitions.
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


    if manifest.empty:

        raise SystemExit(
            "ERROR: sensitivity scenario manifest is empty."
        )


    # --------------------------------------------------------
    # Scenario names
    # --------------------------------------------------------

    missing_scenario = (
        manifest[
            "scenario"
        ]
        .astype(str)
        .str.strip()
        == ""
    )

    if missing_scenario.any():

        raise SystemExit(
            "ERROR: empty scenario names in sensitivity "
            "scenario manifest."
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
            "ERROR: duplicate sensitivity scenarios:\n"
            + "\n".join(
                sorted(
                    duplicated
                )
            )
        )


    # --------------------------------------------------------
    # Parameter values
    # --------------------------------------------------------

    try:

        manifest[
            "reciprocal_overlap"
        ] = pd.to_numeric(
            manifest[
                "reciprocal_overlap"
            ],
            errors="raise",
        )

        manifest[
            "local_gene_distance_bp"
        ] = pd.to_numeric(
            manifest[
                "local_gene_distance_bp"
            ],
            errors="raise",
        )

    except ValueError as exc:

        raise SystemExit(
            "ERROR: invalid numeric values in sensitivity "
            "scenario manifest."
        ) from exc


    invalid_overlap = (
        (
            manifest[
                "reciprocal_overlap"
            ]
            <= 0
        )
        |
        (
            manifest[
                "reciprocal_overlap"
            ]
            > 1
        )
    )

    if invalid_overlap.any():

        raise SystemExit(
            "ERROR: reciprocal-overlap thresholds in scenario "
            "manifest must be > 0 and <= 1."
        )


    non_integer_distance = (
        manifest[
            "local_gene_distance_bp"
        ]
        % 1
        != 0
    )

    if non_integer_distance.any():

        raise SystemExit(
            "ERROR: local-gene-distance thresholds must be "
            "integer base-pair values."
        )


    manifest[
        "local_gene_distance_bp"
    ] = manifest[
        "local_gene_distance_bp"
    ].astype(
        int
    )


    if (
        manifest[
            "local_gene_distance_bp"
        ]
        < 0
    ).any():

        raise SystemExit(
            "ERROR: local-gene-distance thresholds must be "
            ">= 0."
        )


    # --------------------------------------------------------
    # Primary scenario
    # --------------------------------------------------------

    primary_matches = (
        manifest[
            "scenario"
        ]
        == primary_scenario
    )

    if int(
        primary_matches.sum()
    ) != 1:

        raise SystemExit(
            "ERROR: primary sensitivity scenario must occur "
            "exactly once in scenario manifest.\n"
            f"Scenario: {primary_scenario}\n"
            f"Matches:  {int(primary_matches.sum())}"
        )


    return manifest.reset_index(
        drop=True
    )


def load_species_state_table(
    path,
    species,
    expected_reference_cres,
):
    """
    Load and validate one scenario x species CRE-state table.
    """

    require_file(
        path,
        f"CRE classification for {species}",
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
        str(
            path
        ),
    )


    # --------------------------------------------------------
    # Required values
    # --------------------------------------------------------

    missing_cre = (
        df[
            "dmel_cre_id"
        ]
        .astype(str)
        .str.strip()
        == ""
    )

    if missing_cre.any():

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
    # Exactly one row per reference CRE
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

    unknown = sorted(
        set(
            df[
                "class"
            ]
        )
        - VALID_STATES
    )

    if unknown:

        raise SystemExit(
            "ERROR: unexpected CRE-state values in:\n"
            f"{path}\n\n"
            + "\n".join(
                unknown
            )
        )


    out = (
        df[
            [
                "dmel_cre_id",
                "class",
            ]
        ]
        .copy()
    )

    out[
        "species"
    ] = species


    return out


def load_scenario(
    sensitivity_dir,
    scenario,
    species,
    expected_reference_cres,
):
    """
    Load the complete CRE x species table for one scenario.
    """

    scenario_dir = (
        sensitivity_dir
        / scenario
    )


    require_directory(
        scenario_dir,
        f"sensitivity scenario directory '{scenario}'",
    )


    parts = []


    for sp in species:

        path = (
            scenario_dir
            / f"dmel_to_{sp}_cre_turnover.tsv"
        )


        part = load_species_state_table(
            path,
            sp,
            expected_reference_cres,
        )

        parts.append(
            part
        )


    out = pd.concat(
        parts,
        ignore_index=True,
    )


    # --------------------------------------------------------
    # Unique CRE x species cells
    # --------------------------------------------------------

    if out.duplicated(
        subset=[
            "dmel_cre_id",
            "species",
        ]
    ).any():

        duplicated = (
            out.loc[
                out.duplicated(
                    subset=[
                        "dmel_cre_id",
                        "species",
                    ],
                    keep=False,
                ),
                [
                    "dmel_cre_id",
                    "species",
                ],
            ]
            .drop_duplicates()
        )

        raise SystemExit(
            "ERROR: duplicate CRE x species rows in "
            f"scenario {scenario}:\n\n"
            + duplicated.to_string(
                index=False
            )
        )


    expected_rows = (
        len(
            species
        )
        * expected_reference_cres
    )


    if len(out) != expected_rows:

        raise SystemExit(
            "ERROR: unexpected CRE x species count for "
            f"scenario {scenario}.\n"
            f"Expected: {expected_rows}\n"
            f"Observed: {len(out)}"
        )


    return out


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
    # Validate inputs
    # ========================================================

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


    scenario_manifest = load_scenario_manifest(
        args.scenario_manifest,
        args.primary_scenario,
    )


    scenario_names = (
        scenario_manifest[
            "scenario"
        ]
        .tolist()
    )


    scenario_parameters = {
        row[
            "scenario"
        ]: {
            "reciprocal_overlap":
                float(
                    row[
                        "reciprocal_overlap"
                    ]
                ),

            "local_gene_distance_bp":
                int(
                    row[
                        "local_gene_distance_bp"
                    ]
                ),
        }
        for _, row in scenario_manifest.iterrows()
    }


    # ========================================================
    # Prepare outputs
    # ========================================================

    for path in [
        args.out_global,
        args.out_species,
        args.out_transitions,
        args.out_cre,
        args.metadata_out,
    ]:

        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )


    # ========================================================
    # Load all sensitivity scenarios
    # ========================================================

    print()
    print("=" * 72)
    print("LOADING SENSITIVITY SCENARIOS")
    print("=" * 72)
    print()


    cache = {}


    for index, scenario in enumerate(
        scenario_names,
        start=1,
    ):

        print(
            f"[{index:02d}/{len(scenario_names):02d}] "
            f"{scenario}"
        )


        cache[
            scenario
        ] = load_scenario(
            args.sensitivity_dir,
            scenario,
            species,
            args.expected_reference_cres,
        )


    # ========================================================
    # Primary scenario
    # ========================================================

    primary = (
        cache[
            args.primary_scenario
        ]
        .rename(
            columns={
                "class":
                    "primary_class",
            }
        )
    )


    expected_cells = (
        len(
            species
        )
        * args.expected_reference_cres
    )


    if len(primary) != expected_cells:

        raise SystemExit(
            "ERROR: unexpected primary-scenario size.\n"
            f"Expected: {expected_cells}\n"
            f"Observed: {len(primary)}"
        )


    # ========================================================
    # Sensitivity comparisons
    # ========================================================

    summary_rows = []
    transition_rows = []
    species_rows = []
    cre_rows = []


    for scenario in scenario_names:

        parameters = scenario_parameters[
            scenario
        ]

        overlap = parameters[
            "reciprocal_overlap"
        ]

        distance = parameters[
            "local_gene_distance_bp"
        ]


        test = cache[
            scenario
        ]


        # ----------------------------------------------------
        # Match scenario to primary CRE x species universe
        # ----------------------------------------------------

        x = primary.merge(
            test,
            on=[
                "dmel_cre_id",
                "species",
            ],
            how="outer",
            validate="one_to_one",
            indicator=True,
        )


        if not x[
            "_merge"
        ].eq(
            "both"
        ).all():

            missing_primary = int(
                (
                    x[
                        "_merge"
                    ]
                    == "right_only"
                ).sum()
            )

            missing_test = int(
                (
                    x[
                        "_merge"
                    ]
                    == "left_only"
                ).sum()
            )


            raise SystemExit(
                "ERROR: CRE x species row mismatch for "
                f"scenario {scenario}.\n"
                f"Missing from primary:  {missing_primary}\n"
                f"Missing from scenario: {missing_test}"
            )


        x = x.drop(
            columns=[
                "_merge",
            ]
        )


        if len(x) != expected_cells:

            raise SystemExit(
                "ERROR: unexpected comparison size for "
                f"scenario {scenario}.\n"
                f"Expected: {expected_cells}\n"
                f"Observed: {len(x)}"
            )


        # ----------------------------------------------------
        # Global stability
        # ----------------------------------------------------

        same = (
            x[
                "primary_class"
            ]
            == x[
                "class"
            ]
        )


        n_same = int(
            same.sum()
        )

        n_changed = int(
            (
                ~same
            ).sum()
        )


        counts = (
            test[
                "class"
            ]
            .value_counts()
        )


        summary_rows.append({
            "scenario":
                scenario,

            "reciprocal_overlap":
                overlap,

            "local_gene_distance_bp":
                distance,

            "is_primary_scenario":
                (
                    "yes"
                    if scenario
                    == args.primary_scenario
                    else "no"
                ),

            "n_cells":
                len(
                    x
                ),

            "n_same_as_primary":
                n_same,

            "n_changed":
                n_changed,

            "percent_stable":
                100.0
                * n_same
                / len(
                    x
                ),

            "present":
                int(
                    counts.get(
                        "present",
                        0,
                    )
                ),

            "turnover_candidate":
                int(
                    counts.get(
                        "turnover_candidate",
                        0,
                    )
                ),

            "no_detected_CRE":
                int(
                    counts.get(
                        "no_detected_CRE",
                        0,
                    )
                ),

            "uncertain":
                int(
                    counts.get(
                        "uncertain",
                        0,
                    )
                ),
        })


        # ====================================================
        # State transitions
        # ====================================================

        changed = (
            x.loc[
                ~same
            ]
            .copy()
        )


        if not changed.empty:

            transitions = (
                changed
                .groupby(
                    [
                        "primary_class",
                        "class",
                    ],
                    sort=True,
                )
                .size()
                .reset_index(
                    name="n",
                )
            )


            transitions[
                "percent_of_changed"
            ] = (
                100.0
                * transitions[
                    "n"
                ]
                / n_changed
            )


            transitions.insert(
                0,
                "local_gene_distance_bp",
                distance,
            )

            transitions.insert(
                0,
                "reciprocal_overlap",
                overlap,
            )

            transitions.insert(
                0,
                "scenario",
                scenario,
            )


            transition_rows.append(
                transitions
            )


        # ====================================================
        # Species-level stability
        # ====================================================

        for sp in species:

            g = x.loc[
                x[
                    "species"
                ]
                == sp
            ]


            if len(g) != args.expected_reference_cres:

                raise SystemExit(
                    "ERROR: unexpected number of CREs for "
                    f"{scenario} / {sp}.\n"
                    f"Expected: "
                    f"{args.expected_reference_cres}\n"
                    f"Observed: {len(g)}"
                )


            stable = (
                g[
                    "primary_class"
                ]
                == g[
                    "class"
                ]
            )


            n_species_changed = int(
                (
                    ~stable
                ).sum()
            )


            species_rows.append({
                "scenario":
                    scenario,

                "reciprocal_overlap":
                    overlap,

                "local_gene_distance_bp":
                    distance,

                "species":
                    sp,

                "n_cells":
                    len(
                        g
                    ),

                "n_same_as_primary":
                    int(
                        stable.sum()
                    ),

                "n_changed":
                    n_species_changed,

                "percent_stable":
                    100.0
                    * stable.mean(),
            })


        # ====================================================
        # CRE-level stability
        # ====================================================

        for cre, g in x.groupby(
            "dmel_cre_id",
            sort=True,
        ):

            if len(g) != len(species):

                raise SystemExit(
                    "ERROR: unexpected species count for CRE "
                    f"{cre} in scenario {scenario}.\n"
                    f"Expected: {len(species)}\n"
                    f"Observed: {len(g)}"
                )


            stable = (
                g[
                    "primary_class"
                ]
                == g[
                    "class"
                ]
            )


            n_cre_changed = int(
                (
                    ~stable
                ).sum()
            )


            cre_rows.append({
                "scenario":
                    scenario,

                "reciprocal_overlap":
                    overlap,

                "local_gene_distance_bp":
                    distance,

                "dmel_cre_id":
                    cre,

                "n_species":
                    len(
                        g
                    ),

                "n_same_as_primary":
                    int(
                        stable.sum()
                    ),

                "n_changed":
                    n_cre_changed,

                "percent_stable":
                    100.0
                    * stable.mean(),
            })


    # ========================================================
    # Build output tables
    # ========================================================

    summary = pd.DataFrame(
        summary_rows
    )


    species_summary = pd.DataFrame(
        species_rows
    )


    cre_summary = pd.DataFrame(
        cre_rows
    )


    if transition_rows:

        transitions = pd.concat(
            transition_rows,
            ignore_index=True,
        )

    else:

        transitions = pd.DataFrame(
            columns=[
                "scenario",
                "reciprocal_overlap",
                "local_gene_distance_bp",
                "primary_class",
                "class",
                "n",
                "percent_of_changed",
            ]
        )


    # ========================================================
    # Final output QC
    # ========================================================

    if len(summary) != len(scenario_names):

        raise SystemExit(
            "ERROR: global sensitivity summary contains an "
            "unexpected number of scenarios."
        )


    expected_species_rows = (
        len(
            scenario_names
        )
        * len(
            species
        )
    )


    if len(species_summary) != expected_species_rows:

        raise SystemExit(
            "ERROR: unexpected number of species-stability "
            "rows.\n"
            f"Expected: {expected_species_rows}\n"
            f"Observed: {len(species_summary)}"
        )


    expected_cre_rows = (
        len(
            scenario_names
        )
        * args.expected_reference_cres
    )


    if len(cre_summary) != expected_cre_rows:

        raise SystemExit(
            "ERROR: unexpected number of CRE-stability rows.\n"
            f"Expected: {expected_cre_rows}\n"
            f"Observed: {len(cre_summary)}"
        )


    # --------------------------------------------------------
    # Primary scenario must be completely stable by definition
    # --------------------------------------------------------

    primary_summary = summary.loc[
        summary[
            "scenario"
        ]
        == args.primary_scenario
    ]


    if len(primary_summary) != 1:

        raise SystemExit(
            "ERROR: primary scenario missing or duplicated in "
            "global sensitivity summary."
        )


    if int(
        primary_summary.iloc[
            0
        ][
            "n_changed"
        ]
    ) != 0:

        raise SystemExit(
            "ERROR: primary scenario differs from itself."
        )


    # ========================================================
    # Deterministic ordering
    # ========================================================

    scenario_order = {
        scenario:
            index
        for index, scenario in enumerate(
            scenario_names
        )
    }


    species_order = {
        sp:
            index
        for index, sp in enumerate(
            species
        )
    }


    summary[
        "_scenario_rank"
    ] = summary[
        "scenario"
    ].map(
        scenario_order
    )


    summary = (
        summary
        .sort_values(
            [
                "_scenario_rank",
            ],
            kind="mergesort",
        )
        .drop(
            columns=[
                "_scenario_rank",
            ]
        )
        .reset_index(
            drop=True
        )
    )


    species_summary[
        "_scenario_rank"
    ] = species_summary[
        "scenario"
    ].map(
        scenario_order
    )

    species_summary[
        "_species_rank"
    ] = species_summary[
        "species"
    ].map(
        species_order
    )


    species_summary = (
        species_summary
        .sort_values(
            [
                "_scenario_rank",
                "_species_rank",
            ],
            kind="mergesort",
        )
        .drop(
            columns=[
                "_scenario_rank",
                "_species_rank",
            ]
        )
        .reset_index(
            drop=True
        )
    )


    cre_summary[
        "_scenario_rank"
    ] = cre_summary[
        "scenario"
    ].map(
        scenario_order
    )


    cre_summary = (
        cre_summary
        .sort_values(
            [
                "_scenario_rank",
                "dmel_cre_id",
            ],
            kind="mergesort",
        )
        .drop(
            columns=[
                "_scenario_rank",
            ]
        )
        .reset_index(
            drop=True
        )
    )


    if not transitions.empty:

        transitions[
            "_scenario_rank"
        ] = transitions[
            "scenario"
        ].map(
            scenario_order
        )


        transitions = (
            transitions
            .sort_values(
                [
                    "_scenario_rank",
                    "primary_class",
                    "class",
                ],
                kind="mergesort",
            )
            .drop(
                columns=[
                    "_scenario_rank",
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
        args.out_global,
        sep="\t",
        index=False,
    )


    species_summary.to_csv(
        args.out_species,
        sep="\t",
        index=False,
    )


    transitions.to_csv(
        args.out_transitions,
        sep="\t",
        index=False,
    )


    cre_summary.to_csv(
        args.out_cre,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Run metadata
    # ========================================================

    primary_parameters = (
        scenario_manifest.loc[
            scenario_manifest[
                "scenario"
            ]
            == args.primary_scenario
        ]
        .iloc[
            0
        ]
    )


    non_primary = summary.loc[
        summary[
            "scenario"
        ]
        != args.primary_scenario
    ]


    if non_primary.empty:

        minimum_stability = 100.0
        maximum_changed = 0

    else:

        minimum_stability = float(
            non_primary[
                "percent_stable"
            ].min()
        )

        maximum_changed = int(
            non_primary[
                "n_changed"
            ].max()
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

            "sensitivity_directory":
                str(
                    args.sensitivity_dir.resolve()
                ),

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

            "primary_scenario":
                args.primary_scenario,

            "primary_reciprocal_overlap":
                float(
                    primary_parameters[
                        "reciprocal_overlap"
                    ]
                ),

            "primary_local_gene_distance_bp":
                int(
                    primary_parameters[
                        "local_gene_distance_bp"
                    ]
                ),

            "expected_reference_cres":
                args.expected_reference_cres,

            "n_target_species":
                len(
                    species
                ),

            "n_scenarios":
                len(
                    scenario_names
                ),

            "n_cre_species_cells_per_scenario":
                expected_cells,

            "minimum_nonprimary_percent_stable":
                minimum_stability,

            "maximum_nonprimary_changed_cells":
                maximum_changed,

            "global_summary_output":
                str(
                    args.out_global.resolve()
                ),

            "species_stability_output":
                str(
                    args.out_species.resolve()
                ),

            "state_transitions_output":
                str(
                    args.out_transitions.resolve()
                ),

            "cre_stability_output":
                str(
                    args.out_cre.resolve()
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
    print("GLOBAL SENSITIVITY ANALYSIS")
    print("=" * 72)
    print()


    print(
        summary.to_string(
            index=False,
            float_format=lambda value:
                f"{value:.2f}",
        )
    )


    print()

    print(
        f"Primary scenario:      "
        f"{args.primary_scenario}"
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
        f"Sensitivity scenarios: "
        f"{len(scenario_names)}"
    )

    print(
        f"Cells per scenario:    "
        f"{expected_cells}"
    )

    print(
        f"Minimum stability:     "
        f"{minimum_stability:.2f}%"
    )


    print()

    print(
        f"Wrote global summary:\n"
        f"{args.out_global}"
    )

    print()

    print(
        f"Wrote species stability:\n"
        f"{args.out_species}"
    )

    print()

    print(
        f"Wrote state transitions:\n"
        f"{args.out_transitions}"
    )

    print()

    print(
        f"Wrote CRE stability:\n"
        f"{args.out_cre}"
    )

    print()

    print(
        f"Wrote metadata:\n"
        f"{args.metadata_out}"
    )

    print()

    print(
        "PASS: global sensitivity analysis completed."
    )


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    main()
