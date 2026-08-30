#!/usr/bin/env python3

# ============================================================
# Analyze secondary singleton CRE-state contrasts
#
# Purpose:
#   Identify CREs for which exactly one species within a focal
#   clade differs from a common state shared by all remaining
#   species.
#
# Categories:
#   - tier1: present <-> turnover_candidate
#   - tier2: positive CRE state <-> no_detected_CRE
#   - tier3: other singleton contrasts
#   - other_pattern: no valid singleton contrast
#
# Secondary Tier-1 candidates are additionally summarized
# across clades to identify recurrent CREs.
#
# Input/output paths and clade definitions:
#   Supplied by the pipeline wrapper using
#   config/candidate_config.sh.
# ============================================================

import argparse
from datetime import datetime
from pathlib import Path
import platform
import sys

import pandas as pd


# ============================================================
# CRE-state definitions
# ============================================================

REFERENCE_CRE_ID_COLUMN = "dmel_cre_id"

VALID_STATES = {
    "present",
    "turnover_candidate",
    "no_detected_CRE",
    "uncertain",
}

POSITIVE_STATES = {
    "present",
    "turnover_candidate",
}

OUTPUT_TIERS = [
    "tier1",
    "tier2",
    "tier3",
    "other_pattern",
]


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Identify singleton CRE-state contrasts within "
            "predefined phylogenetic clades and summarize "
            "secondary Tier-1 candidates."
        )
    )

    parser.add_argument(
        "--matrix",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--groups",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--out-long",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--out-summary",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--out-group-summary",
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
        default=None,
    )

    parser.add_argument(
        "--recurrence-min-clades",
        type=int,
        default=2,
    )

    return parser.parse_args()


# ============================================================
# Helper functions
# ============================================================

def require_file(path):
    """
    Abort if a required input file is missing.
    """

    if not path.is_file():
        raise SystemExit(
            f"ERROR: required input file not found:\n{path}"
        )


def join_unique(values):
    """
    Join unique non-empty values in deterministic order.
    """

    clean = sorted({
        str(value).strip()
        for value in values
        if str(value).strip()
        not in {
            "",
            "NA",
            "nan",
            "None",
        }
    })

    return "|".join(clean)


def load_groups(path):
    """
    Load clade definitions from focal_clades.tsv.

    The focal/comparison distinction is ignored here because
    every species can act as the discordant singleton.
    """

    groups = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
    ).fillna("")

    required = {
        "group_name",
        "focal_species",
        "comparison_species",
    }

    missing = (
        required
        - set(groups.columns)
    )

    if missing:
        raise SystemExit(
            "ERROR: clade definition file is missing columns:\n"
            + "\n".join(
                sorted(missing)
            )
        )

    for column in required:
        groups[column] = (
            groups[column]
            .astype(str)
            .str.strip()
        )

    duplicated_groups = (
        groups.loc[
            groups["group_name"].duplicated(
                keep=False
            ),
            "group_name",
        ]
        .unique()
        .tolist()
    )

    if duplicated_groups:
        raise SystemExit(
            "ERROR: duplicate clade names:\n"
            + "\n".join(
                sorted(duplicated_groups)
            )
        )

    result = {}

    for row in groups.itertuples(
        index=False
    ):

        comparisons = [
            value.strip()
            for value
            in row.comparison_species.split("|")
            if value.strip()
        ]

        species = [
            row.focal_species,
            *comparisons,
        ]

        if len(species) < 3:
            raise SystemExit(
                f"ERROR: {row.group_name} contains fewer "
                "than three species."
            )

        if len(species) != len(
            set(species)
        ):
            raise SystemExit(
                f"ERROR: duplicate species in {row.group_name}."
            )

        result[
            row.group_name
        ] = species

    return result


def classify_secondary_pattern(
    row,
    species,
):
    """
    Classify one CRE across one phylogenetic clade.

    A singleton contrast requires exactly one species to differ
    from one common state shared by all remaining species.
    """

    states = {
        sp: row[sp]
        for sp in species
    }

    unknown = (
        set(states.values())
        - VALID_STATES
    )

    if unknown:
        return {
            "pattern_class":
                "invalid_state",

            "tier":
                "other_pattern",

            "discordant_species":
                "NA",

            "discordant_state":
                "NA",

            "consensus_state":
                "NA",

            "n_consensus_species":
                0,
        }

    counts = (
        pd.Series(
            list(
                states.values()
            )
        )
        .value_counts()
    )

    # Exactly two states must occur.
    if len(counts) != 2:
        return {
            "pattern_class":
                "not_singleton_contrast",

            "tier":
                "other_pattern",

            "discordant_species":
                "NA",

            "discordant_state":
                "NA",

            "consensus_state":
                "NA",

            "n_consensus_species":
                0,
        }

    singleton_states = [
        state
        for state, n in counts.items()
        if n == 1
    ]

    # Exactly one state must occur once.
    if len(singleton_states) != 1:
        return {
            "pattern_class":
                "not_singleton_contrast",

            "tier":
                "other_pattern",

            "discordant_species":
                "NA",

            "discordant_state":
                "NA",

            "consensus_state":
                "NA",

            "n_consensus_species":
                0,
        }

    discordant_state = (
        singleton_states[0]
    )

    discordant_species = [
        sp
        for sp in species
        if states[sp]
        == discordant_state
    ]

    if len(discordant_species) != 1:
        return {
            "pattern_class":
                "not_singleton_contrast",

            "tier":
                "other_pattern",

            "discordant_species":
                "NA",

            "discordant_state":
                "NA",

            "consensus_state":
                "NA",

            "n_consensus_species":
                0,
        }

    discordant_species = (
        discordant_species[0]
    )

    consensus_species = [
        sp
        for sp in species
        if sp != discordant_species
    ]

    consensus_states = {
        states[sp]
        for sp in consensus_species
    }

    if len(consensus_states) != 1:
        return {
            "pattern_class":
                "not_singleton_contrast",

            "tier":
                "other_pattern",

            "discordant_species":
                "NA",

            "discordant_state":
                "NA",

            "consensus_state":
                "NA",

            "n_consensus_species":
                0,
        }

    consensus_state = next(
        iter(
            consensus_states
        )
    )

    n_consensus_species = len(
        consensus_species
    )

    # Tier 1:
    # present <-> turnover_candidate

    if (
        discordant_state in POSITIVE_STATES
        and
        consensus_state in POSITIVE_STATES
        and
        discordant_state != consensus_state
    ):
        return {
            "pattern_class":
                "present_vs_turnover",

            "tier":
                "tier1",

            "discordant_species":
                discordant_species,

            "discordant_state":
                discordant_state,

            "consensus_state":
                consensus_state,

            "n_consensus_species":
                n_consensus_species,
        }

    # Tier 2:
    # positive state <-> no_detected_CRE

    if (
        (
            discordant_state
            == "no_detected_CRE"
            and
            consensus_state
            in POSITIVE_STATES
        )
        or
        (
            consensus_state
            == "no_detected_CRE"
            and
            discordant_state
            in POSITIVE_STATES
        )
    ):
        return {
            "pattern_class":
                "detection_contrast",

            "tier":
                "tier2",

            "discordant_species":
                discordant_species,

            "discordant_state":
                discordant_state,

            "consensus_state":
                consensus_state,

            "n_consensus_species":
                n_consensus_species,
        }

    # Tier 3:
    # remaining singleton contrasts.

    return {
        "pattern_class":
            "other_singleton_contrast",

        "tier":
            "tier3",

        "discordant_species":
            discordant_species,

        "discordant_state":
            discordant_state,

        "consensus_state":
            consensus_state,

        "n_consensus_species":
            n_consensus_species,
    }


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    require_file(
        args.matrix
    )

    require_file(
        args.groups
    )

    if args.recurrence_min_clades < 1:
        raise SystemExit(
            "ERROR: recurrence-min-clades must be >= 1."
        )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    for path in [
        args.out_long,
        args.out_summary,
        args.out_group_summary,
        args.metadata_out,
    ]:
        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )


    # ========================================================
    # Load clade definitions
    # ========================================================

    groups = load_groups(
        args.groups
    )


    # ========================================================
    # Load CRE-state matrix
    # ========================================================

    df = pd.read_csv(
        args.matrix,
        sep="\t",
        dtype=str,
    ).fillna("")

    if REFERENCE_CRE_ID_COLUMN not in df.columns:
        raise SystemExit(
            "ERROR: CRE-state matrix lacks "
            f"'{REFERENCE_CRE_ID_COLUMN}'."
        )

    if df[
        REFERENCE_CRE_ID_COLUMN
    ].duplicated().any():

        duplicated = (
            df.loc[
                df[
                    REFERENCE_CRE_ID_COLUMN
                ].duplicated(
                    keep=False
                ),
                REFERENCE_CRE_ID_COLUMN,
            ]
            .drop_duplicates()
            .tolist()
        )

        raise SystemExit(
            "ERROR: duplicate D. melanogaster CRE IDs:\n"
            + "\n".join(
                sorted(duplicated)
            )
        )

    if (
        args.expected_reference_cres
        is not None
        and
        len(df)
        != args.expected_reference_cres
    ):
        raise SystemExit(
            "ERROR: unexpected number of reference CREs.\n"
            f"Expected: {args.expected_reference_cres}\n"
            f"Observed: {len(df)}"
        )


    # ========================================================
    # Validate species
    # ========================================================

    required_species = sorted({
        sp
        for species
        in groups.values()
        for sp
        in species
    })

    missing_species = sorted(
        set(
            required_species
        )
        - set(
            df.columns
        )
    )

    if missing_species:
        raise SystemExit(
            "ERROR: required species missing from "
            "CRE-state matrix:\n"
            + "\n".join(
                missing_species
            )
        )


    # ========================================================
    # Validate state vocabulary
    # ========================================================

    observed_states = {
        value
        for value
        in pd.unique(
            df[
                required_species
            ].values.ravel()
        )
        if value != ""
    }

    unknown_states = sorted(
        observed_states
        - VALID_STATES
    )

    if unknown_states:
        raise SystemExit(
            "ERROR: unexpected CRE-state values in matrix:\n"
            + "\n".join(
                unknown_states
            )
        )


    # ========================================================
    # Analyze each clade
    # ========================================================

    group_summary_rows = []
    secondary_tier1_rows = []


    for (
        group_name,
        species,
    ) in groups.items():

        sub = df[
            [
                REFERENCE_CRE_ID_COLUMN,
                *species,
            ]
        ].copy()

        classified = sub.apply(
            lambda row:
            classify_secondary_pattern(
                row,
                species,
            ),
            axis=1,
            result_type="expand",
        )

        sub = pd.concat(
            [
                sub,
                classified,
            ],
            axis=1,
        )

        sub["group_name"] = (
            group_name
        )

        sub["group_species"] = (
            "|".join(
                species
            )
        )

        sub["n_group_species"] = (
            len(
                species
            )
        )

        sub = (
            sub
            .sort_values(
                [
                    "tier",
                    REFERENCE_CRE_ID_COLUMN,
                ],
                kind="mergesort",
            )
            .reset_index(
                drop=True
            )
        )


        # ----------------------------------------------------
        # Write complete and tier-specific tables
        # ----------------------------------------------------

        all_file = (
            args.outdir
            / f"{group_name}_all_patterns.tsv"
        )

        sub.to_csv(
            all_file,
            sep="\t",
            index=False,
        )

        for tier in OUTPUT_TIERS:

            tier_file = (
                args.outdir
                / f"{group_name}_{tier}.tsv"
            )

            sub.loc[
                sub["tier"] == tier
            ].to_csv(
                tier_file,
                sep="\t",
                index=False,
            )


        # ----------------------------------------------------
        # Collect secondary Tier-1 candidates
        # ----------------------------------------------------

        tier1 = sub.loc[
            sub["tier"] == "tier1"
        ].copy()

        for _, row in tier1.iterrows():

            secondary_tier1_rows.append({
                REFERENCE_CRE_ID_COLUMN:
                    row[
                        REFERENCE_CRE_ID_COLUMN
                    ],

                "group_name":
                    group_name,

                "discordant_species":
                    row[
                        "discordant_species"
                    ],

                "discordant_state":
                    row[
                        "discordant_state"
                    ],

                "consensus_state":
                    row[
                        "consensus_state"
                    ],

                "n_group_species":
                    len(
                        species
                    ),

                "n_consensus_species":
                    row[
                        "n_consensus_species"
                    ],

                "group_species":
                    "|".join(
                        species
                    ),

                "pattern_class":
                    row[
                        "pattern_class"
                    ],
            })


        # ----------------------------------------------------
        # Group summary
        # ----------------------------------------------------

        counts = (
            sub["tier"]
            .value_counts()
        )

        group_summary_rows.append({
            "group_name":
                group_name,

            "group_species":
                "|".join(
                    species
                ),

            "n_group_species":
                len(
                    species
                ),

            "n_reference_cres":
                len(
                    sub
                ),

            "n_tier1":
                int(
                    counts.get(
                        "tier1",
                        0,
                    )
                ),

            "n_tier2":
                int(
                    counts.get(
                        "tier2",
                        0,
                    )
                ),

            "n_tier3":
                int(
                    counts.get(
                        "tier3",
                        0,
                    )
                ),

            "n_other_pattern":
                int(
                    counts.get(
                        "other_pattern",
                        0,
                    )
                ),
        })


        print()
        print("=" * 72)
        print(group_name)
        print("=" * 72)

        print(
            "Species: "
            + ", ".join(
                species
            )
        )

        print()

        print(
            "Tier 1 "
            "(singleton present vs turnover): "
            f"{counts.get('tier1', 0)}"
        )

        print(
            "Tier 2 "
            "(singleton detection contrast): "
            f"{counts.get('tier2', 0)}"
        )

        print(
            "Tier 3 "
            "(other singleton contrast): "
            f"{counts.get('tier3', 0)}"
        )

        print(
            "Other/non-singleton patterns: "
            f"{counts.get('other_pattern', 0)}"
        )


    # ========================================================
    # Combined secondary Tier-1 table
    # ========================================================

    secondary_long = pd.DataFrame(
        secondary_tier1_rows
    )

    if not secondary_long.empty:

        secondary_long = (
            secondary_long
            .sort_values(
                [
                    REFERENCE_CRE_ID_COLUMN,
                    "group_name",
                    "discordant_species",
                ],
                kind="mergesort",
            )
            .reset_index(
                drop=True
            )
        )

        if secondary_long.duplicated(
            subset=[
                REFERENCE_CRE_ID_COLUMN,
                "group_name",
            ]
        ).any():

            raise SystemExit(
                "ERROR: duplicate CRE x clade records "
                "in secondary Tier-1 table."
            )

    else:

        secondary_long = pd.DataFrame(
            columns=[
                REFERENCE_CRE_ID_COLUMN,
                "group_name",
                "discordant_species",
                "discordant_state",
                "consensus_state",
                "n_group_species",
                "n_consensus_species",
                "group_species",
                "pattern_class",
            ]
        )

    secondary_long.to_csv(
        args.out_long,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Recurrent secondary Tier-1 summary
    # ========================================================

    if not secondary_long.empty:

        secondary_summary = (
            secondary_long
            .groupby(
                REFERENCE_CRE_ID_COLUMN,
                sort=True,
            )
            .agg(
                n_secondary_tier1_clades=(
                    "group_name",
                    "nunique",
                ),

                secondary_tier1_clades=(
                    "group_name",
                    join_unique,
                ),

                discordant_species=(
                    "discordant_species",
                    join_unique,
                ),

                discordant_states=(
                    "discordant_state",
                    join_unique,
                ),

                consensus_states=(
                    "consensus_state",
                    join_unique,
                ),
            )
            .reset_index()
        )

        secondary_summary[
            "secondary_recurrence"
        ] = secondary_summary[
            "n_secondary_tier1_clades"
        ].apply(
            lambda n:
                "recurrent"
                if int(n)
                >= args.recurrence_min_clades
                else "single_clade"
        )

        secondary_summary = (
            secondary_summary
            .sort_values(
                [
                    "n_secondary_tier1_clades",
                    REFERENCE_CRE_ID_COLUMN,
                ],
                ascending=[
                    False,
                    True,
                ],
                kind="mergesort",
            )
            .reset_index(
                drop=True
            )
        )

    else:

        secondary_summary = pd.DataFrame(
            columns=[
                REFERENCE_CRE_ID_COLUMN,
                "n_secondary_tier1_clades",
                "secondary_tier1_clades",
                "discordant_species",
                "discordant_states",
                "consensus_states",
                "secondary_recurrence",
            ]
        )

    secondary_summary.to_csv(
        args.out_summary,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Group summary
    # ========================================================

    group_summary = (
        pd.DataFrame(
            group_summary_rows
        )
        .sort_values(
            "group_name",
            kind="mergesort",
        )
        .reset_index(
            drop=True
        )
    )

    group_summary.to_csv(
        args.out_group_summary,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Run metadata
    # ========================================================

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

            "matrix_input":
                str(
                    args.matrix.resolve()
                ),

            "clade_definition_input":
                str(
                    args.groups.resolve()
                ),

            "n_reference_cres":
                len(
                    df
                ),

            "n_clades":
                len(
                    groups
                ),

            "recurrence_min_clades":
                args.recurrence_min_clades,

            "n_secondary_tier1_rows":
                len(
                    secondary_long
                ),

            "n_unique_secondary_tier1_cres":
                (
                    secondary_long[
                        REFERENCE_CRE_ID_COLUMN
                    ].nunique()
                    if not secondary_long.empty
                    else 0
                ),

            "n_recurrent_secondary_cres":
                (
                    (
                        secondary_summary[
                            "secondary_recurrence"
                        ]
                        == "recurrent"
                    ).sum()
                    if not secondary_summary.empty
                    else 0
                ),
        }
    ])

    metadata.to_csv(
        args.metadata_out,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Final report
    # ========================================================

    print()
    print("=" * 72)
    print(
        "Secondary Tier-1 analysis complete"
    )
    print("=" * 72)

    print(
        f"Reference CREs: "
        f"{len(df)}"
    )

    print(
        f"Clades analyzed: "
        f"{len(groups)}"
    )

    print(
        f"Secondary Tier-1 group x CRE rows: "
        f"{len(secondary_long)}"
    )

    print(
        "Unique Secondary Tier-1 CREs: "
        f"{secondary_long[REFERENCE_CRE_ID_COLUMN].nunique() "
        f"if not secondary_long.empty else 0}"
    )

    recurrent = (
        secondary_summary.loc[
            secondary_summary[
                "secondary_recurrence"
            ] == "recurrent"
        ]
        if not secondary_summary.empty
        else secondary_summary
    )

    print(
        f"Recurrent Secondary Tier-1 CREs: "
        f"{len(recurrent)}"
    )

    print()
    print(
        f"Wrote per-clade tables to:\n"
        f"{args.outdir}"
    )

    print(
        f"Wrote combined Tier-1 table:\n"
        f"{args.out_long}"
    )

    print(
        f"Wrote recurrent CRE summary:\n"
        f"{args.out_summary}"
    )

    print(
        f"Wrote clade summary:\n"
        f"{args.out_group_summary}"
    )

    print(
        f"Wrote run metadata:\n"
        f"{args.metadata_out}"
    )


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    main()
