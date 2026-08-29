#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime
import argparse
import platform
import sys

import pandas as pd


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = (
    Path.home()
    / "cre_turnover"
    / "project"
)

CLASS_DIR = (
    PROJECT_DIR
    / "cre_classification"
)

CANDIDATE_DIR = (
    PROJECT_DIR
    / "downstream_analyses"
    / "candidate_analysis"
)

RESULTS_DIR = (
    CANDIDATE_DIR
    / "results"
)

SECONDARY_DIR = (
    RESULTS_DIR
    / "secondary_clades"
)

SECONDARY_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


DEFAULT_MATRIX = (
    CLASS_DIR
    / "results"
    / "cre_turnover_matrix.tsv"
)

DEFAULT_OUT_LONG = (
    SECONDARY_DIR
    / "secondary_tier1_all_clades.tsv"
)

DEFAULT_OUT_SUMMARY = (
    SECONDARY_DIR
    / "secondary_tier1_recurrent_cres.tsv"
)

DEFAULT_OUT_GROUP_SUMMARY = (
    SECONDARY_DIR
    / "secondary_clade_summary.tsv"
)

DEFAULT_METADATA = (
    SECONDARY_DIR
    / "secondary_tier1_run_metadata.tsv"
)


# ============================================================
# Clade definitions
#
# These are exploratory phylogenetic groups.
# No species is predefined as focal here.
# Any one species may be the discordant singleton.
# ============================================================

GROUPS = {
    "rufa_group": [
        "druf",
        "dkik",
        "dbun",
        "dbir",
        "d_serrata",
    ],

    "immigrans_group": [
        "dimm",
        "dfor",
        "dsul",
        "dnas",
    ],

    "obscura_group": [
        "dsub",
        "dbif",
        "dobs",
    ],

    "azteca_affinis_miranda_group": [
        "dazt",
        "d_affinis",
        "dmir",
        "dper",
    ],

    "teissieri_group": [
        "dtei",
        "d_simulans",
        "d_lutescens",
        "dsuz",
    ],

    "repleta_group": [
        "d_repleta",
        "d_buzzatii",
        "d_mojavensis",
        "d_arizonae",
    ],
}


# ============================================================
# CRE-state definitions
# ============================================================

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
# Argument parsing
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Identify singleton discordant CRE-state patterns "
            "within predefined phylogenetic clades and summarize "
            "recurrent secondary Tier-1 CRE candidates."
        )
    )

    parser.add_argument(
        "--matrix",
        type=Path,
        default=DEFAULT_MATRIX,
        help=(
            "CRE turnover matrix "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=SECONDARY_DIR,
        help=(
            "Output directory for per-clade tables "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--out-long",
        type=Path,
        default=DEFAULT_OUT_LONG,
        help=(
            "Combined secondary Tier-1 long table "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--out-summary",
        type=Path,
        default=DEFAULT_OUT_SUMMARY,
        help=(
            "Recurrent secondary Tier-1 CRE summary "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--out-group-summary",
        type=Path,
        default=DEFAULT_OUT_GROUP_SUMMARY,
        help=(
            "Per-clade pattern-count summary "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--metadata-out",
        type=Path,
        default=DEFAULT_METADATA,
        help=(
            "Run metadata output "
            "(default: %(default)s)"
        ),
    )

    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_file(path):
    if not path.exists():
        raise SystemExit(
            "ERROR: required input file not found:\n"
            f"{path}"
        )


def join_unique(values):
    """
    Join unique, non-empty values in deterministic order.
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


def classify_secondary_pattern(
    row,
    species,
):
    """
    Classify one CRE across one phylogenetic clade.

    A secondary singleton pattern requires exactly one species
    to differ from all other species in the clade.

    Examples
    --------
    3 species:
        A != B = B

    4 species:
        A != B = B = B

    5 species:
        A != B = B = B = B

    Returns
    -------
    dict
        pattern_class
        tier
        discordant_species
        discordant_state
        consensus_state
        n_consensus_species
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

    # --------------------------------------------------------
    # Count state frequencies
    # --------------------------------------------------------

    counts = (
        pd.Series(
            list(
                states.values()
            )
        )
        .value_counts()
    )

    # --------------------------------------------------------
    # Singleton pattern requires exactly two states
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # Exactly one state must occur exactly once
    # --------------------------------------------------------

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

    if len(
        discordant_species
    ) != 1:

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
        if sp
        != discordant_species
    ]

    consensus_states = {
        states[sp]
        for sp
        in consensus_species
    }

    if (
        len(
            consensus_states
        )
        != 1
    ):

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

    # --------------------------------------------------------
    # Tier 1:
    # present <-> turnover_candidate
    # --------------------------------------------------------

    if (
        discordant_state
        in POSITIVE_STATES
        and
        consensus_state
        in POSITIVE_STATES
        and
        discordant_state
        != consensus_state
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

    # --------------------------------------------------------
    # Tier 2:
    # positive <-> no_detected_CRE
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # Tier 3:
    # singleton contrast exists, but involves uncertain
    # or another non-priority combination
    # --------------------------------------------------------

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

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.out_long.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.out_summary.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.out_group_summary.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.metadata_out.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ========================================================
    # Load CRE-state matrix
    # ========================================================

    df = pd.read_csv(
        args.matrix,
        sep="\t",
        dtype=str,
    ).fillna("")

    if (
        "dmel_cre_id"
        not in df.columns
    ):

        raise SystemExit(
            "ERROR: CRE-state matrix lacks "
            "dmel_cre_id."
        )

    # --------------------------------------------------------
    # Unique CRE IDs
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
            "ERROR: duplicate dmel_cre_id values "
            "in CRE-state matrix:\n"
            + "\n".join(
                sorted(
                    duplicated
                )
            )
        )

    # ========================================================
    # Validate required species
    # ========================================================

    required_species = sorted({
        sp
        for species
        in GROUPS.values()
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
            "ERROR: unexpected CRE-state values "
            "in matrix:\n"
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
    ) in GROUPS.items():

        # ----------------------------------------------------
        # Extract group matrix
        # ----------------------------------------------------

        sub = df[
            [
                "dmel_cre_id",
                *species,
            ]
        ].copy()

        # ----------------------------------------------------
        # Classify all CREs
        # ----------------------------------------------------

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

        sub[
            "group_name"
        ] = group_name

        sub[
            "group_species"
        ] = "|".join(
            species
        )

        sub[
            "n_group_species"
        ] = len(
            species
        )

        # ----------------------------------------------------
        # Stable ordering
        # ----------------------------------------------------

        sub = (
            sub
            .sort_values(
                [
                    "tier",
                    "dmel_cre_id",
                ],
                kind="mergesort",
            )
            .reset_index(
                drop=True
            )
        )

        # ----------------------------------------------------
        # Complete group table
        # ----------------------------------------------------

        all_file = (
            args.outdir
            / (
                f"{group_name}"
                "_all_patterns.tsv"
            )
        )

        sub.to_csv(
            all_file,
            sep="\t",
            index=False,
        )

        # ----------------------------------------------------
        # Tier-specific tables
        # ----------------------------------------------------

        for tier in OUTPUT_TIERS:

            tier_df = (
                sub.loc[
                    sub[
                        "tier"
                    ]
                    == tier
                ]
                .copy()
            )

            tier_file = (
                args.outdir
                / (
                    f"{group_name}_"
                    f"{tier}.tsv"
                )
            )

            tier_df.to_csv(
                tier_file,
                sep="\t",
                index=False,
            )

        # ----------------------------------------------------
        # Extract Secondary Tier-1 rows for combined summary
        # ----------------------------------------------------

        tier1 = (
            sub.loc[
                sub[
                    "tier"
                ]
                == "tier1"
            ]
            .copy()
        )

        if not tier1.empty:

            for _, row in (
                tier1.iterrows()
            ):

                secondary_tier1_rows.append({
                    "dmel_cre_id":
                        row[
                            "dmel_cre_id"
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
        # Group-level summary
        # ----------------------------------------------------

        counts = (
            sub[
                "tier"
            ]
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

        # ----------------------------------------------------
        # Console report
        # ----------------------------------------------------

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

        if not tier1.empty:

            print()
            print(
                "Secondary Tier-1 candidates:"
            )

            print(
                tier1[
                    [
                        "dmel_cre_id",
                        "discordant_species",
                        "discordant_state",
                        "consensus_state",
                        *species,
                    ]
                ].to_string(
                    index=False
                )
            )

    # ========================================================
    # Combined Secondary Tier-1 long table
    # ========================================================

    secondary_long = pd.DataFrame(
        secondary_tier1_rows
    )

    if not secondary_long.empty:

        secondary_long = (
            secondary_long
            .sort_values(
                [
                    "dmel_cre_id",
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
                "dmel_cre_id",
                "group_name",
            ]
        ).any():

            raise SystemExit(
                "ERROR: duplicate CRE x clade records "
                "in secondary Tier-1 table."
            )

    secondary_long.to_csv(
        args.out_long,
        sep="\t",
        index=False,
    )

    # ========================================================
    # Recurrent Secondary Tier-1 summary
    # ========================================================

    if not secondary_long.empty:

        secondary_summary = (
            secondary_long
            .groupby(
                "dmel_cre_id",
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
                if int(n) >= 2
                else "single_clade"
        )

        secondary_summary = (
            secondary_summary
            .sort_values(
                [
                    "n_secondary_tier1_clades",
                    "dmel_cre_id",
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
                "dmel_cre_id",
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

    group_summary = pd.DataFrame(
        group_summary_rows
    )

    group_summary = (
        group_summary
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

            "n_reference_cres":
                len(
                    df
                ),

            "n_clades":
                len(
                    GROUPS
                ),

            "n_secondary_tier1_rows":
                len(
                    secondary_long
                ),

            "n_unique_secondary_tier1_cres":
                (
                    secondary_long[
                        "dmel_cre_id"
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
    # Final console summary
    # ========================================================

    print()
    print("=" * 72)
    print("Secondary Tier analysis complete")
    print("=" * 72)

    print(
        f"Reference CREs: "
        f"{len(df)}"
    )

    print(
        f"Clades analyzed: "
        f"{len(GROUPS)}"
    )

    print(
        f"Secondary Tier-1 group x CRE rows: "
        f"{len(secondary_long)}"
    )

    print(
        "Unique Secondary Tier-1 CREs: "
        f"{secondary_long['dmel_cre_id'].nunique() if not secondary_long.empty else 0}"
    )

    if not secondary_summary.empty:

        recurrent = (
            secondary_summary.loc[
                secondary_summary[
                    "secondary_recurrence"
                ]
                == "recurrent"
            ]
        )

        print(
            f"Recurrent Secondary Tier-1 CREs: "
            f"{len(recurrent)}"
        )

        if not recurrent.empty:

            print()
            print(
                "Recurrent Secondary Tier-1 candidates:"
            )

            print(
                recurrent[
                    [
                        "dmel_cre_id",
                        "n_secondary_tier1_clades",
                        "secondary_tier1_clades",
                        "discordant_species",
                    ]
                ].to_string(
                    index=False
                )
            )

    print()
    print(
        f"Wrote per-clade tables to:\n"
        f"{args.outdir}"
    )

    print()

    print(
        f"Wrote combined Tier-1 table:\n"
        f"{args.out_long}"
    )

    print()

    print(
        f"Wrote recurrent CRE summary:\n"
        f"{args.out_summary}"
    )

    print()

    print(
        f"Wrote clade summary:\n"
        f"{args.out_group_summary}"
    )

    print()

    print(
        f"Wrote run metadata:\n"
        f"{args.metadata_out}"
    )


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    main()
