#!/usr/bin/env python3

from pathlib import Path
import argparse

import pandas as pd


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = (
    Path.home()
    / "cre_turnover"
    / "project"
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

TABLE_DIR = (
    RESULTS_DIR
    / "tables"
)

TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

DEFAULT_INPUT = (
    TABLE_DIR
    / "tier1_candidate_fbgn_exploded.tsv"
)

DEFAULT_OUTPUT = (
    TABLE_DIR
    / "tier1_gene_level_summary.tsv"
)


# ============================================================
# Required columns
# ============================================================

REQUIRED_COLUMNS = {
    "dmel_cre_id",
    "fbgn",

    "gene_role",

    "candidate_priority",
    "downstream_priority",

    "n_total_tier1_clades",
    "n_focal_tier1_clades",
    "n_secondary_only_clades",

    "tier1_clades",

    "is_primary_candidate",
    "is_secondary_candidate",

    "gene_assignment_qc",
}


# ============================================================
# Priority rankings
# ============================================================

CANDIDATE_PRIORITY_ORDER = {
    "focal_recurrent": 1,
    "focal_plus_secondary_recurrent": 2,
    "secondary_recurrent": 3,
    "focal_single": 4,
    "secondary_single": 5,
}

DOWNSTREAM_PRIORITY_ORDER = {
    "high": 1,
    "medium": 2,
    "exploratory": 3,
}

GENE_ROLE_ORDER = {
    "primary_candidate": 1,
    "secondary_candidate": 2,
    "reference_target_only": 3,
}


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Build a gene-level summary of Tier-1 CRE candidates "
            "from the exploded CRE x FBgn table."
        )
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=(
            "Exploded CRE x FBgn candidate table "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=(
            "Gene-level candidate summary output "
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


def require_columns(
    df,
    required,
    label,
):

    missing = sorted(
        set(required)
        - set(df.columns)
    )

    if missing:
        raise SystemExit(
            f"ERROR: required columns missing from {label}:\n"
            + "\n".join(missing)
        )


def normalize_missing(value):

    if pd.isna(value):
        return "NA"

    value = str(value).strip()

    if value.lower() in {
        "",
        "na",
        "nan",
        "none",
    }:
        return "NA"

    return value


def join_unique(values):
    """
    Join unique non-missing values in deterministic order.
    """

    clean = sorted({
        normalize_missing(value)
        for value in values
        if normalize_missing(value) != "NA"
    })

    return "|".join(clean)


def split_pipe_values(values):
    """
    Collect unique values from one or more pipe-separated cells.
    """

    result = set()

    for value in values:

        value = normalize_missing(
            value
        )

        if value == "NA":
            continue

        for item in value.split("|"):

            item = item.strip()

            if item:
                result.add(
                    item
                )

    return sorted(
        result
    )


def best_priority(
    values,
    ranking,
):
    """
    Return strongest priority according to predefined ranking.
    """

    clean = [
        normalize_missing(value)
        for value in values
        if normalize_missing(value) != "NA"
    ]

    if not clean:
        return "NA"

    unknown = sorted({
        value
        for value in clean
        if value not in ranking
    })

    if unknown:
        raise SystemExit(
            "ERROR: unknown priority values encountered:\n"
            + "\n".join(
                unknown
            )
        )

    return min(
        clean,
        key=lambda value:
            ranking[value],
    )


def best_gene_role(values):
    """
    Return strongest observed gene role.
    """

    clean = [
        normalize_missing(value)
        for value in values
        if normalize_missing(value) != "NA"
    ]

    if not clean:
        return "NA"

    unknown = sorted({
        value
        for value in clean
        if value not in GENE_ROLE_ORDER
    })

    if unknown:
        raise SystemExit(
            "ERROR: unknown gene_role values encountered:\n"
            + "\n".join(
                unknown
            )
        )

    return min(
        clean,
        key=lambda value:
            GENE_ROLE_ORDER[value],
    )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    require_file(
        args.input
    )

    args.out.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


    # ========================================================
    # Load input
    # ========================================================

    df = pd.read_csv(
        args.input,
        sep="\t",
        dtype=str,
    ).fillna("")

    require_columns(
        df,
        REQUIRED_COLUMNS,
        "exploded CRE x FBgn table",
    )


    # --------------------------------------------------------
    # Validate CRE x FBgn uniqueness
    # --------------------------------------------------------

    if df.duplicated(
        subset=[
            "dmel_cre_id",
            "fbgn",
        ]
    ).any():

        duplicated = (
            df.loc[
                df.duplicated(
                    subset=[
                        "dmel_cre_id",
                        "fbgn",
                    ],
                    keep=False,
                ),
                [
                    "dmel_cre_id",
                    "fbgn",
                ],
            ]
            .drop_duplicates()
        )

        raise SystemExit(
            "ERROR: duplicate CRE x FBgn pairs in input:\n\n"
            + duplicated.to_string(
                index=False
            )
        )


    # --------------------------------------------------------
    # Validate gene assignment QC
    # --------------------------------------------------------

    non_pass_qc = (
        df[
            "gene_assignment_qc"
        ]
        != "PASS"
    )

    if non_pass_qc.any():

        bad = (
            df.loc[
                non_pass_qc,
                [
                    "dmel_cre_id",
                    "fbgn",
                    "gene_assignment_qc",
                ],
            ]
        )

        raise SystemExit(
            "ERROR: non-PASS gene assignments found.\n"
            "Resolve them before building gene-level summary:\n\n"
            + bad.to_string(
                index=False
            )
        )


    # ========================================================
    # Numeric columns
    # ============================================================

    for column in [
        "n_total_tier1_clades",
        "n_focal_tier1_clades",
        "n_secondary_only_clades",
    ]:

        df[
            column
        ] = pd.to_numeric(
            df[
                column
            ],
            errors="coerce",
        )


    # ========================================================
    # Derived CRE-level flags
    # ============================================================

    df[
        "is_recurrent_cre"
    ] = (
        df[
            "n_total_tier1_clades"
        ]
        >= 2
    )

    df[
        "has_focal_support"
    ] = (
        df[
            "n_focal_tier1_clades"
        ]
        >= 1
    )

    df[
        "is_high_priority_cre"
    ] = (
        df[
            "downstream_priority"
        ]
        == "high"
    )

    df[
        "is_medium_priority_cre"
    ] = (
        df[
            "downstream_priority"
        ]
        == "medium"
    )

    df[
        "is_exploratory_cre"
    ] = (
        df[
            "downstream_priority"
        ]
        == "exploratory"
    )


    # ========================================================
    # Gene-level aggregation
    # ============================================================

    gene_rows = []

    for fbgn, sub in (
        df.groupby(
            "fbgn",
            sort=True,
        )
    ):

        cre_ids = sorted(
            sub[
                "dmel_cre_id"
            ].unique()
        )

        primary_cre_ids = sorted(
            sub.loc[
                sub[
                    "is_primary_candidate"
                ]
                == "yes",
                "dmel_cre_id",
            ].unique()
        )

        secondary_cre_ids = sorted(
            sub.loc[
                sub[
                    "is_secondary_candidate"
                ]
                == "yes",
                "dmel_cre_id",
            ].unique()
        )

        recurrent_cre_ids = sorted(
            sub.loc[
                sub[
                    "is_recurrent_cre"
                ],
                "dmel_cre_id",
            ].unique()
        )

        focal_supported_cre_ids = sorted(
            sub.loc[
                sub[
                    "has_focal_support"
                ],
                "dmel_cre_id",
            ].unique()
        )

        high_priority_cre_ids = sorted(
            sub.loc[
                sub[
                    "is_high_priority_cre"
                ],
                "dmel_cre_id",
            ].unique()
        )

        medium_priority_cre_ids = sorted(
            sub.loc[
                sub[
                    "is_medium_priority_cre"
                ],
                "dmel_cre_id",
            ].unique()
        )

        exploratory_cre_ids = sorted(
            sub.loc[
                sub[
                    "is_exploratory_cre"
                ],
                "dmel_cre_id",
            ].unique()
        )


        # ----------------------------------------------------
        # Aggregate clades across all CREs assigned to gene
        # ----------------------------------------------------

        all_clades = split_pipe_values(
            sub[
                "tier1_clades"
            ]
        )


        # ----------------------------------------------------
        # Aggregate counts
        # ----------------------------------------------------

        n_candidate_cres = len(
            cre_ids
        )

        n_primary_cres = len(
            primary_cre_ids
        )

        n_secondary_cres = len(
            secondary_cre_ids
        )

        n_recurrent_cres = len(
            recurrent_cre_ids
        )

        n_focal_supported_cres = len(
            focal_supported_cre_ids
        )

        n_high_priority_cres = len(
            high_priority_cre_ids
        )

        n_medium_priority_cres = len(
            medium_priority_cre_ids
        )

        n_exploratory_cres = len(
            exploratory_cre_ids
        )


        # ----------------------------------------------------
        # Maximum support across CREs
        # ----------------------------------------------------

        max_total_tier1_clades = int(
            sub[
                "n_total_tier1_clades"
            ]
            .max()
        )

        max_focal_tier1_clades = int(
            sub[
                "n_focal_tier1_clades"
            ]
            .max()
        )

        max_secondary_only_clades = int(
            sub[
                "n_secondary_only_clades"
            ]
            .max()
        )


        # ----------------------------------------------------
        # Gene-level primary support
        # ----------------------------------------------------

        has_primary_assignment = (
            n_primary_cres > 0
        )

        has_secondary_assignment = (
            n_secondary_cres > 0
        )


        # ----------------------------------------------------
        # Summary row
        # ----------------------------------------------------

        gene_rows.append({
            "fbgn":
                fbgn,

            # --------------------------------------------
            # Overall CRE support
            # --------------------------------------------

            "n_candidate_cres":
                n_candidate_cres,

            "candidate_cre_ids":
                "|".join(
                    cre_ids
                ),

            # --------------------------------------------
            # Primary / secondary gene-role support
            # --------------------------------------------

            "n_primary_cres":
                n_primary_cres,

            "primary_cre_ids":
                "|".join(
                    primary_cre_ids
                ),

            "n_secondary_cres":
                n_secondary_cres,

            "secondary_cre_ids":
                "|".join(
                    secondary_cre_ids
                ),

            "best_gene_role":
                best_gene_role(
                    sub[
                        "gene_role"
                    ]
                ),

            "has_primary_assignment":
                (
                    "yes"
                    if has_primary_assignment
                    else "no"
                ),

            "has_secondary_assignment":
                (
                    "yes"
                    if has_secondary_assignment
                    else "no"
                ),

            # --------------------------------------------
            # Recurrence support
            # --------------------------------------------

            "n_recurrent_cres":
                n_recurrent_cres,

            "recurrent_cre_ids":
                "|".join(
                    recurrent_cre_ids
                ),

            # --------------------------------------------
            # Focal support
            # --------------------------------------------

            "n_focal_supported_cres":
                n_focal_supported_cres,

            "focal_supported_cre_ids":
                "|".join(
                    focal_supported_cre_ids
                ),

            # --------------------------------------------
            # Downstream priorities
            # --------------------------------------------

            "n_high_priority_cres":
                n_high_priority_cres,

            "high_priority_cre_ids":
                "|".join(
                    high_priority_cre_ids
                ),

            "n_medium_priority_cres":
                n_medium_priority_cres,

            "medium_priority_cre_ids":
                "|".join(
                    medium_priority_cre_ids
                ),

            "n_exploratory_cres":
                n_exploratory_cres,

            "exploratory_cre_ids":
                "|".join(
                    exploratory_cre_ids
                ),

            # --------------------------------------------
            # Best observed priority
            # --------------------------------------------

            "best_candidate_priority":
                best_priority(
                    sub[
                        "candidate_priority"
                    ],
                    CANDIDATE_PRIORITY_ORDER,
                ),

            "best_downstream_priority":
                best_priority(
                    sub[
                        "downstream_priority"
                    ],
                    DOWNSTREAM_PRIORITY_ORDER,
                ),

            # --------------------------------------------
            # Clade support
            # --------------------------------------------

            "n_unique_tier1_clades":
                len(
                    all_clades
                ),

            "tier1_clades":
                "|".join(
                    all_clades
                ),

            "max_total_tier1_clades_per_cre":
                max_total_tier1_clades,

            "max_focal_tier1_clades_per_cre":
                max_focal_tier1_clades,

            "max_secondary_only_clades_per_cre":
                max_secondary_only_clades,
        })


    # ========================================================
    # Build gene-level table
    # ========================================================

    out = pd.DataFrame(
        gene_rows
    )


    if out.empty:
        raise SystemExit(
            "ERROR: no gene-level rows generated."
        )


    # ========================================================
    # Gene-level prioritization support variables
    #
    # These are descriptive, not a new biological
    # classification.
    # ========================================================

    out[
        "gene_support_score"
    ] = (
        4
        * out[
            "n_high_priority_cres"
        ]
        +
        2
        * out[
            "n_medium_priority_cres"
        ]
        +
        out[
            "n_exploratory_cres"
        ]
        +
        2
        * out[
            "n_recurrent_cres"
        ]
        +
        out[
            "n_focal_supported_cres"
        ]
        +
        out[
            "n_primary_cres"
        ]
    )


    # ========================================================
    # Stable deterministic ordering
    # ========================================================

    out[
        "_best_downstream_rank"
    ] = (
        out[
            "best_downstream_priority"
        ]
        .map(
            DOWNSTREAM_PRIORITY_ORDER
        )
        .fillna(
            99
        )
    )


    out[
        "_best_candidate_rank"
    ] = (
        out[
            "best_candidate_priority"
        ]
        .map(
            CANDIDATE_PRIORITY_ORDER
        )
        .fillna(
            99
        )
    )


    out = (
        out
        .sort_values(
            [
                "_best_downstream_rank",
                "gene_support_score",
                "n_high_priority_cres",
                "n_recurrent_cres",
                "n_focal_supported_cres",
                "n_candidate_cres",
                "_best_candidate_rank",
                "fbgn",
            ],
            ascending=[
                True,
                False,
                False,
                False,
                False,
                False,
                True,
                True,
            ],
            kind="mergesort",
        )
        .drop(
            columns=[
                "_best_downstream_rank",
                "_best_candidate_rank",
            ]
        )
        .reset_index(
            drop=True
        )
    )


    # ========================================================
    # Final QC
    # ========================================================

    if out[
        "fbgn"
    ].duplicated().any():

        raise SystemExit(
            "ERROR: duplicate FBgn rows in gene-level summary."
        )


    # ========================================================
    # Useful column order
    # ========================================================

    front = [
        "fbgn",

        "gene_support_score",

        "best_downstream_priority",
        "best_candidate_priority",
        "best_gene_role",

        "n_candidate_cres",
        "candidate_cre_ids",

        "n_primary_cres",
        "primary_cre_ids",

        "n_secondary_cres",
        "secondary_cre_ids",

        "n_recurrent_cres",
        "recurrent_cre_ids",

        "n_focal_supported_cres",
        "focal_supported_cre_ids",

        "n_high_priority_cres",
        "high_priority_cre_ids",

        "n_medium_priority_cres",
        "medium_priority_cre_ids",

        "n_exploratory_cres",
        "exploratory_cre_ids",

        "n_unique_tier1_clades",
        "tier1_clades",

        "max_total_tier1_clades_per_cre",
        "max_focal_tier1_clades_per_cre",
        "max_secondary_only_clades_per_cre",

        "has_primary_assignment",
        "has_secondary_assignment",
    ]


    remaining = [
        column
        for column in out.columns
        if column not in front
    ]


    out = out[
        front
        + remaining
    ]


    # ========================================================
    # Write
    # ========================================================

    out.to_csv(
        args.out,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Console summary
    # ========================================================

    print()
    print("=" * 72)
    print("Gene-level candidate summary complete")
    print("=" * 72)

    print(
        f"Input CRE x FBgn rows: "
        f"{len(df)}"
    )

    print(
        f"Unique candidate CREs: "
        f"{df['dmel_cre_id'].nunique()}"
    )

    print(
        f"Unique candidate genes: "
        f"{len(out)}"
    )


    print()
    print(
        "Best downstream priority:"
    )

    print(
        out[
            "best_downstream_priority"
        ]
        .value_counts()
        .to_string()
    )


    print()
    print(
        "Genes supported by multiple candidate CREs:"
    )

    multi = (
        out.loc[
            out[
                "n_candidate_cres"
            ]
            > 1
        ]
    )

    if multi.empty:

        print(
            "None"
        )

    else:

        print(
            multi[
                [
                    "fbgn",
                    "gene_support_score",
                    "n_candidate_cres",
                    "n_primary_cres",
                    "n_recurrent_cres",
                    "n_focal_supported_cres",
                    "n_high_priority_cres",
                    "n_unique_tier1_clades",
                    "candidate_cre_ids",
                ]
            ].to_string(
                index=False
            )
        )


    print()
    print(
        "Top gene-level candidates:"
    )

    print(
        out[
            [
                "fbgn",
                "gene_support_score",
                "best_downstream_priority",
                "best_candidate_priority",
                "best_gene_role",
                "n_candidate_cres",
                "n_recurrent_cres",
                "n_focal_supported_cres",
                "n_high_priority_cres",
                "n_unique_tier1_clades",
            ]
        ]
        .head(
            20
        )
        .to_string(
            index=False
        )
    )


    print()
    print(
        f"Wrote gene-level summary:\n"
        f"{args.out}"
    )


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    main()
