#!/usr/bin/env python3

# ============================================================
# Build gene-level Tier-1 candidate summary
#
# Purpose:
#   Aggregate CRE x FBgn candidate evidence into one row per
#   unique D. melanogaster gene.
#
# Summary:
#   Gene-level support is described by:
#   - number of associated Tier-1 CREs
#   - primary and secondary candidate-gene assignments
#   - original SCRMshaw target-gene support
#   - recurrent and focal CRE support
#   - strongest CRE-level candidate priority
#   - Tier-1 clade support
#   - carried-forward gene-assignment QC
#
# Notes:
#   Gene-assignment QC flags are retained as descriptive
#   evidence and are not used as automatic exclusion criteria.
#
#   No weighted gene-level score is calculated. Candidate genes
#   are instead ordered by transparent evidence variables.
#
# Input/output paths:
#   Supplied by the pipeline wrapper using
#   config/candidate_config.sh.
# ============================================================

import argparse
from datetime import datetime
import hashlib
from pathlib import Path
import platform
import sys

import pandas as pd


# ============================================================
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a gene-level summary from the exploded Tier-1 CRE x FBgn table."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--metadata-out", type=Path, required=True)
    return parser.parse_args()


# ============================================================
# Required columns
# ============================================================

REQUIRED_COLUMNS = {
    "dmel_cre_id",
    "fbgn",

    "gene_role",
    "assignment_evidence",

    "candidate_priority",
    "downstream_priority",

    "n_total_tier1_clades",
    "n_focal_tier1_clades",
    "n_secondary_only_clades",

    "tier1_clades",

    "is_reference_target",
    "is_primary_candidate",
    "is_secondary_candidate",

    "gene_assignment_qc",
}


# ============================================================
# Priority definitions
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


RECURRENT_CANDIDATE_PRIORITIES = {
    "focal_recurrent",
    "focal_plus_secondary_recurrent",
    "secondary_recurrent",
}


EXPECTED_DOWNSTREAM_PRIORITY = {
    "focal_recurrent":
        "high",

    "focal_plus_secondary_recurrent":
        "high",

    "secondary_recurrent":
        "medium",

    "focal_single":
        "medium",

    "secondary_single":
        "exploratory",
}


VALID_ASSIGNMENT_EVIDENCE = {
    "primary_and_reference_target",
    "secondary_and_reference_target",
    "primary_not_reference_target",
    "secondary_not_reference_target",
    "reference_target_only",
}


# ============================================================
# Helpers
# ============================================================

def require_file(path):
    """
    Abort if a required input file is missing.
    """

    if not path.is_file():

        raise SystemExit(
            "ERROR: required input file not found:\n"
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


def normalize_missing(value):
    """
    Convert missing or empty textual values to 'NA'.
    """

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

    clean = set()

    for value in values:

        value = normalize_missing(
            value
        )

        if value != "NA":
            clean.add(
                value
            )

    return "|".join(
        sorted(
            clean
        )
    )


def split_delimited_values(
    values,
    delimiter,
):
    """
    Collect unique values from delimited cells.
    """

    result = set()

    for value in values:

        value = normalize_missing(
            value
        )

        if value == "NA":
            continue

        for item in value.split(
            delimiter
        ):

            item = item.strip()

            if (
                item
                and item.lower()
                not in {
                    "na",
                    "nan",
                    "none",
                }
            ):
                result.add(
                    item
                )

    return sorted(
        result
    )


def file_sha256(path):
    """
    Calculate SHA256 checksum for reproducibility metadata.
    """

    sha = hashlib.sha256()

    with open(
        path,
        "rb",
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


def best_priority(
    values,
    ranking,
):
    """
    Return strongest priority according to predefined ranking.
    """

    clean = [
        normalize_missing(
            value
        )
        for value in values
        if normalize_missing(
            value
        )
        != "NA"
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
            ranking[
                value
            ],
    )


def best_gene_role(values):
    """
    Return strongest observed gene role.
    """

    clean = [
        normalize_missing(
            value
        )
        for value in values
        if normalize_missing(
            value
        )
        != "NA"
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
            GENE_ROLE_ORDER[
                value
            ],
    )


def validate_vocabulary(
    df,
    column,
    allowed,
):
    """
    Validate categorical values against an allowed vocabulary.
    """

    observed = {
        normalize_missing(
            value
        )
        for value in df[
            column
        ]
    }

    invalid = sorted(
        observed
        - set(
            allowed
        )
    )

    if invalid:

        raise SystemExit(
            f"ERROR: unexpected values in {column}:\n"
            + "\n".join(
                invalid
            )
        )


def validate_cre_level_consistency(
    df,
    columns,
):
    """
    Ensure CRE-level attributes are identical across all
    exploded FBgn rows belonging to the same CRE.
    """

    for column in columns:

        counts = (
            df.groupby(
                "dmel_cre_id",
                sort=False,
            )[
                column
            ]
            .nunique(
                dropna=False
            )
        )

        inconsistent = counts.loc[
            counts
            > 1
        ]

        if not inconsistent.empty:

            bad_ids = (
                inconsistent.index
                .astype(str)
                .tolist()
            )

            raise SystemExit(
                "ERROR: inconsistent CRE-level values for "
                f"column '{column}' in:\n"
                + "\n".join(
                    sorted(
                        bad_ids
                    )
                )
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

    args.metadata_out.parent.mkdir(
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

    if df.empty:

        raise SystemExit(
            "ERROR: exploded CRE x FBgn table contains no rows."
        )


    # ========================================================
    # Normalize key categorical columns
    # ========================================================

    categorical_columns = [
        "dmel_cre_id",
        "fbgn",
        "gene_role",
        "assignment_evidence",
        "candidate_priority",
        "downstream_priority",
        "is_reference_target",
        "is_primary_candidate",
        "is_secondary_candidate",
        "gene_assignment_qc",
    ]

    for column in categorical_columns:

        df[
            column
        ] = df[
            column
        ].apply(
            normalize_missing
        )


    # ========================================================
    # Basic identifier QC
    # ========================================================

    for column in [
        "dmel_cre_id",
        "fbgn",
    ]:

        invalid = (
            df[
                column
            ]
            == "NA"
        )

        if invalid.any():

            raise SystemExit(
                f"ERROR: missing values detected in {column}."
            )


    # --------------------------------------------------------
    # Each CRE x FBgn pair must be unique
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


    # ========================================================
    # Validate categorical vocabularies
    # ========================================================

    validate_vocabulary(
        df,
        "candidate_priority",
        CANDIDATE_PRIORITY_ORDER,
    )

    validate_vocabulary(
        df,
        "downstream_priority",
        DOWNSTREAM_PRIORITY_ORDER,
    )

    validate_vocabulary(
        df,
        "gene_role",
        GENE_ROLE_ORDER,
    )

    validate_vocabulary(
        df,
        "assignment_evidence",
        VALID_ASSIGNMENT_EVIDENCE,
    )

    for column in [
        "is_reference_target",
        "is_primary_candidate",
        "is_secondary_candidate",
    ]:

        validate_vocabulary(
            df,
            column,
            {
                "yes",
                "no",
            },
        )


    # ========================================================
    # Validate candidate/downstream priority consistency
    # ========================================================

    expected_downstream = df[
        "candidate_priority"
    ].map(
        EXPECTED_DOWNSTREAM_PRIORITY
    )

    inconsistent_priority = (
        expected_downstream
        != df[
            "downstream_priority"
        ]
    )

    if inconsistent_priority.any():

        bad = df.loc[
            inconsistent_priority,
            [
                "dmel_cre_id",
                "fbgn",
                "candidate_priority",
                "downstream_priority",
            ],
        ].drop_duplicates()

        raise SystemExit(
            "ERROR: candidate_priority and downstream_priority "
            "are inconsistent:\n\n"
            + bad.to_string(
                index=False
            )
        )


    # ========================================================
    # Validate gene roles against role flags
    # ========================================================

    both_primary_secondary = (
        (
            df[
                "is_primary_candidate"
            ]
            == "yes"
        )
        &
        (
            df[
                "is_secondary_candidate"
            ]
            == "yes"
        )
    )

    if both_primary_secondary.any():

        bad = df.loc[
            both_primary_secondary,
            [
                "dmel_cre_id",
                "fbgn",
                "gene_role",
            ],
        ]

        raise SystemExit(
            "ERROR: CRE x FBgn rows marked as both primary "
            "and secondary candidates:\n\n"
            + bad.to_string(
                index=False
            )
        )


    invalid_primary_role = (
        (
            df[
                "gene_role"
            ]
            == "primary_candidate"
        )
        !=
        (
            df[
                "is_primary_candidate"
            ]
            == "yes"
        )
    )

    if invalid_primary_role.any():

        bad = df.loc[
            invalid_primary_role,
            [
                "dmel_cre_id",
                "fbgn",
                "gene_role",
                "is_primary_candidate",
            ],
        ]

        raise SystemExit(
            "ERROR: primary gene-role flags are inconsistent:\n\n"
            + bad.to_string(
                index=False
            )
        )


    invalid_secondary_role = (
        (
            df[
                "gene_role"
            ]
            == "secondary_candidate"
        )
        !=
        (
            df[
                "is_secondary_candidate"
            ]
            == "yes"
        )
    )

    if invalid_secondary_role.any():

        bad = df.loc[
            invalid_secondary_role,
            [
                "dmel_cre_id",
                "fbgn",
                "gene_role",
                "is_secondary_candidate",
            ],
        ]

        raise SystemExit(
            "ERROR: secondary gene-role flags are inconsistent:\n\n"
            + bad.to_string(
                index=False
            )
        )


    # ========================================================
    # Numeric CRE-support columns
    # ========================================================

    numeric_columns = [
        "n_total_tier1_clades",
        "n_focal_tier1_clades",
        "n_secondary_only_clades",
    ]

    for column in numeric_columns:

        missing = (
            df[
                column
            ]
            .astype(str)
            .str.strip()
            .isin({
                "",
                "NA",
                "na",
                "nan",
                "None",
            })
        )

        if missing.any():

            bad = (
                df.loc[
                    missing,
                    [
                        "dmel_cre_id",
                        "fbgn",
                        column,
                    ],
                ]
                .drop_duplicates()
            )

            raise SystemExit(
                f"ERROR: missing values in required numeric "
                f"column '{column}':\n\n"
                + bad.to_string(
                    index=False
                )
            )

        df[
            column
        ] = pd.to_numeric(
            df[
                column
            ],
            errors="raise",
        )

        non_integer = (
            df[
                column
            ]
            % 1
            != 0
        )

        if non_integer.any():

            raise SystemExit(
                "ERROR: non-integer values found in "
                f"{column}."
            )

        if (
            df[
                column
            ]
            < 0
        ).any():

            raise SystemExit(
                "ERROR: negative values found in "
                f"{column}."
            )

        df[
            column
        ] = df[
            column
        ].astype(
            int
        )


    # --------------------------------------------------------
    # Total support must equal focal + secondary-only support
    # --------------------------------------------------------

    inconsistent_counts = (
        df[
            "n_total_tier1_clades"
        ]
        !=
        (
            df[
                "n_focal_tier1_clades"
            ]
            +
            df[
                "n_secondary_only_clades"
            ]
        )
    )

    if inconsistent_counts.any():

        bad = (
            df.loc[
                inconsistent_counts,
                [
                    "dmel_cre_id",
                    "n_total_tier1_clades",
                    "n_focal_tier1_clades",
                    "n_secondary_only_clades",
                ],
            ]
            .drop_duplicates()
        )

        raise SystemExit(
            "ERROR: inconsistent Tier-1 clade counts:\n\n"
            + bad.to_string(
                index=False
            )
        )


    # ========================================================
    # Validate repeated CRE-level information
    # ========================================================

    validate_cre_level_consistency(
        df,
        [
            "candidate_priority",
            "downstream_priority",
            "n_total_tier1_clades",
            "n_focal_tier1_clades",
            "n_secondary_only_clades",
            "tier1_clades",
            "gene_assignment_qc",
        ],
    )


    # ========================================================
    # Derived CRE-level support flags
    # ========================================================

    # Recurrence is inherited from the explicit candidate
    # priority assigned in Step 04 rather than redefined here.
    df[
        "is_recurrent_cre"
    ] = df[
        "candidate_priority"
    ].isin(
        RECURRENT_CANDIDATE_PRIORITIES
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
    # ========================================================

    gene_rows = []

    for fbgn, sub in (
        df.groupby(
            "fbgn",
            sort=True,
        )
    ):

        # ----------------------------------------------------
        # CRE sets associated with this gene
        # ----------------------------------------------------

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


        reference_target_cre_ids = sorted(
            sub.loc[
                sub[
                    "is_reference_target"
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
        # Gene-assignment QC support
        # ----------------------------------------------------

        qc_pass_cre_ids = sorted(
            sub.loc[
                sub[
                    "gene_assignment_qc"
                ]
                == "PASS",
                "dmel_cre_id",
            ].unique()
        )


        qc_flagged_cre_ids = sorted(
            sub.loc[
                sub[
                    "gene_assignment_qc"
                ]
                != "PASS",
                "dmel_cre_id",
            ].unique()
        )


        qc_values = join_unique(
            sub[
                "gene_assignment_qc"
            ]
        )


        qc_flags = split_delimited_values(
            sub.loc[
                sub[
                    "gene_assignment_qc"
                ]
                != "PASS",
                "gene_assignment_qc",
            ],
            ";",
        )


        # ----------------------------------------------------
        # Aggregate clade support
        # ----------------------------------------------------

        all_clades = split_delimited_values(
            sub[
                "tier1_clades"
            ],
            "|",
        )


        # ----------------------------------------------------
        # Aggregate evidence classes
        # ----------------------------------------------------

        assignment_evidence_values = join_unique(
            sub[
                "assignment_evidence"
            ]
        )


        # ----------------------------------------------------
        # Counts
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

        n_reference_target_cres = len(
            reference_target_cre_ids
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

        n_qc_pass_cres = len(
            qc_pass_cre_ids
        )

        n_qc_flagged_cres = len(
            qc_flagged_cre_ids
        )


        # ----------------------------------------------------
        # Maximum support across CREs
        # ----------------------------------------------------

        max_total_tier1_clades = int(
            sub[
                "n_total_tier1_clades"
            ].max()
        )


        max_focal_tier1_clades = int(
            sub[
                "n_focal_tier1_clades"
            ].max()
        )


        max_secondary_only_clades = int(
            sub[
                "n_secondary_only_clades"
            ].max()
        )


        # ----------------------------------------------------
        # Gene-level support flags
        # ----------------------------------------------------

        has_primary_assignment = (
            n_primary_cres
            > 0
        )

        has_secondary_assignment = (
            n_secondary_cres
            > 0
        )

        has_reference_target_support = (
            n_reference_target_cres
            > 0
        )


        # ----------------------------------------------------
        # Summary row
        # ----------------------------------------------------

        gene_rows.append({

            "fbgn":
                fbgn,

            # --------------------------------------------
            # Strongest observed evidence
            # --------------------------------------------

            "best_downstream_priority":
                best_priority(
                    sub[
                        "downstream_priority"
                    ],
                    DOWNSTREAM_PRIORITY_ORDER,
                ),

            "best_candidate_priority":
                best_priority(
                    sub[
                        "candidate_priority"
                    ],
                    CANDIDATE_PRIORITY_ORDER,
                ),

            "best_gene_role":
                best_gene_role(
                    sub[
                        "gene_role"
                    ]
                ),

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
            # Primary / secondary support
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

            # --------------------------------------------
            # Original SCRMshaw target support
            # --------------------------------------------

            "n_reference_target_cres":
                n_reference_target_cres,

            "reference_target_cre_ids":
                "|".join(
                    reference_target_cre_ids
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
            # Downstream-priority support
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

            # --------------------------------------------
            # QC support
            # --------------------------------------------

            "n_qc_pass_cres":
                n_qc_pass_cres,

            "qc_pass_cre_ids":
                "|".join(
                    qc_pass_cre_ids
                ),

            "n_qc_flagged_cres":
                n_qc_flagged_cres,

            "qc_flagged_cre_ids":
                "|".join(
                    qc_flagged_cre_ids
                ),

            "gene_assignment_qc_values":
                qc_values,

            "gene_assignment_qc_flags":
                ";".join(
                    qc_flags
                ),

            # --------------------------------------------
            # Assignment evidence
            # --------------------------------------------

            "assignment_evidence":
                assignment_evidence_values,

            # --------------------------------------------
            # Gene-level support flags
            # --------------------------------------------

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

            "has_reference_target_support":
                (
                    "yes"
                    if has_reference_target_support
                    else "no"
                ),
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
    # Stable deterministic ordering
    #
    # No weighted support score is used. Genes are ordered by
    # explicit, interpretable evidence variables.
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
                "n_high_priority_cres",
                "n_recurrent_cres",
                "n_focal_supported_cres",
                "n_primary_cres",
                "n_candidate_cres",
                "n_unique_tier1_clades",
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
    # Explicit final column order
    # ========================================================

    front = [
        "fbgn",

        "best_downstream_priority",
        "best_candidate_priority",
        "best_gene_role",

        "n_candidate_cres",
        "candidate_cre_ids",

        "n_primary_cres",
        "primary_cre_ids",

        "n_secondary_cres",
        "secondary_cre_ids",

        "n_reference_target_cres",
        "reference_target_cre_ids",

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

        "n_qc_pass_cres",
        "qc_pass_cre_ids",

        "n_qc_flagged_cres",
        "qc_flagged_cre_ids",

        "gene_assignment_qc_values",
        "gene_assignment_qc_flags",

        "assignment_evidence",

        "has_primary_assignment",
        "has_secondary_assignment",
        "has_reference_target_support",
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
    # Write gene-level summary
    # ========================================================

    out.to_csv(
        args.out,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Run metadata
    # ========================================================

    unique_input_cres = (
        df[
            "dmel_cre_id"
        ].nunique()
    )


    flagged_input_cres = (
        df.loc[
            df[
                "gene_assignment_qc"
            ]
            != "PASS",
            "dmel_cre_id",
        ]
        .nunique()
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

            "input_file":
                str(
                    args.input.resolve()
                ),

            "input_sha256":
                file_sha256(
                    args.input
                ),

            "output_file":
                str(
                    args.out.resolve()
                ),

            "output_sha256":
                file_sha256(
                    args.out
                ),

            "n_cre_fbgn_input_rows":
                len(
                    df
                ),

            "n_unique_candidate_cres":
                unique_input_cres,

            "n_unique_candidate_genes":
                len(
                    out
                ),

            "n_input_cres_with_qc_flags":
                int(
                    flagged_input_cres
                ),

            "n_genes_with_primary_support":
                int(
                    (
                        out[
                            "has_primary_assignment"
                        ]
                        == "yes"
                    ).sum()
                ),

            "n_genes_with_secondary_support":
                int(
                    (
                        out[
                            "has_secondary_assignment"
                        ]
                        == "yes"
                    ).sum()
                ),

            "n_genes_with_reference_target_support":
                int(
                    (
                        out[
                            "has_reference_target_support"
                        ]
                        == "yes"
                    ).sum()
                ),

            "n_genes_with_flagged_cre_evidence":
                int(
                    (
                        out[
                            "n_qc_flagged_cres"
                        ]
                        > 0
                    ).sum()
                ),

            "n_genes_with_high_priority_support":
                int(
                    (
                        out[
                            "n_high_priority_cres"
                        ]
                        > 0
                    ).sum()
                ),

            "n_genes_with_recurrent_support":
                int(
                    (
                        out[
                            "n_recurrent_cres"
                        ]
                        > 0
                    ).sum()
                ),

            "n_genes_with_focal_support":
                int(
                    (
                        out[
                            "n_focal_supported_cres"
                        ]
                        > 0
                    ).sum()
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
    print("Gene-level candidate summary complete")
    print("=" * 72)

    print(
        f"Input CRE x FBgn rows: "
        f"{len(df)}"
    )

    print(
        f"Unique candidate CREs: "
        f"{unique_input_cres}"
    )

    print(
        f"Unique candidate genes: "
        f"{len(out)}"
    )

    print(
        f"CREs with gene-assignment QC flags: "
        f"{flagged_input_cres}"
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


    # --------------------------------------------------------
    # Genes supported by multiple candidate CREs
    # --------------------------------------------------------

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
                    "best_downstream_priority",
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


    # --------------------------------------------------------
    # Top gene-level candidates
    # --------------------------------------------------------

    print()
    print(
        "Top gene-level candidates:"
    )

    print(
        out[
            [
                "fbgn",
                "best_downstream_priority",
                "best_candidate_priority",
                "best_gene_role",
                "n_candidate_cres",
                "n_primary_cres",
                "n_reference_target_cres",
                "n_recurrent_cres",
                "n_focal_supported_cres",
                "n_high_priority_cres",
                "n_unique_tier1_clades",
                "n_qc_flagged_cres",
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

    print()

    print(
        f"Wrote metadata:\n"
        f"{args.metadata_out}"
    )


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    main()
