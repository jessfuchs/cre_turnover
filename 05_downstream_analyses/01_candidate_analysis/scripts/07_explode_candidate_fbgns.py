#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime
import argparse
import hashlib
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
    / "tier1_candidate_gene_assignments.tsv"
)

DEFAULT_OUTPUT = (
    TABLE_DIR
    / "tier1_candidate_fbgn_exploded.tsv"
)

# ============================================================
# Required columns
# ============================================================

REQUIRED_COLUMNS = {
    "dmel_cre_id",
    "candidate_priority",
    "downstream_priority",
    "fbgn_target_genes",

    "primary_candidate_fbgn",
    "primary_candidate_distance_bp",
    "primary_candidate_relation",

    "secondary_candidate_fbgn",
    "secondary_candidate_distance_bp",
    "secondary_candidate_relation",

    "gene_assignment_qc",
}


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Convert Tier-1 CRE candidate gene assignments "
            "into one row per CRE x FBgn for downstream "
            "gene-level candidate analysis."
        )
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=(
            "Tier-1 candidate gene-assignment table "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=(
            "Long-form CRE x FBgn output "
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


def parse_fbgns(value):
    """
    Parse pipe-separated FBgn identifiers.
    """

    value = normalize_missing(
        value
    )

    if value == "NA":
        return []

    return sorted({
        gene.strip()
        for gene in value.split("|")
        if gene.strip()
        and gene.strip().lower()
        not in {
            "",
            "na",
            "nan",
            "none",
        }
    })


def file_sha256(path):

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


def optional_value(
    row,
    column,
):
    """
    Safely copy an optional source column.
    """

    if column not in row.index:
        return "NA"

    return normalize_missing(
        row[
            column
        ]
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
        "candidate gene-assignment table",
    )


    # --------------------------------------------------------
    # One row per CRE
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
            "ERROR: duplicate CRE IDs in input:\n"
            + "\n".join(
                sorted(
                    duplicated
                )
            )
        )


    # ========================================================
    # Explode to one row per CRE x FBgn
    # ========================================================

    rows = []


    for _, row in (
        df.iterrows()
    ):

        cre_id = row[
            "dmel_cre_id"
        ]

        reference_genes = parse_fbgns(
            row[
                "fbgn_target_genes"
            ]
        )

        primary_fbgn = normalize_missing(
            row[
                "primary_candidate_fbgn"
            ]
        )

        secondary_fbgn = normalize_missing(
            row[
                "secondary_candidate_fbgn"
            ]
        )


        # ----------------------------------------------------
        # Gene universe for this CRE
        #
        # Usually primary/secondary are already members of
        # reference_genes, but unioning them explicitly makes
        # the transformation robust and auditable.
        # ----------------------------------------------------

        gene_set = set(
            reference_genes
        )

        if primary_fbgn != "NA":
            gene_set.add(
                primary_fbgn
            )

        if secondary_fbgn != "NA":
            gene_set.add(
                secondary_fbgn
            )


        # ----------------------------------------------------
        # No gene information at all
        # ----------------------------------------------------

        if not gene_set:

            raise SystemExit(
                "ERROR: candidate CRE has no usable FBgn "
                "assignment:\n"
                f"{cre_id}"
            )


        # ----------------------------------------------------
        # One row per unique FBgn
        # ----------------------------------------------------

        for fbgn in sorted(
            gene_set
        ):

            is_reference_target = (
                fbgn
                in reference_genes
            )

            is_primary = (
                fbgn
                == primary_fbgn
            )

            is_secondary = (
                fbgn
                == secondary_fbgn
            )


            # ------------------------------------------------
            # Assignment role
            # ------------------------------------------------

            if is_primary:
                gene_role = (
                    "primary_candidate"
                )

                candidate_distance = (
                    normalize_missing(
                        row[
                            "primary_candidate_distance_bp"
                        ]
                    )
                )

                candidate_relation = (
                    normalize_missing(
                        row[
                            "primary_candidate_relation"
                        ]
                    )
                )

            elif is_secondary:
                gene_role = (
                    "secondary_candidate"
                )

                candidate_distance = (
                    normalize_missing(
                        row[
                            "secondary_candidate_distance_bp"
                        ]
                    )
                )

                candidate_relation = (
                    normalize_missing(
                        row[
                            "secondary_candidate_relation"
                        ]
                    )
                )

            else:
                gene_role = (
                    "reference_target_only"
                )

                candidate_distance = (
                    "NA"
                )

                candidate_relation = (
                    "NA"
                )


            # ------------------------------------------------
            # Assignment evidence class
            # ------------------------------------------------

            if (
                is_primary
                and is_reference_target
            ):
                assignment_evidence = (
                    "primary_and_reference_target"
                )

            elif (
                is_secondary
                and is_reference_target
            ):
                assignment_evidence = (
                    "secondary_and_reference_target"
                )

            elif (
                is_primary
                and not is_reference_target
            ):
                assignment_evidence = (
                    "primary_not_reference_target"
                )

            elif (
                is_secondary
                and not is_reference_target
            ):
                assignment_evidence = (
                    "secondary_not_reference_target"
                )

            else:
                assignment_evidence = (
                    "reference_target_only"
                )


            # ------------------------------------------------
            # Long-form row
            # ------------------------------------------------

            rows.append({
                "dmel_cre_id":
                    cre_id,

                "fbgn":
                    fbgn,

                "gene_role":
                    gene_role,

                "assignment_evidence":
                    assignment_evidence,

                "is_reference_target":
                    (
                        "yes"
                        if is_reference_target
                        else "no"
                    ),

                "is_primary_candidate":
                    (
                        "yes"
                        if is_primary
                        else "no"
                    ),

                "is_secondary_candidate":
                    (
                        "yes"
                        if is_secondary
                        else "no"
                    ),

                "candidate_gene_distance_bp":
                    candidate_distance,

                "candidate_gene_relation":
                    candidate_relation,

                # --------------------------------------------
                # CRE-level candidate evidence
                # --------------------------------------------

                "candidate_priority":
                    row[
                        "candidate_priority"
                    ],

                "downstream_priority":
                    row[
                        "downstream_priority"
                    ],

                "n_total_tier1_clades":
                    optional_value(
                        row,
                        "n_total_tier1_clades",
                    ),

                "n_focal_tier1_clades":
                    optional_value(
                        row,
                        "n_focal_tier1_clades",
                    ),

                "n_secondary_only_clades":
                    optional_value(
                        row,
                        "n_secondary_only_clades",
                    ),

                "tier1_clades":
                    optional_value(
                        row,
                        "tier1_clades",
                    ),

                # --------------------------------------------
                # Original CRE gene assignment
                # --------------------------------------------

                "fbgn_target_genes":
                    normalize_missing(
                        row[
                            "fbgn_target_genes"
                        ]
                    ),

                "primary_candidate_fbgn":
                    primary_fbgn,

                "secondary_candidate_fbgn":
                    secondary_fbgn,

                "gene_assignment_qc":
                    normalize_missing(
                        row[
                            "gene_assignment_qc"
                        ]
                    ),
            })


    # ========================================================
    # Build output
    # ========================================================

    out = pd.DataFrame(
        rows
    )


    # ========================================================
    # QC
    # ========================================================

    if out.empty:
        raise SystemExit(
            "ERROR: no CRE x FBgn rows generated."
        )


    # --------------------------------------------------------
    # Each CRE x FBgn pair must be unique
    # --------------------------------------------------------

    if out.duplicated(
        subset=[
            "dmel_cre_id",
            "fbgn",
        ]
    ).any():

        duplicated = (
            out.loc[
                out.duplicated(
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
            "ERROR: duplicate CRE x FBgn pairs generated:\n\n"
            + duplicated.to_string(
                index=False
            )
        )


    # --------------------------------------------------------
    # Exactly one primary gene per CRE
    # --------------------------------------------------------

    primary_counts = (
        out.loc[
            out[
                "is_primary_candidate"
            ]
            == "yes"
        ]
        .groupby(
            "dmel_cre_id"
        )
        .size()
    )


    missing_primary = sorted(
        set(
            df[
                "dmel_cre_id"
            ]
        )
        - set(
            primary_counts.index
        )
    )


    if missing_primary:

        raise SystemExit(
            "ERROR: CREs without primary candidate "
            "in exploded table:\n"
            + "\n".join(
                missing_primary
            )
        )


    if (
        primary_counts
        != 1
    ).any():

        bad = (
            primary_counts.loc[
                primary_counts
                != 1
            ]
        )

        raise SystemExit(
            "ERROR: CREs with !=1 primary candidate:\n\n"
            + bad.to_string()
        )


    # --------------------------------------------------------
    # At most one secondary gene per CRE
    # --------------------------------------------------------

    secondary_counts = (
        out.loc[
            out[
                "is_secondary_candidate"
            ]
            == "yes"
        ]
        .groupby(
            "dmel_cre_id"
        )
        .size()
    )


    if (
        secondary_counts
        > 1
    ).any():

        bad = (
            secondary_counts.loc[
                secondary_counts
                > 1
            ]
        )

        raise SystemExit(
            "ERROR: CREs with >1 secondary candidate:\n\n"
            + bad.to_string()
        )


    # ========================================================
    # Stable deterministic ordering
    # ========================================================

    downstream_order = {
        "high": 1,
        "medium": 2,
        "exploratory": 3,
    }

    role_order = {
        "primary_candidate": 1,
        "secondary_candidate": 2,
        "reference_target_only": 3,
    }


    out[
        "_downstream_rank"
    ] = (
        out[
            "downstream_priority"
        ]
        .map(
            downstream_order
        )
        .fillna(
            99
        )
    )


    out[
        "_role_rank"
    ] = (
        out[
            "gene_role"
        ]
        .map(
            role_order
        )
        .fillna(
            99
        )
    )


    out = (
        out
        .sort_values(
            [
                "_downstream_rank",
                "dmel_cre_id",
                "_role_rank",
                "fbgn",
            ],
            kind="mergesort",
        )
        .drop(
            columns=[
                "_downstream_rank",
                "_role_rank",
            ]
        )
        .reset_index(
            drop=True
        )
    )


    # ========================================================
    # Write output
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
    print("Candidate FBgn explosion complete")
    print("=" * 72)

    print(
        f"Candidate CREs: "
        f"{df['dmel_cre_id'].nunique()}"
    )

    print(
        f"CRE x FBgn rows: "
        f"{len(out)}"
    )

    print(
        f"Unique FBgns: "
        f"{out['fbgn'].nunique()}"
    )

    print()

    print(
        "Gene roles:"
    )

    print(
        out[
            "gene_role"
        ]
        .value_counts()
        .to_string()
    )

    print()

    print(
        "Assignment evidence:"
    )

    print(
        out[
            "assignment_evidence"
        ]
        .value_counts()
        .to_string()
    )

    print()
    print(
        f"Wrote long-form table:\n"
        f"{args.out}"
    )

# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    main()
