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
    / "tier1_candidate_gene_distances.tsv"
)

DEFAULT_OUTPUT = (
    TABLE_DIR
    / "tier1_candidate_gene_assignments.tsv"
)

DEFAULT_METADATA = (
    TABLE_DIR
    / "tier1_candidate_gene_assignments_metadata.tsv"
)


# ============================================================
# Required columns
# ============================================================

REQUIRED_COLUMNS = {
    "dmel_cre_id",
    "candidate_priority",
    "downstream_priority",
    "fbgn_target_genes",

    "dmel_fbgn_flanking_gene",
    "distance_flanking_gene",

    "dmel_fbgn_next_flanking_gene",
    "distance_next_gene",

    "gene_distance_qc",
}


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Assign primary and secondary candidate genes "
            "to Tier-1 CRE candidates using D. melanogaster "
            "flanking-gene distances."
        )
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=(
            "Candidate gene-distance table "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=(
            "Candidate gene-assignment output "
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


def parse_target_genes(value):
    """
    Convert pipe-separated FBgn target genes into a set.
    """

    raw = normalize_missing(
        value
    )

    if raw == "NA":
        return set()

    return {
        gene.strip()
        for gene in raw.split("|")
        if gene.strip()
        and gene.strip().lower()
        not in {
            "na",
            "nan",
            "none",
        }
    }


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


def classify_relation(
    distance,
):
    """
    Describe CRE-gene relationship based on distance.

    distance == 0
        overlapping

    distance > 0
        nearest_nonoverlapping

    missing
        unknown
    """

    if pd.isna(
        distance
    ):
        return "unknown"

    if float(
        distance
    ) == 0:
        return "overlapping"

    return "nearest_nonoverlapping"


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
    # Load
    # ========================================================

    df = pd.read_csv(
        args.input,
        sep="\t",
        dtype=str,
    ).fillna("")

    require_columns(
        df,
        REQUIRED_COLUMNS,
        "candidate gene-distance table",
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
    # Numeric distances
    # ========================================================

    for column in [
        "distance_flanking_gene",
        "distance_next_gene",
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
    # Assign genes
    # ========================================================

    assignment_rows = []


    for _, row in (
        df.iterrows()
    ):

        cre_id = row[
            "dmel_cre_id"
        ]

        target_genes = parse_target_genes(
            row[
                "fbgn_target_genes"
            ]
        )


        # ----------------------------------------------------
        # Raw distance-derived genes
        # ----------------------------------------------------

        fbgn1 = normalize_missing(
            row[
                "dmel_fbgn_flanking_gene"
            ]
        )

        fbgn2 = normalize_missing(
            row[
                "dmel_fbgn_next_flanking_gene"
            ]
        )

        d1 = row[
            "distance_flanking_gene"
        ]

        d2 = row[
            "distance_next_gene"
        ]


        # ----------------------------------------------------
        # Normalize impossible gene/distance combinations
        # ----------------------------------------------------

        if fbgn1 == "NA":
            d1 = pd.NA

        if fbgn2 == "NA":
            d2 = pd.NA


        # ====================================================
        # Determine primary and secondary candidate gene
        #
        # Primary rule:
        #   smallest genomic distance to CRE
        #
        # Tie rule:
        #   deterministic lexical FBgn ordering
        # ====================================================

        candidates = []

        if (
            fbgn1 != "NA"
            and not pd.isna(
                d1
            )
        ):

            candidates.append({
                "fbgn":
                    fbgn1,

                "distance":
                    float(
                        d1
                    ),

                "source":
                    "flanking_gene",
            })


        if (
            fbgn2 != "NA"
            and not pd.isna(
                d2
            )
        ):

            candidates.append({
                "fbgn":
                    fbgn2,

                "distance":
                    float(
                        d2
                    ),

                "source":
                    "next_flanking_gene",
            })


        # ----------------------------------------------------
        # Deduplicate identical FBgn assignments
        #
        # If the same FBgn appears twice, retain its shortest
        # distance and one deterministic source label.
        # ----------------------------------------------------

        collapsed = {}

        for candidate in candidates:

            fbgn = candidate[
                "fbgn"
            ]

            if (
                fbgn
                not in collapsed
            ):

                collapsed[
                    fbgn
                ] = candidate

            else:

                existing = collapsed[
                    fbgn
                ]

                if (
                    candidate[
                        "distance"
                    ]
                    <
                    existing[
                        "distance"
                    ]
                ):

                    collapsed[
                        fbgn
                    ] = candidate


        candidates = list(
            collapsed.values()
        )


        # ----------------------------------------------------
        # Stable deterministic ranking
        # ----------------------------------------------------

        candidates = sorted(
            candidates,
            key=lambda x: (
                x[
                    "distance"
                ],
                x[
                    "fbgn"
                ],
                x[
                    "source"
                ],
            ),
        )


        # ====================================================
        # Primary / secondary assignment
        # ====================================================

        if len(
            candidates
        ) == 0:

            primary_fbgn = "NA"
            primary_distance = pd.NA
            primary_source = "NA"

            secondary_fbgn = "NA"
            secondary_distance = pd.NA
            secondary_source = "NA"


        elif len(
            candidates
        ) == 1:

            primary_fbgn = (
                candidates[
                    0
                ][
                    "fbgn"
                ]
            )

            primary_distance = (
                candidates[
                    0
                ][
                    "distance"
                ]
            )

            primary_source = (
                candidates[
                    0
                ][
                    "source"
                ]
            )

            secondary_fbgn = "NA"
            secondary_distance = pd.NA
            secondary_source = "NA"


        else:

            primary_fbgn = (
                candidates[
                    0
                ][
                    "fbgn"
                ]
            )

            primary_distance = (
                candidates[
                    0
                ][
                    "distance"
                ]
            )

            primary_source = (
                candidates[
                    0
                ][
                    "source"
                ]
            )

            secondary_fbgn = (
                candidates[
                    1
                ][
                    "fbgn"
                ]
            )

            secondary_distance = (
                candidates[
                    1
                ][
                    "distance"
                ]
            )

            secondary_source = (
                candidates[
                    1
                ][
                    "source"
                ]
            )


        # ====================================================
        # Relations
        # ====================================================

        primary_relation = (
            classify_relation(
                primary_distance
            )
        )

        secondary_relation = (
            classify_relation(
                secondary_distance
            )
        )


        # ====================================================
        # Gene assignment QC
        # ====================================================

        qc_flags = []


        # ----------------------------------------------------
        # Carry forward gene-distance QC
        # ----------------------------------------------------

        distance_qc = normalize_missing(
            row[
                "gene_distance_qc"
            ]
        )

        if (
            distance_qc
            != "PASS"
        ):

            qc_flags.append(
                f"GENE_DISTANCE_QC_{distance_qc}"
            )


        # ----------------------------------------------------
        # No distance-derived gene
        # ----------------------------------------------------

        if (
            primary_fbgn
            == "NA"
        ):

            qc_flags.append(
                "NO_DISTANCE_DERIVED_GENE"
            )


        # ----------------------------------------------------
        # Check against original SCRMshaw FBgn targets
        # ----------------------------------------------------

        if (
            primary_fbgn
            != "NA"
            and
            primary_fbgn
            not in target_genes
        ):

            qc_flags.append(
                "PRIMARY_FBGN_NOT_IN_REFERENCE_TARGETS"
            )


        if (
            secondary_fbgn
            != "NA"
            and
            secondary_fbgn
            not in target_genes
        ):

            qc_flags.append(
                "SECONDARY_FBGN_NOT_IN_REFERENCE_TARGETS"
            )


        # ----------------------------------------------------
        # If reference target genes are missing entirely
        # ----------------------------------------------------

        if not target_genes:

            qc_flags.append(
                "NO_REFERENCE_TARGET_GENES"
            )


        # ----------------------------------------------------
        # Equal distance tie
        # ----------------------------------------------------

        distance_tie = (
            primary_fbgn != "NA"
            and
            secondary_fbgn != "NA"
            and
            not pd.isna(
                primary_distance
            )
            and
            not pd.isna(
                secondary_distance
            )
            and
            float(
                primary_distance
            )
            ==
            float(
                secondary_distance
            )
        )


        if distance_tie:

            qc_flags.append(
                "PRIMARY_SECONDARY_DISTANCE_TIE"
            )


        # ----------------------------------------------------
        # Determine whether primary gene is one of original
        # SCRMshaw target-gene assignments
        # ----------------------------------------------------

        primary_in_reference = (
            primary_fbgn != "NA"
            and
            primary_fbgn
            in target_genes
        )


        secondary_in_reference = (
            secondary_fbgn != "NA"
            and
            secondary_fbgn
            in target_genes
        )


        # ----------------------------------------------------
        # Final QC
        # ----------------------------------------------------

        qc_flags = sorted(
            set(
                qc_flags
            )
        )


        assignment_qc = (
            "PASS"
            if not qc_flags
            else ";".join(
                qc_flags
            )
        )


        # ====================================================
        # Output row
        # ====================================================

        assignment_rows.append({
            "dmel_cre_id":
                cre_id,

            "candidate_priority":
                row[
                    "candidate_priority"
                ],

            "downstream_priority":
                row[
                    "downstream_priority"
                ],

            "n_total_tier1_clades":
                normalize_missing(
                    row.get(
                        "n_total_tier1_clades",
                        "NA",
                    )
                ),

            "n_focal_tier1_clades":
                normalize_missing(
                    row.get(
                        "n_focal_tier1_clades",
                        "NA",
                    )
                ),

            "n_secondary_only_clades":
                normalize_missing(
                    row.get(
                        "n_secondary_only_clades",
                        "NA",
                    )
                ),

            "tier1_clades":
                normalize_missing(
                    row.get(
                        "tier1_clades",
                        "NA",
                    )
                ),

            "fbgn_target_genes":
                normalize_missing(
                    row[
                        "fbgn_target_genes"
                    ]
                ),

            "n_reference_target_genes":
                len(
                    target_genes
                ),

            "primary_candidate_fbgn":
                primary_fbgn,

            "primary_candidate_distance_bp":
                primary_distance,

            "primary_candidate_relation":
                primary_relation,

            "primary_candidate_source":
                primary_source,

            "primary_in_reference_targets":
                (
                    "yes"
                    if primary_in_reference
                    else "no"
                ),

            "secondary_candidate_fbgn":
                secondary_fbgn,

            "secondary_candidate_distance_bp":
                secondary_distance,

            "secondary_candidate_relation":
                secondary_relation,

            "secondary_candidate_source":
                secondary_source,

            "secondary_in_reference_targets":
                (
                    "yes"
                    if secondary_in_reference
                    else "no"
                ),

            "primary_secondary_distance_tie":
                (
                    "yes"
                    if distance_tie
                    else "no"
                ),

            # Raw provenance
            "dmel_fbgn_flanking_gene":
                fbgn1,

            "distance_flanking_gene":
                d1,

            "dmel_fbgn_next_flanking_gene":
                fbgn2,

            "distance_next_gene":
                d2,

            "gene_distance_qc":
                distance_qc,

            "gene_assignment_qc":
                assignment_qc,
        })


    # ========================================================
    # Build output
    # ========================================================

    out = pd.DataFrame(
        assignment_rows
    )


    # --------------------------------------------------------
    # Stable ordering
    # --------------------------------------------------------

    priority_order = {
        "high": 1,
        "medium": 2,
        "exploratory": 3,
    }


    out[
        "_priority_rank"
    ] = (
        out[
            "downstream_priority"
        ]
        .map(
            priority_order
        )
        .fillna(
            99
        )
    )


    out = (
        out
        .sort_values(
            [
                "_priority_rank",
                "dmel_cre_id",
            ],
            kind="mergesort",
        )
        .drop(
            columns=[
                "_priority_rank",
            ]
        )
        .reset_index(
            drop=True
        )
    )


    # --------------------------------------------------------
    # Final uniqueness QC
    # --------------------------------------------------------

    if out[
        "dmel_cre_id"
    ].duplicated().any():

        raise SystemExit(
            "ERROR: duplicate CRE IDs in gene-assignment output."
        )


    # ========================================================
    # Write
    # ========================================================

    out.to_csv(
        args.out,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Metadata
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

            "input_file":
                str(
                    args.input.resolve()
                ),

            "input_sha256":
                file_sha256(
                    args.input
                ),

            "n_candidate_cres":
                len(
                    out
                ),

            "n_primary_assigned":
                int(
                    (
                        out[
                            "primary_candidate_fbgn"
                        ]
                        != "NA"
                    ).sum()
                ),

            "n_secondary_assigned":
                int(
                    (
                        out[
                            "secondary_candidate_fbgn"
                        ]
                        != "NA"
                    ).sum()
                ),

            "n_assignment_pass":
                int(
                    (
                        out[
                            "gene_assignment_qc"
                        ]
                        == "PASS"
                    ).sum()
                ),

            "n_primary_in_reference_targets":
                int(
                    (
                        out[
                            "primary_in_reference_targets"
                        ]
                        == "yes"
                    ).sum()
                ),

            "n_secondary_in_reference_targets":
                int(
                    (
                        out[
                            "secondary_in_reference_targets"
                        ]
                        == "yes"
                    ).sum()
                ),

            "n_distance_ties":
                int(
                    (
                        out[
                            "primary_secondary_distance_tie"
                        ]
                        == "yes"
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
    print("Candidate gene assignment complete")
    print("=" * 72)

    print(
        f"Candidate CREs: "
        f"{len(out)}"
    )

    print(
        "Primary candidate genes assigned: "
        f"{(out['primary_candidate_fbgn'] != 'NA').sum()}"
    )

    print(
        "Secondary candidate genes assigned: "
        f"{(out['secondary_candidate_fbgn'] != 'NA').sum()}"
    )

    print()

    print(
        "Gene assignment QC:"
    )

    print(
        out[
            "gene_assignment_qc"
        ]
        .value_counts()
        .to_string()
    )


    # --------------------------------------------------------
    # Show non-PASS assignments
    # --------------------------------------------------------

    flagged = out.loc[
        out[
            "gene_assignment_qc"
        ]
        != "PASS"
    ]


    if not flagged.empty:

        print()
        print(
            "Assignments requiring inspection:"
        )

        print(
            flagged[
                [
                    "dmel_cre_id",
                    "candidate_priority",
                    "fbgn_target_genes",

                    "primary_candidate_fbgn",
                    "primary_candidate_distance_bp",
                    "primary_in_reference_targets",

                    "secondary_candidate_fbgn",
                    "secondary_candidate_distance_bp",
                    "secondary_in_reference_targets",

                    "gene_assignment_qc",
                ]
            ].to_string(
                index=False
            )
        )


    print()
    print(
        f"Wrote gene assignments:\n"
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
