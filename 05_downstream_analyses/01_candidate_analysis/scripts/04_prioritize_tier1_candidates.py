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

MAP_DIR = (
    PROJECT_DIR
    / "mapping_orthologs"
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

QC_DIR = (
    RESULTS_DIR
    / "qc"
)

TABLE_DIR = (
    RESULTS_DIR
    / "tables"
)

TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

DEFAULT_QC = (
    QC_DIR
    / "tier1_candidate_qc_summary.tsv"
)

DEFAULT_REF = (
    MAP_DIR
    / "reference_cres/dmel_reference_cres.tsv"
)


DEFAULT_OUT_EVIDENCE = (
    TABLE_DIR
    / "tier1_candidate_clade_evidence.tsv"
)

DEFAULT_OUT_PRIORITIZED = (
    TABLE_DIR
    / "tier1_candidates_prioritized.tsv"
)

DEFAULT_METADATA = (
    TABLE_DIR
    / "tier1_candidate_prioritization_metadata.tsv"
)


# ============================================================
# Priority definitions
# ============================================================

FOCAL_RECURRENT_MIN_CLADES = 2
SECONDARY_RECURRENT_MIN_CLADES = 2


# ============================================================
# Required columns
# ============================================================

QC_REQUIRED_COLUMNS = {
    "group_name",
    "dmel_cre_id",
    "tier1_scope",
    "is_focal_tier1",
    "is_singleton_tier1",
    "discordant_species",
    "discordant_state",
    "consensus_state",
    "group_species",
    "candidate_qc",
}


REF_COLUMNS = [
    "dmel_cre_id",
    "chrom",
    "start0",
    "end0",
    "peak_amplitude",
    "scrmshaw_score",
    "fbgn_target_genes",
    "n_fbgn_target_genes",
    "training_set",
    "method",
    "rank",
]


# ============================================================
# Argument parser
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Consolidate QC-passing focal and secondary Tier-1 "
            "CRE evidence and derive reproducible candidate "
            "priority categories."
        )
    )

    parser.add_argument(
        "--qc",
        type=Path,
        default=DEFAULT_QC,
        help=(
            "Tier-1 candidate QC summary "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--reference",
        type=Path,
        default=DEFAULT_REF,
        help=(
            "D. melanogaster reference CRE annotation "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--out-evidence",
        type=Path,
        default=DEFAULT_OUT_EVIDENCE,
        help=(
            "Clade-level Tier-1 evidence table "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--out-prioritized",
        type=Path,
        default=DEFAULT_OUT_PRIORITIZED,
        help=(
            "One-row-per-CRE prioritized table "
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
            "ERROR: required file not found:\n"
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
            f"ERROR: missing required columns in {label}:\n"
            + "\n".join(
                missing
            )
        )


def join_unique(values):
    """
    Join unique non-empty values in deterministic order.
    """

    clean = {
        str(value).strip()
        for value in values
        if str(value).strip().lower()
        not in {
            "",
            "na",
            "nan",
            "none",
        }
    }

    return "|".join(
        sorted(clean)
    )


def assign_priority(
    n_focal_clades,
    n_secondary_only_clades,
):
    """
    Assign descriptive, reproducible candidate priority.

    Priority hierarchy
    ------------------
    focal_recurrent
        Focal Tier-1 in >=2 clades.

    focal_plus_secondary_recurrent
        Focal Tier-1 in >=1 clade and secondary-only Tier-1
        in >=1 additional clade.

    focal_single
        Exactly one focal Tier-1 clade and no secondary-only
        support elsewhere.

    secondary_recurrent
        No focal Tier-1 support, but secondary-only Tier-1
        in >=2 clades.

    secondary_single
        One secondary-only Tier-1 clade only.
    """

    n_focal_clades = int(
        n_focal_clades
    )

    n_secondary_only_clades = int(
        n_secondary_only_clades
    )

    if (
        n_focal_clades
        >= FOCAL_RECURRENT_MIN_CLADES
    ):
        return "focal_recurrent"

    if (
        n_focal_clades >= 1
        and
        n_secondary_only_clades >= 1
    ):
        return (
            "focal_plus_secondary_recurrent"
        )

    if (
        n_focal_clades == 1
        and
        n_secondary_only_clades == 0
    ):
        return "focal_single"

    if (
        n_focal_clades == 0
        and
        n_secondary_only_clades
        >= SECONDARY_RECURRENT_MIN_CLADES
    ):
        return "secondary_recurrent"

    if (
        n_focal_clades == 0
        and
        n_secondary_only_clades == 1
    ):
        return "secondary_single"

    return "unclassified"


def priority_rank(priority):
    """
    Deterministic display/sorting rank.

    Lower value = shown earlier.
    """

    ranks = {
        "focal_recurrent": 1,
        "focal_plus_secondary_recurrent": 2,
        "secondary_recurrent": 3,
        "focal_single": 4,
        "secondary_single": 5,
        "unclassified": 99,
    }

    return ranks.get(
        priority,
        99,
    )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    require_file(
        args.qc
    )

    require_file(
        args.reference
    )

    args.out_evidence.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.out_prioritized.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.metadata_out.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ========================================================
    # Load QC summary
    # ========================================================

    qc = pd.read_csv(
        args.qc,
        sep="\t",
        dtype=str,
    ).fillna("")

    require_columns(
        qc,
        QC_REQUIRED_COLUMNS,
        "Tier-1 candidate QC summary",
    )

    # --------------------------------------------------------
    # Each group x CRE should occur once
    # --------------------------------------------------------

    if qc.duplicated(
        subset=[
            "group_name",
            "dmel_cre_id",
        ]
    ).any():

        duplicated = (
            qc.loc[
                qc.duplicated(
                    subset=[
                        "group_name",
                        "dmel_cre_id",
                    ],
                    keep=False,
                ),
                [
                    "group_name",
                    "dmel_cre_id",
                ],
            ]
            .drop_duplicates()
        )

        raise SystemExit(
            "ERROR: duplicate group x CRE rows "
            "in Tier-1 QC table:\n\n"
            + duplicated.to_string(
                index=False
            )
        )

    # ========================================================
    # Validate scope vocabulary
    # ========================================================

    valid_scopes = {
        "focal",
        "secondary_only",
    }

    observed_scopes = set(
        qc[
            "tier1_scope"
        ]
    )

    invalid_scopes = sorted(
        observed_scopes
        - valid_scopes
    )

    if invalid_scopes:

        raise SystemExit(
            "ERROR: unexpected tier1_scope values:\n"
            + "\n".join(
                invalid_scopes
            )
        )

    # ========================================================
    # Retain only QC-passing evidence
    # ========================================================

    evidence = (
        qc.loc[
            qc[
                "candidate_qc"
            ]
            == "PASS"
        ]
        .copy()
    )

    # --------------------------------------------------------
    # Stable ordering
    # --------------------------------------------------------

    if not evidence.empty:

        evidence = (
            evidence
            .sort_values(
                [
                    "dmel_cre_id",
                    "group_name",
                    "tier1_scope",
                ],
                kind="mergesort",
            )
            .reset_index(
                drop=True
            )
        )

    # ========================================================
    # Load canonical Dmel reference annotation
    # ========================================================

    ref = pd.read_csv(
        args.reference,
        sep="\t",
        dtype=str,
    ).fillna("")

    require_columns(
        ref,
        REF_COLUMNS,
        "Dmel reference CRE table",
    )

    if ref[
        "dmel_cre_id"
    ].duplicated().any():

        duplicated = (
            ref.loc[
                ref[
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
            "ERROR: duplicated CRE IDs in "
            "Dmel reference table:\n"
            + "\n".join(
                sorted(
                    duplicated
                )
            )
        )

    ref = ref[
        REF_COLUMNS
    ].copy()

    # ========================================================
    # Ensure every candidate exists in reference
    # ========================================================

    candidate_ids = set(
        evidence[
            "dmel_cre_id"
        ]
    )

    reference_ids = set(
        ref[
            "dmel_cre_id"
        ]
    )

    missing_reference = sorted(
        candidate_ids
        - reference_ids
    )

    if missing_reference:

        raise SystemExit(
            "ERROR: QC-passing Tier-1 CREs missing "
            "from Dmel reference table:\n"
            + "\n".join(
                missing_reference
            )
        )

    # ========================================================
    # Add reference annotation to evidence table
    # ========================================================

    # Avoid duplicating columns already provided by Step 03.
    ref_columns_to_add = [
        column
        for column in REF_COLUMNS
        if (
            column != "dmel_cre_id"
            and column not in evidence.columns
        )
    ]

    evidence = evidence.merge(
        ref[
            [
                "dmel_cre_id",
                *ref_columns_to_add,
            ]
        ],
        on="dmel_cre_id",
        how="left",
        validate="many_to_one",
    )

    # ========================================================
    # Aggregate one row per CRE
    # ========================================================

    prioritized_rows = []

    for cre_id, cre_df in (
        evidence.groupby(
            "dmel_cre_id",
            sort=True,
        )
    ):

        focal_df = cre_df.loc[
            cre_df[
                "tier1_scope"
            ]
            == "focal"
        ]

        secondary_df = cre_df.loc[
            cre_df[
                "tier1_scope"
            ]
            == "secondary_only"
        ]

        # ----------------------------------------------------
        # Counts
        # ----------------------------------------------------

        total_clades = (
            cre_df[
                "group_name"
            ].nunique()
        )

        focal_clades = (
            focal_df[
                "group_name"
            ].nunique()
        )

        secondary_clades = (
            secondary_df[
                "group_name"
            ].nunique()
        )

        # Sanity check:
        # focal and secondary_only must be mutually exclusive
        # at group x CRE level.
        if (
            total_clades
            != focal_clades
            + secondary_clades
        ):

            raise SystemExit(
                "ERROR: focal + secondary clade counts "
                "do not equal total for:\n"
                f"{cre_id}"
            )

        # ----------------------------------------------------
        # Priority
        # ----------------------------------------------------

        priority = assign_priority(
            focal_clades,
            secondary_clades,
        )

        # ----------------------------------------------------
        # Direction of Tier-1 contrasts
        # ----------------------------------------------------

        contrast_directions = (
            cre_df[
                "discordant_state"
            ]
            .astype(str)
            + "->"
            + cre_df[
                "consensus_state"
            ]
            .astype(str)
        )

        # ----------------------------------------------------
        # Reference row
        # ----------------------------------------------------

        ref_row = (
            ref.loc[
                ref[
                    "dmel_cre_id"
                ]
                == cre_id
            ]
            .iloc[0]
        )

        # ----------------------------------------------------
        # Build row
        # ----------------------------------------------------

        prioritized_rows.append({
            "dmel_cre_id":
                cre_id,

            "candidate_priority":
                priority,

            "priority_rank":
                priority_rank(
                    priority
                ),

            "n_total_tier1_clades":
                total_clades,

            "tier1_clades":
                join_unique(
                    cre_df[
                        "group_name"
                    ]
                ),

            "n_focal_tier1_clades":
                focal_clades,

            "focal_tier1_clades":
                join_unique(
                    focal_df[
                        "group_name"
                    ]
                ),

            "n_secondary_only_clades":
                secondary_clades,

            "secondary_only_clades":
                join_unique(
                    secondary_df[
                        "group_name"
                    ]
                ),

            "has_focal_support":
                (
                    "yes"
                    if focal_clades > 0
                    else "no"
                ),

            "has_secondary_support":
                (
                    "yes"
                    if secondary_clades > 0
                    else "no"
                ),

            "discordant_species":
                join_unique(
                    cre_df[
                        "discordant_species"
                    ]
                ),

            "focal_discordant_species":
                join_unique(
                    focal_df[
                        "discordant_species"
                    ]
                ),

            "secondary_discordant_species":
                join_unique(
                    secondary_df[
                        "discordant_species"
                    ]
                ),

            "discordant_states":
                join_unique(
                    cre_df[
                        "discordant_state"
                    ]
                ),

            "consensus_states":
                join_unique(
                    cre_df[
                        "consensus_state"
                    ]
                ),

            "tier1_state_directions":
                join_unique(
                    contrast_directions
                ),

            "fbgn_target_genes":
                ref_row[
                    "fbgn_target_genes"
                ],

            "n_fbgn_target_genes":
                ref_row[
                    "n_fbgn_target_genes"
                ],

            "chrom":
                ref_row[
                    "chrom"
                ],

            "start0":
                ref_row[
                    "start0"
                ],

            "end0":
                ref_row[
                    "end0"
                ],

            "peak_amplitude":
                ref_row[
                    "peak_amplitude"
                ],

            "scrmshaw_score":
                ref_row[
                    "scrmshaw_score"
                ],

            "rank":
                ref_row[
                    "rank"
                ],

            "training_set":
                ref_row[
                    "training_set"
                ],

            "method":
                ref_row[
                    "method"
                ],
        })

    # ========================================================
    # Prioritized DataFrame
    # ========================================================

    prioritized = pd.DataFrame(
        prioritized_rows
    )

    # ========================================================
    # Stable deterministic sorting
    # ========================================================

    if not prioritized.empty:

        numeric_columns = [
            "priority_rank",
            "n_total_tier1_clades",
            "n_focal_tier1_clades",
            "n_secondary_only_clades",
        ]

        for column in numeric_columns:

            prioritized[
                column
            ] = pd.to_numeric(
                prioritized[
                    column
                ],
                errors="raise",
            )

        prioritized = (
            prioritized
            .sort_values(
                [
                    "priority_rank",
                    "n_total_tier1_clades",
                    "n_focal_tier1_clades",
                    "n_secondary_only_clades",
                    "dmel_cre_id",
                ],
                ascending=[
                    True,
                    False,
                    False,
                    False,
                    True,
                ],
                kind="mergesort",
            )
            .reset_index(
                drop=True
            )
        )

    # ========================================================
    # Final QC
    # ========================================================

    if not prioritized.empty:

        if prioritized[
            "dmel_cre_id"
        ].duplicated().any():

            raise SystemExit(
                "ERROR: duplicate CRE IDs in prioritized table."
            )

        if (
            prioritized[
                "candidate_priority"
            ]
            == "unclassified"
        ).any():

            bad = prioritized.loc[
                prioritized[
                    "candidate_priority"
                ]
                == "unclassified",
                [
                    "dmel_cre_id",
                    "n_focal_tier1_clades",
                    "n_secondary_only_clades",
                ],
            ]

            raise SystemExit(
                "ERROR: candidates could not be assigned "
                "a reproducible priority:\n\n"
                + bad.to_string(
                    index=False
                )
            )

    # ========================================================
    # Final evidence ordering
    # ========================================================

    if not evidence.empty:

        evidence = (
            evidence
            .sort_values(
                [
                    "dmel_cre_id",
                    "tier1_scope",
                    "group_name",
                ],
                kind="mergesort",
            )
            .reset_index(
                drop=True
            )
        )

    # ========================================================
    # Write outputs
    # ========================================================

    evidence.to_csv(
        args.out_evidence,
        sep="\t",
        index=False,
    )

    prioritized.to_csv(
        args.out_prioritized,
        sep="\t",
        index=False,
    )

    # ========================================================
    # Run metadata
    # ========================================================

    priority_counts = (
        prioritized[
            "candidate_priority"
        ]
        .value_counts()
        if not prioritized.empty
        else pd.Series(
            dtype=int
        )
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

            "focal_recurrent_min_clades":
                FOCAL_RECURRENT_MIN_CLADES,

            "secondary_recurrent_min_clades":
                SECONDARY_RECURRENT_MIN_CLADES,

            "qc_input":
                str(
                    args.qc.resolve()
                ),

            "reference_input":
                str(
                    args.reference.resolve()
                ),

            "n_qc_rows":
                len(
                    qc
                ),

            "n_qc_pass_group_cre_rows":
                len(
                    evidence
                ),

            "n_unique_qc_pass_cres":
                (
                    evidence[
                        "dmel_cre_id"
                    ].nunique()
                    if not evidence.empty
                    else 0
                ),

            "n_focal_recurrent":
                int(
                    priority_counts.get(
                        "focal_recurrent",
                        0,
                    )
                ),

            "n_focal_plus_secondary_recurrent":
                int(
                    priority_counts.get(
                        "focal_plus_secondary_recurrent",
                        0,
                    )
                ),

            "n_secondary_recurrent":
                int(
                    priority_counts.get(
                        "secondary_recurrent",
                        0,
                    )
                ),

            "n_focal_single":
                int(
                    priority_counts.get(
                        "focal_single",
                        0,
                    )
                ),

            "n_secondary_single":
                int(
                    priority_counts.get(
                        "secondary_single",
                        0,
                    )
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
    print("Tier-1 candidate prioritization complete")
    print("=" * 72)

    print(
        f"QC input rows: "
        f"{len(qc)}"
    )

    print(
        f"QC-passing group x CRE rows: "
        f"{len(evidence)}"
    )

    print(
        "Unique QC-passing CREs: "
        f"{evidence['dmel_cre_id'].nunique() if not evidence.empty else 0}"
    )

    if not prioritized.empty:

        print()
        print(
            "Candidate priorities:"
        )

        print(
            prioritized[
                "candidate_priority"
            ]
            .value_counts()
            .to_string()
        )

        recurrent = (
            prioritized.loc[
                prioritized[
                    "n_total_tier1_clades"
                ]
                > 1
            ]
        )

        if not recurrent.empty:

            print()
            print(
                "Recurrent Tier-1 candidates:"
            )

            print(
                recurrent[
                    [
                        "dmel_cre_id",
                        "candidate_priority",
                        "n_total_tier1_clades",
                        "n_focal_tier1_clades",
                        "n_secondary_only_clades",
                        "tier1_clades",
                        "discordant_species",
                        "fbgn_target_genes",
                    ]
                ].to_string(
                    index=False
                )
            )

    print()
    print(
        f"Wrote clade-level evidence:\n"
        f"{args.out_evidence}"
    )

    print()

    print(
        f"Wrote prioritized candidate table:\n"
        f"{args.out_prioritized}"
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
