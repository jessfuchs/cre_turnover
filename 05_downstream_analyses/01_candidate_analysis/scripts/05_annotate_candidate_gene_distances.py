#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime
import argparse
import hashlib
import platform
import sys

import pandas as pd


# ============================================================
# Project paths
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

TABLE_DIR = (
    RESULTS_DIR
    / "tables"
)

TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


DEFAULT_CANDIDATES = (
    TABLE_DIR
    / "tier1_candidates_prioritized.tsv"
)

DEFAULT_SO = (
    MAP_DIR
    / "ortholog_results/SO_all_species_fbgn.tsv"
)


DEFAULT_OUT = (
    TABLE_DIR
    / "tier1_candidate_gene_distances.tsv"
)

DEFAULT_METADATA = (
    TABLE_DIR
    / "tier1_candidate_gene_distances_metadata.tsv"
)


# ============================================================
# Required columns
# ============================================================

CANDIDATE_REQUIRED_COLUMNS = {
    "dmel_cre_id",
    "candidate_priority",
    "chrom",
    "start0",
    "end0",
    "fbgn_target_genes",
}


SO_REQUIRED_COLUMNS = {
    "chrom",
    "start0",
    "end0",
}


# ============================================================
# Candidate downstream priority
# ============================================================

def downstream_priority(candidate_priority):
    """
    Broad downstream-analysis priority.

    high
        Recurrent candidate with focal support.

    medium
        Recurrent secondary candidate or a single focal
        candidate.

    exploratory
        Secondary-only single-clade candidate.
    """

    if candidate_priority in {
        "focal_recurrent",
        "focal_plus_secondary_recurrent",
    }:
        return "high"

    if candidate_priority in {
        "secondary_recurrent",
        "focal_single",
    }:
        return "medium"

    if candidate_priority == "secondary_single":
        return "exploratory"

    return "unclassified"


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Annotate prioritized D. melanogaster Tier-1 CRE "
            "candidates with flanking-gene information from "
            "the Dmel rows of the combined SCRMshaw orthology "
            "table."
        )
    )

    parser.add_argument(
        "--candidates",
        type=Path,
        default=DEFAULT_CANDIDATES,
        help=(
            "Prioritized Tier-1 candidate table "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--so-table",
        type=Path,
        default=DEFAULT_SO,
        help=(
            "Combined SO/SCRMshaw orthology table "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT,
        help=(
            "Gene-distance annotation output "
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
    """
    Convert missing/empty textual values to 'NA'.
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


def file_sha256(path):
    """
    Calculate SHA256 checksum for reproducibility metadata.
    """

    sha = hashlib.sha256()

    with open(path, "rb") as handle:

        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            sha.update(block)

    return sha.hexdigest()


def get_optional(
    row,
    column,
):
    """
    Safely retrieve an optional column value.
    """

    if column not in row.index:
        return "NA"

    return normalize_missing(
        row[column]
    )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    require_file(
        args.candidates
    )

    require_file(
        args.so_table
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
    # Load candidate table
    # ========================================================

    candidates = pd.read_csv(
        args.candidates,
        sep="\t",
        dtype=str,
    ).fillna("")

    require_columns(
        candidates,
        CANDIDATE_REQUIRED_COLUMNS,
        "prioritized Tier-1 candidate table",
    )


    # --------------------------------------------------------
    # One row per Dmel CRE is required
    # --------------------------------------------------------

    if candidates[
        "dmel_cre_id"
    ].duplicated().any():

        duplicated = (
            candidates.loc[
                candidates[
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
            "ERROR: duplicate CRE IDs in candidate table:\n"
            + "\n".join(
                sorted(duplicated)
            )
        )


    # --------------------------------------------------------
    # Numeric candidate coordinates
    # --------------------------------------------------------

    for column in [
        "start0",
        "end0",
    ]:

        candidates[
            column
        ] = pd.to_numeric(
            candidates[
                column
            ],
            errors="raise",
        )


    # --------------------------------------------------------
    # Validate candidate intervals
    # --------------------------------------------------------

    invalid_candidate_coords = (
        candidates[
            "end0"
        ]
        <= candidates[
            "start0"
        ]
    )

    if invalid_candidate_coords.any():

        bad = candidates.loc[
            invalid_candidate_coords,
            [
                "dmel_cre_id",
                "chrom",
                "start0",
                "end0",
            ],
        ]

        raise SystemExit(
            "ERROR: invalid candidate CRE coordinates:\n\n"
            + bad.to_string(
                index=False
            )
        )


    # --------------------------------------------------------
    # Add downstream priority
    # --------------------------------------------------------

    candidates[
        "downstream_priority"
    ] = candidates[
        "candidate_priority"
    ].apply(
        downstream_priority
    )

    if (
        candidates[
            "downstream_priority"
        ]
        == "unclassified"
    ).any():

        bad_priorities = sorted(
            candidates.loc[
                candidates[
                    "downstream_priority"
                ]
                == "unclassified",
                "candidate_priority",
            ]
            .unique()
        )

        raise SystemExit(
            "ERROR: unknown candidate_priority values:\n"
            + "\n".join(
                bad_priorities
            )
        )


    # ========================================================
    # Load combined SO table
    # ========================================================

    so = pd.read_csv(
        args.so_table,
        sep="\t",
        dtype=str,
    ).fillna("")

    require_columns(
        so,
        SO_REQUIRED_COLUMNS,
        "SO table",
    )


    if (
        "species_key"
        not in so.columns
        and
        "species"
        not in so.columns
    ):

        raise SystemExit(
            "ERROR: SO table must contain at least one "
            "species identifier column:\n"
            "species_key or species"
        )


    # ========================================================
    # Identify D. melanogaster records
    # ========================================================

    dmel_mask = pd.Series(
        False,
        index=so.index,
    )


    if "species_key" in so.columns:

        values = (
            so[
                "species_key"
            ]
            .astype(str)
            .str.strip()
        )

        dmel_mask |= (
            values.str.lower()
            .isin({
                "dmel",
                "d_melanogaster",
                "drosophila_melanogaster",
                "drosophila melanogaster",
            })
        )

        dmel_mask |= (
            values.str.contains(
                "melanogaster",
                case=False,
                regex=False,
                na=False,
            )
        )


    if "species" in so.columns:

        values = (
            so[
                "species"
            ]
            .astype(str)
            .str.strip()
        )

        dmel_mask |= (
            values.str.lower()
            .isin({
                "dmel",
                "d_melanogaster",
                "drosophila_melanogaster",
                "drosophila melanogaster",
            })
        )

        dmel_mask |= (
            values.str.contains(
                "melanogaster",
                case=False,
                regex=False,
                na=False,
            )
        )


    dmel = (
        so.loc[
            dmel_mask
        ]
        .copy()
    )


    if dmel.empty:

        # Print potentially useful species names before aborting.
        print(
            "Species values containing 'mel':"
        )

        for column in [
            "species_key",
            "species",
        ]:

            if column in so.columns:

                vals = (
                    so.loc[
                        so[column]
                        .astype(str)
                        .str.contains(
                            "mel",
                            case=False,
                            na=False,
                        ),
                        column,
                    ]
                    .drop_duplicates()
                    .tolist()
                )

                print(
                    f"{column}: {vals}"
                )

        raise SystemExit(
            "\nERROR: no D. melanogaster records "
            "identified in SO table."
        )


    print(
        f"Dmel rows identified in SO table: "
        f"{len(dmel)}"
    )


    # ========================================================
    # Numeric SO coordinates
    # ========================================================

    for column in [
        "start0",
        "end0",
        "distance_flanking_gene",
        "distance_next_gene",
    ]:

        if column in dmel.columns:

            dmel[
                column
            ] = pd.to_numeric(
                dmel[
                    column
                ],
                errors="coerce",
            )


    # --------------------------------------------------------
    # Remove SO rows with invalid genomic coordinates
    # --------------------------------------------------------

    invalid_so_coords = (
        dmel[
            "start0"
        ].isna()
        |
        dmel[
            "end0"
        ].isna()
        |
        (
            dmel[
                "end0"
            ]
            <= dmel[
                "start0"
            ]
        )
    )

    n_invalid_so_coords = int(
        invalid_so_coords.sum()
    )

    if n_invalid_so_coords:

        print(
            "WARNING: dropping "
            f"{n_invalid_so_coords} Dmel SO rows "
            "with invalid coordinates."
        )

    dmel = (
        dmel.loc[
            ~invalid_so_coords
        ]
        .copy()
    )


    # ========================================================
    # Stable deterministic Dmel ordering
    # ========================================================

    dmel = (
        dmel
        .sort_values(
            [
                "chrom",
                "start0",
                "end0",
            ],
            kind="mergesort",
        )
        .reset_index(
            drop=True
        )
    )


    # ========================================================
    # Match each Tier-1 candidate CRE to Dmel prediction
    # ========================================================

    output_rows = []


    candidates = (
        candidates
        .sort_values(
            "dmel_cre_id",
            kind="mergesort",
        )
        .reset_index(
            drop=True
        )
    )


    for _, candidate in (
        candidates.iterrows()
    ):

        cre_id = candidate[
            "dmel_cre_id"
        ]

        chrom = candidate[
            "chrom"
        ]

        cre_start = int(
            candidate[
                "start0"
            ]
        )

        cre_end = int(
            candidate[
                "end0"
            ]
        )

        cre_length = (
            cre_end
            - cre_start
        )


        # ----------------------------------------------------
        # Find all Dmel predictions overlapping the CRE
        # ----------------------------------------------------

        hits = (
            dmel.loc[
                (
                    dmel[
                        "chrom"
                    ]
                    == chrom
                )
                &
                (
                    dmel[
                        "start0"
                    ]
                    < cre_end
                )
                &
                (
                    dmel[
                        "end0"
                    ]
                    > cre_start
                )
            ]
            .copy()
        )


        # ----------------------------------------------------
        # No overlapping prediction
        # ----------------------------------------------------

        if hits.empty:

            output_rows.append({
                "dmel_cre_id":
                    cre_id,

                "candidate_priority":
                    candidate[
                        "candidate_priority"
                    ],

                "downstream_priority":
                    candidate[
                        "downstream_priority"
                    ],

                "n_total_tier1_clades":
                    get_optional(
                        candidate,
                        "n_total_tier1_clades",
                    ),

                "n_focal_tier1_clades":
                    get_optional(
                        candidate,
                        "n_focal_tier1_clades",
                    ),

                "n_secondary_only_clades":
                    get_optional(
                        candidate,
                        "n_secondary_only_clades",
                    ),

                "tier1_clades":
                    get_optional(
                        candidate,
                        "tier1_clades",
                    ),

                "fbgn_target_genes":
                    candidate[
                        "fbgn_target_genes"
                    ],

                "cre_chrom":
                    chrom,

                "cre_start0":
                    cre_start,

                "cre_end0":
                    cre_end,

                "cre_length_bp":
                    cre_length,

                "matched_dmel_prediction":
                    "NO",

                "n_overlapping_predictions":
                    0,

                "best_match_tied":
                    "no",

                "prediction_species_key":
                    "NA",

                "prediction_species":
                    "NA",

                "prediction_chrom":
                    "NA",

                "prediction_start0":
                    "NA",

                "prediction_end0":
                    "NA",

                "prediction_length_bp":
                    "NA",

                "overlap_bp":
                    "NA",

                "overlap_fraction_cre":
                    "NA",

                "overlap_fraction_prediction":
                    "NA",

                "flanking_gene":
                    "NA",

                "dmel_ortholog_flanking_gene":
                    "NA",

                "dmel_fbgn_flanking_gene":
                    "NA",

                "distance_flanking_gene":
                    "NA",

                "next_flanking_gene":
                    "NA",

                "dmel_ortholog_next_flanking_gene":
                    "NA",

                "dmel_fbgn_next_flanking_gene":
                    "NA",

                "distance_next_gene":
                    "NA",

                "gene_distance_qc":
                    "NO_OVERLAPPING_DMel_PREDICTION",
            })

            continue


        # ====================================================
        # Calculate interval overlaps
        # ====================================================

        hits[
            "overlap_bp"
        ] = (
            hits[
                "end0"
            ]
            .clip(
                upper=cre_end
            )
            -
            hits[
                "start0"
            ]
            .clip(
                lower=cre_start
            )
        )


        hits[
            "prediction_length_bp"
        ] = (
            hits[
                "end0"
            ]
            -
            hits[
                "start0"
            ]
        )


        hits[
            "overlap_fraction_cre"
        ] = (
            hits[
                "overlap_bp"
            ]
            / cre_length
        )


        hits[
            "overlap_fraction_prediction"
        ] = (
            hits[
                "overlap_bp"
            ]
            /
            hits[
                "prediction_length_bp"
            ]
        )


        # ====================================================
        # Deterministic best-match selection
        #
        # Priority:
        # 1. largest overlapping bp
        # 2. largest fraction of reference CRE
        # 3. largest fraction of prediction
        # 4. shortest coordinate distance between starts
        # 5. shortest prediction interval
        # 6. smallest genomic start
        # 7. smallest genomic end
        # ====================================================

        hits[
            "start_distance"
        ] = (
            hits[
                "start0"
            ]
            - cre_start
        ).abs()


        hits = (
            hits
            .sort_values(
                [
                    "overlap_bp",
                    "overlap_fraction_cre",
                    "overlap_fraction_prediction",
                    "start_distance",
                    "prediction_length_bp",
                    "start0",
                    "end0",
                ],
                ascending=[
                    False,
                    False,
                    False,
                    True,
                    True,
                    True,
                    True,
                ],
                kind="mergesort",
            )
            .reset_index(
                drop=True
            )
        )


        best = hits.iloc[
            0
        ]


        # ====================================================
        # Check whether another prediction had exactly the
        # same main overlap statistics
        # ====================================================

        tied = hits.loc[
            (
                hits[
                    "overlap_bp"
                ]
                == best[
                    "overlap_bp"
                ]
            )
            &
            (
                hits[
                    "overlap_fraction_cre"
                ]
                == best[
                    "overlap_fraction_cre"
                ]
            )
            &
            (
                hits[
                    "overlap_fraction_prediction"
                ]
                == best[
                    "overlap_fraction_prediction"
                ]
            )
        ]


        best_match_tied = (
            len(
                tied
            )
            > 1
        )


        # ====================================================
        # Match QC
        # ====================================================

        qc_flags = []


        if len(
            hits
        ) > 1:

            qc_flags.append(
                "MULTIPLE_OVERLAPPING_DMel_PREDICTIONS"
            )


        if best_match_tied:

            qc_flags.append(
                "TIED_PRIMARY_OVERLAP"
            )


        # Candidate CRE and matched Dmel prediction are expected
        # to represent the same original Dmel SCRMshaw region.
        # Partial overlap is therefore worth flagging.
        if (
            float(
                best[
                    "overlap_fraction_cre"
                ]
            )
            < 1.0
        ):

            qc_flags.append(
                "PARTIAL_CRE_OVERLAP"
            )


        # Exact coordinates are the strongest possible match.
        exact_coordinate_match = (
            int(
                best[
                    "start0"
                ]
            )
            == cre_start
            and
            int(
                best[
                    "end0"
                ]
            )
            == cre_end
        )


        if not exact_coordinate_match:

            qc_flags.append(
                "NON_EXACT_COORDINATE_MATCH"
            )


        gene_distance_qc = (
            "PASS"
            if not qc_flags
            else ";".join(
                sorted(
                    set(
                        qc_flags
                    )
                )
            )
        )


        # ====================================================
        # Output record
        # ====================================================

        output_rows.append({
            "dmel_cre_id":
                cre_id,

            "candidate_priority":
                candidate[
                    "candidate_priority"
                ],

            "downstream_priority":
                candidate[
                    "downstream_priority"
                ],

            "n_total_tier1_clades":
                get_optional(
                    candidate,
                    "n_total_tier1_clades",
                ),

            "n_focal_tier1_clades":
                get_optional(
                    candidate,
                    "n_focal_tier1_clades",
                ),

            "n_secondary_only_clades":
                get_optional(
                    candidate,
                    "n_secondary_only_clades",
                ),

            "tier1_clades":
                get_optional(
                    candidate,
                    "tier1_clades",
                ),

            "fbgn_target_genes":
                candidate[
                    "fbgn_target_genes"
                ],

            "cre_chrom":
                chrom,

            "cre_start0":
                cre_start,

            "cre_end0":
                cre_end,

            "cre_length_bp":
                cre_length,

            "matched_dmel_prediction":
                "YES",

            "n_overlapping_predictions":
                len(
                    hits
                ),

            "best_match_tied":
                (
                    "yes"
                    if best_match_tied
                    else "no"
                ),

            "exact_coordinate_match":
                (
                    "yes"
                    if exact_coordinate_match
                    else "no"
                ),

            "prediction_species_key":
                get_optional(
                    best,
                    "species_key",
                ),

            "prediction_species":
                get_optional(
                    best,
                    "species",
                ),

            "prediction_chrom":
                get_optional(
                    best,
                    "chrom",
                ),

            "prediction_start0":
                int(
                    best[
                        "start0"
                    ]
                ),

            "prediction_end0":
                int(
                    best[
                        "end0"
                    ]
                ),

            "prediction_length_bp":
                int(
                    best[
                        "prediction_length_bp"
                    ]
                ),

            "overlap_bp":
                int(
                    best[
                        "overlap_bp"
                    ]
                ),

            "overlap_fraction_cre":
                float(
                    best[
                        "overlap_fraction_cre"
                    ]
                ),

            "overlap_fraction_prediction":
                float(
                    best[
                        "overlap_fraction_prediction"
                    ]
                ),

            "flanking_gene":
                get_optional(
                    best,
                    "flanking_gene",
                ),

            "dmel_ortholog_flanking_gene":
                get_optional(
                    best,
                    "dmel_ortholog_flanking_gene",
                ),

            "dmel_fbgn_flanking_gene":
                get_optional(
                    best,
                    "dmel_fbgn_flanking_gene",
                ),

            "distance_flanking_gene":
                get_optional(
                    best,
                    "distance_flanking_gene",
                ),

            "next_flanking_gene":
                get_optional(
                    best,
                    "next_flanking_gene",
                ),

            "dmel_ortholog_next_flanking_gene":
                get_optional(
                    best,
                    "dmel_ortholog_next_flanking_gene",
                ),

            "dmel_fbgn_next_flanking_gene":
                get_optional(
                    best,
                    "dmel_fbgn_next_flanking_gene",
                ),

            "distance_next_gene":
                get_optional(
                    best,
                    "distance_next_gene",
                ),

            "gene_distance_qc":
                gene_distance_qc,
        })


    # ========================================================
    # Build result table
    # ========================================================

    out = pd.DataFrame(
        output_rows
    )


    if not out.empty:

        # ----------------------------------------------------
        # One result per candidate CRE
        # ----------------------------------------------------

        if out[
            "dmel_cre_id"
        ].duplicated().any():

            raise SystemExit(
                "ERROR: duplicate dmel_cre_id values "
                "in gene-distance output."
            )


        # ----------------------------------------------------
        # Stable sorting
        # ----------------------------------------------------

        downstream_order = {
            "high": 1,
            "medium": 2,
            "exploratory": 3,
        }

        candidate_order = {
            "focal_recurrent": 1,
            "focal_plus_secondary_recurrent": 2,
            "secondary_recurrent": 3,
            "focal_single": 4,
            "secondary_single": 5,
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
            "_candidate_rank"
        ] = (
            out[
                "candidate_priority"
            ]
            .map(
                candidate_order
            )
            .fillna(
                99
            )
        )


        if (
            "n_total_tier1_clades"
            in out.columns
        ):

            out[
                "_n_clades"
            ] = pd.to_numeric(
                out[
                    "n_total_tier1_clades"
                ],
                errors="coerce",
            ).fillna(
                0
            )

        else:

            out[
                "_n_clades"
            ] = 0


        out = (
            out
            .sort_values(
                [
                    "_downstream_rank",
                    "_candidate_rank",
                    "_n_clades",
                    "dmel_cre_id",
                ],
                ascending=[
                    True,
                    True,
                    False,
                    True,
                ],
                kind="mergesort",
            )
            .drop(
                columns=[
                    "_downstream_rank",
                    "_candidate_rank",
                    "_n_clades",
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
    # Run metadata
    # ========================================================

    n_matched = (
        int(
            (
                out[
                    "matched_dmel_prediction"
                ]
                == "YES"
            ).sum()
        )
        if not out.empty
        else 0
    )


    n_exact = (
        int(
            (
                out.get(
                    "exact_coordinate_match",
                    pd.Series(
                        dtype=str
                    ),
                )
                == "yes"
            ).sum()
        )
        if not out.empty
        else 0
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

            "candidate_input":
                str(
                    args.candidates.resolve()
                ),

            "candidate_input_sha256":
                file_sha256(
                    args.candidates
                ),

            "so_input":
                str(
                    args.so_table.resolve()
                ),

            "so_input_sha256":
                file_sha256(
                    args.so_table
                ),

            "n_candidate_cres":
                len(
                    candidates
                ),

            "n_dmel_so_rows_initial":
                int(
                    dmel_mask.sum()
                ),

            "n_invalid_dmel_so_coordinates":
                n_invalid_so_coords,

            "n_dmel_so_rows_used":
                len(
                    dmel
                ),

            "n_matched_candidates":
                n_matched,

            "n_unmatched_candidates":
                len(
                    candidates
                )
                - n_matched,

            "n_exact_coordinate_matches":
                n_exact,

            "n_high_priority":
                int(
                    (
                        candidates[
                            "downstream_priority"
                        ]
                        == "high"
                    ).sum()
                ),

            "n_medium_priority":
                int(
                    (
                        candidates[
                            "downstream_priority"
                        ]
                        == "medium"
                    ).sum()
                ),

            "n_exploratory_priority":
                int(
                    (
                        candidates[
                            "downstream_priority"
                        ]
                        == "exploratory"
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
    # Console report
    # ========================================================

    print()
    print("=" * 72)
    print("Candidate gene-distance annotation complete")
    print("=" * 72)

    print(
        f"Candidate CREs: "
        f"{len(candidates)}"
    )

    print(
        f"Dmel SO rows used: "
        f"{len(dmel)}"
    )

    print(
        f"Matched candidates: "
        f"{n_matched}/{len(candidates)}"
    )

    print(
        f"Exact coordinate matches: "
        f"{n_exact}/{len(candidates)}"
    )


    print()
    print(
        "Downstream priorities:"
    )

    print(
        candidates[
            "downstream_priority"
        ]
        .value_counts()
        .to_string()
    )


    if not out.empty:

        print()
        print(
            "Gene-distance QC:"
        )

        print(
            out[
                "gene_distance_qc"
            ]
            .value_counts()
            .to_string()
        )


        # ----------------------------------------------------
        # Print non-PASS rows for immediate inspection
        # ----------------------------------------------------

        flagged = out.loc[
            out[
                "gene_distance_qc"
            ]
            != "PASS"
        ]

        if not flagged.empty:

            print()
            print(
                "Candidates requiring inspection:"
            )

            print(
                flagged[
                    [
                        "dmel_cre_id",
                        "candidate_priority",
                        "matched_dmel_prediction",
                        "n_overlapping_predictions",
                        "overlap_bp",
                        "overlap_fraction_cre",
                        "overlap_fraction_prediction",
                        "gene_distance_qc",
                    ]
                ].to_string(
                    index=False
                )
            )


    print()
    print(
        f"Wrote annotation:\n"
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
