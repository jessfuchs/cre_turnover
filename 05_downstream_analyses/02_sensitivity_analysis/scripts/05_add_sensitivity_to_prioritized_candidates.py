#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT_DIR = (
    Path.home()
    / "cre_turnover"
    / "project"
)

CANDIDATE_TABLE = (
    PROJECT_DIR
    / "downstream_analyses"
    / "candidate_analysis"
    / "results"
    / "tables"
    / "tier1_candidates_prioritized.tsv"
)

SENSITIVITY_TABLE = (
    PROJECT_DIR
    / "downstream_analyses"
    / "sensitivity_analysis"
    / "results"
    / "candidate_sensitivity_summary.tsv"
)

OUT = (
    PROJECT_DIR
    / "downstream_analyses"
    / "candidate_analysis"
    / "results"
    / "tables"
    / "tier1_candidates_prioritized_with_sensitivity.tsv"
)


def main():

    candidates = pd.read_csv(
        CANDIDATE_TABLE,
        sep="\t",
        dtype=str,
    ).fillna("")

    sensitivity = pd.read_csv(
        SENSITIVITY_TABLE,
        sep="\t",
        dtype=str,
    ).fillna("")


    # --------------------------------------------------------
    # Checks
    # --------------------------------------------------------

    if "dmel_cre_id" not in candidates.columns:
        raise SystemExit(
            "ERROR: candidate table lacks dmel_cre_id"
        )

    if "dmel_cre_id" not in sensitivity.columns:
        raise SystemExit(
            "ERROR: sensitivity table lacks dmel_cre_id"
        )

    if sensitivity["dmel_cre_id"].duplicated().any():
        raise SystemExit(
            "ERROR: duplicate dmel_cre_id values "
            "in sensitivity summary."
        )


    # --------------------------------------------------------
    # Convert percentages to numeric
    # --------------------------------------------------------

    sensitivity[
        "percent_candidate_retained"
    ] = pd.to_numeric(
        sensitivity[
            "percent_candidate_retained"
        ],
        errors="raise",
    )

    sensitivity[
        "percent_same_priority_as_primary"
    ] = pd.to_numeric(
        sensitivity[
            "percent_same_priority_as_primary"
        ],
        errors="raise",
    )


    # --------------------------------------------------------
    # Robustness classes
    # --------------------------------------------------------

    def robustness_class(value):

        if value == 100:
            return "robust"

        if value >= 66.6:
            return "moderate"

        return "sensitive"


    sensitivity[
        "robustness_class"
    ] = sensitivity[
        "percent_candidate_retained"
    ].apply(
        robustness_class
    )


    sensitivity[
        "priority_stability"
    ] = sensitivity[
        "percent_same_priority_as_primary"
    ].apply(
        lambda value:
            "stable"
            if value == 100
            else "variable"
    )


    # --------------------------------------------------------
    # Keep only useful sensitivity columns
    # --------------------------------------------------------

    sensitivity_keep = sensitivity[
        [
            "dmel_cre_id",

            "percent_candidate_retained",
            "robust_candidate_all_9",
            "robustness_class",

            "percent_same_priority_as_primary",
            "robust_priority_all_9",
            "priority_stability",

            "lost_in_scenarios",
            "changed_priority_in_scenarios",
        ]
    ].copy()


    # --------------------------------------------------------
    # Merge
    # --------------------------------------------------------

    out = candidates.merge(
        sensitivity_keep,
        on="dmel_cre_id",
        how="left",
        validate="one_to_one",
    )


    # --------------------------------------------------------
    # QC
    # --------------------------------------------------------

    missing = out[
        "percent_candidate_retained"
    ].isna()

    if missing.any():

        print(
            "ERROR: candidate(s) missing "
            "sensitivity information:"
        )

        print(
            out.loc[
                missing,
                "dmel_cre_id",
            ].to_string(
                index=False
            )
        )

        raise SystemExit(1)


    # --------------------------------------------------------
    # Optional robustness rank for sorting
    # --------------------------------------------------------

    robustness_rank = {
        "robust": 1,
        "moderate": 2,
        "sensitive": 3,
    }

    out[
        "_robustness_rank"
    ] = out[
        "robustness_class"
    ].map(
        robustness_rank
    )


    # Keep your original candidate ranking first,
    # then use robustness as additional evidence.
    if "priority_rank" in out.columns:

        sort_cols = [
            "priority_rank",
            "_robustness_rank",
            "percent_candidate_retained",
            "dmel_cre_id",
        ]

        ascending = [
            True,
            True,
            False,
            True,
        ]

    else:

        sort_cols = [
            "_robustness_rank",
            "percent_candidate_retained",
            "dmel_cre_id",
        ]

        ascending = [
            True,
            False,
            True,
        ]


    out = (
        out
        .sort_values(
            sort_cols,
            ascending=ascending,
            kind="mergesort",
        )
        .drop(
            columns="_robustness_rank"
        )
        .reset_index(
            drop=True
        )
    )


    # --------------------------------------------------------
    # Write
    # --------------------------------------------------------

    OUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    out.to_csv(
        OUT,
        sep="\t",
        index=False,
    )


    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------

    print()
    print("=" * 72)
    print("Candidate sensitivity added to prioritized table")
    print("=" * 72)

    print(
        f"Candidates: "
        f"{len(out)}"
    )

    print()

    print(
        "Robustness classes:"
    )

    print(
        out[
            "robustness_class"
        ]
        .value_counts()
        .to_string()
    )

    print()

    print(
        "Priority stability:"
    )

    print(
        out[
            "priority_stability"
        ]
        .value_counts()
        .to_string()
    )

    print()

    print(
        f"Wrote:\n{OUT}"
    )


if __name__ == "__main__":
    main()
