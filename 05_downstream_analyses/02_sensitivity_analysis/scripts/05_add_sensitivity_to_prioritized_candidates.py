#!/usr/bin/env python3

# ============================================================
# Add sensitivity evidence to prioritized Tier-1 candidates
#
# Purpose:
#   Merge candidate-level sensitivity results with the primary
#   prioritized Tier-1 CRE table.
#
# Added annotations:
#   - candidate retention across sensitivity scenarios
#   - robustness class
#   - candidate-priority stability
#   - scenarios in which candidates are lost or reprioritized
#
# The original candidate priority remains the primary ranking;
# sensitivity robustness is added as supporting evidence.
#
# Input/output paths and robustness threshold are supplied by
# config/sensitivity_config.sh via the pipeline wrapper.
# ============================================================

from datetime import datetime
from pathlib import Path
import argparse

import pandas as pd


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Add candidate-sensitivity evidence to the "
            "prioritized Tier-1 candidate table."
        )
    )

    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--sensitivity", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--metadata-out", type=Path, required=True)
    parser.add_argument( "--moderate-robustness-min", type=float, required=True)
    
    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_file(path):
    """Require an existing regular input file."""
    if not path.is_file():
        raise SystemExit(f"ERROR: required file not found:\n{path}")


def require_columns(df, required, label):
    """Require the columns needed for sensitivity annotation."""
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(f"ERROR: {label} missing columns:\n" + "\n".join(missing))

def robustness_class(
    percent_retained,
    moderate_min,
):
    """
    Classify candidate robustness from sensitivity retention.

    robust:
        retained in all sensitivity scenarios

    moderate:
        retained in at least the configured percentage

    sensitive:
        retained below the moderate threshold
    """
    if percent_retained == 100:
        return "robust"
    if percent_retained >= moderate_min:
        return "moderate"
    return "sensitive"


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    if not 0 <= args.moderate_robustness_min <= 100:
        raise SystemExit("ERROR: moderate-robustness-min must be between 0 and 100.")

    require_file(args.candidates)
    require_file(args.sensitivity)
    
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_out.parent.mkdir(parents=True, exist_ok=True)

    # ========================================================
    # Load inputs
    # ========================================================

    candidates = pd.read_csv(
        args.candidates,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    sensitivity = pd.read_csv(
        args.sensitivity,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    require_columns(
        candidates,
        {
            "dmel_cre_id",
        },
        "prioritized candidate table",
    )

    require_columns(
        sensitivity,
        {
            "dmel_cre_id",
            "percent_candidate_retained",
            "robust_candidate_all_9",
            "percent_same_priority_as_primary",
            "robust_priority_all_9",
            "lost_in_scenarios",
            "changed_priority_in_scenarios",
        },
        "candidate sensitivity summary",
    )


    # ========================================================
    # Input QC
    # ========================================================

    if candidates.empty:
        raise SystemExit("ERROR: prioritized candidate table is empty.")
    if candidates["dmel_cre_id"].duplicated().any():
        raise SystemExit("ERROR: duplicate CRE IDs in prioritized candidate table.")
    if sensitivity["dmel_cre_id"].duplicated().any():
        raise SystemExit("ERROR: duplicate CRE IDs in candidate sensitivity summary.")

    # ========================================================
    # Sensitivity annotations
    # ========================================================

    # Percentages are required numeric evidence and should not
    # silently become missing values.
    for column in [
        "percent_candidate_retained",
        "percent_same_priority_as_primary",
    ]:
        sensitivity[column] = pd.to_numeric(
            sensitivity[column],
            errors="raise",
        )

        if (
            (sensitivity[column] < 0)
            | (sensitivity[column] > 100)
        ).any():
            raise SystemExit(f"ERROR: {column} contains values outside 0-100.")

    # Candidate robustness describes whether the CRE remains
    # a Tier-1 candidate across classification thresholds.
    sensitivity["robustness_class"] = (
        sensitivity[
            "percent_candidate_retained"
        ]
        .apply(
            lambda value:
                robustness_class(
                    value,
                    args.moderate_robustness_min,
                )
        )
    )

    # Priority stability is stricter than candidate retention:
    # the exact candidate-priority category must remain unchanged.
    sensitivity["priority_stability"] = (
        sensitivity[
            "percent_same_priority_as_primary"
        ]
        .apply(
            lambda value:
                "stable"
                if value == 100
                else "variable"
        )
    )


    # ========================================================
    # Keep sensitivity columns relevant downstream
    # ========================================================

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


    # ========================================================
    # Merge with primary candidate table
    # ========================================================

    out = candidates.merge(
        sensitivity_keep,
        on="dmel_cre_id",
        how="left",
        validate="one_to_one",
    )

    # Every primary Tier-1 candidate must have been evaluated
    # by the candidate-sensitivity analysis.
    missing = out[
        "percent_candidate_retained"
    ].isna()

    if missing.any():

        bad = out.loc[
            missing,
            "dmel_cre_id",
        ].tolist()

        raise SystemExit(
            "ERROR: candidate(s) missing sensitivity "
            "information:\n"
            + "\n".join(bad)
        )


    # ========================================================
    # Candidate ordering
    # ========================================================

    # Sensitivity evidence supplements rather than replaces the
    # original candidate-priority ranking.
    robustness_rank = {
        "robust": 1,
        "moderate": 2,
        "sensitive": 3,
    }

    out["_robustness_rank"] = (out["robustness_class"].map(robustness_rank))

    if "priority_rank" in out.columns:
        sort_columns = ["priority_rank", "_robustness_rank", "percent_candidate_retained", "dmel_cre_id"]
        ascending = [True, True, False, True]
    else:
        sort_columns = ["_robustness_rank", "percent_candidate_retained", "dmel_cre_id"]
        ascending = [True, False, True,]

    out = (out.sort_values(sort_columns, ascending=ascending, kind="mergesort")
           .drop(columns=["_robustness_rank"]).reset_index(drop=True))

    # ========================================================
    # Write annotated candidate table
    # ========================================================

    out.to_csv(args.out, sep="\t", index=False)

    # ========================================================
    # Run metadata
    # ========================================================

    robustness_counts = (out["robustness_class"].value_counts())

    metadata = pd.DataFrame([{
        "script": Path(__file__).name,
        "run_timestamp": datetime.now().astimezone().isoformat(),
        "n_candidates": len(out),
        "moderate_robustness_min": args.moderate_robustness_min,
        "n_robust": int(robustness_counts.get("robust", 0)),
        "n_moderate": int(robustness_counts.get("moderate", 0)),
        "n_sensitive": int(robustness_counts.get("sensitive", 0)),
        "n_priority_stable": int((out["priority_stability"] == "stable").sum()),
    }])

    metadata.to_csv(args.metadata_out, sep="\t", index=False)

    # ========================================================
    # Console summary
    # ========================================================

    print()
    print("=" * 72)
    print("Candidate sensitivity added to prioritized table")
    print("=" * 72)
    print(
        f"Candidates: "
        f"{len(out)}"
    )
    print()
    print("Robustness classes:")
    print(out["robustness_class"].value_counts().to_string())
    print()
    print("Priority stability:")
    print(out["priority_stability"].value_counts().to_string())
    print()
    print(
        f"Wrote annotated candidates:\n"
        f"{args.out}"
    )
    print(
        f"Wrote metadata:\n"
        f"{args.metadata_out}"
    )

if __name__ == "__main__":
    main()
