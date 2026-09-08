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


def parse_args():
    p = argparse.ArgumentParser(
        description="Add sensitivity robustness to the primary prioritized Tier-1 candidate table."
    )
    p.add_argument("--candidates", type=Path, required=True)
    p.add_argument("--sensitivity", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--metadata-out", type=Path, required=True)
    p.add_argument("--moderate-robustness-min", type=float, required=True)
    return p.parse_args()


def require_file(path):
    if not path.is_file():
        raise SystemExit(f"ERROR: required file not found:\n{path}")


def require_columns(df, required, label):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(f"ERROR: {label} missing columns:\n" + "\n".join(missing))


def robustness_class(percent_retained, moderate_min):
    if percent_retained == 100:
        return "robust"
    if percent_retained >= moderate_min:
        return "moderate"
    return "sensitive"


def main():
    args = parse_args()
    if not 0 <= args.moderate_robustness_min < 100:
        raise SystemExit("ERROR: moderate-robustness-min must be between 0 and <100.")

    require_file(args.candidates)
    require_file(args.sensitivity)

    candidates = pd.read_csv(args.candidates, sep="\t", dtype=str, keep_default_na=False)
    sensitivity = pd.read_csv(args.sensitivity, sep="\t", dtype=str, keep_default_na=False)

    require_columns(candidates, {"dmel_cre_id", "candidate_priority"}, "prioritized candidate table")
    require_columns(
        sensitivity,
        {
            "dmel_cre_id",
            "primary_candidate_priority",
            "percent_candidate_retained",
            "robust_candidate_all_scenarios",
            "percent_same_priority_as_primary",
            "robust_priority_all_scenarios",
            "lost_in_scenarios",
            "changed_priority_in_scenarios",
        },
        "candidate sensitivity summary",
    )

    if candidates.empty:
        raise SystemExit("ERROR: prioritized candidate table is empty.")
    if candidates["dmel_cre_id"].duplicated().any() or sensitivity["dmel_cre_id"].duplicated().any():
        raise SystemExit("ERROR: duplicate CRE IDs in candidate or sensitivity table.")

    candidate_ids = set(candidates["dmel_cre_id"])
    sensitivity_ids = set(sensitivity["dmel_cre_id"])
    if candidate_ids != sensitivity_ids:
        missing = sorted(candidate_ids - sensitivity_ids)
        extra = sorted(sensitivity_ids - candidate_ids)
        raise SystemExit(
            "ERROR: candidate and sensitivity CRE sets differ.\n"
            + ("Missing sensitivity: " + ", ".join(missing) + "\n" if missing else "")
            + ("Unexpected sensitivity: " + ", ".join(extra) if extra else "")
        )

    # The sensitivity summary must be anchored to the original candidate priority.
    priority_check = candidates[["dmel_cre_id", "candidate_priority"]].merge(
        sensitivity[["dmel_cre_id", "primary_candidate_priority"]],
        on="dmel_cre_id",
        validate="one_to_one",
    )
    mismatch = priority_check["candidate_priority"] != priority_check["primary_candidate_priority"]
    if mismatch.any():
        raise SystemExit(
            "ERROR: primary candidate priority differs between candidate and sensitivity tables.\n\n"
            + priority_check.loc[mismatch].to_string(index=False)
        )

    for column in ["percent_candidate_retained", "percent_same_priority_as_primary"]:
        sensitivity[column] = pd.to_numeric(sensitivity[column], errors="raise")
        if ((sensitivity[column] < 0) | (sensitivity[column] > 100)).any():
            raise SystemExit(f"ERROR: {column} contains values outside 0-100.")

    sensitivity["robustness_class"] = sensitivity["percent_candidate_retained"].apply(
        lambda value: robustness_class(value, args.moderate_robustness_min)
    )
    sensitivity["priority_stability"] = sensitivity["percent_same_priority_as_primary"].apply(
        lambda value: "stable" if value == 100 else "variable"
    )

    keep = sensitivity[[
        "dmel_cre_id",
        "percent_candidate_retained",
        "robust_candidate_all_scenarios",
        "robustness_class",
        "percent_same_priority_as_primary",
        "robust_priority_all_scenarios",
        "priority_stability",
        "lost_in_scenarios",
        "changed_priority_in_scenarios",
    ]].copy()

    out = candidates.merge(keep, on="dmel_cre_id", how="left", validate="one_to_one")
    if out["percent_candidate_retained"].isna().any():
        raise SystemExit("ERROR: sensitivity information missing after candidate merge.")

    robustness_rank = {"robust": 1, "moderate": 2, "sensitive": 3}
    out["_robustness_rank"] = out["robustness_class"].map(robustness_rank)
    if "priority_rank" in out.columns:
        sort_cols = ["priority_rank", "_robustness_rank", "percent_candidate_retained", "dmel_cre_id"]
        ascending = [True, True, False, True]
    else:
        sort_cols = ["_robustness_rank", "percent_candidate_retained", "dmel_cre_id"]
        ascending = [True, False, True]
    out = out.sort_values(sort_cols, ascending=ascending, kind="mergesort").drop(columns="_robustness_rank").reset_index(drop=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    counts = out["robustness_class"].value_counts()
    metadata = pd.DataFrame([{
        "script": Path(__file__).name,
        "run_timestamp": datetime.now().astimezone().isoformat(),
        "n_candidates": len(out),
        "moderate_robustness_min": args.moderate_robustness_min,
        "n_robust": int(counts.get("robust", 0)),
        "n_moderate": int(counts.get("moderate", 0)),
        "n_sensitive": int(counts.get("sensitive", 0)),
        "n_priority_stable": int(out["priority_stability"].eq("stable").sum()),
    }])
    metadata.to_csv(args.metadata_out, sep="\t", index=False)

    print("\n" + "=" * 72)
    print("CANDIDATE SENSITIVITY ANNOTATION")
    print("=" * 72)
    print(f"Candidates: {len(out)}")
    print(out["robustness_class"].value_counts().to_string())
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
