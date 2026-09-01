#!/usr/bin/env python3

# ============================================================
# Export Azteca focal Tier-1 candidates
#
# Purpose:
#   Export a compact thesis-ready table of prioritized CREs
#   showing focal Tier-1 support in the Azteca clade.
#
# Candidate robustness is taken from the sensitivity-annotated
# candidate table. The Azteca-specific state direction is taken
# from the focal Tier-1 analysis rather than inferred from
# aggregated cross-clade annotations.
#
# Input/output paths and the focal-clade name are supplied by
# the climate-analysis configuration via the wrapper.
# ============================================================

from pathlib import Path
import argparse

import pandas as pd


# ============================================================
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Export thesis-ready Azteca focal Tier-1 CRE candidates."
    )
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--focal-events", type=Path, required=True)
    parser.add_argument("--clade", required=True)
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_file(path):
    if not path.is_file():
        raise SystemExit(f"ERROR: required file not found:\n{path}")


def require_columns(df, required, label):
    missing = sorted(set(required) - set(df.columns))

    if missing:
        raise SystemExit(
            f"ERROR: {label} missing columns:\n" + "\n".join(missing)
        )


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    require_file(args.candidates)
    require_file(args.focal_events)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    # ========================================================
    # Load prioritized candidates
    # ========================================================

    candidates = pd.read_csv(
        args.candidates,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    require_columns(
        candidates,
        {
            "dmel_cre_id",
            "fbgn_target_genes",
            "n_fbgn_target_genes",
            "n_total_tier1_clades",
            "tier1_clades",
            "percent_candidate_retained",
            "robust_candidate_all_9",
            "robustness_class",
            "percent_same_priority_as_primary",
            "robust_priority_all_9",
            "priority_stability",
        },
        "sensitivity-annotated candidate table",
    )

    if candidates["dmel_cre_id"].duplicated().any():
        raise SystemExit("ERROR: duplicate CRE IDs in candidate table.")

    # ========================================================
    # Load focal Tier-1 events
    # ========================================================

    focal = pd.read_csv(
        args.focal_events,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    require_columns(
        focal,
        {
            "clade",
            "dmel_cre_id",
            "is_focal_tier1",
            "focal_direction",
            "discordant_state",
            "consensus_state",
        },
        "focal Tier-1 event table",
    )

    # Keep only true focal Tier-1 events from the selected clade.
    focal = focal.loc[
        (focal["clade"] == args.clade)
        & (focal["is_focal_tier1"].astype(str).str.lower().isin(["true", "yes", "1"]))
    ].copy()

    if focal.empty:
        raise SystemExit(
            f"ERROR: no focal Tier-1 events found for {args.clade}."
        )

    if focal["dmel_cre_id"].duplicated().any():
        raise SystemExit(
            f"ERROR: duplicate CRE IDs among focal events for {args.clade}."
        )

    # ========================================================
    # Merge focal direction with prioritized candidates
    # ========================================================

    focal = focal[
        [
            "dmel_cre_id",
            "focal_direction",
            "discordant_state",
            "consensus_state",
        ]
    ]

    azteca = candidates.merge(
        focal,
        on="dmel_cre_id",
        how="inner",
        validate="one_to_one",
    )

    if azteca.empty:
        raise SystemExit(
            "ERROR: no Azteca focal Tier-1 events overlap "
            "the prioritized candidate table."
        )

    unexpected_directions = sorted(
        set(azteca["focal_direction"])
        - {"focal_turnover", "focal_present"}
    )

    if unexpected_directions:
        raise SystemExit(
            "ERROR: unexpected focal-direction values:\n"
            + "\n".join(unexpected_directions)
        )

    # ========================================================
    # Biological interpretation
    # ========================================================

    direction_labels = {
        "focal_turnover":
            "Azteca turnover candidate; comparison species present",
        "focal_present":
            "Azteca present; comparison species turnover candidates",
    }

    azteca["azteca_pattern"] = (
        azteca["focal_direction"].map(direction_labels)
    )

    azteca = azteca.rename(
        columns={"focal_direction": "azteca_focal_direction"}
    )

    # ========================================================
    # Thesis-ready output
    # ========================================================

    output = azteca[
        [
            "dmel_cre_id",
            "azteca_focal_direction",
            "azteca_pattern",
            "fbgn_target_genes",
            "n_fbgn_target_genes",
            "n_total_tier1_clades",
            "tier1_clades",
            "percent_candidate_retained",
            "robust_candidate_all_9",
            "robustness_class",
            "percent_same_priority_as_primary",
            "robust_priority_all_9",
            "priority_stability",
        ]
    ].copy()

    for column in [
        "n_fbgn_target_genes",
        "n_total_tier1_clades",
        "percent_candidate_retained",
        "percent_same_priority_as_primary",
    ]:
        output[column] = pd.to_numeric(output[column], errors="raise")

    # ========================================================
    # Candidate ordering
    # ========================================================

    # Robustness and priority stability are used only to order
    # the thesis table, not to redefine candidate status.
    robustness_rank = {
        "robust": 1,
        "moderate": 2,
        "sensitive": 3,
    }

    stability_rank = {
        "stable": 1,
        "variable": 2,
    }

    output["_robust_rank"] = output["robustness_class"].map(robustness_rank)
    output["_stability_rank"] = output["priority_stability"].map(stability_rank)

    output = (
        output
        .sort_values(
            [
                "_robust_rank",
                "_stability_rank",
                "n_total_tier1_clades",
                "dmel_cre_id",
            ],
            ascending=[True, True, False, True],
            kind="mergesort",
        )
        .drop(columns=["_robust_rank", "_stability_rank"])
        .reset_index(drop=True)
    )

    # ========================================================
    # Write
    # ========================================================

    output.to_csv(args.out, sep="\t", index=False)

    print()
    print("=" * 72)
    print("AZTECA FOCAL TIER-1 CANDIDATES")
    print("=" * 72)

    display_columns = [
        "dmel_cre_id",
        "azteca_focal_direction",
        "fbgn_target_genes",
        "n_total_tier1_clades",
        "percent_candidate_retained",
        "robustness_class",
        "priority_stability",
    ]

    print()
    print(output[display_columns].to_string(index=False))

    print()
    print(f"Number of candidates: {len(output)}")
    print()
    print(f"Wrote:\n{args.out}")


if __name__ == "__main__":
    main()
