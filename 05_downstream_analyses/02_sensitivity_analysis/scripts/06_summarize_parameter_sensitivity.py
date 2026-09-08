#!/usr/bin/env python3

# ============================================================
# Summarize parameter-specific Tier-1 candidate sensitivity
#
# Purpose:
#   Separate the effects of the two CRE-state classification
#   parameters on Tier-1 candidate retention:
#
#   1. reciprocal-overlap threshold
#      varied at fixed local distance = 24 kb
#
#   2. local-distance threshold
#      varied at fixed reciprocal overlap = 0.50
#
# Input:
#   candidate_sensitivity_retention_matrix.tsv
#   sensitivity_global_summary.tsv
#
# Output:
#   candidate_parameter_sensitivity.tsv
#   candidate_parameter_sensitivity_summary.tsv
#   candidate_parameter_sensitivity_class_summary.tsv
#
# This script does not rerun classification.
# ============================================================


from pathlib import Path
import argparse
import pandas as pd


# ============================================================
# Scenario definitions
# ============================================================

PRIMARY_SCENARIO = "ov050_dist24000"

OVERLAP_SCENARIOS = [
    "ov025_dist24000",
    "ov050_dist24000",
    "ov075_dist24000",
]

OVERLAP_NONPRIMARY = [
    "ov025_dist24000",
    "ov075_dist24000",
]

DISTANCE_SCENARIOS = [
    "ov050_dist12000",
    "ov050_dist24000",
    "ov050_dist48000",
]

DISTANCE_NONPRIMARY = [
    "ov050_dist12000",
    "ov050_dist48000",
]

VALID_RETENTION_VALUES = {"yes", "no"}


# ============================================================
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize parameter-specific sensitivity of Tier-1 CRE candidate retention."
    )
    parser.add_argument("--retention-matrix", type=Path, required=True)
    parser.add_argument("--global-summary", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_columns(df, required, table_name):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(
            f"ERROR: {table_name} missing columns:\n" + "\n".join(missing)
        )


def normalize_retention(series):
    values = series.astype(str).str.strip().str.lower()
    unknown = sorted(set(values) - VALID_RETENTION_VALUES)

    if unknown:
        raise SystemExit(
            "ERROR: unexpected candidate-retention values:\n" + "\n".join(unknown)
        )

    return values == "yes"


def get_lost_scenarios(row, scenarios):
    lost = [scenario for scenario in scenarios if not bool(row[scenario])]
    return "|".join(lost) if lost else "NA"


def classify_parameter_sensitivity(overlap_retained, distance_retained):
    overlap_sensitive = overlap_retained < 3
    distance_sensitive = distance_retained < 3

    if not overlap_sensitive and not distance_sensitive:
        return "robust_to_both"
    if overlap_sensitive and not distance_sensitive:
        return "overlap_sensitive_only"
    if not overlap_sensitive and distance_sensitive:
        return "distance_sensitive_only"
    return "sensitive_to_both"


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    if not args.retention_matrix.is_file():
        raise SystemExit(f"ERROR: retention matrix not found:\n{args.retention_matrix}")

    if not args.global_summary.is_file():
        raise SystemExit(f"ERROR: global sensitivity summary not found:\n{args.global_summary}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Candidate-retention matrix
    # --------------------------------------------------------

    retention = pd.read_csv(
        args.retention_matrix,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required_retention = {
        "dmel_cre_id",
        *OVERLAP_SCENARIOS,
        *DISTANCE_SCENARIOS,
    }

    require_columns(retention, required_retention, "candidate retention matrix")

    if retention.empty:
        raise SystemExit("ERROR: candidate retention matrix is empty.")

    if retention["dmel_cre_id"].duplicated().any():
        duplicates = (
            retention.loc[
                retention["dmel_cre_id"].duplicated(keep=False),
                "dmel_cre_id",
            ]
            .drop_duplicates()
            .tolist()
        )
        raise SystemExit(
            "ERROR: duplicate CRE IDs in retention matrix:\n" + "\n".join(duplicates)
        )

    scenarios_needed = sorted(set(OVERLAP_SCENARIOS + DISTANCE_SCENARIOS))

    bool_matrix = retention[["dmel_cre_id", *scenarios_needed]].copy()

    for scenario in scenarios_needed:
        bool_matrix[scenario] = normalize_retention(bool_matrix[scenario])

    if not bool_matrix[PRIMARY_SCENARIO].all():
        failed = bool_matrix.loc[
            ~bool_matrix[PRIMARY_SCENARIO],
            "dmel_cre_id",
        ].tolist()
        raise SystemExit(
            "ERROR: primary Tier-1 candidates were not retained in the primary scenario:\n"
            + "\n".join(failed)
        )

    # --------------------------------------------------------
    # Candidate-level parameter sensitivity
    # --------------------------------------------------------

    rows = []

    for _, row in bool_matrix.iterrows():
        overlap_retained_n = int(row[OVERLAP_SCENARIOS].sum())
        distance_retained_n = int(row[DISTANCE_SCENARIOS].sum())

        overlap_retention_pct = 100.0 * overlap_retained_n / len(OVERLAP_SCENARIOS)
        distance_retention_pct = 100.0 * distance_retained_n / len(DISTANCE_SCENARIOS)

        rows.append({
            "dmel_cre_id": row["dmel_cre_id"],
            "overlap_retained_n": overlap_retained_n,
            "overlap_retention_pct": overlap_retention_pct,
            "overlap_sensitive": "yes" if overlap_retained_n < 3 else "no",
            "overlap_lost_scenarios": get_lost_scenarios(row, OVERLAP_NONPRIMARY),
            "distance_retained_n": distance_retained_n,
            "distance_retention_pct": distance_retention_pct,
            "distance_sensitive": "yes" if distance_retained_n < 3 else "no",
            "distance_lost_scenarios": get_lost_scenarios(row, DISTANCE_NONPRIMARY),
            "parameter_sensitivity_class": classify_parameter_sensitivity(
                overlap_retained_n,
                distance_retained_n,
            ),
        })

    candidate = pd.DataFrame(rows)

    class_order = {
        "robust_to_both": 0,
        "overlap_sensitive_only": 1,
        "distance_sensitive_only": 2,
        "sensitive_to_both": 3,
    }

    candidate["_class_order"] = candidate["parameter_sensitivity_class"].map(class_order)
    candidate = (
        candidate
        .sort_values(["_class_order", "dmel_cre_id"], kind="mergesort")
        .drop(columns="_class_order")
        .reset_index(drop=True)
    )

    # --------------------------------------------------------
    # Candidate losses per non-primary scenario
    # --------------------------------------------------------

    candidate_losses = {
        scenario: int((~bool_matrix[scenario]).sum())
        for scenario in OVERLAP_NONPRIMARY + DISTANCE_NONPRIMARY
    }

    # --------------------------------------------------------
    # Global CRE-state sensitivity
    # --------------------------------------------------------

    global_summary = pd.read_csv(args.global_summary, sep="\t")
    require_columns(global_summary, {"scenario", "n_changed"}, "global sensitivity summary")

    global_summary["n_changed"] = pd.to_numeric(
        global_summary["n_changed"],
        errors="raise",
    )

    global_by_scenario = global_summary.set_index("scenario")["n_changed"].to_dict()

    required_global = set(OVERLAP_SCENARIOS + DISTANCE_SCENARIOS)
    missing_global = sorted(required_global - set(global_by_scenario))

    if missing_global:
        raise SystemExit(
            "ERROR: global summary missing scenarios:\n" + "\n".join(missing_global)
        )

    if int(global_by_scenario[PRIMARY_SCENARIO]) != 0:
        raise SystemExit(
            "ERROR: primary global sensitivity scenario does not contain zero changed states."
        )

    # --------------------------------------------------------
    # Parameter-level summary
    # --------------------------------------------------------

    summary_rows = []

    parameter_defs = [
        {
            "parameter_dimension": "reciprocal_overlap",
            "fixed_parameter": "local_distance=24000",
            "scenarios": OVERLAP_SCENARIOS,
            "nonprimary": OVERLAP_NONPRIMARY,
            "retention_col": "overlap_retention_pct",
            "retained_n_col": "overlap_retained_n",
        },
        {
            "parameter_dimension": "local_distance",
            "fixed_parameter": "reciprocal_overlap=0.50",
            "scenarios": DISTANCE_SCENARIOS,
            "nonprimary": DISTANCE_NONPRIMARY,
            "retention_col": "distance_retention_pct",
            "retained_n_col": "distance_retained_n",
        },
    ]

    for definition in parameter_defs:
        retained_n = candidate[definition["retained_n_col"]]
        retention_pct = candidate[definition["retention_col"]]

        nonprimary_candidate_losses = [
            candidate_losses[s] for s in definition["nonprimary"]
        ]
        nonprimary_global_changes = [
            int(global_by_scenario[s]) for s in definition["nonprimary"]
        ]

        summary_rows.append({
            "parameter_dimension": definition["parameter_dimension"],
            "fixed_parameter": definition["fixed_parameter"],
            "scenarios": "|".join(definition["scenarios"]),
            "n_candidates": len(candidate),
            "n_retained_3_of_3": int((retained_n == 3).sum()),
            "pct_retained_3_of_3": 100.0 * (retained_n == 3).mean(),
            "n_retained_2_of_3": int((retained_n == 2).sum()),
            "pct_retained_2_of_3": 100.0 * (retained_n == 2).mean(),
            "n_retained_1_of_3": int((retained_n == 1).sum()),
            "pct_retained_1_of_3": 100.0 * (retained_n == 1).mean(),
            "median_retention_pct": retention_pct.median(),
            "mean_retention_pct": retention_pct.mean(),
            "n_sensitive_to_parameter": int((retained_n < 3).sum()),
            "pct_sensitive_to_parameter": 100.0 * (retained_n < 3).mean(),
            "candidate_losses_lower_setting": nonprimary_candidate_losses[0],
            "candidate_losses_upper_setting": nonprimary_candidate_losses[1],
            "mean_candidate_losses_nonprimary": sum(nonprimary_candidate_losses) / 2,
            "median_candidate_losses_nonprimary": pd.Series(
                nonprimary_candidate_losses
            ).median(),
            "global_state_changes_lower_setting": nonprimary_global_changes[0],
            "global_state_changes_upper_setting": nonprimary_global_changes[1],
            "mean_global_state_changes_nonprimary": sum(nonprimary_global_changes) / 2,
            "median_global_state_changes_nonprimary": pd.Series(
                nonprimary_global_changes
            ).median(),
        })

    parameter_summary = pd.DataFrame(summary_rows)

    # --------------------------------------------------------
    # Sensitivity-class summary
    # --------------------------------------------------------

    class_summary = (
        candidate["parameter_sensitivity_class"]
        .value_counts()
        .rename_axis("parameter_sensitivity_class")
        .reset_index(name="n_candidates")
    )

    all_classes = pd.DataFrame({
        "parameter_sensitivity_class": list(class_order.keys())
    })

    class_summary = all_classes.merge(
        class_summary,
        on="parameter_sensitivity_class",
        how="left",
    )

    class_summary["n_candidates"] = class_summary["n_candidates"].fillna(0).astype(int)
    class_summary["percent_candidates"] = (
        100.0 * class_summary["n_candidates"] / len(candidate)
    )

    class_summary["_order"] = class_summary["parameter_sensitivity_class"].map(class_order)
    class_summary = (
        class_summary
        .sort_values("_order")
        .drop(columns="_order")
        .reset_index(drop=True)
    )

    # --------------------------------------------------------
    # Round
    # --------------------------------------------------------

    candidate["overlap_retention_pct"] = candidate["overlap_retention_pct"].round(1)
    candidate["distance_retention_pct"] = candidate["distance_retention_pct"].round(1)

    for col in [
        "pct_retained_3_of_3",
        "pct_retained_2_of_3",
        "pct_retained_1_of_3",
        "median_retention_pct",
        "mean_retention_pct",
        "pct_sensitive_to_parameter",
        "mean_candidate_losses_nonprimary",
        "median_candidate_losses_nonprimary",
        "mean_global_state_changes_nonprimary",
        "median_global_state_changes_nonprimary",
    ]:
        parameter_summary[col] = parameter_summary[col].round(1)

    class_summary["percent_candidates"] = class_summary[
        "percent_candidates"
    ].round(1)

    # --------------------------------------------------------
    # Write
    # --------------------------------------------------------

    out_candidates = args.out_dir / "candidate_parameter_sensitivity.tsv"
    out_summary = args.out_dir / "candidate_parameter_sensitivity_summary.tsv"
    out_classes = args.out_dir / "candidate_parameter_sensitivity_class_summary.tsv"

    candidate.to_csv(out_candidates, sep="\t", index=False)
    parameter_summary.to_csv(out_summary, sep="\t", index=False)
    class_summary.to_csv(out_classes, sep="\t", index=False)

    # --------------------------------------------------------
    # Console report
    # --------------------------------------------------------

    overlap = parameter_summary.loc[
        parameter_summary["parameter_dimension"] == "reciprocal_overlap"
    ].iloc[0]

    distance = parameter_summary.loc[
        parameter_summary["parameter_dimension"] == "local_distance"
    ].iloc[0]

    print()
    print("=" * 72)
    print("PARAMETER-SPECIFIC TIER-1 CANDIDATE SENSITIVITY")
    print("=" * 72)
    print(f"Candidates evaluated: {len(candidate)}")
    print(f"Primary scenario: {PRIMARY_SCENARIO}")
    print()

    print("OVERLAP THRESHOLD")
    print("-" * 72)
    print(f"Retained in 3/3 settings: {int(overlap['n_retained_3_of_3'])}/{len(candidate)} "
          f"({overlap['pct_retained_3_of_3']:.1f}%)")
    print(f"Retained in 2/3 settings: {int(overlap['n_retained_2_of_3'])}/{len(candidate)} "
          f"({overlap['pct_retained_2_of_3']:.1f}%)")
    print(f"Retained in 1/3 settings: {int(overlap['n_retained_1_of_3'])}/{len(candidate)} "
          f"({overlap['pct_retained_1_of_3']:.1f}%)")
    print(f"Candidate losses at overlap 0.25: {candidate_losses['ov025_dist24000']}")
    print(f"Candidate losses at overlap 0.75: {candidate_losses['ov075_dist24000']}")
    print(f"Global CRE-state changes at overlap 0.25: "
          f"{int(global_by_scenario['ov025_dist24000'])}")
    print(f"Global CRE-state changes at overlap 0.75: "
          f"{int(global_by_scenario['ov075_dist24000'])}")
    print()

    print("LOCAL-DISTANCE THRESHOLD")
    print("-" * 72)
    print(f"Retained in 3/3 settings: {int(distance['n_retained_3_of_3'])}/{len(candidate)} "
          f"({distance['pct_retained_3_of_3']:.1f}%)")
    print(f"Retained in 2/3 settings: {int(distance['n_retained_2_of_3'])}/{len(candidate)} "
          f"({distance['pct_retained_2_of_3']:.1f}%)")
    print(f"Retained in 1/3 settings: {int(distance['n_retained_1_of_3'])}/{len(candidate)} "
          f"({distance['pct_retained_1_of_3']:.1f}%)")
    print(f"Candidate losses at 12 kb: {candidate_losses['ov050_dist12000']}")
    print(f"Candidate losses at 48 kb: {candidate_losses['ov050_dist48000']}")
    print(f"Global CRE-state changes at 12 kb: "
          f"{int(global_by_scenario['ov050_dist12000'])}")
    print(f"Global CRE-state changes at 48 kb: "
          f"{int(global_by_scenario['ov050_dist48000'])}")
    print()

    print("CANDIDATE SENSITIVITY CLASSES")
    print("-" * 72)
    print(class_summary.to_string(index=False))
    print()

    print(f"Wrote:\n{out_candidates}")
    print(f"Wrote:\n{out_summary}")
    print(f"Wrote:\n{out_classes}")


if __name__ == "__main__":
    main()