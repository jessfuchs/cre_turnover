#!/usr/bin/env python3

# ============================================================
# Summarize Tier-1 candidate sensitivity
#
# Purpose:
#   Re-evaluate prioritized Tier-1 CRE candidates across all
#   sensitivity scenarios and quantify candidate robustness.
#
# For each candidate and scenario, the script determines:
#   - focal Tier-1 clade support
#   - secondary-only Tier-1 clade support
#   - resulting candidate priority
#   - retention relative to the primary scenario
#   - priority stability relative to the primary scenario
#
# Input/output paths and recurrence thresholds are supplied by
# config/sensitivity_config.sh via the pipeline wrapper.
# ============================================================

from collections import Counter
from datetime import datetime
from pathlib import Path
import argparse

import pandas as pd


# ============================================================
# CRE-state definitions
# ============================================================

POSITIVE_STATES = {
    "present",
    "turnover_candidate",
}

VALID_STATES = {
    "present",
    "turnover_candidate",
    "no_detected_CRE",
    "uncertain",
}


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Summarize Tier-1 candidate robustness across "
            "CRE-classification sensitivity scenarios."
        )
    )

    parser.add_argument("--sensitivity-dir", type=Path, required=True)
    parser.add_argument("--scenario-manifest", type=Path, required=True)
    parser.add_argument("--groups", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--primary-scenario", required=True)

    parser.add_argument("--out-long", type=Path, required=True)
    parser.add_argument("--out-summary", type=Path, required=True)
    parser.add_argument("--out-priority-matrix", type=Path, required=True)
    parser.add_argument("--out-retention-matrix", type=Path, required=True)
    parser.add_argument("--metadata-out", type=Path, required=True)

    parser.add_argument("--focal-recurrent-min-clades", type=int, required=True)
    parser.add_argument("--secondary-recurrent-min-clades", type=int, required=True)
    parser.add_argument("--expected-reference-cres", type=int, required=True)

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
        raise SystemExit(f"ERROR: {label} missing columns:\n" + "\n".join(missing))

def load_scenarios(path, primary_scenario):
    """Read scenario order from the Step-01 scenario manifest."""
    
    require_file(path)
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    require_columns(
        df,
        {"scenario"},
        "scenario manifest",
    )

    scenarios = df["scenario"].tolist()

    if not scenarios:
        raise SystemExit("ERROR: no sensitivity scenarios found.")
    if len(scenarios) != len(set(scenarios)):
        raise SystemExit("ERROR: duplicate scenarios in scenario manifest.")
    if primary_scenario not in scenarios:
        raise SystemExit("ERROR: primary scenario not found in scenario manifest.")

    return scenarios


def load_groups(path):
    """
    Load the same focal-clade definitions used in the primary
    candidate analysis.
    """

    require_file(path)
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    require_columns(
        df,
        {
            "group_name",
            "focal_species",
            "comparison_species",
        },
        "focal-clade table",
    )

    if df["group_name"].duplicated().any():
        raise SystemExit("ERROR: duplicate focal-clade names.")

    groups = {}

    for _, row in df.iterrows():
        focal = row["focal_species"].strip()
        comparisons = [
            x.strip()
            for x in row["comparison_species"].split("|")
            if x.strip()
        ]
        if not focal or not comparisons:
            raise SystemExit(
                f"ERROR: incomplete focal-clade definition "
                f"for {row['group_name']}."
            )
        groups[row["group_name"]] = {
            "focal": focal,
            "comparison": comparisons,
        }

    return groups


def load_states(
    sensitivity_dir,
    scenario,
    species,
    expected_reference_cres,
):
    """
    Load one species-level scenario table as:
        dmel_cre_id -> CRE state
    """

    path = (
        sensitivity_dir
        / scenario
        / f"dmel_to_{species}_cre_turnover.tsv"
    )

    require_file(path)
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    require_columns(
        df,
        {
            "dmel_cre_id",
            "class",
        },
        str(path),
    )

    if df["dmel_cre_id"].duplicated().any():
        raise SystemExit(f"ERROR: duplicate CRE IDs in:\n{path}")
    if len(df) != expected_reference_cres:
        raise SystemExit(
            f"ERROR: unexpected CRE count in:\n{path}\n"
            f"Expected: {expected_reference_cres}\n"
            f"Observed: {len(df)}"
        )

    unknown = sorted(set(df["class"]) - VALID_STATES)
    if unknown:
        raise SystemExit(f"ERROR: unknown CRE states in:\n{path}\n" + "\n".join(unknown))
    return dict(zip(df["dmel_cre_id"], df["class"]))

def is_focal_tier1(
    focal_state,
    comparison_states,
):
    """
    Focal Tier 1 requires unanimous comparison species and a
    present <-> turnover_candidate contrast with the focal species.
    """

    if len(set(comparison_states)) != 1:
        return False

    consensus = comparison_states[0]
    return (focal_state in POSITIVE_STATES and consensus in POSITIVE_STATES and focal_state != consensus)


def is_singleton_tier1(states):
    """
    Singleton Tier 1 requires exactly one discordant species,
    with both observed states being positive CRE states.
    """

    values = list(states.values())

    if any(state not in POSITIVE_STATES for state in values
    ):
        return False

    counts = Counter(values)

    return (len(counts) == 2 and sorted(counts.values()) == [1, len(values) - 1])


def assign_priority(
    n_focal,
    n_secondary,
    focal_recurrent_min,
    secondary_recurrent_min,
):
    """Apply the same priority rules as the primary candidate analysis."""

    if n_focal >= focal_recurrent_min:
        return "focal_recurrent"
    if n_focal >= 1 and n_secondary >= 1:
        return "focal_plus_secondary_recurrent"
    if n_focal == 1:
        return "focal_single"
    if n_secondary >= secondary_recurrent_min:
        return "secondary_recurrent"
    if n_secondary == 1:
        return "secondary_single"
    return "not_tier1_candidate"


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    require_file(args.candidates)
    scenarios = load_scenarios(args.scenario_manifest,args.primary_scenario)
    groups = load_groups(args.groups)

    # --------------------------------------------------------
    # Candidate set from the primary analysis
    # --------------------------------------------------------

    candidates = pd.read_csv(
        args.candidates,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    require_columns(
        candidates,
        {"dmel_cre_id"},
        "prioritized candidate table",
    )

    if candidates.empty:
        raise SystemExit("ERROR: prioritized candidate table is empty.")

    if candidates["dmel_cre_id"].duplicated().any():
        raise SystemExit("ERROR: duplicate CRE IDs in candidate table.")

    candidate_ids = candidates["dmel_cre_id"].tolist()

    # Only species belonging to the focal clades are required.
    species_needed = sorted({
        species
        for group in groups.values()
        for species in [
            group["focal"],
            *group["comparison"],
        ]
    })


    # ========================================================
    # Re-evaluate candidates across sensitivity scenarios
    # ========================================================

    rows = []

    for scenario in scenarios:

        # Load each required species once for this scenario.
        state_tables = {
            species: load_states(
                args.sensitivity_dir,
                scenario,
                species,
                args.expected_reference_cres,
            )
            for species in species_needed
        }

        for cre in candidate_ids:

            focal_groups = []
            secondary_groups = []

            for group_name, group in groups.items():

                focal = group["focal"]
                comparisons = group["comparison"]

                states = {species:state_tables[species].get(cre,"MISSING")
                    for species in [focal, *comparisons]}

                if "MISSING" in states.values():
                    raise SystemExit(
                        f"ERROR: missing state for {cre} "
                        f"in {scenario}."
                    )

                focal_tier1 = is_focal_tier1(
                    states[focal], [states[species] for species in comparisons])

                singleton_tier1 = (is_singleton_tier1(states))

                # Singleton support is secondary only when the
                # discordant species is not the predefined focal.
                secondary_only = (singleton_tier1 and not focal_tier1)

                if focal_tier1:
                    focal_groups.append(group_name)

                if secondary_only:
                    secondary_groups.append(group_name)

            n_focal = len(focal_groups)

            n_secondary = len(secondary_groups)

            priority = assign_priority(
                n_focal,
                n_secondary,
                args.focal_recurrent_min_clades,
                args.secondary_recurrent_min_clades,
            )

            rows.append({
                "dmel_cre_id": cre,
                "scenario": scenario,
                "n_focal_tier1_clades": n_focal,
                "focal_tier1_clades": "|".join(focal_groups),
                "n_secondary_only_tier1_clades": n_secondary,
                "secondary_only_tier1_clades": "|".join(secondary_groups),
                "n_total_tier1_clades": n_focal + n_secondary,
                "candidate_priority": priority,
                "candidate_retained":("yes" if priority != "not_tier1_candidate" else "no"),
            })

    long = pd.DataFrame(rows)

    # ========================================================
    # Primary-scenario priority
    # ========================================================

    primary = (
        long.loc[
            long["scenario"] == args.primary_scenario,
            [
                "dmel_cre_id",
                "candidate_priority",
            ],
        ]
        .rename(
            columns={
                "candidate_priority":
                    "primary_candidate_priority",
            }
        )
    )

    long = long.merge(
        primary,
        on="dmel_cre_id",
        how="left",
        validate="many_to_one",
    )

    long["same_priority_as_primary"] = (
        long["candidate_priority"]
        == long["primary_candidate_priority"]
    ).map({
        True: "yes",
        False: "no",
    })


    # ========================================================
    # Candidate-level robustness summary
    # ========================================================

    summary_rows = []

    for cre, group in long.groupby(
        "dmel_cre_id",
        sort=True,
    ):

        retained = (
            group["candidate_retained"]
            == "yes"
        )

        same_priority = (
            group["same_priority_as_primary"]
            == "yes"
        )

        lost = group.loc[
            ~retained,
            "scenario",
        ].tolist()

        changed = group.loc[
            ~same_priority,
            [
                "scenario",
                "candidate_priority",
            ],
        ]

        summary_rows.append({
            "dmel_cre_id": cre,
            "primary_candidate_priority":
                group[
                    "primary_candidate_priority"
                ].iloc[0],
            "n_scenarios": len(group),
            "n_candidate_retained": int(retained.sum()),
            "percent_candidate_retained": 100 * retained.mean(),

            # Names retained for compatibility with existing
            # downstream plotting scripts.
            "robust_candidate_all_9": "yes" if retained.all() else "no",
            "n_same_priority_as_primary": int(same_priority.sum()),
            "percent_same_priority_as_primary": 100 * same_priority.mean(),
            "robust_priority_all_9": "yes" if same_priority.all() else "no",
            "lost_in_scenarios":( "|".join(lost) if lost else "NA"),
            "changed_priority_in_scenarios":
                (
                    "|".join(
                        f"{row.scenario}:{row.candidate_priority}"
                        for row in changed.itertuples()
                    )
                    if not changed.empty
                    else "NA"
                ),
        })

    summary = (
        pd.DataFrame(summary_rows).sort_values(
            [
                "percent_candidate_retained",
                "percent_same_priority_as_primary",
                "dmel_cre_id",
            ],
            ascending=[
                False,
                False,
                True,
            ],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )


    # ========================================================
    # Scenario matrices
    # ========================================================

    priority_matrix = (
        long
        .pivot(
            index="dmel_cre_id",
            columns="scenario",
            values="candidate_priority",
        )
        .reindex(
            columns=scenarios
        )
        .reset_index()
    )

    retention_matrix = (
        long
        .pivot(
            index="dmel_cre_id",
            columns="scenario",
            values="candidate_retained",
        )
        .reindex(
            columns=scenarios
        )
        .reset_index()
    )


    # ========================================================
    # Write outputs
    # ========================================================

    for path in [
        args.out_long,
        args.out_summary,
        args.out_priority_matrix,
        args.out_retention_matrix,
        args.metadata_out,
    ]:
        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    long.to_csv(
        args.out_long,
        sep="\t",
        index=False,
    )

    summary.to_csv(
        args.out_summary,
        sep="\t",
        index=False,
    )

    priority_matrix.to_csv(
        args.out_priority_matrix,
        sep="\t",
        index=False,
    )

    retention_matrix.to_csv(
        args.out_retention_matrix,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Run metadata
    # ========================================================

    metadata = pd.DataFrame([{
        "script": Path(__file__).name,
        "run_timestamp": datetime.now().astimezone().isoformat(),
        "primary_scenario": args.primary_scenario,
        "n_candidates": len(candidate_ids),
        "n_scenarios": len(scenarios),
        "n_focal_clades": len(groups),
        "focal_recurrent_min_clades": args.focal_recurrent_min_clades,
        "secondary_recurrent_min_clades": args.secondary_recurrent_min_clades,
        "n_candidates_retained_all_scenarios": int((summary["robust_candidate_all_9"] == "yes").sum()),
        "n_candidates_same_priority_all_scenarios": int((summary["robust_priority_all_9"] == "yes").sum()),
    }])

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
    print("FINAL TIER-1 CANDIDATE SENSITIVITY")
    print("=" * 72)

    print(
        f"Candidates evaluated: "
        f"{len(candidate_ids)}"
    )

    print(
        f"Sensitivity scenarios: "
        f"{len(scenarios)}"
    )

    print(
        f"Focal clades: "
        f"{len(groups)}"
    )

    print()
    print("Candidate retained in all scenarios:")

    print(summary["robust_candidate_all_9"].value_counts().to_string())

    print()
    print("Exact priority retained in all scenarios:")

    print(summary["robust_priority_all_9"].value_counts().to_string())

    print()
    print(summary[
            [
                "dmel_cre_id",
                "primary_candidate_priority",
                "percent_candidate_retained",
                "percent_same_priority_as_primary",
                "robust_candidate_all_9",
                "robust_priority_all_9",
            ]
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.1f}",
        )
    )

    print()
    print(
        f"Wrote candidate summary:\n"
        f"{args.out_summary}"
    )


if __name__ == "__main__":
    main()
