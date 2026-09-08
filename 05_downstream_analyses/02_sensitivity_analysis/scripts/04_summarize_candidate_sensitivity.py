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

POSITIVE_STATES = {"positional_match", "turnover_candidate"}
VALID_STATES = POSITIVE_STATES | {"no_detected_CRE", "uncertain"}
VALID_PRIORITIES = {
    "focal_recurrent",
    "focal_plus_secondary_recurrent",
    "secondary_recurrent",
    "focal_single",
    "secondary_single",
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Re-evaluate prioritized Tier-1 CRE candidates across sensitivity scenarios."
    )
    p.add_argument("--sensitivity-dir", type=Path, required=True)
    p.add_argument("--scenario-manifest", type=Path, required=True)
    p.add_argument("--groups", type=Path, required=True)
    p.add_argument("--candidates", type=Path, required=True)
    p.add_argument("--candidate-metadata", type=Path, required=True)
    p.add_argument("--primary-scenario", required=True)
    p.add_argument("--out-long", type=Path, required=True)
    p.add_argument("--out-summary", type=Path, required=True)
    p.add_argument("--out-priority-matrix", type=Path, required=True)
    p.add_argument("--out-retention-matrix", type=Path, required=True)
    p.add_argument("--metadata-out", type=Path, required=True)
    p.add_argument("--focal-recurrent-min-clades", type=int, required=True)
    p.add_argument("--secondary-recurrent-min-clades", type=int, required=True)
    p.add_argument("--expected-reference-cres", type=int, required=True)
    return p.parse_args()


def require_file(path):
    if not path.is_file():
        raise SystemExit(f"ERROR: required file not found:\n{path}")


def require_columns(df, required, label):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(f"ERROR: {label} missing columns:\n" + "\n".join(missing))


def load_scenarios(path, primary):
    require_file(path)
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(df, {"scenario"}, "scenario manifest")
    scenarios = df["scenario"].tolist()
    if not scenarios or len(scenarios) != len(set(scenarios)):
        raise SystemExit("ERROR: scenario manifest is empty or contains duplicate scenarios.")
    if primary not in scenarios:
        raise SystemExit("ERROR: primary scenario not found in scenario manifest.")
    return scenarios


def load_groups(path):
    require_file(path)
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(df, {"group_name", "focal_species", "comparison_species"}, "focal-clade table")
    if df["group_name"].duplicated().any():
        raise SystemExit("ERROR: duplicate focal-clade names.")

    groups = {}
    for _, row in df.iterrows():
        focal = row["focal_species"].strip()
        comparisons = [x.strip() for x in row["comparison_species"].split("|") if x.strip()]
        if not focal or not comparisons:
            raise SystemExit(f"ERROR: incomplete focal-clade definition for {row['group_name']}.")
        species = [focal, *comparisons]
        if len(species) != len(set(species)):
            raise SystemExit(f"ERROR: duplicate species in focal clade {row['group_name']}.")
        groups[row["group_name"]] = {"focal": focal, "comparison": comparisons}
    return groups


def validate_candidate_metadata(path, focal_min, secondary_min):
    require_file(path)
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(
        df,
        {"focal_recurrent_min_clades", "secondary_recurrent_min_clades"},
        "candidate prioritization metadata",
    )
    if len(df) != 1:
        raise SystemExit("ERROR: candidate prioritization metadata must contain exactly one row.")

    observed_focal = int(df.iloc[0]["focal_recurrent_min_clades"])
    observed_secondary = int(df.iloc[0]["secondary_recurrent_min_clades"])
    if (observed_focal, observed_secondary) != (focal_min, secondary_min):
        raise SystemExit(
            "ERROR: sensitivity recurrence thresholds do not match the primary candidate analysis.\n"
            f"Primary: focal={observed_focal}, secondary={observed_secondary}\n"
            f"Sensitivity: focal={focal_min}, secondary={secondary_min}"
        )


def load_states(root, scenario, species, expected_n):
    path = root / scenario / f"dmel_to_{species}_cre_turnover.tsv"
    require_file(path)
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(df, {"dmel_cre_id", "class"}, str(path))
    if len(df) != expected_n or df["dmel_cre_id"].duplicated().any():
        raise SystemExit(f"ERROR: invalid CRE count or duplicate CRE IDs in:\n{path}")
    unknown = sorted(set(df["class"]) - VALID_STATES)
    if unknown:
        raise SystemExit(f"ERROR: unknown CRE states in:\n{path}\n" + "\n".join(unknown))
    return dict(zip(df["dmel_cre_id"], df["class"]))


def is_focal_tier1(focal_state, comparison_states):
    if len(set(comparison_states)) != 1:
        return False
    consensus = comparison_states[0]
    return focal_state in POSITIVE_STATES and consensus in POSITIVE_STATES and focal_state != consensus


def is_singleton_tier1(states):
    values = list(states.values())
    if any(state not in POSITIVE_STATES for state in values):
        return False
    counts = Counter(values)
    return len(counts) == 2 and sorted(counts.values()) == [1, len(values) - 1]


def assign_priority(n_focal, n_secondary, focal_min, secondary_min):
    if n_focal >= focal_min:
        return "focal_recurrent"
    if n_focal >= 1 and n_secondary >= 1:
        return "focal_plus_secondary_recurrent"
    if n_focal == 1 and n_secondary == 0:
        return "focal_single"
    if n_focal == 0 and n_secondary >= secondary_min:
        return "secondary_recurrent"
    if n_focal == 0 and n_secondary == 1:
        return "secondary_single"
    return "not_tier1_candidate"


def main():
    args = parse_args()
    if args.focal_recurrent_min_clades < 2 or args.secondary_recurrent_min_clades < 2:
        raise SystemExit("ERROR: recurrence thresholds must be >= 2.")
    if args.expected_reference_cres < 1:
        raise SystemExit("ERROR: expected-reference-cres must be >= 1.")

    scenarios = load_scenarios(args.scenario_manifest, args.primary_scenario)
    groups = load_groups(args.groups)
    validate_candidate_metadata(
        args.candidate_metadata,
        args.focal_recurrent_min_clades,
        args.secondary_recurrent_min_clades,
    )

    require_file(args.candidates)
    candidates = pd.read_csv(args.candidates, sep="\t", dtype=str, keep_default_na=False)
    require_columns(candidates, {"dmel_cre_id", "candidate_priority"}, "prioritized candidate table")
    if candidates.empty or candidates["dmel_cre_id"].duplicated().any():
        raise SystemExit("ERROR: candidate table is empty or contains duplicate CRE IDs.")
    unknown = sorted(set(candidates["candidate_priority"]) - VALID_PRIORITIES)
    if unknown:
        raise SystemExit("ERROR: unknown primary candidate priorities:\n" + "\n".join(unknown))

    candidate_ids = candidates["dmel_cre_id"].tolist()
    primary_priority = dict(zip(candidates["dmel_cre_id"], candidates["candidate_priority"]))
    species_needed = sorted({sp for group in groups.values() for sp in [group["focal"], *group["comparison"]]})

    rows = []
    for scenario in scenarios:
        state_tables = {
            sp: load_states(args.sensitivity_dir, scenario, sp, args.expected_reference_cres)
            for sp in species_needed
        }
        for cre in candidate_ids:
            focal_groups, secondary_groups = [], []
            for group_name, group in groups.items():
                focal = group["focal"]
                species = [focal, *group["comparison"]]
                states = {sp: state_tables[sp].get(cre, "MISSING") for sp in species}
                if "MISSING" in states.values():
                    raise SystemExit(f"ERROR: missing state for {cre} in {scenario}.")

                focal_tier1 = is_focal_tier1(states[focal], [states[sp] for sp in group["comparison"]])
                singleton_tier1 = is_singleton_tier1(states)
                if focal_tier1:
                    focal_groups.append(group_name)
                elif singleton_tier1:
                    secondary_groups.append(group_name)

            n_focal, n_secondary = len(focal_groups), len(secondary_groups)
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
                "focal_tier1_clades": "|".join(focal_groups) or "NA",
                "n_secondary_only_tier1_clades": n_secondary,
                "secondary_only_tier1_clades": "|".join(secondary_groups) or "NA",
                "n_total_tier1_clades": n_focal + n_secondary,
                "candidate_priority": priority,
                "candidate_retained": "yes" if priority != "not_tier1_candidate" else "no",
            })

    long = pd.DataFrame(rows)

    # The primary sensitivity scenario must reproduce the actual primary
    # candidate priorities, not merely a newly calculated internal baseline.
    primary_calc = long.loc[
        long["scenario"] == args.primary_scenario,
        ["dmel_cre_id", "candidate_priority"],
    ].set_index("dmel_cre_id")["candidate_priority"]
    mismatches = [
        cre for cre in candidate_ids
        if primary_calc.get(cre, "MISSING") != primary_priority[cre]
    ]
    if mismatches:
        detail = pd.DataFrame({
            "dmel_cre_id": mismatches,
            "primary_candidate_priority": [primary_priority[x] for x in mismatches],
            "recalculated_primary_priority": [primary_calc.get(x, "MISSING") for x in mismatches],
        })
        raise SystemExit(
            "ERROR: primary sensitivity scenario does not reproduce the primary candidate priorities.\n\n"
            + detail.to_string(index=False)
        )

    long["primary_candidate_priority"] = long["dmel_cre_id"].map(primary_priority)
    long["same_priority_as_primary"] = (
        long["candidate_priority"] == long["primary_candidate_priority"]
    ).map({True: "yes", False: "no"})

    summary_rows = []
    for cre, group in long.groupby("dmel_cre_id", sort=True):
        retained = group["candidate_retained"].eq("yes")
        same_priority = group["same_priority_as_primary"].eq("yes")
        lost = group.loc[~retained, "scenario"].tolist()
        changed = group.loc[~same_priority, ["scenario", "candidate_priority"]]
        summary_rows.append({
            "dmel_cre_id": cre,
            "primary_candidate_priority": primary_priority[cre],
            "n_scenarios": len(group),
            "n_candidate_retained": int(retained.sum()),
            "percent_candidate_retained": 100.0 * retained.mean(),
            "robust_candidate_all_scenarios": "yes" if retained.all() else "no",
            "n_same_priority_as_primary": int(same_priority.sum()),
            "percent_same_priority_as_primary": 100.0 * same_priority.mean(),
            "robust_priority_all_scenarios": "yes" if same_priority.all() else "no",
            "lost_in_scenarios": "|".join(lost) if lost else "NA",
            "changed_priority_in_scenarios": (
                "|".join(f"{row.scenario}:{row.candidate_priority}" for row in changed.itertuples())
                if not changed.empty else "NA"
            ),
        })

    summary = pd.DataFrame(summary_rows).sort_values(
        ["percent_candidate_retained", "percent_same_priority_as_primary", "dmel_cre_id"],
        ascending=[False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)

    priority_matrix = long.pivot(index="dmel_cre_id", columns="scenario", values="candidate_priority")
    retention_matrix = long.pivot(index="dmel_cre_id", columns="scenario", values="candidate_retained")
    priority_matrix = priority_matrix.reindex(index=candidate_ids, columns=scenarios).reset_index()
    retention_matrix = retention_matrix.reindex(index=candidate_ids, columns=scenarios).reset_index()

    for path in [
        args.out_long,
        args.out_summary,
        args.out_priority_matrix,
        args.out_retention_matrix,
        args.metadata_out,
    ]:
        path.parent.mkdir(parents=True, exist_ok=True)

    long.to_csv(args.out_long, sep="\t", index=False)
    summary.to_csv(args.out_summary, sep="\t", index=False)
    priority_matrix.to_csv(args.out_priority_matrix, sep="\t", index=False)
    retention_matrix.to_csv(args.out_retention_matrix, sep="\t", index=False)

    metadata = pd.DataFrame([{
        "script": Path(__file__).name,
        "run_timestamp": datetime.now().astimezone().isoformat(),
        "primary_scenario": args.primary_scenario,
        "n_candidates": len(candidate_ids),
        "n_scenarios": len(scenarios),
        "n_focal_clades": len(groups),
        "focal_recurrent_min_clades": args.focal_recurrent_min_clades,
        "secondary_recurrent_min_clades": args.secondary_recurrent_min_clades,
        "primary_priority_regression": "PASS",
        "n_candidates_retained_all_scenarios": int(summary["robust_candidate_all_scenarios"].eq("yes").sum()),
        "n_candidates_same_priority_all_scenarios": int(summary["robust_priority_all_scenarios"].eq("yes").sum()),
    }])
    metadata.to_csv(args.metadata_out, sep="\t", index=False)

    print("\n" + "=" * 72)
    print("TIER-1 CANDIDATE SENSITIVITY")
    print("=" * 72)
    print(f"Candidates: {len(candidate_ids)}")
    print(f"Scenarios:  {len(scenarios)}")
    print("Primary candidate-priority regression: PASS")
    print(f"Retained in all scenarios: {int(summary['robust_candidate_all_scenarios'].eq('yes').sum())}")
    print(f"Same priority in all scenarios: {int(summary['robust_priority_all_scenarios'].eq('yes').sum())}")
    print(f"Wrote: {args.out_summary}")


if __name__ == "__main__":
    main()

