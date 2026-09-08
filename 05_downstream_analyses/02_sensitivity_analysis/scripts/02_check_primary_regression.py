#!/usr/bin/env python3

# ============================================================
# Check primary sensitivity-scenario regression
#
# Purpose:
#   Verify that the primary sensitivity scenario exactly
#   reproduces the baseline CRE-state classification.
#
# Input/output paths are supplied by
# config/sensitivity_config.sh via the pipeline wrapper.
# ============================================================

from pathlib import Path
from datetime import datetime
import argparse

import pandas as pd


VALID_STATES = {
    "positional_match",
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
            "Check whether the primary sensitivity scenario "
            "reproduces the baseline classification."
        )
    )

    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--sensitivity-dir", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)

    parser.add_argument("--scenario-manifest", type=Path, required=True)
    parser.add_argument("--primary-scenario", required=True)
    parser.add_argument("--out-summary", type=Path, required=True)
    parser.add_argument("--out-changes", type=Path, required=True)
    
    parser.add_argument("--metadata-out", type=Path, required=True)
    parser.add_argument("--expected-reference-cres", type=int, required=True)

    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_file(path):

    if not path.is_file():
        raise SystemExit(
            f"ERROR: required file not found:\n{path}"
        )


def load_species(path):

    require_file(path)

    species = [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip()
        and not line.lstrip().startswith("#")
    ]

    if not species:
        raise SystemExit(
            "ERROR: no target species found."
        )

    if len(species) != len(set(species)):
        raise SystemExit(
            "ERROR: duplicate target species."
        )

    return species


def load_states(
    path,
    expected_reference_cres,
):

    require_file(path)

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "dmel_cre_id",
        "class",
    }

    missing = required - set(df.columns)

    if missing:
        raise SystemExit(
            f"ERROR: {path} missing columns: "
            + ", ".join(sorted(missing))
        )

    if df["dmel_cre_id"].duplicated().any():
        raise SystemExit(
            f"ERROR: duplicate CRE IDs in:\n{path}"
        )

    if len(df) != expected_reference_cres:
        raise SystemExit(
            f"ERROR: unexpected number of CREs in:\n{path}\n"
            f"Expected: {expected_reference_cres}\n"
            f"Observed: {len(df)}"
        )

    unknown = sorted(
        set(df["class"])
        - VALID_STATES
    )

    if unknown:
        raise SystemExit(
            f"ERROR: unknown CRE state(s) in:\n{path}\n"
            + "\n".join(unknown)
        )

    return df[
        [
            "dmel_cre_id",
            "class",
        ]
    ].copy()


def load_primary_parameters(
    manifest_path,
    primary_scenario,
):

    require_file(manifest_path)

    manifest = pd.read_csv(
        manifest_path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "scenario",
        "reciprocal_overlap",
        "local_gene_distance_bp",
    }

    missing = required - set(manifest.columns)

    if missing:
        raise SystemExit(
            "ERROR: scenario manifest missing columns: "
            + ", ".join(sorted(missing))
        )

    match = manifest.loc[
        manifest["scenario"]
        == primary_scenario
    ]

    if len(match) != 1:
        raise SystemExit(
            "ERROR: primary scenario must occur exactly once "
            "in scenario manifest."
        )

    row = match.iloc[0]

    return (
        float(row["reciprocal_overlap"]),
        int(row["local_gene_distance_bp"]),
    )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    species = load_species(
        args.targets
    )

    overlap, distance = load_primary_parameters(
        args.scenario_manifest,
        args.primary_scenario,
    )

    scenario_dir = (
        args.sensitivity_dir
        / args.primary_scenario
    )

    if not scenario_dir.is_dir():
        raise SystemExit(
            f"ERROR: primary scenario directory not found:\n"
            f"{scenario_dir}"
        )

    for path in [
        args.out_summary,
        args.out_changes,
        args.metadata_out,
    ]:
        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    summary_rows = []
    changed_parts = []

    # --------------------------------------------------------
    # Compare baseline vs primary scenario
    # --------------------------------------------------------

    for sp in species:

        baseline_file = (
            args.baseline_dir
            / f"dmel_to_{sp}_cre_turnover.tsv"
        )

        scenario_file = (
            scenario_dir
            / f"dmel_to_{sp}_cre_turnover.tsv"
        )

        baseline = load_states(
            baseline_file,
            args.expected_reference_cres,
        )

        scenario = load_states(
            scenario_file,
            args.expected_reference_cres,
        )

        merged = baseline.merge(
            scenario,
            on="dmel_cre_id",
            how="outer",
            suffixes=(
                "_baseline",
                "_scenario",
            ),
            indicator=True,
            validate="one_to_one",
        )

        if not merged["_merge"].eq("both").all():
            raise SystemExit(
                f"ERROR: CRE-ID mismatch for species {sp}"
            )

        merged = merged.drop(
            columns="_merge"
        )

        changed = merged.loc[
            merged["class_baseline"]
            != merged["class_scenario"]
        ].copy()

        n_cre = len(merged)
        n_changed = len(changed)
        n_identical = n_cre - n_changed

        summary_rows.append({
            "species": sp,
            "n_cre": n_cre,
            "n_identical": n_identical,
            "n_changed": n_changed,
            "percent_identical":
                100 * n_identical / n_cre,
        })

        if not changed.empty:

            changed.insert(
                1,
                "species",
                sp,
            )

            changed["state_transition"] = (
                changed["class_baseline"]
                + "->"
                + changed["class_scenario"]
            )

            changed_parts.append(
                changed[
                    [
                        "dmel_cre_id",
                        "species",
                        "class_baseline",
                        "class_scenario",
                        "state_transition",
                    ]
                ]
            )

    # --------------------------------------------------------
    # Build outputs
    # --------------------------------------------------------

    summary = pd.DataFrame(
        summary_rows
    )

    if changed_parts:

        changes = pd.concat(
            changed_parts,
            ignore_index=True,
        )

    else:

        changes = pd.DataFrame(
            columns=[
                "dmel_cre_id",
                "species",
                "class_baseline",
                "class_scenario",
                "state_transition",
            ]
        )

    summary.to_csv(
        args.out_summary,
        sep="\t",
        index=False,
    )

    changes.to_csv(
        args.out_changes,
        sep="\t",
        index=False,
    )

    total = int(
        summary["n_cre"].sum()
    )

    n_changed = int(
        summary["n_changed"].sum()
    )

    expected_total = (
        len(species)
        * args.expected_reference_cres
    )

    if total != expected_total:
        raise SystemExit(
            "ERROR: unexpected total CRE x species count.\n"
            f"Expected: {expected_total}\n"
            f"Observed: {total}"
        )

    percent_identical = (
        100
        * (total - n_changed)
        / total
    )

    # --------------------------------------------------------
    # Metadata
    # --------------------------------------------------------

    metadata = pd.DataFrame([{
        "script":
            Path(__file__).name,

        "run_timestamp":
            datetime.now().astimezone().isoformat(),

        "primary_scenario":
            args.primary_scenario,

        "reciprocal_overlap":
            overlap,

        "local_gene_distance_bp":
            distance,

        "n_target_species":
            len(species),

        "n_reference_cres":
            args.expected_reference_cres,

        "n_cre_species":
            total,

        "n_changed":
            n_changed,

        "percent_identical":
            percent_identical,

        "status":
            "PASS" if n_changed == 0 else "FAIL",
    }])

    metadata.to_csv(
        args.metadata_out,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Console summary
    # --------------------------------------------------------

    print()
    print("=" * 72)
    print("PRIMARY REGRESSION CHECK")
    print("=" * 72)
    print()

    print(
        summary.to_string(
            index=False
        )
    )

    print()
    print(f"Primary scenario:    {args.primary_scenario}")
    print(f"Overlap:             {overlap}")
    print(f"Distance:            {distance} bp")
    print(f"Total CRE x species: {total}")
    print(f"Changed states:      {n_changed}")
    print(f"Percent identical:   {percent_identical:.4f}%")
    print()

    if n_changed:
        raise SystemExit(
            "ERROR: primary sensitivity scenario does not "
            "reproduce the baseline classification."
        )

    print(
        f"PASS: {args.primary_scenario} reproduces "
        "all baseline states."
    )


if __name__ == "__main__":
    main()
