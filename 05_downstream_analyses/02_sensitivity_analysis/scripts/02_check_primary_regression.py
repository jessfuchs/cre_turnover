#!/usr/bin/env python3

from pathlib import Path
import argparse
import pandas as pd


PROJECT_DIR = Path.home() / "cre_turnover" / "project"

BASELINE_DIR = (
    PROJECT_DIR
    / "cre_classification"
    / "results"
    / "turnover_by_species"
)

SENSITIVITY_DIR = (
    PROJECT_DIR
    / "cre_classification"
    / "sensitivity"
)

TARGETS = (
    PROJECT_DIR
    / "pairwise_wga"
    / "target_species.txt"
)

OUT_DIR = (
    PROJECT_DIR
    / "downstream_analyses"
    / "sensitivity_analysis"
    / "results"
)

PRIMARY_SCENARIO = "ov050_dist24000"


def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Verify that the primary sensitivity scenario "
            "reproduces the baseline CRE-state classification."
        )
    )

    parser.add_argument(
        "--baseline-dir",
        type=Path,
        default=BASELINE_DIR,
    )

    parser.add_argument(
        "--sensitivity-dir",
        type=Path,
        default=SENSITIVITY_DIR,
    )

    parser.add_argument(
        "--targets",
        type=Path,
        default=TARGETS,
    )

    parser.add_argument(
        "--out-dir",
        type=Path,
        default=OUT_DIR,
    )

    parser.add_argument(
        "--primary-scenario",
        default=PRIMARY_SCENARIO,
    )

    return parser.parse_args()


def load_species(path):

    if not path.exists():
        raise SystemExit(
            f"ERROR: targets file not found:\n{path}"
        )

    species = [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip()
        and not line.lstrip().startswith("#")
    ]

    if not species:
        raise SystemExit(
            f"ERROR: no species found in:\n{path}"
        )

    return species


def load_states(path):

    if not path.exists():
        raise SystemExit(
            f"ERROR: missing file:\n{path}"
        )

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
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

    return df[
        [
            "dmel_cre_id",
            "class",
        ]
    ].copy()


def main():

    args = parse_args()

    args.out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    species = load_species(
        args.targets
    )

    scenario_dir = (
        args.sensitivity_dir
        / args.primary_scenario
    )

    summary_rows = []
    changed_parts = []

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
            baseline_file
        )

        scenario = load_states(
            scenario_file
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

        changed = merged[
            merged["class_baseline"]
            != merged["class_scenario"]
        ].copy()

        summary_rows.append({
            "species": sp,
            "n_cre": len(merged),
            "n_identical":
                len(merged) - len(changed),
            "n_changed":
                len(changed),
            "percent_identical":
                100
                * (len(merged) - len(changed))
                / len(merged),
        })

        if not changed.empty:

            changed.insert(
                1,
                "species",
                sp,
            )

            changed_parts.append(
                changed[
                    [
                        "dmel_cre_id",
                        "species",
                        "class_baseline",
                        "class_scenario",
                    ]
                ]
            )

    summary = pd.DataFrame(
        summary_rows
    )

    summary_path = (
        args.out_dir
        / "primary_regression_check.tsv"
    )

    changes_path = (
        args.out_dir
        / "primary_regression_state_changes.tsv"
    )

    summary.to_csv(
        summary_path,
        sep="\t",
        index=False,
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
            ]
        )

    changes.to_csv(
        changes_path,
        sep="\t",
        index=False,
    )

    total = int(
        summary["n_cre"].sum()
    )

    n_changed = int(
        summary["n_changed"].sum()
    )

    percent_identical = (
        100
        * (total - n_changed)
        / total
    )

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
    print(
        f"Primary scenario:     "
        f"{args.primary_scenario}"
    )

    print(
        f"Total CRE x species: "
        f"{total}"
    )

    print(
        f"Changed states:      "
        f"{n_changed}"
    )

    print(
        f"Percent identical:   "
        f"{percent_identical:.4f}%"
    )

    print()
    print(f"Wrote: {summary_path}")
    print(f"Wrote: {changes_path}")
    print()

    if n_changed:

        raise SystemExit(
            "ERROR: primary sensitivity setting does not "
            "reproduce the baseline classification."
        )

    print(
        f"PASS: {args.primary_scenario} "
        "reproduces all baseline states."
    )


if __name__ == "__main__":
    main()
