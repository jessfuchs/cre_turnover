#!/usr/bin/env python3

from pathlib import Path
import argparse
import pandas as pd


PROJECT_DIR = Path.home() / "cre_turnover" / "project"

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


SCENARIOS = [
    "ov025_dist12000",
    "ov025_dist24000",
    "ov025_dist48000",

    "ov050_dist12000",
    "ov050_dist24000",
    "ov050_dist48000",

    "ov075_dist12000",
    "ov075_dist24000",
    "ov075_dist48000",
]


PRIMARY_SCENARIO = (
    "ov050_dist24000"
)


VALID_STATES = {
    "present",
    "turnover_candidate",
    "no_detected_CRE",
    "uncertain",
}


def parse_args():

    parser = argparse.ArgumentParser()

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

    return parser.parse_args()


def load_species(path):

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

    return species


def load_scenario(
    sensitivity_dir,
    scenario,
    species,
):

    parts = []

    for sp in species:

        path = (
            sensitivity_dir
            / scenario
            / f"dmel_to_{sp}_cre_turnover.tsv"
        )

        if not path.exists():

            raise SystemExit(
                f"ERROR: missing:\n{path}"
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

        missing = (
            required
            - set(df.columns)
        )

        if missing:

            raise SystemExit(
                f"ERROR: {path} missing columns: "
                + ", ".join(
                    sorted(missing)
                )
            )

        if df[
            "dmel_cre_id"
        ].duplicated().any():

            raise SystemExit(
                f"ERROR: duplicate CRE IDs in {path}"
            )

        unknown = sorted(
            set(df["class"].dropna())
            - VALID_STATES
        )

        if unknown:

            raise SystemExit(
                f"ERROR: unknown state(s) in {path}:\n"
                + "\n".join(unknown)
            )

        x = df[
            [
                "dmel_cre_id",
                "class",
            ]
        ].copy()

        x["species"] = sp

        parts.append(x)

    out = pd.concat(
        parts,
        ignore_index=True,
    )

    if out[
        [
            "dmel_cre_id",
            "species",
        ]
    ].duplicated().any():

        raise SystemExit(
            f"ERROR: duplicate CRE x species rows "
            f"in {scenario}"
        )

    return out


def main():

    args = parse_args()

    args.out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    species = load_species(
        args.targets
    )

    cache = {
        scenario:
            load_scenario(
                args.sensitivity_dir,
                scenario,
                species,
            )
        for scenario
        in SCENARIOS
    }

    primary = (
        cache[
            PRIMARY_SCENARIO
        ]
        .rename(
            columns={
                "class":
                    "primary_class"
            }
        )
    )

    summary_rows = []
    transition_rows = []
    species_rows = []
    cre_rows = []

    for scenario in SCENARIOS:

        test = cache[
            scenario
        ]

        x = primary.merge(
            test,
            on=[
                "dmel_cre_id",
                "species",
            ],
            how="outer",
            validate="one_to_one",
            indicator=True,
        )

        if not x[
            "_merge"
        ].eq(
            "both"
        ).all():

            raise SystemExit(
                f"ERROR: row mismatch for {scenario}"
            )

        x = x.drop(
            columns="_merge"
        )

        same = (
            x["primary_class"]
            == x["class"]
        )

        counts = (
            test["class"]
            .value_counts()
        )

        summary_rows.append({
            "scenario": scenario,
            "n_cells": len(x),
            "n_same_as_primary":
                int(same.sum()),
            "n_changed":
                int((~same).sum()),
            "percent_stable":
                100 * same.mean(),

            "present":
                int(
                    counts.get(
                        "present",
                        0,
                    )
                ),

            "turnover_candidate":
                int(
                    counts.get(
                        "turnover_candidate",
                        0,
                    )
                ),

            "no_detected_CRE":
                int(
                    counts.get(
                        "no_detected_CRE",
                        0,
                    )
                ),

            "uncertain":
                int(
                    counts.get(
                        "uncertain",
                        0,
                    )
                ),
        })


        # State transitions
        changed = x[
            ~same
        ].copy()

        if not changed.empty:

            trans = (
                changed
                .groupby(
                    [
                        "primary_class",
                        "class",
                    ]
                )
                .size()
                .reset_index(
                    name="n"
                )
            )

            trans.insert(
                0,
                "scenario",
                scenario,
            )

            transition_rows.append(
                trans
            )


        # Species stability
        for sp, g in x.groupby(
            "species"
        ):

            stable = (
                g["primary_class"]
                == g["class"]
            )

            species_rows.append({
                "scenario": scenario,
                "species": sp,
                "n_cells": len(g),
                "n_changed":
                    int(
                        (~stable).sum()
                    ),
                "percent_stable":
                    100
                    * stable.mean(),
            })


        # CRE stability
        for cre, g in x.groupby(
            "dmel_cre_id"
        ):

            stable = (
                g["primary_class"]
                == g["class"]
            )

            cre_rows.append({
                "scenario":
                    scenario,

                "dmel_cre_id":
                    cre,

                "n_species":
                    len(g),

                "n_changed":
                    int(
                        (~stable).sum()
                    ),

                "percent_stable":
                    100
                    * stable.mean(),
            })


    summary = pd.DataFrame(
        summary_rows
    )

    species_summary = pd.DataFrame(
        species_rows
    )

    cre_summary = pd.DataFrame(
        cre_rows
    )

    if transition_rows:

        transitions = pd.concat(
            transition_rows,
            ignore_index=True,
        )

    else:

        transitions = pd.DataFrame(
            columns=[
                "scenario",
                "primary_class",
                "class",
                "n",
            ]
        )


    summary.to_csv(
        args.out_dir
        / "sensitivity_global_summary.tsv",
        sep="\t",
        index=False,
    )


    species_summary.to_csv(
        args.out_dir
        / "sensitivity_species_stability.tsv",
        sep="\t",
        index=False,
    )


    transitions.to_csv(
        args.out_dir
        / "sensitivity_state_transitions.tsv",
        sep="\t",
        index=False,
    )


    cre_summary.to_csv(
        args.out_dir
        / "sensitivity_cre_stability.tsv",
        sep="\t",
        index=False,
    )


    print()
    print("=" * 72)
    print("GLOBAL SENSITIVITY ANALYSIS")
    print("=" * 72)
    print()

    print(
        summary.to_string(
            index=False,
            float_format=
                lambda x:
                    f"{x:.2f}",
        )
    )

    print()
    print(
        "Sensitivity analysis completed."
    )


if __name__ == "__main__":
    main()
