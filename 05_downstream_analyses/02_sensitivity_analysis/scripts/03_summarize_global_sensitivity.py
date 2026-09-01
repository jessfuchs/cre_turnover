#!/usr/bin/env python3

# ============================================================
# Summarize global CRE-classification sensitivity
#
# Purpose:
#   Quantify CRE-state stability across all sensitivity
#   scenarios relative to the primary scenario.
#
# Outputs:
#   - global scenario-level summary
#   - species-level stability
#   - CRE-level stability
#   - state-transition counts
#   - run metadata
#
# Input/output paths are supplied by
# config/sensitivity_config.sh via the pipeline wrapper.
# ============================================================

from pathlib import Path
from datetime import datetime
import argparse

import pandas as pd


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
            "Summarize CRE-state stability across sensitivity "
            "scenarios relative to the primary scenario."
        )
    )

    parser.add_argument("--sensitivity-dir", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--scenario-manifest", type=Path, required=True)
    parser.add_argument("--primary-scenario", required=True)
    
    parser.add_argument("--out-global", type=Path, required=True)
    parser.add_argument("--out-species", type=Path, required=True)
    parser.add_argument("--out-transitions", type=Path, required=True)
    parser.add_argument("--out-cre", type=Path, required=True)
    
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


def load_manifest(
    path,
    primary_scenario,
):

    require_file(path)

    df = pd.read_csv(
        path,
        sep="\t",
        keep_default_na=False,
    )

    required = {
        "scenario",
        "reciprocal_overlap",
        "local_gene_distance_bp",
    }

    missing = required - set(df.columns)

    if missing:
        raise SystemExit(
            "ERROR: scenario manifest missing columns: "
            + ", ".join(sorted(missing))
        )

    if df["scenario"].duplicated().any():
        raise SystemExit(
            "ERROR: duplicate scenarios in manifest."
        )

    if primary_scenario not in set(df["scenario"]):
        raise SystemExit(
            "ERROR: primary scenario not found in manifest."
        )

    return df


def load_scenario(
    sensitivity_dir,
    scenario,
    species,
    expected_reference_cres,
):

    parts = []

    for sp in species:

        path = (
            sensitivity_dir
            / scenario
            / f"dmel_to_{sp}_cre_turnover.tsv"
        )

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
                f"ERROR: unexpected CRE count in:\n{path}\n"
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

        part = df[
            [
                "dmel_cre_id",
                "class",
            ]
        ].copy()

        part["species"] = sp

        parts.append(part)

    return pd.concat(
        parts,
        ignore_index=True,
    )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    species = load_species(
        args.targets
    )

    manifest = load_manifest(
        args.scenario_manifest,
        args.primary_scenario,
    )

    scenarios = manifest[
        "scenario"
    ].tolist()

    parameters = {
        row["scenario"]: (
            float(row["reciprocal_overlap"]),
            int(row["local_gene_distance_bp"]),
        )
        for _, row in manifest.iterrows()
    }

    for path in [
        args.out_global,
        args.out_species,
        args.out_transitions,
        args.out_cre,
        args.metadata_out,
    ]:
        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    # --------------------------------------------------------
    # Load primary scenario once
    # --------------------------------------------------------

    primary = load_scenario(
        args.sensitivity_dir,
        args.primary_scenario,
        species,
        args.expected_reference_cres,
    ).rename(
        columns={
            "class":
                "primary_class",
        }
    )

    expected_cells = (
        len(species)
        * args.expected_reference_cres
    )

    if len(primary) != expected_cells:
        raise SystemExit(
            "ERROR: unexpected primary scenario size."
        )

    # --------------------------------------------------------
    # Compare all scenarios with primary
    # --------------------------------------------------------

    summary_rows = []
    species_rows = []
    cre_rows = []
    transition_rows = []

    for scenario in scenarios:

        overlap, distance = parameters[
            scenario
        ]

        test = load_scenario(
            args.sensitivity_dir,
            scenario,
            species,
            args.expected_reference_cres,
        )

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

        if not x["_merge"].eq("both").all():
            raise SystemExit(
                f"ERROR: CRE x species mismatch for {scenario}"
            )

        x = x.drop(
            columns="_merge"
        )

        same = (
            x["primary_class"]
            == x["class"]
        )

        n_same = int(same.sum())
        n_changed = int((~same).sum())

        counts = (
            test["class"]
            .value_counts()
        )

        # ----------------------------------------------------
        # Global summary
        # ----------------------------------------------------

        summary_rows.append({
            "scenario": scenario,
            "reciprocal_overlap": overlap,
            "local_gene_distance_bp": distance,
            "is_primary_scenario":
                "yes"
                if scenario == args.primary_scenario
                else "no",
            "n_cells": len(x),
            "n_same_as_primary": n_same,
            "n_changed": n_changed,
            "percent_stable":
                100 * n_same / len(x),
            "present":
                int(counts.get("present", 0)),
            "turnover_candidate":
                int(counts.get("turnover_candidate", 0)),
            "no_detected_CRE":
                int(counts.get("no_detected_CRE", 0)),
            "uncertain":
                int(counts.get("uncertain", 0)),
        })

        # ----------------------------------------------------
        # State transitions
        # ----------------------------------------------------

        changed = x.loc[
            ~same
        ]

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
                "local_gene_distance_bp",
                distance,
            )

            trans.insert(
                0,
                "reciprocal_overlap",
                overlap,
            )

            trans.insert(
                0,
                "scenario",
                scenario,
            )

            transition_rows.append(
                trans
            )

        # ----------------------------------------------------
        # Species-level stability
        # ----------------------------------------------------

        for sp in species:

            g = x.loc[
                x["species"] == sp
            ]

            stable = (
                g["primary_class"]
                == g["class"]
            )

            species_rows.append({
                "scenario": scenario,
                "reciprocal_overlap": overlap,
                "local_gene_distance_bp": distance,
                "species": sp,
                "n_cells": len(g),
                "n_changed":
                    int((~stable).sum()),
                "percent_stable":
                    100 * stable.mean(),
            })

        # ----------------------------------------------------
        # CRE-level stability
        # ----------------------------------------------------

        for cre, g in x.groupby(
            "dmel_cre_id",
            sort=True,
        ):

            stable = (
                g["primary_class"]
                == g["class"]
            )

            cre_rows.append({
                "scenario": scenario,
                "reciprocal_overlap": overlap,
                "local_gene_distance_bp": distance,
                "dmel_cre_id": cre,
                "n_species": len(g),
                "n_changed":
                    int((~stable).sum()),
                "percent_stable":
                    100 * stable.mean(),
            })

    # --------------------------------------------------------
    # Build output tables
    # --------------------------------------------------------

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
                "reciprocal_overlap",
                "local_gene_distance_bp",
                "primary_class",
                "class",
                "n",
            ]
        )

    # --------------------------------------------------------
    # Minimal final QC
    # --------------------------------------------------------

    if len(summary) != len(scenarios):
        raise SystemExit(
            "ERROR: incomplete scenario summary."
        )

    primary_row = summary.loc[
        summary["scenario"]
        == args.primary_scenario
    ]

    if (
        len(primary_row) != 1
        or int(primary_row.iloc[0]["n_changed"]) != 0
    ):
        raise SystemExit(
            "ERROR: primary scenario is not self-consistent."
        )

    # --------------------------------------------------------
    # Write outputs
    # --------------------------------------------------------

    summary.to_csv(
        args.out_global,
        sep="\t",
        index=False,
    )

    species_summary.to_csv(
        args.out_species,
        sep="\t",
        index=False,
    )

    transitions.to_csv(
        args.out_transitions,
        sep="\t",
        index=False,
    )

    cre_summary.to_csv(
        args.out_cre,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Metadata
    # --------------------------------------------------------

    non_primary = summary.loc[
        summary["scenario"]
        != args.primary_scenario
    ]

    metadata = pd.DataFrame([{
        "script":
            Path(__file__).name,

        "run_timestamp":
            datetime.now().astimezone().isoformat(),

        "primary_scenario":
            args.primary_scenario,

        "n_target_species":
            len(species),

        "n_reference_cres":
            args.expected_reference_cres,

        "n_scenarios":
            len(scenarios),

        "n_cells_per_scenario":
            expected_cells,

        "minimum_nonprimary_percent_stable":
            (
                float(
                    non_primary["percent_stable"].min()
                )
                if not non_primary.empty
                else 100.0
            ),
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
    print("GLOBAL SENSITIVITY ANALYSIS")
    print("=" * 72)
    print()

    print(
        summary.to_string(
            index=False,
            float_format=lambda x:
                f"{x:.2f}",
        )
    )

    print()
    print(
        "PASS: global sensitivity analysis completed."
    )
    print()


if __name__ == "__main__":
    main()
