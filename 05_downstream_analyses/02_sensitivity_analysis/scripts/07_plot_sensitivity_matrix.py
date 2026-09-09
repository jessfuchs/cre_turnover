#!/usr/bin/env python3

# ============================================================
# Plot Tier-1 candidate sensitivity matrix
#
# Purpose:
#   Visualize whether prioritized Tier-1 CRE candidates are
#   retained across CRE-classification sensitivity scenarios.
#
# Encoding:
#   - rows: prioritized Tier-1 CRE candidates
#   - columns: sensitivity scenarios
#   - green: candidate retained
#   - red: candidate not retained
#   - black outline: primary sensitivity scenario
#
# Scenario order, overlap thresholds, distance thresholds and
# the primary scenario are supplied by the sensitivity pipeline.
# ============================================================

from pathlib import Path
import argparse

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import pandas as pd


# ============================================================
# Candidate priority
# ============================================================

PRIORITY_ORDER = {
    "focal_recurrent": 1,
    "focal_plus_secondary_recurrent": 2,
    "secondary_recurrent": 3,
    "focal_single": 4,
    "secondary_single": 5,
}


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Plot Tier-1 candidate retention across "
            "sensitivity scenarios."
        )
    )

    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--scenario-manifest", type=Path, required=True)
    parser.add_argument("--primary-scenario", required=True)
    parser.add_argument("--out-png", type=Path, required=True)
    parser.add_argument("--out-pdf", type=Path, required=True)

    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_file(path):
    if not path.is_file():
        raise SystemExit(
            f"ERROR: required file not found:\n{path}"
        )


def require_columns(df, required, label):
    missing = sorted(set(required) - set(df.columns))

    if missing:
        raise SystemExit(
            f"ERROR: {label} missing columns:\n"
            + "\n".join(missing)
        )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()
    for path in [
        args.matrix,
        args.summary,
        args.scenario_manifest,
    ]:
        require_file(path)
    args.out_png.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.out_pdf.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


    # ========================================================
    # Load scenario definitions
    # ========================================================

    scenarios = pd.read_csv(
        args.scenario_manifest,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    require_columns(
        scenarios,
        {
            "scenario",
            "reciprocal_overlap",
            "local_gene_distance_bp",
        },
        "scenario manifest",
    )

    if scenarios["scenario"].duplicated().any():
        raise SystemExit(
            "ERROR: duplicate scenarios in scenario manifest."
        )

    scenario_order = scenarios[
        "scenario"
    ].tolist()

    if args.primary_scenario not in scenario_order:
        raise SystemExit(
            "ERROR: primary scenario not found "
            "in scenario manifest."
        )

    scenarios["reciprocal_overlap"] = pd.to_numeric(
        scenarios["reciprocal_overlap"],
        errors="raise",
    )

    scenarios["local_gene_distance_bp"] = pd.to_numeric(
        scenarios["local_gene_distance_bp"],
        errors="raise",
    )

    # Labels are derived directly from the scenario parameters
    # rather than maintained separately in the plotting script.
    scenario_labels = [
        (
            f"{row.reciprocal_overlap:.2f} / "
            f"{row.local_gene_distance_bp / 1000:g} kb"
        )
        for row in scenarios.itertuples()
    ]


    # ========================================================
    # Load candidate retention and ranking information
    # ========================================================

    matrix = pd.read_csv(
        args.matrix,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    summary = pd.read_csv(
        args.summary,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    require_columns(
        matrix,
        {
            "dmel_cre_id",
            *scenario_order,
        },
        "candidate retention matrix",
    )

    require_columns(
        summary,
        {
            "dmel_cre_id",
            "percent_candidate_retained",
            "primary_candidate_priority",
        },
        "candidate sensitivity summary",
    )

    if matrix["dmel_cre_id"].duplicated().any():
        raise SystemExit(
            "ERROR: duplicate CRE IDs in retention matrix."
        )

    if summary["dmel_cre_id"].duplicated().any():
        raise SystemExit(
            "ERROR: duplicate CRE IDs in sensitivity summary."
        )


    # ========================================================
    # Merge ranking information
    # ========================================================

    df = matrix.merge(
        summary[
            [
                "dmel_cre_id",
                "percent_candidate_retained",
                "primary_candidate_priority",
            ]
        ],
        on="dmel_cre_id",
        how="left",
        validate="one_to_one",
    )

    if df["percent_candidate_retained"].eq("").any():
        raise SystemExit(
            "ERROR: retention summary missing for one or "
            "more candidate CREs."
        )

    df["percent_candidate_retained"] = pd.to_numeric(
        df["percent_candidate_retained"],
        errors="raise",
    )

    df["priority_rank"] = (
        df["primary_candidate_priority"]
        .map(PRIORITY_ORDER)
    )

    if df["priority_rank"].isna().any():

        unknown = sorted(
            df.loc[
                df["priority_rank"].isna(),
                "primary_candidate_priority",
            ].unique()
        )

        raise SystemExit(
            "ERROR: unknown candidate priorities:\n"
            + "\n".join(unknown)
        )


    # ========================================================
    # Candidate ordering
    # ========================================================

    # Most robust candidates are shown first, with primary
    # candidate priority used as a secondary ranking criterion.
    df = (
        df
        .sort_values(
            [
                "percent_candidate_retained",
                "priority_rank",
                "dmel_cre_id",
            ],
            ascending=[
                False,
                True,
                True,
            ],
            kind="mergesort",
        )
        .reset_index(
            drop=True
        )
    )


    # ========================================================
    # Convert retention states to numeric matrix
    # ========================================================

    retention_values = set(
        df[scenario_order]
        .stack()
        .unique()
    )

    unexpected = sorted(
        retention_values - {"yes", "no"}
    )

    if unexpected:
        raise SystemExit(
            "ERROR: unexpected candidate-retention values:\n"
            + "\n".join(unexpected)
        )

    values = (
        df[scenario_order]
        .replace({
            "yes": 1,
            "no": 0,
        })
        .astype(int)
        .to_numpy()
    )


    # ========================================================
    # Plot
    # ========================================================

    fig_height = max(
        7,
        len(df) * 0.23,
    )

    fig, ax = plt.subplots(
        figsize=(
            9,
            fig_height,
        )
    )

    cmap = ListedColormap([
        "#C62828",
        "#2E7D32",
    ])

    ax.imshow(
        values,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=1,
    )


    # ========================================================
    # Axis labels
    # ========================================================

    ax.set_xticks(
        range(len(scenario_order))
    )

    ax.set_xticklabels(
        scenario_labels,
        rotation=45,
        ha="right",
    )

    ax.set_yticks(
        range(len(df))
    )

    ax.set_yticklabels(
        df["dmel_cre_id"],
        fontsize=8,
    )

    ax.set_xlabel(
        "Reciprocal overlap threshold / local distance threshold"
    )

    ax.set_ylabel(
        "Tier-1 CRE candidate"
    )

    ax.set_title(
        "Retention of Tier-1 candidates across sensitivity scenarios"
    )


    # ========================================================
    # Highlight primary scenario
    # ========================================================

    primary_index = scenario_order.index(
        args.primary_scenario
    )

    ax.axvline(
        primary_index - 0.5,
        color="black",
        linewidth=1.5,
    )

    ax.axvline(
        primary_index + 0.5,
        color="black",
        linewidth=1.5,
    )


    # ========================================================
    # Separate reciprocal-overlap groups
    # ========================================================

    overlaps = scenarios[
        "reciprocal_overlap"
    ].tolist()

    for index in range(
        1,
        len(overlaps),
    ):

        if overlaps[index] != overlaps[index - 1]:

            ax.axvline(
                index - 0.5,
                color="white",
                linewidth=2,
            )


    # ========================================================
    # Legend
    # ========================================================

    legend_handles = [
        Patch(
            facecolor="#2E7D32",
            label="Tier-1 candidate retained",
        ),
        Patch(
            facecolor="#C62828",
            label="Candidate not retained",
        ),
    ]

    ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(
            1.02,
            1.0,
        ),
        frameon=False,
    )

    fig.tight_layout()


    # ========================================================
    # Save
    # ========================================================

    fig.savefig(
        args.out_png,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        args.out_pdf,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(
        f"Wrote PNG:\n"
        f"{args.out_png}"
    )

    print(
        f"Wrote PDF:\n"
        f"{args.out_pdf}"
    )


if __name__ == "__main__":
    main()
