#!/usr/bin/env python3

# ============================================================
# Plot ranked Tier-1 candidate robustness
#
# Purpose:
#   Visualize the sensitivity robustness of prioritized Tier-1
#   CRE candidates.
#
# Encoding:
#   - x position: candidate retention across scenarios
#   - marker shape: primary candidate priority
#   - marker color: robustness class
#   - black outline: candidate priority stable across all scenarios
#
# Candidate categories and robustness classes are taken directly
# from the sensitivity-annotated candidate table and are not
# recalculated by this plotting script.
#
# Input/output paths and the moderate-robustness threshold are
# supplied by config/sensitivity_config.sh via the wrapper.
# ============================================================

from pathlib import Path
import argparse

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


# ============================================================
# Candidate-priority plotting setup
# ============================================================

PRIORITY_ORDER = {
    "focal_recurrent": 1,
    "focal_plus_secondary_recurrent": 2,
    "secondary_recurrent": 3,
    "focal_single": 4,
    "secondary_single": 5,
}

MARKERS = {
    "focal_recurrent": "*",
    "focal_plus_secondary_recurrent": "D",
    "secondary_recurrent": "s",
    "focal_single": "o",
    "secondary_single": "^",
}

PRIORITY_LABELS = {
    "focal_recurrent": "Focal recurrent",
    "focal_plus_secondary_recurrent": "Focal + secondary recurrent",
    "secondary_recurrent": "Secondary recurrent",
    "focal_single": "Focal single",
    "secondary_single": "Secondary single",
}


# ============================================================
# Robustness plotting setup
# ============================================================

ROBUSTNESS_COLORS = {
    "robust": "#2E7D32",
    "moderate": "#E6A700",
    "sensitive": "#C62828",
}

BACKGROUND_COLORS = {
    "robust": "#E6F4EA",
    "moderate": "#FFF5D6",
    "sensitive": "#FCE8E6",
}


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Plot ranked robustness of Tier-1 CRE candidates "
            "across sensitivity scenarios."
        )
    )
    
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-png", type=Path, required=True)
    parser.add_argument("--out-pdf", type=Path, required=True)
    parser.add_argument("--moderate-robustness-min", type=float, required=True)
    
    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_file(path):
    if not path.is_file():
        raise SystemExit(
            f"ERROR: required file not found:\n{path}"
        )


def require_columns(df, required):
    missing = sorted(
        set(required)
        - set(df.columns)
    )

    if missing:
        raise SystemExit(
            "ERROR: input table missing columns:\n"
            + "\n".join(missing)
        )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    require_file(args.input)
    
    if not 0 <= args.moderate_robustness_min < 100:
        raise SystemExit(
            "ERROR: moderate-robustness-min must be "
            "between 0 and <100."
        )
    args.out_png.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    args.out_pdf.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


    # ========================================================
    # Load candidate table
    # ========================================================

    df = pd.read_csv(
        args.input,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    require_columns(
        df,
        {
            "dmel_cre_id",
            "candidate_priority",
            "percent_candidate_retained",
            "robustness_class",
            "percent_same_priority_as_primary",
            "priority_stability",
        },
    )

    if df.empty:
        raise SystemExit(
            "ERROR: candidate table is empty."
        )

    if df["dmel_cre_id"].duplicated().any():
        raise SystemExit(
            "ERROR: duplicate CRE IDs in candidate table."
        )


    # ========================================================
    # Prepare plotting variables
    # ========================================================

    for column in [
        "percent_candidate_retained",
        "percent_same_priority_as_primary",
    ]:
        df[column] = pd.to_numeric(
            df[column],
            errors="raise",
        )

    unknown_priorities = sorted(
        set(df["candidate_priority"])
        - set(PRIORITY_ORDER)
    )

    if unknown_priorities:
        raise SystemExit(
            "ERROR: unknown candidate priorities:\n"
            + "\n".join(unknown_priorities)
        )

    unknown_robustness = sorted(
        set(df["robustness_class"])
        - set(ROBUSTNESS_COLORS)
    )

    if unknown_robustness:
        raise SystemExit(
            "ERROR: unknown robustness classes:\n"
            + "\n".join(unknown_robustness)
        )

    unknown_stability = sorted(
        set(df["priority_stability"])
        - {"stable", "variable"}
    )

    if unknown_stability:
        raise SystemExit(
            "ERROR: unknown priority-stability values:\n"
            + "\n".join(unknown_stability)
        )

    df["priority_rank"] = (
        df["candidate_priority"]
        .map(PRIORITY_ORDER)
    )


    # ========================================================
    # Candidate ordering
    # ========================================================

    # Candidates with the strongest robustness evidence are
    # displayed first. Candidate priority and exact priority
    # stability are used as secondary ordering criteria.
    df = (
        df
        .sort_values(
            [
                "percent_candidate_retained",
                "priority_rank",
                "percent_same_priority_as_primary",
                "dmel_cre_id",
            ],
            ascending=[
                False,
                True,
                False,
                True,
            ],
            kind="mergesort",
        )
        .reset_index(
            drop=True
        )
    )

    # Reverse y coordinates so the strongest candidate appears
    # at the top of the figure.
    df["y"] = range(
        len(df),
        0,
        -1,
    )


    # ========================================================
    # Figure layout
    # ========================================================

    height = max(
        7,
        len(df) * 0.26,
    )

    fig = plt.figure(
        figsize=(
            13.5,
            height,
        )
    )

    gs = fig.add_gridspec(
        nrows=1,
        ncols=2,
        width_ratios=[
            4.8,
            1.8,
        ],
        wspace=0.03,
    )

    ax = fig.add_subplot(
        gs[0, 0]
    )

    ax_leg = fig.add_subplot(
        gs[0, 1]
    )

    ax_leg.axis(
        "off"
    )


    # ========================================================
    # Robustness zones
    # ========================================================

    moderate_min = (
        args.moderate_robustness_min
    )

    # Sensitive candidates.
    ax.axvspan(
        0,
        moderate_min,
        color=BACKGROUND_COLORS["sensitive"],
        alpha=0.45,
        zorder=0,
    )

    # Moderately robust candidates.
    ax.axvspan(
        moderate_min,
        100,
        color=BACKGROUND_COLORS["moderate"],
        alpha=0.45,
        zorder=0,
    )

    # Fully robust candidates occur exactly at 100% retention.
    # A narrow region beyond 100 visually separates this class.
    ax.axvspan(
        100,
        103,
        color=BACKGROUND_COLORS["robust"],
        alpha=0.55,
        zorder=0,
    )

    ax.axvline(
        moderate_min,
        linestyle="--",
        linewidth=1,
        color="grey",
        alpha=0.8,
    )

    ax.axvline(
        100,
        linestyle="--",
        linewidth=1,
        color="grey",
        alpha=0.8,
    )


    # ========================================================
    # Candidate points
    # ========================================================

    # Marker shape represents the original candidate priority.
    # Color represents sensitivity robustness.
    # A black outline indicates an unchanged priority category
    # across all sensitivity scenarios.
    for row in df.itertuples():

        priority_stable = (
            row.priority_stability
            == "stable"
        )

        ax.scatter(
            row.percent_candidate_retained,
            row.y,
            s=95,
            marker=MARKERS[
                row.candidate_priority
            ],
            facecolor=ROBUSTNESS_COLORS[
                row.robustness_class
            ],
            edgecolor=(
                "black"
                if priority_stable
                else "none"
            ),
            linewidth=(
                1.5
                if priority_stable
                else 0
            ),
            zorder=3,
        )


    # ========================================================
    # Axes
    # ========================================================

    ax.set_yticks(
        df["y"]
    )

    ax.set_yticklabels(
        df["dmel_cre_id"],
        fontsize=8,
    )

    ax.set_xlim(
        0,
        103,
    )

    ax.set_xlabel(
        "Candidate retention across sensitivity scenarios (%)",
        fontsize=11,
    )

    ax.set_ylabel(
        "Tier-1 CRE candidate",
        fontsize=11,
    )

    ax.set_title(
        "Robustness of Tier-1 CRE candidates",
        fontsize=14,
        pad=10,
    )

    ax.grid(
        axis="x",
        linestyle=":",
        linewidth=0.8,
        alpha=0.35,
    )


    # ========================================================
    # Robustness-zone labels
    # ========================================================

    top = (
        df["y"].max()
        + 1.6
    )

    ax.text(
        moderate_min / 2,
        top,
        "sensitive",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
        color=ROBUSTNESS_COLORS["sensitive"],
    )

    ax.text(
        (
            moderate_min
            + 100
        ) / 2,
        top,
        "moderate",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
        color=ROBUSTNESS_COLORS["moderate"],
    )

    ax.text(
        101.2,
        top,
        "robust",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
        color=ROBUSTNESS_COLORS["robust"],
    )

    ax.set_ylim(
        0,
        top + 0.8,
    )


    # ========================================================
    # Legends
    # ========================================================

    priority_handles = [
        Line2D(
            [0],
            [0],
            marker=MARKERS[key],
            color="none",
            markerfacecolor="grey",
            markeredgecolor="grey",
            markersize=8,
            label=PRIORITY_LABELS[key],
        )
        for key in PRIORITY_ORDER
    ]

    legend1 = ax_leg.legend(
        handles=priority_handles,
        title="Primary candidate priority",
        loc="upper left",
        bbox_to_anchor=(
            0.0,
            1.0,
        ),
        frameon=False,
        fontsize=9,
        title_fontsize=10,
    )

    ax_leg.add_artist(
        legend1
    )


    robustness_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=color,
            markeredgecolor="none",
            markersize=8,
            label=label.capitalize(),
        )
        for label, color
        in ROBUSTNESS_COLORS.items()
    ]

    legend2 = ax_leg.legend(
        handles=robustness_handles,
        title="Candidate robustness",
        loc="upper left",
        bbox_to_anchor=(
            0.0,
            0.62,
        ),
        frameon=False,
        fontsize=9,
        title_fontsize=10,
    )

    ax_leg.add_artist(
        legend2
    )


    stability_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="lightgrey",
            markeredgecolor="black",
            markeredgewidth=1.5,
            markersize=8,
            label="Same priority in all scenarios",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="lightgrey",
            markeredgecolor="none",
            markersize=8,
            label="Priority changes",
        ),
    ]

    legend3 = ax_leg.legend(
        handles=stability_handles,
        title="Priority stability",
        loc="upper left",
        bbox_to_anchor=(
            0.0,
            0.34,
        ),
        frameon=False,
        fontsize=9,
        title_fontsize=10,
    )

    ax_leg.add_artist(
        legend3
    )

    for legend in [
        legend1,
        legend2,
        legend3,
    ]:
        legend._legend_box.align = "left"


    # ========================================================
    # Save
    # ========================================================

    fig.subplots_adjust(
        left=0.18,
        right=0.98,
        top=0.95,
        bottom=0.07,
    )

    fig.savefig(
        args.out_png,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        args.out_pdf,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )

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
