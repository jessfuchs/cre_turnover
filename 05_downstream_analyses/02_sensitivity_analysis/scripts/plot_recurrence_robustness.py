#!/usr/bin/env python3

# ============================================================
# Plot recurrence and sensitivity robustness of Tier-1 candidates
#
# Purpose:
#   Relate phylogenetic recurrence of Tier-1 support to candidate
#   robustness across CRE-classification sensitivity scenarios.
#
# Encoding:
#   - x: number of clades showing Tier-1 support
#   - y: candidate retention across sensitivity scenarios
#   - orange: candidate with focal Tier-1 support
#   - grey: secondary-only Tier-1 support
#
# Recurrent and sensitivity-sensitive candidates are labelled.
# Bold labels indicate stable candidate priority across scenarios.
# ============================================================

from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


# ============================================================
# Plot settings
# ============================================================

FOCAL_COLOR = "#E69F00"
SECONDARY_COLOR = "#BDBDBD"
EDGE_COLOR = "#222222"
LEADER_COLOR = "#777777"

LABEL_OFFSETS = [
    (8, 8),
    (8, -10),
    (-8, 8),
    (-8, -10),
    (12, 4),
    (-12, 4),
]


# ============================================================
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot Tier-1 candidate recurrence versus sensitivity robustness."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-png", type=Path, required=True)
    parser.add_argument("--out-pdf", type=Path, required=True)
    parser.add_argument("--moderate-robustness-min", type=float, required=True)
    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_columns(df, required):
    missing = sorted(set(required) - set(df.columns))

    if missing:
        raise SystemExit(
            "ERROR: candidate table missing columns:\n"
            + "\n".join(missing)
        )


def short_id(cre_id):
    return str(cre_id).replace("DMEL_CRE_", "")


def as_bool(value):
    return str(value).strip().lower() in {"yes", "true", "1"}


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    if not args.input.is_file():
        raise SystemExit(f"ERROR: input file not found:\n{args.input}")

    if not 0 <= args.moderate_robustness_min < 100:
        raise SystemExit(
            "ERROR: moderate-robustness-min must be between 0 and <100."
        )

    args.out_png.parent.mkdir(parents=True, exist_ok=True)
    args.out_pdf.parent.mkdir(parents=True, exist_ok=True)

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
            "n_total_tier1_clades",
            "has_focal_support",
            "percent_candidate_retained",
            "robustness_class",
            "priority_stability",
        },
    )

    if df.empty:
        raise SystemExit("ERROR: candidate table is empty.")

    if df["dmel_cre_id"].duplicated().any():
        raise SystemExit("ERROR: duplicate CRE IDs in candidate table.")

    df["n_total_tier1_clades"] = pd.to_numeric(
        df["n_total_tier1_clades"], errors="raise"
    )
    df["percent_candidate_retained"] = pd.to_numeric(
        df["percent_candidate_retained"], errors="raise"
    )

    unknown_robustness = sorted(
        set(df["robustness_class"]) - {"robust", "moderate", "sensitive"}
    )
    if unknown_robustness:
        raise SystemExit(
            "ERROR: unknown robustness classes:\n"
            + "\n".join(unknown_robustness)
        )

    unknown_stability = sorted(
        set(df["priority_stability"]) - {"stable", "variable"}
    )
    if unknown_stability:
        raise SystemExit(
            "ERROR: unknown priority-stability values:\n"
            + "\n".join(unknown_stability)
        )

    df["is_focal"] = df["has_focal_support"].apply(as_bool)

    # ========================================================
    # Deterministic x offsets
    # ========================================================

    df = df.sort_values(
        ["n_total_tier1_clades", "percent_candidate_retained", "dmel_cre_id"]
    ).copy()

    df["x_plot"] = df["n_total_tier1_clades"].astype(float)

    # Slight horizontal offsets separate candidates sharing the
    # same recurrence count without introducing random jitter.
    for n_clades in sorted(df["n_total_tier1_clades"].unique()):
        idx = df.index[df["n_total_tier1_clades"] == n_clades]
        offsets = (
            np.array([0.0])
            if len(idx) == 1
            else np.linspace(-0.13, 0.13, len(idx))
        )
        df.loc[idx, "x_plot"] = float(n_clades) + offsets

    # ========================================================
    # Candidate labels
    # ========================================================

    # Labels focus on candidates that are recurrent across clades
    # or sensitive to the classification parameters.
    df["label_this"] = (
        (df["n_total_tier1_clades"] >= 2)
        | (df["robustness_class"] == "sensitive")
    )

    secondary = df.loc[~df["is_focal"]]
    focal = df.loc[df["is_focal"]]

    # ========================================================
    # Plot
    # ========================================================

    fig, ax = plt.subplots(figsize=(10.2, 6.4))

    ax.scatter(
        secondary["x_plot"],
        secondary["percent_candidate_retained"],
        s=62,
        facecolor=SECONDARY_COLOR,
        edgecolor=EDGE_COLOR,
        linewidth=0.65,
        alpha=0.82,
        marker="o",
        zorder=2,
    )

    ax.scatter(
        focal["x_plot"],
        focal["percent_candidate_retained"],
        s=88,
        facecolor=FOCAL_COLOR,
        edgecolor=EDGE_COLOR,
        linewidth=0.75,
        alpha=0.96,
        marker="o",
        zorder=3,
    )

    # Robustness thresholds correspond to the classes assigned
    # during candidate sensitivity analysis.
    ax.axhline(
        args.moderate_robustness_min,
        color="#AAAAAA",
        linestyle=":",
        linewidth=1.0,
        zorder=1,
    )

    ax.axhline(
        100,
        color="#666666",
        linestyle=":",
        linewidth=1.0,
        zorder=1,
    )

    # ========================================================
    # Labels
    # ========================================================

    label_df = df.loc[df["label_this"]].copy()

    # Focal candidates are labelled first, followed by secondary
    # candidates, to keep label placement deterministic.
    label_df["_label_priority"] = np.where(label_df["is_focal"], 0, 1)

    label_df = label_df.sort_values(
        [
            "_label_priority",
            "n_total_tier1_clades",
            "percent_candidate_retained",
            "dmel_cre_id",
        ],
        ascending=[True, False, False, True],
    )

    for i, row in enumerate(label_df.itertuples()):
        dx, dy = LABEL_OFFSETS[i % len(LABEL_OFFSETS)]

        ax.annotate(
            short_id(row.dmel_cre_id),
            xy=(row.x_plot, row.percent_candidate_retained),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=8.5,
            ha="left" if dx > 0 else "right",
            va="bottom" if dy > 0 else "top",
            fontweight="bold" if row.priority_stability == "stable" else "normal",
            arrowprops={
                "arrowstyle": "-",
                "linewidth": 0.55,
                "color": LEADER_COLOR,
                "shrinkA": 2,
                "shrinkB": 2,
            },
            zorder=5,
        )

    # ========================================================
    # Formatting
    # ========================================================

    max_clades = int(df["n_total_tier1_clades"].max())

    ax.set_xticks(range(1, max_clades + 1))
    ax.set_xlim(0.65, max_clades + 0.4)
    ax.set_ylim(0, 104)

    ax.set_xlabel(
        "Number of clades showing Tier-1 support",
        fontsize=11,
    )

    ax.set_ylabel(
        "Candidate retention across sensitivity scenarios (%)",
        fontsize=11,
    )

    ax.set_title(
        "Recurrence and robustness of Tier-1 CRE candidates",
        fontsize=13,
        pad=10,
    )

    ax.grid(
        axis="both",
        linestyle=":",
        linewidth=0.6,
        alpha=0.35,
    )

    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_handles = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="none",
            markerfacecolor=FOCAL_COLOR,
            markeredgecolor=EDGE_COLOR,
            markersize=8,
            label="Focal Tier-1 support",
        ),
        Line2D(
            [0], [0],
            marker="o",
            linestyle="none",
            markerfacecolor=SECONDARY_COLOR,
            markeredgecolor=EDGE_COLOR,
            markersize=8,
            label="Secondary-only Tier-1 support",
        ),
    ]

    ax.legend(
        handles=legend_handles,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0,
        fontsize=9,
    )

    fig.tight_layout()
    fig.subplots_adjust(right=0.78)

    # ========================================================
    # Save
    # ========================================================

    fig.savefig(args.out_png, dpi=300, bbox_inches="tight")
    fig.savefig(args.out_pdf, bbox_inches="tight")
    plt.close(fig)

    print()
    print("=" * 72)
    print("TIER-1 RECURRENCE AND ROBUSTNESS")
    print("=" * 72)
    print(f"Candidates plotted: {len(df)}")
    print(f"Candidates labelled: {int(df['label_this'].sum())}")
    print(f"Candidates with focal support: {int(df['is_focal'].sum())}")
    print()
    print(f"Wrote PNG:\n{args.out_png}")
    print(f"Wrote PDF:\n{args.out_pdf}")


if __name__ == "__main__":
    main()
