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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


# ============================================================
# Default paths
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
SENSITIVITY_DIR = SCRIPT_DIR.parent
DOWNSTREAM_DIR = SENSITIVITY_DIR.parent
CANDIDATE_DIR = DOWNSTREAM_DIR / "01_candidate_analysis"
FIG_DIR = SENSITIVITY_DIR / "results" / "figures"

DEFAULT_INPUT = (
    CANDIDATE_DIR / "results" / "tables" /
    "tier1_candidates_prioritized_with_sensitivity.tsv"
)
DEFAULT_OUT_PNG = FIG_DIR / "candidate_recurrence_robustness.png"
DEFAULT_OUT_PDF = FIG_DIR / "candidate_recurrence_robustness.pdf"
DEFAULT_MODERATE_ROBUSTNESS_MIN = 66.6

RECURRENT_PRIORITIES = {
    "focal_recurrent",
    "focal_plus_secondary_recurrent",
    "secondary_recurrent",
}
VALID_PRIORITIES = RECURRENT_PRIORITIES | {"focal_single", "secondary_single"}


# ============================================================
# Plot settings
# ============================================================

FOCAL_COLOR = "#E69F00"
SECONDARY_COLOR = "#BDBDBD"
STABLE_EDGE_COLOR = "#111111"
VARIABLE_EDGE_COLOR = "#777777"
LEADER_COLOR = "#777777"

FOCAL_MARKER = "o"
SECONDARY_MARKER = "s"
NONRECURRENT_SIZE = 48
RECURRENT_SIZE = 88
NONRECURRENT_ALPHA = 0.55
RECURRENT_ALPHA = 0.98

AXIS_LABEL_SIZE = 12
TICK_SIZE = 11
LEGEND_SIZE = 10
LEGEND_TITLE_SIZE = 11
ANNOT_SIZE = 9.5


# ============================================================
# Deterministic label placement
# ============================================================

LABEL_SPECS = {
    "DMEL_CRE_00065": {"dx": -18, "dy": 18, "ha": "right", "va": "bottom"},
    "DMEL_CRE_00078": {"dx": -18, "dy": -18, "ha": "right", "va": "top"},
    "DMEL_CRE_00080": {"dx": 0, "dy": 26, "ha": "center", "va": "bottom"},
    "DMEL_CRE_00275": {"dx": 0, "dy": -28, "ha": "center", "va": "top"},
    "DMEL_CRE_00310": {"dx": 18, "dy": -18, "ha": "left", "va": "top"},
    "DMEL_CRE_00152": {"dx": -18, "dy": 14, "ha": "right", "va": "bottom"},
    "DMEL_CRE_00309": {"dx": 18, "dy": 14, "ha": "left", "va": "bottom"},
    "DMEL_CRE_00052": {"dx": -18, "dy": 16, "ha": "right", "va": "bottom"},
    "DMEL_CRE_00330": {"dx": 18, "dy": 16, "ha": "left", "va": "bottom"},
    "DMEL_CRE_00151": {"dx": 18, "dy": 14, "ha": "left", "va": "bottom"},
    "DMEL_CRE_00331": {"dx": 18, "dy": 14, "ha": "left", "va": "bottom"},
}


# ============================================================
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot cross-clade recurrence and sensitivity robustness of Tier-1 CRE candidates."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-png", type=Path, default=DEFAULT_OUT_PNG)
    parser.add_argument("--out-pdf", type=Path, default=DEFAULT_OUT_PDF)
    parser.add_argument(
        "--moderate-robustness-min",
        type=float,
        default=DEFAULT_MODERATE_ROBUSTNESS_MIN,
        help="Minimum candidate-retention percentage defining moderate robustness.",
    )
    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_columns(df, required):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(
            "ERROR: candidate table missing columns:\n" + "\n".join(missing)
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

    df = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    require_columns(
        df,
        {
            "dmel_cre_id",
            "n_total_tier1_clades",
            "candidate_priority",
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
    df["_retention_group"] = df["percent_candidate_retained"].round(4)

    unknown_priorities = sorted(set(df["candidate_priority"]) - VALID_PRIORITIES)
    if unknown_priorities:
        raise SystemExit(
            "ERROR: unknown candidate priorities:\n" + "\n".join(unknown_priorities)
        )

    unknown_robustness = sorted(
        set(df["robustness_class"]) - {"robust", "moderate", "sensitive"}
    )
    if unknown_robustness:
        raise SystemExit(
            "ERROR: unknown robustness classes:\n" + "\n".join(unknown_robustness)
        )

    unknown_stability = sorted(
        set(df["priority_stability"]) - {"stable", "variable"}
    )
    if unknown_stability:
        raise SystemExit(
            "ERROR: unknown priority-stability values:\n" + "\n".join(unknown_stability)
        )

    df["is_focal"] = df["has_focal_support"].apply(as_bool)
    df["is_recurrent"] = df["candidate_priority"].isin(RECURRENT_PRIORITIES)
    df["is_stable"] = df["priority_stability"] == "stable"

    # Deterministic horizontal offsets for overlapping points.
    df = df.sort_values(
        ["n_total_tier1_clades", "percent_candidate_retained", "dmel_cre_id"],
        kind="mergesort",
    ).copy()
    df["x_plot"] = df["n_total_tier1_clades"].astype(float)

    for _, idx in df.groupby(
        ["n_total_tier1_clades", "_retention_group"], sort=True
    ).groups.items():
        idx = list(idx)
        if len(idx) == 1:
            offsets = np.array([0.0])
        else:
            n_clades = int(df.loc[idx[0], "n_total_tier1_clades"])
            width = 0.13 if n_clades == 1 else 0.10
            offsets = np.linspace(-width, width, len(idx))
        n_clades = float(df.loc[idx[0], "n_total_tier1_clades"])
        df.loc[idx, "x_plot"] = n_clades + offsets

    df["label_this"] = df["is_recurrent"]

    # ========================================================
    # Plot
    # ========================================================

    fig, ax = plt.subplots(figsize=(12.6, 7.2))

    def plot_candidates(subset, marker, facecolor, recurrent):
        if subset.empty:
            return
        edgecolors = np.where(
            subset["is_stable"], STABLE_EDGE_COLOR, VARIABLE_EDGE_COLOR
        )
        linewidths = np.where(subset["is_stable"], 1.8, 0.7)
        ax.scatter(
            subset["x_plot"],
            subset["percent_candidate_retained"],
            s=RECURRENT_SIZE if recurrent else NONRECURRENT_SIZE,
            facecolor=facecolor,
            edgecolor=edgecolors,
            linewidths=linewidths,
            alpha=RECURRENT_ALPHA if recurrent else NONRECURRENT_ALPHA,
            marker=marker,
            zorder=4 if recurrent else 2,
        )

    plot_candidates(
        df.loc[df["is_focal"] & ~df["is_recurrent"]],
        FOCAL_MARKER,
        FOCAL_COLOR,
        False,
    )
    plot_candidates(
        df.loc[~df["is_focal"] & ~df["is_recurrent"]],
        SECONDARY_MARKER,
        SECONDARY_COLOR,
        False,
    )
    plot_candidates(
        df.loc[df["is_focal"] & df["is_recurrent"]],
        FOCAL_MARKER,
        FOCAL_COLOR,
        True,
    )
    plot_candidates(
        df.loc[~df["is_focal"] & df["is_recurrent"]],
        SECONDARY_MARKER,
        SECONDARY_COLOR,
        True,
    )

    ax.axhline(
        args.moderate_robustness_min,
        color="#AAAAAA",
        linestyle=":",
        linewidth=0.9,
        zorder=1,
    )
    ax.axhline(100, color="#777777", linestyle=":", linewidth=0.9, zorder=1)

    # Candidate labels: recurrent candidates only.
    label_df = df.loc[df["label_this"]].sort_values(
        ["n_total_tier1_clades", "percent_candidate_retained", "dmel_cre_id"],
        ascending=[False, False, True],
    )

    default_spec = {"dx": 10, "dy": 10, "ha": "left", "va": "bottom"}
    for row in label_df.itertuples():
        spec = LABEL_SPECS.get(row.dmel_cre_id, default_spec)
        ax.annotate(
            short_id(row.dmel_cre_id),
            xy=(row.x_plot, row.percent_candidate_retained),
            xytext=(spec["dx"], spec["dy"]),
            textcoords="offset points",
            fontsize=ANNOT_SIZE,
            ha=spec["ha"],
            va=spec["va"],
            fontweight="normal",
            arrowprops={
                "arrowstyle": "-",
                "linewidth": 0.55,
                "color": LEADER_COLOR,
                "shrinkA": 2,
                "shrinkB": 2,
            },
            zorder=6,
        )

    # ========================================================
    # Axes and formatting
    # ========================================================

    max_clades = int(df["n_total_tier1_clades"].max())
    ax.set_xticks(range(1, max_clades + 1))
    ax.set_xlim(0.65, max_clades + 0.4)

    min_retention = float(df["percent_candidate_retained"].min())
    ymin = max(0, np.floor(min_retention / 10) * 10 - 5)
    ax.set_ylim(ymin, 105)

    ax.set_xlabel(
        "Number of analysed clades with Tier-1 support", fontsize=AXIS_LABEL_SIZE
    )
    ax.set_ylabel(
        "Candidate retention across sensitivity scenarios (%)",
        fontsize=AXIS_LABEL_SIZE,
    )
    ax.tick_params(axis="both", labelsize=TICK_SIZE)
    ax.grid(axis="both", linestyle=":", linewidth=0.6, alpha=0.35)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ========================================================
    # Legends
    # ========================================================

    support_handles = [
        Line2D(
            [0], [0], marker=FOCAL_MARKER, linestyle="none",
            markerfacecolor=FOCAL_COLOR, markeredgecolor=VARIABLE_EDGE_COLOR,
            markeredgewidth=0.7, markersize=8, label="Focal Tier-1 support",
        ),
        Line2D(
            [0], [0], marker=SECONDARY_MARKER, linestyle="none",
            markerfacecolor=SECONDARY_COLOR, markeredgecolor=VARIABLE_EDGE_COLOR,
            markeredgewidth=0.7, markersize=8, label="Secondary-only Tier-1 support",
        ),
    ]

    support_legend = ax.legend(
        handles=support_handles,
        title="Candidate support",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.03, 1.00),
        borderaxespad=0,
        fontsize=LEGEND_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
        handletextpad=0.7,
        labelspacing=0.4,
    )
    support_legend._legend_box.align = "left"
    support_legend.get_title().set_ha("left")
    ax.add_artist(support_legend)

    stability_handles = [
        Line2D(
            [0], [0], marker="o", linestyle="none", markerfacecolor="white",
            markeredgecolor=STABLE_EDGE_COLOR, markeredgewidth=1.8,
            markersize=8, label="Stable priority",
        ),
        Line2D(
            [0], [0], marker="o", linestyle="none", markerfacecolor="white",
            markeredgecolor=VARIABLE_EDGE_COLOR, markeredgewidth=0.7,
            markersize=8, label="Variable priority",
        ),
    ]

    stability_legend = ax.legend(
        handles=stability_handles,
        title="Priority stability",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.03, 0.79),
        borderaxespad=0,
        fontsize=LEGEND_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
        handletextpad=0.7,
        labelspacing=0.4,
    )
    stability_legend._legend_box.align = "left"
    stability_legend.get_title().set_ha("left")

    # No internal figure title; use the available vertical space.
    fig.subplots_adjust(left=0.10, right=0.72, bottom=0.12, top=0.94)

    # ========================================================
    # Save
    # ========================================================

    extra_artists = (support_legend, stability_legend)
    fig.savefig(
        args.out_png,
        dpi=300,
        bbox_inches="tight",
        bbox_extra_artists=extra_artists,
        pad_inches=0.25,
    )
    fig.savefig(
        args.out_pdf,
        bbox_inches="tight",
        bbox_extra_artists=extra_artists,
        pad_inches=0.25,
    )
    plt.close(fig)

    print()
    print("=" * 72)
    print("TIER-1 RECURRENCE AND ROBUSTNESS")
    print("=" * 72)
    print(f"Candidates plotted: {len(df)}")
    print(f"Candidates labelled: {int(df['label_this'].sum())}")
    print(f"Candidates with focal support: {int(df['is_focal'].sum())}")
    print(f"Recurrent candidates: {int(df['is_recurrent'].sum())}")
    print(f"Stable-priority candidates: {int(df['is_stable'].sum())}")
    print()
    print(f"Wrote PNG:\n{args.out_png}")
    print(f"Wrote PDF:\n{args.out_pdf}")


if __name__ == "__main__":
    main()
