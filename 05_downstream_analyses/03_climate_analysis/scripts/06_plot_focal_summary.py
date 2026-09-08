#!/usr/bin/env python3

# ============================================================
# Plot focal Tier-1 summary
#
# Purpose:
#   Summarize clade-wide Tier-1 singleton contrasts by showing
#   the fraction assigned to the predefined focal lineage.
#
# Focal events are separated into:
#   - focal turnover_candidate
#   - focal positional_match
#
# Remaining singleton events are assigned to other species.
# The vertical marker shows the expected focal share under
# equal probability across all species in the clade.
# ============================================================

from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


# ============================================================
# Default paths
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
CLIMATE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = CLIMATE_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"

DEFAULT_INPUT = RESULTS_DIR / "focal_tier1_enrichment.tsv"
DEFAULT_OUT_PNG = FIG_DIR / "focal_tier1_summary.png"
DEFAULT_OUT_PDF = FIG_DIR / "focal_tier1_summary.pdf"


# ============================================================
# Plot settings
# ============================================================

# Colorblind-friendly palette used throughout the CRE figures.
TURNOVER_COLOR = "#E69F00"
PRESENT_COLOR = "#0072B2"
OTHER_COLOR = "#D9D9D9"
EDGE_COLOR = "#333333"


# ============================================================
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot focal Tier-1 singleton composition by clade."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-png", type=Path, default=DEFAULT_OUT_PNG)
    parser.add_argument("--out-pdf", type=Path, default=DEFAULT_OUT_PDF)
    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_columns(df, required):
    missing = sorted(set(required) - set(df.columns))

    if missing:
        raise SystemExit(
            "ERROR: focal Tier-1 enrichment table missing columns:\n"
            + "\n".join(missing)
        )


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    if not args.input.is_file():
        raise SystemExit(f"ERROR: input file not found:\n{args.input}")

    args.out_png.parent.mkdir(parents=True, exist_ok=True)
    args.out_pdf.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")

    require_columns(
        df,
        {
            "clade_display",
            "n_tier1_singletons",
            "n_focal_tier1",
            "n_focal_turnover",
            "n_focal_positional_match",
            "observed_focal_share",
            "expected_focal_share",
        },
    )

    if df.empty:
        raise SystemExit("ERROR: focal Tier-1 enrichment table is empty.")

    numeric_columns = [
        "n_tier1_singletons",
        "n_focal_tier1",
        "n_focal_turnover",
        "n_focal_positional_match",
        "observed_focal_share",
        "expected_focal_share",
    ]

    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="raise")

    # The focal count must equal the sum of both focal directions.
    if not (
        df["n_focal_tier1"]
        == df["n_focal_turnover"] + df["n_focal_positional_match"]
    ).all():
        raise SystemExit(
            "ERROR: n_focal_tier1 does not equal "
            "n_focal_turnover + n_focal_positional_match."
        )

    # Clades without singleton events have no defined composition.
    plot_df = df.loc[df["n_tier1_singletons"] > 0].copy()

    if plot_df.empty:
        raise SystemExit("ERROR: no Tier-1 singleton events available to plot.")

    plot_df["focal_turnover_share"] = (
        plot_df["n_focal_turnover"] / plot_df["n_tier1_singletons"]
    )

    plot_df["focal_positional_match_share"] = (
        plot_df["n_focal_positional_match"] / plot_df["n_tier1_singletons"]
    )

    plot_df["other_share"] = (
        1
        - plot_df["focal_turnover_share"]
        - plot_df["focal_positional_match_share"]
    )

    # ========================================================
    # Plot
    # ========================================================

    fig_height = max(4.8, len(plot_df) * 0.75)
    fig, ax = plt.subplots(figsize=(10.2, fig_height))

    y = range(len(plot_df))
    height = 0.68

    turnover_pct = plot_df["focal_turnover_share"] * 100
    positional_match_pct = plot_df["focal_positional_match_share"] * 100
    other_pct = plot_df["other_share"] * 100
    observed_pct = plot_df["observed_focal_share"] * 100
    expected_pct = plot_df["expected_focal_share"] * 100

    ax.barh(
        y,
        turnover_pct,
        height=height,
        color=TURNOVER_COLOR,
        edgecolor="white",
        linewidth=1.0,
    )

    ax.barh(
        y,
        positional_match_pct,
        left=turnover_pct,
        height=height,
        color=PRESENT_COLOR,
        edgecolor="white",
        linewidth=1.0,
    )

    ax.barh(
        y,
        other_pct,
        left=observed_pct,
        height=height,
        color=OTHER_COLOR,
        edgecolor="white",
        linewidth=1.0,
    )

    # Expected focal share under equal lineage probability.
    for i, value in enumerate(expected_pct):
        ax.plot(
            [value, value],
            [i - height / 2 - 0.08, i + height / 2 + 0.08],
            color=EDGE_COLOR,
            linewidth=2.0,
            zorder=5,
        )

    # Show focal / total singleton counts.
    for i, row in enumerate(plot_df.itertuples()):
        label = f"{int(row.n_focal_tier1)}/{int(row.n_tier1_singletons)}"

        ax.text(
            101.2,
            i,
            label,
            va="center",
            ha="left",
            fontsize=9,
        )

    # ========================================================
    # Formatting
    # ========================================================

    ax.set_yticks(list(y))
    ax.set_yticklabels(plot_df["clade_display"], fontsize=10)
    ax.invert_yaxis()

    ax.set_xlim(0, 108)

    ax.set_xlabel(
        "Share of clade-wide singleton Tier-1 contrasts (%)",
        fontsize=11,
        labelpad=12,
    )

    ax.grid(axis="x", linestyle=":", linewidth=0.6, alpha=0.35)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [
        Line2D(
            [0], [0],
            color=TURNOVER_COLOR,
            linewidth=10,
            label="Focal turnover",
        ),
        Line2D(
            [0], [0],
            color=PRESENT_COLOR,
            linewidth=10,
            label="Focal positional_match",
        ),
        Line2D(
            [0], [0],
            color=OTHER_COLOR,
            linewidth=10,
            label="Other singleton",
        ),
        Line2D(
            [0], [0],
            color=EDGE_COLOR,
            linewidth=2.5,
            label="Expected focal share",
        ),
    ]

    ax.legend(
        handles=handles,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0,
        fontsize=10,
        labelspacing=0.8,
        handletextpad=0.8,
        handlelength=1.0,
    )

    fig.tight_layout()
    fig.subplots_adjust(right=0.79)

    # ========================================================
    # Save
    # ========================================================

    fig.savefig(args.out_png, dpi=300, bbox_inches="tight")
    fig.savefig(args.out_pdf, bbox_inches="tight")
    plt.close(fig)

    print()
    print("=" * 72)
    print("FOCAL TIER-1 SUMMARY PLOT")
    print("=" * 72)
    print(f"Clades plotted: {len(plot_df)}")
    print()
    print(f"Wrote PNG:\n{args.out_png}")
    print(f"Wrote PDF:\n{args.out_pdf}")


if __name__ == "__main__":
    main()
