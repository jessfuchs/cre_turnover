#!/usr/bin/env python3

# ============================================================
# Plot focal Tier-1 enrichment
#
# Purpose:
#   Compare the observed fraction of Tier-1 singleton contrasts
#   assigned to the predefined focal lineage with the null
#   expectation of equal probability across clade members.
#
# The observed focal share is shown with an exact 95% binomial
# confidence interval calculated in the enrichment analysis.
#
# Input/output paths are supplied by the climate-analysis
# configuration via the pipeline wrapper.
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
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot observed versus expected focal Tier-1 enrichment."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-png", type=Path, required=True)
    parser.add_argument("--out-pdf", type=Path, required=True)
    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_file(path):
    if not path.is_file():
        raise SystemExit(f"ERROR: required file not found:\n{path}")


def require_columns(df, required):
    missing = sorted(set(required) - set(df.columns))

    if missing:
        raise SystemExit(
            "ERROR: enrichment table missing columns:\n"
            + "\n".join(missing)
        )


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    require_file(args.input)
    args.out_png.parent.mkdir(parents=True, exist_ok=True)
    args.out_pdf.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")

    require_columns(
        df,
        {
            "clade_display",
            "observed_focal_share",
            "expected_focal_share",
            "ci95_low",
            "ci95_high",
            "n_focal_tier1",
            "n_tier1_singletons",
        },
    )

    if df.empty:
        raise SystemExit("ERROR: focal Tier-1 enrichment table is empty.")

    numeric_columns = [
        "observed_focal_share",
        "expected_focal_share",
        "ci95_low",
        "ci95_high",
        "n_focal_tier1",
        "n_tier1_singletons",
    ]

    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    # Clades without any Tier-1 singleton event have no defined
    # observed focal proportion and cannot be displayed here.
    plot_df = df.loc[
        df["observed_focal_share"].notna()
        & df["ci95_low"].notna()
        & df["ci95_high"].notna()
    ].copy()

    if plot_df.empty:
        raise SystemExit(
            "ERROR: no clades with defined focal Tier-1 enrichment."
        )

    observed = plot_df["observed_focal_share"].to_numpy() * 100
    expected = plot_df["expected_focal_share"].to_numpy() * 100
    ci_low = plot_df["ci95_low"].to_numpy() * 100
    ci_high = plot_df["ci95_high"].to_numpy() * 100
    y = np.arange(len(plot_df))

    # ========================================================
    # Plot
    # ========================================================

    fig, ax = plt.subplots(figsize=(8.2, 5.6))

    # Connecting lines emphasize deviation from the clade-specific
    # null expectation.
    for i in range(len(plot_df)):
        ax.plot(
            [expected[i], observed[i]],
            [y[i], y[i]],
            color="#B0B0B0",
            linewidth=1.2,
            zorder=1,
        )

    xerr = np.vstack([
        observed - ci_low,
        ci_high - observed,
    ])

    ax.errorbar(
        observed,
        y,
        xerr=xerr,
        fmt="o",
        markersize=7,
        color="#222222",
        markerfacecolor="#E69F00",
        markeredgecolor="#222222",
        linewidth=1.1,
        capsize=3,
        zorder=3,
    )

    ax.scatter(
        expected,
        y,
        s=55,
        facecolor="white",
        edgecolor="#222222",
        linewidth=1.2,
        zorder=4,
    )

    # Show the underlying focal / total singleton counts next to
    # each observed proportion.
    for i, row in enumerate(plot_df.itertuples()):
        label = f"{int(row.n_focal_tier1)}/{int(row.n_tier1_singletons)}"

        if observed[i] > 92:
            x_text = observed[i] - 2
            ha = "right"
        else:
            x_text = observed[i] + 2
            ha = "left"

        ax.text(
            x_text,
            y[i],
            label,
            va="center",
            ha=ha,
            fontsize=8,
        )

    # ========================================================
    # Formatting
    # ========================================================

    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["clade_display"], fontsize=10)
    ax.invert_yaxis()

    ax.set_xlim(0, 100)

    ax.set_xlabel(
        "Tier-1 singleton contrasts assigned to focal lineage (%)",
        fontsize=11,
    )

    ax.set_title(
        "Focal Tier-1 contrasts relative to clade-wide singleton background",
        fontsize=13,
        pad=10,
    )

    ax.grid(
        axis="x",
        linestyle=":",
        linewidth=0.6,
        alpha=0.5,
    )

    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_handles = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#E69F00",
            markeredgecolor="#222222",
            markersize=7,
            label="Observed focal share",
        ),
        Line2D(
            [0], [0],
            marker="o",
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor="#222222",
            markersize=7,
            label="Expected under equal lineage probability",
        ),
    ]

    ax.legend(
        handles=legend_handles,
        frameon=False,
        loc="lower right",
    )

    fig.tight_layout()

    # ========================================================
    # Save
    # ========================================================

    fig.savefig(args.out_png, dpi=300, bbox_inches="tight")
    fig.savefig(args.out_pdf, bbox_inches="tight")
    plt.close(fig)

    print()
    print("=" * 72)
    print("FOCAL TIER-1 ENRICHMENT PLOT")
    print("=" * 72)
    print(f"Clades plotted: {len(plot_df)}")
    print()
    print(f"Wrote PNG:\n{args.out_png}")
    print(f"Wrote PDF:\n{args.out_pdf}")


if __name__ == "__main__":
    main()
