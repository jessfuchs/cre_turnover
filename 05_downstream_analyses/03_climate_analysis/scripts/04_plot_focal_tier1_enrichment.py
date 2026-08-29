#!/usr/bin/env python3

"""
Plot observed focal Tier-1 share versus the clade-specific null expectation.

Observed:
    n_focal_tier1 / n_tier1_singletons

Expected:
    1 / n_species

The observed proportion is shown with an exact 95% binomial confidence
interval.

Output
-------
results/figures/focal_tier1_enrichment.png
"""

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# ============================================================
# Paths
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPT_DIR.parent

RESULTS_DIR = ANALYSIS_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


INPUT_FILE = (
    RESULTS_DIR
    / "focal_tier1_enrichment.tsv"
)

OUT_PNG = (
    FIG_DIR
    / "focal_tier1_enrichment.png"
)


# ============================================================
# Load
# ============================================================

if not INPUT_FILE.exists():

    raise SystemExit(
        f"ERROR: input file missing:\n{INPUT_FILE}"
    )


df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
)


required = {
    "clade_display",
    "observed_focal_share",
    "expected_focal_share",
    "ci95_low",
    "ci95_high",
    "n_focal_tier1",
    "n_tier1_singletons",
}


missing = (
    required
    - set(df.columns)
)


if missing:
    raise SystemExit(
        "ERROR: enrichment table missing columns:\n"
        + "\n".join(
            sorted(missing)
        )
    )


# ============================================================
# Values
# ============================================================

observed = (
    df[
        "observed_focal_share"
    ]
    .astype(float)
    .to_numpy()
)


expected = (
    df[
        "expected_focal_share"
    ]
    .astype(float)
    .to_numpy()
)


ci_low = (
    df[
        "ci95_low"
    ]
    .astype(float)
    .to_numpy()
)


ci_high = (
    df[
        "ci95_high"
    ]
    .astype(float)
    .to_numpy()
)


y = np.arange(
    len(df)
)


# ============================================================
# Plot
# ============================================================

fig, ax = plt.subplots(
    figsize=(
        8.2,
        5.6,
    )
)


# ------------------------------------------------------------
# Connecting line from expected to observed
# ------------------------------------------------------------

for i in range(
    len(df)
):

    ax.plot(
        [
            expected[i] * 100,
            observed[i] * 100,
        ],
        [
            y[i],
            y[i],
        ],
        color="#B0B0B0",
        linewidth=1.2,
        zorder=1,
    )


# ------------------------------------------------------------
# Confidence intervals
# ------------------------------------------------------------

xerr = np.vstack([
    (
        observed
        - ci_low
    ) * 100,

    (
        ci_high
        - observed
    ) * 100,
])


ax.errorbar(
    observed * 100,
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


# ------------------------------------------------------------
# Expected share
# ------------------------------------------------------------

ax.scatter(
    expected * 100,
    y,
    s=55,
    facecolor="white",
    edgecolor="#222222",
    linewidth=1.2,
    marker="o",
    zorder=4,
)


# ============================================================
# Counts
# ============================================================

for i, row in df.iterrows():

    label = (
        f"{int(row['n_focal_tier1'])}/"
        f"{int(row['n_tier1_singletons'])}"
    )


    x_text = (
        row[
            "observed_focal_share"
        ]
        * 100
        + 2
    )


    ax.text(
        x_text,
        i,
        label,
        va="center",
        ha="left",
        fontsize=8,
    )


# ============================================================
# Formatting
# ============================================================

ax.set_yticks(
    y
)


ax.set_yticklabels(
    df[
        "clade_display"
    ],
    fontsize=10,
)


ax.invert_yaxis()


ax.set_xlim(
    0,
    100,
)


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


ax.set_axisbelow(
    True
)


# ============================================================
# Legend
# ============================================================

legend_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        markerfacecolor="#E69F00",
        markeredgecolor="#222222",
        markersize=7,
        label="Observed focal share",
    ),

    Line2D(
        [0],
        [0],
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


# ============================================================
# Save
# ============================================================

fig.savefig(
    OUT_PNG,
    dpi=300,
    bbox_inches="tight",
)

plt.close(
    fig
)


print()
print("=" * 72)
print("FOCAL TIER-1 ENRICHMENT PLOT")
print("=" * 72)

print(
    f"Wrote:\n{OUT_PNG}"
)
