#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


# ============================================================
# Paths
# ============================================================

PROJECT = (
    Path.home()
    / "cre_turnover"
    / "project"
)

SENS_DIR = (
    PROJECT
    / "downstream_analyses"
    / "sensitivity_analysis"
)

INPUT = (
    SENS_DIR
    / "results"
    / "candidate_sensitivity_summary.tsv"
)

FIG_DIR = (
    SENS_DIR
    / "results"
    / "figures"
)

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Candidate priority setup
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
# Robustness colors
# ============================================================

ROBUSTNESS_COLORS = {
    "robust": "#2E7D32",
    "moderate": "#E6A700",
    "sensitive": "#C62828",
}


# ============================================================
# Load input
# ============================================================

df = pd.read_csv(
    INPUT,
    sep="\t",
)

required = {
    "dmel_cre_id",
    "primary_candidate_priority",
    "percent_candidate_retained",
    "percent_same_priority_as_primary",
    "robust_priority_all_9",
}

missing = required - set(df.columns)

if missing:
    raise SystemExit(
        "ERROR: missing required columns:\n"
        + "\n".join(sorted(missing))
    )


# ============================================================
# Robustness class
# ============================================================

def classify_robustness(value):
    if value >= 99.999:
        return "robust"
    if value >= 66.6:
        return "moderate"
    return "sensitive"


df["robustness_class"] = df[
    "percent_candidate_retained"
].apply(classify_robustness)


# ============================================================
# Priority rank
# ============================================================

df["priority_rank"] = df[
    "primary_candidate_priority"
].map(PRIORITY_ORDER)

if df["priority_rank"].isna().any():
    unknown = sorted(
        df.loc[
            df["priority_rank"].isna(),
            "primary_candidate_priority",
        ].unique()
    )
    raise SystemExit(
        "ERROR: unknown priority categories:\n"
        + "\n".join(unknown)
    )


# ============================================================
# Sort candidates
#
# Order:
# 1. highest retention first
# 2. strongest priority class first
# 3. most priority-stable first
# 4. CRE ID
# ============================================================

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
    )
    .reset_index(drop=True)
)

# strongest / most robust at top
df["y"] = list(
    range(len(df), 0, -1)
)


# ============================================================
# Figure and axes
# Use a dedicated legend axis so legends never get clipped.
# ============================================================

height = max(
    9,
    len(df) * 0.26,
)

fig = plt.figure(
    figsize=(13.5, height)
)

gs = fig.add_gridspec(
    nrows=1,
    ncols=2,
    width_ratios=[4.8, 1.8],
    wspace=0.03,
)

ax = fig.add_subplot(gs[0, 0])
ax_leg = fig.add_subplot(gs[0, 1])

ax_leg.axis("off")


# ============================================================
# Background robustness zones
# ============================================================

ax.axvspan(
    25,
    66.6,
    color="#FCE8E6",
    alpha=0.45,
    zorder=0,
)

ax.axvspan(
    66.6,
    99.999,
    color="#FFF5D6",
    alpha=0.45,
    zorder=0,
)

ax.axvspan(
    99.999,
    103,
    color="#E6F4EA",
    alpha=0.55,
    zorder=0,
)

ax.axvline(
    66.6,
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


# ============================================================
# Plot points
#
# Shape = candidate priority
# Fill color = robustness class
# Black outline = same priority in all 9 scenarios
# ============================================================

MARKER_SIZE = 95

for row in df.itertuples():
    priority = row.primary_candidate_priority
    robustness = row.robustness_class
    priority_stable = (
        str(row.robust_priority_all_9).lower() == "yes"
    )

    ax.scatter(
        row.percent_candidate_retained,
        row.y,
        s=MARKER_SIZE,
        marker=MARKERS[priority],
        facecolor=ROBUSTNESS_COLORS[robustness],
        edgecolor="black" if priority_stable else "none",
        linewidth=1.5 if priority_stable else 0,
        zorder=3,
    )


# ============================================================
# Axis formatting
# ============================================================

ax.set_yticks(df["y"])
ax.set_yticklabels(
    df["dmel_cre_id"],
    fontsize=8,
)

ax.set_xlim(25, 103)

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


# ============================================================
# Zone labels
# ============================================================

top = df["y"].max() + 1.6

ax.text(
    45,
    top,
    "sensitive",
    ha="center",
    va="center",
    fontsize=10,
    fontweight="bold",
    color=ROBUSTNESS_COLORS["sensitive"],
)

ax.text(
    83,
    top,
    "moderate",
    ha="center",
    va="center",
    fontsize=10,
    fontweight="bold",
    color=ROBUSTNESS_COLORS["moderate"],
)

ax.text(
    100.6,
    top,
    "robust",
    ha="center",
    va="center",
    fontsize=10,
    fontweight="bold",
    color=ROBUSTNESS_COLORS["robust"],
)

ax.set_ylim(0, top + 0.8)


# ============================================================
# Legends in dedicated legend axis
# ============================================================

# ---- priority legend
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
    bbox_to_anchor=(0.0, 1.0),
    frameon=False,
    borderaxespad=0.0,
    labelspacing=0.55,
    handletextpad=0.7,
    fontsize=9,
    title_fontsize=10,
)
ax_leg.add_artist(legend1)


# ---- robustness legend
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
    for label, color in ROBUSTNESS_COLORS.items()
]

legend2 = ax_leg.legend(
    handles=robustness_handles,
    title="Candidate robustness",
    loc="upper left",
    bbox_to_anchor=(0.0, 0.62),
    frameon=False,
    borderaxespad=0.0,
    labelspacing=0.55,
    handletextpad=0.7,
    fontsize=9,
    title_fontsize=10,
)
ax_leg.add_artist(legend2)


# ---- priority stability legend
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
        label="Same priority in all 9",
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
    bbox_to_anchor=(0.0, 0.34),
    frameon=False,
    borderaxespad=0.0,
    labelspacing=0.55,
    handletextpad=0.7,
    fontsize=9,
    title_fontsize=10,
)
ax_leg.add_artist(legend3)


# left-align legend contents
for lg in [legend1, legend2, legend3]:
    lg._legend_box.align = "left"


# ============================================================
# Final layout
# ============================================================

fig.subplots_adjust(
    left=0.18,
    right=0.98,
    top=0.95,
    bottom=0.06,
)


# ============================================================
# Save
# Since the legends are inside the figure canvas, no clipping.
# ============================================================

for extension in ["png"]:
    out = (
        FIG_DIR
        / f"tier1_candidate_robustness_ranked.{extension}"
    )

    fig.savefig(
        out,
        dpi=300,
    )

    print(f"Wrote: {out}")

plt.close(fig)
