#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
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

RESULTS = (
    SENS_DIR
    / "results"
)

MATRIX_FILE = (
    RESULTS
    / "candidate_sensitivity_retention_matrix.tsv"
)

SUMMARY_FILE = (
    RESULTS
    / "candidate_sensitivity_summary.tsv"
)

FIG_DIR = (
    RESULTS
    / "figures"
)

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Scenario order
# ============================================================

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


SCENARIO_LABELS = [

    "0.25 / 12 kb",
    "0.25 / 24 kb",
    "0.25 / 48 kb",

    "0.50 / 12 kb",
    "0.50 / 24 kb",
    "0.50 / 48 kb",

    "0.75 / 12 kb",
    "0.75 / 24 kb",
    "0.75 / 48 kb",
]


# ============================================================
# Load
# ============================================================

matrix = pd.read_csv(
    MATRIX_FILE,
    sep="\t",
)


summary = pd.read_csv(
    SUMMARY_FILE,
    sep="\t",
)


# ============================================================
# Merge ranking information
# ============================================================

summary_keep = summary[
    [
        "dmel_cre_id",
        "percent_candidate_retained",
        "primary_candidate_priority",
    ]
].copy()


df = matrix.merge(
    summary_keep,
    on="dmel_cre_id",
    how="left",
    validate="one_to_one",
)


# ============================================================
# Priority rank
# ============================================================

priority_order = {

    "focal_recurrent":
        1,

    "focal_plus_secondary_recurrent":
        2,

    "secondary_recurrent":
        3,

    "focal_single":
        4,

    "secondary_single":
        5,
}


df[
    "priority_rank"
] = df[
    "primary_candidate_priority"
].map(
    priority_order
)


# ============================================================
# Sort candidates
# ============================================================

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
    )
    .reset_index(
        drop=True
    )
)


# ============================================================
# yes/no -> numeric
# ============================================================

values = df[
    SCENARIOS
].replace(
    {
        "yes": 1,
        "no": 0,
    }
).astype(
    int
).to_numpy()


# ============================================================
# Plot
# ============================================================

fig_height = max(
    8,
    len(df) * 0.23,
)


fig, ax = plt.subplots(
    figsize=(
        9,
        fig_height,
    )
)


cmap = ListedColormap(
    [
        "#C62828",
        "#2E7D32",
    ]
)


ax.imshow(
    values,
    aspect="auto",
    interpolation="nearest",
    cmap=cmap,
    vmin=0,
    vmax=1,
)


# ============================================================
# Labels
# ============================================================

ax.set_xticks(
    range(
        len(SCENARIOS)
    )
)

ax.set_xticklabels(
    SCENARIO_LABELS,
    rotation=45,
    ha="right",
)


ax.set_yticks(
    range(
        len(df)
    )
)

ax.set_yticklabels(
    df[
        "dmel_cre_id"
    ],
    fontsize=8,
)


ax.set_xlabel(
    "Overlap threshold / local distance threshold"
)

ax.set_ylabel(
    "Tier-1 CRE candidate"
)


ax.set_title(
    "Retention of Tier-1 candidates across sensitivity scenarios"
)


# ============================================================
# Highlight primary scenario
# ============================================================

primary_index = SCENARIOS.index(
    "ov050_dist24000"
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


# ============================================================
# Group separators
# ============================================================

for position in [
    2.5,
    5.5,
]:

    ax.axvline(
        position,
        color="white",
        linewidth=2,
    )


# ============================================================
# Legend
# ============================================================

from matplotlib.patches import Patch


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


# ============================================================
# Save
# ============================================================

for extension in [
    "png",
    "pdf",
    "svg",
]:

    out = (
        FIG_DIR
        / f"tier1_sensitivity_matrix.{extension}"
    )

    fig.savefig(
        out,
        dpi=300,
        bbox_inches="tight",
    )

    print(
        f"Wrote: {out}"
    )


plt.close(fig)
