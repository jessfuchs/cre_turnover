#!/usr/bin/env python3

"""
Recurrence and robustness of Tier-1 CRE candidates.

x:
    number of analysed clades showing Tier-1 support

y:
    percentage of sensitivity scenarios in which the candidate
    remains classified as Tier 1

Symbols:
    triangle = Azteca focal Tier-1 candidate
    orange circle = other focal Tier-1-supported candidate
    grey circle = secondary-only Tier-1-supported candidate

Labels:
    - all Azteca focal candidates are labeled
    - all focal-supported candidates around 66.7% retention are labeled
    - selected recurrent and isolated candidates are labeled
    - selected dense clusters are manually arranged
    - thin leader lines connect labels to points
    - bold labels indicate stable candidate priority

Input
-----
downstream_analyses/candidate_analysis/results/tables/
    tier1_candidates_prioritized_with_sensitivity.tsv

Output
------
downstream_analyses/climate_analysis/results/figures/
    candidate_recurrence_robustness.png
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
DOWNSTREAM_DIR = ANALYSIS_DIR.parent
PROJECT_DIR = DOWNSTREAM_DIR.parent

INPUT_FILE = (
    PROJECT_DIR
    / "downstream_analyses"
    / "candidate_analysis"
    / "results"
    / "tables"
    / "tier1_candidates_prioritized_with_sensitivity.tsv"
)

FIG_DIR = (
    ANALYSIS_DIR
    / "results"
    / "figures"
)

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUT_PNG = (
    FIG_DIR
    / "candidate_recurrence_robustness.png"
)


# ============================================================
# Plot settings
# ============================================================

FOCAL_COLOR = "#E69F00"
SECONDARY_COLOR = "#BDBDBD"

EDGE_COLOR = "#222222"
LEADER_COLOR = "#777777"

AZTECA_CLADE = "azteca_affinis_miranda_group"


# ============================================================
# Manual label positions
#
# (-x, +y) = upper left
# (+x, +y) = upper right
# (-x, -y) = lower left
# (+x, -y) = lower right
# (0, -y)  = vertically below
# ============================================================


MANUAL_LABEL_OFFSETS = {

    # Azteca focal candidates
    "DMEL_CRE_00080": (12, 12),
    "DMEL_CRE_00309": (10, 12),
    "DMEL_CRE_00310": (10, -16),
    "DMEL_CRE_00330": (10, 12),
    "DMEL_CRE_00331": (10, 12),

    # 100% cluster around x ~ 2
    "DMEL_CRE_00052": (-10, 12),
    "DMEL_CRE_00065": (-12, 12),
    "DMEL_CRE_00078": (-10, -14),
    "DMEL_CRE_00275": (0, -18),

    # 66.7% focal-supported cluster: four-way fan
    "DMEL_CRE_00017": (-14, 10),   # upper left
    "DMEL_CRE_00121": (14, 10),    # upper right
    "DMEL_CRE_00122": (-14, -10),  # lower left
    "DMEL_CRE_00306": (14, -10),   # lower right

    # recurrent candidates
    "DMEL_CRE_00151": (-12, 12),
    "DMEL_CRE_00152": (-12, 12),

    # isolated candidates
    "DMEL_CRE_00255": (10, 10),
    "DMEL_CRE_00175": (10, 12),

    # ~33% cluster
    "DMEL_CRE_00154": (-12, 12),
    "DMEL_CRE_00277": (12, 12),
    "DMEL_CRE_00186": (10, -12),
}


FORCE_LABEL_IDS = set(
    MANUAL_LABEL_OFFSETS.keys()
)


# ============================================================
# Fallback offsets
# ============================================================

DEFAULT_OFFSETS = [
    (8, 8),
    (8, -10),
    (-8, 8),
    (-8, -10),
    (12, 4),
    (-12, 4),
]


# ============================================================
# Dense-cluster detection
# ============================================================

DENSE_X_RADIUS = 0.16
DENSE_Y_RADIUS = 2.0
DENSE_MIN_POINTS = 3


# ============================================================
# Helpers
# ============================================================

def short_id(cre_id):
    return str(cre_id).replace(
        "DMEL_CRE_",
        "",
    )


def contains_azteca(value):

    entries = [
        item.strip()
        for item in str(value).split("|")
        if item.strip()
    ]

    return AZTECA_CLADE in entries


def local_neighbor_count(
    row,
    data,
):

    dx = np.abs(
        data["x_plot"]
        - row["x_plot"]
    )

    dy = np.abs(
        data["percent_candidate_retained"]
        - row["percent_candidate_retained"]
    )

    return int(
        (
            (dx <= DENSE_X_RADIUS)
            & (dy <= DENSE_Y_RADIUS)
        ).sum()
    )


def choose_default_offset(index):

    return DEFAULT_OFFSETS[
        index % len(DEFAULT_OFFSETS)
    ]


# ============================================================
# Input check
# ============================================================

if not INPUT_FILE.exists():

    raise SystemExit(
        f"ERROR: input file not found:\n{INPUT_FILE}"
    )


# ============================================================
# Load data
# ============================================================

df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
)


required_columns = {
    "dmel_cre_id",
    "n_total_tier1_clades",
    "has_focal_support",
    "percent_candidate_retained",
    "priority_stability",
    "focal_tier1_clades",
}


missing = (
    required_columns
    - set(df.columns)
)


if missing:

    raise SystemExit(
        "ERROR: input table missing columns:\n"
        + "\n".join(
            sorted(missing)
        )
    )


# ============================================================
# Numeric conversion
# ============================================================

df["n_total_tier1_clades"] = pd.to_numeric(
    df["n_total_tier1_clades"],
    errors="coerce",
)

df["percent_candidate_retained"] = pd.to_numeric(
    df["percent_candidate_retained"],
    errors="coerce",
)

df = df.dropna(
    subset=[
        "n_total_tier1_clades",
        "percent_candidate_retained",
    ]
).copy()


# ============================================================
# Candidate classes
# ============================================================

df["is_focal"] = (
    df["has_focal_support"]
    .astype(str)
    .str.lower()
    == "yes"
)

df["is_azteca_focal"] = (
    df["focal_tier1_clades"]
    .apply(
        contains_azteca
    )
)


# ============================================================
# Deterministic x jitter
# ============================================================

df = df.sort_values(
    [
        "n_total_tier1_clades",
        "percent_candidate_retained",
        "dmel_cre_id",
    ]
).copy()

df["x_plot"] = (
    df["n_total_tier1_clades"]
    .astype(float)
)


for n_clades in sorted(
    df["n_total_tier1_clades"]
    .unique()
):

    idx = df.index[
        df["n_total_tier1_clades"]
        == n_clades
    ]

    n = len(idx)

    if n == 1:

        offsets = np.array([
            0.0
        ])

    else:

        offsets = np.linspace(
            -0.13,
            0.13,
            n,
        )

    df.loc[
        idx,
        "x_plot",
    ] = (
        float(n_clades)
        + offsets
    )


# ============================================================
# Dense-cluster detection
# ============================================================

df["local_neighbor_count"] = df.apply(
    lambda row:
        local_neighbor_count(
            row,
            df,
        ),
    axis=1,
)

df["is_dense_cluster"] = (
    df["local_neighbor_count"]
    >= DENSE_MIN_POINTS
)


# ============================================================
# Decide which candidates receive labels
# ============================================================

def should_label(row):

    cre_id = row["dmel_cre_id"]

    # Explicitly selected labels
    if cre_id in FORCE_LABEL_IDS:
        return True

    # All Azteca focal candidates
    if row["is_azteca_focal"]:
        return True

    # All focal-supported candidates around 66.7% retention
    if (
        row["is_focal"]
        and abs(
            row["percent_candidate_retained"]
            - 66.6667
        ) < 0.5
    ):
        return True

    # Recurrent candidates outside dense clusters
    if (
        row["n_total_tier1_clades"] >= 2
        and not row["is_dense_cluster"]
    ):
        return True

    # Isolated lower-retention candidates
    if (
        row["percent_candidate_retained"] < 80
        and not row["is_dense_cluster"]
    ):
        return True

    return False


df["label_this"] = df.apply(
    should_label,
    axis=1,
)


# ============================================================
# Split plotting groups
# ============================================================

secondary = df[
    ~df["is_focal"]
].copy()

focal_non_azteca = df[
    df["is_focal"]
    & ~df["is_azteca_focal"]
].copy()

azteca_focal = df[
    df["is_azteca_focal"]
].copy()


# ============================================================
# Plot
# ============================================================

fig, ax = plt.subplots(
    figsize=(
        10.8,
        6.6,
    )
)


# ------------------------------------------------------------
# Secondary candidates
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Other focal-supported candidates
# ------------------------------------------------------------

ax.scatter(
    focal_non_azteca["x_plot"],
    focal_non_azteca["percent_candidate_retained"],
    s=88,
    facecolor=FOCAL_COLOR,
    edgecolor=EDGE_COLOR,
    linewidth=0.75,
    alpha=0.96,
    marker="o",
    zorder=3,
)


# ------------------------------------------------------------
# Azteca focal candidates
# ------------------------------------------------------------

ax.scatter(
    azteca_focal["x_plot"],
    azteca_focal["percent_candidate_retained"],
    s=125,
    facecolor=FOCAL_COLOR,
    edgecolor=EDGE_COLOR,
    linewidth=0.9,
    alpha=0.98,
    marker="^",
    zorder=4,
)


# ============================================================
# Reference lines
# ============================================================

ax.axhline(
    100,
    color="#666666",
    linestyle=":",
    linewidth=1.0,
    zorder=1,
)

ax.axhline(
    66.6667,
    color="#AAAAAA",
    linestyle=":",
    linewidth=0.9,
    zorder=1,
)


# ============================================================
# Labels
# ============================================================

label_df = df.loc[
    df["label_this"]
].copy()


# priority:
# Azteca -> focal -> secondary
label_df["_label_priority"] = 2

label_df.loc[
    label_df["is_focal"],
    "_label_priority",
] = 1

label_df.loc[
    label_df["is_azteca_focal"],
    "_label_priority",
] = 0


label_df = label_df.sort_values(
    [
        "_label_priority",
        "n_total_tier1_clades",
        "percent_candidate_retained",
        "dmel_cre_id",
    ],
    ascending=[
        True,
        True,
        False,
        True,
    ],
)


generic_counter = 0


for _, row in label_df.iterrows():

    cre_id = row["dmel_cre_id"]


    # --------------------------------------------------------
    # Offset
    # --------------------------------------------------------

    if cre_id in MANUAL_LABEL_OFFSETS:

        dx, dy = (
            MANUAL_LABEL_OFFSETS[
                cre_id
            ]
        )

    else:

        dx, dy = (
            choose_default_offset(
                generic_counter
            )
        )

        generic_counter += 1


    # --------------------------------------------------------
    # Alignment
    # --------------------------------------------------------

    if dx > 0:
        ha = "left"

    elif dx < 0:
        ha = "right"

    else:
        ha = "center"


    if dy > 0:
        va = "bottom"

    elif dy < 0:
        va = "top"

    else:
        va = "center"


    # --------------------------------------------------------
    # Annotation
    # --------------------------------------------------------

    ax.annotate(
        short_id(cre_id),
        xy=(
            row["x_plot"],
            row["percent_candidate_retained"],
        ),
        xytext=(
            dx,
            dy,
        ),
        textcoords="offset points",
        fontsize=8.7,
        ha=ha,
        va=va,
        fontweight=(
            "bold"
            if str(
                row["priority_stability"]
            ).lower()
            == "stable"
            else "normal"
        ),
        arrowprops=dict(
            arrowstyle="-",
            linewidth=0.55,
            color=LEADER_COLOR,
            shrinkA=2,
            shrinkB=2,
        ),
        zorder=5,
    )


# ============================================================
# Axes
# ============================================================

max_clades = int(
    df["n_total_tier1_clades"]
    .max()
)

ax.set_xticks(
    range(
        1,
        max_clades + 1,
    )
)

ax.set_xlim(
    0.65,
    max_clades + 0.42,
)

ax.set_ylim(
    25,
    104,
)

ax.set_xlabel(
    "Number of analysed clades showing Tier-1 support",
    fontsize=11,
)

ax.set_ylabel(
    "Candidate retention across sensitivity scenarios (%)",
    fontsize=11,
)

ax.set_title(
    "Recurrence and robustness of Tier-1 CRE candidates",
    fontsize=14,
    pad=12,
)


# ============================================================
# Grid / spines
# ============================================================

ax.grid(
    axis="both",
    linestyle=":",
    linewidth=0.6,
    alpha=0.35,
)

ax.set_axisbelow(
    True
)

ax.spines[
    "top"
].set_visible(
    False
)

ax.spines[
    "right"
].set_visible(
    False
)


# ============================================================
# Legend
# ============================================================

legend_handles = [
    Line2D(
        [0],
        [0],
        marker="^",
        linestyle="none",
        markerfacecolor=FOCAL_COLOR,
        markeredgecolor=EDGE_COLOR,
        markersize=8,
        label="Azteca focal Tier-1 candidate",
    ),

    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        markerfacecolor=FOCAL_COLOR,
        markeredgecolor=EDGE_COLOR,
        markersize=8,
        label="Other focal Tier-1 support",
    ),

    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        markerfacecolor=SECONDARY_COLOR,
        markeredgecolor=EDGE_COLOR,
        markersize=8,
        label="Secondary Tier-1 support",
    ),
]


ax.legend(
    handles=legend_handles,
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(
        1.01,
        1.00,
    ),
    borderaxespad=0.0,
    fontsize=10,
)


# ============================================================
# Layout
# ============================================================

plt.tight_layout()

plt.subplots_adjust(
    right=0.77
)


# ============================================================
# Save
# ============================================================

fig.savefig(
    OUT_PNG,
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)


# ============================================================
# Console summary
# ============================================================

print(
    f"Candidates plotted: {len(df)}"
)

print(
    f"Candidates labelled: "
    f"{int(df['label_this'].sum())}"
)

print(
    f"Azteca focal candidates: "
    f"{int(df['is_azteca_focal'].sum())}"
)

print()

print(
    f"Wrote:\n{OUT_PNG}"
)
