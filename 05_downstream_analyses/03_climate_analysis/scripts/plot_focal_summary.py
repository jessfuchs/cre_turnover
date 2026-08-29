#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPT_DIR.parent

INPUT_FILE = ANALYSIS_DIR / "results" / "focal_tier1_enrichment.tsv"
FIG_DIR = ANALYSIS_DIR / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_PNG = FIG_DIR / "focal_tier1_summary.png"

TURNOVER_COLOR = "#E69F00"
PRESENT_COLOR = "#0072B2"
OTHER_COLOR = "#D9D9D9"
EDGE_COLOR = "#333333"

df = pd.read_csv(INPUT_FILE, sep="\t")

required = {
    "clade_display",
    "n_tier1_singletons",
    "n_focal_tier1",
    "n_focal_turnover",
    "n_focal_present",
    "observed_focal_share",
    "expected_focal_share",
}
missing = required - set(df.columns)
if missing:
    raise SystemExit(
        "Missing columns:\n" + "\n".join(sorted(missing))
    )

df["focal_turnover_share"] = (
    df["n_focal_turnover"] / df["n_tier1_singletons"]
)
df["focal_present_share"] = (
    df["n_focal_present"] / df["n_tier1_singletons"]
)
df["other_share"] = (
    1.0
    - df["focal_turnover_share"]
    - df["focal_present_share"]
)

fig, ax = plt.subplots(figsize=(10.2, 5.8))

y = range(len(df))
height = 0.68

turnover_pct = df["focal_turnover_share"] * 100
present_pct = df["focal_present_share"] * 100
other_pct = df["other_share"] * 100
observed_pct = df["observed_focal_share"] * 100
expected_pct = df["expected_focal_share"] * 100

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
    present_pct,
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

for i, value in enumerate(expected_pct):
    ax.plot(
        [value, value],
        [i - height / 2 - 0.08, i + height / 2 + 0.08],
        color=EDGE_COLOR,
        linewidth=2.0,
        zorder=5,
    )

for i, row in df.iterrows():
    label = f"{int(row['n_focal_tier1'])}/{int(row['n_tier1_singletons'])}"
    ax.text(
        101.2,
        i,
        label,
        va="center",
        ha="left",
        fontsize=10,
    )

ax.set_yticks(list(y))
ax.set_yticklabels(df["clade_display"], fontsize=11)
ax.invert_yaxis()
ax.set_xlim(0, 108)

ax.set_xlabel(
    "Share of clade-wide singleton Tier-1 contrasts (%)",
    fontsize=11,
)

ax.set_title(
    "Concentration of lineage-specific Tier-1 contrasts in focal species",
    fontsize=14,
    pad=12,
)

ax.grid(axis="x", linestyle=":", linewidth=0.6, alpha=0.35)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

handles = [
    Line2D([0], [0], color=TURNOVER_COLOR, linewidth=10, label="Focal turnover"),
    Line2D([0], [0], color=PRESENT_COLOR, linewidth=10, label="Focal present"),
    Line2D([0], [0], color=OTHER_COLOR, linewidth=10, label="Other singleton"),
    Line2D([0], [0], color=EDGE_COLOR, linewidth=2.5, label="Expected focal share"),
]

ax.legend(
    handles=handles,
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(1.01, 1.00),
    borderaxespad=0.0,
    fontsize=10,
)

plt.tight_layout()
plt.subplots_adjust(right=0.79)

fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
plt.close(fig)

print("Wrote:")
print(OUT_PNG)
