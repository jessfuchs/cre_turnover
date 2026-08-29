#!/usr/bin/env python3

"""
Phylogeny + species labels + CRE-state heatmap + climate strip.

Layout:
    phylogeny | species labels | climate | CRE states | legends

Important:
- Branch lengths are taken unchanged from the Newick tree.
- Species labels are drawn in their own panel and therefore do not
  influence the phylogenetic x-axis.
- Tree tips, labels, CRE rows, and climate rows share exactly the
  same species order and y coordinates.
- D. melanogaster is shown only as a reference row and is NOT
  included as a target species in CRE statistics.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch, Rectangle

from Bio import Phylo


# ============================================================
# Paths
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
CLASS_DIR = SCRIPT_DIR.parent
PROJECT_DIR = CLASS_DIR.parent

PHYLO_DIR = CLASS_DIR / "phylogeny"

RESULTS_DIR = CLASS_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

# Pruned 40-species phylogeny
TREE_FILE = (
    PHYLO_DIR
    / "results"
    / "301Fly_HOG_UCLDtree_40species.nw"
)

# Central species manifest
MANIFEST_FILE = (
    PROJECT_DIR
    / "external_scrmshaw"
    / "combined_manifest.tsv"
)

# Climatic-zone annotations
TRAITS_FILE = (
    PHYLO_DIR
    / "data"
    / "species_traits.tsv"
)

# Final CRE-state matrix
MATRIX_FILE = (
    RESULTS_DIR
    / "cre_turnover_matrix.tsv"
)

OUT_PNG = (
    FIG_DIR
    / "phylogeny_CRE_heatmap_climate.png"
)

# ============================================================
# Constants
# ============================================================

DMEL_SLUG = "d_melanogaster"
DMEL_TREE_NAME = "DROSOPHILA_MELANOGASTER"


# ============================================================
# CRE-state palette
# ============================================================

STATE_COLORS = {
    "present": "#0072B2",
    "turnover_candidate": "#E69F00",
    "no_detected_CRE": "#E5E5E5",
    "uncertain": "#9E9E9E",
    "reference": "#F5F5F5",
}

STATE_LABELS = {
    "present": "Positional match",
    "turnover_candidate": "Turnover candidate",
    "no_detected_CRE": "No CRE detected",
    "uncertain": "Uncertain",
    "reference": "Reference species",
}

STATE_TO_NUM = {
    "no_detected_CRE": 0,
    "turnover_candidate": 1,
    "present": 2,
    "uncertain": 3,
    "reference": 4,
}


# ============================================================
# Climate palette
# ============================================================

CLIMATE_COLORS = {
    "TROP": "#E41A1C",
    "ARID": "#FBC02D",
    "TEMP": "#009E73",
    "BORE": "#0072B2",
}

CLIMATE_LABELS = {
    "TROP": "Tropical",
    "ARID": "Arid",
    "TEMP": "Temperate",
    "BORE": "Boreal",
}

CLIMATE_TO_NUM = {
    "TROP": 0,
    "ARID": 1,
    "TEMP": 2,
    "BORE": 3,
}


# ============================================================
# Helper functions
# ============================================================

def require_file(path):
    if not path.exists():
        raise SystemExit(
            f"ERROR: required file not found:\n{path}"
        )


def short_species_name(
    slug,
    manifest,
):
    """
    Return abbreviated scientific species name.

    Species with externally obtained SCRMshaw predictions
    are marked with an asterisk. D. melanogaster is excluded
    because it is the reference species.
    """

    hit = manifest.loc[
        manifest["slug"] == slug
    ]

    if len(hit) != 1:
        return slug

    row = hit.iloc[0]

    name = row["species"]
    source = row["source"]

    if name.startswith("Drosophila "):
        label = (
            "D. "
            + name.split(" ", 1)[1]
        )

    elif name.startswith("Zaprionus "):
        label = (
            "Z. "
            + name.split(" ", 1)[1]
        )

    else:
        label = name

    if (
        source == "external"
        and slug != DMEL_SLUG
    ):
        label += "*"

    return label


def tree_name_from_species(species_name):
    """
    Convert scientific species name to the naming convention
    used in the published phylogeny.

    Example
    -------
    Drosophila virilis
        -> DROSOPHILA_VIRILIS
    """

    return (
        str(species_name)
        .strip()
        .upper()
        .replace(" ", "_")
    )

# ============================================================
# Input checks
# ============================================================

for path in [
    TREE_FILE,
    MANIFEST_FILE,
    TRAITS_FILE,
    MATRIX_FILE,
]:
    require_file(path)


# ============================================================
# Load species manifest
# ============================================================

manifest = pd.read_csv(
    MANIFEST_FILE,
    sep="\t",
    dtype=str,
).fillna("")

required = {
    "slug",
    "species",
    "source",
}

missing = (
    required
    - set(manifest.columns)
)

if missing:
    raise SystemExit(
        "ERROR: combined_manifest.tsv missing columns:\n"
        + "\n".join(
            sorted(missing)
        )
    )


for column in [
    "slug",
    "species",
    "source",
]:
    manifest[column] = (
        manifest[column]
        .astype(str)
        .str.strip()
    )


# Construct tree names deterministically.
manifest["tree_name"] = (
    manifest["species"]
    .map(tree_name_from_species)
)


# Slugs must be unique.
duplicated_slugs = (
    manifest.loc[
        manifest["slug"].duplicated(
            keep=False
        ),
        "slug",
    ]
    .unique()
    .tolist()
)

if duplicated_slugs:
    raise SystemExit(
        "ERROR: duplicate slugs in combined_manifest.tsv:\n"
        + "\n".join(
            sorted(duplicated_slugs)
        )
    )


# Inferred tree names must also be unique.
duplicated_tree_names = (
    manifest.loc[
        manifest["tree_name"].duplicated(
            keep=False
        ),
        "tree_name",
    ]
    .unique()
    .tolist()
)

if duplicated_tree_names:
    raise SystemExit(
        "ERROR: duplicate inferred tree names "
        "in combined_manifest.tsv:\n"
        + "\n".join(
            sorted(duplicated_tree_names)
        )
    )


tree_to_slug = dict(
    zip(
        manifest["tree_name"],
        manifest["slug"],
    )
)

slug_to_tree = dict(
    zip(
        manifest["slug"],
        manifest["tree_name"],
    )
)


# ============================================================
# Load climate traits
# ============================================================

traits = pd.read_csv(
    TRAITS_FILE,
    sep="\t",
    dtype=str,
).fillna("")

required = {
    "tree_name",
    "climatic_zone",
}

missing = required - set(traits.columns)

if missing:
    raise SystemExit(
        "ERROR: species_traits.tsv missing columns:\n"
        + "\n".join(sorted(missing))
    )


tree_to_climate = dict(
    zip(
        traits["tree_name"],
        traits["climatic_zone"],
    )
)


unknown_climate = sorted(
    set(traits["climatic_zone"])
    - set(CLIMATE_TO_NUM)
    - {""}
)

if unknown_climate:
    raise SystemExit(
        "ERROR: unknown climatic-zone values:\n"
        + "\n".join(unknown_climate)
    )


# ============================================================
# Load CRE matrix
# ============================================================

matrix = pd.read_csv(
    MATRIX_FILE,
    sep="\t",
    index_col=0,
)


if matrix.index.duplicated().any():

    duplicates = (
        matrix.index[
            matrix.index.duplicated()
        ]
        .unique()
        .tolist()
    )

    raise SystemExit(
        "ERROR: duplicated CRE IDs:\n"
        + "\n".join(duplicates)
    )


target_species = list(
    matrix.columns
)


if DMEL_SLUG in target_species:
    raise SystemExit(
        "ERROR: D. melanogaster occurs as a target species "
        "in cre_turnover_matrix.tsv."
        )


allowed_states = {
    "present",
    "turnover_candidate",
    "no_detected_CRE",
    "uncertain",
}

observed_states = set(
    pd.unique(
        matrix.values.ravel()
    )
)

unknown_states = sorted(
    observed_states
    - allowed_states
)

if unknown_states:
    raise SystemExit(
        "ERROR: unknown CRE states:\n"
        + "\n".join(unknown_states)
    )


# ============================================================
# Validate species metadata
# ============================================================

missing_species = [
    sp
    for sp in target_species
    if sp not in slug_to_tree
]

if missing_species:
    raise SystemExit(
        "ERROR: target species missing from combined_manifest.tsv:\n"
        + "\n".join(missing_species)
    )


# ============================================================
# Sort CREs
# ============================================================

n_present = (
    matrix == "present"
).sum(axis=1)

n_turnover = (
    matrix == "turnover_candidate"
).sum(axis=1)

n_uncertain = (
    matrix == "uncertain"
).sum(axis=1)

n_evaluable = (
    len(target_species)
    - n_uncertain
)


cre_sort = pd.DataFrame({
    "n_present": n_present,
    "n_turnover": n_turnover,
    "n_evaluable": n_evaluable,
})


cre_order = (
    cre_sort
    .sort_values(
        [
            "n_present",
            "n_turnover",
            "n_evaluable",
        ],
        ascending=[
            False,
            False,
            False,
        ],
        kind="mergesort",
    )
    .index
)


matrix = matrix.loc[
    cre_order
]


# ============================================================
# Load and prune tree
# ============================================================

tree = Phylo.read(
    TREE_FILE,
    "newick",
)


wanted_tree_names = {
    slug_to_tree[sp]
    for sp in target_species
}

wanted_tree_names.add(
    DMEL_TREE_NAME
)


all_tree_tips = {
    tip.name
    for tip in tree.get_terminals()
}


missing_from_tree = sorted(
    wanted_tree_names
    - all_tree_tips
)

if missing_from_tree:
    raise SystemExit(
        "ERROR: species missing from tree:\n"
        + "\n".join(missing_from_tree)
    )


for tip in list(
    tree.get_terminals()
):
    if tip.name not in wanted_tree_names:
        tree.prune(tip)


tips = tree.get_terminals()


if len(tips) != (
    len(target_species) + 1
):
    raise SystemExit(
        "ERROR: incorrect tip count after pruning."
    )


tree_tip_names = [
    tip.name
    for tip in tips
]


# ============================================================
# Species order
# ============================================================

species_order = [
    tree_to_slug[name]
    for name in tree_tip_names
]


if DMEL_SLUG not in species_order:
    raise SystemExit(
        "ERROR: Dmel missing from pruned tree."
    )


dmel_row = species_order.index(
    DMEL_SLUG
)


# ============================================================
# Species display labels
# ============================================================

species_labels = []

for sp in species_order:

    label = short_species_name(
        sp,
        manifest,
    )

    species_labels.append(
        label
    )

# ============================================================
# Climate order
# ============================================================

missing_climate = [
    name
    for name in tree_tip_names
    if (
        name not in tree_to_climate
        or not tree_to_climate[name]
    )
]

if missing_climate:
    raise SystemExit(
        "ERROR: missing climate annotation for:\n"
        + "\n".join(missing_climate)
    )


climate_order = [
    tree_to_climate[name]
    for name in tree_tip_names
]


climate_numeric = np.asarray(
    [
        CLIMATE_TO_NUM[zone]
        for zone in climate_order
    ],
    dtype=float,
).reshape(
    -1,
    1,
)


# ============================================================
# Build species x CRE matrix
# ============================================================

heatmap_rows = []


for sp in species_order:

    if sp == DMEL_SLUG:

        row = pd.Series(
            "reference",
            index=matrix.index,
            name=sp,
        )

    else:

        row = matrix[sp].copy()
        row.name = sp

    heatmap_rows.append(
        row
    )


heatmap = pd.DataFrame(
    heatmap_rows
)


numeric = heatmap.replace(
    STATE_TO_NUM
).astype(float)


# ============================================================
# Alignment QC
# ============================================================

print()
print("ROW ALIGNMENT QC")
print("----------------")

for i, (
    sp,
    climate,
) in enumerate(
    zip(
        species_order,
        climate_order,
    ),
    start=1,
):

    if sp == DMEL_SLUG:

        state_info = "REFERENCE"

    else:

        counts = (
            heatmap.loc[sp]
            .value_counts()
        )

        state_info = (
            f"present={counts.get('present', 0):3d}  "
            f"turnover={counts.get('turnover_candidate', 0):3d}  "
            f"no_detected={counts.get('no_detected_CRE', 0):3d}  "
            f"uncertain={counts.get('uncertain', 0):3d}"
        )

    print(
        f"{i:2d}  "
        f"{sp:18s}  "
        f"{climate:4s}  "
        f"{state_info}"
    )


# ============================================================
# Remove tree labels
#
# Species names are now drawn in a dedicated panel.
# ============================================================

for tip in tree.get_terminals():
    tip.name = ""


# ============================================================
# Colormaps
# ============================================================

state_cmap = ListedColormap([
    STATE_COLORS["no_detected_CRE"],
    STATE_COLORS["turnover_candidate"],
    STATE_COLORS["present"],
    STATE_COLORS["uncertain"],
    STATE_COLORS["reference"],
])


state_norm = BoundaryNorm(
    [
        -0.5,
        0.5,
        1.5,
        2.5,
        3.5,
        4.5,
    ],
    state_cmap.N,
)


climate_cmap = ListedColormap([
    CLIMATE_COLORS["TROP"],
    CLIMATE_COLORS["ARID"],
    CLIMATE_COLORS["TEMP"],
    CLIMATE_COLORS["BORE"],
])


climate_norm = BoundaryNorm(
    [
        -0.5,
        0.5,
        1.5,
        2.5,
        3.5,
    ],
    climate_cmap.N,
)


# ============================================================
# Figure layout
#
# tree | labels | heatmap | climate | legends
# ============================================================

fig = plt.figure(
    figsize=(21, 10.2)
)

gs = fig.add_gridspec(
    nrows=1,
    ncols=5,
    width_ratios=[
        1.30,   # phylogeny
        0.5,   # species labels
        0.10,   # climate strip
        4.25,   # CRE heatmap
        0.75,   # legends
    ],
    wspace=0.01,
)

ax_tree = fig.add_subplot(gs[0, 0])
ax_labels = fig.add_subplot(gs[0, 1])
ax_climate = fig.add_subplot(gs[0, 2])
ax_heat = fig.add_subplot(gs[0, 3])
ax_legend = fig.add_subplot(gs[0, 4])

# ============================================================
# Shared y coordinates
# ============================================================

n_species = len(
    species_order
)

n_cre = numeric.shape[1]


shared_ylim = (
    n_species + 0.5,
    0.5,
)


# ============================================================
# Draw phylogeny
# ============================================================

Phylo.draw(
    tree,
    axes=ax_tree,
    do_show=False,
    show_confidence=False,
)


ax_tree.set_ylim(
    *shared_ylim
)


ax_tree.set_title(
    "",
    fontsize=12,
    pad=10,
)


ax_tree.set_ylabel("")


# ------------------------------------------------------------
# Preserve original tree branch-length scale.
#
# Do NOT use the labels to determine x limits.
# ------------------------------------------------------------

depths = tree.depths()

terminal_depths = [
    depths[tip]
    for tip in tree.get_terminals()
]

max_tree_depth = max(
    terminal_depths
)


# Small visual margin only.
tree_xlim_max = (
    max_tree_depth * 1.03
)


ax_tree.set_xlim(
    0,
    tree_xlim_max,
)


# Major ticks up to 60 where applicable.
candidate_ticks = [
    0,
    20,
    40,
    60,
]

tree_ticks = [
    x
    for x in candidate_ticks
    if x <= tree_xlim_max
]

ax_tree.set_xticks(
    tree_ticks
)


ax_tree.set_xlabel(
    "Branch length",
    fontsize=10,
)


for spine in [
    "top",
    "right",
    "left",
]:
    ax_tree.spines[
        spine
    ].set_visible(False)


ax_tree.tick_params(
    axis="y",
    left=False,
    labelleft=False,
)


# ============================================================
# Dedicated species-label panel
# ============================================================

ax_labels.set_xlim(
    0,
    1,
)

ax_labels.set_ylim(
    *shared_ylim
)


for row_index, (
    sp,
    label,
) in enumerate(
    zip(
        species_order,
        species_labels,
    ),
    start=1,
):

    ax_labels.text(
        0.02,
        row_index,
        label,
        ha="left",
        va="center",
        fontsize=9,
        fontweight=(
            "bold"
            if sp == DMEL_SLUG
            else "normal"
        ),
    )

ax_labels.set_xticks([])
ax_labels.set_yticks([])

for spine in ax_labels.spines.values():
    spine.set_visible(False)


# ============================================================
# CRE heatmap
# ============================================================

ax_heat.imshow(
    numeric.values,
    aspect="auto",
    interpolation="nearest",
    cmap=state_cmap,
    norm=state_norm,
    origin="upper",
    extent=[
        -0.5,
        n_cre - 0.5,
        n_species + 0.5,
        0.5,
    ],
)


ax_heat.set_ylim(
    *shared_ylim
)


ax_heat.set_title(
    "",
    fontsize=12,
    pad=10,
)


ax_heat.set_yticks([])
ax_heat.set_ylabel("")


# ============================================================
# Sparse CRE labels
# ============================================================

MAX_CRE_LABELS = 113


if n_cre <= MAX_CRE_LABELS:

    tick_indices = np.arange(
        n_cre
    )

else:

    tick_indices = np.unique(
        np.linspace(
            0,
            n_cre - 1,
            MAX_CRE_LABELS,
            dtype=int,
        )
    )


ax_heat.set_xticks(
    tick_indices
)


ax_heat.set_xticklabels(
    [
        numeric.columns[i]
        for i in tick_indices
    ],
    rotation=90,
    fontsize=7,
)


ax_heat.set_xlabel(
    "D. melanogaster reference CREs (ordered by conservation)",
    fontsize=10,
    labelpad=16,
)


# ============================================================
# Dmel reference row
# ============================================================

dmel_y = (
    dmel_row + 1
)


reference_rectangle = Rectangle(
    (
        -0.5,
        dmel_y - 0.5,
    ),
    n_cre,
    1.0,
    fill=False,
    edgecolor="#555555",
    linewidth=1.2,
    zorder=10,
)


ax_heat.add_patch(
    reference_rectangle
)


# ============================================================
# Climate strip
# ============================================================

ax_climate.imshow(
    climate_numeric,
    aspect="auto",
    interpolation="nearest",
    cmap=climate_cmap,
    norm=climate_norm,
    origin="upper",
    extent=[
        -0.5,
        0.5,
        n_species + 0.5,
        0.5,
    ],
)


ax_climate.set_ylim(
    *shared_ylim
)


ax_climate.set_yticks([])


ax_climate.set_xticks(
    [0]
)


ax_climate.set_xticklabels(
    [""],
    rotation=90,
    fontsize=9,
)


ax_climate.tick_params(
    axis="both",
    length=0,
)


climate_reference_rectangle = Rectangle(
    (
        -0.5,
        dmel_y - 0.5,
    ),
    1.0,
    1.0,
    fill=False,
    edgecolor="#555555",
    linewidth=1.2,
    zorder=10,
)


ax_climate.add_patch(
    climate_reference_rectangle
)


# ============================================================
# Dedicated legend panel
# ============================================================

ax_legend.axis("off")


cre_handles = [
    Patch(
        facecolor=STATE_COLORS["present"],
        label=STATE_LABELS["present"],
    ),
    Patch(
        facecolor=STATE_COLORS["turnover_candidate"],
        label=STATE_LABELS["turnover_candidate"],
    ),
    Patch(
        facecolor=STATE_COLORS["no_detected_CRE"],
        label=STATE_LABELS["no_detected_CRE"],
    ),
    Patch(
        facecolor=STATE_COLORS["uncertain"],
        label=STATE_LABELS["uncertain"],
    ),
    Patch(
        facecolor=STATE_COLORS["reference"],
        edgecolor="#555555",
        label=STATE_LABELS["reference"],
    ),
]


climate_handles = [
    Patch(
        facecolor=CLIMATE_COLORS["TROP"],
        label=CLIMATE_LABELS["TROP"],
    ),
    Patch(
        facecolor=CLIMATE_COLORS["ARID"],
        label=CLIMATE_LABELS["ARID"],
    ),
    Patch(
        facecolor=CLIMATE_COLORS["TEMP"],
        label=CLIMATE_LABELS["TEMP"],
    ),
    Patch(
        facecolor=CLIMATE_COLORS["BORE"],
        label=CLIMATE_LABELS["BORE"],
    ),
]


cre_legend = ax_legend.legend(
    handles=cre_handles,
    title="CRE state",
    loc="upper left",
    bbox_to_anchor=(0.0, 1.0),
    frameon=False,
    fontsize=10,
    title_fontsize=11,
    borderaxespad=0,
    labelspacing=0.4,
    handlelength=1.6,
    handletextpad=0.5,
)

# Left-align legend title and entries.
cre_legend._legend_box.align = "left"

ax_legend.add_artist(
    cre_legend
)


climate_legend = ax_legend.legend(
    handles=climate_handles,
    title="Climatic zone",
    loc="upper left",

    # Smaller vertical gap below CRE-state legend.
    bbox_to_anchor=(0.0, 0.78),

    frameon=False,
    fontsize=10,
    title_fontsize=11,
    borderaxespad=0,
    labelspacing=0.4,
    handlelength=1.6,
    handletextpad=0.5,
)

climate_legend._legend_box.align = "left"


# ============================================================
# Main title
# ============================================================

fig.suptitle(
    "Phylogenetic distribution of predicted CRE "
    "conservation and turnover",
    fontsize=15,
    y=0.985,
)


# ============================================================
# Margins
# ============================================================

fig.subplots_adjust(
    left=0.035,
    right=0.985,
    bottom=0.14,
    top=0.92,
)


# ============================================================
# Save
# ============================================================

fig.savefig(
    OUT_PNG,
    dpi=300,
)


plt.close(
    fig
)


# ============================================================
# Console summary
# ============================================================

print()
print("PLOT SUMMARY")
print("------------")

print(
    f"Wrote: {OUT_PNG}"
)

print()

print(
    f"Target species: "
    f"{len(target_species)}"
)

print(
    f"Displayed species incl. Dmel: "
    f"{len(species_order)}"
)

print(
    f"Reference CREs: "
    f"{n_cre}"
)

print(
    f"Maximum root-to-tip branch length: "
    f"{max_tree_depth:.3f}"
)

print(
    f"Tree x-axis maximum used: "
    f"{tree_xlim_max:.3f}"
)

print()

print(
    "Branch lengths modified: NO"
)

print(
    "Species labels affect branch-length scaling: NO"
)

print(
    "Tree / labels / heatmap / climate synchronized: YES"
)

print(
    "Dmel included in target statistics: NO"
)

