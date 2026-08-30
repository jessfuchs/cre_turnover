#!/usr/bin/env python3

# ============================================================
# Plot focal-clade CRE-state heatmaps
#
# Purpose:
# Identify stringent lineage-specific CRE-state contrasts across six focal clades.
#
# Analysis:
# Strict singleton contrasts between `present` and `turnover_candidate` states are identified. 
# Focal Tier-1 and Secondary Tier-1 candidates are kept separate and ordered deterministically.
#
# Figure structure:
# phylogeny | species labels | climate strip | gap | CRE-state heatmap
#
# Notes:
# - D. melanogaster is not displayed as a species row.
# - Heatmap columns represent D. melanogaster reference CREs.
#
# Output:
# - one PNG and PDF figure per focal clade
# - focal_clade_heatmap_CREs.tsv
# ============================================================


# ============================================================
# Imports
# ============================================================

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle


# ============================================================
# Repository-relative default paths
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
CLASS_DIR = SCRIPT_DIR.parent
PROJECT_ROOT = CLASS_DIR.parent

DEFAULT_MATRIX_FILE = CLASS_DIR / "results" / "cre_turnover_matrix.tsv"
DEFAULT_TREE_FILE = (
    CLASS_DIR
    / "phylogeny"
    / "results"
    / "301Fly_HOG_UCLDtree_40species.nw"
)
DEFAULT_MANIFEST_FILE = (
    PROJECT_ROOT
    / "01_scrmshaw"
    / "external"
    / "combined_manifest.tsv"
)
DEFAULT_TRAITS_FILE = CLASS_DIR / "phylogeny" / "data" / "species_traits.tsv"
DEFAULT_FIG_DIR = CLASS_DIR / "results" / "figures"
DEFAULT_OUT_TSV = CLASS_DIR / "results" / "focal_clade_heatmap_CREs.tsv"


# ============================================================
# Figure geometry
# ============================================================

# All six figures use identical physical dimensions.
FIG_WIDTH = 13.2
FIG_HEIGHT = 6.6
FIGSIZE = (FIG_WIDTH, FIG_HEIGHT)

# Horizontal structure:
# tree | species labels | climate | gap | heatmap
TREE_WIDTH = 1.60
HEATMAP_SPECIES_LABEL_WIDTH = 0.95
CLIMATE_STRIP_WIDTH = 0.085
CLIMATE_HEATMAP_GAP_WIDTH = 0.0055
HEATMAP_WIDTH = 4.44

GRID_WIDTH_RATIOS = [
    TREE_WIDTH,
    HEATMAP_SPECIES_LABEL_WIDTH,
    CLIMATE_STRIP_WIDTH,
    CLIMATE_HEATMAP_GAP_WIDTH,
    HEATMAP_WIDTH,
]

GRID_WSPACE = 0.006

FIG_LEFT = 0.028
FIG_RIGHT = 0.985
FIG_TOP = 0.855
HEATMAP_BOTTOM_MARGIN = 0.345


# ============================================================
# Font sizes
# ============================================================

TITLE_FONTSIZE = 15
TREE_AXIS_FONTSIZE = 9
TREE_TICK_FONTSIZE = 8

HEATMAP_SPECIES_FONTSIZE = 11
HEATMAP_CRE_FONTSIZE = 9

CRE_AXIS_FONTSIZE = 10
BLOCK_LABEL_FONTSIZE = 10


# ============================================================
# Plot styling
# ============================================================

# Species labels remain right-aligned toward the climate strip.
HEATMAP_SPECIES_LABEL_X = 0.91

# Preserve branch-length ratios while using almost the entire tree axis.
TREE_RIGHT_PADDING_FACTOR = 1.01

HEATMAP_BORDER_WIDTH = 0.55
FOCAL_BLOCK_LINEWIDTH = 1.5
OUTPUT_DPI = 300


# ============================================================
# Focal-clade definitions
# ============================================================

CLADES = {
    "rufa_group": {
        "title": "Rufa group",
        "focal": "druf",
        "species": [
            "druf",
            "dkik",
            "dbun",
            "dbir",
            "d_serrata",
        ],
    },
    "immigrans_group": {
        "title": "Immigrans group",
        "focal": "dimm",
        "species": [
            "dimm",
            "dfor",
            "dsul",
            "dnas",
        ],
    },
    "obscura_group": {
        "title": "Subobscura group",
        "focal": "dsub",
        "species": [
            "dsub",
            "dbif",
            "dobs",
        ],
    },
    "azteca_affinis_miranda_group": {
        "title": "Azteca group",
        "focal": "dazt",
        "species": [
            "dazt",
            "d_affinis",
            "dmir",
            "dper",
        ],
    },
    "teissieri_group": {
        "title": "Teissieri group",
        "focal": "dtei",
        "species": [
            "dtei",
            "d_simulans",
            "d_lutescens",
            "dsuz",
        ],
    },
    "repleta_group": {
        "title": "Repleta group",
        "focal": "d_repleta",
        "species": [
            "d_repleta",
            "d_buzzatii",
            "d_mojavensis",
            "d_arizonae",
        ],
    },
}


# ============================================================
# CRE-state palette
# ============================================================

STATE_COLORS = {
    "present": "#0072B2",
    "turnover_candidate": "#E69F00",
}

STATE_TO_NUM = {
    "present": 0,
    "turnover_candidate": 1,
}

STATE_CMAP = ListedColormap([
    STATE_COLORS["present"],
    STATE_COLORS["turnover_candidate"],
])

STATE_NORM = BoundaryNorm(
    [-0.5, 0.5, 1.5],
    STATE_CMAP.N,
)


# ============================================================
# Climate palette
# ============================================================

CLIMATE_COLORS = {
    "TROP": "#E41A1C",
    "ARID": "#FBC02D",
    "TEMP": "#009E73",
    "BORE": "#0072B2",
}

CLIMATE_TO_NUM = {
    "TROP": 0,
    "ARID": 1,
    "TEMP": 2,
    "BORE": 3,
}

CLIMATE_CMAP = ListedColormap([
    CLIMATE_COLORS["TROP"],
    CLIMATE_COLORS["ARID"],
    CLIMATE_COLORS["TEMP"],
    CLIMATE_COLORS["BORE"],
])

CLIMATE_NORM = BoundaryNorm(
    [-0.5, 0.5, 1.5, 2.5, 3.5],
    CLIMATE_CMAP.N,
)


# ============================================================
# Command-line arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot strict focal-clade CRE-state contrasts using a "
            "phylogeny, climate annotations, and the CRE-state matrix."
        )
    )

    parser.add_argument(
        "--matrix",
        type=Path,
        default=DEFAULT_MATRIX_FILE,
        help=f"CRE-state matrix (default: {DEFAULT_MATRIX_FILE})",
    )
    parser.add_argument(
        "--tree",
        type=Path,
        default=DEFAULT_TREE_FILE,
        help=f"Pruned 40-species Newick tree (default: {DEFAULT_TREE_FILE})",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST_FILE,
        help=f"Combined SCRMshaw manifest (default: {DEFAULT_MANIFEST_FILE})",
    )
    parser.add_argument(
        "--traits",
        type=Path,
        default=DEFAULT_TRAITS_FILE,
        help=f"Species trait table (default: {DEFAULT_TRAITS_FILE})",
    )
    parser.add_argument(
        "--fig-dir",
        type=Path,
        default=DEFAULT_FIG_DIR,
        help=f"Figure output directory (default: {DEFAULT_FIG_DIR})",
    )
    parser.add_argument(
        "--out-tsv",
        type=Path,
        default=DEFAULT_OUT_TSV,
        help=f"Plotted-CRE summary table (default: {DEFAULT_OUT_TSV})",
    )

    return parser.parse_args()


# ============================================================
# General helpers
# ============================================================

def require_file(path):
    if not path.is_file():
        raise SystemExit(
            f"ERROR: required file not found or not a regular file:\n{path}"
        )


def require_columns(df, required, source_name):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(
            f"ERROR: {source_name} is missing required columns:\n"
            + "\n".join(missing)
        )


def tree_name_from_species(species_name):
    """Convert a manifest species name to the tree naming convention."""
    return str(species_name).strip().upper().replace(" ", "_")


def short_species_name(slug, manifest_by_slug):
    """
    Return an abbreviated species name.

    Species using externally generated SCRMshaw predictions receive '*'.
    """
    if slug not in manifest_by_slug:
        return slug

    row = manifest_by_slug[slug]
    name = row["species"]
    source = row["source"]

    if name.startswith("Drosophila "):
        label = "D. " + name.split(" ", 1)[1]
    elif name.startswith("Zaprionus "):
        label = "Z. " + name.split(" ", 1)[1]
    else:
        label = name

    if source == "external":
        label += "*"

    return label


def count_turnover_states(row):
    return int((row == "turnover_candidate").sum())


def contrast_direction(discordant_state):
    if discordant_state == "turnover_candidate":
        return "discordant_turnover"
    if discordant_state == "present":
        return "discordant_present"
    raise ValueError(f"Unexpected discordant state: {discordant_state}")


# ============================================================
# Candidate-selection helpers
# ============================================================

def is_focal_tier1(row, focal, comparison_species):
    """
    Test whether the predefined focal species differs from every comparison
    species in a strict present <-> turnover_candidate contrast.
    """
    focal_state = row[focal]
    comparison_states = [row[species] for species in comparison_species]

    focal_turnover = (
        focal_state == "turnover_candidate"
        and all(state == "present" for state in comparison_states)
    )

    focal_present = (
        focal_state == "present"
        and all(
            state == "turnover_candidate"
            for state in comparison_states
        )
    )

    return focal_turnover or focal_present


def get_singleton_tier1_info(row):
    """
    Identify a strict singleton present <-> turnover_candidate contrast.

    Returns
    -------
    tuple or None
        (discordant_species, discordant_state, consensus_state)
    """
    allowed_states = {"present", "turnover_candidate"}

    if not all(state in allowed_states for state in row.values):
        return None

    counts = row.value_counts()
    if len(counts) != 2:
        return None

    singleton_states = counts[counts == 1].index.tolist()
    consensus_states = counts[counts == len(row) - 1].index.tolist()

    if len(singleton_states) != 1 or len(consensus_states) != 1:
        return None

    discordant_state = singleton_states[0]
    consensus_state = consensus_states[0]
    discordant_species = row.index[row == discordant_state][0]

    return discordant_species, discordant_state, consensus_state


# ============================================================
# Phylogeny helpers
# ============================================================

def prune_tree_to_species(tree_file, species, slug_to_tree):
    tree = Phylo.read(tree_file, "newick")

    wanted_tree_names = {
        slug_to_tree[species_name]
        for species_name in species
    }

    available_tree_names = {
        tip.name
        for tip in tree.get_terminals()
    }

    missing = sorted(wanted_tree_names - available_tree_names)
    if missing:
        raise SystemExit(
            "ERROR: focal-clade species missing from phylogenetic tree:\n"
            + "\n".join(missing)
        )

    for tip in list(tree.get_terminals()):
        if tip.name not in wanted_tree_names:
            tree.prune(tip)

    return tree


# ============================================================
# Input loading and validation
# ============================================================

def load_manifest(path):
    manifest = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    require_columns(
        manifest,
        {"slug", "species", "source"},
        "combined_manifest.tsv",
    )

    for column in ("slug", "species", "source"):
        manifest[column] = manifest[column].astype(str).str.strip()

    if manifest["slug"].eq("").any():
        raise SystemExit("ERROR: empty species slug in combined_manifest.tsv")

    if manifest["species"].eq("").any():
        raise SystemExit("ERROR: empty species name in combined_manifest.tsv")

    if manifest["slug"].duplicated().any():
        duplicates = sorted(
            manifest.loc[
                manifest["slug"].duplicated(keep=False),
                "slug",
            ].unique()
        )
        raise SystemExit(
            "ERROR: duplicate species slugs in combined_manifest.tsv:\n"
            + "\n".join(duplicates)
        )

    manifest["tree_name"] = manifest["species"].map(tree_name_from_species)

    if manifest["tree_name"].duplicated().any():
        duplicates = sorted(
            manifest.loc[
                manifest["tree_name"].duplicated(keep=False),
                "tree_name",
            ].unique()
        )
        raise SystemExit(
            "ERROR: duplicate inferred tree names in combined_manifest.tsv:\n"
            + "\n".join(duplicates)
        )

    return manifest


def load_traits(path):
    traits = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    require_columns(
        traits,
        {"tree_name", "climatic_zone"},
        "species_traits.tsv",
    )

    for column in ("tree_name", "climatic_zone"):
        traits[column] = traits[column].astype(str).str.strip()

    if traits["tree_name"].duplicated().any():
        duplicates = sorted(
            traits.loc[
                traits["tree_name"].duplicated(keep=False),
                "tree_name",
            ].unique()
        )
        raise SystemExit(
            "ERROR: duplicate tree_name values in species_traits.tsv:\n"
            + "\n".join(duplicates)
        )

    return traits


def load_matrix(path):
    matrix = pd.read_csv(path, sep="\t", index_col=0, dtype=str).fillna("")

    if matrix.index.duplicated().any():
        duplicates = sorted(
            matrix.index[matrix.index.duplicated()].unique().tolist()
        )
        raise SystemExit(
            "ERROR: duplicated CRE IDs in cre_turnover_matrix.tsv:\n"
            + "\n".join(duplicates)
        )

    return matrix


def validate_clades(matrix, slug_to_tree):
    all_clade_species = sorted({
        species
        for config in CLADES.values()
        for species in config["species"]
    })

    missing_from_matrix = [
        species
        for species in all_clade_species
        if species not in matrix.columns
    ]
    if missing_from_matrix:
        raise SystemExit(
            "ERROR: focal-clade species missing from cre_turnover_matrix.tsv:\n"
            + "\n".join(missing_from_matrix)
        )

    missing_from_manifest = [
        species
        for species in all_clade_species
        if species not in slug_to_tree
    ]
    if missing_from_manifest:
        raise SystemExit(
            "ERROR: focal-clade species missing from combined_manifest.tsv:\n"
            + "\n".join(missing_from_manifest)
        )

    for clade_name, config in CLADES.items():
        declared_species = config["species"]
        focal = config["focal"]

        if focal not in declared_species:
            raise SystemExit(
                f"ERROR: focal species {focal!r} is not listed in "
                f"clade {clade_name!r}."
            )

        if len(declared_species) != len(set(declared_species)):
            raise SystemExit(
                f"ERROR: duplicate species in clade definition {clade_name!r}."
            )


# ============================================================
# Candidate selection
# ============================================================

def build_clade_data(matrix, tree_file, slug_to_tree, tree_to_slug):
    clade_data = {}
    summary_rows = []

    for clade_name, config in CLADES.items():
        focal = config["focal"]
        declared_species = config["species"]

        clade_tree = prune_tree_to_species(
            tree_file,
            declared_species,
            slug_to_tree,
        )

        tree_tip_names = [
            tip.name
            for tip in clade_tree.get_terminals()
        ]

        missing_reverse_mapping = [
            name
            for name in tree_tip_names
            if name not in tree_to_slug
        ]
        if missing_reverse_mapping:
            raise SystemExit(
                "ERROR: tree tips could not be mapped back to species slugs:\n"
                + "\n".join(sorted(missing_reverse_mapping))
            )

        species_order = [
            tree_to_slug[tree_name]
            for tree_name in tree_tip_names
        ]

        comparison_species = [
            species
            for species in declared_species
            if species != focal
        ]

        sub = matrix[declared_species].copy()

        # Focal Tier 1.
        focal_mask = sub.apply(
            lambda row: is_focal_tier1(
                row,
                focal,
                comparison_species,
            ),
            axis=1,
        )

        focal_ids = sub.index[focal_mask].tolist()
        focal_ids = sorted(
            focal_ids,
            key=lambda cre_id: (
                count_turnover_states(sub.loc[cre_id, declared_species]),
                cre_id,
            ),
        )

        # All strict singleton Tier-1 contrasts.
        singleton_info = {}
        for cre_id, row in sub.iterrows():
            info = get_singleton_tier1_info(row)
            if info is not None:
                singleton_info[cre_id] = info

        # Secondary Tier 1 excludes focal Tier-1 candidates.
        secondary_ids = [
            cre_id
            for cre_id in singleton_info
            if cre_id not in focal_ids
        ]

        species_rank = {
            species: index
            for index, species in enumerate(species_order)
        }

        def secondary_sort_key(cre_id):
            discordant_species, _, _ = singleton_info[cre_id]
            return (
                count_turnover_states(sub.loc[cre_id, declared_species]),
                species_rank[discordant_species],
                cre_id,
            )

        secondary_ids = sorted(secondary_ids, key=secondary_sort_key)
        cre_order = focal_ids + secondary_ids

        if not cre_order:
            raise SystemExit(
                f"ERROR: no Tier-1 CREs found for {clade_name}."
            )

        selected = sub.loc[cre_order].copy()

        clade_data[clade_name] = {
            "matrix": selected,
            "focal_ids": focal_ids,
            "secondary_ids": secondary_ids,
            "singleton_info": singleton_info,
            "species_order": species_order,
        }

        for plot_order, cre_id in enumerate(cre_order, start=1):
            discordant_species, discordant_state, consensus_state = (
                singleton_info[cre_id]
            )
            n_turnover = count_turnover_states(
                sub.loc[cre_id, declared_species]
            )

            summary_rows.append({
                "clade": clade_name,
                "clade_title": config["title"],
                "plot_order": plot_order,
                "cre_id": cre_id,
                "category": (
                    "focal_tier1"
                    if cre_id in focal_ids
                    else "secondary_tier1"
                ),
                "focal_species": focal,
                "discordant_species": discordant_species,
                "discordant_state": discordant_state,
                "consensus_state": consensus_state,
                "contrast_direction": contrast_direction(discordant_state),
                "n_turnover_states": n_turnover,
                "n_present_states": len(declared_species) - n_turnover,
            })

    return clade_data, pd.DataFrame(summary_rows)


# ============================================================
# Plot helpers
# ============================================================

def climate_array(species_order, slug_to_tree, tree_to_climate):
    values = []

    for species in species_order:
        tree_name = slug_to_tree[species]
        climate = tree_to_climate.get(tree_name, "")

        if not climate:
            raise SystemExit(
                f"ERROR: missing climate annotation for {species}."
            )

        if climate not in CLIMATE_TO_NUM:
            raise SystemExit(
                f"ERROR: unknown climatic zone {climate!r} for {species}."
            )

        values.append(CLIMATE_TO_NUM[climate])

    return np.asarray(values, dtype=float).reshape(-1, 1)


def style_tree_axis(ax_tree, tree, shared_ylim):
    ax_tree.set_ylim(*shared_ylim)

    depths = tree.depths()
    terminal_depths = [depths[tip] for tip in tree.get_terminals()]
    max_tree_depth = max(terminal_depths)

    ax_tree.set_xlim(
        0,
        max_tree_depth * TREE_RIGHT_PADDING_FACTOR,
    )

    ax_tree.set_ylabel("")
    ax_tree.set_xlabel(
        "Branch length",
        fontsize=TREE_AXIS_FONTSIZE,
        labelpad=4,
    )
    ax_tree.tick_params(
        axis="x",
        labelsize=TREE_TICK_FONTSIZE,
    )
    ax_tree.tick_params(
        axis="y",
        left=False,
        labelleft=False,
    )

    for spine in ("top", "right", "left"):
        ax_tree.spines[spine].set_visible(False)


def draw_species_labels(
    ax_labels,
    species_order,
    species_labels,
    focal,
    shared_ylim,
):
    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(*shared_ylim)

    for row_index, (species, label) in enumerate(
        zip(species_order, species_labels),
        start=1,
    ):
        ax_labels.text(
            HEATMAP_SPECIES_LABEL_X,
            row_index,
            label,
            ha="right",
            va="center",
            fontsize=HEATMAP_SPECIES_FONTSIZE,
            fontweight="bold" if species == focal else "normal",
        )

    ax_labels.set_xticks([])
    ax_labels.set_yticks([])

    for spine in ax_labels.spines.values():
        spine.set_visible(False)


def draw_climate_strip(ax_climate, climate_numeric, y_edges, shared_ylim):
    climate_x_edges = np.array([-0.5, 0.5])

    ax_climate.pcolormesh(
        climate_x_edges,
        y_edges,
        climate_numeric,
        cmap=CLIMATE_CMAP,
        norm=CLIMATE_NORM,
        shading="flat",
        edgecolors="none",
        antialiased=False,
    )

    ax_climate.set_xlim(-0.5, 0.5)
    ax_climate.set_ylim(*shared_ylim)
    ax_climate.set_xticks([])
    ax_climate.set_yticks([])

    for spine in ax_climate.spines.values():
        spine.set_visible(False)


def draw_heatmap(
    ax_heat,
    numeric,
    n_species,
    n_focal,
    n_secondary,
    shared_ylim,
):
    n_cre = numeric.shape[1]
    y_edges = np.arange(n_species + 1) + 0.5
    cre_x_edges = np.arange(n_cre + 1) - 0.5

    ax_heat.pcolormesh(
        cre_x_edges,
        y_edges,
        numeric.values,
        cmap=STATE_CMAP,
        norm=STATE_NORM,
        shading="flat",
        edgecolors="white",
        linewidth=HEATMAP_BORDER_WIDTH,
        antialiased=False,
    )

    ax_heat.set_xlim(-0.5, n_cre - 0.5)
    ax_heat.set_ylim(*shared_ylim)
    ax_heat.set_yticks([])
    ax_heat.set_ylabel("")

    ax_heat.set_xticks(np.arange(n_cre))
    ax_heat.set_xticklabels(
        numeric.columns,
        rotation=90,
        fontsize=HEATMAP_CRE_FONTSIZE,
        ha="center",
        va="top",
    )
    ax_heat.tick_params(
        axis="x",
        length=3,
        pad=3,
    )

    ax_heat.set_xlabel(
        r"$D.\ melanogaster$ reference CRE",
        fontsize=CRE_AXIS_FONTSIZE,
        labelpad=10,
    )

    if n_focal > 0:
        focal_rectangle = Rectangle(
            (-0.5, 0.5),
            n_focal,
            n_species,
            fill=False,
            edgecolor="#222222",
            linewidth=FOCAL_BLOCK_LINEWIDTH,
            zorder=10,
        )
        ax_heat.add_patch(focal_rectangle)

    if n_focal > 0 and n_secondary > 0:
        ax_heat.axvline(
            n_focal - 0.5,
            color="#222222",
            linewidth=FOCAL_BLOCK_LINEWIDTH,
            zorder=11,
        )

    if n_focal > 0:
        focal_center = (n_focal - 1) / 2
        ax_heat.text(
            focal_center,
            1.045,
            "Focal Tier 1",
            transform=ax_heat.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=BLOCK_LABEL_FONTSIZE,
            fontweight="bold",
        )

    if n_secondary > 0:
        secondary_center = (n_focal + n_cre - 1) / 2
        ax_heat.text(
            secondary_center,
            1.045,
            "Secondary Tier 1",
            transform=ax_heat.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=BLOCK_LABEL_FONTSIZE,
        )


def plot_clade(
    clade_name,
    config,
    clade_info,
    tree_file,
    slug_to_tree,
    manifest_by_slug,
    tree_to_climate,
    fig_dir,
):
    focal = config["focal"]
    declared_species = config["species"]

    selected = clade_info["matrix"]
    focal_ids = clade_info["focal_ids"]
    secondary_ids = clade_info["secondary_ids"]
    species_order = clade_info["species_order"]

    n_focal = len(focal_ids)
    n_secondary = len(secondary_ids)

    heatmap = selected[species_order].T

    invalid_states = sorted(
        set(pd.unique(heatmap.values.ravel()))
        - set(STATE_TO_NUM)
    )
    if invalid_states:
        raise SystemExit(
            f"ERROR: invalid states in Tier-1 heatmap for {clade_name}:\n"
            + "\n".join(str(x) for x in invalid_states)
        )

    numeric = heatmap.apply(
        lambda column: column.map(STATE_TO_NUM)
    ).astype(float)
    n_species, n_cre = numeric.shape

    shared_ylim = (n_species + 0.5, 0.5)
    y_edges = np.arange(n_species + 1) + 0.5

    species_labels = [
        short_species_name(species, manifest_by_slug)
        for species in species_order
    ]

    climate_numeric = climate_array(
        species_order,
        slug_to_tree,
        tree_to_climate,
    )

    fig = plt.figure(figsize=FIGSIZE)
    grid = fig.add_gridspec(
        nrows=1,
        ncols=5,
        width_ratios=GRID_WIDTH_RATIOS,
        wspace=GRID_WSPACE,
    )

    ax_tree = fig.add_subplot(grid[0, 0])
    ax_labels = fig.add_subplot(grid[0, 1])
    ax_climate = fig.add_subplot(grid[0, 2])
    ax_spacer = fig.add_subplot(grid[0, 3])
    ax_heat = fig.add_subplot(grid[0, 4])

    ax_spacer.axis("off")

    tree = prune_tree_to_species(
        tree_file,
        declared_species,
        slug_to_tree,
    )

    # Tip labels are drawn in the dedicated species-label axis.
    for tip in tree.get_terminals():
        tip.name = ""

    Phylo.draw(
        tree,
        axes=ax_tree,
        do_show=False,
        show_confidence=False,
    )

    style_tree_axis(ax_tree, tree, shared_ylim)
    draw_species_labels(
        ax_labels,
        species_order,
        species_labels,
        focal,
        shared_ylim,
    )
    draw_climate_strip(
        ax_climate,
        climate_numeric,
        y_edges,
        shared_ylim,
    )
    draw_heatmap(
        ax_heat,
        numeric,
        n_species,
        n_focal,
        n_secondary,
        shared_ylim,
    )

    fig.suptitle(
        config["title"],
        fontsize=TITLE_FONTSIZE,
        fontweight="bold",
        y=0.975,
    )

    fig.subplots_adjust(
        left=FIG_LEFT,
        right=FIG_RIGHT,
        top=FIG_TOP,
        bottom=HEATMAP_BOTTOM_MARGIN,
    )

    stem = f"focal_clade_tier1_heatmap_{clade_name}"
    out_png = fig_dir / f"{stem}.png"
    out_pdf = fig_dir / f"{stem}.pdf"

    # Deliberately no bbox_inches="tight": all figures retain identical size.
    fig.savefig(
        out_png,
        dpi=OUTPUT_DPI,
        facecolor="white",
    )
    fig.savefig(
        out_pdf,
        facecolor="white",
    )
    plt.close(fig)

    return {
        "png": out_png,
        "pdf": out_pdf,
        "n_cre": n_cre,
        "n_focal": n_focal,
        "n_secondary": n_secondary,
        "species_order": species_order,
        "selected": selected,
    }


# ============================================================
# Reporting
# ============================================================

def report_clade(config, result, declared_species):
    print()
    print("=" * 72)
    print(config["title"])
    print("=" * 72)
    print("Tree order: " + ", ".join(result["species_order"]))
    print(
        f"Displayed CREs: {result['n_cre']} "
        f"({result['n_focal']} focal Tier 1 + "
        f"{result['n_secondary']} secondary Tier 1)"
    )
    print("CRE display order:")

    selected = result["selected"]
    focal_ids = set(selected.index[: result["n_focal"]])

    for display_index, cre_id in enumerate(selected.index, start=1):
        n_turnover = count_turnover_states(
            selected.loc[cre_id, declared_species]
        )
        category = "focal" if cre_id in focal_ids else "secondary"

        print(
            f"  {display_index:>2}. "
            f"{cre_id:<18} "
            f"{category:<10} "
            f"turnover states = {n_turnover}"
        )

    print()
    print(f"Wrote PNG:\n{result['png']}")
    print(f"Wrote PDF:\n{result['pdf']}")


def report_final(written_pngs, written_pdfs, summary_df, out_tsv):
    print()
    print("=" * 72)
    print("FOCAL-CLADE HEATMAPS COMPLETE")
    print("=" * 72)
    print(f"Focal clades: {len(CLADES)}")
    print(f"PNG figures written: {len(written_pngs)}")
    print(f"PDF figures written: {len(written_pdfs)}")
    print()
    print(f"Fixed figure dimensions: {FIG_WIDTH} x {FIG_HEIGHT} inches")
    print(f"Species-label fontsize: {HEATMAP_SPECIES_FONTSIZE}")
    print(f"CRE-label fontsize: {HEATMAP_CRE_FONTSIZE}")
    print(f"Species-label width ratio: {HEATMAP_SPECIES_LABEL_WIDTH}")
    print(f"Species-label x position: {HEATMAP_SPECIES_LABEL_X}")
    print(f"Bottom margin: {HEATMAP_BOTTOM_MARGIN}")
    print("Heatmap renderer: pcolormesh")
    print("PDF-safe categorical rendering: YES")
    print("Continuous climate strip: YES")
    print("Dedicated climate-to-heatmap spacer: YES")
    print(
        "Climate-to-heatmap gap width ratio: "
        f"{CLIMATE_HEATMAP_GAP_WIDTH}"
    )
    print("Species labels right-aligned before climate strip: YES")
    print("Tree branch lengths modified: NO")
    print("Panel letters A-F: NO")
    print("Legends in individual figures: NO")
    print()
    print("CRE ordering:")
    print("  Focal Tier 1 -> Secondary Tier 1")
    print(
        "  within blocks: increasing number of "
        "turnover_candidate states"
    )
    print(
        "  Secondary tie-break: discordant species in displayed "
        "tree order -> CRE ID"
    )
    print()
    print(f"Wrote plotted-CRE table:\n{out_tsv}")
    print(
        "Total clade-specific CRE occurrences: "
        f"{len(summary_df)}"
    )


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    matrix_file = args.matrix.resolve()
    tree_file = args.tree.resolve()
    manifest_file = args.manifest.resolve()
    traits_file = args.traits.resolve()
    fig_dir = args.fig_dir.resolve()
    out_tsv = args.out_tsv.resolve()

    for path in (
        matrix_file,
        tree_file,
        manifest_file,
        traits_file,
    ):
        require_file(path)

    fig_dir.mkdir(parents=True, exist_ok=True)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    manifest = load_manifest(manifest_file)
    traits = load_traits(traits_file)
    matrix = load_matrix(matrix_file)

    manifest_by_slug = (
        manifest
        .set_index("slug")[["species", "source"]]
        .to_dict(orient="index")
    )

    slug_to_tree = dict(zip(manifest["slug"], manifest["tree_name"]))
    tree_to_slug = dict(zip(manifest["tree_name"], manifest["slug"]))
    tree_to_climate = dict(zip(traits["tree_name"], traits["climatic_zone"]))

    validate_clades(matrix, slug_to_tree)

    clade_data, summary_df = build_clade_data(
        matrix,
        tree_file,
        slug_to_tree,
        tree_to_slug,
    )

    summary_df.to_csv(
        out_tsv,
        sep="\t",
        index=False,
    )

    written_pngs = []
    written_pdfs = []

    for clade_name, config in CLADES.items():
        result = plot_clade(
            clade_name,
            config,
            clade_data[clade_name],
            tree_file,
            slug_to_tree,
            manifest_by_slug,
            tree_to_climate,
            fig_dir,
        )

        written_pngs.append(result["png"])
        written_pdfs.append(result["pdf"])

        report_clade(
            config,
            result,
            config["species"],
        )

    report_final(
        written_pngs,
        written_pdfs,
        summary_df,
        out_tsv,
    )


if __name__ == "__main__":
    main()
