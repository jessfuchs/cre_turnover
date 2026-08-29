#!/usr/bin/env python3

"""
Plot stringent lineage-specific CRE-state contrasts for six focal clades.

One separate figure is generated for each clade.

Figure structure
----------------
Each figure contains:

    phylogeny | species labels | climate strip | gap | CRE-state heatmap

All six figures use identical physical dimensions and identical subplot
proportions so that they can later be combined consistently in LaTeX.

Candidate categories
--------------------
1. Focal Tier 1

   The predefined focal species differs from all comparison species,
   restricted to the contrast:

       present <-> turnover_candidate

   Both directions are retained:

       focal = turnover_candidate
       comparisons = present

   OR

       focal = present
       comparisons = turnover_candidate


2. Secondary Tier 1

   Exactly one species within the clade differs from all remaining
   species, again restricted to:

       present <-> turnover_candidate

   Focal Tier-1 CREs are excluded from the secondary set to avoid
   duplicate display.


CRE ordering
------------
Focal and Secondary Tier-1 CREs remain separate blocks.

Within each block, CREs are ordered by increasing number of
turnover_candidate states:

    predominantly present
        ->
    predominantly turnover_candidate

For singleton contrasts this corresponds to:

    one turnover_candidate + all other species present
        ->
    one present + all other species turnover_candidate

Within the Secondary Tier-1 block, candidates with the same number of
turnover states are additionally ordered by the discordant species
according to the displayed phylogenetic order and then by CRE ID.


Display
-------
- No A-F panel letters are added in Python.
- The clade name is retained as the figure title.
- The focal species is shown in bold.
- Species using externally generated SCRMshaw prediction sets receive '*'.
- Climatic zone is shown as one continuous vertical strip.
- A dedicated empty spacer separates the climate strip from the heatmap.
- PDF-safe categorical rendering uses pcolormesh().
- PNG and PDF versions are written.
- No bbox_inches="tight" is used, ensuring identical output dimensions.

D. melanogaster is not displayed as a species row. CRE identifiers
refer to the D. melanogaster reference CRE set.
"""


# ============================================================
# Imports
# ============================================================

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle

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


MATRIX_FILE = (
    RESULTS_DIR
    / "cre_turnover_matrix.tsv"
)

TREE_FILE = (
    PHYLO_DIR
    / "results"
    / "301Fly_HOG_UCLDtree_40species.nw"
)

MANIFEST_FILE = (
    PROJECT_DIR
    / "external_scrmshaw"
    / "combined_manifest.tsv"
)

TRAITS_FILE = (
    PHYLO_DIR
    / "data"
    / "species_traits.tsv"
)

OUT_TSV = (
    RESULTS_DIR
    / "focal_clade_heatmap_CREs.tsv"
)


# ============================================================
# Figure geometry
# ============================================================

# All six figures have exactly the same outer dimensions.
FIG_WIDTH = 13.2
FIG_HEIGHT = 6.6

FIGSIZE = (
    FIG_WIDTH,
    FIG_HEIGHT,
)


# ------------------------------------------------------------
# Horizontal structure:
#
# tree | species | climate | gap | heatmap
#
# These are the main values to adjust if the spacing should
# later be fine-tuned.
# ------------------------------------------------------------

TREE_WIDTH = 1.60

# Smaller value brings species names + climate strip closer
# to the tree while retaining their mutual spacing.
SPECIES_LABEL_WIDTH = 0.75

CLIMATE_STRIP_WIDTH = 0.085

# Controls ONLY the additional space between climate strip
# and CRE-state heatmap.
CLIMATE_HEATMAP_GAP_WIDTH = 0.0055

HEATMAP_WIDTH = 4.44


GRID_WIDTH_RATIOS = [
    TREE_WIDTH,
    SPECIES_LABEL_WIDTH,
    CLIMATE_STRIP_WIDTH,
    CLIMATE_HEATMAP_GAP_WIDTH,
    HEATMAP_WIDTH,
]


# Small general spacing between all GridSpec columns.
GRID_WSPACE = 0.006


# Fixed outer margins.
FIG_LEFT = 0.028
FIG_RIGHT = 0.985
FIG_TOP = 0.855
FIG_BOTTOM = 0.305


# ============================================================
# Font sizes
# ============================================================

TITLE_FONTSIZE = 15

TREE_AXIS_FONTSIZE = 9
TREE_TICK_FONTSIZE = 8

SPECIES_FONTSIZE = 10.5

CRE_LABEL_FONTSIZE = 7.5
CRE_AXIS_FONTSIZE = 10

BLOCK_LABEL_FONTSIZE = 10


# ============================================================
# Plot styling
# ============================================================

# Species labels remain right-aligned toward the climate strip.
# This keeps the species-name <-> climate-strip spacing uniform.
SPECIES_LABEL_X = 0.89

# Tree uses nearly all of its available x-axis width.
TREE_RIGHT_PADDING_FACTOR = 1.01

HEATMAP_BORDER_WIDTH = 0.55

FOCAL_BLOCK_LINEWIDTH = 1.5


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
    [
        -0.5,
        0.5,
        1.5,
    ],
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
    [
        -0.5,
        0.5,
        1.5,
        2.5,
        3.5,
    ],
    CLIMATE_CMAP.N,
)


# ============================================================
# Helper functions
# ============================================================

def require_file(path):
    """
    Stop execution if a required input file does not exist.
    """

    if not path.exists():
        raise SystemExit(
            f"ERROR: required file not found:\n{path}"
        )


def tree_name_from_species(species_name):
    """
    Convert manifest species name into the naming convention
    used by the phylogenetic tree.
    """

    return (
        str(species_name)
        .strip()
        .upper()
        .replace(
            " ",
            "_",
        )
    )


def short_species_name(
    slug,
    manifest,
):
    """
    Return abbreviated species name.

    Species using externally generated SCRMshaw predictions
    receive an asterisk.
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
            + name.split(
                " ",
                1,
            )[1]
        )


    elif name.startswith("Zaprionus "):

        label = (
            "Z. "
            + name.split(
                " ",
                1,
            )[1]
        )


    else:

        label = name


    if source == "external":
        label += "*"


    return label


def is_focal_tier1(
    row,
    focal,
    comparison_species,
):
    """
    Test whether the focal species differs from all comparison
    species in a strict present <-> turnover_candidate contrast.
    """

    focal_state = row[focal]

    comparison_states = [
        row[species]
        for species
        in comparison_species
    ]


    # Focal turnover candidate, all others positional matches.
    if (
        focal_state == "turnover_candidate"

        and all(
            state == "present"
            for state
            in comparison_states
        )
    ):
        return True


    # Focal positional match, all others turnover candidates.
    if (
        focal_state == "present"

        and all(
            state == "turnover_candidate"
            for state
            in comparison_states
        )
    ):
        return True


    return False


def get_singleton_tier1_info(row):
    """
    Identify a strict singleton present <-> turnover_candidate contrast.

    Returns
    -------
    tuple or None

        (
            discordant_species,
            discordant_state,
            consensus_state,
        )
    """

    allowed_states = {
        "present",
        "turnover_candidate",
    }


    if not all(
        state in allowed_states
        for state
        in row.values
    ):
        return None


    counts = row.value_counts()


    if len(counts) != 2:
        return None


    singleton_states = (
        counts[
            counts == 1
        ]
        .index
        .tolist()
    )


    consensus_states = (
        counts[
            counts == len(row) - 1
        ]
        .index
        .tolist()
    )


    if len(singleton_states) != 1:
        return None


    if len(consensus_states) != 1:
        return None


    discordant_state = (
        singleton_states[0]
    )

    consensus_state = (
        consensus_states[0]
    )


    discordant_species = (
        row.index[
            row == discordant_state
        ][0]
    )


    return (
        discordant_species,
        discordant_state,
        consensus_state,
    )


def count_turnover_states(row):
    """
    Count turnover_candidate states across one clade.
    """

    return int(
        (
            row
            == "turnover_candidate"
        ).sum()
    )


def contrast_direction(
    discordant_state,
):
    """
    Return compact direction label for the output table.
    """

    if (
        discordant_state
        == "turnover_candidate"
    ):
        return "discordant_turnover"


    if (
        discordant_state
        == "present"
    ):
        return "discordant_present"


    raise ValueError(
        f"Unexpected discordant state: "
        f"{discordant_state}"
    )


def prune_tree_to_species(
    tree_file,
    species,
    slug_to_tree,
):
    """
    Read the full tree and prune it to one focal clade.
    """

    tree = Phylo.read(
        tree_file,
        "newick",
    )


    wanted_tree_names = {
        slug_to_tree[species_name]
        for species_name
        in species
    }


    available_tree_names = {
        tip.name
        for tip
        in tree.get_terminals()
    }


    missing = sorted(
        wanted_tree_names
        - available_tree_names
    )


    if missing:

        raise SystemExit(
            "ERROR: focal-clade species missing "
            "from phylogenetic tree:\n"
            + "\n".join(
                missing
            )
        )


    for tip in list(
        tree.get_terminals()
    ):

        if (
            tip.name
            not in wanted_tree_names
        ):
            tree.prune(
                tip
            )


    return tree


# ============================================================
# Input checks
# ============================================================

for path in [
    MATRIX_FILE,
    TREE_FILE,
    MANIFEST_FILE,
    TRAITS_FILE,
]:
    require_file(
        path
    )


# ============================================================
# Load manifest
# ============================================================

manifest = pd.read_csv(
    MANIFEST_FILE,
    sep="\t",
    dtype=str,
).fillna("")


required_manifest_columns = {
    "slug",
    "species",
    "source",
}


missing_manifest_columns = (
    required_manifest_columns
    - set(
        manifest.columns
    )
)


if missing_manifest_columns:

    raise SystemExit(
        "ERROR: combined_manifest.tsv "
        "missing required columns:\n"
        + "\n".join(
            sorted(
                missing_manifest_columns
            )
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


manifest["tree_name"] = (
    manifest["species"]
    .map(
        tree_name_from_species
    )
)


# ============================================================
# Validate manifest uniqueness
# ============================================================

if manifest["slug"].duplicated().any():

    duplicates = sorted(
        manifest.loc[
            manifest["slug"].duplicated(
                keep=False
            ),
            "slug",
        ]
        .unique()
        .tolist()
    )


    raise SystemExit(
        "ERROR: duplicate species slugs "
        "in combined_manifest.tsv:\n"
        + "\n".join(
            duplicates
        )
    )


if manifest["tree_name"].duplicated().any():

    duplicates = sorted(
        manifest.loc[
            manifest["tree_name"].duplicated(
                keep=False
            ),
            "tree_name",
        ]
        .unique()
        .tolist()
    )


    raise SystemExit(
        "ERROR: duplicate inferred tree names "
        "in combined_manifest.tsv:\n"
        + "\n".join(
            duplicates
        )
    )


slug_to_tree = dict(
    zip(
        manifest["slug"],
        manifest["tree_name"],
    )
)


tree_to_slug = dict(
    zip(
        manifest["tree_name"],
        manifest["slug"],
    )
)


# ============================================================
# Load climate annotations
# ============================================================

traits = pd.read_csv(
    TRAITS_FILE,
    sep="\t",
    dtype=str,
).fillna("")


required_trait_columns = {
    "tree_name",
    "climatic_zone",
}


missing_trait_columns = (
    required_trait_columns
    - set(
        traits.columns
    )
)


if missing_trait_columns:

    raise SystemExit(
        "ERROR: species_traits.tsv "
        "missing required columns:\n"
        + "\n".join(
            sorted(
                missing_trait_columns
            )
        )
    )


for column in [
    "tree_name",
    "climatic_zone",
]:

    traits[column] = (
        traits[column]
        .astype(str)
        .str.strip()
    )


tree_to_climate = dict(
    zip(
        traits["tree_name"],
        traits["climatic_zone"],
    )
)


# ============================================================
# Load CRE-state matrix
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
        + "\n".join(
            duplicates
        )
    )


# ============================================================
# Validate focal-clade species
# ============================================================

all_clade_species = sorted({
    species
    for config
    in CLADES.values()
    for species
    in config["species"]
})


missing_from_matrix = [
    species
    for species
    in all_clade_species
    if species
    not in matrix.columns
]


if missing_from_matrix:

    raise SystemExit(
        "ERROR: focal-clade species missing "
        "from cre_turnover_matrix.tsv:\n"
        + "\n".join(
            missing_from_matrix
        )
    )


missing_from_manifest = [
    species
    for species
    in all_clade_species
    if species
    not in slug_to_tree
]


if missing_from_manifest:

    raise SystemExit(
        "ERROR: focal-clade species missing "
        "from combined_manifest.tsv:\n"
        + "\n".join(
            missing_from_manifest
        )
    )


# ============================================================
# Candidate selection and sorting
# ============================================================

clade_data = {}
summary_rows = []


for clade_name, config in CLADES.items():

    focal = config["focal"]
    declared_species = config["species"]


    # --------------------------------------------------------
    # Determine phylogenetic display order first.
    #
    # This is also used as the secondary sorting order for
    # discordant species.
    # --------------------------------------------------------

    clade_tree = prune_tree_to_species(
        TREE_FILE,
        declared_species,
        slug_to_tree,
    )


    tree_tip_names = [
        tip.name
        for tip
        in clade_tree.get_terminals()
    ]


    species_order = [
        tree_to_slug[tree_name]
        for tree_name
        in tree_tip_names
    ]


    comparison_species = [
        species
        for species
        in declared_species
        if species != focal
    ]


    sub = (
        matrix[
            declared_species
        ]
        .copy()
    )


    # ========================================================
    # Focal Tier 1
    # ========================================================

    focal_mask = sub.apply(
        lambda row: is_focal_tier1(
            row,
            focal,
            comparison_species,
        ),
        axis=1,
    )


    focal_ids = (
        sub.index[
            focal_mask
        ]
        .tolist()
    )


    # Predominantly positional -> predominantly turnover.
    focal_ids = sorted(
        focal_ids,
        key=lambda cre_id: (
            count_turnover_states(
                sub.loc[
                    cre_id,
                    declared_species,
                ]
            ),
            cre_id,
        ),
    )


    # ========================================================
    # All singleton Tier-1 contrasts
    # ========================================================

    singleton_info = {}


    for cre_id, row in sub.iterrows():

        info = (
            get_singleton_tier1_info(
                row
            )
        )


        if info is not None:

            singleton_info[
                cre_id
            ] = info


    # ========================================================
    # Secondary Tier 1
    # ========================================================

    secondary_ids = [
        cre_id
        for cre_id
        in singleton_info
        if cre_id
        not in focal_ids
    ]


    def secondary_sort_key(
        cre_id,
    ):

        (
            discordant_species,
            discordant_state,
            consensus_state,
        ) = singleton_info[
            cre_id
        ]


        n_turnover = (
            count_turnover_states(
                sub.loc[
                    cre_id,
                    declared_species,
                ]
            )
        )


        return (
            # Positional -> turnover.
            n_turnover,

            # Follow displayed phylogenetic species order.
            species_order.index(
                discordant_species
            ),

            # Stable deterministic tie-break.
            cre_id,
        )


    secondary_ids = sorted(
        secondary_ids,
        key=secondary_sort_key,
    )


    # ========================================================
    # Final CRE order
    # ========================================================

    cre_order = (
        focal_ids
        + secondary_ids
    )


    if not cre_order:

        raise SystemExit(
            f"ERROR: no Tier-1 CREs found "
            f"for {clade_name}."
        )


    selected = (
        sub.loc[
            cre_order
        ]
        .copy()
    )


    clade_data[
        clade_name
    ] = {
        "matrix": selected,
        "focal_ids": focal_ids,
        "secondary_ids": secondary_ids,
        "singleton_info": singleton_info,
        "species_order": species_order,
    }


    # ========================================================
    # Summary rows
    # ========================================================

    for plot_order, cre_id in enumerate(
        cre_order,
        start=1,
    ):


        category = (
            "focal_tier1"
            if cre_id in focal_ids
            else "secondary_tier1"
        )


        (
            discordant_species,
            discordant_state,
            consensus_state,
        ) = singleton_info[
            cre_id
        ]


        n_turnover = (
            count_turnover_states(
                sub.loc[
                    cre_id,
                    declared_species,
                ]
            )
        )


        summary_rows.append({

            "clade":
                clade_name,

            "clade_title":
                config["title"],

            "plot_order":
                plot_order,

            "cre_id":
                cre_id,

            "category":
                category,

            "focal_species":
                focal,

            "discordant_species":
                discordant_species,

            "discordant_state":
                discordant_state,

            "consensus_state":
                consensus_state,

            "contrast_direction":
                contrast_direction(
                    discordant_state
                ),

            "n_turnover_states":
                n_turnover,

            "n_present_states":
                (
                    len(
                        declared_species
                    )
                    - n_turnover
                ),
        })


# ============================================================
# Write summary table
# ============================================================

summary_df = pd.DataFrame(
    summary_rows
)


summary_df.to_csv(
    OUT_TSV,
    sep="\t",
    index=False,
)


# ============================================================
# Plot each clade separately
# ============================================================

written_pngs = []
written_pdfs = []


for clade_name, config in CLADES.items():

    focal = config["focal"]
    declared_species = config["species"]


    selected = (
        clade_data[
            clade_name
        ]["matrix"]
    )


    focal_ids = (
        clade_data[
            clade_name
        ]["focal_ids"]
    )


    secondary_ids = (
        clade_data[
            clade_name
        ]["secondary_ids"]
    )


    species_order = (
        clade_data[
            clade_name
        ]["species_order"]
    )


    n_focal = len(
        focal_ids
    )


    n_secondary = len(
        secondary_ids
    )


    # ========================================================
    # Heatmap matrix
    # ========================================================

    heatmap = (
        selected[
            species_order
        ]
        .T
    )


    invalid_states = sorted(
        set(
            pd.unique(
                heatmap
                .values
                .ravel()
            )
        )
        - set(
            STATE_TO_NUM
        )
    )


    if invalid_states:

        raise SystemExit(
            f"ERROR: invalid states in Tier-1 heatmap "
            f"for {clade_name}:\n"
            + "\n".join(
                invalid_states
            )
        )


    numeric = (
        heatmap
        .replace(
            STATE_TO_NUM
        )
        .astype(float)
    )


    n_species = len(
        species_order
    )


    n_cre = (
        numeric
        .shape[1]
    )


    # ========================================================
    # IMPORTANT: shared y-axis coordinates
    #
    # This is defined BEFORE any plotting code that uses it.
    # ========================================================

    shared_ylim = (
        n_species + 0.5,
        0.5,
    )


    y_edges = (
        np.arange(
            n_species + 1
        )
        + 0.5
    )


    # ========================================================
    # Species labels
    # ========================================================

    species_labels = [
        short_species_name(
            species,
            manifest,
        )
        for species
        in species_order
    ]


    # ========================================================
    # Climate values
    # ========================================================

    climate_numeric = []


    for species in species_order:

        tree_name = (
            slug_to_tree[
                species
            ]
        )


        climate = (
            tree_to_climate.get(
                tree_name,
                "",
            )
        )


        if not climate:

            raise SystemExit(
                f"ERROR: missing climate annotation "
                f"for {species}."
            )


        if climate not in CLIMATE_TO_NUM:

            raise SystemExit(
                f"ERROR: unknown climatic zone "
                f"'{climate}' for {species}."
            )


        climate_numeric.append(
            CLIMATE_TO_NUM[
                climate
            ]
        )


    climate_numeric = (
        np.asarray(
            climate_numeric,
            dtype=float,
        )
        .reshape(
            -1,
            1,
        )
    )


    # ========================================================
    # Figure and axes
    # ========================================================

    fig = plt.figure(
        figsize=FIGSIZE
    )


    grid = fig.add_gridspec(
        nrows=1,
        ncols=5,
        width_ratios=GRID_WIDTH_RATIOS,
        wspace=GRID_WSPACE,
    )


    ax_tree = fig.add_subplot(
        grid[
            0,
            0,
        ]
    )


    ax_labels = fig.add_subplot(
        grid[
            0,
            1,
        ]
    )


    ax_climate = fig.add_subplot(
        grid[
            0,
            2,
        ]
    )


    # Dedicated blank column between climatic strip and heatmap.
    ax_spacer = fig.add_subplot(
        grid[
            0,
            3,
        ]
    )

    ax_spacer.axis(
        "off"
    )


    ax_heat = fig.add_subplot(
        grid[
            0,
            4,
        ]
    )


    # ========================================================
    # Draw phylogeny
    # ========================================================

    tree = prune_tree_to_species(
        TREE_FILE,
        declared_species,
        slug_to_tree,
    )


    # Species labels are drawn in the dedicated neighbouring axis.
    for tip in tree.get_terminals():
        tip.name = ""


    Phylo.draw(
        tree,
        axes=ax_tree,
        do_show=False,
        show_confidence=False,
    )


    ax_tree.set_ylim(
        *shared_ylim
    )


    depths = tree.depths()


    terminal_depths = [
        depths[
            tip
        ]
        for tip
        in tree.get_terminals()
    ]


    max_tree_depth = max(
        terminal_depths
    )


    # Use almost the complete x-axis width.
    # Branch-length ratios themselves remain unchanged.
    ax_tree.set_xlim(
        0,
        max_tree_depth
        * TREE_RIGHT_PADDING_FACTOR,
    )


    ax_tree.set_ylabel(
        ""
    )


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


    for spine in [
        "top",
        "right",
        "left",
    ]:

        ax_tree.spines[
            spine
        ].set_visible(
            False
        )


    # ========================================================
    # Species labels
    # ========================================================

    ax_labels.set_xlim(
        0,
        1,
    )


    ax_labels.set_ylim(
        *shared_ylim
    )


    for row_index, (
        species,
        label,
    ) in enumerate(
        zip(
            species_order,
            species_labels,
        ),
        start=1,
    ):


        # Right-aligned toward climate strip:
        # preserves species-name <-> climate-strip spacing.
        ax_labels.text(
            SPECIES_LABEL_X,
            row_index,
            label,
            ha="right",
            va="center",
            fontsize=SPECIES_FONTSIZE,
            fontweight=(
                "bold"
                if species == focal
                else "normal"
            ),
        )


    ax_labels.set_xticks(
        []
    )


    ax_labels.set_yticks(
        []
    )


    for spine in (
        ax_labels
        .spines
        .values()
    ):

        spine.set_visible(
            False
        )


    # ========================================================
    # Continuous climate strip
    # ========================================================

    climate_x_edges = np.array([
        -0.5,
        0.5,
    ])


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


    ax_climate.set_xlim(
        -0.5,
        0.5,
    )


    ax_climate.set_ylim(
        *shared_ylim
    )


    ax_climate.set_xticks(
        []
    )


    ax_climate.set_yticks(
        []
    )


    for spine in (
        ax_climate
        .spines
        .values()
    ):

        spine.set_visible(
            False
        )


    # ========================================================
    # CRE-state heatmap
    # ========================================================

    cre_x_edges = (
        np.arange(
            n_cre + 1
        )
        - 0.5
    )


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


    ax_heat.set_xlim(
        -0.5,
        n_cre - 0.5,
    )


    ax_heat.set_ylim(
        *shared_ylim
    )


    ax_heat.set_yticks(
        []
    )


    ax_heat.set_ylabel(
        ""
    )


    # ========================================================
    # CRE labels
    # ========================================================

    ax_heat.set_xticks(
        np.arange(
            n_cre
        )
    )


    ax_heat.set_xticklabels(
        numeric.columns,
        rotation=90,
        fontsize=CRE_LABEL_FONTSIZE,
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


    # ========================================================
    # Focal Tier-1 block outline
    # ========================================================

    if n_focal > 0:

        focal_rectangle = Rectangle(
            (
                -0.5,
                0.5,
            ),
            n_focal,
            n_species,
            fill=False,
            edgecolor="#222222",
            linewidth=FOCAL_BLOCK_LINEWIDTH,
            zorder=10,
        )


        ax_heat.add_patch(
            focal_rectangle
        )


    # ========================================================
    # Focal / Secondary block boundary
    # ========================================================

    if (
        n_focal > 0
        and n_secondary > 0
    ):

        ax_heat.axvline(
            n_focal - 0.5,
            color="#222222",
            linewidth=FOCAL_BLOCK_LINEWIDTH,
            zorder=11,
        )


    # ========================================================
    # Block labels
    # ========================================================

    if n_focal > 0:

        focal_center = (
            n_focal - 1
        ) / 2


        ax_heat.text(
            focal_center,
            1.045,
            "Focal Tier 1",
            transform=(
                ax_heat
                .get_xaxis_transform()
            ),
            ha="center",
            va="bottom",
            fontsize=BLOCK_LABEL_FONTSIZE,
            fontweight="bold",
        )


    if n_secondary > 0:

        secondary_start = n_focal
        secondary_end = n_cre - 1

        secondary_center = (
            secondary_start
            + secondary_end
        ) / 2


        ax_heat.text(
            secondary_center,
            1.045,
            "Secondary Tier 1",
            transform=(
                ax_heat
                .get_xaxis_transform()
            ),
            ha="center",
            va="bottom",
            fontsize=BLOCK_LABEL_FONTSIZE,
        )


    # ========================================================
    # Clade title
    # ========================================================

    fig.suptitle(
        config["title"],
        fontsize=TITLE_FONTSIZE,
        fontweight="bold",
        y=0.975,
    )


    # ========================================================
    # Fixed margins
    # ========================================================

    fig.subplots_adjust(
        left=FIG_LEFT,
        right=FIG_RIGHT,
        top=FIG_TOP,
        bottom=FIG_BOTTOM,
    )


    # ========================================================
    # Output files
    # ========================================================

    stem = (
        "focal_clade_tier1_heatmap_"
        f"{clade_name}"
    )


    out_png = (
        FIG_DIR
        / f"{stem}.png"
    )


    out_pdf = (
        FIG_DIR
        / f"{stem}.pdf"
    )


    # Deliberately no bbox_inches="tight".
    # This keeps all exported figures exactly the same size.

    fig.savefig(
        out_png,
        dpi=300,
        facecolor="white",
    )


    fig.savefig(
        out_pdf,
        facecolor="white",
    )


    plt.close(
        fig
    )


    written_pngs.append(
        out_png
    )


    written_pdfs.append(
        out_pdf
    )


    # ========================================================
    # Per-clade QC
    # ========================================================

    print()
    print("=" * 72)
    print(
        config["title"]
    )
    print("=" * 72)


    print(
        "Tree order: "
        + ", ".join(
            species_order
        )
    )


    print(
        f"Displayed CREs: "
        f"{n_cre} "
        f"({n_focal} focal Tier 1 + "
        f"{n_secondary} secondary Tier 1)"
    )


    print(
        "CRE display order:"
    )


    for display_index, cre_id in enumerate(
        selected.index,
        start=1,
    ):

        n_turnover = (
            count_turnover_states(
                selected.loc[
                    cre_id,
                    declared_species,
                ]
            )
        )


        category = (
            "focal"
            if cre_id in focal_ids
            else "secondary"
        )


        print(
            f"  {display_index:>2}. "
            f"{cre_id:<18} "
            f"{category:<10} "
            f"turnover states = {n_turnover}"
        )


    print()


    print(
        f"Wrote PNG:\n"
        f"{out_png}"
    )


    print(
        f"Wrote PDF:\n"
        f"{out_pdf}"
    )


# ============================================================
# Final summary
# ============================================================

print()
print("=" * 72)
print("FOCAL-CLADE HEATMAPS COMPLETE")
print("=" * 72)


print(
    f"Focal clades: "
    f"{len(CLADES)}"
)


print(
    f"PNG figures written: "
    f"{len(written_pngs)}"
)


print(
    f"PDF figures written: "
    f"{len(written_pdfs)}"
)


print()


print(
    f"Fixed figure dimensions: "
    f"{FIG_WIDTH} x "
    f"{FIG_HEIGHT} inches"
)


print(
    "Heatmap renderer: "
    "pcolormesh"
)


print(
    "PDF-safe categorical rendering: YES"
)


print(
    "Continuous climate strip: YES"
)


print(
    "Dedicated climate-to-heatmap spacer: YES"
)


print(
    f"Climate-to-heatmap gap width ratio: "
    f"{CLIMATE_HEATMAP_GAP_WIDTH}"
)


print(
    "Species labels right-aligned before climate strip: YES"
)


print(
    "Tree horizontally expanded: YES"
)


print(
    "Tree branch lengths modified: NO"
)


print(
    "Panel letters A-F: NO"
)


print(
    "Legends in individual figures: NO"
)


print()


print(
    "CRE ordering:"
)


print(
    "  Focal Tier 1 -> Secondary Tier 1"
)


print(
    "  within blocks: "
    "increasing number of turnover_candidate states"
)


print(
    "  Secondary tie-break: "
    "discordant species in displayed tree order -> CRE ID"
)


print()


print(
    f"Wrote plotted-CRE table:\n"
    f"{OUT_TSV}"
)


print(
    f"Total clade-specific CRE occurrences: "
    f"{len(summary_df)}"
)
