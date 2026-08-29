#!/usr/bin/env python3

"""
Plot species-level CRE turnover rates by climatic zone.

Input
-----
cre_classification/results/species_summary.tsv

This table already contains the species-level CRE-state counts and rates,
including:

    turnover_candidate
    mapped
    uncertain
    turnover_rate_evaluable

The column called "species" in species_summary.tsv contains species slugs
(e.g. dazt, dbif, d_repleta), not scientific names.

Additional metadata are joined from:

    external_scrmshaw/combined_manifest.tsv
    cre_classification/phylogeny/data/species_traits.tsv

Primary metric
--------------
turnover_rate_evaluable

This corresponds to:

    turnover_candidate / mapped

where "mapped" is the number of evaluable reference CREs.

Outputs
-------
downstream_analyses/climate_analysis/results/
    climate_turnover_summary.tsv
    climate_turnover_descriptive_statistics.tsv
    climate_turnover_kruskal.tsv

downstream_analyses/climate_analysis/results/figures/
    turnover_rate_by_climate.png
    turnover_rate_by_climate.pdf

The Kruskal-Wallis test is exploratory only because it does not account
for phylogenetic non-independence. The phylogenetically controlled
inference is performed separately using PGLS.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

from scipy.stats import kruskal


# ============================================================
# Paths
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPT_DIR.parent
DOWNSTREAM_DIR = ANALYSIS_DIR.parent
PROJECT_DIR = DOWNSTREAM_DIR.parent

CLASS_DIR = PROJECT_DIR / "cre_classification"

SPECIES_SUMMARY_FILE = (
    CLASS_DIR
    / "results"
    / "species_summary.tsv"
)

MANIFEST_FILE = (
    PROJECT_DIR
    / "external_scrmshaw"
    / "combined_manifest.tsv"
)

TRAITS_FILE = (
    CLASS_DIR
    / "phylogeny"
    / "data"
    / "species_traits.tsv"
)

RESULTS_DIR = (
    ANALYSIS_DIR
    / "results"
)

FIG_DIR = (
    RESULTS_DIR
    / "figures"
)

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


OUT_SUMMARY = (
    RESULTS_DIR
    / "climate_turnover_summary.tsv"
)

OUT_DESCRIPTIVE = (
    RESULTS_DIR
    / "climate_turnover_descriptive_statistics.tsv"
)

OUT_KRUSKAL = (
    RESULTS_DIR
    / "climate_turnover_kruskal.tsv"
)

OUT_PNG = (
    FIG_DIR
    / "turnover_rate_by_climate.png"
)


# ============================================================
# Constants
# ============================================================

METRIC = "turnover_rate_evaluable"

CLIMATE_ORDER = [
    "TROP",
    "ARID",
    "TEMP",
    "BORE",
]

CLIMATE_LABELS = {
    "TROP": "Tropical",
    "ARID": "Arid",
    "TEMP": "Temperate",
    "BORE": "Boreal",
}

CLIMATE_COLORS = {
    "TROP": "#E41A1C",
    "ARID": "#FBC02D",
    "TEMP": "#009E73",
    "BORE": "#0072B2",
}


# ============================================================
# Helpers
# ============================================================

def require_file(path):

    if not path.exists():

        raise SystemExit(
            f"ERROR: required file not found:\n{path}"
        )


def tree_name_from_species(
    scientific_name,
):

    return (
        str(scientific_name)
        .strip()
        .upper()
        .replace(" ", "_")
    )


# ============================================================
# Input checks
# ============================================================

for path in [
    SPECIES_SUMMARY_FILE,
    MANIFEST_FILE,
    TRAITS_FILE,
]:

    require_file(
        path
    )


# ============================================================
# Load species summary
#
# Important:
# The column "species" contains SLUGS.
# ============================================================

summary = pd.read_csv(
    SPECIES_SUMMARY_FILE,
    sep="\t",
)


required_summary_columns = {
    "species",
    "n_reference_cres",
    "mapped",
    "unmapped",
    "mapping_rate",
    "present",
    "turnover_candidate",
    "no_detected_CRE",
    "uncertain",
    "turnover_rate_evaluable",
}


missing = (
    required_summary_columns
    - set(summary.columns)
)


if missing:

    raise SystemExit(
        "ERROR: species_summary.tsv missing columns:\n"
        + "\n".join(
            sorted(
                missing
            )
        )
    )


summary["species"] = (
    summary["species"]
    .astype(str)
    .str.strip()
)


if summary["species"].duplicated().any():

    duplicates = (
        summary.loc[
            summary["species"].duplicated(
                keep=False
            ),
            "species",
        ]
        .unique()
        .tolist()
    )

    raise SystemExit(
        "ERROR: duplicate species slugs "
        "in species_summary.tsv:\n"
        + "\n".join(
            sorted(
                duplicates
            )
        )
    )


# ============================================================
# Validate species-summary arithmetic
# ============================================================

expected_mapped = (
    summary["present"]
    + summary["turnover_candidate"]
    + summary["no_detected_CRE"]
)


bad_mapped = (
    expected_mapped
    != summary["mapped"]
)


if bad_mapped.any():

    bad = summary.loc[
        bad_mapped,
        [
            "species",
            "mapped",
            "present",
            "turnover_candidate",
            "no_detected_CRE",
        ],
    ]

    raise SystemExit(
        "ERROR: mapped does not equal "
        "present + turnover_candidate + no_detected_CRE:\n"
        + bad.to_string(
            index=False
        )
    )


expected_total = (
    summary["mapped"]
    + summary["uncertain"]
)


bad_total = (
    expected_total
    != summary["n_reference_cres"]
)


if bad_total.any():

    bad = summary.loc[
        bad_total,
        [
            "species",
            "n_reference_cres",
            "mapped",
            "uncertain",
        ],
    ]

    raise SystemExit(
        "ERROR: mapped + uncertain does not equal "
        "n_reference_cres:\n"
        + bad.to_string(
            index=False
        )
    )


# ============================================================
# Validate turnover_rate_evaluable
# ============================================================

recomputed_turnover_rate = (
    summary["turnover_candidate"]
    / summary["mapped"]
)


rate_difference = (
    recomputed_turnover_rate
    - summary[METRIC]
).abs()


bad_rate = (
    rate_difference
    > 1e-6
)


if bad_rate.any():

    bad = summary.loc[
        bad_rate,
        [
            "species",
            "turnover_candidate",
            "mapped",
            METRIC,
        ],
    ].copy()


    bad[
        "recomputed_turnover_rate"
    ] = recomputed_turnover_rate.loc[
        bad_rate
    ].values


    raise SystemExit(
        "ERROR: turnover_rate_evaluable does not match "
        "turnover_candidate / mapped:\n"
        + bad.to_string(
            index=False
        )
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


missing = (
    required_manifest_columns
    - set(manifest.columns)
)


if missing:

    raise SystemExit(
        "ERROR: combined_manifest.tsv missing columns:\n"
        + "\n".join(
            sorted(
                missing
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


# Scientific species name -> tree-style name.
manifest["tree_name"] = (
    manifest["species"]
    .map(
        tree_name_from_species
    )
)


# Rename scientific-name column to avoid ambiguity with the
# "species" slug column in species_summary.tsv.
manifest_small = (
    manifest[
        [
            "slug",
            "species",
            "source",
            "tree_name",
        ]
    ]
    .rename(
        columns={
            "species": "scientific_name"
        }
    )
    .copy()
)


if manifest_small["slug"].duplicated().any():

    duplicates = (
        manifest_small.loc[
            manifest_small["slug"].duplicated(
                keep=False
            ),
            "slug",
        ]
        .unique()
        .tolist()
    )

    raise SystemExit(
        "ERROR: duplicate slugs "
        "in combined_manifest.tsv:\n"
        + "\n".join(
            sorted(
                duplicates
            )
        )
    )


if manifest_small["tree_name"].duplicated().any():

    duplicates = (
        manifest_small.loc[
            manifest_small[
                "tree_name"
            ].duplicated(
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
            sorted(
                duplicates
            )
        )
    )


# ============================================================
# Merge species summary with manifest
#
# species_summary.species = slug
# manifest_small.slug      = slug
# ============================================================

df = summary.merge(
    manifest_small,
    left_on="species",
    right_on="slug",
    how="left",
    validate="one_to_one",
)


missing_manifest = (
    df["tree_name"].isna()
    | (
        df["tree_name"]
        .astype(str)
        .str.strip()
        == ""
    )
)


if missing_manifest.any():

    missing_species = (
        df.loc[
            missing_manifest,
            "species",
        ]
        .tolist()
    )

    raise SystemExit(
        "ERROR: species slugs missing from "
        "combined_manifest.tsv:\n"
        + "\n".join(
            missing_species
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


required_trait_columns = {
    "tree_name",
    "climatic_zone",
}


missing = (
    required_trait_columns
    - set(traits.columns)
)


if missing:

    raise SystemExit(
        "ERROR: species_traits.tsv missing columns:\n"
        + "\n".join(
            sorted(
                missing
            )
        )
    )


traits["tree_name"] = (
    traits["tree_name"]
    .astype(str)
    .str.strip()
)


traits["climatic_zone"] = (
    traits["climatic_zone"]
    .astype(str)
    .str.strip()
)


if traits["tree_name"].duplicated().any():

    duplicates = (
        traits.loc[
            traits["tree_name"].duplicated(
                keep=False
            ),
            "tree_name",
        ]
        .unique()
        .tolist()
    )

    raise SystemExit(
        "ERROR: duplicate tree names "
        "in species_traits.tsv:\n"
        + "\n".join(
            sorted(
                duplicates
            )
        )
    )


# ============================================================
# Merge climate annotation
# ============================================================

df = df.merge(
    traits[
        [
            "tree_name",
            "climatic_zone",
        ]
    ],
    on="tree_name",
    how="left",
    validate="many_to_one",
)


missing_climate = (
    df["climatic_zone"].isna()
    | (
        df["climatic_zone"]
        .astype(str)
        .str.strip()
        == ""
    )
)


if missing_climate.any():

    missing_species = (
        df.loc[
            missing_climate,
            [
                "species",
                "scientific_name",
                "tree_name",
            ],
        ]
    )

    raise SystemExit(
        "ERROR: missing climate annotation for:\n"
        + missing_species.to_string(
            index=False
        )
    )


unknown_climate = sorted(
    set(
        df["climatic_zone"]
    )
    - set(
        CLIMATE_ORDER
    )
)


if unknown_climate:

    raise SystemExit(
        "ERROR: unknown climatic-zone values:\n"
        + "\n".join(
            unknown_climate
        )
    )


# ============================================================
# Sort output
# ============================================================

climate_rank = {
    climate: rank
    for rank, climate
    in enumerate(
        CLIMATE_ORDER
    )
}


df["_climate_rank"] = (
    df["climatic_zone"]
    .map(
        climate_rank
    )
)


df = (
    df
    .sort_values(
        [
            "_climate_rank",
            "scientific_name",
        ],
        kind="mergesort",
    )
    .drop(
        columns="_climate_rank"
    )
    .reset_index(
        drop=True
    )
)


# ============================================================
# Write joined species-level table
# ============================================================

df.to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False,
)


# ============================================================
# Descriptive statistics by climate
# ============================================================

descriptive_rows = []


for climate in CLIMATE_ORDER:

    values = (
        df.loc[
            df["climatic_zone"]
            == climate,
            METRIC,
        ]
        .astype(float)
        .to_numpy()
    )


    if len(values) == 0:
        continue


    descriptive_rows.append({
        "climatic_zone": climate,
        "climate_label": CLIMATE_LABELS[
            climate
        ],
        "n_species": len(
            values
        ),
        "mean": float(
            np.mean(
                values
            )
        ),
        "median": float(
            np.median(
                values
            )
        ),
        "sd": (
            float(
                np.std(
                    values,
                    ddof=1,
                )
            )
            if len(values) > 1
            else np.nan
        ),
        "min": float(
            np.min(
                values
            )
        ),
        "max": float(
            np.max(
                values
            )
        ),
    })


descriptive_df = pd.DataFrame(
    descriptive_rows
)


descriptive_df.to_csv(
    OUT_DESCRIPTIVE,
    sep="\t",
    index=False,
)


# ============================================================
# Exploratory Kruskal-Wallis test
# ============================================================

groups = []

used_climates = []


for climate in CLIMATE_ORDER:

    values = (
        df.loc[
            df["climatic_zone"]
            == climate,
            METRIC,
        ]
        .dropna()
        .astype(float)
        .to_numpy()
    )


    if len(values) > 0:

        groups.append(
            values
        )

        used_climates.append(
            climate
        )


if len(groups) >= 2:

    test = kruskal(
        *groups
    )


    kruskal_df = pd.DataFrame([
        {
            "test": "Kruskal-Wallis",
            "metric": METRIC,
            "groups": ",".join(
                used_climates
            ),
            "statistic": float(
                test.statistic
            ),
            "p_value": float(
                test.pvalue
            ),
            "phylogenetically_corrected": False,
        }
    ])

else:

    kruskal_df = pd.DataFrame([
        {
            "test": "Kruskal-Wallis",
            "metric": METRIC,
            "groups": ",".join(
                used_climates
            ),
            "statistic": np.nan,
            "p_value": np.nan,
            "phylogenetically_corrected": False,
        }
    ])


kruskal_df.to_csv(
    OUT_KRUSKAL,
    sep="\t",
    index=False,
)


# ============================================================
# Plot
# ============================================================

fig, ax = plt.subplots(
    figsize=(
        7.8,
        5.8,
    )
)


x_positions = {
    climate: i
    for i, climate
    in enumerate(
        CLIMATE_ORDER
    )
}


# ------------------------------------------------------------
# Plot individual species
#
# Deterministic offsets are used instead of random jitter so
# the plot is exactly reproducible.
# ------------------------------------------------------------

for climate in CLIMATE_ORDER:

    sub = (
        df.loc[
            df["climatic_zone"]
            == climate
        ]
        .sort_values(
            "scientific_name"
        )
    )


    n_species = len(
        sub
    )


    if n_species == 0:
        continue


    if n_species == 1:

        offsets = np.array([
            0.0
        ])

    else:

        offsets = np.linspace(
            -0.16,
            0.16,
            n_species,
        )


    x = (
        x_positions[
            climate
        ]
        + offsets
    )


    y = (
        sub[
            METRIC
        ]
        .astype(float)
        .to_numpy()
        * 100
    )


    ax.scatter(
        x,
        y,
        s=60,
        facecolor=CLIMATE_COLORS[
            climate
        ],
        edgecolor="#222222",
        linewidth=0.6,
        alpha=0.90,
        zorder=3,
    )


    # --------------------------------------------------------
    # Median line
    # --------------------------------------------------------

    median = float(
        np.median(
            y
        )
    )


    ax.plot(
        [
            x_positions[
                climate
            ] - 0.22,

            x_positions[
                climate
            ] + 0.22,
        ],
        [
            median,
            median,
        ],
        color="#222222",
        linewidth=2.0,
        zorder=4,
    )


# ============================================================
# Formatting
# ============================================================

ax.set_xticks(
    range(
        len(
            CLIMATE_ORDER
        )
    )
)


ax.set_xticklabels(
    [
        CLIMATE_LABELS[
            climate
        ]
        for climate in CLIMATE_ORDER
    ],
    fontsize=10,
)


ax.set_ylabel(
    "Turnover candidates among evaluable CRE loci (%)",
    fontsize=11,
)


ax.set_xlabel(
    "Climatic zone",
    fontsize=11,
)


ax.set_title(
    "Species-level CRE turnover across climatic zones",
    fontsize=13,
    pad=10,
)


ax.grid(
    axis="y",
    linestyle=":",
    linewidth=0.6,
    alpha=0.5,
)


ax.set_axisbelow(
    True
)


# Remove unnecessary spines.
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


# ============================================================
# Console summary
# ============================================================

print()
print("=" * 72)
print("SPECIES-LEVEL TURNOVER BY CLIMATE")
print("=" * 72)
print()


print(
    f"Species analysed: "
    f"{len(df)}"
)

print()


for _, row in descriptive_df.iterrows():

    print(
        f"{row['climatic_zone']:4s}  "
        f"n={int(row['n_species']):2d}  "
        f"mean={row['mean']:.4f}  "
        f"median={row['median']:.4f}  "
        f"range="
        f"{row['min']:.4f}-"
        f"{row['max']:.4f}"
    )


print()
print("Exploratory test:")
print(
    kruskal_df.to_string(
        index=False
    )
)


print()
print(
    "NOTE: Kruskal-Wallis does not account for "
    "phylogenetic non-independence."
)

print(
    "Use the PGLS analysis for the primary "
    "phylogenetically corrected inference."
)


print()
print("Wrote:")
print(
    OUT_SUMMARY
)

print(
    OUT_DESCRIPTIVE
)

print(
    OUT_KRUSKAL
)

print(
    OUT_PNG
)
