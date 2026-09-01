#!/usr/bin/env python3

# ============================================================
# Plot species-level CRE turnover by climatic zone
#
# Purpose:
#   Compare species-level CRE turnover rates among climatic
#   zones and generate descriptive statistics and an
#   exploratory Kruskal-Wallis test.
#
# Primary metric:
#   turnover_rate_evaluable =
#       turnover_candidate / mapped reference CREs
#
# Notes:
#   The Kruskal-Wallis test is exploratory because species are
#   not phylogenetically independent. Phylogenetically
#   controlled inference is performed separately using PGLS.
#
# Input/output paths are supplied by the climate-analysis
# configuration via the pipeline wrapper.
# ============================================================

from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal


# ============================================================
# Analysis setup
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
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Analyze and plot species-level CRE turnover "
            "across climatic zones."
        )
    )

    parser.add_argument("--species-summary", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--traits", type=Path, required=True)
    parser.add_argument("--out-summary", type=Path, required=True)
    parser.add_argument("--out-descriptive", type=Path, required=True)
    parser.add_argument("--out-kruskal", type=Path, required=True)
    parser.add_argument("--out-png", type=Path, required=True)
    parser.add_argument("--out-pdf", type=Path, required=True)

    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_file(path):
    if not path.is_file():
        raise SystemExit(
            f"ERROR: required file not found:\n{path}"
        )


def require_columns(df, required, label):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(
            f"ERROR: {label} missing columns:\n"
            + "\n".join(missing)
        )


def tree_name_from_species(scientific_name):
    """Convert a scientific species name to the tree-name format."""

    return (
        str(scientific_name)
        .strip()
        .upper()
        .replace(" ", "_")
    )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    for path in [
        args.species_summary,
        args.manifest,
        args.traits,
    ]:
        require_file(path)

    for path in [
        args.out_summary,
        args.out_descriptive,
        args.out_kruskal,
        args.out_png,
        args.out_pdf,
    ]:
        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )


    # ========================================================
    # Species-level CRE turnover
    # ========================================================

    summary = pd.read_csv(args.species_summary, sep="\t")

    require_columns(
        summary,
        {
            "species",
            "mapped",
            "turnover_candidate",
            METRIC,
        },
        "species summary",
    )

    if summary["species"].duplicated().any():
        raise SystemExit(
            "ERROR: duplicate species in species summary."
        )

    summary["species"] = (
        summary["species"]
        .astype(str)
        .str.strip()
    )

    summary["mapped"] = pd.to_numeric(
        summary["mapped"],
        errors="raise",
    )

    summary["turnover_candidate"] = pd.to_numeric(
        summary["turnover_candidate"],
        errors="raise",
    )

    summary[METRIC] = pd.to_numeric(
        summary[METRIC],
        errors="raise",
    )

    # Verify the metric used for all downstream climate analyses.
    expected_rate = (
        summary["turnover_candidate"]
        / summary["mapped"]
    )

    if not np.allclose(
        summary[METRIC],
        expected_rate,
        atol=1e-6,
        rtol=0,
    ):
        raise SystemExit(
            "ERROR: turnover_rate_evaluable does not match "
            "turnover_candidate / mapped."
        )


    # ========================================================
    # Species manifest
    # ========================================================

    manifest = pd.read_csv(
        args.manifest,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    require_columns(
        manifest,
        {
            "slug",
            "species",
            "source",
        },
        "combined manifest",
    )

    if manifest["slug"].duplicated().any():
        raise SystemExit(
            "ERROR: duplicate species slugs in manifest."
        )

    manifest["tree_name"] = (
        manifest["species"]
        .map(tree_name_from_species)
    )

    manifest = (
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
                "species":
                    "scientific_name",
            }
        )
    )


    # ========================================================
    # Climatic-zone annotation
    # ========================================================

    traits = pd.read_csv(
        args.traits,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    require_columns(
        traits,
        {
            "tree_name",
            "climatic_zone",
        },
        "species traits",
    )

    if traits["tree_name"].duplicated().any():
        raise SystemExit(
            "ERROR: duplicate tree names in species traits."
        )

    # species_summary uses species slugs, while climatic-zone
    # annotations are linked through scientific/tree names.
    df = summary.merge(
        manifest,
        left_on="species",
        right_on="slug",
        how="left",
        validate="one_to_one",
    )

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

    if df["climatic_zone"].isna().any():
        bad = df.loc[
            df["climatic_zone"].isna(),
            "species",
        ]

        raise SystemExit(
            "ERROR: missing climatic-zone annotation for:\n"
            + "\n".join(bad)
        )

    unknown_climates = sorted(
        set(df["climatic_zone"])
        - set(CLIMATE_ORDER)
    )

    if unknown_climates:
        raise SystemExit(
            "ERROR: unknown climatic-zone values:\n"
            + "\n".join(unknown_climates)
        )


    # ========================================================
    # Sort and write joined species table
    # ========================================================

    climate_rank = {
        climate: rank
        for rank, climate
        in enumerate(CLIMATE_ORDER)
    }

    df["_climate_rank"] = (
        df["climatic_zone"]
        .map(climate_rank)
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

    df.to_csv(
        args.out_summary,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Descriptive statistics
    # ========================================================

    descriptive_rows = []

    for climate in CLIMATE_ORDER:

        values = (
            df.loc[
                df["climatic_zone"] == climate,
                METRIC,
            ]
            .dropna()
            .to_numpy(
                dtype=float
            )
        )

        if len(values) == 0:
            continue

        descriptive_rows.append({
            "climatic_zone": climate,
            "climate_label": CLIMATE_LABELS[climate],
            "n_species": len(values),
            "mean": float(np.mean(values)),
            "median": float(np.median(values)),
            "sd":(float(np.std(values, ddof=1)) if len(values) > 1 else np.nan),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
        })

    descriptive = pd.DataFrame(descriptive_rows)

    descriptive.to_csv(
        args.out_descriptive,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Exploratory Kruskal-Wallis test
    # ========================================================

    groups = []
    used_climates = []

    for climate in CLIMATE_ORDER:

        values = (
            df.loc[
                df["climatic_zone"] == climate,
                METRIC,
            ]
            .dropna()
            .to_numpy(
                dtype=float
            )
        )

        if len(values) > 0:
            groups.append(values)
            used_climates.append(climate)

    if len(groups) >= 2:

        test = kruskal(*groups)
        statistic = float(test.statistic)
        p_value = float(test.pvalue)

    else:

        statistic = np.nan
        p_value = np.nan

    kruskal_df = pd.DataFrame([{
        "test": "Kruskal-Wallis",
        "metric": METRIC,
        "groups": ",".join(used_climates),
        "statistic": statistic,
        "p_value": p_value,
        "phylogenetically_corrected": False,
    }])

    kruskal_df.to_csv(
        args.out_kruskal,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Plot
    # ========================================================

    fig, ax = plt.subplots(figsize=(7.8, 5.8))

    x_positions = {
        climate: index
        for index, climate
        in enumerate(CLIMATE_ORDER)
    }

    for climate in CLIMATE_ORDER:

        sub = (
            df.loc[
                df["climatic_zone"] == climate
            ]
            .sort_values(
                "scientific_name"
            )
        )

        if sub.empty:
            continue

        # Deterministic horizontal offsets avoid random jitter
        # and keep the figure exactly reproducible.
        if len(sub) == 1:

            offsets = np.array([0.0])

        else:

            offsets = np.linspace(-0.16, 0.16, len(sub))

        x = (x_positions[climate] + offsets)

        y = (sub[METRIC].to_numpy(dtype=float) * 100)

        ax.scatter(
            x,
            y,
            s=60,
            facecolor=CLIMATE_COLORS[climate],
            edgecolor="#222222",
            linewidth=0.6,
            alpha=0.90,
            zorder=3,
        )

        # Horizontal line indicates the climatic-zone median.
        median = float(
            np.median(y)
        )

        ax.plot(
            [
                x_positions[climate] - 0.22,
                x_positions[climate] + 0.22,
            ],
            [
                median,
                median,
            ],
            color="#222222",
            linewidth=2.0,
            zorder=4,
        )


    # ========================================================
    # Figure formatting
    # ========================================================

    ax.set_xticks(
        range(
            len(CLIMATE_ORDER)
        )
    )

    ax.set_xticklabels(
        [
            CLIMATE_LABELS[climate]
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


    # ========================================================
    # Save
    # ========================================================

    fig.savefig(
        args.out_png,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        args.out_pdf,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )


    # ========================================================
    # Console summary
    # ========================================================

    print()
    print("=" * 72)
    print("SPECIES-LEVEL TURNOVER BY CLIMATIC ZONE")
    print("=" * 72)
    print(
        f"Species analyzed: "
        f"{len(df)}"
    )

    print()
    for row in descriptive.itertuples():
        print(
            f"{row.climatic_zone:4s}  "
            f"n={row.n_species:2d}  "
            f"mean={row.mean:.4f}  "
            f"median={row.median:.4f}  "
            f"range={row.min:.4f}-{row.max:.4f}"
        )
    print()
    print("Exploratory Kruskal-Wallis test:")
    print(kruskal_df.to_string(index=False))
    print()
    print(
        "NOTE: Kruskal-Wallis does not account for "
        "phylogenetic non-independence."
    )
    print(
        "Use the PGLS analysis for phylogenetically "
        "controlled inference."
    )
    print()
    print(
        f"Wrote PNG:\n"
        f"{args.out_png}"
    )
    print(
        f"Wrote PDF:\n"
        f"{args.out_pdf}"
    )


if __name__ == "__main__":
    main()
