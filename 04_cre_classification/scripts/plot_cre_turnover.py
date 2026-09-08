#!/usr/bin/env python3

# ============================================================
# Generate CRE-turnover summary figures and QC outputs
#
# Purpose:
#   Visualize species- and CRE-level conservation patterns,
#   mapping quality, and CRE-state composition.
#
# Input/output paths:
#   Supplied by the pipeline wrapper using
#   config/classification_config.sh.
# ============================================================

import argparse
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# ============================================================
# Repository-relative default paths
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
CLASS_DIR = SCRIPT_DIR.parent
PROJECT_ROOT = CLASS_DIR.parent
RESULTS_DIR = CLASS_DIR / "results"

DEFAULT_SUMMARY_FILE = RESULTS_DIR / "species_summary.tsv"
DEFAULT_MATRIX_FILE = RESULTS_DIR / "cre_turnover_matrix.tsv"
DEFAULT_LONG_FILE = RESULTS_DIR / "cre_turnover_all_species.tsv"
DEFAULT_QC_FILE = RESULTS_DIR / "species_qc_summary.tsv"
DEFAULT_TARGETS_FILE = PROJECT_ROOT / "03_pairwise_wga" / "target_species.txt"
DEFAULT_MANIFEST_FILE = PROJECT_ROOT / "01_scrmshaw" / "external" / "combined_manifest.tsv"
DEFAULT_TREE_ORDER_FILE = CLASS_DIR / "phylogeny" / "results" / "species_order_40_tree_names.txt"
DEFAULT_TRAITS_FILE = CLASS_DIR / "phylogeny" / "data" / "species_traits.tsv"
DEFAULT_FIG_DIR = RESULTS_DIR / "figures"
DEFAULT_CRE_SUMMARY_OUT = RESULTS_DIR / "CRE_conservation_summary.tsv"
DEFAULT_PLOT_SUMMARY_OUT = RESULTS_DIR / "species_plot_qc_summary.tsv"


# ============================================================
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate species- and CRE-level CRE-turnover summary figures and QC outputs."
    )
    parser.add_argument("--species-summary", type=Path, default=DEFAULT_SUMMARY_FILE)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX_FILE)
    parser.add_argument("--long-table", type=Path, default=DEFAULT_LONG_FILE)
    parser.add_argument("--species-qc", type=Path, default=DEFAULT_QC_FILE)
    parser.add_argument("--targets", type=Path, default=DEFAULT_TARGETS_FILE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST_FILE)
    parser.add_argument("--tree-order", type=Path, default=DEFAULT_TREE_ORDER_FILE)
    parser.add_argument("--traits", type=Path, default=DEFAULT_TRAITS_FILE)
    parser.add_argument("--fig-dir", type=Path, default=DEFAULT_FIG_DIR)
    parser.add_argument("--cre-summary-out", type=Path, default=DEFAULT_CRE_SUMMARY_OUT)
    parser.add_argument("--plot-summary-out", type=Path, default=DEFAULT_PLOT_SUMMARY_OUT)
    return parser.parse_args()


args = parse_args()
SUMMARY_FILE = Path(args.species_summary)
MATRIX_FILE = Path(args.matrix)
LONG_FILE = Path(args.long_table)
QC_FILE = Path(args.species_qc)
TARGETS_FILE = Path(args.targets)
MANIFEST_FILE = Path(args.manifest)
TREE_ORDER_FILE = Path(args.tree_order)
TRAITS_FILE = Path(args.traits)
FIG_DIR = Path(args.fig_dir)
CRE_SUMMARY_OUT = Path(args.cre_summary_out)
PLOT_SUMMARY_OUT = Path(args.plot_summary_out)

FIG_DIR.mkdir(parents=True, exist_ok=True)
CRE_SUMMARY_OUT.parent.mkdir(parents=True, exist_ok=True)
PLOT_SUMMARY_OUT.parent.mkdir(parents=True, exist_ok=True)


# ============================================================
# Original analysis and plotting logic
# ============================================================
DMEL_SLUG = 'd_melanogaster'
STATE_COLORS = {
    "positional_match": "#0072B2", "turnover_candidate": "#E69F00", "no_detected_CRE": "#E5E5E5",
    "uncertain": "#9E9E9E", "reference": "#FFFFFF",
}
SOURCE_MARKERS = {'generated': 'o', 'external': '^'}
SOURCE_LABELS = {'generated': 'Generated', 'external': 'External'}
EVALUABLE_COLOR = '#0072B2'
CLIMATE_COLORS = {"TROP": "#D55E00", "ARID": "#F0E442", "TEMP": "#009E73",  "BORE": "#56B4E9"}
CLIMATE_LABELS = {'TROP': 'Tropical', 'ARID': 'Arid', 'TEMP': 'Temperate', 'BORE': 'Boreal'}
VALID_CLIMATE_ZONES = set(CLIMATE_COLORS)
STATE_LABELS = {
    "positional_match": "Positional match", "turnover_candidate": "Turnover candidate",
    "no_detected_CRE": "No CRE detected", "uncertain": "Uncertain",
    "reference": r"$D.\ melanogaster$ reference",
}
VALID_STATES = {'positional_match', 'turnover_candidate', 'no_detected_CRE', 'uncertain'}

# ============================================================
# Plot typography and legend layout
# ============================================================

AXIS_LABEL_FONTSIZE = 12
LEGEND_TITLE_FONTSIZE = 12
LEGEND_FONTSIZE = 11

# Right-side legend anchors. Increasing the vertical distance between
# these y-values increases the whitespace between separate legend blocks.
AXES_LEGEND_X = 1.02
AXES_REFERENCE_Y = 1.00
AXES_CLIMATE_Y = 0.76

FIGURE_LEGEND_X = 0.735
FIGURE_REFERENCE_Y = 0.92
FIGURE_SECONDARY_Y = 0.78
FIGURE_CLIMATE_Y_LONG = 0.48   # after a 3-entry secondary legend
FIGURE_CLIMATE_Y_SHORT = 0.55  # after a 2-entry secondary legend

SCATTER_SOURCE_Y = 1.00
SCATTER_CLIMATE_Y = 0.66

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.labelsize": AXIS_LABEL_FONTSIZE,
    "legend.title_fontsize": LEGEND_TITLE_FONTSIZE,
    "legend.fontsize": LEGEND_FONTSIZE,
})

# ============================================================
# Utility functions
# ============================================================

def require_file(path):
    """
    Abort if a required file is missing.
    """
    if not path.exists():
        raise SystemExit(f'ERROR: required file does not exist:\n{path}')

def save_figure(fig, basename):
    """
    Save one figure as PNG and PDF.
    """
    out_png = FIG_DIR / f'{basename}.png'
    out_pdf = FIG_DIR / f'{basename}.pdf'
    fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Wrote: {out_png}')
    print(f'Wrote: {out_pdf}')

def tree_name_from_species(species_name):
    """
    Convert a scientific species name to the naming convention
    used in the published Newick tree.

    Example
    -------
    Drosophila virilis
        -> DROSOPHILA_VIRILIS
    """
    return str(species_name).strip().upper().replace(' ', '_')

def identify_scatter_outliers(df, x_col, y_col, iqr_factor=1.5):
    """
    Identify scatter-plot outliers using the 1.5 x IQR rule
    independently for both axes.

    A species is labelled if its x OR y value lies outside
    the corresponding IQR bounds.
    """
    x = df[x_col]
    y = df[y_col]
    x_q1 = x.quantile(0.25)
    x_q3 = x.quantile(0.75)
    x_iqr = x_q3 - x_q1
    y_q1 = y.quantile(0.25)
    y_q3 = y.quantile(0.75)
    y_iqr = y_q3 - y_q1
    x_low = x_q1 - iqr_factor * x_iqr
    x_high = x_q3 + iqr_factor * x_iqr
    y_low = y_q1 - iqr_factor * y_iqr
    y_high = y_q3 + iqr_factor * y_iqr
    return (x < x_low) | (x > x_high) | (y < y_low) | (y > y_high)

# ============================================================
# Species information
# ============================================================

def load_manifest():
    """
    Load combined_manifest.tsv and construct a temporary
    tree_name column.

    The input manifest itself is NOT modified.
    """
    require_file(MANIFEST_FILE)
    df = pd.read_csv(MANIFEST_FILE, sep='\t', dtype=str).fillna('')
    required = {'slug', 'species', 'source'}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit('ERROR: combined_manifest.tsv lacks columns:\n' + '\n'.join(sorted(missing)))
    for column in ['slug', 'species', 'source']:
        df[column] = df[column].astype(str).str.strip()
    empty_slug = df['slug'] == ''
    empty_species = df['species'] == ''
    if empty_slug.any():
        raise SystemExit(
            "ERROR: combined_manifest.tsv contains empty slugs:\n"
            + df.loc[empty_slug, "species"].to_string(index=False)
        )
    if empty_species.any():
        raise SystemExit('ERROR: combined_manifest.tsv contains empty species names.')
    df['tree_name'] = df['species'].map(tree_name_from_species)
    duplicated_slugs = df.loc[df['slug'].duplicated(keep=False), 'slug'].unique().tolist()
    if duplicated_slugs:
        raise SystemExit(
            "ERROR: duplicate slugs in combined_manifest.tsv:\n" + "\n".join(sorted(duplicated_slugs))
        )
    duplicated_tree_names = df.loc[df['tree_name'].duplicated(keep=False), 'tree_name'].unique().tolist()
    if duplicated_tree_names:
        raise SystemExit(
            "ERROR: duplicate inferred tree names in combined_manifest.tsv:\n"
            + "\n".join(sorted(duplicated_tree_names))
        )
    reference_rows = (df['slug'] == DMEL_SLUG).sum()
    if reference_rows != 1:
        raise SystemExit(
            f"ERROR: expected exactly one {DMEL_SLUG} row in combined_manifest.tsv; found {reference_rows}."
        )
    return df

def load_species_traits(manifest):
    """
    Load climatic-zone information and map tree names to
    pipeline species slugs using combined_manifest.tsv.

    Returns
    -------
    slug_to_climate : dict
        Mapping:
            pipeline slug -> climatic zone

        Example:
            dvir -> TEMP
            dazt -> ARID
    """
    require_file(TRAITS_FILE)
    traits = pd.read_csv(TRAITS_FILE, sep='\t', dtype=str).fillna('')
    required = {'tree_name', 'climatic_zone'}
    missing = required - set(traits.columns)
    if missing:
        raise SystemExit('ERROR: species_traits.tsv lacks columns:\n' + '\n'.join(sorted(missing)))
    traits['tree_name'] = traits['tree_name'].astype(str).str.strip()
    traits['climatic_zone'] = traits['climatic_zone'].astype(str).str.strip().str.upper()
    duplicated = traits.loc[traits['tree_name'].duplicated(keep=False), 'tree_name'].unique().tolist()
    if duplicated:
        raise SystemExit(
            "ERROR: duplicate tree names in species_traits.tsv:\n" + "\n".join(sorted(duplicated))
        )
    unknown_zones = sorted(set(traits['climatic_zone']) - VALID_CLIMATE_ZONES - {''})
    if unknown_zones:
        raise SystemExit('ERROR: unknown climatic zones in species_traits.tsv:\n' + '\n'.join(unknown_zones))
    tree_to_climate = dict(zip(traits['tree_name'], traits['climatic_zone']))
    slug_to_tree = dict(zip(manifest['slug'], manifest['tree_name']))
    slug_to_climate = {}
    for slug, tree_name in slug_to_tree.items():
        if tree_name in tree_to_climate:
            slug_to_climate[slug] = tree_to_climate[tree_name]
    return slug_to_climate

def validate_climate_data(plot_order, slug_to_climate):
    """
    Require climatic-zone information for every species
    displayed in phylogenetic species plots.
    """
    missing = [sp for sp in plot_order if sp not in slug_to_climate or not slug_to_climate[sp]]
    if missing:
        raise SystemExit('ERROR: climatic-zone information missing for:\n' + '\n'.join(missing))

def load_tree_order():
    """
    Read the 40 species in the order in which they occur
    in the pruned phylogenetic tree.
    """
    require_file(TREE_ORDER_FILE)
    tree_order = [x.strip() for x in TREE_ORDER_FILE.read_text().splitlines() if x.strip()]
    if not tree_order:
        raise SystemExit('ERROR: species_order_40_tree_names.txt is empty.')
    if len(tree_order) != 40:
        raise SystemExit(
            f"ERROR: expected 40 species in species_order_40_tree_names.txt, found {len(tree_order)}."
        )
    if len(tree_order) != len(set(tree_order)):
        raise SystemExit('ERROR: duplicate species in species_order_40_tree_names.txt.')
    return tree_order

def load_targets():
    """
    Read the species that currently have CRE-classification
    results.
    """
    require_file(TARGETS_FILE)
    targets = [x.strip() for x in TARGETS_FILE.read_text().splitlines() if x.strip()]
    if not targets:
        raise SystemExit('ERROR: target_species.txt is empty.')
    if len(targets) != len(set(targets)):
        raise SystemExit('ERROR: duplicate slugs in target_species.txt.')
    if DMEL_SLUG in targets:
        raise SystemExit('ERROR: D. melanogaster must not occur in target_species.txt.')
    return targets

def build_species_orders(manifest, tree_order, targets):
    """
    Map tree names to pipeline slugs.

    Returns
    -------
    target_order
        Analysed target species in phylogenetic order.

    plot_order
        Same list with D. melanogaster inserted at its
        phylogenetic position.
    """
    tree_to_slug = dict(zip(manifest['tree_name'], manifest['slug']))
    missing_manifest = [tree_name for tree_name in tree_order if tree_name not in tree_to_slug]
    if missing_manifest:
        raise SystemExit(
            "ERROR: tree species missing from combined_manifest.tsv:\n" + "\n".join(missing_manifest)
        )
    manifest_slugs = set(manifest['slug'])
    missing_targets_manifest = sorted(set(targets) - manifest_slugs)
    if missing_targets_manifest:
        raise SystemExit(
            "ERROR: target species missing from combined_manifest.tsv:\n"
            + "\n".join(missing_targets_manifest)
        )
    full_slug_order = [tree_to_slug[tree_name] for tree_name in tree_order]
    tree_slugs = set(full_slug_order)
    missing_targets_tree = sorted(set(targets) - tree_slugs)
    if missing_targets_tree:
        raise SystemExit(
            "ERROR: target species missing from 40-species phylogenetic tree:\n"
            + "\n".join(missing_targets_tree)
        )
    target_set = set(targets)
    target_order = [slug for slug in full_slug_order if slug in target_set]
    if len(target_order) != len(targets):
        raise SystemExit('ERROR: unexpected number of target species after phylogenetic ordering.')
    wanted_for_plot = target_set | {DMEL_SLUG}
    plot_order = [slug for slug in full_slug_order if slug in wanted_for_plot]
    if DMEL_SLUG not in plot_order:
        raise SystemExit('ERROR: D. melanogaster is missing from the 40-species phylogenetic tree.')
    if len(plot_order) != len(target_order) + 1:
        raise SystemExit('ERROR: unexpected number of species in plot order.')
    return (target_order, plot_order)

def short_species_name(slug, manifest):
    """
    Return abbreviated scientific name for plot labels.

    External SCRMshaw species and the D. melanogaster reference are marked with an asterisk.

    Examples
    --------
    dvir
        -> D. virilis

    d_mojavensis
        -> D. mojavensis*
    """
    hit = manifest.loc[manifest['slug'] == slug]
    if len(hit) != 1:
        return slug
    row = hit.iloc[0]
    name = row['species']
    source = row['source']
    if name.startswith('Drosophila '):
        label = 'D. ' + name.split(' ', 1)[1]
    elif name.startswith('Zaprionus '):
        label = 'Z. ' + name.split(' ', 1)[1]
    else:
        label = name
    if source == 'external' or slug == DMEL_SLUG:
        label += '*'
    return label

# ============================================================
# Plot helpers
# ============================================================

def insert_reference_placeholder(values_by_species, plot_order):
    """
    Insert NaN at the Dmel position.

    Dmel therefore receives no numerical value and does not
    contribute to any target-species calculation.
    """
    values = []
    for sp in plot_order:
        if sp == DMEL_SLUG:
            values.append(np.nan)
        else:
            values.append(values_by_species[sp])
    return np.asarray(values, dtype=float)

def add_reference_marker(ax, dmel_position):
    """
    Mark D. melanogaster as a visual reference.
    """
    ax.axvspan(
        dmel_position - 0.48, dmel_position + 0.48, facecolor="#F5F5F5", edgecolor="#555555",
        linewidth=1.0, linestyle="--", zorder=0,
    )

def set_species_axis(ax, plot_order, manifest):
    """
    Add species labels in phylogenetic order.

    Species names are italicized. External targets and D. melanogaster receive an asterisk.
    D. melanogaster is additionally shown in bold as the reference species.
    """
    labels = [short_species_name(sp, manifest) for sp in plot_order]
    x = np.arange(len(plot_order))
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=60, ha='right', fontstyle='italic')
    ax.tick_params(axis='x', pad=13)
    ticklabels = ax.get_xticklabels()
    ticklabels[plot_order.index(DMEL_SLUG)].set_fontweight('bold')
    return x

def add_climate_strip(ax, plot_order, slug_to_climate, y=-0.045, height=0.025):
    """
    Add a narrow climatic-zone strip below the x axis.

    Each species receives one rectangle spanning exactly
    one x-axis category.

    Coordinates along x use data coordinates, while y uses
    axes coordinates. This keeps the strip aligned with the
    species bars/columns independently of the y-axis scale.
    """
    from matplotlib.patches import Rectangle
    for i, sp in enumerate(plot_order):
        zone = slug_to_climate[sp]
        rectangle = Rectangle(
            (i - 0.5, y), width=1.0, height=height, transform=ax.get_xaxis_transform(),
            facecolor=CLIMATE_COLORS[zone], edgecolor="none", clip_on=False, zorder=10,
        )
        ax.add_patch(rectangle)

def left_align_legend(legend):
    """
    Left-align legend title and entries.
    """
    legend._legend_box.align = 'left'
    if legend.get_title() is not None:
        legend.get_title().set_ha('left')

def add_reference_legend(fig, anchor=(0.735, 0.94)):
    legend = fig.legend(
        handles=[reference_legend_handle()], frameon=False, bbox_to_anchor=anchor,
        loc="upper left", borderaxespad=0.0,
    )
    left_align_legend(legend)
    return legend

def add_climate_legend(fig, anchor=(0.735, 0.6)):
    legend = fig.legend(
        handles=climate_legend_handles(), title="Climate zone", frameon=False,
        bbox_to_anchor=anchor, loc="upper left", borderaxespad=0.0,
    )
    left_align_legend(legend)
    return legend

def climate_legend_handles():
    """
    Return legend patches for climatic zones.
    """
    return [
        Patch(facecolor=CLIMATE_COLORS[zone], edgecolor="none", label=CLIMATE_LABELS[zone])
        for zone in ["TROP", "ARID", "TEMP", "BORE"]
    ]

def reference_legend_handle():
    """
    Return legend handle for the D. melanogaster
    reference-species marker.
    """
    return Patch(
        facecolor="#F5F5F5", edgecolor="#555555", linewidth=1.2,
        linestyle="--", label="Reference species",
    )

def finalize_species_plot(fig, right=0.72, bottom=0.3):
    fig.subplots_adjust(right=right, bottom=bottom, top=0.92)

def validate_species(observed, expected, dataset_name):
    """
    Verify that a result table contains every analysed target
    species.

    Additional species are tolerated but reported.
    """
    observed = set(observed)
    expected = set(expected)
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing:
        raise SystemExit(f'ERROR: {dataset_name} is missing species:\n' + '\n'.join(missing))
    if extra:
        print(f'WARNING: {dataset_name} contains additional species that are not plotted:')
        for sp in extra:
            print(f'  {sp}')

def get_label_position(x, y, x_min, x_max, y_min, y_max):
    """
    Determine a deterministic label position from the point's
    relative location within the plotting area.
    """
    x_rel = (x - x_min) / (x_max - x_min)
    y_rel = (y - y_min) / (y_max - y_min)
    if x_rel >= 0.85:
        return {'xytext': (-8, 4), 'ha': 'right', 'va': 'bottom'}
    if x_rel <= 0.15:
        return {'xytext': (8, 4), 'ha': 'left', 'va': 'bottom'}
    if y_rel >= 0.8:
        return {'xytext': (6, -5), 'ha': 'left', 'va': 'top'}
    if y_rel <= 0.2:
        return {'xytext': (6, 5), 'ha': 'left', 'va': 'bottom'}
    return {'xytext': (6, 5), 'ha': 'left', 'va': 'bottom'}
for path in [SUMMARY_FILE, MATRIX_FILE, LONG_FILE, TARGETS_FILE, MANIFEST_FILE, TREE_ORDER_FILE, TRAITS_FILE]:
    require_file(path)
manifest = load_manifest()
tree_order = load_tree_order()
targets = load_targets()
target_order, plot_order = build_species_orders(manifest, tree_order, targets)
dmel_position = plot_order.index(DMEL_SLUG)
slug_to_climate = load_species_traits(manifest)
validate_climate_data(plot_order, slug_to_climate)
print()
print('Species / phylogeny setup')
print('-------------------------')
print(f'Species in combined manifest: {len(manifest)}')
print(f'Species in pruned tree:       {len(tree_order)}')
print(f'Target species analysed:      {len(target_order)}')
print(f'Plot positions incl. Dmel:    {len(plot_order)}')
print(f'Dmel plot position:           {dmel_position + 1}')
print()
summary = pd.read_csv(SUMMARY_FILE, sep='\t')
long_df = pd.read_csv(LONG_FILE, sep='\t')
matrix = pd.read_csv(MATRIX_FILE, sep='\t', index_col=0)
validate_species(summary['species'], target_order, 'species_summary.tsv')
validate_species(long_df['species'], target_order, 'cre_turnover_all_species.tsv')
missing_matrix = [sp for sp in target_order if sp not in matrix.columns]
if missing_matrix:
    raise SystemExit('ERROR: species missing from cre_turnover_matrix.tsv:\n' + '\n'.join(missing_matrix))
summary = summary.set_index('species').loc[target_order].reset_index()
matrix = matrix.loc[:, target_order]
long_df = long_df[long_df['species'].isin(target_order)].copy()
n_reference_cres = matrix.shape[0]
n_target_species = len(target_order)
expected_long_rows = n_reference_cres * n_target_species
if len(long_df) != expected_long_rows:
    raise SystemExit(
        f"ERROR: unexpected long-table size.\nExpected: {expected_long_rows}\nObserved: {len(long_df)}"
    )
matrix_states = set(pd.unique(matrix.values.ravel()))
unknown_states = matrix_states - VALID_STATES
if unknown_states:
    raise SystemExit(
        "ERROR: unknown CRE states in cre_turnover_matrix.tsv:\n" + "\n".join(sorted(unknown_states))
    )
print('CRE-classification input')
print('------------------------')
print(f'Reference CREs: {n_reference_cres}')
print(f'Target species: {n_target_species}')
print(f'Long-table rows: {len(long_df)}')
print()
qc = None
if QC_FILE.exists():
    qc = pd.read_csv(QC_FILE, sep='\t')
    validate_species(qc['species'], target_order, 'species_qc_summary.tsv')
    qc = qc.set_index('species').loc[target_order].reset_index()
else:
    print('WARNING: species_qc_summary.tsv was not found.')
    print('Figures requiring QC metrics will be skipped.')
    print()
# ============================================================
# Species-level figures
# ============================================================

mapping_dict = dict(zip(summary['species'], summary['mapping_rate'] * 100))
mapping_values = insert_reference_placeholder(mapping_dict, plot_order)
fig, ax = plt.subplots(figsize=(14, 6.2))
x = set_species_axis(ax, plot_order, manifest)
ax.bar(x, mapping_values)
add_climate_strip(ax, plot_order, slug_to_climate)
add_reference_marker(ax, dmel_position)
finalize_species_plot(fig)
ref_legend = add_reference_legend(ax, anchor=(AXES_LEGEND_X, AXES_REFERENCE_Y))
ax.add_artist(ref_legend)
add_climate_legend(ax, anchor=(AXES_LEGEND_X, AXES_CLIMATE_Y))
ax.set_ylabel('Reference CREs successfully projected (%)')
ax.set_xlabel('Species in phylogenetic order')
ax.set_ylim(0, 106)
save_figure(fig, 'mapping_rate_phylogenetic')
if qc is not None:
    if 'n_scrmshaw_peaks' not in qc.columns:
        print('WARNING: n_scrmshaw_peaks absent from QC table; skipping SCRMshaw prediction-count figure.')
    else:
        peaks_dict = dict(zip(qc['species'], qc['n_scrmshaw_peaks']))
        peak_values = insert_reference_placeholder(peaks_dict, plot_order)
        fig, ax = plt.subplots(figsize=(14, 6.2))
        x = set_species_axis(ax, plot_order, manifest)
        ax.bar(x, peak_values)
        add_climate_strip(ax, plot_order, slug_to_climate)
        add_reference_marker(ax, dmel_position)
        finalize_species_plot(fig)
        ref_legend = add_reference_legend(ax, anchor=(AXES_LEGEND_X, AXES_REFERENCE_Y))
        ax.add_artist(ref_legend)
        add_climate_legend(ax, anchor=(AXES_LEGEND_X, AXES_CLIMATE_Y))
        ax.set_ylabel('Number of SCRMshaw predictions')
        ax.set_xlabel('Species in phylogenetic order')
        save_figure(fig, 'scrmshaw_peak_count_phylogenetic')
evaluable = summary['positional_match'] + summary['turnover_candidate'] + summary['no_detected_CRE']
if (evaluable == 0).any():
    bad = summary.loc[evaluable == 0, 'species'].tolist()
    raise SystemExit('ERROR: species with zero evaluable CREs:\n' + '\n'.join(bad))
positional_match_pct = summary['positional_match'] / evaluable * 100
turnover_pct = summary['turnover_candidate'] / evaluable * 100
nodetect_pct = summary['no_detected_CRE'] / evaluable * 100
positional_match_dict = dict(zip(summary['species'], positional_match_pct))
turnover_dict = dict(zip(summary['species'], turnover_pct))
nodetect_dict = dict(zip(summary['species'], nodetect_pct))
positional_match_values = insert_reference_placeholder(positional_match_dict, plot_order)
turnover_values = insert_reference_placeholder(turnover_dict, plot_order)
nodetect_values = insert_reference_placeholder(nodetect_dict, plot_order)
p_stack = np.nan_to_num(positional_match_values, nan=0.0)
t_stack = np.nan_to_num(turnover_values, nan=0.0)
n_stack = np.nan_to_num(nodetect_values, nan=0.0)
fig, ax = plt.subplots(figsize=(14, 6.2))
x = set_species_axis(ax, plot_order, manifest)
bars_positional_match = ax.bar(x, p_stack, label=STATE_LABELS['positional_match'], color=STATE_COLORS['positional_match'])
bars_turnover = ax.bar(
    x, t_stack, bottom=p_stack, label=STATE_LABELS["turnover_candidate"],
    color=STATE_COLORS["turnover_candidate"],
)
bars_nodetect = ax.bar(
    x, n_stack, bottom=p_stack + t_stack, label=STATE_LABELS["no_detected_CRE"],
    color=STATE_COLORS["no_detected_CRE"],
)
add_climate_strip(ax, plot_order, slug_to_climate)
for container in [bars_positional_match, bars_turnover, bars_nodetect]:
    container[dmel_position].set_visible(False)
add_reference_marker(ax, dmel_position)
ax.set_ylabel('CRE state among evaluable loci (%)')
ax.set_xlabel('Species in phylogenetic order')
ax.set_ylim(0, 100)
finalize_species_plot(fig)
add_reference_legend(fig, anchor=(FIGURE_LEGEND_X, FIGURE_REFERENCE_Y))
state_legend = fig.legend(
    handles=[
        Patch(facecolor=STATE_COLORS["positional_match"], label=STATE_LABELS["positional_match"]),
        Patch(facecolor=STATE_COLORS["turnover_candidate"], label=STATE_LABELS["turnover_candidate"]),
        Patch(facecolor=STATE_COLORS["no_detected_CRE"], label=STATE_LABELS["no_detected_CRE"]),
    ],
    title="CRE state", frameon=False, bbox_to_anchor=(FIGURE_LEGEND_X, FIGURE_SECONDARY_Y),
    loc="upper left", borderaxespad=0.0,
)
left_align_legend(state_legend)
add_climate_legend(fig, anchor=(FIGURE_LEGEND_X, FIGURE_CLIMATE_Y_LONG))
save_figure(fig, 'cre_state_composition_evaluable')
mapped_pct = summary['mapped'] / summary['n_reference_cres'] * 100
uncertain_pct = summary['uncertain'] / summary['n_reference_cres'] * 100
mapped_dict = dict(zip(summary['species'], mapped_pct))
uncertain_dict = dict(zip(summary['species'], uncertain_pct))
mapped_values = insert_reference_placeholder(mapped_dict, plot_order)
uncertain_values = insert_reference_placeholder(uncertain_dict, plot_order)
mapped_stack = np.nan_to_num(mapped_values, nan=0.0)
uncertain_stack = np.nan_to_num(uncertain_values, nan=0.0)
fig, ax = plt.subplots(figsize=(14, 6.2))
x = set_species_axis(ax, plot_order, manifest)
bars_mapped = ax.bar(x, mapped_stack, label='Evaluable', color=EVALUABLE_COLOR)
bars_uncertain = ax.bar(
    x, uncertain_stack, bottom=mapped_stack, label="Uncertain", color=STATE_COLORS["uncertain"]
)
add_climate_strip(ax, plot_order, slug_to_climate)
bars_mapped[dmel_position].set_visible(False)
bars_uncertain[dmel_position].set_visible(False)
add_reference_marker(ax, dmel_position)
ax.set_ylabel(r'$D.\ melanogaster$ reference CREs (%)')
ax.set_xlabel('Species in phylogenetic order')
ax.set_ylim(0, 100)
finalize_species_plot(fig)
add_reference_legend(fig, anchor=(FIGURE_LEGEND_X, FIGURE_REFERENCE_Y))
mapping_legend = fig.legend(
    handles=[
        Patch(facecolor=EVALUABLE_COLOR, label="Evaluable"),
        Patch(facecolor=STATE_COLORS["uncertain"], label="Uncertain"),
    ],
    title="Mapping status", frameon=False, bbox_to_anchor=(FIGURE_LEGEND_X, FIGURE_SECONDARY_Y),
    loc="upper left", borderaxespad=0.0,
)
left_align_legend(mapping_legend)
add_climate_legend(fig, anchor=(FIGURE_LEGEND_X, FIGURE_CLIMATE_Y_SHORT))
save_figure(fig, 'mapping_evaluability')
# ============================================================
# CRE-level summaries and figures
# ============================================================

cre_summary = pd.DataFrame(index=matrix.index)
cre_summary.index.name = 'dmel_cre_id'
cre_summary['n_positional_match'] = (matrix == 'positional_match').sum(axis=1)
cre_summary['n_turnover_candidate'] = (matrix == 'turnover_candidate').sum(axis=1)
cre_summary['n_no_detected_CRE'] = (matrix == 'no_detected_CRE').sum(axis=1)
cre_summary['n_uncertain'] = (matrix == 'uncertain').sum(axis=1)
cre_summary['n_evaluable'] = n_target_species - cre_summary['n_uncertain']
cre_summary["positional_conservation_rate_evaluable"] = np.where(
    cre_summary["n_evaluable"] > 0, cre_summary["n_positional_match"] / cre_summary["n_evaluable"], np.nan
)
cre_summary["turnover_candidate_rate_evaluable"] = np.where(
    cre_summary["n_evaluable"] > 0,
    cre_summary["n_turnover_candidate"] / cre_summary["n_evaluable"], np.nan,
)
cre_summary.reset_index().to_csv(CRE_SUMMARY_OUT, sep='\t', index=False)
print(f'Wrote: {CRE_SUMMARY_OUT}')
positional_match_counts = Counter(cre_summary['n_positional_match'])
xs = np.arange(0, n_target_species + 1)
ys = np.asarray([positional_match_counts.get(i, 0) for i in xs])
fig, ax = plt.subplots(figsize=(8.5, 5.5))
ax.bar(xs, ys)
ax.set_xlabel('Number of target species with a positional CRE match')
ax.set_ylabel(r'Number of $D.\ melanogaster$ reference CREs')
ax.set_xticks(xs)
fig.tight_layout()
save_figure(fig, 'CRE_positional_conservation_breadth')
rates = cre_summary['positional_conservation_rate_evaluable'].dropna()
fig, ax = plt.subplots(figsize=(8.5, 5.5))
ax.hist(rates, bins=np.linspace(0, 1, 21))
ax.set_xlabel('Fraction of evaluable target species with a positional CRE match')
ax.set_ylabel(r'Number of $D.\ melanogaster$ reference CREs')
ax.set_xlim(0, 1)
fig.tight_layout()
save_figure(fig, 'CRE_conservation_rate_evaluable')
if qc is not None:
    required = {'n_scrmshaw_peaks', 'positional_match_rate_evaluable'}
    if required.issubset(qc.columns):
        qc_plot = qc.copy()
        qc_plot['climatic_zone'] = qc_plot['species'].map(slug_to_climate)
        slug_to_source = dict(zip(manifest['slug'], manifest['source']))
        qc_plot['source'] = qc_plot['species'].map(slug_to_source)
        qc_plot['n_scrmshaw_peaks'] = pd.to_numeric(qc_plot['n_scrmshaw_peaks'], errors='coerce')
        qc_plot['positional_match_rate_pct'] = pd.to_numeric(qc_plot['positional_match_rate_evaluable'], errors='coerce') * 100
        qc_plot = qc_plot.dropna(
            subset=["n_scrmshaw_peaks", "positional_match_rate_pct", "climatic_zone", "source"]
        ).copy()
        unknown_sources = sorted(set(qc_plot['source']) - set(SOURCE_MARKERS))
        if unknown_sources:
            raise SystemExit('ERROR: unknown SCRMshaw source categories:\n' + '\n'.join(unknown_sources))
        qc_plot["is_outlier"] = identify_scatter_outliers(
            qc_plot, x_col="n_scrmshaw_peaks", y_col="positional_match_rate_pct", iqr_factor=1.5
        )
        top_x_species = set(qc_plot.nlargest(2, 'n_scrmshaw_peaks')['species'])
        top_y_species = set(qc_plot.nlargest(2, 'positional_match_rate_pct')['species'])
        qc_plot["label_point"] = (
            qc_plot["is_outlier"] | qc_plot["species"].isin(top_x_species)
            | qc_plot["species"].isin(top_y_species)
        )
        outliers = qc_plot[qc_plot['label_point']].copy()
        x_min = qc_plot['n_scrmshaw_peaks'].min()
        x_max = qc_plot['n_scrmshaw_peaks'].max()
        y_min = qc_plot['positional_match_rate_pct'].min()
        y_max = qc_plot['positional_match_rate_pct'].max()
        fig, ax = plt.subplots(figsize=(9.0, 6.5))
        for zone in ['TROP', 'ARID', 'TEMP', 'BORE']:
            for source in ['generated', 'external']:
                subset = qc_plot[(qc_plot['climatic_zone'] == zone) & (qc_plot['source'] == source)]
                if subset.empty:
                    continue
                ax.scatter(
                    subset["n_scrmshaw_peaks"], subset["positional_match_rate_pct"], color=CLIMATE_COLORS[zone],
                    marker=SOURCE_MARKERS[source], s=65, alpha=0.9, edgecolors="black", linewidths=0.4,
                )
        labelled = outliers.sort_values('positional_match_rate_pct', ascending=False).copy()
        previous_y = None
        for _, row in labelled.iterrows():
            x_val = row['n_scrmshaw_peaks']
            y_val = row['positional_match_rate_pct']
            style = get_label_position(x_val, y_val, x_min, x_max, y_min, y_max)
            if previous_y is not None:
                if abs(y_val - previous_y) < 2.0:
                    dx, dy = style['xytext']
                    style['xytext'] = (dx, dy - 8)
                    style['va'] = 'top'
            ax.annotate(
                short_species_name(row["species"], manifest), (x_val, y_val), xytext=style["xytext"],
                textcoords="offset points", ha=style["ha"], va=style["va"], fontsize=9,
                fontstyle="italic",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.75),
            )
            previous_y = y_val
        climate_handles = [
            Line2D(
                [0], [0], marker="o", linestyle="none", markersize=8,
                markerfacecolor=CLIMATE_COLORS[zone], markeredgecolor="none", label=CLIMATE_LABELS[zone],
            )
            for zone in ["TROP", "ARID", "TEMP", "BORE"]
        ]
        climate_legend = ax.legend(
            handles=climate_handles, title="Climate zone", frameon=False,
            bbox_to_anchor=(AXES_LEGEND_X, SCATTER_CLIMATE_Y), loc="upper left", borderaxespad=0.0,
        )
        left_align_legend(climate_legend)
        ax.add_artist(climate_legend)
        source_handles = [
            Line2D(
                [0], [0], marker=SOURCE_MARKERS[source], linestyle="none", markersize=8,
                markerfacecolor="#666666", markeredgecolor="black", markeredgewidth=0.4,
                label=SOURCE_LABELS[source],
            )
            for source in ["generated", "external"]
        ]
        source_legend = ax.legend(
            handles=source_handles, title="SCRMshaw predictions", frameon=False,
            bbox_to_anchor=(AXES_LEGEND_X, SCATTER_SOURCE_Y), loc="upper left", borderaxespad=0.0,
        )
        left_align_legend(source_legend)
        ax.set_xlabel('Number of SCRMshaw predictions')
        ax.set_ylabel('Positional CRE matches among evaluable loci (%)')
        fig.subplots_adjust(right=0.76, bottom=0.14, top=0.9)
        save_figure(fig, 'scrmshaw_peaks_vs_positional_match_rate')
        print()
        print('Scatterplot outliers labelled:')
        if outliers.empty:
            print('  none')
        else:
            for _, row in outliers.iterrows():
                print(
                    "  "
                    + short_species_name(row["species"], manifest)
                    + " "
                    + f"(peaks={row['n_scrmshaw_peaks']:.0f}, "
                    + f"positional_match={row['positional_match_rate_pct']:.1f}%)"
                )
if qc is not None:
    required = {'turnover_candidate', 'unique_turnover_target_peaks'}
    if required.issubset(qc.columns):
        turnover_pair_dict = dict(zip(qc['species'], qc['turnover_candidate']))
        unique_peak_dict = dict(zip(qc['species'], qc['unique_turnover_target_peaks']))
        pair_values = insert_reference_placeholder(turnover_pair_dict, plot_order)
        unique_values = insert_reference_placeholder(unique_peak_dict, plot_order)
        fig, ax = plt.subplots(figsize=(14, 6.2))
        x = set_species_axis(ax, plot_order, manifest)
        width = 0.36
        pair_bars = ax.bar(
            x - width / 2, np.nan_to_num(pair_values, nan=0.0), width=width,
            label="CRE-level candidates",
        )
        unique_bars = ax.bar(
            x + width / 2, np.nan_to_num(unique_values, nan=0.0), width=width,
            label="Unique target peaks",
        )
        pair_bars[dmel_position].set_visible(False)
        unique_bars[dmel_position].set_visible(False)
        add_climate_strip(ax, plot_order, slug_to_climate)
        add_reference_marker(ax, dmel_position)
        ax.set_ylabel('Count')
        ax.set_xlabel('Species in phylogenetic order')
        finalize_species_plot(fig)
        add_reference_legend(fig, anchor=(FIGURE_LEGEND_X, FIGURE_REFERENCE_Y))
        turnover_legend = fig.legend(
            handles=[pair_bars, unique_bars], title="Turnover metric", frameon=False,
            bbox_to_anchor=(FIGURE_LEGEND_X, FIGURE_SECONDARY_Y), loc="upper left", borderaxespad=0.0,
        )
        left_align_legend(turnover_legend)
        add_climate_legend(fig, anchor=(FIGURE_LEGEND_X, FIGURE_CLIMATE_Y_SHORT))
        save_figure(fig, 'turnover_candidates_vs_unique_peaks')
plot_summary = summary.copy()
if qc is not None:
    qc_cols = [
        col for col in [
            "species", "n_scrmshaw_peaks", "n_prediction_seqids", "n_lifted_seqids", "n_common_seqids",
            "any_positional_overlap", "reciprocal_overlap_ge_50pct", "unique_turnover_target_peaks",
            "qc_flags",
        ]
        if col in qc.columns
    ]
    plot_summary = plot_summary.merge(qc[qc_cols], on='species', how='left', validate='one_to_one')
plot_summary['phylogenetic_target_order'] = np.arange(1, len(plot_summary) + 1)
slug_to_name = dict(zip(manifest['slug'], manifest['species']))
plot_summary.insert(1, 'species_name', plot_summary['species'].map(slug_to_name))
plot_summary.to_csv(PLOT_SUMMARY_OUT, sep='\t', index=False)
print(f'Wrote: {PLOT_SUMMARY_OUT}')
print()
print('Plot generation complete.')
print('-------------------------')
print(f'Target species analysed: {n_target_species}')
print('D. melanogaster included in calculations: NO')
print('D. melanogaster shown as visual reference: YES')
print()
print(f'Figures directory:\n{FIG_DIR}')
print()
print('Species displayed in phylogenetic order:')
for i, sp in enumerate(plot_order, start=1):
    if sp == DMEL_SLUG:
        print(f'{i:2d}\tD. melanogaster* [REFERENCE; no data included]')
    else:
        print(f'{i:2d}\t{short_species_name(sp, manifest)} ({sp})')
