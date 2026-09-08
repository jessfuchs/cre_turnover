#!/usr/bin/env python3

# ============================================================
# Test focal Tier-1 enrichment
#
# Purpose:
#   Test whether the predefined focal species is enriched as
#   the discordant species among Tier-1 singleton contrasts.
#
# Tier-1 singleton:
#   - all species are positional_match or turnover_candidate
#   - exactly one species differs from all remaining species
#
# Null model:
#   Each species in the clade is equally likely to carry the
#   singleton difference, giving an expected focal share of
#   1 / number of species.
#
# A one-sided exact binomial test evaluates focal enrichment.
# P-values across focal clades are adjusted using Benjamini-
# Hochberg correction.
#
# The test is exploratory because CRE events are not guaranteed
# to repositional_match independent evolutionary events.
#
# Input/output paths and focal-clade definitions are supplied
# by the climate-analysis configuration via the wrapper.
# ============================================================

from pathlib import Path
import argparse

import numpy as np
import pandas as pd
from scipy.stats import binomtest


POSITIVE_STATES = {"positional_match", "turnover_candidate"}

CLADE_LABELS = {
    "rufa_group": "Rufa group",
    "immigrans_group": "Immigrans group",
    "obscura_group": "Obscura group",
    "azteca_affinis_miranda_group": "Azteca group",
    "teissieri_group": "Teissieri group",
    "repleta_group": "Repleta group",
}


# ============================================================
# Repository-relative default paths
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
CLIMATE_DIR = SCRIPT_DIR.parent
DOWNSTREAM_DIR = CLIMATE_DIR.parent
PROJECT_ROOT = DOWNSTREAM_DIR.parent
RESULTS_DIR = CLIMATE_DIR / "results"

DEFAULT_MATRIX = PROJECT_ROOT / "04_cre_classification" / "results" / "cre_turnover_matrix.tsv"
DEFAULT_GROUPS = DOWNSTREAM_DIR / "01_candidate_analysis" / "config" / "focal_clades.tsv"
DEFAULT_OUT_CANDIDATES = RESULTS_DIR / "focal_tier1_candidates.tsv"
DEFAULT_OUT_ENRICHMENT = RESULTS_DIR / "focal_tier1_enrichment.tsv"


# ============================================================
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Test enrichment of focal Tier-1 singleton CRE contrasts."
    )
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--groups", type=Path, default=DEFAULT_GROUPS)
    parser.add_argument("--out-candidates", type=Path, default=DEFAULT_OUT_CANDIDATES)
    parser.add_argument("--out-enrichment", type=Path, default=DEFAULT_OUT_ENRICHMENT)
    return parser.parse_args()


# ============================================================
# Helpers
# ============================================================

def require_file(path):
    if not path.is_file():
        raise SystemExit(f"ERROR: required file not found:\n{path}")


def require_columns(df, required, label):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(
            f"ERROR: {label} missing columns:\n" + "\n".join(missing)
        )


def load_groups(path):
    """Load the shared focal-clade definitions."""

    groups = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    require_columns(
        groups,
        {"group_name", "focal_species", "comparison_species"},
        "focal-clade table",
    )

    if groups["group_name"].duplicated().any():
        raise SystemExit("ERROR: duplicate group names in focal-clade table.")

    result = {}

    for row in groups.itertuples(index=False):
        comparisons = [
            x.strip()
            for x in row.comparison_species.split("|")
            if x.strip()
        ]

        if not row.focal_species or not comparisons:
            raise SystemExit(
                f"ERROR: incomplete focal-clade definition for {row.group_name}."
            )

        result[row.group_name] = {
            "focal": row.focal_species,
            "species": [row.focal_species, *comparisons],
        }

    return result


def get_singleton_tier1(row):
    """
    Return discordant species and state direction for a strict
    positional_match <-> turnover_candidate singleton contrast.
    """

    if not all(state in POSITIVE_STATES for state in row.values):
        return None

    counts = row.value_counts()

    if len(counts) != 2:
        return None

    singleton = counts[counts == 1]
    consensus = counts[counts == len(row) - 1]

    if len(singleton) != 1 or len(consensus) != 1:
        return None

    discordant_state = singleton.index[0]
    consensus_state = consensus.index[0]
    discordant_species = row.index[row == discordant_state][0]

    return discordant_species, discordant_state, consensus_state


def benjamini_hochberg(p_values):
    """Benjamini-Hochberg adjustment preserving original order."""

    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    ranked = p_values[order]

    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    result = np.empty(len(adjusted))
    result[order] = adjusted

    return result


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    require_file(args.matrix)
    require_file(args.groups)

    args.out_candidates.parent.mkdir(parents=True, exist_ok=True)
    args.out_enrichment.parent.mkdir(parents=True, exist_ok=True)

    groups = load_groups(args.groups)

    matrix = pd.read_csv(
        args.matrix,
        sep="\t",
        index_col=0,
        dtype=str,
        keep_default_na=False,
    )

    if matrix.index.duplicated().any():
        raise SystemExit("ERROR: duplicate CRE IDs in CRE-state matrix.")

    candidate_rows = []
    enrichment_rows = []

    # ========================================================
    # Analyze each focal clade
    # ========================================================

    for clade_name, config in groups.items():
        species = config["species"]
        focal = config["focal"]

        missing_species = [sp for sp in species if sp not in matrix.columns]

        if missing_species:
            raise SystemExit(
                f"ERROR: species missing for {clade_name}:\n"
                + "\n".join(missing_species)
            )

        events = []

        for cre_id, row in matrix[species].iterrows():
            result = get_singleton_tier1(row)

            if result is None:
                continue

            discordant_species, discordant_state, consensus_state = result
            is_focal = discordant_species == focal

            focal_direction = ""

            if is_focal:
                if (
                    discordant_state == "turnover_candidate"
                    and consensus_state == "positional_match"
                ):
                    focal_direction = "focal_turnover"

                elif (
                    discordant_state == "positional_match"
                    and consensus_state == "turnover_candidate"
                ):
                    focal_direction = "focal_positional_match"

            event = {
                "dmel_cre_id": cre_id,
                "discordant_species": discordant_species,
                "discordant_state": discordant_state,
                "consensus_state": consensus_state,
                "is_focal_tier1": is_focal,
                "focal_direction": focal_direction,
            }

            events.append(event)

            # The output retains all singleton events because both
            # focal and non-focal events define the enrichment test.
            candidate_rows.append({
                "clade": clade_name,
                "clade_display": CLADE_LABELS.get(
                    clade_name,
                    clade_name.replace("_", " ").title(),
                ),
                "dmel_cre_id": cre_id,
                "focal_species": focal,
                **event,
            })

        events_df = pd.DataFrame(events)

        n_all = len(events_df)

        if n_all:
            focal_events = events_df.loc[events_df["is_focal_tier1"]]
            n_focal = len(focal_events)

            n_focal_turnover = int(
                (focal_events["focal_direction"] == "focal_turnover").sum()
            )
            n_focal_positional_match = int(
                (focal_events["focal_direction"] == "focal_positional_match").sum()
            )

        else:
            n_focal = 0
            n_focal_turnover = 0
            n_focal_positional_match = 0

        n_species = len(species)
        expected_share = 1 / n_species

        observed_share = n_focal / n_all if n_all else np.nan
        enrichment_ratio = observed_share / expected_share if n_all else np.nan

        if n_all:
            test = binomtest(
                k=n_focal,
                n=n_all,
                p=expected_share,
                alternative="greater",
            )

            ci = test.proportion_ci(confidence_level=0.95, method="exact")
            p_value = test.pvalue
            ci_low = ci.low
            ci_high = ci.high

        else:
            p_value = np.nan
            ci_low = np.nan
            ci_high = np.nan

        enrichment_rows.append({
            "clade": clade_name,
            "clade_display": CLADE_LABELS.get(
                clade_name,
                clade_name.replace("_", " ").title(),
            ),
            "focal_species": focal,
            "n_species": n_species,
            "n_tier1_singletons": n_all,
            "n_focal_tier1": n_focal,
            "n_focal_turnover": n_focal_turnover,
            "n_focal_positional_match": n_focal_positional_match,
            "observed_focal_share": observed_share,
            "expected_focal_share": expected_share,
            "enrichment_ratio": enrichment_ratio,
            "ci95_low": ci_low,
            "ci95_high": ci_high,
            "p_value": p_value,
        })

    # ========================================================
    # Multiple-testing correction and output
    # ========================================================

    candidate_df = pd.DataFrame(candidate_rows)
    enrichment_df = pd.DataFrame(enrichment_rows)

    valid = enrichment_df["p_value"].notna()

    enrichment_df.loc[valid, "p_adj_BH"] = benjamini_hochberg(
        enrichment_df.loc[valid, "p_value"]
    )

    candidate_df.to_csv(args.out_candidates, sep="\t", index=False)
    enrichment_df.to_csv(args.out_enrichment, sep="\t", index=False)

    # ========================================================
    # Console summary
    # ========================================================

    print()
    print("=" * 72)
    print("FOCAL TIER-1 ENRICHMENT")
    print("=" * 72)

    columns = [
        "clade_display",
        "n_tier1_singletons",
        "n_focal_tier1",
        "n_focal_turnover",
        "n_focal_positional_match",
        "observed_focal_share",
        "expected_focal_share",
        "enrichment_ratio",
        "p_value",
        "p_adj_BH",
    ]

    print()
    print(enrichment_df[columns].to_string(index=False))

    print()
    print(
        "NOTE: exact binomial enrichment is exploratory because "
        "CRE events may not repositional_match independent evolutionary events."
    )

    print()
    print(f"Wrote enrichment:\n{args.out_enrichment}")
    print(f"Wrote singleton events:\n{args.out_candidates}")


if __name__ == "__main__":
    main()
