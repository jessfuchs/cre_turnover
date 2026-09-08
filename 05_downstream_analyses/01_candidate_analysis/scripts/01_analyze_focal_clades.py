#!/usr/bin/env python3

# ============================================================
# Analyze focal-clade CRE-state contrasts
#
# Purpose:
#   Identify lineage-specific CRE-state contrasts between one
#   focal species and phylogenetically close comparison species.
#
# Categories:
#   - tier1: positional_match <-> turnover_candidate
#   - tier2: positive CRE state <-> no_detected_CRE
#   - tier3: other focal contrasts with comparison consensus
#   - comparison_mixed: comparison species do not agree
#   - all_same: focal and comparison species share one state
#   - invalid: unexpected CRE-state value
#
# Input/output paths and focal-clade definitions:
#   Supplied by the pipeline wrapper using
#   config/candidate_config.sh.
# ============================================================

from pathlib import Path
import argparse
import pandas as pd

# ============================================================
# CRE-state definitions
# ============================================================

POSITIVE_STATES = {"positional_match", "turnover_candidate"}
VALID_STATES = {"positional_match", "turnover_candidate", "no_detected_CRE", "uncertain"}
OUTPUT_CATEGORIES = ["tier1", "tier2", "tier3", "comparison_mixed", "all_same", "invalid"]

# ============================================================
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze predefined focal-clade CRE-state contrasts."
    )
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--groups", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--traits", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    parser.add_argument("--reference-species", required=True)
    parser.add_argument("--expected-reference-cres", type=int, required=True)
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
            f"ERROR: required columns missing from {label}:\n" + "\n".join(missing)
        )


def tree_name_from_species(species_name):
    return str(species_name).strip().upper().replace(" ", "_")


def load_groups(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(df, {"group_name", "focal_species", "comparison_species"}, "focal-clade table")
    if df["group_name"].duplicated().any():
        raise SystemExit("ERROR: duplicate group_name values in focal-clade table.")

    groups = []
    for _, row in df.iterrows():
        comparisons = [x.strip() for x in row["comparison_species"].split("|") if x.strip()]
        if not row["group_name"].strip() or not row["focal_species"].strip() or not comparisons:
            raise SystemExit("ERROR: incomplete focal-clade definition.")
        if row["focal_species"] in comparisons or len(comparisons) != len(set(comparisons)):
            raise SystemExit(f"ERROR: invalid species definition for {row['group_name']}.")
        groups.append((row["group_name"], row["focal_species"], comparisons))
    return groups


def load_climates(manifest_path, traits_path):
    """
    Map pipeline species slugs to phylogeny names and climatic zones.
    The manifest links slugs to species names; species_traits.tsv
    links phylogeny-compatible tree names to climatic zones.
    """
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
    traits = pd.read_csv(traits_path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(manifest, {"slug", "species", "source"}, "combined manifest")
    require_columns(traits, {"tree_name", "climatic_zone"}, "species traits")

    if manifest["slug"].duplicated().any():
        raise SystemExit("ERROR: duplicate species slugs in combined manifest.")
    if traits["tree_name"].duplicated().any():
        raise SystemExit("ERROR: duplicate tree_name values in species traits.")

    tree_to_climate = dict(zip(traits["tree_name"], traits["climatic_zone"]))
    slug_to_tree = dict(zip(manifest["slug"], manifest["species"].map(tree_name_from_species)))
    slug_to_climate = {
        slug: tree_to_climate.get(tree_name, "")
        for slug, tree_name in slug_to_tree.items()
    }
    return slug_to_tree, slug_to_climate


def classify_focal(row, focal, comparisons):
    """
    Classify one reference CRE within one focal clade.
    Tier 1 requires unanimous comparison species and a strict
    positional_match <-> turnover_candidate contrast with the focal species.
    """
    focal_state = row[focal]
    comparison_states = [row[sp] for sp in comparisons]
    all_states = [focal_state, *comparison_states]

    if not set(all_states).issubset(VALID_STATES):
        return "invalid"
    if len(set(comparison_states)) != 1:
        return "comparison_mixed"

    consensus = comparison_states[0]
    if focal_state == consensus:
        return "all_same"
    if focal_state in POSITIVE_STATES and consensus in POSITIVE_STATES:
        return "tier1"
    if (
        focal_state in POSITIVE_STATES and consensus == "no_detected_CRE"
    ) or (
        focal_state == "no_detected_CRE" and consensus in POSITIVE_STATES
    ):
        return "tier2"
    return "tier3"


def comparison_consensus(row, comparisons):
    values = [row[sp] for sp in comparisons]
    return values[0] if len(set(values)) == 1 else "mixed"

# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()
    for path in [args.matrix, args.groups, args.manifest, args.traits]:
        require_file(path)

    groups = load_groups(args.groups)
    slug_to_tree, slug_to_climate = load_climates(args.manifest, args.traits)

    if args.reference_species not in slug_to_tree:
        raise SystemExit(
            f"ERROR: reference species not found in combined manifest: {args.reference_species}"
        )

    df = pd.read_csv(args.matrix, sep="\t", dtype=str, keep_default_na=False)
    require_columns(df, {"dmel_cre_id"}, "CRE-state matrix")
    if df["dmel_cre_id"].duplicated().any():
        raise SystemExit("ERROR: duplicate dmel_cre_id values in CRE-state matrix.")
    if len(df) != args.expected_reference_cres:
        raise SystemExit(
            "ERROR: unexpected number of reference CREs.\n"
            f"Expected: {args.expected_reference_cres}\nObserved: {len(df)}"
        )

    # --------------------------------------------------------
    # Validate focal-clade species and annotations
    # --------------------------------------------------------

    required_species = sorted({sp for _, focal, comps in groups for sp in [focal, *comps]})
    missing = sorted(set(required_species) - set(df.columns))
    if missing:
        raise SystemExit("ERROR: focal-clade species missing from CRE-state matrix:\n" + "\n".join(missing))
    missing_manifest = sorted(set(required_species) - set(slug_to_tree))
    if missing_manifest:
        raise SystemExit("ERROR: focal-clade species missing from combined manifest:\n" + "\n".join(missing_manifest))
    missing_climate = sorted(sp for sp in required_species if not slug_to_climate.get(sp, ""))
    if missing_climate:
        raise SystemExit("ERROR: climatic-zone data missing for:\n" + "\n".join(missing_climate))
    if args.reference_species in required_species:
        raise SystemExit("ERROR: reference species must not be included in focal-clade target species.")

    observed = {x for x in pd.unique(df[required_species].values.ravel()) if x != ""}
    unknown = sorted(observed - VALID_STATES)
    if unknown:
        raise SystemExit("ERROR: unexpected CRE-state values:\n" + "\n".join(unknown))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.summary_out.parent.mkdir(parents=True, exist_ok=True)
    summary_rows = []

    # ========================================================
    # Analyze each predefined focal clade independently
    # ========================================================

    for group_name, focal, comparisons in groups:
        species = [focal, *comparisons]
        comparison_climates = {slug_to_climate[sp] for sp in comparisons}
        if len(comparison_climates) != 1:
            raise SystemExit(
                f"ERROR: comparison species do not share one climatic zone in {group_name}."
            )

        sub = df[["dmel_cre_id", *species]].copy()
        sub["category"] = sub.apply(lambda row: classify_focal(row, focal, comparisons), axis=1)
        sub["group_name"] = group_name
        sub["focal_species"] = focal
        sub["focal_climate"] = slug_to_climate[focal]
        sub["focal_state"] = sub[focal]
        sub["comparison_species"] = "|".join(comparisons)
        sub["comparison_climate"] = next(iter(comparison_climates))
        sub["comparison_consensus_state"] = sub.apply(
            lambda row: comparison_consensus(row, comparisons), axis=1
        )

        front = [
            "group_name", "dmel_cre_id", "category", "focal_species", "focal_climate",
            "focal_state", "comparison_species", "comparison_climate",
            "comparison_consensus_state", *species,
        ]
        sub = sub[front + [c for c in sub.columns if c not in front]]
        sub.to_csv(args.out_dir / f"{group_name}_all.tsv", sep="\t", index=False)

        for category in OUTPUT_CATEGORIES:
            sub.loc[sub["category"] == category].to_csv(
                args.out_dir / f"{group_name}_{category}.tsv", sep="\t", index=False
            )

        # ----------------------------------------------------
        # Group-level summary
        # ----------------------------------------------------

        counts = sub["category"].value_counts()
        row = {
            "group_name": group_name,
            "focal_species": focal,
            "focal_climate": slug_to_climate[focal],
            "comparison_species": "|".join(comparisons),
            "comparison_climate": next(iter(comparison_climates)),
            "n_comparison_species": len(comparisons),
            "n_reference_cres": len(sub),
            **{f"n_{cat}": int(counts.get(cat, 0)) for cat in OUTPUT_CATEGORIES},
        }
        summary_rows.append(row)

        print()
        print("=" * 72)
        print(group_name)
        print("=" * 72)
        print(f"Focal species: {focal} ({slug_to_climate[focal]})")
        print(f"Comparison species: {', '.join(comparisons)}")
        print(f"Comparison climate: {next(iter(comparison_climates))}")
        print(f"Tier 1: {row['n_tier1']}")
        print(f"Tier 2: {row['n_tier2']}")
        print(f"Tier 3: {row['n_tier3']}")
        print(f"Comparison mixed: {row['n_comparison_mixed']}")
        print(f"All same: {row['n_all_same']}")
        print(f"Invalid: {row['n_invalid']}")

    
    # ========================================================
    # Cross-clade summary
    # ========================================================

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(args.summary_out, sep="\t", index=False)

    print()
    print("=" * 72)
    print("Focal-clade candidate analysis complete")
    print("=" * 72)
    print(f"Reference CREs: {len(df)}")
    print(f"Focal clades analyzed: {len(groups)}")
    print(f"Total Tier 1 rows: {summary['n_tier1'].sum()}")
    print(f"Total Tier 2 rows: {summary['n_tier2'].sum()}")
    print(f"Wrote focal-clade tables to:\n{args.out_dir}")
    print(f"Wrote clade summary:\n{args.summary_out}")


if __name__ == "__main__":
    main()
