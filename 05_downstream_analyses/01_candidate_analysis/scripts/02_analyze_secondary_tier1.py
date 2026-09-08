#!/usr/bin/env python3

# ============================================================
# Analyze secondary singleton CRE-state contrasts
#
# Purpose:
#   Identify CREs for which exactly one species within a focal
#   clade differs from a common state shared by all remaining
#   species.
#
# Categories:
#   - tier1: positional_match <-> turnover_candidate
#   - tier2: positive CRE state <-> no_detected_CRE
#   - tier3: other singleton contrasts
#   - other_pattern: no valid singleton contrast
#
# Secondary Tier-1 candidates are additionally summarized
# across clades to identify recurrent CREs.
#
# Input/output paths and clade definitions:
#   Supplied by the pipeline wrapper using
#   config/candidate_config.sh.
# ============================================================

from pathlib import Path
from datetime import datetime
import argparse
import pandas as pd

# ============================================================
# CRE-state definitions
# ============================================================

VALID_STATES = {"positional_match", "turnover_candidate", "no_detected_CRE", "uncertain"}
POSITIVE_STATES = {"positional_match", "turnover_candidate"}
OUTPUT_TIERS = ["tier1", "tier2", "tier3", "other_pattern"]

# ============================================================
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify singleton CRE-state contrasts within predefined focal clades."
    )
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--groups", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--out-long", type=Path, required=True)
    parser.add_argument("--out-summary", type=Path, required=True)
    parser.add_argument("--out-group-summary", type=Path, required=True)
    parser.add_argument("--metadata-out", type=Path, required=True)
    parser.add_argument("--expected-reference-cres", type=int, required=True)
    parser.add_argument("--recurrence-min-clades", type=int, required=True)
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


def load_groups(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(df, {"group_name", "focal_species", "comparison_species"}, "focal-clade table")
    if df["group_name"].duplicated().any():
        raise SystemExit("ERROR: duplicate group_name values in focal-clade table.")

    groups = []
    for _, row in df.iterrows():
        comparisons = [x.strip() for x in row["comparison_species"].split("|") if x.strip()]
        species = [row["focal_species"].strip(), *comparisons]
        if not species[0] or not comparisons or len(species) != len(set(species)):
            raise SystemExit(f"ERROR: invalid species definition for {row['group_name']}.")
        groups.append((row["group_name"], species))
    return groups


def join_unique(values):
    clean = {
        str(value).strip()
        for value in values
        if str(value).strip().lower() not in {"", "na", "nan", "none"}
    }
    return "|".join(sorted(clean))


def empty_pattern(pattern_class="not_singleton_contrast"):
    return {
        "pattern_class": pattern_class,
        "tier": "other_pattern",
        "discordant_species": "NA",
        "discordant_state": "NA",
        "consensus_state": "NA",
        "n_consensus_species": 0,
    }


def classify_secondary_pattern(row, species):
    """
    Classify one CRE across one phylogenetic clade.
    
    A secondary singleton pattern requires exactly one species
    to differ from all other species in the clade.
    
    Examples
    --------
    3 species:
        A != B = B
    
    4 species:
        A != B = B = B
    
    5 species:
        A != B = B = B = B
    """
    states = {sp: row[sp] for sp in species}
    if set(states.values()) - VALID_STATES:
        return empty_pattern("invalid_state")

    counts = pd.Series(list(states.values())).value_counts()
    if len(counts) != 2:
        return empty_pattern()

    singleton_states = [state for state, n in counts.items() if n == 1]
    if len(singleton_states) != 1:
        return empty_pattern()

    discordant_state = singleton_states[0]
    discordant = [sp for sp in species if states[sp] == discordant_state]
    if len(discordant) != 1:
        return empty_pattern()

    discordant_species = discordant[0]
    consensus_states = {states[sp] for sp in species if sp != discordant_species}
    if len(consensus_states) != 1:
        return empty_pattern()

    consensus_state = next(iter(consensus_states))
    n_consensus = len(species) - 1

    # --------------------------------------------------------
    # Tier 1:
    # positional_match <-> turnover_candidate
    # --------------------------------------------------------
    if (
        discordant_state in POSITIVE_STATES
        and consensus_state in POSITIVE_STATES
        and discordant_state != consensus_state
    ):
        pattern_class, tier = "positional_match_vs_turnover", "tier1"

    # --------------------------------------------------------
    # Tier 2:
    # positive <-> no_detected_CRE
    # --------------------------------------------------------
    elif (
        discordant_state == "no_detected_CRE" and consensus_state in POSITIVE_STATES
    ) or (
        consensus_state == "no_detected_CRE" and discordant_state in POSITIVE_STATES
    ):
        pattern_class, tier = "detection_contrast", "tier2"

    # --------------------------------------------------------
    # Tier 3:
    # singleton contrast exists, but involves uncertain
    # or another non-priority combination
    # --------------------------------------------------------
    else:
        pattern_class, tier = "other_singleton_contrast", "tier3"

    return {
        "pattern_class": pattern_class,
        "tier": tier,
        "discordant_species": discordant_species,
        "discordant_state": discordant_state,
        "consensus_state": consensus_state,
        "n_consensus_species": n_consensus,
    }

# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()
    require_file(args.matrix)
    require_file(args.groups)
    if args.recurrence_min_clades < 1:
        raise SystemExit("ERROR: recurrence-min-clades must be >= 1.")

    groups = load_groups(args.groups)
    df = pd.read_csv(args.matrix, sep="\t", dtype=str, keep_default_na=False)
    require_columns(df, {"dmel_cre_id"}, "CRE-state matrix")
    if df["dmel_cre_id"].duplicated().any():
        raise SystemExit("ERROR: duplicate dmel_cre_id values in CRE-state matrix.")
    if len(df) != args.expected_reference_cres:
        raise SystemExit(
            "ERROR: unexpected number of reference CREs.\n"
            f"Expected: {args.expected_reference_cres}\nObserved: {len(df)}"
        )

    # ========================================================
    # Validate required species
    # ========================================================

    required_species = sorted({sp for _, species in groups for sp in species})
    missing = sorted(set(required_species) - set(df.columns))
    if missing:
        raise SystemExit("ERROR: required species missing from CRE-state matrix:\n" + "\n".join(missing))
    observed = {x for x in pd.unique(df[required_species].values.ravel()) if x != ""}
    unknown = sorted(observed - VALID_STATES)
    if unknown:
        raise SystemExit("ERROR: unexpected CRE-state values:\n" + "\n".join(unknown))

    for path in [args.outdir, args.out_long.parent, args.out_summary.parent,
                 args.out_group_summary.parent, args.metadata_out.parent]:
        path.mkdir(parents=True, exist_ok=True)

    # ========================================================
    # Analyze each clade
    # ========================================================

    group_summary_rows = []
    secondary_tier1_rows = []

    for group_name, species in groups:
        sub = df[["dmel_cre_id", *species]].copy()
        classified = sub.apply(
            lambda row: classify_secondary_pattern(row, species),
            axis=1,
            result_type="expand",
        )
        sub = pd.concat([sub, classified], axis=1)
        sub["group_name"] = group_name
        sub["group_species"] = "|".join(species)
        sub["n_group_species"] = len(species)
        sub = sub.sort_values(["tier", "dmel_cre_id"], kind="mergesort").reset_index(drop=True)

        sub.to_csv(args.outdir / f"{group_name}_all_patterns.tsv", sep="\t", index=False)
        for tier in OUTPUT_TIERS:
            sub.loc[sub["tier"] == tier].to_csv(
                args.outdir / f"{group_name}_{tier}.tsv", sep="\t", index=False
            )

        tier1 = sub.loc[sub["tier"] == "tier1"]
        for _, row in tier1.iterrows():
            secondary_tier1_rows.append({
                "dmel_cre_id": row["dmel_cre_id"],
                "group_name": group_name,
                "discordant_species": row["discordant_species"],
                "discordant_state": row["discordant_state"],
                "consensus_state": row["consensus_state"],
                "n_group_species": len(species),
                "n_consensus_species": row["n_consensus_species"],
                "group_species": "|".join(species),
                "pattern_class": row["pattern_class"],
            })

        # ----------------------------------------------------
        # Group-level summary
        # ----------------------------------------------------

        counts = sub["tier"].value_counts()
        group_summary_rows.append({
            "group_name": group_name,
            "group_species": "|".join(species),
            "n_group_species": len(species),
            "n_reference_cres": len(sub),
            "n_tier1": int(counts.get("tier1", 0)),
            "n_tier2": int(counts.get("tier2", 0)),
            "n_tier3": int(counts.get("tier3", 0)),
            "n_other_pattern": int(counts.get("other_pattern", 0)),
        })

        print()
        print("=" * 72)
        print(group_name)
        print("=" * 72)
        print(f"Species: {', '.join(species)}")
        print(f"Tier 1: {counts.get('tier1', 0)}")
        print(f"Tier 2: {counts.get('tier2', 0)}")
        print(f"Tier 3: {counts.get('tier3', 0)}")
        print(f"Other patterns: {counts.get('other_pattern', 0)}")

    # ========================================================
    # Combined Secondary Tier-1 long table
    # ========================================================

    secondary_long = pd.DataFrame(secondary_tier1_rows)
    if not secondary_long.empty:
        secondary_long = secondary_long.sort_values(
            ["dmel_cre_id", "group_name", "discordant_species"], kind="mergesort"
        ).reset_index(drop=True)
        if secondary_long.duplicated(subset=["dmel_cre_id", "group_name"]).any():
            raise SystemExit("ERROR: duplicate CRE x clade rows in secondary Tier-1 table.")
    secondary_long.to_csv(args.out_long, sep="\t", index=False)

    if secondary_long.empty:
        secondary_summary = pd.DataFrame(columns=[
            "dmel_cre_id", "n_secondary_tier1_clades", "secondary_tier1_clades",
            "discordant_species", "discordant_states", "consensus_states",
            "secondary_recurrence",
        ])
    else:
        secondary_summary = secondary_long.groupby("dmel_cre_id", sort=True).agg(
            n_secondary_tier1_clades=("group_name", "nunique"),
            secondary_tier1_clades=("group_name", join_unique),
            discordant_species=("discordant_species", join_unique),
            discordant_states=("discordant_state", join_unique),
            consensus_states=("consensus_state", join_unique),
        ).reset_index()
        secondary_summary["secondary_recurrence"] = secondary_summary[
            "n_secondary_tier1_clades"
        ].apply(lambda n: "recurrent" if int(n) >= args.recurrence_min_clades else "single_clade")
        secondary_summary = secondary_summary.sort_values(
            ["n_secondary_tier1_clades", "dmel_cre_id"],
            ascending=[False, True], kind="mergesort"
        ).reset_index(drop=True)

    secondary_summary.to_csv(args.out_summary, sep="\t", index=False)
    pd.DataFrame(group_summary_rows).sort_values("group_name", kind="mergesort").to_csv(
        args.out_group_summary, sep="\t", index=False
    )

    # ========================================================
    # Run metadata
    # ========================================================

    metadata = pd.DataFrame([{
        "script": Path(__file__).name,
        "run_timestamp": datetime.now().astimezone().isoformat(),
        "n_reference_cres": len(df),
        "n_clades": len(groups),
        "recurrence_min_clades": args.recurrence_min_clades,
        "n_secondary_tier1_rows": len(secondary_long),
        "n_unique_secondary_tier1_cres": secondary_long["dmel_cre_id"].nunique() if not secondary_long.empty else 0,
        "n_recurrent_secondary_cres": int((secondary_summary["secondary_recurrence"] == "recurrent").sum()) if not secondary_summary.empty else 0,
    }])
    metadata.to_csv(args.metadata_out, sep="\t", index=False)

    # ========================================================
    # Final console summary
    # ========================================================

    print()
    print("=" * 72)
    print("Secondary Tier-1 analysis complete")
    print("=" * 72)
    print(f"Reference CREs: {len(df)}")
    print(f"Clades analyzed: {len(groups)}")
    print(f"Secondary Tier-1 group x CRE rows: {len(secondary_long)}")
    print(f"Wrote combined Tier-1 table:\n{args.out_long}")
    print(f"Wrote recurrent CRE summary:\n{args.out_summary}")


if __name__ == "__main__":
    main()
