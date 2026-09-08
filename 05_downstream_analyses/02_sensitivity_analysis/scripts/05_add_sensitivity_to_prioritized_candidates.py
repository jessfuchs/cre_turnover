#!/usr/bin/env python3

# ============================================================
# Add sensitivity evidence to prioritized Tier-1 candidates
#
# Purpose:
#   Merge candidate-level sensitivity results with the primary
#   prioritized Tier-1 CRE table.
#
# Added annotations:
#   - candidate retention across sensitivity scenarios
#   - robustness class
#   - candidate-priority stability
#   - scenarios in which candidates are lost or reprioritized
#
# The original candidate priority remains the primary ranking;
# sensitivity robustness is added as supporting evidence.
#
# Input/output paths and robustness threshold are supplied by
# config/sensitivity_config.sh via the pipeline wrapper.
# ============================================================

from datetime import datetime
from pathlib import Path
import argparse
import math

import pandas as pd

VALID_RETENTION_VALUES = {"yes", "no"}


def parse_args():
    p = argparse.ArgumentParser(
        description="Separate Tier-1 candidate sensitivity to overlap and local-distance thresholds."
    )
    p.add_argument("--retention-matrix", type=Path, required=True)
    p.add_argument("--global-summary", type=Path, required=True)
    p.add_argument("--scenario-manifest", type=Path, required=True)
    p.add_argument("--primary-scenario", required=True)
    p.add_argument("--out-candidates", type=Path, required=True)
    p.add_argument("--out-summary", type=Path, required=True)
    p.add_argument("--out-classes", type=Path, required=True)
    p.add_argument("--metadata-out", type=Path, required=True)
    return p.parse_args()


def require_file(path):
    if not path.is_file():
        raise SystemExit(f"ERROR: required file not found:\n{path}")


def require_columns(df, required, label):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(f"ERROR: {label} missing columns:\n" + "\n".join(missing))


def normalize_retention(series):
    values = series.astype(str).str.strip().str.lower()
    unknown = sorted(set(values) - VALID_RETENTION_VALUES)
    if unknown:
        raise SystemExit("ERROR: unexpected retention values:\n" + "\n".join(unknown))
    return values.eq("yes")


def classify_parameter_sensitivity(overlap_retained, overlap_n, distance_retained, distance_n):
    overlap_sensitive = overlap_retained < overlap_n
    distance_sensitive = distance_retained < distance_n
    if not overlap_sensitive and not distance_sensitive:
        return "robust_to_both"
    if overlap_sensitive and not distance_sensitive:
        return "overlap_sensitive_only"
    if not overlap_sensitive and distance_sensitive:
        return "distance_sensitive_only"
    return "sensitive_to_both"


def scenario_list(manifest, primary_overlap, primary_distance, dimension):
    if dimension == "overlap":
        mask = manifest["local_gene_distance_bp"].eq(primary_distance)
        return manifest.loc[mask].sort_values("reciprocal_overlap")["scenario"].tolist()
    mask = manifest["reciprocal_overlap"].apply(lambda x: math.isclose(x, primary_overlap, rel_tol=0, abs_tol=1e-12))
    return manifest.loc[mask].sort_values("local_gene_distance_bp")["scenario"].tolist()


def main():
    args = parse_args()
    for path in [args.retention_matrix, args.global_summary, args.scenario_manifest]:
        require_file(path)

    manifest = pd.read_csv(args.scenario_manifest, sep="\t", dtype=str, keep_default_na=False)
    require_columns(manifest, {"scenario", "reciprocal_overlap", "local_gene_distance_bp"}, "scenario manifest")
    if manifest["scenario"].duplicated().any():
        raise SystemExit("ERROR: duplicate scenarios in scenario manifest.")
    manifest["reciprocal_overlap"] = pd.to_numeric(manifest["reciprocal_overlap"], errors="raise")
    manifest["local_gene_distance_bp"] = pd.to_numeric(manifest["local_gene_distance_bp"], errors="raise").astype(int)

    primary = manifest.loc[manifest["scenario"].eq(args.primary_scenario)]
    if len(primary) != 1:
        raise SystemExit("ERROR: primary scenario must occur exactly once in scenario manifest.")
    primary_overlap = float(primary.iloc[0]["reciprocal_overlap"])
    primary_distance = int(primary.iloc[0]["local_gene_distance_bp"])

    overlap_scenarios = scenario_list(manifest, primary_overlap, primary_distance, "overlap")
    distance_scenarios = scenario_list(manifest, primary_overlap, primary_distance, "distance")
    if args.primary_scenario not in overlap_scenarios or args.primary_scenario not in distance_scenarios:
        raise SystemExit("ERROR: primary scenario is missing from one parameter axis.")
    if len(overlap_scenarios) < 2 or len(distance_scenarios) < 2:
        raise SystemExit("ERROR: each parameter-specific analysis requires the primary plus at least one alternative setting.")

    retention = pd.read_csv(args.retention_matrix, sep="\t", dtype=str, keep_default_na=False)
    required = {"dmel_cre_id", *overlap_scenarios, *distance_scenarios}
    require_columns(retention, required, "candidate retention matrix")
    if retention.empty or retention["dmel_cre_id"].duplicated().any():
        raise SystemExit("ERROR: retention matrix is empty or contains duplicate CRE IDs.")

    needed = list(dict.fromkeys([*overlap_scenarios, *distance_scenarios]))
    bool_matrix = retention[["dmel_cre_id", *needed]].copy()
    for scenario in needed:
        bool_matrix[scenario] = normalize_retention(bool_matrix[scenario])
    if not bool_matrix[args.primary_scenario].all():
        failed = bool_matrix.loc[~bool_matrix[args.primary_scenario], "dmel_cre_id"].tolist()
        raise SystemExit(
            "ERROR: primary candidates are not retained in the primary scenario:\n" + "\n".join(failed)
        )

    overlap_nonprimary = [s for s in overlap_scenarios if s != args.primary_scenario]
    distance_nonprimary = [s for s in distance_scenarios if s != args.primary_scenario]

    rows = []
    for _, row in bool_matrix.iterrows():
        n_overlap = int(row[overlap_scenarios].sum())
        n_distance = int(row[distance_scenarios].sum())
        rows.append({
            "dmel_cre_id": row["dmel_cre_id"],
            "overlap_n_settings": len(overlap_scenarios),
            "overlap_retained_n": n_overlap,
            "overlap_retention_pct": 100.0 * n_overlap / len(overlap_scenarios),
            "overlap_sensitive": "yes" if n_overlap < len(overlap_scenarios) else "no",
            "overlap_lost_scenarios": "|".join(s for s in overlap_nonprimary if not row[s]) or "NA",
            "distance_n_settings": len(distance_scenarios),
            "distance_retained_n": n_distance,
            "distance_retention_pct": 100.0 * n_distance / len(distance_scenarios),
            "distance_sensitive": "yes" if n_distance < len(distance_scenarios) else "no",
            "distance_lost_scenarios": "|".join(s for s in distance_nonprimary if not row[s]) or "NA",
            "parameter_sensitivity_class": classify_parameter_sensitivity(
                n_overlap, len(overlap_scenarios), n_distance, len(distance_scenarios)
            ),
        })

    candidate = pd.DataFrame(rows)
    class_order = {
        "robust_to_both": 0,
        "overlap_sensitive_only": 1,
        "distance_sensitive_only": 2,
        "sensitive_to_both": 3,
    }
    candidate["_order"] = candidate["parameter_sensitivity_class"].map(class_order)
    candidate = candidate.sort_values(["_order", "dmel_cre_id"], kind="mergesort").drop(columns="_order").reset_index(drop=True)
    candidate["overlap_retention_pct"] = candidate["overlap_retention_pct"].round(1)
    candidate["distance_retention_pct"] = candidate["distance_retention_pct"].round(1)

    global_summary = pd.read_csv(args.global_summary, sep="\t", dtype=str, keep_default_na=False)
    require_columns(global_summary, {"scenario", "n_changed"}, "global sensitivity summary")
    global_summary["n_changed"] = pd.to_numeric(global_summary["n_changed"], errors="raise").astype(int)
    if global_summary["scenario"].duplicated().any():
        raise SystemExit("ERROR: duplicate scenarios in global sensitivity summary.")
    global_by_scenario = dict(zip(global_summary["scenario"], global_summary["n_changed"]))
    missing_global = sorted(set(needed) - set(global_by_scenario))
    if missing_global:
        raise SystemExit("ERROR: global sensitivity summary missing scenarios:\n" + "\n".join(missing_global))
    if global_by_scenario[args.primary_scenario] != 0:
        raise SystemExit("ERROR: primary scenario has non-zero global state changes.")

    def dimension_summary(name, scenarios, nonprimary, retained_col, pct_col, fixed_parameter):
        retained = candidate[retained_col]
        losses = {s: int((~bool_matrix[s]).sum()) for s in nonprimary}
        global_changes = {s: int(global_by_scenario[s]) for s in nonprimary}
        return {
            "parameter_dimension": name,
            "fixed_parameter": fixed_parameter,
            "scenarios": "|".join(scenarios),
            "n_settings": len(scenarios),
            "n_candidates": len(candidate),
            "n_retained_all_settings": int(retained.eq(len(scenarios)).sum()),
            "pct_retained_all_settings": round(100.0 * retained.eq(len(scenarios)).mean(), 1),
            "mean_retention_pct": round(float(candidate[pct_col].mean()), 1),
            "median_retention_pct": round(float(candidate[pct_col].median()), 1),
            "n_sensitive_to_parameter": int(retained.lt(len(scenarios)).sum()),
            "pct_sensitive_to_parameter": round(100.0 * retained.lt(len(scenarios)).mean(), 1),
            "candidate_losses_nonprimary": "|".join(f"{s}:{losses[s]}" for s in nonprimary) or "NA",
            "global_state_changes_nonprimary": "|".join(f"{s}:{global_changes[s]}" for s in nonprimary) or "NA",
        }

    parameter_summary = pd.DataFrame([
        dimension_summary(
            "reciprocal_overlap",
            overlap_scenarios,
            overlap_nonprimary,
            "overlap_retained_n",
            "overlap_retention_pct",
            f"local_gene_distance_bp={primary_distance}",
        ),
        dimension_summary(
            "local_gene_distance",
            distance_scenarios,
            distance_nonprimary,
            "distance_retained_n",
            "distance_retention_pct",
            f"reciprocal_overlap={primary_overlap:g}",
        ),
    ])

    class_summary = candidate["parameter_sensitivity_class"].value_counts().rename_axis(
        "parameter_sensitivity_class"
    ).reset_index(name="n_candidates")
    class_summary = pd.DataFrame({"parameter_sensitivity_class": list(class_order)}).merge(
        class_summary, on="parameter_sensitivity_class", how="left"
    )
    class_summary["n_candidates"] = class_summary["n_candidates"].fillna(0).astype(int)
    class_summary["percent_candidates"] = (100.0 * class_summary["n_candidates"] / len(candidate)).round(1)

    for path in [args.out_candidates, args.out_summary, args.out_classes, args.metadata_out]:
        path.parent.mkdir(parents=True, exist_ok=True)
    candidate.to_csv(args.out_candidates, sep="\t", index=False)
    parameter_summary.to_csv(args.out_summary, sep="\t", index=False)
    class_summary.to_csv(args.out_classes, sep="\t", index=False)

    metadata = pd.DataFrame([{
        "script": Path(__file__).name,
        "run_timestamp": datetime.now().astimezone().isoformat(),
        "primary_scenario": args.primary_scenario,
        "primary_reciprocal_overlap": primary_overlap,
        "primary_local_gene_distance_bp": primary_distance,
        "n_overlap_settings": len(overlap_scenarios),
        "n_distance_settings": len(distance_scenarios),
        "n_candidates": len(candidate),
    }])
    metadata.to_csv(args.metadata_out, sep="\t", index=False)

    print("\n" + "=" * 72)
    print("PARAMETER-SPECIFIC TIER-1 SENSITIVITY")
    print("=" * 72)
    print(f"Primary scenario: {args.primary_scenario}")
    print(f"Overlap settings: {len(overlap_scenarios)}")
    print(f"Distance settings: {len(distance_scenarios)}")
    print(class_summary.to_string(index=False))
    print(f"Wrote: {args.out_candidates}")


if __name__ == "__main__":
    main()
