#!/usr/bin/env python3

# ============================================================
# Build gene-level Tier-1 candidate summary
#
# Purpose:
#   Aggregate CRE x FBgn candidate evidence into one row per
#   unique D. melanogaster gene.
#
# Summary:
#   Gene-level support is described by:
#   - number of associated Tier-1 CREs
#   - primary and secondary candidate-gene assignments
#   - original SCRMshaw target-gene support
#   - recurrent and focal CRE support
#   - strongest CRE-level candidate priority
#   - Tier-1 clade support
#   - carried-forward gene-assignment QC
#
# Notes:
#   Gene-assignment QC flags are retained as descriptive
#   evidence and are not used as automatic exclusion criteria.
#
#   No weighted gene-level score is calculated. Candidate genes
#   are instead ordered by transparent evidence variables.
#
# Input/output paths:
#   Supplied by the pipeline wrapper using
#   config/candidate_config.sh.
# ============================================================

import argparse
from datetime import datetime
import hashlib
from pathlib import Path
import platform
import sys

import pandas as pd


REQUIRED_COLUMNS = {
    "dmel_cre_id", "fbgn", "gene_role", "assignment_evidence",
    "candidate_priority", "downstream_priority",
    "n_total_tier1_clades", "n_focal_tier1_clades",
    "n_secondary_only_clades", "tier1_clades",
    "is_reference_target", "is_primary_candidate",
    "is_secondary_candidate", "gene_assignment_qc",
}

CANDIDATE_PRIORITY_ORDER = {
    "focal_recurrent": 1,
    "focal_plus_secondary_recurrent": 2,
    "secondary_recurrent": 3,
    "focal_single": 4,
    "secondary_single": 5,
}

DOWNSTREAM_PRIORITY_ORDER = {"high": 1, "medium": 2, "exploratory": 3}

GENE_ROLE_ORDER = {"primary_candidate": 1, "secondary_candidate": 2, "reference_target_only": 3}

RECURRENT_PRIORITIES = {
    "focal_recurrent", "focal_plus_secondary_recurrent", "secondary_recurrent"
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a gene-level summary from the exploded Tier-1 CRE x FBgn table."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--metadata-out", type=Path, required=True)
    return parser.parse_args()


def require_file(path):
    if not path.is_file():
        raise SystemExit(f"ERROR: required file not found:\n{path}")


def require_columns(df, required, label):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(
            f"ERROR: required columns missing from {label}:\n" + "\n".join(missing)
        )


def normalize_missing(value):
    if pd.isna(value):
        return "NA"
    value = str(value).strip()
    return "NA" if value.lower() in {"", "na", "nan", "none"} else value


def join_unique(values, delimiter="|"):
    clean = {
        normalize_missing(value)
        for value in values
        if normalize_missing(value) != "NA"
    }
    return delimiter.join(sorted(clean))


def split_pipe_values(values):
    result = set()
    for value in values:
        value = normalize_missing(value)
        if value != "NA":
            result.update(x.strip() for x in value.split("|") if x.strip())
    return sorted(result)


def best_value(values, ranking, label):
    clean = [normalize_missing(value) for value in values if normalize_missing(value) != "NA"]
    unknown = sorted(set(clean) - set(ranking))
    if unknown:
        raise SystemExit(f"ERROR: unknown {label} values:\n" + "\n".join(unknown))
    return min(clean, key=ranking.get) if clean else "NA"


def main():
    args = parse_args()
    require_file(args.input)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_out.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    require_columns(df, REQUIRED_COLUMNS, "exploded CRE x FBgn table")
    if df.empty:
        raise SystemExit("ERROR: exploded CRE x FBgn table is empty.")
    if df.duplicated(subset=["dmel_cre_id", "fbgn"]).any():
        bad = df.loc[
            df.duplicated(subset=["dmel_cre_id", "fbgn"], keep=False),
            ["dmel_cre_id", "fbgn"],
        ].drop_duplicates()
        raise SystemExit("ERROR: duplicate CRE x FBgn pairs:\n\n" + bad.to_string(index=False))

    for column in ["n_total_tier1_clades", "n_focal_tier1_clades", "n_secondary_only_clades"]:
        df[column] = pd.to_numeric(df[column], errors="raise").astype(int)

    df["is_recurrent_cre"] = df["candidate_priority"].isin(RECURRENT_PRIORITIES)
    df["has_focal_support"] = df["n_focal_tier1_clades"] >= 1
    df["is_high_priority_cre"] = df["downstream_priority"] == "high"
    df["is_medium_priority_cre"] = df["downstream_priority"] == "medium"
    df["is_exploratory_cre"] = df["downstream_priority"] == "exploratory"

    gene_rows = []
    for fbgn, sub in df.groupby("fbgn", sort=True):
        cre_ids = sorted(sub["dmel_cre_id"].unique())
        ids = {
            "primary": sorted(sub.loc[sub["is_primary_candidate"] == "yes", "dmel_cre_id"].unique()),
            "secondary": sorted(sub.loc[sub["is_secondary_candidate"] == "yes", "dmel_cre_id"].unique()),
            "reference": sorted(sub.loc[sub["is_reference_target"] == "yes", "dmel_cre_id"].unique()),
            "recurrent": sorted(sub.loc[sub["is_recurrent_cre"], "dmel_cre_id"].unique()),
            "focal": sorted(sub.loc[sub["has_focal_support"], "dmel_cre_id"].unique()),
            "high": sorted(sub.loc[sub["is_high_priority_cre"], "dmel_cre_id"].unique()),
            "medium": sorted(sub.loc[sub["is_medium_priority_cre"], "dmel_cre_id"].unique()),
            "exploratory": sorted(sub.loc[sub["is_exploratory_cre"], "dmel_cre_id"].unique()),
            "qc_pass": sorted(sub.loc[sub["gene_assignment_qc"] == "PASS", "dmel_cre_id"].unique()),
            "qc_flagged": sorted(sub.loc[sub["gene_assignment_qc"] != "PASS", "dmel_cre_id"].unique()),
        }
        clades = split_pipe_values(sub["tier1_clades"])

        gene_rows.append({
            "fbgn": fbgn,
            "best_downstream_priority": best_value(sub["downstream_priority"], DOWNSTREAM_PRIORITY_ORDER, "downstream_priority"),
            "best_candidate_priority": best_value(sub["candidate_priority"], CANDIDATE_PRIORITY_ORDER, "candidate_priority"),
            "best_gene_role": best_value(sub["gene_role"], GENE_ROLE_ORDER, "gene_role"),
            "n_candidate_cres": len(cre_ids),
            "candidate_cre_ids": "|".join(cre_ids),
            "n_primary_cres": len(ids["primary"]),
            "primary_cre_ids": "|".join(ids["primary"]),
            "n_secondary_cres": len(ids["secondary"]),
            "secondary_cre_ids": "|".join(ids["secondary"]),
            "n_reference_target_cres": len(ids["reference"]),
            "reference_target_cre_ids": "|".join(ids["reference"]),
            "n_recurrent_cres": len(ids["recurrent"]),
            "recurrent_cre_ids": "|".join(ids["recurrent"]),
            "n_focal_supported_cres": len(ids["focal"]),
            "focal_supported_cre_ids": "|".join(ids["focal"]),
            "n_high_priority_cres": len(ids["high"]),
            "high_priority_cre_ids": "|".join(ids["high"]),
            "n_medium_priority_cres": len(ids["medium"]),
            "medium_priority_cre_ids": "|".join(ids["medium"]),
            "n_exploratory_cres": len(ids["exploratory"]),
            "exploratory_cre_ids": "|".join(ids["exploratory"]),
            "n_unique_tier1_clades": len(clades),
            "tier1_clades": "|".join(clades),
            "max_total_tier1_clades_per_cre": int(sub["n_total_tier1_clades"].max()),
            "max_focal_tier1_clades_per_cre": int(sub["n_focal_tier1_clades"].max()),
            "max_secondary_only_clades_per_cre": int(sub["n_secondary_only_clades"].max()),
            "n_qc_pass_cres": len(ids["qc_pass"]),
            "qc_pass_cre_ids": "|".join(ids["qc_pass"]),
            "n_qc_flagged_cres": len(ids["qc_flagged"]),
            "qc_flagged_cre_ids": "|".join(ids["qc_flagged"]),
            "gene_assignment_qc_values": join_unique(sub["gene_assignment_qc"]),
            "assignment_evidence": join_unique(sub["assignment_evidence"]),
            "has_primary_assignment": "yes" if ids["primary"] else "no",
            "has_secondary_assignment": "yes" if ids["secondary"] else "no",
            "has_reference_target_support": "yes" if ids["reference"] else "no",
        })

    out = pd.DataFrame(gene_rows)
    if out.empty:
        raise SystemExit("ERROR: no gene-level rows generated.")
    if out["fbgn"].duplicated().any():
        raise SystemExit("ERROR: duplicate FBgn rows in gene-level summary.")

    out["_downstream_rank"] = out["best_downstream_priority"].map(DOWNSTREAM_PRIORITY_ORDER).fillna(99)
    out["_candidate_rank"] = out["best_candidate_priority"].map(CANDIDATE_PRIORITY_ORDER).fillna(99)
    out = out.sort_values(
        [
            "_downstream_rank", "n_high_priority_cres", "n_recurrent_cres",
            "n_focal_supported_cres", "n_primary_cres", "n_candidate_cres",
            "n_unique_tier1_clades", "_candidate_rank", "fbgn",
        ],
        ascending=[True, False, False, False, False, False, False, True, True],
        kind="mergesort",
    ).drop(columns=["_downstream_rank", "_candidate_rank"]).reset_index(drop=True)

    out.to_csv(args.out, sep="\t", index=False)

    metadata = pd.DataFrame([{
        "script": Path(__file__).name,
        "run_timestamp": datetime.now().astimezone().isoformat(),
        "n_cre_fbgn_input_rows": len(df),
        "n_unique_candidate_cres": df["dmel_cre_id"].nunique(),
        "n_unique_candidate_genes": len(out),
        "n_genes_with_primary_support": int((out["has_primary_assignment"] == "yes").sum()),
        "n_genes_with_reference_target_support": int((out["has_reference_target_support"] == "yes").sum()),
        "n_genes_with_flagged_cre_evidence": int((out["n_qc_flagged_cres"] > 0).sum()),
    }])
    metadata.to_csv(args.metadata_out, sep="\t", index=False)

    print()
    print("=" * 72)
    print("Gene-level candidate summary complete")
    print("=" * 72)
    print(f"Input CRE x FBgn rows: {len(df)}")
    print(f"Unique candidate CREs: {df['dmel_cre_id'].nunique()}")
    print(f"Unique candidate genes: {len(out)}")
    print()
    print("Best downstream priority:")
    print(out["best_downstream_priority"].value_counts().to_string())
    print()
    print("Top gene-level candidates:")
    print(out[[
        "fbgn", "best_downstream_priority", "best_candidate_priority", "best_gene_role",
        "n_candidate_cres", "n_primary_cres", "n_recurrent_cres",
        "n_focal_supported_cres", "n_high_priority_cres", "n_unique_tier1_clades",
        "n_qc_flagged_cres",
    ]].head(20).to_string(index=False))
    print()
    print(f"Wrote gene-level summary:\n{args.out}")
    print(f"Wrote metadata:\n{args.metadata_out}")


if __name__ == "__main__":
    main()
