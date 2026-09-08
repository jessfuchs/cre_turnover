#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize final SCRMshaw prediction counts separately "
            "for generated and external species."
        )
    )
    parser.add_argument(
        "--qc",
        required=True,
        type=Path,
        help="species_qc_summary.tsv",
    )
    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help="combined_manifest.tsv containing slug and source",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output TSV with source-specific summary statistics",
    )
    parser.add_argument(
        "--out-species",
        type=Path,
        default=None,
        help="Optional merged per-species QC table",
    )

    args = parser.parse_args()

    qc = pd.read_csv(args.qc, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")

    # ------------------------------------------------------------
    # Required columns
    # ------------------------------------------------------------

    required_qc = {
        "species",
        "n_scrmshaw_peaks",
        "positional_match_rate_evaluable",
    }

    required_manifest = {
        "slug",
        "source",
    }

    missing_qc = required_qc - set(qc.columns)
    missing_manifest = required_manifest - set(manifest.columns)

    if missing_qc:
        raise SystemExit(
            "Missing QC columns: "
            + ", ".join(sorted(missing_qc))
        )

    if missing_manifest:
        raise SystemExit(
            "Missing manifest columns: "
            + ", ".join(sorted(missing_manifest))
        )

    # ------------------------------------------------------------
    # Join prediction source onto species QC
    # ------------------------------------------------------------

    manifest_small = (
        manifest[["slug", "source"]]
        .drop_duplicates()
        .rename(columns={"slug": "species"})
    )

    merged = qc.merge(
        manifest_small,
        on="species",
        how="left",
        validate="one_to_one",
    )

    missing_source = merged.loc[
        merged["source"].isna(),
        "species",
    ].tolist()

    if missing_source:
        raise SystemExit(
            "No prediction source found for:\n"
            + "\n".join(missing_source)
        )

    unexpected_sources = sorted(
        set(merged["source"]) - {"generated", "external"}
    )

    if unexpected_sources:
        raise SystemExit(
            "Unexpected source values: "
            + ", ".join(unexpected_sources)
        )

    # ------------------------------------------------------------
    # Source-specific prediction-depth summary
    # ------------------------------------------------------------

    rows = []

    for source, group in merged.groupby("source", sort=False):
        peaks = group["n_scrmshaw_peaks"]

        rows.append(
            {
                "source": source,
                "n_species": len(group),
                "median_peaks": peaks.median(),
                "q1_peaks": peaks.quantile(0.25),
                "q3_peaks": peaks.quantile(0.75),
                "min_peaks": peaks.min(),
                "max_peaks": peaks.max(),
            }
        )

    summary = pd.DataFrame(rows)

    # Make source ordering deterministic
    summary["source"] = pd.Categorical(
        summary["source"],
        categories=["generated", "external"],
        ordered=True,
    )

    summary = (
        summary
        .sort_values("source")
        .reset_index(drop=True)
    )

    # ------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------

    args.out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.out, sep="\t", index=False)

    if args.out_species is not None:
        args.out_species.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        keep = [
            "species",
            "source",
            "n_scrmshaw_peaks",
            "mapping_rate",
            "positional_match",
            "turnover_candidate",
            "no_detected_CRE",
            "uncertain",
            "positional_match_rate_evaluable",
            "turnover_rate_evaluable",
            "qc_flags",
        ]

        keep = [x for x in keep if x in merged.columns]

        merged[keep].to_csv(
            args.out_species,
            sep="\t",
            index=False,
        )

    # ------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------

    print("\nSCRMshaw prediction depth by source")
    print(summary.to_string(index=False))

    print("\nThesis-ready values:")

    for _, row in summary.iterrows():
        print(
            f"{row['source']}: "
            f"n={int(row['n_species'])}, "
            f"median={row['median_peaks']:g}, "
            f"IQR={row['q1_peaks']:g}-{row['q3_peaks']:g}, "
            f"range={row['min_peaks']:g}-{row['max_peaks']:g}"
        )


if __name__ == "__main__":
    main()