#!/usr/bin/env python3

# ============================================================
# 03 - Build species-level CRE classification QC summary
#
# Purpose:
#   Summarize mapping quality, sequence-ID compatibility,
#   positional overlap, CRE-state counts, and potential
#   technical QC issues for each target species.
#
# Input:
#   - SO_all_species_fbgn.tsv
#   - target_species.txt
#   - species-specific CRE classifications
#   - lifted reference CRE BED files
#
# Output:
#   - species_qc_summary.tsv
#
# Configuration:
#   Paths and QC thresholds are supplied by the pipeline wrapper
#   using config/classification_config.sh.
# ============================================================


import argparse
import csv
from collections import defaultdict
from pathlib import Path


# ============================================================
# Arguments
# ============================================================

parser = argparse.ArgumentParser()

parser.add_argument("--predictions", required=True)
parser.add_argument("--targets", required=True)
parser.add_argument("--turnover-dir", required=True)
parser.add_argument("--lifted-dir", required=True)
parser.add_argument("--out", required=True)

parser.add_argument(
    "--overlap-threshold",
    type=float,
    default=0.50,
)

parser.add_argument(
    "--min-mapping-rate",
    type=float,
    default=0.50,
)

parser.add_argument(
    "--min-scrmshaw-peaks",
    type=int,
    default=50,
)

parser.add_argument(
    "--min-seqid-overlap",
    type=float,
    default=0.25,
)

parser.add_argument(
    "--min-mapped-for-overlap-check",
    type=int,
    default=100,
)

parser.add_argument(
    "--high-mapping-rate",
    type=float,
    default=0.90,
)

args = parser.parse_args()


PREDICTIONS = Path(args.predictions)
TARGETS = Path(args.targets)
TURNOVER_DIR = Path(args.turnover_dir)
LIFTED_DIR = Path(args.lifted_dir)
OUT = Path(args.out)


# ============================================================
# Input checks
# ============================================================

for path in (
    PREDICTIONS,
    TARGETS,
):
    if not path.is_file():
        raise SystemExit(
            f"ERROR: required input file does not exist: {path}"
        )

for path in (
    TURNOVER_DIR,
    LIFTED_DIR,
):
    if not path.is_dir():
        raise SystemExit(
            f"ERROR: required input directory does not exist: {path}"
        )

OUT.parent.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Read target species
# ============================================================

species = [
    line.strip()
    for line in TARGETS.read_text().splitlines()
    if line.strip()
    and not line.lstrip().startswith("#")
]


# ============================================================
# Count SCRMshaw peaks and prediction sequence IDs
# ============================================================

n_peaks = defaultdict(int)
prediction_chroms = defaultdict(set)

with PREDICTIONS.open(
    encoding="utf-8-sig",
    newline=""
) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )

    for row in reader:

        sp = row["species_key"]

        if sp not in species:
            continue

        n_peaks[sp] += 1
        prediction_chroms[sp].add(
            row["chrom"]
        )


# ============================================================
# Build species-level QC summary
# ============================================================

rows_out = []

for sp in species:

    turnover_file = (
        TURNOVER_DIR
        / f"dmel_to_{sp}_cre_turnover.tsv"
    )

    lifted_file = (
        LIFTED_DIR
        / sp
        / f"dmel_reference_cres.{sp}.bed"
    )


    if not turnover_file.is_file():
        raise SystemExit(
            f"ERROR: missing turnover table for {sp}: "
            f"{turnover_file}"
        )

    if not lifted_file.is_file():
        raise SystemExit(
            f"ERROR: missing lifted BED for {sp}: "
            f"{lifted_file}"
        )


    # --------------------------------------------------------
    # Sequence-ID compatibility
    # --------------------------------------------------------

    lifted_chroms = set()

    with lifted_file.open() as handle:

        for line in handle:

            if not line.strip() or line.startswith("#"):
                continue

            lifted_chroms.add(
                line.split("\t", 1)[0]
            )


    pred_chroms = prediction_chroms[sp]

    common_chroms = (
        pred_chroms
        & lifted_chroms
    )

    n_pred_chroms = len(pred_chroms)
    n_lifted_chroms = len(lifted_chroms)
    n_common_chroms = len(common_chroms)

    pred_chrom_overlap_fraction = (
        n_common_chroms / n_pred_chroms
        if n_pred_chroms
        else 0
    )

    lifted_chrom_overlap_fraction = (
        n_common_chroms / n_lifted_chroms
        if n_lifted_chroms
        else 0
    )


    # --------------------------------------------------------
    # Classification results
    # --------------------------------------------------------

    mapped = 0
    unmapped = 0

    any_overlap = 0
    lifted50 = 0
    peak50 = 0
    reciprocal50 = 0

    positional_match = 0
    turnover_candidate = 0
    no_detected = 0
    uncertain = 0

    same_fbgn_anywhere = 0
    same_fbgn_local = 0

    unique_turnover_peaks = set()


    with turnover_file.open(
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:

            if row["alignment_status"] == "mapped":
                mapped += 1

            elif row["alignment_status"] == "unmapped":
                unmapped += 1


            try:
                overlap_bp = float(
                    row["best_overlap_bp"]
                )
            except ValueError:
                overlap_bp = 0


            try:
                frac_lift = float(
                    row["best_overlap_fraction_lifted"]
                )
            except ValueError:
                frac_lift = 0


            try:
                frac_peak = float(
                    row["best_overlap_fraction_peak"]
                )
            except ValueError:
                frac_peak = 0


            if (
                row["alignment_status"] == "mapped"
                and overlap_bp > 0
            ):
                any_overlap += 1


            if (
                row["alignment_status"] == "mapped"
                and frac_lift >= args.overlap_threshold
            ):
                lifted50 += 1


            if (
                row["alignment_status"] == "mapped"
                and frac_peak >= args.overlap_threshold
            ):
                peak50 += 1


            if (
                row["alignment_status"] == "mapped"
                and frac_lift >= args.overlap_threshold
                and frac_peak >= args.overlap_threshold
            ):
                reciprocal50 += 1


            cls = row["class"]

            if cls == "positional_match":
                positional_match += 1

            elif cls == "turnover_candidate":
                turnover_candidate += 1

            elif cls == "no_detected_CRE":
                no_detected += 1

            elif cls == "uncertain":
                uncertain += 1


            if row["same_fbgn_peak_anywhere"] == "yes":
                same_fbgn_anywhere += 1

            if row["same_fbgn_peak_local"] == "yes":
                same_fbgn_local += 1


            if (
                cls == "turnover_candidate"
                and row["local_same_fbgn_peak_id"] != "NA"
            ):
                unique_turnover_peaks.add(
                    row["local_same_fbgn_peak_id"]
                )


    # --------------------------------------------------------
    # Rates
    # --------------------------------------------------------

    n_ref = mapped + unmapped

    mapping_rate = (
        mapped / n_ref
        if n_ref
        else 0
    )

    any_overlap_rate_evaluable = (
        any_overlap / mapped
        if mapped
        else 0
    )

    reciprocal50_rate_evaluable = (
        reciprocal50 / mapped
        if mapped
        else 0
    )

    positional_match_rate_evaluable = (
        positional_match / mapped
        if mapped
        else 0
    )

    turnover_rate_evaluable = (
        turnover_candidate / mapped
        if mapped
        else 0
    )


    # --------------------------------------------------------
    # QC flags
    # --------------------------------------------------------

    flags = []

    if mapping_rate < args.min_mapping_rate:
        flags.append("LOW_MAPPING")

    if n_peaks[sp] < args.min_scrmshaw_peaks:
        flags.append(
            "LOW_SCRMSHAW_PEAK_COUNT"
        )

    if n_common_chroms == 0:
        flags.append("NO_COMMON_SEQIDS")

    elif (
        pred_chrom_overlap_fraction
        < args.min_seqid_overlap
        and
        lifted_chrom_overlap_fraction
        < args.min_seqid_overlap
    ):
        flags.append("LOW_SEQID_OVERLAP")


    if (
        mapped >= args.min_mapped_for_overlap_check
        and any_overlap == 0
    ):
        flags.append(
            "NO_POSITIONAL_OVERLAPS"
        )


    if (
        mapping_rate >= args.high_mapping_rate
        and reciprocal50 == 0
    ):
        flags.append(
            "HIGH_MAPPING_ZERO_PRESENT"
        )


    if (
        n_peaks[sp] < args.min_scrmshaw_peaks
        and mapping_rate >= args.high_mapping_rate
        and reciprocal50 == 0
    ):
        flags.append(
            "LIKELY_PREDICTION_LIMITED"
        )


    if not flags:
        flags = ["OK"]


    rows_out.append({
        "species": sp,
        "n_reference_cres": n_ref,

        "n_scrmshaw_peaks": n_peaks[sp],

        "mapped": mapped,
        "unmapped": unmapped,
        "mapping_rate": mapping_rate,

        "n_prediction_seqids": n_pred_chroms,
        "n_lifted_seqids": n_lifted_chroms,
        "n_common_seqids": n_common_chroms,

        "prediction_seqid_overlap_fraction":
            pred_chrom_overlap_fraction,

        "lifted_seqid_overlap_fraction":
            lifted_chrom_overlap_fraction,

        "any_positional_overlap": any_overlap,

        "any_overlap_rate_evaluable":
            any_overlap_rate_evaluable,

        "lifted_overlap_ge_50pct": lifted50,
        "peak_overlap_ge_50pct": peak50,

        "reciprocal_overlap_ge_50pct":
            reciprocal50,

        "reciprocal_overlap_rate_evaluable":
            reciprocal50_rate_evaluable,

        "positional_match": positional_match,

        "turnover_candidate":
            turnover_candidate,

        "unique_turnover_target_peaks":
            len(unique_turnover_peaks),

        "no_detected_CRE": no_detected,
        "uncertain": uncertain,

        "same_fbgn_peak_anywhere":
            same_fbgn_anywhere,

        "same_fbgn_peak_local":
            same_fbgn_local,

        "positional_match_rate_evaluable":
            positional_match_rate_evaluable,

        "turnover_rate_evaluable":
            turnover_rate_evaluable,

        "qc_flags": ";".join(flags),
    })


# ============================================================
# Write summary
# ============================================================

fields = [
    "species",
    "n_reference_cres",
    "n_scrmshaw_peaks",

    "mapped",
    "unmapped",
    "mapping_rate",

    "n_prediction_seqids",
    "n_lifted_seqids",
    "n_common_seqids",

    "prediction_seqid_overlap_fraction",
    "lifted_seqid_overlap_fraction",

    "any_positional_overlap",
    "any_overlap_rate_evaluable",

    "lifted_overlap_ge_50pct",
    "peak_overlap_ge_50pct",
    "reciprocal_overlap_ge_50pct",
    "reciprocal_overlap_rate_evaluable",

    "positional_match",
    "turnover_candidate",
    "unique_turnover_target_peaks",
    "no_detected_CRE",
    "uncertain",

    "same_fbgn_peak_anywhere",
    "same_fbgn_peak_local",

    "positional_match_rate_evaluable",
    "turnover_rate_evaluable",

    "qc_flags",
]


with OUT.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=fields,
        lineterminator="\n",
    )

    writer.writeheader()

    for row in rows_out:

        output_row = dict(row)

        for key in (
            "mapping_rate",
            "prediction_seqid_overlap_fraction",
            "lifted_seqid_overlap_fraction",
            "any_overlap_rate_evaluable",
            "reciprocal_overlap_rate_evaluable",
            "positional_match_rate_evaluable",
            "turnover_rate_evaluable",
        ):

            output_row[key] = (
                f"{output_row[key]:.6f}"
            )

        writer.writerow(
            output_row
        )


# ============================================================
# Report
# ============================================================

print(f"Wrote: {OUT}")
print()

for row in rows_out:

    print(
        f"{row['species']:20s} "
        f"peaks={row['n_scrmshaw_peaks']:5d} "
        f"mapped={row['mapped']:3d} "
        f"any={row['any_positional_overlap']:3d} "
        f"positional_match={row['positional_match']:3d} "
        f"turnover={row['turnover_candidate']:3d} "
        f"unique_turnover="
        f"{row['unique_turnover_target_peaks']:3d} "
        f"QC={row['qc_flags']}"
    )
