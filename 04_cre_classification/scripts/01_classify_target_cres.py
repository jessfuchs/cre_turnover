#!/usr/bin/env python3

#!/usr/bin/env python3

# ============================================================
# 01 - Classify D. melanogaster reference CREs
#
# Purpose:
#   Classify each reference CRE in one target species based on
#   positional conservation and local same-FBgn CRE predictions.
#
# Input:
#   - dmel_reference_cres.tsv
#   - lifted reference CRE BED
#   - SO_all_species_fbgn.tsv
#
# Output:
#   - species-specific CRE classification TSV
#
# Configuration:
#   Parameters and paths are supplied by the pipeline wrapper
#   using config/classification_config.sh.
# ============================================================

import argparse
import csv
from pathlib import Path
from collections import defaultdict


# ============================================================
# Helper functions
# ============================================================

def split_fbgn(value):
    """
    Parse a | separated FBgn field.
    """
    if value is None:
        return set()

    value = value.strip()

    if not value or value in {"NA", ".", "None"}:
        return set()

    return {
        x.strip()
        for x in value.split("|")
        if x.strip().startswith("FBgn")
    }


def parse_distance(value):
    """
    Parse a CRE-to-gene distance.

    Absolute distance is used because the local-neighborhood
    criterion concerns genomic distance, not direction.
    """
    if value is None:
        return None

    value = value.strip()

    if not value or value in {"NA", ".", "None"}:
        return None

    try:
        return abs(float(value))
    except ValueError:
        return None


def fbgn_distance_map(row):
    """
    Build:
        FBgn -> minimum CRE-to-associated-gene distance

    for one target-species SCRMshaw prediction.

    The flanking FBgn field is paired with
    distance_flanking_gene.

    The next-flanking FBgn field is paired with
    distance_next_gene.
    """

    result = {}

    associations = [
        (
            row.get("dmel_fbgn_flanking_gene", ""),
            row.get("distance_flanking_gene", ""),
        ),
        (
            row.get("dmel_fbgn_next_flanking_gene", ""),
            row.get("distance_next_gene", ""),
        ),
    ]

    for fbgn_field, distance_field in associations:

        fbgns = split_fbgn(fbgn_field)
        distance = parse_distance(distance_field)

        if distance is None:
            continue

        for fbgn in fbgns:
            if fbgn not in result:
                result[fbgn] = distance
            else:
                result[fbgn] = min(
                    result[fbgn],
                    distance,
                )

    return result


# ============================================================
# Arguments
# ============================================================

parser = argparse.ArgumentParser(
    description=(
        "Classify Dmel reference CREs in a target species "
        "using positional homology, SCRMshaw predictions, "
        "and Dmel FBgn target-gene associations."
    )
)

parser.add_argument(
    "--species",
    required=True,
)

parser.add_argument(
    "--reference",
    required=True,
    help="dmel_reference_cres.tsv",
)

parser.add_argument(
    "--lifted",
    required=True,
    help="Lifted Dmel CRE BED for the target species",
)

parser.add_argument(
    "--predictions",
    required=True,
    help="SO_all_species_fbgn.tsv",
)

parser.add_argument(
    "--out",
    required=True,
)

parser.add_argument(
    "--min-lifted-overlap",
    type=float,
    default=0.50,
    help=(
        "Minimum reciprocal overlap fraction required for "
        "a positional CRE match. Applied to both the lifted "
        "CRE and the target prediction. Default: 0.50"
    ),
)

parser.add_argument(
    "--max-local-gene-distance",
    type=float,
    default=24000,
    help=(
        "Maximum CRE-to-associated-gene distance in bp for "
        "a same-FBgn prediction to count as a local turnover "
        "candidate. Default: 24000 bp."
    ),
)

args = parser.parse_args()


# ============================================================
# 1. Load Dmel reference CREs
# ============================================================

reference = {}

with open(args.reference) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        cre_id = row["dmel_cre_id"]

        if cre_id in reference:
            raise SystemExit(
                f"Duplicate reference CRE ID: {cre_id}"
            )

        reference[cre_id] = {
            **row,
            "fbgn_set": split_fbgn(
                row.get("fbgn_target_genes", "")
            ),
        }


# ============================================================
# 2. Load lifted homologous CRE intervals
#
# BED6:
# chrom start end dmel_cre_id score strand
# ============================================================

lifted = {}

with open(args.lifted) as f:

    for line in f:

        if not line.strip() or line.startswith("#"):
            continue

        cols = line.rstrip("\n").split("\t")

        if len(cols) < 4:
            raise SystemExit(
                f"Invalid lifted BED line:\n{line}"
            )

        cre_id = cols[3]

        if cre_id in lifted:
            raise SystemExit(
                f"Multiple lifted mappings found for {cre_id}. "
                "This workflow expects one primary mapping."
            )

        lifted[cre_id] = {
            "chrom": cols[0],
            "start": int(cols[1]),
            "end": int(cols[2]),
            "strand": (
                cols[5]
                if len(cols) >= 6
                else "."
            ),
        }


# ============================================================
# 3. Load target-species SCRMshaw predictions
# ============================================================

pred_by_chrom = defaultdict(list)

# FBgn -> list of target predictions associated with that FBgn
pred_by_fbgn = defaultdict(list)

all_target_fbgn = set()

species_peak_number = 0

with open(args.predictions) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for source_row_number, row in enumerate(reader, start=2):

        if row["species_key"] != args.species:
            continue

        species_peak_number += 1

        try:
            start = int(row["start0"])
            end = int(row["end0"])
        except ValueError:
            raise SystemExit(
                f"Invalid target coordinates at input line "
                f"{source_row_number}"
            )

        if end <= start:
            raise SystemExit(
                f"Invalid target interval at input line "
                f"{source_row_number}: "
                f"{row['chrom']}:{start}-{end}"
            )

        distances = fbgn_distance_map(row)
        fbgn_set = set(distances)

        # Retain FBgn IDs even if distance information is missing,
        # because they still count for the genome-wide QC field.
        all_fbgn = (
            split_fbgn(
                row.get(
                    "dmel_fbgn_flanking_gene",
                    "",
                )
            )
            |
            split_fbgn(
                row.get(
                    "dmel_fbgn_next_flanking_gene",
                    "",
                )
            )
        )

        all_target_fbgn |= all_fbgn

        pred = {
            "id": (
                f"{args.species}_PEAK_"
                f"{species_peak_number:05d}"
            ),
            "chrom": row["chrom"],
            "start": start,
            "end": end,
            "score": row.get(
                "scrmshaw_score",
                "NA",
            ),
            "amplitude": row.get(
                "peak_amplitude",
                "NA",
            ),
            "fbgn": all_fbgn,
            "fbgn_distance": distances,
            "training_set": row.get(
                "training_set",
                "NA",
            ),
            "method": row.get(
                "method",
                "NA",
            ),
            "rank": row.get(
                "rank",
                "NA",
            ),
        }

        pred_by_chrom[pred["chrom"]].append(pred)

        for fbgn in all_fbgn:
            pred_by_fbgn[fbgn].append(pred)


# ============================================================
# 4. Classification
# ============================================================

output_rows = []


for cre_id, ref in reference.items():

    ref_genes = ref["fbgn_set"]

    # --------------------------------------------------------
    # A. Same-FBgn prediction anywhere in target genome
    # --------------------------------------------------------

    same_fbgn_anywhere = bool(
        ref_genes & all_target_fbgn
    )


    # --------------------------------------------------------
    # B. Find best LOCAL same-FBgn prediction
    #
    # Local means:
    # same Dmel FBgn AND CRE-to-associated-gene distance
    # <= empirical threshold (default 24 kb).
    # --------------------------------------------------------

    local_candidates = []

    for fbgn in ref_genes:

        for pred in pred_by_fbgn.get(fbgn, []):

            distance = pred["fbgn_distance"].get(fbgn)

            # Same-FBgn information may exist without a usable
            # distance. Such cases count as "anywhere" but not
            # as local turnover candidates.
            if distance is None:
                continue

            if distance <= args.max_local_gene_distance:

                local_candidates.append({
                    "pred": pred,
                    "fbgn": fbgn,
                    "distance": distance,
                })


    # Remove duplicate peak/FBgn combinations if necessary.
    unique_local = {}

    for candidate in local_candidates:

        key = (
            candidate["pred"]["id"],
            candidate["fbgn"],
        )

        if key not in unique_local:
            unique_local[key] = candidate
        else:
            if (
                candidate["distance"]
                < unique_local[key]["distance"]
            ):
                unique_local[key] = candidate

    local_candidates = list(
        unique_local.values()
    )

    same_fbgn_local = bool(local_candidates)


    # Choose the nearest local same-FBgn candidate for reporting.
    best_local = None

    if local_candidates:
        best_local = min(
            local_candidates,
            key=lambda x: (
                x["distance"],
                x["pred"]["id"],
                x["fbgn"],
            ),
        )


    # --------------------------------------------------------
    # C. Unmapped homologous interval -> uncertain
    # --------------------------------------------------------

    if cre_id not in lifted:

        output_rows.append({
            "dmel_cre_id": cre_id,
            "species": args.species,
            "alignment_status": "unmapped",

            "target_chrom": "NA",
            "target_start0": "NA",
            "target_end0": "NA",
            "lifted_length": "NA",

            "best_peak_id": "NA",
            "best_overlap_bp": 0,
            "best_overlap_fraction_lifted": 0,
            "best_overlap_fraction_peak": 0,

            "positional_peak": "no",

            "same_fbgn_peak_anywhere": (
                "yes"
                if same_fbgn_anywhere
                else "no"
            ),

            "same_fbgn_peak_local": (
                "yes"
                if same_fbgn_local
                else "no"
            ),

            "local_same_fbgn_peak_id": (
                best_local["pred"]["id"]
                if best_local
                else "NA"
            ),

            "local_same_fbgn": (
                best_local["fbgn"]
                if best_local
                else "NA"
            ),

            "best_local_same_fbgn_distance_bp": (
                f"{best_local['distance']:.1f}"
                if best_local
                else "NA"
            ),

            "shared_fbgn_with_best_peak": "NA",

            "dmel_target_fbgn": (
                "|".join(sorted(ref_genes))
                if ref_genes
                else "NA"
            ),

            "gene_support_at_best_peak": "NA",

            "class": "uncertain",
        })

        continue


    # --------------------------------------------------------
    # D. Find best positional overlap
    # --------------------------------------------------------

    loc = lifted[cre_id]

    lifted_len = (
        loc["end"] - loc["start"]
    )

    if lifted_len <= 0:
        raise SystemExit(
            f"Invalid lifted interval for {cre_id}"
        )

    best = None

    for pred in pred_by_chrom.get(
        loc["chrom"],
        [],
    ):

        overlap = max(
            0,
            min(
                loc["end"],
                pred["end"],
            )
            -
            max(
                loc["start"],
                pred["start"],
            ),
        )

        if overlap <= 0:
            continue

        pred_len = (
            pred["end"] - pred["start"]
        )

        frac_lifted = (
            overlap / lifted_len
        )

        frac_peak = (
            overlap / pred_len
        )

        candidate = {
            "pred": pred,
            "overlap": overlap,
            "frac_lifted": frac_lifted,
            "frac_peak": frac_peak,
        }

        # Best positional candidate is selected primarily
        # by fraction of the lifted CRE covered, then by
        # fraction of target peak covered, then bp overlap.
        if best is None:
            best = candidate
        else:
            old_key = (
                best["frac_lifted"],
                best["frac_peak"],
                best["overlap"],
            )

            new_key = (
                candidate["frac_lifted"],
                candidate["frac_peak"],
                candidate["overlap"],
            )

            if new_key > old_key:
                best = candidate


    # --------------------------------------------------------
    # E. Reciprocal 50/50 positional criterion
    # --------------------------------------------------------

    positional = (
        best is not None
        and
        best["frac_lifted"]
        >= args.min_lifted_overlap
        and
        best["frac_peak"]
        >= args.min_lifted_overlap
    )


    # --------------------------------------------------------
    # F. Biological classification
    # --------------------------------------------------------

    if positional:

        cls = "present"

    elif same_fbgn_local:

        cls = "turnover_candidate"

    else:

        cls = "no_detected_CRE"


    # --------------------------------------------------------
    # G. Information about best positional prediction
    # --------------------------------------------------------

    if best is not None:

        shared = (
            ref_genes
            & best["pred"]["fbgn"]
        )

        best_peak_id = (
            best["pred"]["id"]
        )

        overlap_bp = (
            best["overlap"]
        )

        frac_lifted = (
            best["frac_lifted"]
        )

        frac_peak = (
            best["frac_peak"]
        )

        shared_fbgn = (
            "|".join(sorted(shared))
            if shared
            else "NA"
        )

        if not ref_genes:

            gene_support = (
                "no_fbgn_information"
            )

        elif shared:

            gene_support = (
                "shared_fbgn"
            )

        else:

            gene_support = (
                "discordant_fbgn"
            )

    else:

        best_peak_id = "NA"
        overlap_bp = 0
        frac_lifted = 0
        frac_peak = 0
        shared_fbgn = "NA"

        if not ref_genes:
            gene_support = (
                "no_fbgn_information"
            )
        else:
            gene_support = (
                "no_peak_at_locus"
            )


    # --------------------------------------------------------
    # H. Output
    # --------------------------------------------------------

    output_rows.append({
        "dmel_cre_id": cre_id,
        "species": args.species,
        "alignment_status": "mapped",

        "target_chrom": loc["chrom"],
        "target_start0": loc["start"],
        "target_end0": loc["end"],
        "lifted_length": lifted_len,

        "best_peak_id": best_peak_id,
        "best_overlap_bp": overlap_bp,

        "best_overlap_fraction_lifted": (
            f"{frac_lifted:.6f}"
        ),

        "best_overlap_fraction_peak": (
            f"{frac_peak:.6f}"
        ),

        "positional_peak": (
            "yes"
            if positional
            else "no"
        ),

        "same_fbgn_peak_anywhere": (
            "yes"
            if same_fbgn_anywhere
            else "no"
        ),

        "same_fbgn_peak_local": (
            "yes"
            if same_fbgn_local
            else "no"
        ),

        "local_same_fbgn_peak_id": (
            best_local["pred"]["id"]
            if best_local
            else "NA"
        ),

        "local_same_fbgn": (
            best_local["fbgn"]
            if best_local
            else "NA"
        ),

        "best_local_same_fbgn_distance_bp": (
            f"{best_local['distance']:.1f}"
            if best_local
            else "NA"
        ),

        "shared_fbgn_with_best_peak": (
            shared_fbgn
        ),

        "dmel_target_fbgn": (
            "|".join(sorted(ref_genes))
            if ref_genes
            else "NA"
        ),

        "gene_support_at_best_peak": (
            gene_support
        ),

        "class": cls,
    })


# ============================================================
# 5. Write result
# ============================================================

fields = [
    "dmel_cre_id",
    "species",
    "alignment_status",

    "target_chrom",
    "target_start0",
    "target_end0",
    "lifted_length",

    "best_peak_id",
    "best_overlap_bp",
    "best_overlap_fraction_lifted",
    "best_overlap_fraction_peak",

    "positional_peak",

    "same_fbgn_peak_anywhere",
    "same_fbgn_peak_local",

    "local_same_fbgn_peak_id",
    "local_same_fbgn",
    "best_local_same_fbgn_distance_bp",

    "shared_fbgn_with_best_peak",
    "dmel_target_fbgn",
    "gene_support_at_best_peak",

    "class",
]


with open(
    args.out,
    "w",
    newline="",
) as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(output_rows)


# ============================================================
# 6. Summary
# ============================================================

counts = defaultdict(int)

for row in output_rows:
    counts[row["class"]] += 1


mapped_n = sum(
    row["alignment_status"] == "mapped"
    for row in output_rows
)

unmapped_n = sum(
    row["alignment_status"] == "unmapped"
    for row in output_rows
)

same_any_n = sum(
    row["same_fbgn_peak_anywhere"] == "yes"
    for row in output_rows
)

same_local_n = sum(
    row["same_fbgn_peak_local"] == "yes"
    for row in output_rows
)


print(f"Species: {args.species}")
print(f"Reference CREs: {len(output_rows)}")
print(f"Mapped: {mapped_n}")
print(f"Unmapped: {unmapped_n}")

print(
    "Local same-FBgn distance threshold: "
    f"{args.max_local_gene_distance:.0f} bp"
)

print(
    f"Same-FBgn peak anywhere: {same_any_n}"
)

print(
    f"Same-FBgn peak local: {same_local_n}"
)

for cls in [
    "present",
    "turnover_candidate",
    "no_detected_CRE",
    "uncertain",
]:
    print(
        f"{cls}: {counts[cls]}"
    )

print(f"Wrote: {args.out}")
