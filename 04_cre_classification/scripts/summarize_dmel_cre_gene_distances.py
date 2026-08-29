#!/usr/bin/env python3

import csv
import math
from pathlib import Path

INFILE = Path("SO_all_species_fbgn.tsv")
OUTFILE = Path("dmel_cre_gene_distance_summary.tsv")

pairs = []

with INFILE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        if row["species_key"] != "dmelanogaster":
            continue

        candidates = [
            (
                row.get("dmel_fbgn_flanking_gene", ""),
                row.get("distance_flanking_gene", ""),
                "flanking",
            ),
            (
                row.get("dmel_fbgn_next_flanking_gene", ""),
                row.get("distance_next_gene", ""),
                "next",
            ),
        ]

        for fbgn_field, dist_field, which in candidates:
            fbgns = [
                x.strip()
                for x in fbgn_field.split("|")
                if x.strip().startswith("FBgn")
            ]

            if not fbgns:
                continue

            try:
                dist = abs(float(dist_field))
            except (TypeError, ValueError):
                continue

            for fbgn in fbgns:
                pairs.append({
                    "fbgn": fbgn,
                    "distance_bp": dist,
                    "gene_relation": which,
                })


distances = sorted(x["distance_bp"] for x in pairs)


def percentile(values, p):
    if not values:
        return float("nan")

    k = (len(values) - 1) * p / 100
    lo = math.floor(k)
    hi = math.ceil(k)

    if lo == hi:
        return values[lo]

    return values[lo] * (hi - k) + values[hi] * (k - lo)


with OUTFILE.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
    writer.writerow(["percentile", "distance_bp"])

    for p in [25, 50, 75, 90, 95, 97.5, 99, 100]:
        writer.writerow([p, f"{percentile(distances, p):.1f}"])


print(f"N CRE-gene associations: {len(distances)}")

for p in [25, 50, 75, 90, 95, 97.5, 99, 100]:
    print(f"{p:5g}th percentile: {percentile(distances, p):.1f} bp")

print(f"Wrote: {OUTFILE}")
