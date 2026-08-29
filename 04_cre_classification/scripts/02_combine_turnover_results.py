#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter, defaultdict

BASE = Path.home() / "cre_turnover/project/cre_classification"
RESULTS = BASE / "results"

IN_DIR = RESULTS / "turnover_by_species"
TARGETS = Path.home() / "cre_turnover/project/pairwise_wga/target_species.txt"

LONG_OUT = RESULTS / "cre_turnover_all_species.tsv"
MATRIX_OUT = RESULTS / "cre_turnover_matrix.tsv"
SUMMARY_OUT = RESULTS / "species_summary.tsv"

RESULTS.mkdir(parents=True, exist_ok=True)
IN_DIR.mkdir(parents=True, exist_ok=True)

N_REF_CRES = 337

VALID_CLASSES = {
    "present",
    "turnover_candidate",
    "no_detected_CRE",
    "uncertain",
}

species = [
    x.strip()
    for x in TARGETS.read_text().splitlines()
    if x.strip()
]

all_rows = []
by_cre = defaultdict(dict)
summary = []


# ---------------------------------------------------------
# Read all per-species results
# ---------------------------------------------------------

expected_header = None

for sp in species:

    infile = IN_DIR / f"dmel_to_{sp}_cre_turnover.tsv"

    if not infile.exists():
        raise SystemExit(f"Missing result file: {infile}")

    with infile.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        if expected_header is None:
            expected_header = reader.fieldnames
        elif reader.fieldnames != expected_header:
            raise SystemExit(
                f"Header mismatch in {infile}"
            )

        rows = list(reader)

    if len(rows) != N_REF_CRES:
        raise SystemExit(
            f"{sp}: expected {N_REF_CRES} CRE rows, found {len(rows)}"
        )

    cre_ids = [r["dmel_cre_id"] for r in rows]

    if len(set(cre_ids)) != 337:
        raise SystemExit(
            f"{sp}: duplicate dmel_cre_id values detected"
        )

    counts = Counter()

    mapped = 0
    unmapped = 0

    present_shared_fbgn = 0
    present_discordant_fbgn = 0
    present_no_fbgn = 0

    for r in rows:

        cls = r["class"]

        if cls not in VALID_CLASSES:
            raise SystemExit(
                f"{sp}: unexpected class '{cls}' "
                f"for {r['dmel_cre_id']}"
            )

        counts[cls] += 1

        if r["alignment_status"] == "mapped":
            mapped += 1
        elif r["alignment_status"] == "unmapped":
            unmapped += 1
        else:
            raise SystemExit(
                f"{sp}: unexpected alignment_status "
                f"'{r['alignment_status']}'"
            )

        by_cre[r["dmel_cre_id"]][sp] = cls

        row_out = dict(r)
        all_rows.append(row_out)

        if cls == "present":
            gs = r.get("gene_support_at_best_peak", "NA")

            if gs == "shared_fbgn":
                present_shared_fbgn += 1
            elif gs == "discordant_fbgn":
                present_discordant_fbgn += 1
            else:
                present_no_fbgn += 1

    evaluable = mapped

    present = counts["present"]
    turnover = counts["turnover_candidate"]
    no_detected = counts["no_detected_CRE"]
    uncertain = counts["uncertain"]

    if mapped + unmapped != N_REF_CRES:
        raise SystemExit(
            f"{sp}: mapped + unmapped != 337"
        )

    if sum(counts.values()) != N_REF_CRES:
        raise SystemExit(
            f"{sp}: class counts != 337"
        )

    summary.append({
        "species": sp,
        "n_reference_cres": N_REF_CRES,
        "mapped": mapped,
        "unmapped": unmapped,
        "mapping_rate": mapped / N_REF_CRES,
        "present": present,
        "turnover_candidate": turnover,
        "no_detected_CRE": no_detected,
        "uncertain": uncertain,
        "present_rate_all": present / N_REF_CRES,
        "present_rate_evaluable": (
            present / evaluable if evaluable else 0
        ),
        "turnover_rate_all": turnover / N_REF_CRES,
        "turnover_rate_evaluable": (
            turnover / evaluable if evaluable else 0
        ),
        "no_detected_rate_all": no_detected / N_REF_CRES,
        "no_detected_rate_evaluable": (
            no_detected / evaluable if evaluable else 0
        ),
        "present_shared_fbgn": present_shared_fbgn,
        "present_discordant_fbgn": present_discordant_fbgn,
        "present_other_gene_support": present_no_fbgn,
    })


# ---------------------------------------------------------
# Check that every CRE exists for every species
# ---------------------------------------------------------

if len(by_cre) != N_REF_CRES:
    raise SystemExit(
        f"{sp}: expected {N_REF_CRES} CRE rows, found {len(rows)}"
    )

for cre_id, sp_map in by_cre.items():
    missing = set(species) - set(sp_map)

    if missing:
        raise SystemExit(
            f"{cre_id}: missing species: "
            + ",".join(sorted(missing))
        )


# ---------------------------------------------------------
# 1. Long-format table
# ---------------------------------------------------------

with LONG_OUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=expected_header,
        lineterminator="\n",
    )

    writer.writeheader()

    all_rows.sort(
        key=lambda r: (
            r["dmel_cre_id"],
            species.index(r["species"])
        )
    )

    writer.writerows(all_rows)


# ---------------------------------------------------------
# 2. 337 x 22 class matrix
# ---------------------------------------------------------

with MATRIX_OUT.open("w", newline="") as f:
    writer = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow(
        ["dmel_cre_id"] + species
    )

    for cre_id in sorted(by_cre):
        writer.writerow(
            [cre_id] +
            [by_cre[cre_id][sp] for sp in species]
        )


# ---------------------------------------------------------
# 3. Species summary
# ---------------------------------------------------------

summary_fields = [
    "species",
    "n_reference_cres",
    "mapped",
    "unmapped",
    "mapping_rate",
    "present",
    "turnover_candidate",
    "no_detected_CRE",
    "uncertain",
    "present_rate_all",
    "present_rate_evaluable",
    "turnover_rate_all",
    "turnover_rate_evaluable",
    "no_detected_rate_all",
    "no_detected_rate_evaluable",
    "present_shared_fbgn",
    "present_discordant_fbgn",
    "present_other_gene_support",
]

with SUMMARY_OUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=summary_fields,
        lineterminator="\n",
    )

    writer.writeheader()

    for r in summary:
        rr = dict(r)

        for key in (
            "mapping_rate",
            "present_rate_all",
            "present_rate_evaluable",
            "turnover_rate_all",
            "turnover_rate_evaluable",
            "no_detected_rate_all",
            "no_detected_rate_evaluable",
        ):
            rr[key] = f"{rr[key]:.6f}"

        writer.writerow(rr)


# ---------------------------------------------------------
# Report
# ---------------------------------------------------------

print("Combination complete")
print(f"Species: {len(species)}")
print(f"Dmel CREs: {len(by_cre)}")
print(f"Long-table rows: {len(all_rows)}")
print()
print(f"Wrote: {LONG_OUT}")
print(f"Wrote: {MATRIX_OUT}")
print(f"Wrote: {SUMMARY_OUT}")
