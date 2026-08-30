#!/usr/bin/env python3

# ============================================================
# 02 - Combine species-specific CRE classifications
#
# Purpose:
#   Combine per-species CRE classifications into long-format,
#   matrix-format, and species-level summary tables.
#
# Input:
#   - turnover_by_species/dmel_to_<species>_cre_turnover.tsv
#   - target_species.txt
#
# Output:
#   - cre_turnover_all_species.tsv
#   - cre_turnover_matrix.tsv
#   - species_summary.tsv
#
# Configuration:
#   Paths and expected reference CRE count are supplied by the
#   pipeline wrapper using config/classification_config.sh.
# ============================================================


import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


# ============================================================
# Command-line arguments
# ============================================================

if len(sys.argv) != 7:
    sys.exit(
        "Usage: 02_combine_turnover_results.py "
        "TURNOVER_BY_SPECIES_DIR TARGET_SPECIES_FILE "
        "LONG_OUT MATRIX_OUT SUMMARY_OUT N_REFERENCE_CRES"
    )

input_dir = Path(sys.argv[1])
target_file = Path(sys.argv[2])

long_out = Path(sys.argv[3])
matrix_out = Path(sys.argv[4])
summary_out = Path(sys.argv[5])

n_ref_cres = int(sys.argv[6])


# ============================================================
# Input validation
# ============================================================

if not input_dir.is_dir():
    sys.exit(
        f"ERROR: turnover-results directory does not exist: "
        f"{input_dir}"
    )

if not target_file.is_file():
    sys.exit(
        f"ERROR: target-species file does not exist: "
        f"{target_file}"
    )

for path in (
    long_out,
    matrix_out,
    summary_out,
):
    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )


# ============================================================
# Classification states
# ============================================================

VALID_CLASSES = {
    "present",
    "turnover_candidate",
    "no_detected_CRE",
    "uncertain",
}


# ============================================================
# Read target species
# ============================================================

species = [
    line.strip()
    for line in target_file.read_text().splitlines()
    if line.strip()
    and not line.lstrip().startswith("#")
]

if not species:
    sys.exit(
        f"ERROR: no target species found in {target_file}"
    )

if len(species) != len(set(species)):
    sys.exit(
        "ERROR: duplicate target species detected."
    )


# ============================================================
# Combine species-specific classifications
# ============================================================

all_rows = []
by_cre = defaultdict(dict)
summary = []

expected_header = None


for sp in species:

    infile = (
        input_dir
        / f"dmel_to_{sp}_cre_turnover.tsv"
    )

    if not infile.is_file():
        sys.exit(
            f"ERROR: missing result file: {infile}"
        )

    with infile.open(
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        if reader.fieldnames is None:
            sys.exit(
                f"ERROR: file has no header: {infile}"
            )

        if expected_header is None:
            expected_header = reader.fieldnames

        elif reader.fieldnames != expected_header:
            sys.exit(
                f"ERROR: header mismatch in {infile}"
            )

        rows = list(reader)


    # --------------------------------------------------------
    # Per-species QC
    # --------------------------------------------------------

    if len(rows) != n_ref_cres:
        sys.exit(
            f"ERROR: {sp}: expected {n_ref_cres} CRE rows, "
            f"found {len(rows)}"
        )

    cre_ids = [
        row["dmel_cre_id"]
        for row in rows
    ]

    if len(set(cre_ids)) != n_ref_cres:
        sys.exit(
            f"ERROR: {sp}: duplicate or missing "
            "dmel_cre_id values detected"
        )


    # --------------------------------------------------------
    # Counters
    # --------------------------------------------------------

    counts = Counter()

    mapped = 0
    unmapped = 0

    present_shared_fbgn = 0
    present_discordant_fbgn = 0
    present_no_fbgn = 0


    # --------------------------------------------------------
    # Process CRE classifications
    # --------------------------------------------------------

    for row in rows:

        cls = row["class"]

        if cls not in VALID_CLASSES:
            sys.exit(
                f"ERROR: {sp}: unexpected class '{cls}' "
                f"for {row['dmel_cre_id']}"
            )

        counts[cls] += 1

        alignment_status = row["alignment_status"]

        if alignment_status == "mapped":
            mapped += 1

        elif alignment_status == "unmapped":
            unmapped += 1

        else:
            sys.exit(
                f"ERROR: {sp}: unexpected alignment_status "
                f"'{alignment_status}'"
            )

        by_cre[
            row["dmel_cre_id"]
        ][sp] = cls

        all_rows.append(
            dict(row)
        )


        # ----------------------------------------------------
        # Gene support among positionally conserved CREs
        # ----------------------------------------------------

        if cls == "present":

            gene_support = row.get(
                "gene_support_at_best_peak",
                "NA"
            )

            if gene_support == "shared_fbgn":
                present_shared_fbgn += 1

            elif gene_support == "discordant_fbgn":
                present_discordant_fbgn += 1

            else:
                present_no_fbgn += 1


    # --------------------------------------------------------
    # Per-species summary
    # --------------------------------------------------------

    evaluable = mapped

    present = counts["present"]
    turnover = counts["turnover_candidate"]
    no_detected = counts["no_detected_CRE"]
    uncertain = counts["uncertain"]


    if mapped + unmapped != n_ref_cres:
        sys.exit(
            f"ERROR: {sp}: mapped + unmapped "
            f"!= {n_ref_cres}"
        )

    if sum(counts.values()) != n_ref_cres:
        sys.exit(
            f"ERROR: {sp}: class counts "
            f"!= {n_ref_cres}"
        )


    summary.append({
        "species": sp,
        "n_reference_cres": n_ref_cres,
        "mapped": mapped,
        "unmapped": unmapped,
        "mapping_rate": mapped / n_ref_cres,

        "present": present,
        "turnover_candidate": turnover,
        "no_detected_CRE": no_detected,
        "uncertain": uncertain,

        "present_rate_all": (
            present / n_ref_cres
        ),

        "present_rate_evaluable": (
            present / evaluable
            if evaluable
            else 0
        ),

        "turnover_rate_all": (
            turnover / n_ref_cres
        ),

        "turnover_rate_evaluable": (
            turnover / evaluable
            if evaluable
            else 0
        ),

        "no_detected_rate_all": (
            no_detected / n_ref_cres
        ),

        "no_detected_rate_evaluable": (
            no_detected / evaluable
            if evaluable
            else 0
        ),

        "present_shared_fbgn": (
            present_shared_fbgn
        ),

        "present_discordant_fbgn": (
            present_discordant_fbgn
        ),

        "present_other_gene_support": (
            present_no_fbgn
        ),
    })


# ============================================================
# Cross-species QC
# ============================================================

if len(by_cre) != n_ref_cres:
    sys.exit(
        f"ERROR: expected {n_ref_cres} unique reference CREs "
        f"across combined results, found {len(by_cre)}"
    )

for cre_id, species_map in by_cre.items():

    missing_species = (
        set(species)
        - set(species_map)
    )

    if missing_species:
        sys.exit(
            f"ERROR: {cre_id}: missing species: "
            + ",".join(
                sorted(missing_species)
            )
        )


# ============================================================
# 1. Long-format classification table
# ============================================================

species_order = {
    sp: i
    for i, sp in enumerate(species)
}

all_rows.sort(
    key=lambda row: (
        row["dmel_cre_id"],
        species_order[row["species"]],
    )
)

with long_out.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=expected_header,
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(all_rows)


# ============================================================
# 2. CRE-by-species classification matrix
# ============================================================

with matrix_out.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow(
        ["dmel_cre_id"]
        + species
    )

    for cre_id in sorted(by_cre):

        writer.writerow(
            [cre_id]
            + [
                by_cre[cre_id][sp]
                for sp in species
            ]
        )


# ============================================================
# 3. Species-level summary
# ============================================================

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


with summary_out.open(
    "w",
    encoding="utf-8",
    newline=""
) as handle:

    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=summary_fields,
        lineterminator="\n",
    )

    writer.writeheader()

    for row in summary:

        output_row = dict(row)

        for key in (
            "mapping_rate",
            "present_rate_all",
            "present_rate_evaluable",
            "turnover_rate_all",
            "turnover_rate_evaluable",
            "no_detected_rate_all",
            "no_detected_rate_evaluable",
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

print("Combination complete")
print(f"Target species: {len(species)}")
print(f"Dmel CREs: {len(by_cre)}")
print(f"Long-table rows: {len(all_rows)}")
print()
print(f"Wrote: {long_out}")
print(f"Wrote: {matrix_out}")
print(f"Wrote: {summary_out}")
