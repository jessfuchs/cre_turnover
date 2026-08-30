#!/usr/bin/env bash
set -euo pipefail


# ============================================================
# Complete ortholog-mapping pipeline
#
# Steps:
#
#   01  QC Dmel ortholog annotations in input GFFs
#   02  Map SCRMshaw peak-associated genes to Dmel orthologs
#   03  Combine species-specific SO BED files into one TSV
#   04  Normalize Dmel identifiers to FlyBase FBgn IDs
#   05  QC unresolved Dmel identifiers
#   06  Build final Dmel reference CRE set
#
# Main outputs:
#
#   ortholog_results/SO_all_species.tsv
#   ortholog_results/SO_all_species_fbgn.tsv
#   reference_cres/dmel_reference_cres.tsv
#   reference_cres/dmel_reference_cres.bed
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

source "$PIPELINE_ROOT/config/ortholog_config.sh"


# ============================================================
# Prepare directories
# ============================================================

mkdir -p \
    "$ORTHOLOG_RESULTS_DIR" \
    "$REFERENCE_CRES_DIR"


# ============================================================
# Header
# ============================================================

echo "============================================================"
echo "Drosophila ortholog-mapping pipeline"
echo "============================================================"
echo "Started      : $(date)"
echo "Project root : $PROJECT_ROOT"
echo "Pipeline root: $PIPELINE_ROOT"
echo "Manifest     : $COMBINED_MANIFEST"
echo "Dmel GFF     : $DMEL_GFF"
echo


# ============================================================
# Check required scripts and inputs
# ============================================================

echo "[PRE] Checking required files..."

required_files=(

    "$COMBINED_MANIFEST"
    "$DMEL_GFF"

    "$SCRIPTS_DIR/01_qc_dmel_orthologs.sh"
    "$SCRIPTS_DIR/02_map_peaks_to_dmel_orthologs.py"
    "$SCRIPTS_DIR/03_SO_bed_to_tsv.py"
    "$SCRIPTS_DIR/04_map_dmel_ids_to_fbgn.py"
    "$SCRIPTS_DIR/05_qc_unresolved_dmel_ids.sh"
    "$SCRIPTS_DIR/06_make_dmel_reference_cres.py"
)

for f in "${required_files[@]}"; do

    if [[ ! -s "$f" ]]; then

        echo "ERROR: required file missing or empty:" >&2
        echo "  $f" >&2
        exit 1

    fi

done

echo "[PRE] PASSED"
echo


# ============================================================
# 01. QC input Dmel ortholog annotations
# ============================================================

echo "============================================================"
echo "[01/06] QC of Dmel ortholog annotations"
echo "============================================================"
echo

bash "$SCRIPTS_DIR/01_qc_dmel_orthologs.sh"

echo
echo "[01/06] PASSED"
echo


# ============================================================
# 02. Map SCRMshaw peak genes to Dmel orthologs
# ============================================================

echo "============================================================"
echo "[02/06] Mapping peak-associated genes to Dmel orthologs"
echo "============================================================"
echo

python3 "$SCRIPTS_DIR/02_map_peaks_to_dmel_orthologs.py" \
    > "$ORTHOLOG_RESULTS_DIR/ortholog_mapping.out" \
    2> "$ORTHOLOG_RESULTS_DIR/ortholog_mapping_qc.tsv"


# ------------------------------------------------------------
# QC: one SO_all_peaks.bed per species
# ------------------------------------------------------------

N_SO=$(
    find "$ORTHOLOG_RESULTS_DIR" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name 'SO_all_peaks.bed' \
        | wc -l
)

echo "SO_all_peaks.bed files: $N_SO"

if [[ "$N_SO" -ne "$EXPECTED_SPECIES" ]]; then

    echo "ERROR: expected $EXPECTED_SPECIES species-specific" >&2
    echo "SO_all_peaks.bed files, found $N_SO." >&2
    exit 1

fi


# ------------------------------------------------------------
# QC: species-specific files must not be empty
# ------------------------------------------------------------

EMPTY_SO=$(
    find "$ORTHOLOG_RESULTS_DIR" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name 'SO_all_peaks.bed' \
        -size 0 \
        | wc -l
)

if [[ "$EMPTY_SO" -ne 0 ]]; then

    echo "ERROR: $EMPTY_SO SO_all_peaks.bed files are empty." >&2
    exit 1

fi

echo
echo "[02/06] PASSED"
echo


# ============================================================
# 03. Combine all species
# ============================================================

echo "============================================================"
echo "[03/06] Combining species-specific SO BED files"
echo "============================================================"
echo

python3 "$SCRIPTS_DIR/03_SO_bed_to_tsv.py" \
    "$ORTHOLOG_RESULTS_DIR" \
    "$COMBINED_MANIFEST" \
    > "$SO_ALL_SPECIES" \
    2> "$ORTHOLOG_RESULTS_DIR/SO_bed_to_tsv_qc.log"


[[ -s "$SO_ALL_SPECIES" ]] || {

    echo "ERROR: SO_all_species.tsv was not created." >&2
    exit 1

}


# ------------------------------------------------------------
# QC species representation
# ------------------------------------------------------------

N_COMBINED_SPECIES=$(
    python3 - "$SO_ALL_SPECIES" <<'PY'
import csv
import sys

path = sys.argv[1]

with open(
    path,
    encoding="utf-8-sig",
    newline=""
) as handle:

    reader = csv.DictReader(
        handle,
        delimiter="\t"
    )

    species = {
        row["species_key"].strip()
        for row in reader
        if row.get("species_key", "").strip()
    }

print(len(species))
PY
)

echo "Species in SO_all_species.tsv: $N_COMBINED_SPECIES"

if [[ "$N_COMBINED_SPECIES" -ne "$EXPECTED_SPECIES" ]]; then

    echo "ERROR: expected $EXPECTED_SPECIES species in combined TSV." >&2
    exit 1

fi

echo
echo "[03/06] PASSED"
echo


# ============================================================
# 04. Map Dmel identifiers to FBgn
# ============================================================

echo "============================================================"
echo "[04/06] Mapping Dmel identifiers to FBgn"
echo "============================================================"
echo

python3 "$SCRIPTS_DIR/04_map_dmel_ids_to_fbgn.py" \
    "$SO_ALL_SPECIES" \
    "$DMEL_GFF" \
    "$SO_ALL_SPECIES_FBGN" \
    2> "$ORTHOLOG_RESULTS_DIR/fbgn_mapping_qc.log"


[[ -s "$SO_ALL_SPECIES_FBGN" ]] || {

    echo "ERROR: SO_all_species_fbgn.tsv was not created." >&2
    exit 1

}


[[ -f "$UNRESOLVED_DMEL" ]] || {

    echo "ERROR: unresolved_dmel_identifiers.tsv was not created." >&2
    exit 1

}

echo
echo "[04/06] PASSED"
echo


# ============================================================
# 05. QC unresolved Dmel identifiers
# ============================================================

echo "============================================================"
echo "[05/06] QC of unresolved Dmel identifiers"
echo "============================================================"
echo

bash "$SCRIPTS_DIR/05_qc_unresolved_dmel_ids.sh"


[[ -s "$UNRESOLVED_QC_SUMMARY" ]] || {

    echo "ERROR: unresolved Dmel QC summary was not created." >&2
    exit 1

}

echo
echo "[05/06] PASSED"
echo


# ============================================================
# 06. Build final Dmel reference CRE set
# ============================================================

echo "============================================================"
echo "[06/06] Building final Dmel reference CRE set"
echo "============================================================"
echo

python3 "$SCRIPTS_DIR/06_make_dmel_reference_cres.py" \
    > "$REFERENCE_CRES_DIR/reference_qc.log"


[[ -s "$REFERENCE_CRES_TSV" ]] || {

    echo "ERROR: reference CRE TSV was not created:" >&2
    echo "  $REFERENCE_CRES_TSV" >&2
    exit 1

}


[[ -s "$REFERENCE_CRES_BED" ]] || {

    echo "ERROR: reference CRE BED was not created:" >&2
    echo "  $REFERENCE_CRES_BED" >&2
    exit 1

}


# ------------------------------------------------------------
# Final reference CRE QC
# ------------------------------------------------------------

N_REF=$(
    grep -vcE '^(#|$)' "$REFERENCE_CRES_BED"
)

N_UNIQUE_IDS=$(
    cut -f4 "$REFERENCE_CRES_BED" \
        | sort -u \
        | wc -l
)

echo "Reference CREs       : $N_REF"
echo "Unique reference IDs : $N_UNIQUE_IDS"


if [[ "$N_REF" -ne "$N_UNIQUE_IDS" ]]; then

    echo "ERROR: duplicated Dmel CRE identifiers detected." >&2
    exit 1

fi


if [[ "$N_REF" -ne "$EXPECTED_REFERENCE_CRES" ]]; then

    echo "WARNING: expected currently $EXPECTED_REFERENCE_CRES reference CREs,"
    echo "but found $N_REF."
    echo "This is not necessarily an error if the upstream"
    echo "reference definition has intentionally changed."

fi

echo
echo "[06/06] PASSED"
echo


# ============================================================
# Final summary
# ============================================================

echo "============================================================"
echo "ORTHOLOG-MAPPING PIPELINE COMPLETE"
echo "============================================================"
echo "Finished: $(date)"
echo

echo "Main outputs:"
echo

echo "Annotation QC:"
echo "  $DMEL_ORTHOLOG_QC"
echo

echo "Species-specific ortholog mapping:"
echo "  $ORTHOLOG_RESULTS_DIR/<species>/SO_all_peaks.bed"
echo

echo "Combined mapping:"
echo "  $SO_ALL_SPECIES"
echo

echo "FBgn-normalized mapping:"
echo "  $SO_ALL_SPECIES_FBGN"
echo

echo "Unresolved-ID QC:"
echo "  $UNRESOLVED_DMEL"
echo "  $UNRESOLVED_QC_SUMMARY"
echo

echo "Final Dmel reference CRE set:"
echo "  $REFERENCE_CRES_TSV"
echo "  $REFERENCE_CRES_BED"
echo

echo "Reference CREs: $N_REF"
echo "============================================================"
