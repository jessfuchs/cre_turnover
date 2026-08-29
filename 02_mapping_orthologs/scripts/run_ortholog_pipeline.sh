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
# Paths
# ------------------------------------------------------------

PROJECT="$HOME/cre_turnover/project"

ROOT="$PROJECT/mapping_orthologs"
SCRIPTS="$ROOT/scripts"

ORTHO="$ROOT/ortholog_results"
REFERENCE="$ROOT/reference_cres"

EXTERNAL="$PROJECT/external_scrmshaw"

MANIFEST="$EXTERNAL/combined_manifest.tsv"

DMEL_GFF="$EXTERNAL/external_data/selected/d_melanogaster/annotation.gff3"


# ------------------------------------------------------------
# Expected final project state
# ------------------------------------------------------------

EXPECTED_SPECIES=40
EXPECTED_REFERENCE_CRES=337


# ============================================================
# Prepare directories
# ============================================================

mkdir -p "$ORTHO"
mkdir -p "$REFERENCE"


# ============================================================
# Header
# ============================================================

echo "============================================================"
echo "Drosophila ortholog-mapping pipeline"
echo "============================================================"
echo "Started : $(date)"
echo "Project : $PROJECT"
echo "Root    : $ROOT"
echo


# ============================================================
# Check required scripts and inputs
# ============================================================

echo "[PRE] Checking required files..."

required_files=(

    "$MANIFEST"
    "$DMEL_GFF"

    "$SCRIPTS/01_qc_dmel_orthologs.sh"
    "$SCRIPTS/02_map_peaks_to_dmel_orthologs.py"
    "$SCRIPTS/03_SO_bed_to_tsv.py"
    "$SCRIPTS/04_map_dmel_ids_to_fbgn.py"
    "$SCRIPTS/05_qc_unresolved_dmel_ids.sh"
    "$SCRIPTS/06_make_dmel_reference_cres.py"

)

for f in "${required_files[@]}"; do

    if [[ ! -s "$f" ]]; then

        echo "ERROR: required file missing or empty:" >&2
        echo "  $f" >&2
        exit 1

    fi

done


# ============================================================
# 01. QC input Dmel ortholog annotations
# ============================================================

echo "============================================================"
echo "[01/06] QC of Dmel ortholog annotations"
echo "============================================================"
echo

bash "$SCRIPTS/01_qc_dmel_orthologs.sh"

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

python3 "$SCRIPTS/02_map_peaks_to_dmel_orthologs.py" \
    > "$ORTHO/ortholog_mapping.out" \
    2> "$ORTHO/ortholog_mapping_qc.tsv"


# ------------------------------------------------------------
# QC: one SO_all_peaks.bed per species
# ------------------------------------------------------------

N_SO=$(
    find "$ORTHO" \
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
# QC: files must not be empty
# ------------------------------------------------------------

EMPTY_SO=$(
    find "$ORTHO" \
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

python3 "$SCRIPTS/03_SO_bed_to_tsv.py" \
    "$ORTHO" \
    "$MANIFEST" \
    > "$ORTHO/SO_all_species.tsv" \
    2> "$ORTHO/SO_bed_to_tsv_qc.log"


[[ -s "$ORTHO/SO_all_species.tsv" ]] || {

    echo "ERROR: SO_all_species.tsv was not created." >&2
    exit 1

}


# ------------------------------------------------------------
# QC species representation
# ------------------------------------------------------------

N_COMBINED_SPECIES=$(
    python3 - "$ORTHO/SO_all_species.tsv" <<'PY'
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
echo "[02/05] PASSED"
echo


# ============================================================
# 04. Map Dmel identifiers to FBgn
# ============================================================

echo "============================================================"
echo "[04/06] Mapping Dmel identifiers to FBgn"
echo "============================================================"
echo

python3 "$SCRIPTS/04_map_dmel_ids_to_fbgn.py" \
    "$ORTHO/SO_all_species.tsv" \
    "$DMEL_GFF" \
    "$ORTHO/SO_all_species_fbgn.tsv" \
    2> "$ORTHO/fbgn_mapping_qc.log"


[[ -s "$ORTHO/SO_all_species_fbgn.tsv" ]] || {

    echo "ERROR: SO_all_species_fbgn.tsv was not created." >&2
    exit 1

}


[[ -f "$ORTHO/unresolved_dmel_identifiers.tsv" ]] || {

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

bash "$SCRIPTS/05_qc_unresolved_dmel_ids.sh"


[[ -s "$ORTHO/unresolved_dmel_qc_summary.tsv" ]] || {

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

python3 "$SCRIPTS/06_make_dmel_reference_cres.py" \
    > "$REFERENCE/reference_qc.log"


REF_TSV="$REFERENCE/dmel_reference_cres.tsv"
REF_BED="$REFERENCE/dmel_reference_cres.bed"


[[ -s "$REF_TSV" ]] || {

    echo "ERROR: reference CRE TSV was not created:" >&2
    echo "  $REF_TSV" >&2
    exit 1

}


[[ -s "$REF_BED" ]] || {

    echo "ERROR: reference CRE BED was not created:" >&2
    echo "  $REF_BED" >&2
    exit 1

}


# ------------------------------------------------------------
# Final reference CRE QC
# ------------------------------------------------------------

N_REF=$(
    grep -vcE '^(#|$)' "$REF_BED"
)

N_UNIQUE_IDS=$(
    cut -f4 "$REF_BED" \
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
echo "  $ORTHO/dmel_ortholog_annotation_qc.tsv"
echo

echo "Species-specific ortholog mapping:"
echo "  $ORTHO/<species>/SO_all_peaks.bed"
echo

echo "Combined mapping:"
echo "  $ORTHO/SO_all_species.tsv"
echo

echo "FBgn-normalized mapping:"
echo "  $ORTHO/SO_all_species_fbgn.tsv"
echo

echo "Unresolved-ID QC:"
echo "  $ORTHO/unresolved_dmel_identifiers.tsv"
echo "  $ORTHO/unresolved_dmel_qc_summary.tsv"
echo

echo "Final Dmel reference CRE set:"
echo "  $REF_TSV"
echo "  $REF_BED"
echo
echo "Reference CREs: $N_REF"
echo "============================================================"
