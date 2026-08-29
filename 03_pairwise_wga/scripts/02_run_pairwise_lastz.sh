#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 TARGET_SPECIES" >&2
    exit 1
fi

SP="$1"

ROOT="$HOME/cre_turnover/project/pairwise_wga"

REF="d_melanogaster"

REF_2BIT="$ROOT/twobit/${REF}.2bit"
QRY_2BIT="$ROOT/twobit/${SP}.2bit"

REF_SIZES="$ROOT/chrom_sizes/${REF}.sizes"
QRY_SIZES="$ROOT/chrom_sizes/${SP}.sizes"

OUT="$ROOT/alignments_dmel/${SP}"

mkdir -p "$OUT"

AXT="$OUT/${REF}.${SP}.axt"

CHAIN_RAW="$OUT/${REF}.${SP}.raw.chain"
CHAIN_SORTED="$OUT/${REF}.${SP}.all.chain"

PRENET="$OUT/${REF}.${SP}.prenet.chain"

TNET="$OUT/${REF}.${SP}.target.net"
QNET="$OUT/${REF}.${SP}.query.net"

SYNTENIC_NET="$OUT/${REF}.${SP}.syntenic.net"

LIFTOVER_CHAIN="$OUT/${REF}.${SP}.liftover.chain"

echo "[$(date)] START $SP"
echo "Reference: $REF_2BIT"
echo "Query:     $QRY_2BIT"
echo "Output:    $OUT"

# ----------------------------------------------------------
# 1. Pairwise genome alignment
# ----------------------------------------------------------

lastz \
    "$REF_2BIT[multiple]" \
    "$QRY_2BIT[multiple]" \
    --format=axt \
    --ambiguous=iupac \
    --notransition \
    --step=20 \
    --seed=12of19 \
    --hspthresh=2200 \
    --gappedthresh=4000 \
    --ydrop=3400 \
    --inner=2000 \
    --output="$AXT"

echo "[$(date)] LASTZ done"

# ----------------------------------------------------------
# 2. AXT -> chain
# ----------------------------------------------------------

axtChain \
    -linearGap=medium \
    -minScore=3000 \
    "$AXT" \
    "$REF_2BIT" \
    "$QRY_2BIT" \
    "$CHAIN_RAW"

chainSort \
    "$CHAIN_RAW" \
    "$CHAIN_SORTED"

echo "[$(date)] chaining done"

# ----------------------------------------------------------
# 3. Prepare for netting
# ----------------------------------------------------------

chainPreNet \
    "$CHAIN_SORTED" \
    "$REF_SIZES" \
    "$QRY_SIZES" \
    "$PRENET"

# ----------------------------------------------------------
# 4. chain -> net
# ----------------------------------------------------------

chainNet \
    "$PRENET" \
    "$REF_SIZES" \
    "$QRY_SIZES" \
    "$TNET" \
    "$QNET"

netSyntenic \
    "$TNET" \
    "$SYNTENIC_NET"

echo "[$(date)] netting done"

# ----------------------------------------------------------
# 5. Extract chains represented in the syntenic net
# ----------------------------------------------------------

netChainSubset \
    "$SYNTENIC_NET" \
    "$CHAIN_SORTED" \
    "$LIFTOVER_CHAIN"

gzip -f "$LIFTOVER_CHAIN"

echo "[$(date)] DONE $SP"
echo "Final chain:"
echo "$LIFTOVER_CHAIN.gz"
