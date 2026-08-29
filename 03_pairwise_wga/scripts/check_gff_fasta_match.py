#!/usr/bin/env python3

import sys
import gzip
from collections import defaultdict

if len(sys.argv) != 3:
    sys.exit(
        "Usage: check_gff_fasta_match.py annotation.gff3 genome.fa[.gz]"
    )

gff = sys.argv[1]
fasta = sys.argv[2]


def opener(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


# -------------------------
# Read FASTA sequence sizes
# -------------------------

lengths = {}
name = None
n = 0

with opener(fasta) as f:
    for line in f:
        if line.startswith(">"):
            if name is not None:
                lengths[name] = n

            name = line[1:].split()[0]
            n = 0
        else:
            n += len(line.strip())

    if name is not None:
        lengths[name] = n


# -------------------------
# Read maximum GFF coordinate
# -------------------------

gff_max = defaultdict(int)

with opener(gff) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue

        cols = line.rstrip("\n").split("\t")

        if len(cols) < 5:
            continue

        try:
            end = int(cols[4])
        except ValueError:
            continue

        seqid = cols[0]
        gff_max[seqid] = max(gff_max[seqid], end)


# -------------------------
# Compare
# -------------------------

missing = []
too_short = []

for seqid, max_end in gff_max.items():

    if seqid not in lengths:
        missing.append(seqid)

    elif lengths[seqid] < max_end:
        too_short.append(
            (seqid, max_end, lengths[seqid])
        )


print(f"GFF seqIDs: {len(gff_max)}")
print(f"FASTA seqIDs: {len(lengths)}")
print(f"Missing in FASTA: {len(missing)}")
print(f"Too short: {len(too_short)}")


if missing:
    print("\nFirst missing IDs:")
    for x in missing[:20]:
        print(x)


if too_short:
    print("\nFirst too-short sequences:")
    for seqid, need, have in too_short[:20]:
        print(
            f"{seqid}\tGFF_end={need}\tFASTA_len={have}"
        )


if not missing and not too_short:
    print("\nMATCH_OK")
