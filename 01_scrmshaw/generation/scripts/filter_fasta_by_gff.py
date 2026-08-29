#!/usr/bin/env python3
"""Behält nur FASTA-Sequenzen, deren seqid in der GFF vorkommt."""
from __future__ import annotations
import argparse
from pathlib import Path
from Bio import SeqIO


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--gff", required=True, type=Path)
    p.add_argument("--fasta", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()

    seqids: set[str] = set()
    with args.gff.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2] in {"gene", "exon", "CDS", "mRNA", "transcript"}:
                seqids.add(fields[0])

    if not seqids:
        raise SystemExit(f"Keine annotierten seqids in {args.gff} gefunden.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    found = 0
    total = 0
    with args.output.open("w") as out:
        for record in SeqIO.parse(str(args.fasta), "fasta"):
            total += 1
            if record.id in seqids:
                SeqIO.write(record, out, "fasta")
                found += 1

    missing = len(seqids) - found
    print(f"FASTA-Sequenzen gesamt: {total}")
    print(f"Behalten: {found}")
    print(f"Annotierte seqids ohne FASTA-Treffer: {missing}")
    if found == 0 or missing > 0:
        raise SystemExit("GFF/FASTA-seqid-Mismatch; Preflight-Log prüfen.")


if __name__ == "__main__":
    main()
