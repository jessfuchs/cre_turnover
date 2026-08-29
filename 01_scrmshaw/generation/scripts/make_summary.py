#!/usr/bin/env python3
from __future__ import annotations
import argparse
import csv
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--results", required=True, type=Path)
    args = p.parse_args()

    summary_rows = []
    with args.manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    for row in rows:
        slug = row["slug"]
        species = row["species"]
        bed = args.results / slug / "peaks_AllSets.bed"
        count = 0
        top = []
        if bed.exists():
            with bed.open() as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    count += 1
                    fields = line.rstrip("\n").split("\t")
                    try:
                        score = float(fields[4])
                    except (IndexError, ValueError):
                        score = float("-inf")
                    top.append((score, fields))
            top.sort(key=lambda x: x[0], reverse=True)
            top_path = args.results / slug / "top10_peaks.tsv"
            with top_path.open("w") as out:
                out.write("chrom\tstart\tend\tamplitude\tscore\tgene_left\tgene_right\ttraining\tmethod\trank\n")
                for _, f in top[:10]:
                    selected = [
                        f[0] if len(f) > 0 else "",
                        f[1] if len(f) > 1 else "",
                        f[2] if len(f) > 2 else "",
                        f[3] if len(f) > 3 else "",
                        f[4] if len(f) > 4 else "",
                        f[5] if len(f) > 5 else "",
                        f[10] if len(f) > 10 else "",
                        f[15] if len(f) > 15 else "",
                        f[16] if len(f) > 16 else "",
                        f[17] if len(f) > 17 else "",
                    ]
                    out.write("\t".join(selected) + "\n")
        summary_rows.append((slug, species, count, str(bed) if bed.exists() else "MISSING"))

    args.results.mkdir(parents=True, exist_ok=True)
    summary = args.results / "summary.tsv"
    with summary.open("w") as out:
        out.write("slug\tspecies\tn_peaks\tpeaks_file\n")
        for row in summary_rows:
            out.write("\t".join(map(str, row)) + "\n")
    print(summary)


if __name__ == "__main__":
    main()
