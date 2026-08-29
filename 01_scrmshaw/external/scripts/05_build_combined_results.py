#!/usr/bin/env python3
from pathlib import Path
import argparse,csv,shutil
ap=argparse.ArgumentParser()
for n in ["generated_manifest","external_manifest","generated_results","external_results","combined_manifest","combined_results"]:
    ap.add_argument("--"+n.replace("_","-"),dest=n,type=Path,required=True)
a=ap.parse_args()
rows=[]
with a.generated_manifest.open() as h:
    for r in csv.DictReader(h,delimiter="\t"):
        rows.append((r["slug"].strip(),r["species"].strip(),"generated"))
with a.external_manifest.open() as h:
    for r in csv.DictReader(h,delimiter="\t"):
        rows.append((r["slug"].strip(),r["species"].strip(),"external"))
slugs=[r[0] for r in rows]
if len(slugs)!=len(set(slugs)): raise SystemExit("FEHLER: Doppelte Slugs")
a.combined_manifest.parent.mkdir(parents=True,exist_ok=True)
with a.combined_manifest.open("w") as h:
    h.write("index\tslug\tspecies\tsource\n")
    for i,(s,sp,src) in enumerate(rows): h.write(f"{i}\t{s}\t{sp}\t{src}\n")
if a.combined_results.exists(): shutil.rmtree(a.combined_results)
a.combined_results.mkdir(parents=True)
missing=[]
for slug,species,src in rows:
    source=(a.generated_results if src=="generated" else a.external_results)/slug/"peaks_AllSets.bed"
    if not source.is_file() or source.stat().st_size==0:
        missing.append(source); continue
    d=a.combined_results/slug; d.mkdir()
    (d/"peaks_AllSets.bed").symlink_to(source.resolve())
print(f"Manifest-Arten={len(rows)} Symlinks={len(rows)-len(missing)}")
if missing:
    print("Fehlend:"); [print(x) for x in missing]; raise SystemExit(1)
