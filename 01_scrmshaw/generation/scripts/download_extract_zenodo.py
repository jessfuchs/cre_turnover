#!/usr/bin/env python3
"""
Lädt die Genome- und Annotationsarchive des Zenodo-Datensatzes herunter,
löst die in species.txt genannten Spezies gegen die Archivpfade auf und
extrahiert nur die ausgewählten Dateien.

species.txt:
    slug<TAB>Suchname
oder:
    Suchname

Bei uneindeutigen Treffern bricht das Skript bewusst ab und schreibt
Kandidaten in data/resolution_report.txt.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


GENOME_EXTS = (".fa", ".fasta", ".fna", ".fas")
ANNOT_EXTS = (".gff3", ".gff", ".gtf")


@dataclass(frozen=True)
class SpeciesRequest:
    index: int
    slug: str
    query: str


def die(message: str) -> "NoReturn":
    raise SystemExit(f"FEHLER: {message}")


def slugify(text: str) -> str:
    value = re.sub(r"[^A-Za-z0-9]+", "_", text.strip()).strip("_").lower()
    if not value:
        die(f"Kein gültiger slug aus {text!r}")
    return value


def normalize(text: str) -> str:
    text = text.lower()
    text = re.sub(r"\.(tar\.gz|gz|gff3|gff|gtf|fasta|fna|fa|fas)$", "", text)
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_")


def read_species(path: Path) -> list[SpeciesRequest]:
    requests: list[SpeciesRequest] = []
    seen: set[str] = set()

    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()

        if not line or line.startswith("#"):
            continue

        # Erwartetes Format:
        # slug    Scientific name
        #
        # Funktioniert sowohl mit Tabs als auch mit Leerzeichen.
        fields = line.split(maxsplit=1)

        if len(fields) != 2:
            die(
                f"Ungültige Zeile in {path}: {raw!r}\n"
                "Erwartet: slug<TAB/SPACE>Scientific name"
            )

        slug = slugify(fields[0])
        query = fields[1].strip()

        if not slug:
            die(f"Leerer slug in species.txt: {raw!r}")

        if not query:
            die(f"Leerer Species-Name für slug {slug}")

        if slug in seen:
            die(f"Doppelter slug in species.txt: {slug}")

        seen.add(slug)

        requests.append(
            SpeciesRequest(
                index=len(requests),
                slug=slug,
                query=query,
            )
        )

    if not requests:
        die("species.txt enthält keine Spezies.")

    return requests


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, target: Path, expected_md5: str | None = None) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        if expected_md5 and md5sum(target) == expected_md5:
            print(f"Vorhanden und MD5 korrekt: {target}")
            return
        if not expected_md5:
            print(f"Vorhanden: {target}")
            return
        print(f"Vorhandene Datei hat falsche MD5, lade erneut: {target}")
        target.unlink()

    cmd = [
        "wget", "-c", "--tries=10", "--timeout=60",
        "--output-document", str(target), url
    ]
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)
    if expected_md5:
        observed = md5sum(target)
        if observed != expected_md5:
            die(f"MD5 für {target.name}: erwartet {expected_md5}, erhalten {observed}")


def list_members(archive: Path) -> list[str]:
    print(f"Lese Inhaltsverzeichnis: {archive}")
    with tarfile.open(archive, "r:gz") as tar:
        return [m.name for m in tar.getmembers() if m.isfile()]


def strip_gz(name: str) -> str:
    return name[:-3] if name.lower().endswith(".gz") else name


def has_extension(name: str, exts: tuple[str, ...]) -> bool:
    base = strip_gz(name).lower()
    return base.endswith(exts)


def aliases(query: str) -> list[str]:
    q = normalize(query)
    parts = q.split("_")
    values = {q}
    if len(parts) >= 2:
        values.add("_".join(parts[-2:]))
        values.add(parts[-1])
        values.add(parts[0][0] + "_" + parts[-1])
    return sorted(values, key=len, reverse=True)


def score_candidate(member: str, request: SpeciesRequest, kind: str) -> int:
    path_norm = normalize(member)
    base_norm = normalize(Path(member).name)
    score = 0

    # Exakten Taxon-Dateistamm gegenüber längeren Unterartnamen bevorzugen.
    exact = re.sub(
        r"[^A-Za-z0-9]+",
        "_",
        request.query.strip()
    ).strip("_")

    basename = Path(member).name

    if kind == "genome" and re.match(
        rf"^{re.escape(exact)}\\.",
        basename,
        flags=re.IGNORECASE,
    ):
        score += 1000

    elif kind == "annotation" and re.match(
        rf"^{re.escape(exact)}_final\\.gff3?(?:\\.gz)?$",
        basename,
        flags=re.IGNORECASE,
    ):
        score += 1000

    for alias in aliases(request.query):
        if path_norm == alias or base_norm == alias:
            score = max(score, 200)
        elif f"_{alias}_" in f"_{path_norm}_":
            score = max(score, 150 + min(len(alias), 40))
        elif alias in path_norm:
            score = max(score, 80 + min(len(alias), 40))

    # Der explizite slug darf ebenfalls als Archivmuster dienen.
    slug_alias = normalize(request.slug)
    if slug_alias and slug_alias in path_norm:
        score = max(score, 70 + len(slug_alias))

    lower = member.lower()
    if kind == "genome":
        if any(x in lower for x in ("protein", "peptide", "pep.", "cds", "transcript", "rna.")):
            score -= 200
        if any(x in lower for x in ("genome", "assembly", "softmask", "masked")):
            score += 10
        if strip_gz(lower).endswith((".fna", ".fa", ".fasta")):
            score += 5
    else:
        if strip_gz(lower).endswith(".gff3"):
            score += 30
        elif strip_gz(lower).endswith(".gff"):
            score += 20
        elif strip_gz(lower).endswith(".gtf"):
            score += 10
        if any(x in lower for x in ("gene", "annotation", "cat", "braker")):
            score += 5
    return score


def resolve_one(
    members: list[str], request: SpeciesRequest, kind: str
) -> tuple[str | None, list[tuple[int, str]]]:
    valid = [
        m for m in members
        if has_extension(m, GENOME_EXTS if kind == "genome" else ANNOT_EXTS)
    ]

    query_norm = normalize(request.query)
    exact_matches = []

    for member in valid:
        basename = Path(member).name

        if kind == "genome":
            # Beispiel:
            # Drosophila_pseudoananassae.GCA_021223845.1.rm.fna.gz
            taxon_name = basename.split(".", 1)[0]

        else:
            # Beispiel:
            # DROSOPHILA_PSEUDOANANASSAE_final.gff.gz
            taxon_name = basename

            if taxon_name.lower().endswith(".gz"):
                taxon_name = taxon_name[:-3]

            for suffix in (".gff3", ".gff", ".gtf"):
                if taxon_name.lower().endswith(suffix):
                    taxon_name = taxon_name[:-len(suffix)]
                    break

            if taxon_name.lower().endswith("_final"):
                taxon_name = taxon_name[:-6]

        if normalize(taxon_name) == query_norm:
            exact_matches.append(member)

    if len(exact_matches) == 1:
        return exact_matches[0], [(10000, exact_matches[0])]

    ranked = sorted(
        ((score_candidate(m, request, kind), m) for m in valid),
        key=lambda x: (-x[0], x[1]),
    )
    ranked = [item for item in ranked if item[0] > 0]

    if not ranked:
        return None, []

    best_score = ranked[0][0]
    best = [m for score, m in ranked if score == best_score]

    if len(best) == 1:
        return best[0], ranked[:15]

    return None, ranked[:15]

def extract_members(archive: Path, members: list[str], target: Path) -> None:
    target.mkdir(parents=True, exist_ok=True)
    # Ein tar-Aufruf pro Archiv: das gzip-Archiv wird nur einmal durchlaufen.
    cmd = ["tar", "-xzf", str(archive), "-C", str(target), "--"] + members
    print(f"Extrahiere {len(members)} Dateien aus {archive.name}")
    subprocess.run(cmd, check=True)


def copy_decompress(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.name.lower().endswith(".gz"):
        with gzip.open(source, "rb") as src, target.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copy2(source, target)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True, type=Path)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--record-id", default="18453526")
    parser.add_argument("--genomes-md5", default="")
    parser.add_argument("--annotations-md5", default="")
    args = parser.parse_args()

    requests = read_species(args.species)
    download_dir = args.data_dir / "downloads"
    staging = args.data_dir / "archive_staging"
    selected = args.data_dir / "selected"
    report = args.data_dir / "resolution_report.txt"
    manifest = args.data_dir / "manifest.tsv"
    args.data_dir.mkdir(parents=True, exist_ok=True)

    base = f"https://zenodo.org/records/{args.record_id}/files"
    genomes_archive = download_dir / "genomes.tar.gz"
    annotations_archive = download_dir / "annotations.tar.gz"

    download(
        f"{base}/genomes.tar.gz?download=1",
        genomes_archive,
        args.genomes_md5 or None,
    )
    download(
        f"{base}/annotations.tar.gz?download=1",
        annotations_archive,
        args.annotations_md5 or None,
    )
    # Kleine Metadatendateien sind hilfreich zur manuellen Kontrolle.
    download(
        f"{base}/README.md?download=1",
        download_dir / "README.md",
        None,
    )
    download(
        f"{base}/Species_summary_301Fly.xlsx?download=1",
        download_dir / "Species_summary_301Fly.xlsx",
        None,
    )

    genome_members = list_members(genomes_archive)
    annotation_members = list_members(annotations_archive)

    resolved: list[tuple[SpeciesRequest, str, str]] = []
    problems: list[str] = []
    for req in requests:
        genome, genome_ranked = resolve_one(genome_members, req, "genome")
        annot, annot_ranked = resolve_one(annotation_members, req, "annotation")
        if not genome or not annot:
            problems.append(f"\n[{req.slug}] {req.query}")
            problems.append("  Genome-Kandidaten:")
            problems.extend(f"    {s:4d}  {m}" for s, m in genome_ranked)
            problems.append("  Annotation-Kandidaten:")
            problems.extend(f"    {s:4d}  {m}" for s, m in annot_ranked)
        else:
            resolved.append((req, genome, annot))

    report.write_text(
        "\n".join(problems) if problems else "Alle Spezies eindeutig aufgelöst.\n",
        encoding="utf-8",
    )
    if problems:
        die(
            f"{len(problems)} Auflösungsprobleme. Siehe {report}. "
            "Nutze in species.txt als zweite Spalte ein eindeutigeres Archiv-Muster."
        )

    if staging.exists():
        shutil.rmtree(staging)
    extract_members(genomes_archive, [g for _, g, _ in resolved], staging / "genomes")
    extract_members(
        annotations_archive, [a for _, _, a in resolved], staging / "annotations"
    )

    rows = ["index\tslug\tspecies\tgenome\tannotation\tannotation_format"]
    for req, genome_member, annotation_member in resolved:
        species_dir = selected / req.slug
        species_dir.mkdir(parents=True, exist_ok=True)

        genome_source = staging / "genomes" / genome_member
        annotation_source = staging / "annotations" / annotation_member

        genome_target = species_dir / "genome.fa"
        ann_no_gz = strip_gz(annotation_member).lower()
        if ann_no_gz.endswith(".gtf"):
            annotation_target = species_dir / "annotation.gtf"
            annotation_format = "gtf"
        else:
            annotation_target = species_dir / "annotation.gff3"
            annotation_format = "gff3"

        copy_decompress(genome_source, genome_target)
        copy_decompress(annotation_source, annotation_target)

        rows.append(
            "\t".join(
                [
                    str(req.index),
                    req.slug,
                    req.query,
                    str(genome_target.absolute()),
                    str(annotation_target.absolute()),
                    annotation_format,
                ]
            )
        )

    manifest.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"Manifest geschrieben: {manifest}")
    print(f"Aufgelöste Spezies: {len(resolved)}")


if __name__ == "__main__":
    main()
