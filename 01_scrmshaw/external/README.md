# Externe SCRMshaw-Dateien integrieren

Das Paket automatisiert:

1. GFFs der externen Arten aus `annotations.tar.gz` extrahieren.
2. Betreuer-BEDs auf `adult_muscle + imm` und 17 Spalten filtern.
3. BED/GFF-SeqIDs validieren.
4. `D_ananassae` und `D_pseudoobscura` ausschließen.
5. Externes Postprocessing als Slurm-Array auf `abacus-2`, 5 GB, maximal ein Task.
6. Ursprüngliche 22 + externe 12 Arten als Symlinks in `combined_results/` vereinen.
7. Gemeinsame `combined_results/summary.tsv` erzeugen.

## Verwendung

Kopiere/entpacke den Inhalt in:

```bash
~/projects/cre_turnover/scrmshaw_22species_slurm
```

Pfade bei Bedarf in `config/external_config.sh` ändern.

Start:

```bash
bash run_external_import_pipeline.sh
```

Erwartet:

```bash
find results -name peaks_AllSets.bed -type f -size +0c | wc -l
# 22
find external_results -name peaks_AllSets.bed -type f -size +0c | wc -l
# 12
find combined_results -name peaks_AllSets.bed -type l | wc -l
# 34
wc -l combined_results/summary.tsv
# 35
```
