# SCRMshaw-HD + Postprocessing für viele Drosophilidae-Spezies auf Slurm

Diese Pipeline:

1. liest die gewünschten Spezies aus `species.txt`,
2. lädt `genomes.tar.gz` und `annotations.tar.gz` aus dem Zenodo-Datensatz
   zum Paper *Comparative gene annotation and orthology assignments across
   301 species of Drosophilidae*,
3. extrahiert nur die ausgewählten Genome und Annotationen,
4. prüft GFF3/FASTA, entfernt unannotierte Sequenzen und führt optional TRF-Masking durch,
5. startet 25 SCRMshaw-HD-Offsets pro Spezies,
6. führt das HalfonLab-Postprocessing pro Spezies aus,
7. sammelt `peaks_AllSets.bed`, Peak-Zahlen und Top-10-Tabellen.

## Wichtige Designentscheidung

Jeder neue `task_offset_*`-Ordner wird mit `--step 123` gestartet.
`--step 23` reicht in einem frischen Outdir nicht aus, weil dort unter anderem
`gff/genes` und die offset-spezifischen Fenster fehlen.

Standardmäßig wird nur `--imm` verwendet. Das ist für 22 Spezies wesentlich
schneller. In `config/config.sh` kann auf

```bash
SCORING_FLAGS="--imm --hexmcd --pac"
```

umgestellt werden.

## 1. Conda-Umgebungen

```bash
conda env create -f envs/scrmshaw.yml
conda env create -f envs/postprocess.yml
```

Bei abweichenden Environment-Namen `SCRM_ENV` und `POST_ENV` in
`config/config.sh` ändern.

## 2. Software holen

```bash
bash scripts/setup_software.sh
```

Dadurch werden SCRMshaw_HD, UtilityPrograms, das Postprocessing-Repo,
die Trainingssets und TRF heruntergeladen. Verwendet wird:

```text
final_combined_48Tsets/adult_muscle
```

## 3. Speziesliste ausfüllen

`species.txt`:

```text
dana	Drosophila ananassae
dpse	Drosophila pseudoobscura
```

Erste Spalte: kurzer eindeutiger Ordnername.  
Zweite Spalte: wissenschaftlicher Name oder ein eindeutiges Muster aus den
Zenodo-Archivpfaden.

Die Pipeline funktioniert für beliebig viele Spezies, nicht nur 22. Bei
uneindeutigen Dateinamen stoppt der Download-Job und schreibt Kandidaten nach:

```text
data/resolution_report.txt
```

Dann wird die zweite Spalte in `species.txt` präzisiert.

## 4. Cluster-Konfiguration

In `config/config.sh` bei Bedarf setzen:

- `SLURM_PARTITION`
- `SLURM_QOS`
- `SLURM_ACCOUNT`
- `MAX_OFFSET_JOBS`
- Laufzeit- und RAM-Limits
- `CONDA_BASE`

Die Standardwerte für einen Offset-Task sind 1 CPU, 48 GB RAM und 72 Stunden.

## 5. Pipeline starten

```bash
bash submit_pipeline.sh
```

Für 22 Spezies entstehen 550 Offset-Tasks. Die Abhängigkeiten werden
automatisch gesetzt:

```text
Download → 22× Preprocessing → 550× SCRMshaw → 22× Postprocessing → Summary
```

## Outputs

Pro Spezies:

```text
results/<slug>/peaks_AllSets.bed
results/<slug>/top10_peaks.tsv
results/<slug>/peak_count.txt
```

Gesamtübersicht:

```text
results/summary.tsv
```

## Reruns

Fertige Offset-Tasks werden anhand von `SCRM_FINISHED.ok` übersprungen.
Das Postprocessing räumt nur seine eigenen Zwischenoutputs auf und löscht keine
`task_offset_*`-Ordner.

## Speicherplatz

Die beiden Zenodo-Archive belegen zusammen rund 21.4 GB komprimiert. Die
SCRMshaw-HD-Zwischenverzeichnisse können wesentlich größer werden. Vor einem
22-Spezies-Lauf sollte großzügig freier Speicher eingeplant werden; mehrere
hundert GB sind je nach Genomgrößen und Scoring-Methoden realistisch.
