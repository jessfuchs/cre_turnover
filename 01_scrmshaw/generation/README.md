# Multi-species SCRMshaw-HD Prediction Pipeline

This workflow generates and post-processes SCRMshaw-HD cis-regulatory element (CRE) predictions for the species listed in `species.txt`.

The pipeline standardizes genome acquisition, genome and annotation preparation, SCRMshaw-HD execution, and species-level peak calling across newly analyzed Drosophilidae species. The resulting `peaks_AllSets.bed` files are subsequently integrated with externally generated SCRMshaw predictions by the workflow in [`../external/`](../external/).

The final analysis generated predictions for 25 species. The pipeline itself is species-list driven and can be applied to a different number of species without changing the worker scripts.

---

## Workflow overview

For each species, the pipeline performs five main stages:

1. download the source genome and annotation archives and resolve the requested species;
2. validate and prepare the genome and annotation for SCRMshaw;
3. execute 25 independent SCRMshaw-HD offsets;
4. combine offset-level results and perform SCRMshaw peak calling;
5. generate species-level and cross-species result summaries.

Conceptually:

```text
species.txt
    │
    ▼
Zenodo genome + annotation archives
    │
    ▼
Species resolution and extraction
    │
    ▼
Genome / annotation preparation
    │
    ├─ preflight validation
    ├─ FASTA filtering
    └─ TRF masking
    │
    ▼
25 SCRMshaw-HD offsets per species
    │
    ▼
Top 5,000 hits per offset
    │
    ▼
SCRMshaw post-processing + peak calling
    │
    ▼
results/<species>/peaks_AllSets.bed
    │
    ▼
Cross-species summary
```

The complete workflow is submitted through SLURM using:

```bash
bash submit_pipeline.sh
```

---

## Repository structure

```text
01_scrmshaw/generation/
├── config/
│   └── config.sh
├── envs/
│   ├── scrmshaw.yml
│   └── postprocess.yml
├── scripts/
│   ├── 00_setup_software.sh
│   ├── 01_download_extract_zenodo.py
│   ├── 02_prepare_one_species.sh
│   ├── 03_run_one_offset.sh
│   ├── 04_postprocess_one_species.sh
│   └── 05_make_summary.py
├── slurm/
├── bin/
├── species.txt
├── qc_peaks_outputs.sh
└── submit_pipeline.sh
```

Runtime directories such as `data/`, `runs/`, `results/`, `logs/`, and `software/` are created or populated by the workflow and are not intended to contain manually maintained analysis code.

Paths, software locations, SCRMshaw settings, and SLURM resources are defined centrally in:

```text
config/config.sh
```

---

## Core analysis defaults

The final analysis used the following primary SCRMshaw settings:

| Parameter | Primary setting |
|---|---:|
| Training set | `adult_muscle` |
| Scoring method | `imm` |
| Number of offsets | 25 |
| Offset range | 0–240 bp |
| Offset step size | 10 bp |
| SCRMshaw `--thitw` | 10,000 |
| SCRMshaw `--step` | `123` |
| Retained hits per offset | 5,000 |
| Post-processing `-num` | 5,000 |
| Post-processing `-topN` | `Median` |
| TRF masking | enabled |

For the final 25-species dataset, this corresponds to:

```text
25 species × 25 offsets = 625 SCRMshaw offset jobs
```

For a different species list, the number of offset jobs is calculated automatically as:

```text
number of species × 25
```

The scoring method, SCRMshaw parameters, environment names, and compute resources are controlled through `config/config.sh`.

---

## Input data

### Species list

The species to be analyzed are defined in:

```text
species.txt
```

Each non-empty, non-comment line must contain:

```text
<slug> <species search name>
```

For example:

```text
d_affinis Drosophila affinis
d_pseudoobscura Drosophila pseudoobscura
```

The first field is the workflow slug used for directories and filenames. The remaining text is the species query used to resolve the corresponding files within the source archives.

Species slugs must be unique.

---

### Genome and annotation source data

The default download stage retrieves the genome and annotation archives associated with the configured Zenodo record.

The default record ID is:

```text
18453526
```

The workflow downloads:

```text
genomes.tar.gz
annotations.tar.gz
README.md
Species_summary_301Fly.xlsx
```

into the configured download directory.

Genome files are recognized using:

```text
.fa
.fasta
.fna
.fas
```

and annotation files using:

```text
.gff3
.gff
.gtf
```

Compressed `.gz` members are supported and are decompressed during preparation of the selected species dataset.

Optional MD5 checksums can be supplied for the genome and annotation archives.

---

## Species resolution and extracted inputs

The download script resolves each entry in `species.txt` against the contents of the genome and annotation archives.

Exact taxon-name matches are preferred. If an exact match is unavailable, candidate files are ranked using the supplied species query and slug.

A species is accepted only when both a genome and an annotation can be resolved unambiguously.

If resolution is ambiguous or unsuccessful, the workflow aborts and writes candidate matches to:

```text
data/resolution_report.txt
```

The species query in `species.txt` can then be made more specific before rerunning the pipeline.

Successfully resolved files are extracted only for the requested species and standardized as:

```text
data/selected/<slug>/genome.fa
data/selected/<slug>/annotation.gff3
```

or, when the source annotation is GTF:

```text
data/selected/<slug>/annotation.gtf
```

The resulting species manifest is written to:

```text
data/manifest.tsv
```

with columns:

```text
index
slug
species
genome
annotation
annotation_format
```

This manifest provides the authoritative species-level input description for all later SLURM stages.

---

## Software setup

### Conda environments

Separate environments are used for SCRMshaw execution and post-processing.

Create them using:

```bash
conda env create -f envs/scrmshaw.yml
conda env create -f envs/postprocess.yml
```

The corresponding environment names are configured through:

```text
SCRM_ENV
POST_ENV
```

in `config/config.sh`.

---

### Third-party software and training data

Required SCRMshaw software and training data are prepared with:

```bash
bash scripts/00_setup_software.sh
```

The setup script installs or updates:

- SCRMshaw-HD;
- HalfonLab UtilityPrograms;
- the HalfonLab SCRMshaw post-processing pipeline;
- the *D. melanogaster* SCRMshaw training-set repository;
- Tandem Repeat Finder.

The primary analysis uses:

```text
final_combined_48Tsets/adult_muscle
```

as the SCRMshaw training set.

Third-party repositories are stored below the configured `software/` directory and are excluded from version control.

Tandem Repeat Finder is downloaded as version:

```text
4.09.1
```

The setup step also verifies that the selected training set contains:

```text
crms.fasta
neg.fasta
```

and writes the resolved training-set path to the training-list file consumed by SCRMshaw.

> For exact long-term reproduction, the Git commit or release used for each cloned HalfonLab repository should be recorded or pinned. The setup script itself updates existing repositories to their current fast-forward state.

---

## Running the complete workflow

From:

```text
01_scrmshaw/generation/
```

run:

```bash
bash submit_pipeline.sh
```

Before submission, the wrapper verifies that the required SCRMshaw executables, utility scripts, post-processing pipeline, training data, and TRF binary are available.

The pipeline is then submitted as a sequence of dependent SLURM jobs:

```text
Download and extraction
        │
        ▼
Species preparation array
        │
        ▼
SCRMshaw offset array
        │
        ▼
Species post-processing array
        │
        ▼
Summary generation
```

Each downstream stage uses an `afterok` dependency and therefore begins only after the required upstream job has completed successfully.

---

## 1. Download and resolve selected species

The first SLURM stage executes the genome and annotation download/resolution workflow.

For each requested species, it:

1. resolves the genome and annotation archive members;
2. aborts on unresolved or ambiguous matches;
3. extracts only the selected files;
4. decompresses them into standardized species directories;
5. writes `data/manifest.tsv`.

Existing downloads are reused when possible. If an expected MD5 checksum is supplied and an existing archive does not match it, the file is downloaded again.

---

## 2. Prepare each species for SCRMshaw

Species preparation is executed as a SLURM array with one task per manifest entry.

### Annotation normalization

If the selected annotation is already GFF3, it is linked directly into the species run directory.

If the source annotation is GTF, it is converted to GFF3 before further processing.

### SCRMshaw preflight validation

The prepared genome and annotation are checked using the SCRMshaw preflight utility.

The preparation stage requires successful checks for:

- gene records in the annotation;
- exon records in the annotation;
- valid FASTA sequence characters;
- inclusion of all GFF sequence identifiers in the FASTA assembly.

A failed preflight check terminates preparation for that species.

### Genome filtering

The genome FASTA is filtered to retain only genomic sequences represented by the annotation.

The filtered assembly is written below the species-specific run directory and becomes the sequence input for repeat masking and SCRMshaw.

### Tandem Repeat Finder masking

When:

```text
RUN_TRF=1
```

the filtered genome is masked using Tandem Repeat Finder.

A usable `.mask` output is required before the species is considered prepared.

If TRF returns a non-zero exit code but still produces a valid non-empty mask file, the workflow records a warning and continues.

When TRF masking is disabled, the filtered but unmasked genome is used directly.

Successful preparation records the resolved genome and annotation paths and creates:

```text
PREPARED.ok
```

for the species.

---

## 3. Run SCRMshaw-HD offsets

For each species, the workflow executes 25 independent SCRMshaw-HD runs corresponding to lower-bound offsets:

```text
0, 10, 20, ..., 240 bp
```

The SLURM array maps each global task ID to one species and one offset.

Each offset run uses the prepared genome and GFF together with the configured training set and scoring flags.

The SCRMshaw invocation includes:

```text
--thitw <THITW>
--lb <offset>
--step 123
```

as well as the configured genome, GFF, training-list, and scoring-method arguments.

For the primary analysis:

```text
THITW = 10000
method = imm
```

The workflow verifies that the expected ranked-hit file was created for each enabled scoring method.

A successfully completed offset is marked with:

```text
SCRM_FINISHED.ok
```

inside its `task_offset_*` directory.

---

## 4. Post-process species-level predictions

Post-processing begins only after all 25 offset jobs for a species have completed successfully.

The workflow first verifies the presence of all 25 `SCRM_FINISHED.ok` markers.

For each offset, the top 5,000 ranked SCRMshaw hits are then extracted using:

```text
Generate_top_N_SCRMhits.pl -n 5000
```

The 25 resulting files are concatenated into:

```text
scrmshawOutput_offset_0to240.bed
```

and retained in the species run directory under:

```text
scrmsIndividualHits_0to240offset/
```

The combined offset predictions are passed to the HalfonLab post-processing pipeline using:

```text
-num 5000
-topN Median
```

together with the prepared species GFF.

Peak files generated by the post-processing workflow are concatenated into the final species-level result:

```text
peaks_AllSets.bed
```

Previous post-processing intermediates are removed before a rerun, while the original `task_offset_*` SCRMshaw directories are retained.

---

## 5. Generate summaries and quality-control outputs

For each species, the final result directory contains:

```text
results/<slug>/peaks_AllSets.bed
results/<slug>/scrmshawOutput_peaksCalled_*
results/<slug>/peak_count.txt
```

`peak_count.txt` records the number of lines in the final combined peak file.

The summary stage additionally creates:

```text
results/<slug>/top10_peaks.tsv
```

containing the ten highest-scoring final peaks for that species.

A cross-species summary is written to:

```text
results/summary.tsv
```

with:

```text
slug
species
n_peaks
peaks_file
```

---

## Key outputs

### Species manifest

```text
data/manifest.tsv
```

Defines the genome and annotation selected for each requested species and is used by all species-wise SLURM stages.

---

### Final species-level CRE predictions

```text
results/<slug>/peaks_AllSets.bed
```

This is the principal output of the generation workflow and the file subsequently consumed by the external-integration workflow.

---

### Individual post-processing peak sets

```text
results/<slug>/scrmshawOutput_peaksCalled_*
```

Contain the individual peak-call outputs that are concatenated to produce `peaks_AllSets.bed`.

---

### Species-level QC summaries

```text
results/<slug>/peak_count.txt
results/<slug>/top10_peaks.tsv
```

Provide the final peak count and a compact inspection set of the highest-scoring predictions.

---

### Cross-species summary

```text
results/summary.tsv
```

Reports the number and location of final SCRMshaw peaks for every species in the manifest.

---

## Additional peak-output QC

A separate QC script can be used after pipeline completion:

```bash
bash qc_peaks_outputs.sh results
```

For each species, it checks:

- whether `peaks_AllSets.bed` exists and is non-empty;
- whether individual `scrmshawOutput_peaksCalled_*` files exist;
- whether any individual peak file is empty;
- whether the summed number of lines across individual peak files equals the number of lines in `peaks_AllSets.bed`.

Possible QC states include:

```text
OK
NO_PEAKS_ALLSETS
NO_PEAK_FILES
EMPTY_PEAK_FILE
LINECOUNT_MISMATCH
```

This QC is separate from the automatically submitted five-stage SLURM workflow.

---

## SLURM configuration

SLURM resources are defined in:

```text
config/config.sh
```

For the final analysis, an individual SCRMshaw offset job used:

```text
CPUs      : 1
Memory    : 5 GB
Walltime  : 72 h
```

with a default maximum of:

```text
MAX_OFFSET_JOBS=10
```

concurrent offset jobs.

The submission wrapper also supports configured partition, QoS, and account settings.

The current project wrapper submits jobs to:

```text
abacus-2
```

via the SLURM `--nodelist` option.

Resource settings should be adapted to the local compute environment before running the workflow on another cluster.

---

## Restart and resume behavior

The pipeline is designed to reuse successfully completed work where possible.

### Downloads

Existing non-empty source downloads are reused. If MD5 validation is enabled, an archive is reused only when its checksum matches the expected value.

### Species preparation

Successful preparation is marked by:

```text
PREPARED.ok
```

and produces persistent resolved genome/GFF path files for downstream jobs.

### SCRMshaw offsets

When:

```text
RESUME=1
```

an offset containing a valid:

```text
SCRM_FINISHED.ok
```

marker is skipped.

If an unfinished offset must be executed, its existing offset directory is removed and rebuilt before SCRMshaw is run.

### Post-processing

Post-processing verifies that all 25 offsets are complete before proceeding.

Previous post-processing files are removed before regeneration, but the computationally expensive `task_offset_*` directories are retained.

---

## Configuration

Project-specific paths, biological settings, environment names, and SLURM resources are centralized in:

```text
config/config.sh
```

Important configurable values include:

- species-file location;
- data, run, result, log, and software directories;
- SCRMshaw and HalfonLab software paths;
- training-set path;
- SCRMshaw scoring flags;
- `THITW`;
- TRF masking;
- Conda environment names;
- resume behavior;
- SLURM resources and concurrency limits.

Routine adaptations should be made through `config/config.sh` rather than directly in the worker scripts.

---

## Adding or changing species

To generate SCRMshaw predictions for another species:

1. add a unique slug and sufficiently specific species query to `species.txt`;
2. rerun the pipeline;
3. inspect `data/resolution_report.txt` if genome or annotation resolution is ambiguous.

The number of species preparation, SCRMshaw offset, and post-processing tasks is calculated automatically from the contents of `species.txt`.

Changing the species list does not require modification of the SLURM worker scripts.

---

## Storage considerations

SCRMshaw generation is the most storage-intensive stage of the project.

Disk usage is driven primarily by:

- downloaded genome and annotation archives;
- extracted species genomes;
- prepared and masked genomes;
- 25 independent `task_offset_*` directories per species;
- offset-level prediction files;
- post-processing intermediates.

For the final 25-species analysis, 625 SCRMshaw offset working directories are generated.

Large source archives, third-party software, SCRMshaw working directories, and reproducible intermediate data are therefore excluded from version control.

Before running the complete workflow, ensure that sufficient server-side storage is available for the configured number of species.

---

## Reproducibility

The workflow separates source-data acquisition, species resolution, genome preparation, SCRMshaw execution, post-processing, and summary generation into explicit stages connected through manifests and completion markers.

The principal provenance files retained by the workflow include:

```text
data/manifest.tsv
data/resolution_report.txt
runs/<slug>/preflight.log
runs/<slug>/filter_fasta.log
runs/<slug>/PREPARED.ok
runs/<slug>/task_offset_*/SCRM_FINISHED.ok
results/<slug>/peak_count.txt
results/summary.tsv
```

These files allow the final species-level `peaks_AllSets.bed` results to be traced back to the selected genome and annotation, preparation status, completed SCRMshaw offsets, and post-processing output.

The standardized species-level results generated here are subsequently integrated with externally generated predictions by:

```text
01_scrmshaw/external/
```
