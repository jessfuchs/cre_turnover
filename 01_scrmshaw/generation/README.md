# Multi-species SCRMshaw-HD Prediction Pipeline

This workflow generates and post-processes SCRMshaw-HD cis-regulatory element (CRE) predictions for 25 Drosophilidae species using a SLURM-based high-performance computing environment.

The pipeline was designed to standardize genome preparation, SCRMshaw execution, and peak calling across all newly analyzed species so that the resulting CRE predictions can subsequently be integrated with externally generated SCRMshaw predictions and used in comparative analyses.

## Overview

For each species listed in `species.txt`, the workflow performs the following steps:

1. Downloads the genome and annotation archives associated with the Zenodo dataset accompanying *Comparative gene annotation and orthology assignments across 301 species of Drosophilidae*.
2. Extracts only the genome FASTA and GFF3 annotation files required for the selected species.
3. Validates genome–annotation compatibility and prepares the genome for SCRMshaw.
4. Removes genomic sequences that are not represented in the corresponding annotation and applies Tandem Repeat Finder masking.
5. Executes 25 independent SCRMshaw-HD offsets per species.
6. Applies the HalfonLab SCRMshaw post-processing workflow to the combined predictions.
7. Generates standardized species-level peak sets and summary statistics.

## Configuration

The main SCRMshaw settings are defined in:

```text
config/config.sh
```

The final analysis uses:

```text
Training set : adult_muscle
Method       : IMM
Offsets      : 25 per species
THITW        : 10000
TRF masking  : enabled
```

## 1. Create the Conda environments

The SCRMshaw and post-processing environments can be created from the supplied environment definitions:

```bash
conda env create -f envs/scrmshaw.yml
conda env create -f envs/postprocess.yml
```

The default environment names are defined in `config/config.sh` as:

```bash
SCRM_ENV="scrmshaw"
POST_ENV="scrm_postproc"
```

These variables can be modified if different environment names are used.

## 2. Install the required software

Required third-party software and training data are installed using:

```bash
bash scripts/setup_software.sh
```

The setup includes the components required for the SCRMshaw workflow, including:

* SCRMshaw-HD,
* associated utility programs,
* the HalfonLab post-processing workflow,
* D. melanogaster SCRMshaw training sets,
* Tandem Repeat Finder.

The CRE predictions generated in this study use:

```text
final_combined_48Tsets/adult_muscle
```

as the SCRMshaw training set.

Third-party software is installed under:

```text
software/
```

and is not included in the Git repository.

## 3. Define the target species

The 25 newly analyzed species are specified in:

```text
species.txt
```

Each line contains a short species identifier followed by the species name or another sufficiently specific identifier that allows the corresponding genome and annotation files to be resolved within the Zenodo archives.

For example:

```text
dazt    Drosophila azteca
dpse    Drosophila pseudoananassae
```

If a species identifier does not resolve uniquely, the download and extraction step terminates rather than selecting an ambiguous file. Candidate matches are recorded in:

```text
data/resolution_report.txt
```

The corresponding entry in `species.txt` can then be refined before rerunning the pipeline. 

## 4. Genome and annotation preparation

Genome FASTA and GFF3 files are extracted from the source archives only for species included in `species.txt`.

For each species, the preparation step validates consistency between the genome assembly and genome annotation. Sequence identifiers occurring in the GFF3 file are checked against the FASTA assembly, and genomic sequences without corresponding annotation are removed from the SCRMshaw input assembly.

Tandem Repeat Finder masking is enabled by default:

```bash
RUN_TRF=1
```

This provides an additional repeat-masking step on top of the soft masking already present in the source genome assemblies.

Prepared species datasets are stored below:

```text
data/selected/
```

A species manifest describing the resolved genome and annotation files is generated as:

```text
data/manifest.tsv
```

## 5. SCRMshaw-HD execution

SCRMshaw is executed independently for each species using 25 offsets.

Each newly created `task_offset_*` working directory is initialized with:

```text
--step 123
```

because a newly created output directory does not yet contain all files required by the later SCRMshaw stages, including the `gff/genes` structures and offset-specific genomic windows.

For the final analysis, the workflow therefore executes:

```text
25 species × 25 offsets = 625 SCRMshaw offset jobs
```

using the IMM scoring method.

## 6. SLURM execution

Cluster-resource settings are defined in:

```text
config/config.sh
```

The current default resources for an individual SCRMshaw offset job are:

```text
CPU      : 1
Memory   : 5 GB
Walltime : 72 h
```

At most 10 offset jobs are submitted concurrently by default:

```bash
MAX_OFFSET_JOBS=10
```

These limits can be changed in `config/config.sh` if required by the compute environment.

## 7. Run the complete pipeline

From the `01_scrmshaw/generation/` directory, execute:

```bash
bash submit_pipeline.sh
```

The wrapper submits the individual SLURM stages and automatically establishes the required job dependencies.

For the final 25-species dataset, the workflow structure is:

```text
Download and extraction
        ↓
25 × species preparation
        ↓
625 × SCRMshaw offset jobs
        ↓
25 × species post-processing
        ↓
Summary generation
```

Later stages are started only after their required upstream jobs have completed successfully.

## 8. Post-processing

After all SCRMshaw offsets for a species have completed, their predictions are processed using the HalfonLab SCRMshaw post-processing workflow.

The post-processing step combines offset-level predictions and performs peak calling to generate a standardized CRE prediction set for each species.

The principal species-level output is:

```text
results/<species>/peaks_AllSets.bed
```

Additional locally generated quality-control outputs include:

```text
results/<species>/peak_count.txt
results/<species>/top10_peaks.tsv
```

A cross-species summary is generated as:

```text
results/summary.tsv
```

The standardized `peaks_AllSets.bed` files are subsequently combined with the equivalently processed external SCRMshaw predictions by the workflow in:

```text
../external/
```

## Restart and resume behavior

The pipeline is designed to reuse completed work where possible.

Completed SCRMshaw offset jobs are identified through:

```text
SCRM_FINISHED.ok
```

and can therefore be skipped during subsequent executions.

Resume behavior is enabled by default:

```bash
RESUME=1
```

Existing downloads and completed outputs are reused where supported by the corresponding pipeline step.

The post-processing stage removes only its own temporary intermediate files and does not delete the original `task_offset_*` SCRMshaw working directories.

## Storage considerations

The source genome and annotation archives occupy substantial disk space, and the complete SCRMshaw working directories can considerably exceed the size of the compressed input data.

Storage requirements depend on genome size, the number of species, the number of offsets, and the enabled scoring methods. Because 25 independent offsets are generated for each of 25 species, sufficient server-side storage should be available before executing the complete workflow.

Large source files, SCRMshaw working directories, third-party software, and other reproducible intermediate data are intentionally excluded from version control.
