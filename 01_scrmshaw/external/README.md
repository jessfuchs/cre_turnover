# Integration of External SCRMshaw Predictions

This workflow imports previously generated SCRMshaw predictions and harmonizes them with the predictions generated de novo in [`../generation/`](../generation/).

External prediction files are not regenerated. Instead, records matching the retained SCRMshaw training set and scoring method are selected, validated, and subjected to the same downstream SCRMshaw post-processing procedure used for newly generated predictions. The resulting species-level peak sets are then combined into a common prediction collection for downstream comparative analyses.

The complete workflow is executed with:

```bash
bash run_external_pipeline.sh
```

---

## Workflow overview

The external integration workflow consists of five stages:

1. prepare genome annotations for the external species;
2. select the retained SCRMshaw training set and scoring method from the external BED files;
3. validate BED structure, rank organization, genome-coordinate compatibility, and species completeness;
4. apply the standardized SCRMshaw post-processing workflow;
5. combine the processed external results with the predictions generated in `01_scrmshaw/generation/`.

Conceptually:

```text
External SCRMshaw BEDs ──┐
                         ├─> filtering + QC ─> post-processing ─┐
External GFF3 files ─────┘                                      │
                                                                ├─> combined results
Generated SCRMshaw results ─────────────────────────────────────┘
```

> **Important:** This stage harmonizes existing external predictions at the filtering and post-processing level. It does not rerun the original SCRMshaw scans.

---

## Repository structure

```text
01_scrmshaw/external/
├── config/
│   └── external_config.sh
├── scripts/
│   ├── 01_extract_external_gffs.sh
│   ├── 02_filter_external_beds.sh
│   ├── 03_validate_bed_gff_and_manifest.py
│   ├── 04_postprocess_external_species.sh
│   ├── 05_build_combined_results.py
│   └── 06_run_combined_summary.sh
├── species_external.txt
└── run_external_pipeline.sh
```

Input and output locations are defined in:

```text
config/external_config.sh
```

The worker scripts should normally not require direct path modifications.

---

## Core analysis defaults

The external prediction files are required to match the SCRMshaw configuration retained for the primary analysis.

| Parameter | Primary setting |
|---|---:|
| Training set | `adult_muscle` |
| Scoring method | `imm` |
| Expected SCRMshaw offsets | 25 |
| Offset range represented by the post-processing input | 0–240 bp |
| Maximum accepted rank | 5,000 |
| Maximum theoretical BED rows | 125,000 |
| Required BED columns | 17 |
| Post-processing `-num` | 5,000 |
| Post-processing `-topN` | `Median` |

The filtering stage selects records by training set and scoring method. Rank depth is subsequently validated rather than truncated: files containing ranks above 5,000 fail QC.

The workflow does **not** require exactly 125,000 rows. Species with fewer predictions are accepted as long as the expected offset and rank structure remains valid.

---

## Input data

### External species list

External species are defined in:

```text
species_external.txt
```

Each non-comment, non-empty line must contain:

```text
<slug> <Genus> <species>
```

For example:

```text
d_affinis Drosophila affinis
```

Species slugs must be unique.

The slug is used throughout the workflow for file resolution and output directory naming.

---

### External SCRMshaw predictions

For each species, the workflow expects a non-empty raw SCRMshaw BED file at the configured external BED location with the naming pattern:

```text
<EXTERNAL_BED_DIR>/<slug>.bed
```

The raw BED may contain multiple SCRMshaw training sets and scoring methods.

The workflow selects records for the configured:

```text
TRAINING_SET
METHOD
```

using:

- column 15: training set
- column 16: scoring method
- column 17: rank

Every retained row must contain exactly 17 tab-separated fields.

---

### Genome annotations

For each external species, the source GFF3 annotation is expected at:

```text
<SOURCE_GFF_ROOT>/<slug>/annotation.gff3
```

Prepared annotations are stored using the standardized name:

```text
<EXTERNAL_GFF_DIR>/<slug>.gff3
```

A GFF is considered structurally usable when it is non-empty and contains at least one non-comment record with at least nine tab-separated fields.

Existing prepared GFF files are reused when they pass this basic validation.

---

### Generated SCRMshaw results

The final integration step additionally requires:

- the generated-species manifest;
- the generated species-level result directories from `01_scrmshaw/generation/`.

For every generated or external species, the corresponding result directory must contain a non-empty:

```text
peaks_AllSets.bed
```

---

## How to run the workflow

From:

```text
01_scrmshaw/external/
```

run:

```bash
bash run_external_pipeline.sh
```

The wrapper executes all five stages sequentially and terminates if any required input, QC check, species count, post-processing step, or final result is invalid.

---

## 1. Prepare external GFF annotations

The first stage processes every species listed in `species_external.txt`.

For each species, the workflow:

1. checks whether a valid prepared GFF already exists;
2. otherwise locates:

```text
<SOURCE_GFF_ROOT>/<slug>/annotation.gff3
```

3. copies it to:

```text
<EXTERNAL_GFF_DIR>/<slug>.gff3
```

4. validates the copied file.

Valid existing GFFs are reused, avoiding unnecessary copying.

---

## 2. Filter external SCRMshaw BED files

For every external species, the raw BED:

```text
<EXTERNAL_BED_DIR>/<slug>.bed
```

is filtered to records matching the configured training set and scoring method.

The standardized output naming pattern is:

```text
<FILTERED_BED_DIR>/<slug>.<training>.<method>.bed
```

For the primary analysis this corresponds to:

```text
<slug>.adult_muscle.imm.bed
```

The filtering step does not independently select the first 5,000 rows. Instead, the resulting file must satisfy the configured rank constraints during QC.

Existing filtered BED files are reused only when they still pass all required validation checks. Invalid existing files are rebuilt automatically.

---

## 3. Validate external inputs and construct the manifest

Before post-processing, each species is subjected to a second validation stage.

### Species-file validation

The workflow requires:

- one slug per species;
- a two-part scientific name (`Genus species`);
- no duplicate slugs.

### BED validation

Each filtered BED must:

- be non-empty;
- contain exactly 17 columns per data row;
- contain only the configured training set;
- contain only the configured scoring method;
- contain an integer rank in column 17;
- contain no rank below 1;
- have a minimum rank of 1;
- contain exactly 25 records with rank 1;
- contain no rank above 5,000;
- contain no more than 125,000 rows.

The requirement of 25 rank-1 entries is used to verify the expected 25-offset structure.

Fewer than 5,000 candidates for an individual offset are permitted. Consequently, fewer than 125,000 total rows are valid and are not treated as an error.

### BED/GFF coordinate compatibility

The sequence identifiers occurring in the filtered BED are compared with those present in the corresponding GFF3 file.

Every BED SeqID must also occur in the GFF. Any missing genomic sequence identifier causes validation to fail before post-processing.

### External manifest

Only successfully validated species are written to the external manifest.

The manifest contains:

```text
index
slug
species
bed
gff
```

The top-level wrapper additionally verifies that the number of species in the manifest exactly matches the number of species listed in `species_external.txt`.

---

## 4. Apply standardized SCRMshaw post-processing

Each validated external species is processed independently.

The filtered BED is linked into the species-specific run directory as:

```text
scrmshawOutput_offset_0to240.bed
```

The SCRMshaw post-processing pipeline is then executed as:

```text
postProcessingScrmshawPipeline.py
    -num 5000
    -topN Median
    -so scrmshawOutput_offset_0to240.bed
    -gff <species.gff3>
```

Before execution, previous post-processing outputs for that species are removed to prevent stale peak files from entering a rerun.

The generated:

```text
scrmshawOutput_peaksCalled_*
```

files are concatenated into:

```text
peaks_AllSets.bed
```

The final species result directory contains:

```text
peaks_AllSets.bed
scrmshawOutput_peaksCalled_*
peak_count.txt
```

`peak_count.txt` records the number of final peak records for that species.

---

## 5. Build the combined SCRMshaw result set

After all external species have been processed, the workflow combines the generated and external species manifests.

The combined manifest contains:

```text
index
slug
species
source
```

where:

```text
source = generated
```

or:

```text
source = external
```

Species slugs must be unique across both input manifests. A duplicate slug causes the integration step to terminate.

The combined results directory is rebuilt from scratch on every run.

For each species, the workflow creates:

```text
combined_results/<slug>/peaks_AllSets.bed
```

as a symbolic link to the corresponding final peak file in either the generated or external result directory.

This avoids duplicating the underlying peak data while providing a uniform directory structure for all downstream workflow stages.

The workflow fails if any expected species-level `peaks_AllSets.bed` file is missing or empty.

Finally, the combined manifest and combined results directory are passed to the standard SCRMshaw summary script.

---

## Key outputs

### Filtered external predictions

```text
<FILTERED_BED_DIR>/<slug>.<training>.<method>.bed
```

Contains the external SCRMshaw records retained for the selected analysis configuration.

---

### External species manifest

The configured external manifest records the validated input BED and GFF path for every external species:

```text
index    slug    species    bed    gff
```

This manifest is the authoritative input to external post-processing.

---

### External processed predictions

```text
<EXTERNAL_RESULTS_DIR>/<slug>/peaks_AllSets.bed
```

Contains the final post-processed CRE prediction set for each external species.

---

### Combined manifest

```text
combined_manifest.tsv
```

Defines the complete prediction framework and records whether each species originates from the generated or external prediction set.

---

### Combined species results

```text
combined_results/<slug>/peaks_AllSets.bed
```

Provides the common species-level prediction structure consumed by downstream analyses.

These files are symbolic links to the corresponding generated or external post-processing results.

---

### Combined summary

After the combined result structure has been built, the standard SCRMshaw summary workflow is executed using the combined manifest and results directory.

The exact summary output is produced by the shared `make_summary.py` implementation used by the SCRMshaw workflow.

---

## Restart and reuse behavior

The workflow is designed to distinguish reusable validated intermediates from outputs that should be regenerated.

| Stage | Behavior on rerun |
|---|---|
| Prepared GFF | Reused if valid |
| Filtered BED | Reused only if it passes QC; otherwise rebuilt |
| Species post-processing | Previous species-specific post-processing outputs are removed and regenerated |
| Combined results | Entire directory is removed and rebuilt |
| Combined peak files | Recreated as symlinks to the current generated/external results |

This makes reruns reproducible while avoiding unnecessary repetition of already validated preparation steps.

---

## Configuration

Paths and methodological settings are centralized in:

```text
config/external_config.sh
```

The configuration provides the locations and settings required by the workflow, including:

- external raw BED directory;
- source and prepared GFF directories;
- filtered BED directory;
- external run and result directories;
- generated manifest and result directory;
- external and combined manifests;
- combined result directory;
- retained SCRMshaw training set;
- retained scoring method;
- post-processing environment and SCRMshaw post-processing installation.

Routine path or parameter changes should be made in the configuration file rather than directly in the worker scripts.

---

## Adding an external species

To add another species to the external prediction set:

1. add a unique slug and scientific name to `species_external.txt`;
2. provide the raw external SCRMshaw BED as:

```text
<EXTERNAL_BED_DIR>/<slug>.bed
```

3. provide the source GFF as:

```text
<SOURCE_GFF_ROOT>/<slug>/annotation.gff3
```

4. rerun:

```bash
bash run_external_pipeline.sh
```

The workflow will validate the new species before including it in the combined result set.

---

## Relation to the main SCRMshaw workflow

Predictions generated directly within this project are documented in:

```text
01_scrmshaw/generation/
```

This external workflow begins from already existing SCRMshaw prediction files and harmonizes them for use alongside those generated datasets.

After completion, downstream workflow stages should consume the unified species-level prediction collection from:

```text
combined_results/
```

rather than accessing generated and external prediction directories separately.
