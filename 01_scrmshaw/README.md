# SCRMshaw Prediction Workflows

This directory contains the SCRMshaw workflows used to obtain standardized species-level cis-regulatory element (CRE) predictions for the comparative analysis.

Two complementary routes are available:

1. **de novo generation** of SCRMshaw-HD predictions;
2. **integration of existing SCRMshaw predictions**.

Which route is required depends on the availability of prediction data for the species being analyzed.

<p align="center">
  <img src="docs/scrmshaw_workflow.svg"
       alt="Overview of the SCRMshaw prediction workflows"
       width="750">
</p>

---

## Which workflow should I use?

### All SCRMshaw predictions need to be generated

If no suitable SCRMshaw predictions are available for the species in the analysis, use only:

[`generation/`](generation/)

This workflow downloads and prepares the required genome resources, runs SCRMshaw-HD, and generates standardized species-level:

```text
peaks_AllSets.bed
```

files.

In this case, the `external/` workflow is not required.

---

### SCRMshaw predictions already exist

If suitable SCRMshaw prediction files already exist for the species being analyzed, use:

[`external/`](external/)

The external workflow filters and validates the existing prediction files and applies the standardized SCRMshaw post-processing procedure required by the downstream analysis.

If **all** analyzed species are supplied as existing external predictions, the final step that combines external and newly generated results is unnecessary. The processed external species-level results can instead be used directly as the SCRMshaw prediction input for the subsequent workflow stages.

---

### Some predictions exist and others need to be generated

If the analysis contains both:

- species for which SCRMshaw predictions must be generated de novo, and
- species for which suitable predictions already exist,

run both workflows:

```text
generation/
external/
```

The generated species are processed first, after which the external workflow harmonizes the existing predictions and builds the combined species-level prediction collection.

This is the configuration used for the final analysis in this project.

---

## Workflow comparison

| Workflow | Use when | Main function | Principal output |
|---|---|---|---|
| [`generation/`](generation/) | SCRMshaw predictions must be generated from genome data | Genome preparation, SCRMshaw-HD execution, and post-processing | `results/<species>/peaks_AllSets.bed` |
| [`external/`](external/) | SCRMshaw predictions already exist | Filtering, validation, and standardized post-processing of existing predictions | processed external `peaks_AllSets.bed` files |
| `generation/` + `external/` | Generated and existing predictions are both used | Harmonization and integration of both prediction sources | `external/combined_results/<species>/peaks_AllSets.bed` |

---

## `generation/`

The de novo generation workflow creates SCRMshaw-HD predictions from genome FASTA and annotation data.

It includes:

- genome and annotation acquisition;
- species resolution and input validation;
- genome filtering and repeat masking;
- 25 SCRMshaw-HD offsets per species;
- extraction of ranked SCRMshaw hits;
- standardized post-processing and peak calling;
- species-level and cross-species QC summaries.

The principal output is:

```text
generation/results/<species>/peaks_AllSets.bed
```

See [`generation/README.md`](generation/README.md) for detailed input requirements, SCRMshaw parameters, software setup, SLURM execution, QC, and restart behavior.

---

## `external/`

The external integration workflow starts from previously generated SCRMshaw prediction files rather than rerunning SCRMshaw-HD.

It:

- prepares the corresponding genome annotations;
- selects the retained training set and scoring method;
- validates BED structure, offset/rank organization, and BED/GFF coordinate compatibility;
- applies the same downstream SCRMshaw post-processing procedure used for newly generated predictions.

The principal processed result for each external species is:

```text
external/<external-results-directory>/<species>/peaks_AllSets.bed
```

When both generated and external prediction sources are used, the workflow additionally creates:

```text
external/combined_manifest.tsv
external/combined_results/<species>/peaks_AllSets.bed
```

The combined results provide a uniform species-level input structure for downstream analysis.

See [`external/README.md`](external/README.md) for detailed input formats, validation criteria, post-processing, integration logic, and restart behavior.

---

## Common output requirement

Regardless of how the predictions are obtained, downstream stages require one standardized:

```text
peaks_AllSets.bed
```

file per species.

The subsequent workflow should therefore be pointed to one consistent species-level prediction collection:

- `generation/results/` when all predictions were generated de novo;
- the processed external result directory when all predictions were supplied externally;
- `external/combined_results/` when generated and external prediction sources are combined.

The exact input location used by downstream stages is defined by the corresponding configuration files and manifests.

---

## Primary analysis configuration

The final analysis used:

| Parameter | Setting |
|---|---|
| SCRMshaw training set | `adult_muscle` |
| SCRMshaw scoring method | `imm` |
| Species framework | 40 species |

Detailed SCRMshaw execution parameters and external-data validation rules are documented in the corresponding subworkflow READMEs.

---

## Recommended workflow

In summary:

```text
Need to generate every species?
        │
        ├── yes ──> generation/ only
        │
        └── no
             │
             ├── all predictions already exist
             │       └──> external/ only
             │            (skip combined-results step)
             │
             └── mixture of generated and existing predictions
                     └──> generation/ + external/
                          └──> combined_results/
```

This separation allows the project to reuse existing SCRMshaw predictions where available without forcing unnecessary regeneration, while still producing a standardized prediction set for downstream comparative analyses.
