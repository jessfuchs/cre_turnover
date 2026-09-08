# D. melanogaster Ortholog Mapping and Reference CRE Construction

This workflow links species-specific SCRMshaw CRE predictions to their corresponding *Drosophila melanogaster* target-gene identifiers and constructs the fixed *D. melanogaster* reference CRE set used in the downstream comparative analysis.

Ortholog assignments are derived directly from the `dmel_orthologs` attributes contained in the species-specific GFF3 annotations. These identifiers are subsequently normalized to standard FlyBase FBgn identifiers using the *D. melanogaster* GFF3 annotation.

No separate external orthology database is required.

The final outputs are a cross-species FBgn-annotated SCRMshaw dataset and the annotated *D. melanogaster* reference CRE set in TSV and BED format.

---

## Workflow overview

The workflow performs six sequential steps:

1. validate the availability of *D. melanogaster* ortholog annotations;
2. map genes associated with SCRMshaw peaks to *D. melanogaster* ortholog identifiers;
3. combine species-specific mappings into a single cross-species table;
4. normalize *D. melanogaster* identifiers to FlyBase FBgn IDs;
5. investigate identifiers that could not be resolved to FBgn;
6. construct the fixed *D. melanogaster* reference CRE set.

Conceptually:

```text
Standardized SCRMshaw predictions
            +
species-specific GFF3 annotations
            │
            ▼
Peak-associated genes
    → Dmel ortholog identifiers
            │
            ▼
Combined cross-species table
            │
            ▼
Dmel identifiers → FlyBase FBgn
            │
            ├── unresolved-ID QC
            │
            ▼
D. melanogaster CRE subset
            │
            ▼
Reference CRE TSV + BED
```

The complete workflow is executed with:

```bash
bash run_ortholog_pipeline.sh
```

---

## Core analysis defaults

The final analysis uses:

| Parameter | Setting |
|---|---:|
| Species framework | 40 species |
| Reference species | *D. melanogaster* |
| Reference-species slug | `d_melanogaster` |
| Expected reference CREs | 337 |
| Minimum fraction of mRNAs with Dmel ortholog annotation | 0.10 |
| Ortholog source | GFF3 `dmel_orthologs` attribute |
| Standard target-gene identifier | FlyBase FBgn |

The 0.10 ortholog-annotation threshold is a conservative technical QC criterion and not a biological filtering threshold.

The expected reference-set size of 337 is also used as a consistency check. A different number produces a warning rather than terminating the workflow.

---

## Input data

### Standardized SCRMshaw predictions

The workflow uses the unified SCRMshaw prediction framework generated in `01_scrmshaw/`:

```text
01_scrmshaw/external/combined_manifest.tsv
01_scrmshaw/external/combined_results/
```

The manifest identifies each species and whether its SCRMshaw predictions were generated within the project or supplied externally.

For each species, the workflow expects:

```text
combined_results/<species>/peaks_AllSets.bed
```

---

### Species-specific genome annotations

The corresponding GFF3 annotations are obtained from the configured generated or external SCRMshaw data directories.

For non-reference species, the annotations must contain usable:

```text
dmel_orthologs=
```

attributes on `mRNA` features.

These annotations provide the connection between species-specific genes associated with SCRMshaw peaks and their *D. melanogaster* ortholog identifiers.

Multiple ortholog assignments are retained when present.

---

### D. melanogaster annotation

A *D. melanogaster* GFF3 annotation is additionally required to normalize the different Dmel identifier formats to FlyBase FBgn identifiers.

The path is configured through:

```text
DMEL_GFF
```

in:

```text
config/ortholog_config.sh
```

---

## Configuration

All paths and primary QC parameters are defined in:

```text
config/ortholog_config.sh
```

Important configurable values include:

```text
COMBINED_MANIFEST
COMBINED_RESULTS_DIR

GENERATED_GFF_ROOT
EXTERNAL_GFF_ROOT
DMEL_GFF

RESULTS_DIR
REFERENCE_CRES_DIR

EXPECTED_SPECIES
EXPECTED_REFERENCE_CRES
MIN_ORTHOLOG_FRACTION
DMEL_SLUG
```

Routine path or dataset changes should be made in this configuration file rather than directly in the worker scripts.

---

## Run the complete workflow

From:

```text
02_mapping_orthologs/
```

execute:

```bash
bash run_ortholog_pipeline.sh
```

The wrapper performs all six stages sequentially and stops if a required input or mandatory QC step fails.

---

## 1. Validate Dmel ortholog annotations

Before mapping begins, the workflow checks whether each species-specific GFF3 annotation contains sufficient *D. melanogaster* ortholog information.

For non-reference species, the QC verifies:

- presence of a non-empty GFF3 file;
- presence of `mRNA` features;
- presence of the `dmel_orthologs` attribute;
- at least one usable ortholog assignment;
- an ortholog-annotated mRNA fraction of at least 0.10.

*D. melanogaster* is treated separately because the reference species does not require an orthology assignment to itself.

The QC report is written to:

```text
results/dmel_ortholog_annotation_qc.tsv
```

The workflow also verifies that the expected number of species is represented.

---

## 2. Map SCRMshaw-associated genes to Dmel orthologs

For each non-reference species, the workflow builds a mapping between species-specific genes and their *D. melanogaster* ortholog identifiers using the GFF3 `Parent` and `dmel_orthologs` attributes.

The two gene associations stored for each SCRMshaw prediction are mapped independently:

```text
flanking gene
next flanking gene
```

If no Dmel ortholog can be assigned, the corresponding value is recorded as:

```text
NA
```

For *D. melanogaster*, the original gene identifiers are retained directly.

Species-specific mapped predictions are written to:

```text
results/<species>/SO_all_peaks.bed
```

Mapping completeness is summarized in:

```text
results/ortholog_mapping_qc.tsv
```

---

## 3. Combine species-specific mappings

The species-specific `SO_all_peaks.bed` files are combined into:

```text
results/SO_all_species.tsv
```

The resulting table retains the original SCRMshaw information together with:

- species identity;
- prediction source;
- species-specific associated genes;
- corresponding *D. melanogaster* ortholog identifiers.

The workflow verifies that all expected species are represented in the combined dataset.

---

## 4. Normalize Dmel identifiers to FlyBase FBgn

The Dmel ortholog identifiers are converted to standardized FlyBase FBgn identifiers using the *D. melanogaster* GFF3 annotation.

The mapping uses identifier information available in the annotation, including:

```text
ID
transcript_id
locus_tag
dbxref
```

The two FBgn fields are appended independently for the flanking and next-flanking gene associations.

The resulting cross-species dataset is written to:

```text
results/SO_all_species_fbgn.tsv
```

Identifiers that cannot be converted are retained separately in:

```text
results/unresolved_dmel_identifiers.tsv
```

rather than being silently discarded.

---

## 5. QC unresolved Dmel identifiers

Unresolved Dmel identifiers are investigated separately to distinguish genuine mapping failures from identifier-format differences.

In particular, identifiers of the form:

```text
Dmelanogast_M...
```

are tested after normalization to:

```text
Dmelanogast_...
```

against the *D. melanogaster* GFF3 annotation.

The summary is written to:

```text
results/unresolved_dmel_qc_summary.tsv
```

with additional diagnostic files stored under:

```text
results/unresolved_dmel_qc/
```

This step is diagnostic and does not itself remove CRE predictions from the dataset.

---

## 6. Construct the D. melanogaster reference CRE set

The final stage extracts the *D. melanogaster* SCRMshaw predictions from the FBgn-normalized cross-species table.

FBgn assignments from both associated gene fields are combined for each CRE, with duplicate target-gene identifiers removed.

The CREs are sorted deterministically and assigned stable identifiers:

```text
DMEL_CRE_00001
DMEL_CRE_00002
...
```

Two reference files are generated:

```text
reference_cres/dmel_reference_cres.tsv
reference_cres/dmel_reference_cres.bed
```

The TSV retains CRE coordinates, SCRMshaw scores, prediction metadata, and associated FBgn target genes.

The BED6 file contains the genomic reference intervals and stable CRE identifiers required for downstream coordinate projection with `liftOver`.

For the primary analysis, the resulting reference set contains:

```text
337 CREs
```

The workflow additionally verifies that all reference CRE identifiers are unique.

---

## Key outputs

### Cross-species Dmel ortholog table

```text
results/SO_all_species.tsv
```

Contains the combined species-specific SCRMshaw predictions with associated *D. melanogaster* ortholog identifiers.

### FBgn-normalized prediction table

```text
results/SO_all_species_fbgn.tsv
```

Contains the standardized FlyBase FBgn target-gene assignments used by downstream comparative analyses.

### Reference CRE annotation

```text
reference_cres/dmel_reference_cres.tsv
```

Contains the fixed *D. melanogaster* reference CREs together with their SCRMshaw metadata and associated FBgn target genes.

### Reference CRE coordinates

```text
reference_cres/dmel_reference_cres.bed
```

BED6 representation of the reference CRE set used by the downstream pairwise whole-genome alignment and `liftOver` workflow.

---

## Quality control

QC is performed at several points in the workflow:

| Stage | Main check |
|---|---|
| GFF annotations | usable `dmel_orthologs` information |
| Species mapping | one non-empty mapped peak file per species |
| Combined table | all expected species represented |
| FBgn normalization | mapping rates and unresolved identifiers recorded |
| Reference CRE set | unique CRE IDs and expected reference-set size |

Mapping rates and unresolved identifiers are retained as diagnostic information and are not interpreted as biological evidence of CRE loss or turnover.

---

## Adapting the workflow

To analyze a different species framework, update the upstream SCRMshaw manifest and:

```text
EXPECTED_SPECIES
```

in `config/ortholog_config.sh`.

If the reference CRE definition changes, update:

```text
EXPECTED_REFERENCE_CRES
```

as appropriate.

The technical ortholog-annotation QC threshold can be changed through:

```text
MIN_ORTHOLOG_FRACTION
```

and an alternative *D. melanogaster* annotation can be specified through:

```text
DMEL_GFF
```

Because the workflow does not remove the complete `results/` directory before rerunning, obsolete species-specific result directories should be removed when changing the species set.

---

## Reproducibility

The workflow retains the intermediate ortholog-mapping tables, FBgn-normalized outputs, unresolved identifier reports, and QC summaries required to trace the final target-gene assignments back to the underlying SCRMshaw predictions and species-specific genome annotations.

The final *D. melanogaster* reference CRE set is generated deterministically from the FBgn-normalized table using a fixed reference-species identifier and stable sorting criteria, allowing the same reference CRE IDs and associated target-gene annotations to be reproduced from the same upstream inputs and configuration.
