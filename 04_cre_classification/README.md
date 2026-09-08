# CRE-state Classification

This workflow classifies the fixed *Drosophila melanogaster* reference CREs across all target species by combining homologous genomic positions from pairwise `liftOver`, species-specific SCRMshaw predictions, and standardized FlyBase FBgn target-gene associations.

Each reference CRE–species comparison is assigned one of four operational states:

```text
positional_match
turnover_candidate
no_detected_CRE
uncertain
```

These states provide the basis for the downstream comparative, phylogenetic, and climate-related analyses.

`turnover_candidate` and `no_detected_CRE` are deliberately conservative computational labels. They identify patterns consistent with regulatory turnover or lack of a detected CRE, but do not by themselves demonstrate functional turnover or biological CRE loss.

---

## Workflow overview

For every *D. melanogaster* reference CRE and target species, the workflow first determines whether the homologous reference interval could be projected into the target genome.

Mapped intervals are compared with the target-species SCRMshaw predictions. If no positional match is detected, the workflow searches for an alternative prediction associated with the same *D. melanogaster* FBgn target gene within the defined local regulatory neighborhood.

```text
D. melanogaster reference CRE
            │
            ▼
   Homologous interval mapped?
        │               │
       no              yes
        │               │
        ▼               ▼
    uncertain    Positional SCRMshaw
                       match?
                   │           │
                  yes          no
                   │           │
                   ▼           ▼
            positional_    Local prediction
              match        with same FBgn?
                            │           │
                           yes          no
                            │           │
                            ▼           ▼
                      turnover_     no_detected_
                      candidate         CRE
```

The complete core classification workflow is executed with:

```bash
bash run_cre_class_pipeline.sh
```

---

## Directory layout

```text
04_cre_classification/
├── config/
│   └── classification_config.sh
├── scripts/
│   ├── 01_classify_target_cres.py
│   ├── 02_combine_turnover_results.py
│   ├── 03_build_species_qc_summary.py
│   ├── summarize_dmel_cre_gene_distances.py
│   ├── summarize_prediction_source.py
│   ├── plot_cre_turnover.py
│   ├── plot_tree_heatmap.py
│   └── plot_clade_heatmaps.py
├── phylogeny/
│   ├── data/
│   └── results/
├── results/
│   ├── turnover_by_species/
│   └── figures/
└── run_cre_class_pipeline.sh
```

Core paths, classification thresholds, expected dataset dimensions, and species-level QC parameters are defined centrally in:

```text
config/classification_config.sh
```

---

## Core analysis defaults

The primary analysis uses:

| Parameter | Setting |
|---|---:|
| Reference CREs | 337 |
| Target species | 39 |
| Positional match | reciprocal overlap ≥ 0.50 |
| Local same-FBgn distance | ≤ 24 kb |
| Minimum mapping rate for QC | 0.50 |
| Minimum SCRMshaw peaks for QC | 50 |
| Minimum sequence-ID overlap for QC | 0.25 |
| Minimum mapped CREs for positional-overlap QC | 100 |
| High mapping-rate threshold for QC | 0.90 |

A positional CRE match requires reciprocal overlap: at least 50% of both the lifted *D. melanogaster* CRE interval and the overlapping target-species SCRMshaw prediction must be shared.

The 24 kb local-distance threshold was derived empirically from the distribution of *D. melanogaster* SCRMshaw CRE-to-associated-gene distances. The 95th percentile was approximately 23.6 kb and was rounded to 24 kb for the primary classification.

---

## Input data

### D. melanogaster reference CREs

The fixed reference CRE set is obtained from:

```text
02_mapping_orthologs/reference_cres/dmel_reference_cres.tsv
```

Each reference CRE has a stable:

```text
dmel_cre_id
```

together with its genomic coordinates and associated FBgn target gene or genes.

---

### Homologous target-genome intervals

Mapped reference CRE coordinates are obtained from the preceding pairwise whole-genome alignment workflow:

```text
03_pairwise_wga/lifted_cres_dmel/<species>/
    dmel_reference_cres.<species>.bed
```

The target-species list is taken from:

```text
03_pairwise_wga/target_species.txt
```

A reference CRE absent from the mapped BED cannot be evaluated positionally in that species and is classified as `uncertain`.

---

### SCRMshaw predictions and FBgn associations

Target-species SCRMshaw predictions and their standardized *D. melanogaster* FBgn associations are read from:

```text
02_mapping_orthologs/results/SO_all_species_fbgn.tsv
```

The table provides:

- target-species SCRMshaw intervals;
- flanking and next-flanking gene associations;
- standardized FBgn identifiers;
- CRE-to-associated-gene distances;
- SCRMshaw prediction metadata.

These data are used both for positional CRE matching and for identifying alternative same-FBgn predictions.

---

## Configuration

The main classification settings are defined in:

```text
config/classification_config.sh
```

Important variables include:

```text
REFERENCE_CRES_TSV
SO_ALL_SPECIES_FBGN

LIFTED_CRES_DIR
TARGET_SPECIES_FILE

RESULTS_DIR
TURNOVER_BY_SPECIES_DIR

RECIPROCAL_OVERLAP
LOCAL_GENE_DISTANCE

EXPECTED_REFERENCE_CRES

QC_MIN_MAPPING_RATE
QC_MIN_SCRMSHAW_PEAKS
QC_MIN_SEQID_OVERLAP
QC_MIN_MAPPED_FOR_OVERLAP_CHECK
QC_HIGH_MAPPING_RATE
```

Routine changes to thresholds or input locations should be made through the configuration file rather than directly in the classification scripts.

---

## Run the complete classification workflow

From:

```text
04_cre_classification/
```

execute:

```bash
bash run_cre_class_pipeline.sh
```

The wrapper performs four sequential stages:

1. validate the required inputs and reference CRE count;
2. classify all reference CREs independently for every target species;
3. combine the species-specific classifications;
4. construct a species-level QC summary.

For the primary dataset, every target-species output must contain exactly 337 reference CRE rows.

---

## 1. Classify reference CREs in each target species

Each reference CRE is first checked for a successfully lifted homologous interval.

If the interval is available, all SCRMshaw predictions on the corresponding target sequence are evaluated for positional overlap. The best overlapping prediction is retained together with the overlap in base pairs and the overlap fraction relative to both intervals.

A positional match requires:

```text
overlap / lifted CRE length  >= 0.50
AND
overlap / target peak length >= 0.50
```

This reciprocal criterion avoids defining positional conservation from a large asymmetric overlap alone.

---

## CRE-state definitions

Each reference CRE–species comparison is assigned exactly one state:

| State | Operational definition |
|---|---|
| `positional_match` | The homologous interval was successfully mapped and a target-species SCRMshaw prediction satisfies the reciprocal overlap criterion. |
| `turnover_candidate` | The homologous interval was mapped but no positional match was detected, while another SCRMshaw prediction associated with the same FBgn target gene occurs within 24 kb of that associated gene. |
| `no_detected_CRE` | The homologous interval was mapped, but neither a positional CRE match nor a qualifying local same-FBgn prediction was detected. |
| `uncertain` | The reference CRE could not be successfully projected into the target genome and is therefore not positionally evaluable. |

### Positional matches and target-gene support

FBgn agreement between the reference CRE and the best positional target peak is retained separately as:

```text
shared_fbgn
discordant_fbgn
no_fbgn_information
```

This allows positional conservation and target-gene agreement to be evaluated independently rather than forcing them into a single classification criterion.

---

## Local same-FBgn turnover candidates

When no positional match is detected, the workflow searches the target-species predictions for CREs associated with any FBgn target gene assigned to the reference CRE.

For a prediction to support `turnover_candidate`, it must:

```text
share a D. melanogaster FBgn target gene
AND
have a CRE-to-associated-gene distance <= 24 kb
```

The closest qualifying same-FBgn prediction is retained for reporting.

The workflow additionally records whether a same-FBgn prediction exists anywhere in the target dataset, even if it does not satisfy the local-distance criterion. This information is diagnostic and does not by itself define a turnover candidate.

The 24 kb threshold can be independently summarized from the *D. melanogaster* data using:

```text
scripts/summarize_dmel_cre_gene_distances.py
```

which reports CRE-to-gene distance percentiles in:

```text
results/dmel_cre_gene_distance_summary.tsv
```

---

## 2. Species-specific classification outputs

For each target species, the complete reference CRE set is written to:

```text
results/turnover_by_species/
    dmel_to_<species>_cre_turnover.tsv
```

Each row retains the final state together with diagnostic information such as:

```text
alignment status
lifted coordinates
best positional peak
reciprocal overlap fractions
same-FBgn support
local same-FBgn candidate
CRE-to-gene distance
reference target FBgn
gene support at the positional peak
```

This preserves the evidence underlying every final state assignment.

---

## 3. Combine classifications across species

The species-specific results are combined into three main representations.

### Long-format table

```text
results/cre_turnover_all_species.tsv
```

Contains one row for every reference CRE–target species comparison together with the complete classification evidence.

For the primary dataset this corresponds to:

```text
337 reference CREs × 39 target species
```

comparisons.

### CRE-by-species matrix

```text
results/cre_turnover_matrix.tsv
```

Stores one reference CRE per row and one target species per column, with each cell containing one of the four CRE states.

This matrix is the principal input for phylogenetic visualization and downstream candidate analyses.

### Species summary

```text
results/species_summary.tsv
```

Summarizes for each species:

- mapped and unmapped reference CREs;
- counts of all four CRE states;
- positional-match rates;
- turnover-candidate rates;
- no-detected-CRE rates;
- FBgn support among positional matches.

Rates are reported both relative to all reference CREs and relative to evaluable CREs.

For evaluable rates:

```text
evaluable CREs = mapped reference CREs
```

so `uncertain` comparisons are excluded from the denominator.

---

## 4. Species-level quality control

A separate QC table is generated as:

```text
results/species_qc_summary.tsv
```

This integrates information from the SCRMshaw predictions, lifted reference CREs, and final classifications.

The QC evaluates:

- reference-CRE mapping rate;
- number of SCRMshaw predictions;
- sequence-ID compatibility between lifted intervals and predictions;
- occurrence of positional overlaps;
- reciprocal-overlap counts;
- CRE-state composition;
- same-FBgn prediction support.

Potential technical issues are recorded through QC flags rather than silently altering the classifications.

The primary QC thresholds are defined in `classification_config.sh` and are intended to identify potentially unreliable species-level results, not to redefine individual CRE states.

---

## Key outputs

### Species-specific classifications

```text
results/turnover_by_species/
    dmel_to_<species>_cre_turnover.tsv
```

Detailed classification evidence for every reference CRE in every target species.

### Combined long-format table

```text
results/cre_turnover_all_species.tsv
```

Complete reference CRE–species classification dataset.

### CRE-state matrix

```text
results/cre_turnover_matrix.tsv
```

Compact cross-species representation of the four CRE states.

### Species-level summary

```text
results/species_summary.tsv
```

State counts and rates for every target species.

### Species-level QC

```text
results/species_qc_summary.tsv
```

Technical QC metrics covering mapping, SCRMshaw prediction availability, coordinate compatibility, and positional overlap.

---

## Additional summaries and visualization

Several scripts in this stage operate on the completed classification outputs but are not part of the four-step core classification wrapper.

### CRE and species summaries

```text
scripts/plot_cre_turnover.py
```

generates species- and CRE-level summary figures and writes additional tables including:

```text
results/CRE_conservation_summary.tsv
results/species_plot_qc_summary.tsv
```

The CRE-level summary reports, for each reference CRE, its number of positional matches, turnover candidates, non-detected CREs, and uncertain comparisons across the target species.

---

### Phylogeny-wide CRE-state heatmap

```text
scripts/plot_tree_heatmap.py
```

combines:

```text
cre_turnover_matrix.tsv
phylogenetic tree
species manifest
species trait annotations
```

to visualize CRE states across the complete phylogeny together with climatic-zone information.

The resulting figure is written in PNG and PDF format under:

```text
results/figures/
```

*D. melanogaster* is displayed as the reference species but is not included in target-species CRE statistics.

---

### Focal-clade heatmaps

```text
scripts/plot_clade_heatmaps.py
```

visualizes selected focal-clade CRE-state contrasts together with the corresponding pruned phylogeny and climate annotations.

The CREs included in these figures are additionally recorded in:

```text
results/focal_clade_heatmap_CREs.tsv
```

These figures are downstream representations of the completed CRE-state matrix and do not alter the underlying classifications.

---

### Prediction-source summary

```text
scripts/summarize_prediction_source.py
```

compares SCRMshaw prediction depth between species whose predictions were generated within this project and species supplied from the external prediction dataset.

This provides an additional technical check for systematic differences between the two SCRMshaw prediction sources.

---

## Quality control

The workflow performs QC at several levels:

| Stage | Main check |
|---|---|
| Input | expected 337 reference CREs and required upstream files |
| Species classification | exactly one row per reference CRE |
| Combined results | all reference CREs represented for every target species |
| Alignment | mapped versus unmapped reference CREs |
| Sequence compatibility | overlap of sequence IDs between lifted CREs and predictions |
| SCRMshaw coverage | prediction counts per species |
| Positional comparison | availability of genomic overlaps and reciprocal matches |
| Final dataset | only the four defined CRE states are accepted |

The QC metrics are retained separately from the biological classification so that technical limitations can be evaluated without automatically interpreting them as CRE divergence.

---

## Adapting the workflow

The positional conservation threshold can be changed through:

```text
RECIPROCAL_OVERLAP
```

and the regulatory-neighborhood threshold through:

```text
LOCAL_GENE_DISTANCE
```

in `config/classification_config.sh`.

Changes to either value alter the operational CRE-state definitions and should therefore be treated as distinct analysis settings.

Changes to the reference CRE set require:

```text
EXPECTED_REFERENCE_CRES
```

to be updated accordingly.

Alternative QC thresholds can be specified through the corresponding:

```text
QC_*
```

configuration variables without changing the primary CRE-state definitions.

Systematic evaluation of alternative overlap and local-distance thresholds is handled separately in the downstream sensitivity analysis rather than by modifying the primary classification outputs.

---

## Reproducibility

The workflow retains the complete evidence underlying each CRE-state assignment, including homologous target coordinates, positional-overlap statistics, same-FBgn support, local candidate information, and alignment status. Species-specific classifications are subsequently combined without changing the individual assignments, allowing every entry in the final CRE-state matrix to be traced back to its upstream *D. melanogaster* reference CRE, pairwise coordinate projection, and target-species SCRMshaw prediction data.

Classification parameters and QC thresholds are centralized in `config/classification_config.sh`, while diagnostic summaries and visualization scripts operate on the completed classification tables without modifying them. This separates the primary state assignment from downstream interpretation and visualization.
