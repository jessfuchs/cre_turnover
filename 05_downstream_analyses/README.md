# Downstream Analyses

This stage performs the biological interpretation and robustness analyses of the CRE-state classifications generated in `04_cre_classification/`.

The downstream workflow consists of three complementary analyses:

- **Candidate analysis** identifies and prioritizes recurrent Tier-1 CRE contrasts and links candidate CREs to *Drosophila melanogaster* target genes.
- **Sensitivity analysis** evaluates the robustness of CRE-state assignments and Tier-1 candidates to alternative classification thresholds.
- **Climate analysis** tests whether CRE turnover is associated with climatic zone and evaluates focal-lineage enrichment of Tier-1 contrasts.

Detailed methods, parameters, QC procedures, and output descriptions are provided in the README of each subworkflow.

---

## Workflow overview

```text
04_cre_classification/
        │
        ▼
CRE-state matrix and
species-level classifications
        │
        ├──────────────────────────────┐
        │                              │
        ▼                              ▼
01_candidate_analysis/        03_climate_analysis/
        │                              ▲
        ▼                              │
Prioritized Tier-1 CREs                │
        │                              │
        ▼                              │
02_sensitivity_analysis/ ──────────────┘
        │
        ▼
Sensitivity-annotated
candidate evidence
```

The three analyses address different aspects of the same CRE-classification results and should be interpreted together rather than as independent pipelines.

---

## Directory layout

```text
05_downstream_analyses/
├── 01_candidate_analysis/
│   ├── config/
│   ├── results/
│   └── README.md
├── 02_sensitivity_analysis/
│   ├── config/
│   ├── results/
│   └── README.md
├── 03_climate_analysis/
│   ├── config/
│   ├── results/
│   └── README.md
└── README.md
```

---

## Analysis stages

### 1. Candidate analysis

```text
01_candidate_analysis/
```

Identifies focal and singleton Tier-1 CRE contrasts, validates candidate evidence, summarizes recurrence across focal clades, and produces prioritized CRE- and gene-level candidate tables.

Run with:

```bash
cd 01_candidate_analysis
bash run_candidate_analysis.sh
```

See:

```text
01_candidate_analysis/README.md
```

---

### 2. Sensitivity analysis

```text
02_sensitivity_analysis/
```

Re-runs CRE classification across alternative reciprocal-overlap and local target-gene distance thresholds and evaluates both global CRE-state stability and Tier-1 candidate robustness.

The primary parameter setting is:

```text
reciprocal overlap = 0.50
local gene distance = 24 kb
```

with sensitivity analyses around this baseline.

Run with:

```bash
cd 02_sensitivity_analysis
bash run_sensitivity_pipeline.sh
```

See:

```text
02_sensitivity_analysis/README.md
```

---

### 3. Climate analysis

```text
03_climate_analysis/
```

Evaluates species-level CRE turnover across climatic zones using descriptive and phylogenetically controlled analyses and tests whether predefined focal lineages are enriched for Tier-1 singleton contrasts.

Run with:

```bash
cd 03_climate_analysis
bash run_climate_pipeline.sh
```

See:

```text
03_climate_analysis/README.md
```

---

## Main inputs

The downstream analyses primarily depend on outputs from:

```text
04_cre_classification/results/
```

including:

```text
cre_turnover_matrix.tsv
species_summary.tsv
turnover_by_species/
```

Additional shared inputs include:

```text
02_mapping_orthologs/reference_cres/
04_cre_classification/phylogeny/
01_scrmshaw/external/combined_manifest.tsv
```

Focal-clade definitions are maintained in:

```text
01_candidate_analysis/config/focal_clades.tsv
```

and are reused where required by the sensitivity and climate analyses.

---

## Recommended execution order

Run the downstream analyses in the following order:

```bash
cd 01_candidate_analysis
bash run_candidate_analysis.sh

cd ../02_sensitivity_analysis
bash run_sensitivity_pipeline.sh

cd ../03_climate_analysis
bash run_climate_pipeline.sh
```

The candidate analysis should be completed before the sensitivity analysis because candidate robustness is evaluated from the prioritized Tier-1 candidate set.

The climate analysis uses the CRE-classification outputs and shared focal-clade definitions and can be rerun independently once these upstream inputs are available.

---

## Key outputs

The main downstream result directories are:

```text
01_candidate_analysis/results/
02_sensitivity_analysis/results/
03_climate_analysis/results/
```

Together they contain:

- prioritized Tier-1 CRE and candidate-gene evidence;
- robustness summaries across alternative classification thresholds;
- species-level and phylogenetically controlled climate analyses;
- focal-lineage Tier-1 enrichment results;
- downstream figures and QC outputs.

Refer to the corresponding subworkflow README for exact files and interpretation.

---

## Configuration

Each downstream workflow has its own configuration file:

```text
01_candidate_analysis/config/candidate_config.sh
02_sensitivity_analysis/config/sensitivity_config.sh
03_climate_analysis/config/climate_config.sh
```

Shared biological definitions, such as focal clades or primary CRE-classification thresholds, should be kept consistent across the affected analyses.

---

## Adapting the workflow

Changes to focal clades, candidate recurrence thresholds, sensitivity parameters, climatic-zone annotations, or upstream CRE classifications should be made through the relevant configuration or metadata files.

If an upstream CRE-state definition or reference CRE set changes, all affected downstream analyses should be regenerated.

---

## Reproducibility

The downstream analyses operate on fixed upstream CRE-state classifications and retain intermediate evidence, QC summaries, parameter settings, candidate annotations, statistical results, and figures within their respective result directories.

Each subworkflow is independently reproducible from its documented inputs and configuration, while shared identifiers and focal-clade definitions ensure consistent interpretation across candidate, sensitivity, and climate analyses.
