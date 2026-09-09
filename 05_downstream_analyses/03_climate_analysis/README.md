# Climate Association and Focal Tier-1 Enrichment

This workflow evaluates whether species-level CRE turnover is associated with climatic zone and whether predefined focal lineages are enriched for Tier-1 CRE-state contrasts.

The analysis combines descriptive species-level comparisons with phylogenetically controlled regression and a complementary focal-lineage enrichment analysis. Climate associations are evaluated using the proportion of `turnover_candidate` states among evaluable CRE loci, while focal enrichment is based on strict `positional_match` ↔ `turnover_candidate` singleton contrasts.

---

## Workflow overview

```text
Species-level CRE summary
        +
climatic-zone metadata
        │
        ▼
Turnover rate by climate
        │
        ├───────────────┐
        ▼               ▼
 Descriptive /      Phylogenetic
 Kruskal-Wallis        PGLS
                        │
                        ▼
              Climate-model estimates


CRE-state matrix
        +
focal-clade definitions
        │
        ▼
Tier-1 singleton contrasts
        │
        ▼
Observed focal-lineage share
        │
        ▼
Exact binomial enrichment
        +
BH correction
```

---

## Directory layout

```text
03_climate_analysis/
├── config/
│   └── climate_config.sh
├── 01_plot_turnover_by_climate.py
├── 02_test_turnover_by_climate.R
├── 03_plot_climate_pgls.R
├── 04_test_focal_tier1_enrichment.py
├── 05_plot_focal_tier1_enrichment.py
├── 06_plot_focal_summary.py
├── results/
│   └── figures/
└── run_climate_pipeline.sh
```

Input and output paths are centralized in:

```text
config/climate_config.sh
```

---

## Core analysis defaults

| Parameter | Setting |
|---|---|
| Climate metric | `turnover_candidate / mapped CREs` |
| Climatic zones | `TROP`, `ARID`, `TEMP`, `BORE` |
| PGLS response | empirical-logit turnover rate |
| PGLS predictor | climatic zone |
| PGLS phylogenetic structure | Pagel's λ |
| PGLS reference climate | `TROP` |
| PGLS model fitting | maximum likelihood |
| Tier-1 states | `positional_match` ↔ `turnover_candidate` |
| Focal-enrichment null | focal probability = `1 / clade size` |
| Focal-enrichment test | one-sided exact binomial |
| Multiple-testing correction | Benjamini–Hochberg |

The primary species-level metric is:

```text
turnover_rate_evaluable =
    turnover_candidate / mapped reference CREs
```

Only mapped reference CREs contribute to the denominator.

---

## Input data

### Species-level CRE classification

Species-level turnover counts and rates are obtained from:

```text
04_cre_classification/results/species_summary.tsv
```

Species names and prediction-source metadata are linked through:

```text
01_scrmshaw/external/combined_manifest.tsv
```

Climatic-zone annotations are obtained from:

```text
04_cre_classification/phylogeny/data/species_traits.tsv
```

### Phylogeny

Phylogenetically controlled analysis uses:

```text
04_cre_classification/phylogeny/results/
    301Fly_HOG_UCLDtree_40species.nw
```

The tree is pruned to the species represented in the climate-analysis table before model fitting.

### Focal Tier-1 analysis

The focal enrichment analysis uses:

```text
04_cre_classification/results/cre_turnover_matrix.tsv
```

together with the same focal-clade definitions used in the candidate analysis:

```text
01_candidate_analysis/config/focal_clades.tsv
```

---

## Configuration

Important paths in:

```text
config/climate_config.sh
```

include:

```text
SPECIES_SUMMARY
COMBINED_MANIFEST
SPECIES_TRAITS
PHYLOGENY_TREE

CRE_TURNOVER_MATRIX
FOCAL_CLADES_FILE

RESULTS_DIR
FIGURES_DIR
```

Analysis outputs are also defined centrally in the configuration file.

---

## Run the workflow

From:

```text
05_downstream_analyses/03_climate_analysis/
```

execute:

```bash
bash run_climate_pipeline.sh
```

The wrapper runs the species-level climate analysis, phylogenetically controlled PGLS, focal Tier-1 enrichment test, and associated figures.

---

## Species-level climate analysis

Species classifications are combined with climatic-zone metadata and evaluated using:

```text
turnover_rate_evaluable
```

The workflow first summarizes the distribution of turnover rates within each climatic zone and performs an exploratory Kruskal–Wallis comparison.

Because species are not statistically independent observations, the Kruskal–Wallis result is treated as descriptive rather than as the primary phylogenetically controlled inference. :contentReference[oaicite:3]{index=3}

### Phylogenetically controlled model

The primary climate test uses phylogenetic generalized least squares (PGLS).

The response is the empirical-logit transformed turnover proportion:

```text
logit((turnover_candidate + 0.5) / (mapped + 1))
```

and the fitted models are:

```text
null:
    logit turnover rate ~ 1

climate:
    logit turnover rate ~ climatic_zone
```

Phylogenetic covariance is modeled using Pagel's λ, which is estimated from the data. The null and climate models are fitted by maximum likelihood and compared to test the overall climatic-zone effect. :contentReference[oaicite:4]{index=4}

Model-based estimates and 95% confidence intervals are subsequently back-transformed to the turnover-rate scale for visualization. :contentReference[oaicite:5]{index=5}

---

## Focal Tier-1 enrichment

The second analysis asks whether the predefined focal species within each focal clade occurs as the discordant lineage more often than expected among Tier-1 singleton contrasts.

A Tier-1 singleton requires:

```text
all clade members ∈
    {positional_match, turnover_candidate}

AND

exactly one species differs
from all remaining species
```

Both directions are retained:

```text
focal turnover_candidate
vs.
comparison positional_match

and

focal positional_match
vs.
comparison turnover_candidate
```

Under the null model, every species in a clade is assumed to have equal probability of carrying the singleton difference:

```text
expected focal share = 1 / number of clade species
```

Observed focal counts are tested using a one-sided exact binomial test. P-values across focal clades are adjusted using the Benjamini–Hochberg procedure. :contentReference[oaicite:6]{index=6}

This enrichment analysis is exploratory because individual CRE contrasts are not guaranteed to represent statistically independent evolutionary events.

---

## Key outputs

| Output | Purpose |
|---|---|
| `results/climate_turnover_summary.tsv` | species-level turnover rates with climatic-zone annotation |
| `results/climate_turnover_descriptive_statistics.tsv` | descriptive turnover statistics by climate |
| `results/climate_turnover_kruskal.tsv` | exploratory Kruskal–Wallis result |
| `results/pgls_climate_model_comparison.tsv` | null-versus-climate PGLS comparison |
| `results/pgls_climate_coefficients.tsv` | fitted climate-model coefficients |
| `results/pgls_climate_predictions.tsv` | climate-specific PGLS estimates and 95% CI |
| `results/pgls_climate_stats.tsv` | global climate p-value and estimated Pagel's λ |
| `results/focal_tier1_candidates.tsv` | clade-wide Tier-1 singleton events |
| `results/focal_tier1_enrichment.tsv` | focal enrichment statistics and BH-adjusted p-values |

A human-readable PGLS report is written to:

```text
results/pgls_climate_summary.txt
```

Figures are written to:

```text
results/figures/
```

including:

```text
turnover_rate_by_climate.png
climate_PGLS.png
focal_tier1_enrichment.png
focal_tier1_summary.png
```

---

## Quality control

The workflow verifies species identifiers, climatic-zone annotations, turnover counts and rates, and compatibility between the climate table and phylogenetic tree.

For the species-level analysis, `turnover_rate_evaluable` is explicitly checked against:

```text
turnover_candidate / mapped
```

before downstream analyses are performed. :contentReference[oaicite:7]{index=7}

The focal analysis requires unique CRE identifiers, complete focal-clade species representation in the CRE-state matrix, and strict Tier-1 singleton patterns before events contribute to the enrichment test. :contentReference[oaicite:8]{index=8}

---

## Adapting the workflow

Input and output paths can be changed in:

```text
config/climate_config.sh
```

Focal species and comparison groups are inherited from:

```text
../01_candidate_analysis/config/focal_clades.tsv
```

and should therefore be changed there rather than independently within the climate analysis.

The climatic-zone vocabulary and PGLS reference level are currently defined directly in the analysis scripts as:

```text
TROP
ARID
TEMP
BORE
```

with `TROP` as the reference category. Changing these requires corresponding updates to the climate-analysis scripts and trait metadata.

---

## Reproducibility

Species-level climate comparisons, PGLS inference, and focal Tier-1 enrichment are based on deterministic upstream classifications and predefined species metadata. The phylogeny, climatic-zone assignments, focal-clade definitions, statistical outputs, and underlying singleton events are retained so that each climate or enrichment result can be traced back to the corresponding species-level CRE classifications.
