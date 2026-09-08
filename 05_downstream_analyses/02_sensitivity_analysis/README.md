# CRE-state and Tier-1 Candidate Sensitivity Analysis

This workflow evaluates whether the primary CRE-state classifications and prioritized Tier-1 candidates are robust to alternative positional-overlap and local target-gene distance thresholds.

The same CRE classifier used for the primary analysis is rerun across a predefined parameter grid. All scenarios are compared with the primary setting, and candidate robustness is evaluated by re-identifying Tier-1 support across the same focal clades used in `01_candidate_analysis/`.

---

## Workflow overview

```text
Primary classification inputs
          │
          ▼
  Sensitivity parameter grid
          │
          ▼
CRE classification per scenario
          │
          ▼
Primary-scenario regression
          │
          ├───────────────┐
          ▼               ▼
 Global CRE-state     Tier-1 candidate
    stability           stability
          │               │
          └───────┬───────┘
                  ▼
       Parameter-specific
       sensitivity summaries
                  │
                  ▼
        Tables and figures
```

---

## Directory layout

```text
02_sensitivity_analysis/
├── config/
│   └── sensitivity_config.sh
├── scripts/
│   ├── 01_generate_sensitivity_scenarios.py
│   ├── 02_check_primary_regression.py
│   ├── 03_summarize_global_sensitivity.py
│   ├── 04_summarize_candidate_sensitivity.py
│   ├── 05_add_sensitivity_to_prioritized_candidates.py
│   ├── 06_summarize_parameter_sensitivity.py
│   ├── plot_ranked_tier1.py
│   ├── plot_recurrence_robustness.py
│   └── plot_sensitivity_matrix.py
├── results/
│   ├── scenarios/
│   └── figures/
└── run_sensitivity_analysis.sh
```

Paths, parameter values, recurrence thresholds, and output locations are defined in:

```text
config/sensitivity_config.sh
```

---

## Core analysis defaults

The primary sensitivity grid is:

| Parameter | Values |
|---|---|
| Reciprocal positional overlap | 0.25, 0.50, 0.75 |
| Local same-FBgn distance | 12 kb, 24 kb, 48 kb |
| Number of scenarios | 9 |
| Primary scenario | `ov050_dist24000` |
| Expected reference CREs | 337 |
| Focal recurrence minimum | 2 clades |
| Secondary recurrence minimum | 2 clades |
| Moderate candidate robustness | ≥ 66.6% retention |

The primary scenario reproduces the classification thresholds used in the main CRE-state analysis:

```text
reciprocal overlap = 0.50
local gene distance = 24 kb
```

---

## Input data

Sensitivity classification reuses the same upstream inputs as the primary CRE-state analysis:

```text
02_mapping_orthologs/reference_cres/dmel_reference_cres.tsv
02_mapping_orthologs/results/SO_all_species_fbgn.tsv
03_pairwise_wga/lifted_cres_dmel/
03_pairwise_wga/target_species.txt
04_cre_classification/scripts/01_classify_target_cres.py
```

The original species-level classifications are used as the regression baseline:

```text
04_cre_classification/results/turnover_by_species/
```

Candidate-level sensitivity additionally uses:

```text
01_candidate_analysis/results/tables/tier1_candidates_prioritized.tsv
01_candidate_analysis/config/focal_clades.tsv
```

---

## Configuration

Important variables in:

```text
config/sensitivity_config.sh
```

include:

```text
SENSITIVITY_OVERLAPS
SENSITIVITY_DISTANCES
PRIMARY_SENSITIVITY_SCENARIO

EXPECTED_REFERENCE_CRES

FOCAL_RECURRENCE_MIN_CLADES
SECONDARY_RECURRENCE_MIN_CLADES
CANDIDATE_MODERATE_ROBUSTNESS_MIN
```

Input and output paths are also centralized in this configuration file.

---

## Run the workflow

From:

```text
05_downstream_analyses/02_sensitivity_analysis/
```

execute:

```bash
bash run_sensitivity_analysis.sh
```

The workflow generates the sensitivity classifications, verifies the primary scenario, summarizes global CRE-state stability, and evaluates Tier-1 candidate robustness.

---

## Sensitivity analysis

### Classification scenarios

The full factorial combination of overlap and distance thresholds produces nine sensitivity scenarios.

Scenario names encode both parameters, for example:

```text
ov025_dist12000
ov050_dist24000
ov075_dist48000
```

For every scenario and target species, the standard CRE classifier is rerun with only the reciprocal-overlap and local-distance thresholds changed.

Scenario-specific classifications are stored under:

```text
results/scenarios/<scenario>/
```

A manifest records the parameter values associated with each scenario:

```text
results/sensitivity_scenarios.tsv
```

### Primary regression

Before sensitivity results are interpreted, the primary scenario is compared directly with the original species-level classifications.

```text
ov050_dist24000
```

must reproduce every baseline CRE-state assignment exactly. Any difference causes the regression check to fail.

### Global CRE-state stability

Each non-primary scenario is compared with the primary scenario at three levels:

```text
all CRE × species comparisons
individual species
individual reference CREs
```

State transitions are also recorded to show which operational CRE states change when classification thresholds are altered.

### Tier-1 candidate stability

The prioritized CREs from `01_candidate_analysis/` are re-evaluated across all scenarios using the same focal and singleton Tier-1 definitions and recurrence thresholds.

For each CRE, the analysis records:

```text
candidate retention
candidate-priority stability
scenarios in which the candidate is lost
scenarios in which its priority changes
```

Sensitivity is added as supporting evidence to the original prioritized candidate table; it does not replace the primary candidate ranking.

Candidate robustness is classified as:

| Class | Definition |
|---|---|
| `robust` | retained in all scenarios |
| `moderate` | retained in at least 66.6% of scenarios |
| `sensitive` | retained below 66.6% |

### Parameter-specific sensitivity

The effects of the two classification parameters are also considered separately.

Reciprocal-overlap sensitivity is evaluated while holding the local distance at 24 kb, whereas local-distance sensitivity is evaluated while holding reciprocal overlap at 0.50.

This distinguishes candidates that are primarily sensitive to positional-overlap stringency from those affected by the local target-gene distance threshold.

---

## Key outputs

| Output | Purpose |
|---|---|
| `results/sensitivity_scenarios.tsv` | scenario definitions |
| `results/primary_regression_check.tsv` | baseline reproduction check |
| `results/sensitivity_global_summary.tsv` | global state stability across scenarios |
| `results/sensitivity_species_stability.tsv` | species-level stability |
| `results/sensitivity_cre_stability.tsv` | reference-CRE-level stability |
| `results/sensitivity_state_transitions.tsv` | state changes relative to the primary scenario |
| `results/candidate_sensitivity_summary.tsv` | Tier-1 candidate robustness |
| `results/candidate_sensitivity_retention_matrix.tsv` | candidate retention across scenarios |
| `results/candidate_sensitivity_priority_matrix.tsv` | candidate priority across scenarios |

Sensitivity annotations are additionally merged into:

```text
../01_candidate_analysis/results/tables/
    tier1_candidates_prioritized_with_sensitivity.tsv
```

Parameter-specific candidate summaries and sensitivity figures are written under:

```text
results/
results/figures/
```

---

## Quality control

The workflow verifies that every scenario contains the expected 337 reference CREs for each target species, accepts only the four defined CRE states, and requires unique CRE identifiers.

Most importantly, the primary sensitivity scenario must reproduce the original classification exactly before downstream sensitivity summaries are generated.

---

## Adapting the workflow

The sensitivity grid and primary scenario are controlled through:

```text
SENSITIVITY_OVERLAPS
SENSITIVITY_DISTANCES
PRIMARY_SENSITIVITY_SCENARIO
```

Candidate recurrence thresholds should remain synchronized with those used in `01_candidate_analysis/`.

Changing the primary classification thresholds requires the primary scenario and corresponding baseline classification to be updated consistently.

The moderate-robustness threshold can be changed through:

```text
CANDIDATE_MODERATE_ROBUSTNESS_MIN
```

---

## Reproducibility

All sensitivity scenarios are generated with the same classifier and upstream inputs as the primary analysis, with only the tested classification thresholds changed. The primary regression check ensures that the designated reference scenario exactly reproduces the baseline classification.

Scenario definitions, global state changes, candidate retention, priority changes, and run metadata are retained so that sensitivity conclusions can be traced back to the corresponding parameter settings and species-level CRE classifications.
