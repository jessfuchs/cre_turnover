# Tier-1 CRE Candidate Analysis and Gene Prioritization

This workflow identifies, validates, and prioritizes lineage-specific Tier-1 CRE-state contrasts and links the resulting *Drosophila melanogaster* CRE candidates to candidate target genes.

Tier-1 contrasts are restricted to:

```text
positional_match ↔ turnover_candidate
```

Candidate evidence is evaluated across predefined focal clades, validated against the detailed species-level CRE classifications, summarized across clades, and subsequently converted from CRE-level to gene-level evidence.

No arbitrary weighted candidate score is used.

---

## Workflow overview

```text
CRE-state matrix
      +
focal-clade definitions
      │
      ├──────────────┐
      ▼              ▼
 Focal Tier-1    Singleton Tier-1
   contrasts        contrasts
      │              │
      └──────┬───────┘
             ▼
      Candidate QC
             │
             ▼
 Cross-clade recurrence
   and prioritization
             │
             ▼
     Prioritized CREs
             │
             ▼
 Dmel gene-distance annotation
             │
             ▼
 Primary / secondary genes
             │
             ▼
      CRE × FBgn table
             │
             ▼
      Gene-level summary
```

The complete workflow is executed with:

```bash
bash run_candidate_pipeline.sh
```

---

## Directory layout

```text
01_candidate_analysis/
├── config/
│   ├── candidate_config.sh
│   └── focal_clades.tsv
├── scripts/
│   ├── 01_analyze_focal_clades.py
│   ├── 02_analyze_secondary_tier1.py
│   ├── 03_qc_tier1_candidates.py
│   ├── 04_prioritize_tier1_candidates.py
│   ├── 05_annotate_candidate_gene_distances.py
│   ├── 06_assign_candidate_genes.py
│   ├── 07_explode_candidate_fbgns.py
│   └── 08_build_gene_level_summary.py
├── results/
│   ├── focal_clades/
│   ├── secondary_clades/
│   ├── qc/
│   └── tables/
└── run_candidate_pipeline.sh
```

Paths and analysis parameters are defined in:

```text
config/candidate_config.sh
```

Focal and comparison species are defined independently in:

```text
config/focal_clades.tsv
```

---

## Core analysis defaults

| Parameter | Setting |
|---|---:|
| Reference species | *D. melanogaster* |
| Expected reference CREs | 337 |
| Tier-1 states | `positional_match` ↔ `turnover_candidate` |
| Focal recurrence minimum | 2 clades |
| Secondary recurrence minimum | 2 clades |
| Tier-1 positional QC overlap | ≥ 0.50 |
`no_detected_CRE` and `uncertain` are not included in the Tier-1 definition.

---

## Input data

### CRE-state matrix

The main input is:

```text
04_cre_classification/results/cre_turnover_matrix.tsv
```

containing one state per *D. melanogaster* reference CRE and target species.

### Detailed species classifications

Candidate QC uses:

```text
04_cre_classification/results/turnover_by_species/
```

to verify the species-level evidence underlying each Tier-1 state.

### D. melanogaster reference CREs

Reference coordinates and original FBgn target-gene annotations are obtained from:

```text
02_mapping_orthologs/reference_cres/dmel_reference_cres.tsv
```

### SCRMshaw gene annotations

Candidate-gene distances are recovered from:

```text
02_mapping_orthologs/results/SO_all_species_fbgn.tsv
```

### Focal-clade metadata

Focal groups are defined in:

```text
config/focal_clades.tsv
```

using:

```text
group_name
focal_species
comparison_species
```

Species and climatic-zone metadata are obtained from the combined species manifest and the trait annotations used in the classification workflow.

---

## Configuration

Important settings in:

```text
config/candidate_config.sh
```

include:

```text
CRE_TURNOVER_MATRIX
REFERENCE_CRES_TSV
TURNOVER_BY_SPECIES_DIR
SO_ALL_SPECIES_FBGN

FOCAL_CLADES_FILE
COMBINED_MANIFEST
SPECIES_TRAITS

REFERENCE_SPECIES
EXPECTED_REFERENCE_CRES

FOCAL_RECURRENCE_MIN_CLADES
SECONDARY_RECURRENCE_MIN_CLADES
TIER1_QC_MIN_RECIPROCAL_OVERLAP
```

Routine path or threshold changes should be made in the configuration files rather than directly in the worker scripts.

---

## Run the workflow

From:

```text
05_downstream_analyses/01_candidate_analysis/
```

execute:

```bash
bash run_candidate_analysis.sh
```

The wrapper identifies focal and singleton Tier-1 contrasts, validates the candidates, prioritizes recurrent evidence, annotates CRE-to-gene distances, assigns candidate genes, and builds the final gene-level summary.

---

## Tier-1 candidate identification

### Focal Tier-1

For each predefined clade, the comparison species must share one consensus CRE state.

A focal Tier-1 event requires:

```text
focal state != comparison consensus
AND
{focal state, comparison consensus}
=
{positional_match, turnover_candidate}
```

Thus both contrast directions are retained.

The focal analysis also records other CRE-state patterns for descriptive purposes, but only Tier-1 contrasts are carried forward into candidate prioritization.

Outputs are written below:

```text
results/focal_clades/
```

with the overall summary:

```text
results/focal_clade_summary.tsv
```

### Secondary Tier-1

A complementary singleton analysis searches the same clades without requiring the predefined focal species to be the discordant lineage.

A valid Secondary Tier-1 event requires exactly one species to differ from a common state shared by all remaining clade members, again specifically through:

```text
positional_match ↔ turnover_candidate
```

This retains potentially informative lineage-specific contrasts that occur outside the predefined focal species.

Secondary evidence is summarized across clades and considered recurrent when observed in at least two clades under the primary configuration.

---

## Candidate QC and prioritization

Tier-1 candidates are validated against the detailed per-species classification tables before prioritization.

QC checks include:

- agreement with the expected Tier-1 state pattern;
- successful homologous-locus mapping;
- reciprocal positional overlap for `positional_match`;
- local same-FBgn support for `turnover_candidate`;
- consistency between focal and singleton analyses.

Only:

```text
candidate_qc = PASS
```

evidence is used for prioritization.

The candidate categories are:

```markdown
| Candidate priority | Definition | Downstream priority |
|---|---|---|
| `focal_recurrent` | focal support in ≥2 clades | `high` |
| `focal_plus_secondary_recurrent` | focal plus secondary-only support | `high` |
| `secondary_recurrent` | secondary-only support in ≥2 clades | `medium` |
| `focal_single` | focal support in one clade | `medium` |
| `secondary_single` | secondary-only support in one clade | `exploratory` |
The main prioritized CRE table is:

```text
results/tables/tier1_candidates_prioritized.tsv
```

Only QC-passing clade evidence contributes to these priorities. :contentReference[oaicite:0]{index=0}

---

## Candidate-gene assignment

Prioritized CREs are matched to their corresponding *D. melanogaster* SCRMshaw predictions by genomic interval overlap.

When more than one prediction overlaps a CRE, the best match is selected deterministically using overlap strength and coordinate agreement. Ambiguous, partial, or non-exact matches are retained as QC information rather than discarded. :contentReference[oaicite:1]{index=1}

Candidate genes are then derived from the flanking and next-flanking FBgn annotations of the matched prediction.

Genes are ranked by absolute CRE-to-gene distance:

```text
closest unique FBgn       → primary_candidate
second-closest unique FBgn → secondary_candidate
```

Equal-distance cases are resolved deterministically. Agreement with the original SCRMshaw target-gene assignment is retained as separate QC information and does not determine the distance-based ranking. :contentReference[oaicite:2]{index=2}

---

## Gene-level summary

CRE-level assignments are expanded into unique:

```text
CRE × FBgn
```

records using the union of:

- original SCRMshaw target genes;
- primary distance-derived candidate genes;
- secondary distance-derived candidate genes.

Each FBgn is annotated according to whether it represents a primary candidate, secondary candidate, or reference target only. :contentReference[oaicite:3]{index=3}

These records are then aggregated to one row per candidate gene. The final table summarizes the number and type of supporting CREs, recurrent and focal support, candidate priority, clade support, and gene-assignment QC.

No weighted gene-level score is calculated; genes are ordered using transparent evidence variables. :contentReference[oaicite:4]{index=4}

---

## Key outputs

| Output | Purpose |
|---|---|
| `results/focal_clade_summary.tsv` | summary of focal-clade CRE-state contrasts |
| `results/secondary_clades/secondary_tier1_all_clades.tsv` | all secondary Tier-1 singleton evidence |
| `results/qc/tier1_candidate_qc_summary.tsv` | QC status of Tier-1 candidates |
| `results/tables/tier1_candidates_prioritized.tsv` | prioritized unique Tier-1 CRE candidates |
| `results/tables/tier1_candidate_gene_assignments.tsv` | primary and secondary candidate genes |
| `results/tables/tier1_candidate_fbgn_exploded.tsv` | one row per unique CRE × FBgn association |
| `results/tables/tier1_gene_level_summary.tsv` | final one-row-per-gene candidate summary |

---

## Quality control

QC is performed at several stages:

| Stage | Main check |
|---|---|
| CRE-state input | expected 337 unique reference CREs |
| Focal analysis | unanimous comparison-species state |
| Singleton analysis | exactly one discordant species |
| Candidate validation | mapping and state-specific evidence |
| Cross-step validation | focal Tier-1 reproduced by singleton analysis |
| Prioritization | only QC-passing evidence retained |
| Dmel matching | coordinate and overlap ambiguities recorded |
| Gene assignment | deterministic distance ranking and QC |
| Gene tables | unique CRE × FBgn and FBgn records |

Gene-distance and gene-assignment QC flags are descriptive and do not automatically exclude otherwise valid Tier-1 candidates.

---

## Adapting the workflow

Focal groups can be changed in:

```text
config/focal_clades.tsv
```

Recurrence thresholds are controlled through:

```text
FOCAL_RECURRENCE_MIN_CLADES
SECONDARY_RECURRENCE_MIN_CLADES
```

and the positional-match QC threshold through:

```text
TIER1_QC_MIN_RECIPROCAL_OVERLAP
```

The latter should remain consistent with the reciprocal-overlap criterion used for the primary CRE-state classification.

If the reference CRE set or Tier-1 definition changes, the candidate workflow and dependent sensitivity or climate analyses should be regenerated consistently.

---

## Reproducibility

Candidate selection is based on explicit categorical rules rather than manual selection or an opaque weighted score. Focal and singleton evidence is retained at clade level, and only QC-passing evidence contributes to candidate prioritization.

CRE-to-gene matching and candidate-gene assignment use deterministic overlap- and distance-based rules, while ambiguous matches and discrepancies with the original SCRMshaw target-gene annotation are retained as QC information. The final tables preserve CRE IDs, FBgn assignments, clade support, priority classes, and gene-assignment evidence so that each gene-level candidate can be traced back to the underlying CRE-state contrasts.
