# Pairwise Whole-genome Alignment and CRE Liftover

This workflow identifies homologous genomic positions of the fixed *Drosophila melanogaster* reference CREs across the target species.

Rather than using a single multi-species alignment, each target genome is aligned independently to *D. melanogaster*. The resulting pairwise alignments are processed through the UCSC chain/net workflow to obtain synteny-filtered chain files, which are subsequently used to project the *D. melanogaster* reference CRE coordinates into each target genome using `liftOver`.

The principal outputs are the species-specific synteny-filtered chain files and the mapped and unmapped reference CRE intervals used by the downstream CRE-state classification.

---

## Workflow overview

The workflow consists of three main stages:

1. prepare all genome assemblies for pairwise alignment;
2. align each target genome independently to *D. melanogaster* using LASTZ and the UCSC chain/net workflow;
3. project the fixed *D. melanogaster* reference CREs into each target genome using `liftOver`.

Conceptually:

```text
Genome FASTA files
        │
        ▼
FASTA → 2bit
+ chromosome sizes
        │
        ▼
D. melanogaster
        │
        ├── LASTZ → target species 1
        ├── LASTZ → target species 2
        ├── ...
        └── LASTZ → target species N
                │
                ▼
        UCSC chain/net workflow
                │
                ▼
        synteny-filtered chain
                │
                ▼
       Dmel reference CRE BED
                │
                ▼
             liftOver
          ┌─────┴─────┐
          ▼           ▼
       mapped      unmapped
        CREs          CREs
```

The complete workflow is submitted with:

```bash
bash run_wga_pipeline.sh
```

---

## Directory layout

```text
03_pairwise_wga/
├── bin/
├── config/
│   └── wga_config.sh
├── scripts/
│   ├── 01_prepare_wga_genomes.sh
│   ├── 02_run_pairwise_lastz.sh
│   └── 03_run_liftover.sh
├── slurm/
├── twobit/
├── chrom_sizes/
├── alignments_dmel/
├── lifted_cres_dmel/
├── logs/
├── target_species.txt
└── run_wga_pipeline.sh
```

Input locations, alignment parameters, `liftOver` settings, and SLURM resources are defined centrally in:

```text
config/wga_config.sh
```

---

## Core analysis defaults

The final analysis uses:

| Parameter | Setting |
|---|---|
| Reference species | `d_melanogaster` |
| Alignment strategy | independent pairwise alignment |
| LASTZ output format | `axt` |
| LASTZ ambiguous bases | `iupac` |
| LASTZ transitions | disabled |
| LASTZ step | 20 |
| LASTZ seed | `12of19` |
| LASTZ HSP threshold | 2200 |
| LASTZ gapped threshold | 4000 |
| LASTZ Y-drop | 3400 |
| LASTZ inner threshold | 2000 |
| `axtChain -linearGap` | `medium` |
| `axtChain -minScore` | 3000 |
| `liftOver -minMatch` | 0.50 |

The primary workflow uses *D. melanogaster* as the fixed reference genome. All other species present in the upstream species manifest are treated as independent alignment targets.

---

## Input data

### Species manifest and genome assemblies

The species framework is obtained from:

```text
01_scrmshaw/external/combined_manifest.tsv
```

The manifest must contain:

```text
slug
source
```

for each species.

The `source` field determines whether the corresponding genome FASTA is obtained from the generated or external SCRMshaw input collection:

```text
01_scrmshaw/generation/data/selected/<species>/genome.fa
01_scrmshaw/external/external_data/selected/<species>/genome.fa
```

The reference species is included in the manifest but excluded from the target-species list used for pairwise alignment.

---

### D. melanogaster reference CREs

The fixed reference CRE coordinates are taken from:

```text
02_mapping_orthologs/reference_cres/dmel_reference_cres.bed
```

The primary analysis uses the 337 *D. melanogaster* reference CREs constructed in the preceding workflow stage.

The BED file contains stable reference CRE identifiers and provides the coordinates projected into every target genome.

---

## Software requirements

The workflow requires:

- LASTZ;
- UCSC `faToTwoBit`;
- UCSC `twoBitInfo`;
- UCSC `axtChain`;
- UCSC `chainSort`;
- UCSC `chainPreNet`;
- UCSC `chainNet`;
- UCSC `netSyntenic`;
- UCSC `netChainSubset`;
- UCSC `liftOver`;
- GNU Bash;
- SLURM.

The UCSC and LASTZ executables used by the project are made available through the configured:

```text
bin/
```

directory.

---

## Configuration

All project-specific paths, alignment parameters, `liftOver` settings, and compute resources are defined in:

```text
config/wga_config.sh
```

Important configurable values include:

```text
COMBINED_MANIFEST

GENERATED_GENOME_ROOT
EXTERNAL_GENOME_ROOT

REFERENCE_CRES_BED
REFERENCE_SPECIES

LASTZ_*
AXTCHAIN_*
LIFTOVER_MIN_MATCH

SLURM_NODE
WGA_MEM
WGA_CPUS
MAX_WGA_JOBS
```

Routine changes should be made through the configuration file rather than directly in the worker scripts.

---

## Run the complete workflow

From:

```text
03_pairwise_wga/
```

execute:

```bash
bash run_wga_pipeline.sh
```

The wrapper first prepares all genome inputs and then submits one pairwise alignment job per target species as a SLURM array.

The `liftOver` stage is submitted with an `afterok` dependency and therefore starts only if all pairwise alignment jobs complete successfully.

For the primary 40-species framework, this corresponds to one *D. melanogaster* reference genome and 39 independent target-species alignments.

---

## 1. Prepare genomes for pairwise alignment

Genome FASTA files are converted to UCSC 2bit format:

```text
twobit/<species>.2bit
```

and chromosome-size files are generated as:

```text
chrom_sizes/<species>.sizes
```

using:

```text
faToTwoBit
twoBitInfo
```

The workflow automatically constructs:

```text
target_species.txt
```

from the combined species manifest, excluding the configured reference species.

Before alignment, the preparation stage verifies:

- availability of the species manifest;
- availability of the reference CRE BED;
- availability of all required genome FASTA files;
- successful generation of the 2bit and chromosome-size files;
- absence of duplicate target species;
- availability of the reference genome in the prepared files;
- compatibility between reference-CRE sequence identifiers and the *D. melanogaster* reference genome.

These checks ensure that the reference CRE coordinates and the genome assembly use compatible sequence identifiers before pairwise alignment begins.

---

## 2. Run pairwise whole-genome alignments

Each target genome is aligned independently against the *D. melanogaster* reference genome using LASTZ.

The initial pairwise alignment is written in AXT format and subsequently processed through the UCSC chain/net workflow:

```text
LASTZ
  ↓
axtChain
  ↓
chainSort
  ↓
chainPreNet
  ↓
chainNet
  ↓
netSyntenic
  ↓
netChainSubset
```

`axtChain` converts the LASTZ alignment into chain format, after which the chains are sorted and prepared for genome net construction.

`chainNet` generates target- and query-side nets, and `netSyntenic` annotates the target-side net with syntenic relationships.

Finally, `netChainSubset` extracts the chains represented in the syntenic net.

The resulting species-specific chain used for coordinate projection is:

```text
alignments_dmel/<species>/
    d_melanogaster.<species>.liftover.chain.gz
```

Intermediate AXT, chain, and net files are retained in the same species-specific alignment directory for reproducibility and QC.

---

## SLURM execution

Pairwise genome alignments are submitted as a SLURM array with one task per target species.

The primary configuration uses:

```text
Node             : abacus-2
Memory per job   : 5 GB
CPUs per job     : 1
Maximum parallel : 10 alignment jobs
```

These values are controlled through:

```text
config/wga_config.sh
```

and can be adapted to another compute environment.

---

## 3. Project reference CREs with liftOver

After all pairwise alignments have completed successfully, the fixed *D. melanogaster* reference CRE coordinates are projected into each target genome using the corresponding synteny-filtered chain file.

The primary analysis uses:

```text
-minMatch=0.50
```

requiring at least 50% of the reference interval to be mapped by `liftOver`.

For each target species, two files are generated:

```text
lifted_cres_dmel/<species>/
    dmel_reference_cres.<species>.bed

lifted_cres_dmel/<species>/
    dmel_reference_cres.<species>.unmapped.bed
```

The first contains successfully projected reference CRE intervals, while the second retains CREs that could not be projected under the selected mapping criterion.

Unmapped reference CREs are retained explicitly rather than being discarded, allowing downstream analyses to distinguish lack of positional mapping from an observed CRE-state difference.

---

## Key outputs

### Prepared genomes

```text
twobit/<species>.2bit
chrom_sizes/<species>.sizes
```

Standardized genome representations used by LASTZ and the UCSC chain/net workflow.

### Target-species list

```text
target_species.txt
```

Contains all species aligned independently against the *D. melanogaster* reference.

### Pairwise alignment chains

```text
alignments_dmel/<species>/
    d_melanogaster.<species>.liftover.chain.gz
```

Synteny-filtered chain files used for downstream coordinate projection.

### Mapped reference CREs

```text
lifted_cres_dmel/<species>/
    dmel_reference_cres.<species>.bed
```

Contains the homologous target-genome positions of successfully projected *D. melanogaster* reference CREs.

### Unmapped reference CREs

```text
lifted_cres_dmel/<species>/
    dmel_reference_cres.<species>.unmapped.bed
```

Retains reference CREs for which no acceptable coordinate projection was obtained.

---

## Quality control

QC is performed at both genome-preparation and coordinate-projection stages.

| Stage | Main check |
|---|---|
| Genome preparation | required genomes and tools are available |
| 2bit conversion | non-empty 2bit and chromosome-size files |
| Reference coordinates | CRE sequence IDs occur in the reference genome |
| Target list | no duplicate target species |
| Pairwise chains | all target species have a non-empty, valid gzip chain |
| liftOver | mapped + unmapped CREs equal the reference CRE count |
| liftOver | mapped reference CRE identifiers remain unique |

If any pairwise chain is missing or corrupt, `liftOver` is not performed.

The `liftOver` QC additionally requires, for every target species:

```text
mapped CREs + unmapped CREs = total reference CREs
```

and verifies that successfully mapped CRE identifiers are unique.

Any violation causes the liftOver stage to fail.

---

## Adapting the workflow

To change the reference species, update:

```text
REFERENCE_SPECIES
```

in `config/wga_config.sh`.

To modify alignment stringency, adjust the configured:

```text
LASTZ_*
AXTCHAIN_*
```

parameters.

The minimum fraction of a reference CRE required for successful coordinate projection can be changed through:

```text
LIFTOVER_MIN_MATCH
```

When changing the species framework, regenerate the genome-preparation outputs and `target_species.txt` from the updated combined manifest before rerunning the alignments.

Changes to the upstream *D. melanogaster* reference CRE definition require the `liftOver` outputs to be regenerated.

---

## Reproducibility

Each target species is aligned independently against the same *D. melanogaster* reference genome using a consistent LASTZ and UCSC chain/net configuration.

Intermediate alignment files and the final synteny-filtered chain files are retained, allowing individual pairwise alignments and CRE projections to be inspected or recomputed independently.

Mapped and unmapped reference CREs are stored separately for every target species so that downstream CRE-state classification can distinguish successful coordinate projection from non-evaluable reference positions.
