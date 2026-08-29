#!/usr/bin/env python3

"""
Export a compact thesis-ready table of Azteca focal Tier-1 CRE candidates.

Input
-----
downstream_analyses/candidate_analysis/results/tables/
    tier1_candidates_prioritized_with_sensitivity.tsv

Selection
---------
Keep candidates with focal Tier-1 support in:
    azteca_affinis_miranda_group

Output
------
downstream_analyses/climate_analysis/results/
    azteca_focal_tier1_candidates.tsv

The output contains only the most relevant columns for interpretation:
    - CRE ID
    - focal direction
    - FBgn target genes
    - number of Tier-1 clades
    - Tier-1 clades
    - robustness
    - priority stability
"""

from pathlib import Path

import pandas as pd


# ============================================================
# Paths
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPT_DIR.parent
DOWNSTREAM_DIR = ANALYSIS_DIR.parent
PROJECT_DIR = DOWNSTREAM_DIR.parent

INPUT_FILE = (
    PROJECT_DIR
    / "downstream_analyses"
    / "candidate_analysis"
    / "results"
    / "tables"
    / "tier1_candidates_prioritized_with_sensitivity.tsv"
)

RESULTS_DIR = (
    ANALYSIS_DIR
    / "results"
)

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUT_FILE = (
    RESULTS_DIR
    / "azteca_focal_tier1_candidates.tsv"
)


# ============================================================
# Constants
# ============================================================

AZTECA_CLADE = "azteca_affinis_miranda_group"

EXPECTED_CRE_IDS = {
    "DMEL_CRE_00080",
    "DMEL_CRE_00309",
    "DMEL_CRE_00310",
    "DMEL_CRE_00330",
    "DMEL_CRE_00331",
}


# ============================================================
# Input check
# ============================================================

if not INPUT_FILE.exists():

    raise SystemExit(
        f"ERROR: input file not found:\n{INPUT_FILE}"
    )


# ============================================================
# Load
# ============================================================

df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    dtype=str,
).fillna("")


required_columns = {
    "dmel_cre_id",
    "focal_tier1_clades",
    "discordant_states",
    "consensus_states",
    "tier1_state_directions",
    "fbgn_target_genes",
    "n_fbgn_target_genes",
    "n_total_tier1_clades",
    "tier1_clades",
    "percent_candidate_retained",
    "robust_candidate_all_9",
    "robustness_class",
    "percent_same_priority_as_primary",
    "robust_priority_all_9",
    "priority_stability",
}


missing = (
    required_columns
    - set(df.columns)
)


if missing:

    raise SystemExit(
        "ERROR: input table missing required columns:\n"
        + "\n".join(
            sorted(
                missing
            )
        )
    )


# ============================================================
# Select Azteca focal Tier-1 candidates
# ============================================================

def contains_clade(
    value,
    clade,
):

    entries = [
        item.strip()
        for item in str(value).split("|")
        if item.strip()
    ]

    return clade in entries


azteca = df.loc[
    df["focal_tier1_clades"].apply(
        lambda x: contains_clade(
            x,
            AZTECA_CLADE,
        )
    )
].copy()


# ============================================================
# Sanity check expected candidates
# ============================================================

observed_ids = set(
    azteca[
        "dmel_cre_id"
    ]
)


missing_expected = (
    EXPECTED_CRE_IDS
    - observed_ids
)

unexpected = (
    observed_ids
    - EXPECTED_CRE_IDS
)


if missing_expected:

    raise SystemExit(
        "ERROR: expected Azteca focal candidates are missing:\n"
        + "\n".join(
            sorted(
                missing_expected
            )
        )
    )


if unexpected:

    print(
        "WARNING: additional Azteca focal candidates detected:"
    )

    for cre_id in sorted(
        unexpected
    ):
        print(
            f"  {cre_id}"
        )


# ============================================================
# Determine Azteca-specific focal direction
#
# We know the discordant species for Azteca is dazt.
# The source table can contain multiple clades/directions for
# recurrent candidates, so derive the Azteca direction directly
# from discordant/consensus states only when unambiguous.
# ============================================================

def infer_focal_direction(
    row,
):

    discordant = str(
        row[
            "discordant_states"
        ]
    )

    consensus = str(
        row[
            "consensus_states"
        ]
    )


    # Simple single-direction case.
    if (
        discordant == "turnover_candidate"
        and consensus == "present"
    ):
        return "focal_turnover"


    if (
        discordant == "present"
        and consensus == "turnover_candidate"
    ):
        return "focal_present"


    # Recurrent candidates can contain both directions because
    # different clades contribute opposite contrasts.
    # Use the known Azteca focal candidates explicitly.
    explicit = {
        "DMEL_CRE_00080": "focal_present",
        "DMEL_CRE_00309": "focal_present",
        "DMEL_CRE_00310": "focal_turnover",
        "DMEL_CRE_00330": "focal_turnover",
        "DMEL_CRE_00331": "focal_turnover",
    }


    return explicit.get(
        row[
            "dmel_cre_id"
        ],
        "ambiguous",
    )


azteca[
    "azteca_focal_direction"
] = azteca.apply(
    infer_focal_direction,
    axis=1,
)


# ============================================================
# Add simple biological interpretation label
# ============================================================

direction_label = {
    "focal_turnover":
        "Azteca turnover candidate; comparison species present",

    "focal_present":
        "Azteca present; comparison species turnover candidates",

    "ambiguous":
        "Ambiguous",
}


azteca[
    "azteca_pattern"
] = azteca[
    "azteca_focal_direction"
].map(
    direction_label
)


# ============================================================
# Select thesis-relevant columns
# ============================================================

output = azteca[
    [
        "dmel_cre_id",
        "azteca_focal_direction",
        "azteca_pattern",
        "fbgn_target_genes",
        "n_fbgn_target_genes",
        "n_total_tier1_clades",
        "tier1_clades",
        "percent_candidate_retained",
        "robust_candidate_all_9",
        "robustness_class",
        "percent_same_priority_as_primary",
        "robust_priority_all_9",
        "priority_stability",
    ]
].copy()


# ============================================================
# Convert numeric fields
# ============================================================

numeric_columns = [
    "n_fbgn_target_genes",
    "n_total_tier1_clades",
    "percent_candidate_retained",
    "percent_same_priority_as_primary",
]


for column in numeric_columns:

    output[
        column
    ] = pd.to_numeric(
        output[
            column
        ],
        errors="coerce",
    )


# ============================================================
# Sort candidates
#
# Priority:
#   1. robust
#   2. stable priority
#   3. recurrent across more clades
# ============================================================

robust_rank = {
    "robust": 0,
    "moderate": 1,
    "sensitive": 2,
}


stability_rank = {
    "stable": 0,
    "variable": 1,
}


output[
    "_robust_rank"
] = output[
    "robustness_class"
].map(
    robust_rank
).fillna(
    99
)


output[
    "_stability_rank"
] = output[
    "priority_stability"
].map(
    stability_rank
).fillna(
    99
)


output = (
    output
    .sort_values(
        [
            "_robust_rank",
            "_stability_rank",
            "n_total_tier1_clades",
            "dmel_cre_id",
        ],
        ascending=[
            True,
            True,
            False,
            True,
        ],
        kind="mergesort",
    )
    .drop(
        columns=[
            "_robust_rank",
            "_stability_rank",
        ]
    )
    .reset_index(
        drop=True
    )
)


# ============================================================
# Save
# ============================================================

output.to_csv(
    OUT_FILE,
    sep="\t",
    index=False,
)


# ============================================================
# Console output
# ============================================================

print()
print("=" * 88)
print("AZTECA FOCAL TIER-1 CANDIDATES")
print("=" * 88)
print()


display_columns = [
    "dmel_cre_id",
    "azteca_focal_direction",
    "fbgn_target_genes",
    "n_total_tier1_clades",
    "percent_candidate_retained",
    "robustness_class",
    "priority_stability",
]


print(
    output[
        display_columns
    ].to_string(
        index=False
    )
)


print()
print(
    f"Number of candidates: {len(output)}"
)

print()
print(
    f"Wrote:\n{OUT_FILE}"
)
