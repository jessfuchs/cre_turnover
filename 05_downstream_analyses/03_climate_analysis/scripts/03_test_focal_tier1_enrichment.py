#!/usr/bin/env python3

"""
Test whether focal Tier-1 contrasts are enriched within six focal clades.

Tier-1 singleton definition
---------------------------
All species in a focal clade must be either:

    present
    turnover_candidate

and exactly one species differs from all remaining species.

Focal Tier 1
------------
The singleton species is the predefined focal species.

Both directions are retained:

    focal = turnover_candidate
    others = present

or:

    focal = present
    others = turnover_candidate

Null model
----------
If each species in the clade were equally likely to carry the singleton
difference:

    expected focal share = 1 / number of species

A one-sided exact binomial test evaluates whether the focal species is
represented more frequently than expected.

This enrichment test should be interpreted as exploratory because CRE
events are not guaranteed to represent statistically independent
evolutionary events.

Outputs
-------
results/focal_tier1_candidates.tsv
results/focal_tier1_enrichment.tsv
"""

from pathlib import Path

import numpy as np
import pandas as pd

from scipy.stats import binomtest


# ============================================================
# Paths
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = SCRIPT_DIR.parent
DOWNSTREAM_DIR = ANALYSIS_DIR.parent
PROJECT_DIR = DOWNSTREAM_DIR.parent

CLASS_DIR = PROJECT_DIR / "cre_classification"

MATRIX_FILE = (
    CLASS_DIR
    / "results"
    / "cre_turnover_matrix.tsv"
)

RESULTS_DIR = ANALYSIS_DIR / "results"

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


OUT_CANDIDATES = (
    RESULTS_DIR
    / "focal_tier1_candidates.tsv"
)

OUT_ENRICHMENT = (
    RESULTS_DIR
    / "focal_tier1_enrichment.tsv"
)


# ============================================================
# Focal clades
# ============================================================

CLADES = {
    "rufa_group": {
        "display_name": "Rufa group",
        "focal": "druf",
        "species": [
            "druf",
            "dkik",
            "dbun",
            "dbir",
            "d_serrata",
        ],
    },

    "immigrans_group": {
        "display_name": "Immigrans group",
        "focal": "dimm",
        "species": [
            "dimm",
            "dfor",
            "dsul",
            "dnas",
        ],
    },

    "subobscura_group": {
        "display_name": "Subobscura group",
        "focal": "dsub",
        "species": [
            "dsub",
            "dbif",
            "dobs",
        ],
    },

    "azteca_group": {
        "display_name": "Azteca group",
        "focal": "dazt",
        "species": [
            "dazt",
            "d_affinis",
            "dmir",
            "dper",
        ],
    },

    "teissieri_group": {
        "display_name": "Teissieri group",
        "focal": "dtei",
        "species": [
            "dtei",
            "d_simulans",
            "d_lutescens",
            "dsuz",
        ],
    },

    "repleta_group": {
        "display_name": "Repleta group",
        "focal": "d_repleta",
        "species": [
            "d_repleta",
            "d_buzzatii",
            "d_mojavensis",
            "d_arizonae",
        ],
    },
}


# ============================================================
# Helpers
# ============================================================

def get_singleton_tier1(row):

    allowed = {
        "present",
        "turnover_candidate",
    }


    if not all(
        state in allowed
        for state in row.values
    ):
        return None


    counts = (
        row.value_counts()
    )


    if len(counts) != 2:
        return None


    singleton_states = (
        counts[
            counts == 1
        ]
        .index
        .tolist()
    )


    consensus_states = (
        counts[
            counts == len(row) - 1
        ]
        .index
        .tolist()
    )


    if (
        len(singleton_states) != 1
        or len(consensus_states) != 1
    ):
        return None


    discordant_state = (
        singleton_states[0]
    )


    consensus_state = (
        consensus_states[0]
    )


    discordant_species = (
        row.index[
            row == discordant_state
        ][0]
    )


    return (
        discordant_species,
        discordant_state,
        consensus_state,
    )


def benjamini_hochberg(
    p_values,
):

    p_values = np.asarray(
        p_values,
        dtype=float,
    )


    n = len(
        p_values
    )


    order = np.argsort(
        p_values
    )


    ranked = (
        p_values[
            order
        ]
    )


    adjusted_ranked = (
        ranked
        * n
        / np.arange(
            1,
            n + 1,
        )
    )


    adjusted_ranked = np.minimum.accumulate(
        adjusted_ranked[::-1]
    )[::-1]


    adjusted_ranked = np.clip(
        adjusted_ranked,
        0,
        1,
    )


    adjusted = np.empty(
        n,
        dtype=float,
    )


    adjusted[
        order
    ] = adjusted_ranked


    return adjusted


# ============================================================
# Load matrix
# ============================================================

if not MATRIX_FILE.exists():

    raise SystemExit(
        f"ERROR: matrix not found:\n{MATRIX_FILE}"
    )


matrix = pd.read_csv(
    MATRIX_FILE,
    sep="\t",
    index_col=0,
)


# ============================================================
# Analysis
# ============================================================

candidate_rows = []

enrichment_rows = []


for clade_name, config in CLADES.items():

    species = config[
        "species"
    ]

    focal = config[
        "focal"
    ]


    missing_species = [
        sp
        for sp in species
        if sp not in matrix.columns
    ]


    if missing_species:

        raise SystemExit(
            f"ERROR: species missing for {clade_name}:\n"
            + "\n".join(
                missing_species
            )
        )


    sub = matrix[
        species
    ]


    events = []


    for cre_id, row in sub.iterrows():

        result = get_singleton_tier1(
            row
        )


        if result is None:
            continue


        (
            discordant_species,
            discordant_state,
            consensus_state,
        ) = result


        is_focal = (
            discordant_species
            == focal
        )


        focal_direction = ""


        if is_focal:

            if (
                discordant_state
                == "turnover_candidate"
                and consensus_state
                == "present"
            ):

                focal_direction = (
                    "focal_turnover"
                )


            elif (
                discordant_state
                == "present"
                and consensus_state
                == "turnover_candidate"
            ):

                focal_direction = (
                    "focal_present"
                )


            else:

                raise RuntimeError(
                    "Unexpected focal Tier-1 direction."
                )


        event = {
            "cre_id": cre_id,
            "discordant_species": discordant_species,
            "discordant_state": discordant_state,
            "consensus_state": consensus_state,
            "is_focal_tier1": is_focal,
            "focal_direction": focal_direction,
        }


        events.append(
            event
        )


        candidate_rows.append({
            "clade": clade_name,
            "clade_display": config[
                "display_name"
            ],
            "dmel_cre_id": cre_id,
            "focal_species": focal,
            "discordant_species": discordant_species,
            "discordant_state": discordant_state,
            "consensus_state": consensus_state,
            "is_focal_tier1": is_focal,
            "focal_direction": focal_direction,
        })


    events_df = pd.DataFrame(
        events
    )


    n_all = len(
        events_df
    )


    if n_all > 0:

        focal_events = events_df.loc[
            events_df[
                "is_focal_tier1"
            ]
        ]


        n_focal = len(
            focal_events
        )


        n_focal_turnover = int(
            (
                focal_events[
                    "focal_direction"
                ]
                == "focal_turnover"
            ).sum()
        )


        n_focal_present = int(
            (
                focal_events[
                    "focal_direction"
                ]
                == "focal_present"
            ).sum()
        )

    else:

        n_focal = 0
        n_focal_turnover = 0
        n_focal_present = 0


    n_species = len(
        species
    )


    expected_share = (
        1.0
        / n_species
    )


    observed_share = (
        n_focal
        / n_all
        if n_all > 0
        else np.nan
    )


    enrichment_ratio = (
        observed_share
        / expected_share
        if n_all > 0
        else np.nan
    )


    if n_all > 0:

        test = binomtest(
            k=n_focal,
            n=n_all,
            p=expected_share,
            alternative="greater",
        )


        p_value = (
            test.pvalue
        )


        ci = test.proportion_ci(
            confidence_level=0.95,
            method="exact",
        )


        ci_low = (
            ci.low
        )


        ci_high = (
            ci.high
        )

    else:

        p_value = np.nan
        ci_low = np.nan
        ci_high = np.nan


    enrichment_rows.append({
        "clade": clade_name,
        "clade_display": config[
            "display_name"
        ],
        "focal_species": focal,
        "n_species": n_species,
        "n_tier1_singletons": n_all,
        "n_focal_tier1": n_focal,
        "n_focal_turnover": n_focal_turnover,
        "n_focal_present": n_focal_present,
        "observed_focal_share": observed_share,
        "expected_focal_share": expected_share,
        "enrichment_ratio": enrichment_ratio,
        "ci95_low": ci_low,
        "ci95_high": ci_high,
        "p_value": p_value,
    })


# ============================================================
# Output
# ============================================================

candidate_df = pd.DataFrame(
    candidate_rows
)


enrichment_df = pd.DataFrame(
    enrichment_rows
)


valid = (
    enrichment_df[
        "p_value"
    ].notna()
)


enrichment_df.loc[
    valid,
    "p_adj_BH",
] = benjamini_hochberg(
    enrichment_df.loc[
        valid,
        "p_value",
    ]
)


candidate_df.to_csv(
    OUT_CANDIDATES,
    sep="\t",
    index=False,
)


enrichment_df.to_csv(
    OUT_ENRICHMENT,
    sep="\t",
    index=False,
)


# ============================================================
# Console
# ============================================================

print()
print("=" * 72)
print("FOCAL TIER-1 ENRICHMENT")
print("=" * 72)
print()


columns_to_show = [
    "clade_display",
    "n_tier1_singletons",
    "n_focal_tier1",
    "n_focal_turnover",
    "n_focal_present",
    "observed_focal_share",
    "expected_focal_share",
    "enrichment_ratio",
    "p_value",
    "p_adj_BH",
]


print(
    enrichment_df[
        columns_to_show
    ].to_string(
        index=False
    )
)


print()
print(
    "NOTE: exact binomial enrichment is exploratory."
)

print()
print(
    f"Wrote:\n{OUT_ENRICHMENT}"
)

print(
    f"Wrote:\n{OUT_CANDIDATES}"
)
