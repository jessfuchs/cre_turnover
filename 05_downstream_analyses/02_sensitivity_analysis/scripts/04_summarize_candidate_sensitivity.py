#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import argparse

import pandas as pd


PROJECT_DIR = (
    Path.home()
    / "cre_turnover"
    / "project"
)


SENSITIVITY_DIR = (
    PROJECT_DIR
    / "cre_classification"
    / "sensitivity"
)


CANDIDATE_TABLE = (
    PROJECT_DIR
    / "downstream_analyses"
    / "candidate_analysis"
    / "results"
    / "tables"
    / "tier1_candidates_prioritized.tsv"
)


OUT_DIR = (
    PROJECT_DIR
    / "downstream_analyses"
    / "sensitivity_analysis"
    / "results"
)


SCENARIOS = [
    "ov025_dist12000",
    "ov025_dist24000",
    "ov025_dist48000",

    "ov050_dist12000",
    "ov050_dist24000",
    "ov050_dist48000",

    "ov075_dist12000",
    "ov075_dist24000",
    "ov075_dist48000",
]


PRIMARY_SCENARIO = (
    "ov050_dist24000"
)


POSITIVE = {
    "present",
    "turnover_candidate",
}


GROUPS = {

    "rufa_group": {

        "focal":
            "druf",

        "comparison": [
            "dkik",
            "dbun",
            "dbir",
            "d_serrata",
        ],
    },


    "immigrans_group": {

        "focal":
            "dimm",

        "comparison": [
            "dfor",
            "dsul",
            "dnas",
        ],
    },


    "obscura_group": {

        "focal":
            "dsub",

        "comparison": [
            "dbif",
            "dobs",
        ],
    },


    "azteca_affinis_miranda_group": {

        "focal":
            "dazt",

        "comparison": [
            "d_affinis",
            "dmir",
            "dper",
        ],
    },


    "teissieri_group": {

        "focal":
            "dtei",

        "comparison": [
            "d_simulans",
            "d_lutescens",
            "dsuz",
        ],
    },


    "repleta_group": {

        "focal":
            "d_repleta",

        "comparison": [
            "d_buzzatii",
            "d_mojavensis",
            "d_arizonae",
        ],
    },
}


def load_state(
    scenario,
    species,
):

    path = (
        SENSITIVITY_DIR
        / scenario
        / f"dmel_to_{species}_cre_turnover.tsv"
    )

    if not path.exists():

        raise SystemExit(
            f"ERROR: missing:\n{path}"
        )

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
    )

    if df[
        "dmel_cre_id"
    ].duplicated().any():

        raise SystemExit(
            f"ERROR: duplicate CRE IDs in {path}"
        )

    return dict(
        zip(
            df["dmel_cre_id"],
            df["class"],
        )
    )


def is_focal_tier1(
    focal_state,
    comparison_states,
):

    if len(
        set(
            comparison_states
        )
    ) != 1:

        return False

    consensus = (
        comparison_states[0]
    )

    return (
        focal_state in POSITIVE
        and
        consensus in POSITIVE
        and
        focal_state != consensus
    )


def is_singleton_tier1(
    states_by_species,
):

    states = list(
        states_by_species.values()
    )

    if any(
        state not in POSITIVE
        for state in states
    ):

        return False


    counts = Counter(
        states
    )


    if len(counts) != 2:

        return False


    expected = sorted(
        [
            1,
            len(states) - 1,
        ]
    )


    if sorted(
        counts.values()
    ) != expected:

        return False


    return True


def assign_priority(
    n_focal,
    n_secondary,
):

    if n_focal >= 2:

        return (
            "focal_recurrent"
        )


    if (
        n_focal >= 1
        and
        n_secondary >= 1
    ):

        return (
            "focal_plus_secondary_recurrent"
        )


    if (
        n_focal == 1
        and
        n_secondary == 0
    ):

        return (
            "focal_single"
        )


    if (
        n_focal == 0
        and
        n_secondary >= 2
    ):

        return (
            "secondary_recurrent"
        )


    if (
        n_focal == 0
        and
        n_secondary == 1
    ):

        return (
            "secondary_single"
        )


    return (
        "not_tier1_candidate"
    )


def main():

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )


    candidates = pd.read_csv(
        CANDIDATE_TABLE,
        sep="\t",
        dtype=str,
    )


    if "dmel_cre_id" not in candidates.columns:

        raise SystemExit(
            "ERROR: candidate table does not contain "
            "'dmel_cre_id'."
        )


    candidate_ids = (
        candidates[
            "dmel_cre_id"
        ]
        .drop_duplicates()
        .tolist()
    )


    species_needed = sorted({
        species

        for info
        in GROUPS.values()

        for species
        in [
            info["focal"],
            *info["comparison"],
        ]
    })


    cache = {}


    for scenario in SCENARIOS:

        for species in species_needed:

            cache[
                (
                    scenario,
                    species,
                )
            ] = load_state(
                scenario,
                species,
            )


    rows = []


    for cre in candidate_ids:

        for scenario in SCENARIOS:

            focal_groups = []

            secondary_groups = []


            for (
                group_name,
                info,
            ) in GROUPS.items():

                focal = (
                    info["focal"]
                )

                comparisons = (
                    info["comparison"]
                )


                species_states = {

                    species:
                        cache[
                            (
                                scenario,
                                species,
                            )
                        ].get(
                            cre,
                            "MISSING",
                        )

                    for species
                    in [
                        focal,
                        *comparisons,
                    ]
                }


                if (
                    "MISSING"
                    in species_states.values()
                ):

                    raise SystemExit(
                        f"ERROR: missing state for "
                        f"{cre} in {scenario}"
                    )


                focal_state = (
                    species_states[
                        focal
                    ]
                )


                comparison_states = [

                    species_states[
                        species
                    ]

                    for species
                    in comparisons
                ]


                focal_tier1 = (
                    is_focal_tier1(
                        focal_state,
                        comparison_states,
                    )
                )


                singleton = (
                    is_singleton_tier1(
                        species_states
                    )
                )


                secondary_only = (
                    singleton
                    and
                    not focal_tier1
                )


                if focal_tier1:

                    focal_groups.append(
                        group_name
                    )


                if secondary_only:

                    secondary_groups.append(
                        group_name
                    )


            n_focal = (
                len(
                    focal_groups
                )
            )


            n_secondary = (
                len(
                    secondary_groups
                )
            )


            priority = (
                assign_priority(
                    n_focal,
                    n_secondary,
                )
            )


            rows.append({

                "dmel_cre_id":
                    cre,

                "scenario":
                    scenario,

                "n_focal_tier1_clades":
                    n_focal,

                "focal_tier1_clades":
                    "|".join(
                        focal_groups
                    ),

                "n_secondary_only_tier1_clades":
                    n_secondary,

                "secondary_only_tier1_clades":
                    "|".join(
                        secondary_groups
                    ),

                "n_total_tier1_clades":
                    n_focal
                    + n_secondary,

                "candidate_priority":
                    priority,

                "candidate_retained":
                    (
                        "yes"
                        if priority
                        != "not_tier1_candidate"
                        else "no"
                    ),
            })


    long = pd.DataFrame(
        rows
    )


    # --------------------------------------------------------
    # Primary priorities
    # --------------------------------------------------------

    primary = (
        long[
            long["scenario"]
            == PRIMARY_SCENARIO
        ][
            [
                "dmel_cre_id",
                "candidate_priority",
            ]
        ]
        .rename(
            columns={
                "candidate_priority":
                    "primary_candidate_priority"
            }
        )
    )


    long = long.merge(
        primary,
        on="dmel_cre_id",
        how="left",
        validate="many_to_one",
    )


    long[
        "same_priority_as_primary"
    ] = (

        long[
            "candidate_priority"
        ]

        ==

        long[
            "primary_candidate_priority"
        ]

    ).map({
        True: "yes",
        False: "no",
    })


    # --------------------------------------------------------
    # Candidate summary
    # --------------------------------------------------------

    summary_rows = []


    for (
        cre,
        g,
    ) in long.groupby(
        "dmel_cre_id"
    ):

        retained = (
            g[
                "candidate_retained"
            ]
            == "yes"
        )


        same_priority = (
            g[
                "same_priority_as_primary"
            ]
            == "yes"
        )


        lost = (
            g.loc[
                ~retained,
                "scenario",
            ]
            .tolist()
        )


        changed = g.loc[
            ~same_priority,
            [
                "scenario",
                "candidate_priority",
            ],
        ]


        summary_rows.append({

            "dmel_cre_id":
                cre,

            "primary_candidate_priority":
                g[
                    "primary_candidate_priority"
                ].iloc[0],

            "n_scenarios":
                len(g),

            "n_candidate_retained":
                int(
                    retained.sum()
                ),

            "percent_candidate_retained":
                100
                * retained.mean(),

            "robust_candidate_all_9":
                (
                    "yes"
                    if retained.all()
                    else "no"
                ),

            "n_same_priority_as_primary":
                int(
                    same_priority.sum()
                ),

            "percent_same_priority_as_primary":
                100
                * same_priority.mean(),

            "robust_priority_all_9":
                (
                    "yes"
                    if same_priority.all()
                    else "no"
                ),

            "lost_in_scenarios":
                (
                    "|".join(
                        lost
                    )
                    if lost
                    else "NA"
                ),

            "changed_priority_in_scenarios":
                (
                    "|".join(
                        f"{row.scenario}:"
                        f"{row.candidate_priority}"

                        for row
                        in changed.itertuples()
                    )

                    if not changed.empty

                    else "NA"
                ),
        })


    summary = (
        pd.DataFrame(
            summary_rows
        )
        .sort_values(
            [
                "percent_candidate_retained",
                "percent_same_priority_as_primary",
                "dmel_cre_id",
            ],
            ascending=[
                False,
                False,
                True,
            ],
        )
    )


    # --------------------------------------------------------
    # Priority matrix
    # --------------------------------------------------------

    priority_matrix = (
        long
        .pivot(
            index="dmel_cre_id",
            columns="scenario",
            values="candidate_priority",
        )
        .reset_index()
    )


    priority_matrix = (
        priority_matrix[
            [
                "dmel_cre_id",
                *SCENARIOS,
            ]
        ]
    )


    # --------------------------------------------------------
    # Retention matrix
    # --------------------------------------------------------

    retention_matrix = (
        long
        .pivot(
            index="dmel_cre_id",
            columns="scenario",
            values="candidate_retained",
        )
        .reset_index()
    )


    retention_matrix = (
        retention_matrix[
            [
                "dmel_cre_id",
                *SCENARIOS,
            ]
        ]
    )


    # --------------------------------------------------------
    # Write
    # --------------------------------------------------------

    long.to_csv(
        OUT_DIR
        / "candidate_sensitivity_long.tsv",
        sep="\t",
        index=False,
    )


    summary.to_csv(
        OUT_DIR
        / "candidate_sensitivity_summary.tsv",
        sep="\t",
        index=False,
    )


    priority_matrix.to_csv(
        OUT_DIR
        / "candidate_sensitivity_priority_matrix.tsv",
        sep="\t",
        index=False,
    )


    retention_matrix.to_csv(
        OUT_DIR
        / "candidate_sensitivity_retention_matrix.tsv",
        sep="\t",
        index=False,
    )


    # --------------------------------------------------------
    # Console report
    # --------------------------------------------------------

    print()
    print("=" * 72)
    print("FINAL TIER-1 CANDIDATE SENSITIVITY")
    print("=" * 72)
    print()

    print(
        f"Candidates evaluated: "
        f"{len(candidate_ids)}"
    )

    print(
        f"Sensitivity scenarios: "
        f"{len(SCENARIOS)}"
    )

    print()

    print(
        "Primary candidate priorities:"
    )

    print(
        primary[
            "primary_candidate_priority"
        ]
        .value_counts()
        .to_string()
    )

    print()

    print(
        "Candidate retained in all 9 scenarios:"
    )

    print(
        summary[
            "robust_candidate_all_9"
        ]
        .value_counts()
        .to_string()
    )

    print()

    print(
        "Exact candidate priority retained in all 9 scenarios:"
    )

    print(
        summary[
            "robust_priority_all_9"
        ]
        .value_counts()
        .to_string()
    )

    print()

    print(
        summary[
            [
                "dmel_cre_id",
                "primary_candidate_priority",
                "percent_candidate_retained",
                "percent_same_priority_as_primary",
                "robust_candidate_all_9",
                "robust_priority_all_9",
            ]
        ].to_string(
            index=False,
            float_format=
                lambda x:
                    f"{x:.1f}",
        )
    )


if __name__ == "__main__":
    main()
