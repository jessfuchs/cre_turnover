#!/usr/bin/env python3

# ============================================================
# Analyze focal-clade CRE-state contrasts
#
# Purpose:
#   Identify lineage-specific CRE-state contrasts between one
#   focal species and phylogenetically close comparison species.
#
# Categories:
#   - tier1: present <-> turnover_candidate
#   - tier2: positive CRE state <-> no_detected_CRE
#   - tier3: other focal contrasts with comparison consensus
#   - comparison_mixed: comparison species do not agree
#   - all_same: focal and comparison species share one state
#   - invalid: unexpected CRE-state value
#
# Input/output paths and focal-clade definitions:
#   Supplied by the pipeline wrapper using
#   config/candidate_config.sh.
# ============================================================

import argparse
from pathlib import Path

import pandas as pd


# ============================================================
# CRE-state definitions
# ============================================================

REFERENCE_CRE_ID_COLUMN = "dmel_cre_id"

POSITIVE_STATES = {
    "present",
    "turnover_candidate",
}

VALID_STATES = {
    "present",
    "turnover_candidate",
    "no_detected_CRE",
    "uncertain",
}

OUTPUT_CATEGORIES = [
    "tier1",
    "tier2",
    "tier3",
    "comparison_mixed",
    "all_same",
    "invalid",
]


# ============================================================
# Arguments
# ============================================================

parser = argparse.ArgumentParser(
    description=(
        "Analyze focal-clade CRE-state contrasts between one "
        "focal species and phylogenetically close comparison "
        "species."
    )
)

parser.add_argument(
    "--matrix",
    required=True,
)

parser.add_argument(
    "--groups",
    required=True,
)

parser.add_argument(
    "--manifest",
    required=True,
)

parser.add_argument(
    "--traits",
    required=True,
)

parser.add_argument(
    "--out-dir",
    required=True,
)

parser.add_argument(
    "--summary-out",
    required=True,
)

parser.add_argument(
    "--reference-species",
    default="d_melanogaster",
)

parser.add_argument(
    "--expected-reference-cres",
    type=int,
    default=None,
)

args = parser.parse_args()


# ============================================================
# Input and output paths
# ============================================================

MATRIX_FILE = Path(args.matrix)
GROUPS_FILE = Path(args.groups)
MANIFEST_FILE = Path(args.manifest)
TRAITS_FILE = Path(args.traits)

OUTDIR = Path(args.out_dir)
SUMMARY_FILE = Path(args.summary_out)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)

SUMMARY_FILE.parent.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Helper functions
# ============================================================

def require_file(path):
    """
    Abort if a required input file is missing.
    """

    if not path.is_file():
        raise SystemExit(
            f"ERROR: required input file not found:\n{path}"
        )


def tree_name_from_species(species_name):
    """
    Convert a scientific species name to the tree-name
    convention used in the published phylogeny.
    """

    return (
        str(species_name)
        .strip()
        .upper()
        .replace(" ", "_")
    )


def load_groups(path):
    """
    Load focal-clade definitions.
    """

    groups = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
    ).fillna("")

    required = {
        "group_name",
        "focal_species",
        "comparison_species",
    }

    missing = (
        required
        - set(groups.columns)
    )

    if missing:
        raise SystemExit(
            "ERROR: focal-clade definition file is missing "
            "columns:\n"
            + "\n".join(
                sorted(missing)
            )
        )

    for column in required:
        groups[column] = (
            groups[column]
            .astype(str)
            .str.strip()
        )

    if (
        groups[
            list(required)
        ] == ""
    ).any().any():
        raise SystemExit(
            "ERROR: focal-clade definition file contains "
            "empty required values."
        )

    duplicated = (
        groups.loc[
            groups["group_name"].duplicated(
                keep=False
            ),
            "group_name",
        ]
        .unique()
        .tolist()
    )

    if duplicated:
        raise SystemExit(
            "ERROR: duplicate focal-clade names:\n"
            + "\n".join(
                sorted(duplicated)
            )
        )

    result = {}

    for row in groups.itertuples(
        index=False
    ):

        comparisons = [
            value.strip()
            for value
            in row.comparison_species.split("|")
            if value.strip()
        ]

        if not comparisons:
            raise SystemExit(
                "ERROR: no comparison species defined for "
                f"{row.group_name}."
            )

        if len(comparisons) != len(
            set(comparisons)
        ):
            raise SystemExit(
                "ERROR: duplicate comparison species in "
                f"{row.group_name}."
            )

        if row.focal_species in comparisons:
            raise SystemExit(
                "ERROR: focal species also occurs among "
                f"comparisons for {row.group_name}."
            )

        result[row.group_name] = {
            "focal": row.focal_species,
            "comparison": comparisons,
        }

    return result


def load_species_climates(
    manifest_file,
    traits_file,
):
    """
    Map pipeline species slugs to climatic zones.
    """

    manifest = pd.read_csv(
        manifest_file,
        sep="\t",
        dtype=str,
    ).fillna("")

    required_manifest = {
        "slug",
        "species",
    }

    missing = (
        required_manifest
        - set(manifest.columns)
    )

    if missing:
        raise SystemExit(
            "ERROR: combined manifest is missing columns:\n"
            + "\n".join(
                sorted(missing)
            )
        )

    for column in required_manifest:
        manifest[column] = (
            manifest[column]
            .astype(str)
            .str.strip()
        )

    duplicated_slugs = (
        manifest.loc[
            manifest["slug"].duplicated(
                keep=False
            ),
            "slug",
        ]
        .unique()
        .tolist()
    )

    if duplicated_slugs:
        raise SystemExit(
            "ERROR: duplicate species slugs in combined "
            "manifest:\n"
            + "\n".join(
                sorted(duplicated_slugs)
            )
        )

    manifest["tree_name"] = (
        manifest["species"]
        .map(
            tree_name_from_species
        )
    )

    traits = pd.read_csv(
        traits_file,
        sep="\t",
        dtype=str,
    ).fillna("")

    required_traits = {
        "tree_name",
        "climatic_zone",
    }

    missing = (
        required_traits
        - set(traits.columns)
    )

    if missing:
        raise SystemExit(
            "ERROR: species trait table is missing columns:\n"
            + "\n".join(
                sorted(missing)
            )
        )

    for column in required_traits:
        traits[column] = (
            traits[column]
            .astype(str)
            .str.strip()
        )

    duplicated_tree_names = (
        traits.loc[
            traits["tree_name"].duplicated(
                keep=False
            ),
            "tree_name",
        ]
        .unique()
        .tolist()
    )

    if duplicated_tree_names:
        raise SystemExit(
            "ERROR: duplicate tree names in species trait "
            "table:\n"
            + "\n".join(
                sorted(duplicated_tree_names)
            )
        )

    tree_to_climate = dict(
        zip(
            traits["tree_name"],
            traits["climatic_zone"],
        )
    )

    return {
        row.slug: tree_to_climate.get(
            row.tree_name,
            "",
        )
        for row in manifest.itertuples(
            index=False
        )
    }


def classify_focal(
    row,
    focal,
    comparisons,
):
    """
    Classify one focal-clade CRE-state pattern.
    """

    focal_state = row[focal]

    comparison_states = [
        row[species]
        for species in comparisons
    ]

    all_states = [
        focal_state,
        *comparison_states,
    ]

    if not set(
        all_states
    ).issubset(
        VALID_STATES
    ):
        return "invalid"

    if len(
        set(comparison_states)
    ) != 1:
        return "comparison_mixed"

    consensus_state = (
        comparison_states[0]
    )

    if focal_state == consensus_state:
        return "all_same"

    # Tier 1:
    # present <-> turnover_candidate

    if (
        focal_state in POSITIVE_STATES
        and consensus_state in POSITIVE_STATES
    ):
        return "tier1"

    # Tier 2:
    # positive state <-> no_detected_CRE

    if (
        (
            focal_state in POSITIVE_STATES
            and
            consensus_state == "no_detected_CRE"
        )
        or
        (
            focal_state == "no_detected_CRE"
            and
            consensus_state in POSITIVE_STATES
        )
    ):
        return "tier2"

    # Tier 3:
    # remaining focal contrasts, including uncertain.

    return "tier3"


def get_comparison_consensus(
    row,
    comparisons,
):
    """
    Return the common comparison state or 'mixed'.
    """

    states = [
        row[species]
        for species in comparisons
    ]

    unique = set(states)

    if len(unique) == 1:
        return states[0]

    return "mixed"


# ============================================================
# Load inputs
# ============================================================

for path in [
    MATRIX_FILE,
    GROUPS_FILE,
    MANIFEST_FILE,
    TRAITS_FILE,
]:
    require_file(path)


GROUPS = load_groups(
    GROUPS_FILE
)

slug_to_climate = (
    load_species_climates(
        MANIFEST_FILE,
        TRAITS_FILE,
    )
)

df = pd.read_csv(
    MATRIX_FILE,
    sep="\t",
    dtype=str,
).fillna("")


# ============================================================
# Validate CRE-state matrix
# ============================================================

if REFERENCE_CRE_ID_COLUMN not in df.columns:
    raise SystemExit(
        "ERROR: CRE-state matrix does not contain "
        f"'{REFERENCE_CRE_ID_COLUMN}'."
    )


if df[
    REFERENCE_CRE_ID_COLUMN
].duplicated().any():

    duplicated = (
        df.loc[
            df[
                REFERENCE_CRE_ID_COLUMN
            ].duplicated(
                keep=False
            ),
            REFERENCE_CRE_ID_COLUMN,
        ]
        .unique()
        .tolist()
    )

    raise SystemExit(
        "ERROR: duplicated D. melanogaster CRE IDs:\n"
        + "\n".join(
            sorted(duplicated)
        )
    )


if (
    args.expected_reference_cres is not None
    and
    len(df) != args.expected_reference_cres
):
    raise SystemExit(
        "ERROR: unexpected number of reference CREs.\n"
        f"Expected: {args.expected_reference_cres}\n"
        f"Observed: {len(df)}"
    )


# ============================================================
# Validate focal-clade species
# ============================================================

required_species = {
    species
    for info in GROUPS.values()
    for species in [
        info["focal"],
        *info["comparison"],
    ]
}


if args.reference_species in required_species:
    raise SystemExit(
        "ERROR: the reference species must not be included "
        "in a focal-clade definition."
    )


missing_matrix_species = sorted(
    required_species
    - set(df.columns)
)

if missing_matrix_species:
    raise SystemExit(
        "ERROR: focal-clade species missing from "
        "CRE-state matrix:\n"
        + "\n".join(
            missing_matrix_species
        )
    )


missing_climates = sorted(
    species
    for species in required_species
    if not slug_to_climate.get(
        species,
        "",
    )
)

if missing_climates:
    raise SystemExit(
        "ERROR: climatic-zone annotation missing for:\n"
        + "\n".join(
            missing_climates
        )
    )


# ============================================================
# Validate CRE-state vocabulary
# ============================================================

observed_states = {
    value
    for value in pd.unique(
        df[
            sorted(required_species)
        ].values.ravel()
    )
    if value != ""
}

unknown_states = sorted(
    observed_states
    - VALID_STATES
)

if unknown_states:
    raise SystemExit(
        "ERROR: unexpected CRE-state values in focal-clade "
        "species:\n"
        + "\n".join(
            unknown_states
        )
    )


# ============================================================
# Analyze focal clades
# ============================================================

summary_rows = []


for group_name, info in GROUPS.items():

    focal = info["focal"]
    comparisons = info["comparison"]

    focal_climate = (
        slug_to_climate[focal]
    )

    comparison_climates = {
        slug_to_climate[species]
        for species in comparisons
    }

    if len(
        comparison_climates
    ) != 1:
        raise SystemExit(
            "ERROR: comparison species in "
            f"{group_name} do not share one climatic zone:\n"
            + "\n".join(
                f"{species}\t"
                f"{slug_to_climate[species]}"
                for species in comparisons
            )
        )

    comparison_climate = next(
        iter(
            comparison_climates
        )
    )

    if focal_climate == comparison_climate:
        raise SystemExit(
            f"ERROR: focal and comparison species in "
            f"{group_name} share the same climatic zone "
            f"({focal_climate})."
        )

    species = [
        focal,
        *comparisons,
    ]


    # --------------------------------------------------------
    # Extract and classify group-specific CRE states
    # --------------------------------------------------------

    sub = df[
        [
            REFERENCE_CRE_ID_COLUMN,
            *species,
        ]
    ].copy()

    sub["category"] = sub.apply(
        lambda row: classify_focal(
            row,
            focal,
            comparisons,
        ),
        axis=1,
    )


    # --------------------------------------------------------
    # Add analysis metadata
    # --------------------------------------------------------

    sub["group_name"] = group_name
    sub["focal_species"] = focal
    sub["focal_climate"] = focal_climate
    sub["focal_state"] = sub[focal]

    sub[
        "comparison_species"
    ] = "|".join(
        comparisons
    )

    sub[
        "comparison_climate"
    ] = comparison_climate

    sub[
        "comparison_consensus_state"
    ] = sub.apply(
        lambda row:
        get_comparison_consensus(
            row,
            comparisons,
        ),
        axis=1,
    )


    # --------------------------------------------------------
    # Column order
    # --------------------------------------------------------

    columns = [
        "group_name",
        REFERENCE_CRE_ID_COLUMN,
        "category",
        "focal_species",
        "focal_climate",
        "focal_state",
        "comparison_species",
        "comparison_climate",
        "comparison_consensus_state",
        *species,
    ]

    sub = sub[
        columns
    ]


    # --------------------------------------------------------
    # Write group tables
    # --------------------------------------------------------

    all_file = (
        OUTDIR
        / f"{group_name}_all.tsv"
    )

    sub.to_csv(
        all_file,
        sep="\t",
        index=False,
    )


    for category in OUTPUT_CATEGORIES:

        out_file = (
            OUTDIR
            / f"{group_name}_{category}.tsv"
        )

        sub.loc[
            sub["category"] == category
        ].to_csv(
            out_file,
            sep="\t",
            index=False,
        )


    # --------------------------------------------------------
    # Group-level counts
    # --------------------------------------------------------

    counts = (
        sub["category"]
        .value_counts()
    )

    category_counts = {
        category: int(
            counts.get(
                category,
                0,
            )
        )
        for category in OUTPUT_CATEGORIES
    }

    n_reference_cres = len(
        sub
    )

    if (
        sum(
            category_counts.values()
        )
        != n_reference_cres
    ):
        raise SystemExit(
            "ERROR: category counts do not sum to all CREs "
            f"for {group_name}."
        )


    summary_rows.append({
        "group_name":
            group_name,

        "focal_species":
            focal,

        "focal_climate":
            focal_climate,

        "comparison_species":
            "|".join(
                comparisons
            ),

        "comparison_climate":
            comparison_climate,

        "n_comparison_species":
            len(
                comparisons
            ),

        "n_reference_cres":
            n_reference_cres,

        "n_tier1":
            category_counts[
                "tier1"
            ],

        "n_tier2":
            category_counts[
                "tier2"
            ],

        "n_tier3":
            category_counts[
                "tier3"
            ],

        "n_comparison_mixed":
            category_counts[
                "comparison_mixed"
            ],

        "n_all_same":
            category_counts[
                "all_same"
            ],

        "n_invalid":
            category_counts[
                "invalid"
            ],
    })


    # --------------------------------------------------------
    # Console summary
    # --------------------------------------------------------

    print()
    print("=" * 72)
    print(group_name)
    print("=" * 72)

    print(
        f"Focal species: "
        f"{focal} ({focal_climate})"
    )

    print(
        "Comparison species: "
        + ", ".join(
            comparisons
        )
    )

    print(
        f"Comparison climate: "
        f"{comparison_climate}"
    )

    print()

    print(
        "Tier 1 "
        "(present vs turnover_candidate): "
        f"{category_counts['tier1']}"
    )

    print(
        "Tier 2 "
        "(positive vs no_detected_CRE): "
        f"{category_counts['tier2']}"
    )

    print(
        "Tier 3 "
        "(other focal contrast): "
        f"{category_counts['tier3']}"
    )

    print(
        "Comparison mixed: "
        f"{category_counts['comparison_mixed']}"
    )

    print(
        "All same: "
        f"{category_counts['all_same']}"
    )

    print(
        "Invalid: "
        f"{category_counts['invalid']}"
    )


    # --------------------------------------------------------
    # Print Tier-1 candidates
    # --------------------------------------------------------

    tier1 = sub[
        sub["category"]
        == "tier1"
    ]

    if not tier1.empty:

        print()
        print(
            "Tier 1 focal candidates:"
        )

        display_columns = [
            REFERENCE_CRE_ID_COLUMN,
            "focal_state",
            "comparison_consensus_state",
            *species,
        ]

        print(
            tier1[
                display_columns
            ].to_string(
                index=False
            )
        )


# ============================================================
# Write cross-clade summary
# ============================================================

summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    SUMMARY_FILE,
    sep="\t",
    index=False,
)


# ============================================================
# Final QC and report
# ============================================================

if summary[
    "n_invalid"
].sum() != 0:

    print()
    print(
        "WARNING: invalid CRE-state entries were detected."
    )


print()
print("=" * 72)
print(
    "Focal-clade candidate analysis complete"
)
print("=" * 72)

print(
    f"Reference CREs: "
    f"{len(df)}"
)

print(
    f"Focal clades analyzed: "
    f"{len(GROUPS)}"
)

print(
    "Total Tier 1 candidates across clades "
    "(not deduplicated by CRE): "
    f"{summary['n_tier1'].sum()}"
)

print(
    "Total Tier 2 candidates across clades "
    "(not deduplicated by CRE): "
    f"{summary['n_tier2'].sum()}"
)

print()

print(
    f"Wrote focal-clade tables to:\n"
    f"{OUTDIR}"
)

print(
    f"Wrote clade summary:\n"
    f"{SUMMARY_FILE}"
)
