#!/usr/bin/env python3

"""
    Classify the CRE-state pattern for one focal clade.

    Categories
    ----------
    tier1
        Focal and comparison consensus differ between
        present and turnover_candidate.

        Example:
            focal: turnover_candidate
            comparisons: present, present, present

        This is the strongest candidate for lineage-specific
        positional CRE turnover.

    tier2
        Focal and comparison consensus differ between
        a positive CRE state and no_detected_CRE.

        Example:
            focal: no_detected_CRE
            comparisons: present, present, present

        This represents a detection contrast and is interpreted
        more conservatively.

    tier3
        Comparison species agree and focal differs, but the
        contrast involves uncertain or another non-priority
        state combination.

    comparison_mixed
        Comparison species do not share one common state.

    all_same
        Focal and all comparison species share the same state.

    invalid
        At least one value is outside the expected CRE-state
        vocabulary.
"""

from pathlib import Path
import pandas as pd


# ============================================================
# Paths
# ============================================================

PROJECT_DIR = (
    Path.home()
    / "cre_turnover"
    / "project"
)

CLASS_DIR = (
    PROJECT_DIR
    / "cre_classification"
)

CANDIDATE_DIR = (
    PROJECT_DIR
    / "downstream_analyses"
    / "candidate_analysis"
)

RESULTS_DIR = (
    CANDIDATE_DIR
    / "results"
)

OUTDIR = (
    RESULTS_DIR
    / "focal_clades"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)

MATRIX_FILE = (
    CLASS_DIR
    / "results"
    / "cre_turnover_matrix.tsv"
)

SUMMARY_FILE = (
    RESULTS_DIR
    / "focal_clade_summary.tsv"
)


# ============================================================
# Each group contains:
# - one focal species from one climatic zone
# - phylogenetically close comparison species sharing another
#   climatic zone
#
# D. melanogaster is not included because it is the reference
# species and therefore has no target-state column in the matrix
# ============================================================

GROUPS = {
    "rufa_group": {
        "focal": "druf",              # TEMP
        "comparison": [
            "dkik",                   # TROP
            "dbun",                   # TROP
            "dbir",                   # TROP
            "d_serrata",              # TROP
        ],
        "focal_climate": "TEMP",
        "comparison_climate": "TROP",
    },

    "immigrans_group": {
        "focal": "dimm",              # TEMP
        "comparison": [
            "dfor",                   # TROP
            "dsul",                   # TROP
            "dnas",                   # TROP
        ],
        "focal_climate": "TEMP",
        "comparison_climate": "TROP",
    },

    "obscura_group": {
        "focal": "dsub",              # TEMP
        "comparison": [
            "dbif",                   # BORE
            "dobs",                   # BORE
        ],
        "focal_climate": "TEMP",
        "comparison_climate": "BORE",
    },

    "azteca_affinis_miranda_group": {
        "focal": "dazt",              # ARID
        "comparison": [
            "d_affinis",              # TEMP
            "dmir",                   # TEMP
            "dper",                   # TEMP
        ],
        "focal_climate": "ARID",
        "comparison_climate": "TEMP",
    },

    "teissieri_group": {
        "focal": "dtei",              # TROP
        "comparison": [
            "d_simulans",             # TEMP
            "d_lutescens",            # TEMP
            "dsuz",                   # TEMP
        ],
        "focal_climate": "TROP",
        "comparison_climate": "TEMP",
    },

    "repleta_group": {
        "focal": "d_repleta",         # TEMP
        "comparison": [
            "d_buzzatii",             # ARID
            "d_mojavensis",           # ARID
            "d_arizonae",             # ARID
        ],
        "focal_climate": "TEMP",
        "comparison_climate": "ARID",
    },
}


# ============================================================
# CRE-state definitions
# ============================================================

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
# Helper functions
# ============================================================

def require_file(path):
    """
    Abort if a required input file does not exist.
    """

    if not path.exists():
        raise SystemExit(
            "ERROR: required input file not found:\n"
            f"{path}"
        )


def classify_focal(
    row,
    focal,
    comparisons,
):

    focal_state = row[focal]

    comparison_states = [
        row[sp]
        for sp in comparisons
    ]

    all_states = [
        focal_state,
        *comparison_states,
    ]

    # --------------------------------------------------------
    # Validate state vocabulary
    # --------------------------------------------------------

    if not set(all_states).issubset(
        VALID_STATES
    ):
        return "invalid"

    # --------------------------------------------------------
    # Comparison species must agree
    # --------------------------------------------------------

    if len(
        set(comparison_states)
    ) != 1:
        return "comparison_mixed"

    consensus_state = (
        comparison_states[0]
    )

    # --------------------------------------------------------
    # No focal difference
    # --------------------------------------------------------

    if focal_state == consensus_state:
        return "all_same"

    # --------------------------------------------------------
    # Tier 1:
    # present <-> turnover_candidate
    # --------------------------------------------------------

    if (
        focal_state in POSITIVE_STATES
        and
        consensus_state in POSITIVE_STATES
        and
        focal_state != consensus_state
    ):
        return "tier1"

    # --------------------------------------------------------
    # Tier 2:
    # positive state <-> no_detected_CRE
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # Tier 3:
    # remaining focal contrasts, including uncertain
    # --------------------------------------------------------

    return "tier3"


def get_comparison_consensus(
    row,
    comparisons,
):
    """
    Return the common comparison-species state.

    Returns 'mixed' when comparison species disagree.
    """

    values = [
        row[sp]
        for sp in comparisons
    ]

    unique = set(values)

    if len(unique) == 1:
        return values[0]

    return "mixed"


# ============================================================
# Input checks
# ============================================================

require_file(
    MATRIX_FILE
)

df = pd.read_csv(
    MATRIX_FILE,
    sep="\t",
    dtype=str,
).fillna("")


# ============================================================
# Validate matrix structure
# ============================================================

if "dmel_cre_id" not in df.columns:
    raise SystemExit(
        "ERROR: cre_turnover_matrix.tsv does not contain "
        "a dmel_cre_id column."
    )


# ------------------------------------------------------------
# Duplicated reference CRE IDs
# ------------------------------------------------------------

if df["dmel_cre_id"].duplicated().any():

    duplicated = (
        df.loc[
            df[
                "dmel_cre_id"
            ].duplicated(
                keep=False
            ),
            "dmel_cre_id",
        ]
        .unique()
        .tolist()
    )

    raise SystemExit(
        "ERROR: duplicated D. melanogaster CRE IDs "
        "in cre_turnover_matrix.tsv:\n"
        + "\n".join(
            sorted(duplicated)
        )
    )


# ------------------------------------------------------------
# Validate required focal-clade species
# ------------------------------------------------------------

required_species = {
    sp
    for info in GROUPS.values()
    for sp in [
        info["focal"],
        *info["comparison"],
    ]
}

missing_species = sorted(
    required_species
    - set(df.columns)
)

if missing_species:
    raise SystemExit(
        "ERROR: focal-clade species missing from "
        "cre_turnover_matrix.tsv:\n"
        + "\n".join(
            missing_species
        )
    )


# ============================================================
# Validate CRE-state vocabulary
# ============================================================

state_columns = sorted(
    required_species
)

observed_states = {
    value
    for value in pd.unique(
        df[
            state_columns
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

    comparisons = (
        info["comparison"]
    )

    focal_climate = (
        info["focal_climate"]
    )

    comparison_climate = (
        info["comparison_climate"]
    )

    species = [
        focal,
        *comparisons,
    ]

    # --------------------------------------------------------
    # Extract group-specific matrix
    # --------------------------------------------------------

    sub = df[
        [
            "dmel_cre_id",
            *species,
        ]
    ].copy()

    # --------------------------------------------------------
    # Classify each reference CRE
    # --------------------------------------------------------

    sub["category"] = sub.apply(
        lambda row: classify_focal(
            row,
            focal,
            comparisons,
        ),
        axis=1,
    )

    # --------------------------------------------------------
    # Add explicit analysis metadata
    # --------------------------------------------------------

    sub["group_name"] = (
        group_name
    )

    sub["focal_species"] = (
        focal
    )

    sub["focal_climate"] = (
        focal_climate
    )

    sub["focal_state"] = (
        sub[focal]
    )

    sub[
        "comparison_species"
    ] = "|".join(
        comparisons
    )

    sub[
        "comparison_climate"
    ] = (
        comparison_climate
    )

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
    # Useful column order
    # --------------------------------------------------------

    front_columns = [
        "group_name",
        "dmel_cre_id",
        "category",
        "focal_species",
        "focal_climate",
        "focal_state",
        "comparison_species",
        "comparison_climate",
        "comparison_consensus_state",
        *species,
    ]

    remaining_columns = [
        col
        for col in sub.columns
        if col not in front_columns
    ]

    sub = sub[
        front_columns
        + remaining_columns
    ]

    # --------------------------------------------------------
    # Write complete group table
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

    # --------------------------------------------------------
    # Write one table per category
    # --------------------------------------------------------

    for category in OUTPUT_CATEGORIES:

        out = sub[
            sub["category"]
            == category
        ].copy()

        out_file = (
            OUTDIR
            / (
                f"{group_name}_"
                f"{category}.tsv"
            )
        )

        out.to_csv(
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

    n_tier1 = int(
        counts.get(
            "tier1",
            0,
        )
    )

    n_tier2 = int(
        counts.get(
            "tier2",
            0,
        )
    )

    n_tier3 = int(
        counts.get(
            "tier3",
            0,
        )
    )

    n_comparison_mixed = int(
        counts.get(
            "comparison_mixed",
            0,
        )
    )

    n_all_same = int(
        counts.get(
            "all_same",
            0,
        )
    )

    n_invalid = int(
        counts.get(
            "invalid",
            0,
        )
    )

    n_reference_cres = len(
        sub
    )

    # --------------------------------------------------------
    # Sanity check:
    # all categories must sum to all reference CREs
    # --------------------------------------------------------

    n_total_classified = (
        n_tier1
        + n_tier2
        + n_tier3
        + n_comparison_mixed
        + n_all_same
        + n_invalid
    )

    if (
        n_total_classified
        != n_reference_cres
    ):
        raise SystemExit(
            "ERROR: category counts do not sum to "
            f"all CREs for {group_name}.\n"
            f"Expected: {n_reference_cres}\n"
            f"Observed: {n_total_classified}"
        )

    # --------------------------------------------------------
    # Summary row
    # --------------------------------------------------------

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
            n_tier1,

        "n_tier2":
            n_tier2,

        "n_tier3":
            n_tier3,

        "n_comparison_mixed":
            n_comparison_mixed,

        "n_all_same":
            n_all_same,

        "n_invalid":
            n_invalid,
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
        f"{n_tier1}"
    )

    print(
        "Tier 2 "
        "(positive vs no_detected_CRE): "
        f"{n_tier2}"
    )

    print(
        "Tier 3 "
        "(other focal contrast): "
        f"{n_tier3}"
    )

    print(
        "Comparison mixed: "
        f"{n_comparison_mixed}"
    )

    print(
        "All same: "
        f"{n_all_same}"
    )

    print(
        "Invalid: "
        f"{n_invalid}"
    )

    # --------------------------------------------------------
    # Print Tier 1 candidates
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
            "dmel_cre_id",
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
# Cross-clade summary
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
# Final QC
# ============================================================

if summary["n_invalid"].sum() != 0:
    print()
    print(
        "WARNING: invalid CRE-state entries were detected."
    )


# ============================================================
# Final console summary
# ============================================================

print()
print("=" * 72)
print("Focal-clade candidate analysis complete")
print("=" * 72)

print(
    f"Input matrix:\n"
    f"{MATRIX_FILE}"
)

print()

print(
    f"Reference CREs: "
    f"{len(df)}"
)

print(
    f"Focal clades analyzed: "
    f"{len(GROUPS)}"
)

print()

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

print()

print(
    f"Wrote clade summary:\n"
    f"{SUMMARY_FILE}"
)
