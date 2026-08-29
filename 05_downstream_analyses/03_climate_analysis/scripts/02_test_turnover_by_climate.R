#!/usr/bin/env Rscript

# ============================================================
# Phylogenetically controlled climate analysis
#
# Response:
#   empirical-logit transformed turnover_rate_evaluable
#
# Model:
#   logit turnover rate ~ climatic_zone
#
# Phylogenetic covariance:
#   Pagel's lambda
#
# Input:
#   results/climate_turnover_summary.tsv
#
# Outputs:
#   results/pgls_climate_model_comparison.tsv
#   results/pgls_climate_coefficients.tsv
#   results/pgls_climate_summary.txt
# ============================================================


# ============================================================
# Packages
# ============================================================

required_packages <- c(
    "ape",
    "nlme"
)


missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]


if (length(missing_packages) > 0) {

    stop(
        paste0(
            "Missing R packages: ",
            paste(
                missing_packages,
                collapse = ", "
            ),
            "\nInstall with:\n",
            "install.packages(c(",
            paste(
                sprintf(
                    "\"%s\"",
                    missing_packages
                ),
                collapse = ", "
            ),
            "))"
        )
    )
}


library(ape)
library(nlme)


# ============================================================
# Paths
# ============================================================

args <- commandArgs(
    trailingOnly = FALSE
)


file_arg <- grep(
    "^--file=",
    args,
    value = TRUE
)


if (length(file_arg) != 1) {

    stop(
        "Could not determine script path."
    )
}


script_file <- sub(
    "^--file=",
    "",
    file_arg
)


SCRIPT_DIR <- dirname(
    normalizePath(
        script_file
    )
)


ANALYSIS_DIR <- dirname(
    SCRIPT_DIR
)


DOWNSTREAM_DIR <- dirname(
    ANALYSIS_DIR
)


PROJECT_DIR <- dirname(
    DOWNSTREAM_DIR
)


CLASS_DIR <- file.path(
    PROJECT_DIR,
    "cre_classification"
)


RESULTS_DIR <- file.path(
    ANALYSIS_DIR,
    "results"
)


INPUT_FILE <- file.path(
    RESULTS_DIR,
    "climate_turnover_summary.tsv"
)


TREE_FILE <- file.path(
    CLASS_DIR,
    "phylogeny",
    "results",
    "301Fly_HOG_UCLDtree_40species.nw"
)


OUT_MODEL <- file.path(
    RESULTS_DIR,
    "pgls_climate_model_comparison.tsv"
)


OUT_COEFFICIENTS <- file.path(
    RESULTS_DIR,
    "pgls_climate_coefficients.tsv"
)


OUT_SUMMARY <- file.path(
    RESULTS_DIR,
    "pgls_climate_summary.txt"
)


# ============================================================
# Input checks
# ============================================================

if (!file.exists(INPUT_FILE)) {

    stop(
        paste(
            "Missing input file:",
            INPUT_FILE
        )
    )
}


if (!file.exists(TREE_FILE)) {

    stop(
        paste(
            "Missing phylogenetic tree:",
            TREE_FILE
        )
    )
}


# ============================================================
# Load species data
# ============================================================

dat <- read.delim(
    INPUT_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


required_columns <- c(
    "tree_name",
    "climatic_zone",
    "turnover_candidate",
    "mapped"
)


missing_columns <- setdiff(
    required_columns,
    colnames(dat)
)


if (length(missing_columns) > 0) {

    stop(
        paste(
            "Missing columns:",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    )
}


# ============================================================
# Response variable
#
# Empirical logit:
#
#   (turnover + 0.5) / (mapped + 1)
#
# This avoids ±Inf if a species has a turnover proportion
# exactly equal to 0 or 1.
# ============================================================

dat$turnover_rate_adjusted <- (
    dat$turnover_candidate + 0.5
) / (
    dat$mapped + 1
)


dat$logit_turnover_rate <- qlogis(
    dat$turnover_rate_adjusted
)


dat$climatic_zone <- factor(
    dat$climatic_zone,
    levels = c(
        "TROP",
        "ARID",
        "TEMP",
        "BORE"
    )
)


# ============================================================
# Load and prune tree
# ============================================================

tree <- read.tree(
    TREE_FILE
)


data_species <- dat$tree_name


missing_from_tree <- setdiff(
    data_species,
    tree$tip.label
)


if (length(missing_from_tree) > 0) {

    stop(
        paste(
            "Species missing from tree:",
            paste(
                missing_from_tree,
                collapse = ", "
            )
        )
    )
}


tree <- drop.tip(
    tree,
    setdiff(
        tree$tip.label,
        data_species
    )
)


# ============================================================
# Match data order to tree
# ============================================================

dat <- dat[
    match(
        tree$tip.label,
        dat$tree_name
    ),
]


if (!all(
    dat$tree_name == tree$tip.label
)) {

    stop(
        "ERROR: tree and data order do not match."
    )
}


dat$species_id <- dat$tree_name


# ============================================================
# Correlation structures
# ============================================================

cor_null <- corPagel(
    value = 0.5,
    phy = tree,
    fixed = FALSE,
    form = ~ species_id
)


cor_climate <- corPagel(
    value = 0.5,
    phy = tree,
    fixed = FALSE,
    form = ~ species_id
)


# ============================================================
# Models
#
# ML is used because models with different fixed effects
# are compared.
# ============================================================

model_null <- gls(
    logit_turnover_rate ~ 1,
    data = dat,
    correlation = cor_null,
    method = "ML"
)


model_climate <- gls(
    logit_turnover_rate ~ climatic_zone,
    data = dat,
    correlation = cor_climate,
    method = "ML"
)


# ============================================================
# Model comparison
# ============================================================

comparison <- anova(
    model_null,
    model_climate
)


comparison_df <- data.frame(
    model = rownames(
        comparison
    ),
    comparison,
    row.names = NULL,
    check.names = FALSE
)


write.table(
    comparison_df,
    OUT_MODEL,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# Coefficients
# ============================================================

coef_table <- summary(
    model_climate
)$tTable


coef_df <- data.frame(
    term = rownames(
        coef_table
    ),
    coef_table,
    row.names = NULL,
    check.names = FALSE
)


write.table(
    coef_df,
    OUT_COEFFICIENTS,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# Pagel lambda
# ============================================================

lambda_estimate <- coef(
    model_climate$modelStruct$corStruct,
    unconstrained = FALSE
)


# ============================================================
# Human-readable report
# ============================================================

sink(
    OUT_SUMMARY
)


cat(
    "Phylogenetically controlled CRE turnover analysis\n"
)

cat(
    "=================================================\n\n"
)


cat(
    "Response:\n"
)

cat(
    "logit((turnover_candidate + 0.5) / (mapped + 1))\n\n"
)


cat(
    "Predictor:\n"
)

cat(
    "climatic_zone\n\n"
)


cat(
    "Number of species: ",
    nrow(dat),
    "\n",
    sep = ""
)


cat(
    "Estimated Pagel lambda: ",
    lambda_estimate,
    "\n\n",
    sep = ""
)


cat(
    "MODEL COMPARISON\n"
)

cat(
    "----------------\n"
)

print(
    comparison
)


cat(
    "\n\nFULL CLIMATE MODEL\n"
)

cat(
    "------------------\n"
)

print(
    summary(
        model_climate
    )
)


sink()


# ============================================================
# Console summary
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "PGLS CLIMATE ANALYSIS COMPLETE\n"
)

cat(
    "============================================================\n\n"
)


cat(
    "Species: ",
    nrow(dat),
    "\n",
    sep = ""
)


cat(
    "Estimated Pagel lambda: ",
    lambda_estimate,
    "\n\n",
    sep = ""
)


print(
    comparison
)


cat(
    "\nWrote:\n",
    OUT_MODEL,
    "\n",
    OUT_COEFFICIENTS,
    "\n",
    OUT_SUMMARY,
    "\n",
    sep = ""
)
