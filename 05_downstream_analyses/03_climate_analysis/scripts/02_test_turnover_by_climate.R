#!/usr/bin/env Rscript

# ============================================================
# Phylogenetically controlled climate analysis
#
# Purpose:
#   Test whether species-level CRE turnover differs among
#   climatic zones while accounting for phylogenetic
#   non-independence.
#
# Response:
#   Empirical-logit transformed evaluable turnover rate:
#
#       logit((turnover_candidate + 0.5) / (mapped + 1))
#
# Model:
#   logit turnover rate ~ climatic_zone
#
# Phylogenetic covariance:
#   Pagel's lambda estimated by GLS.
#
# Null and climate models are fitted by maximum likelihood
# because models with different fixed effects are compared.
#
# Input/output paths are supplied by the climate-analysis
# configuration via the pipeline wrapper.
# ============================================================


# ============================================================
# Packages
# ============================================================

required_packages <- c("ape", "nlme")

missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
    stop(
        "Missing R packages: ",
        paste(missing_packages, collapse = ", ")
    )
}


# ============================================================
# Arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name) {
    position <- match(name, args)

    if (is.na(position) || position == length(args)) {
        stop("Missing required argument: ", name)
    }

    args[position + 1]
}

INPUT_FILE <- get_arg("--input")
TREE_FILE <- get_arg("--tree")
OUT_PREDICTIONS <- get_arg("--out-predictions")
OUT_STATS <- get_arg("--out-stats")
OUT_MODEL <- get_arg("--out-model")
OUT_COEFFICIENTS <- get_arg("--out-coefficients")
OUT_SUMMARY <- get_arg("--out-summary")


# ============================================================
# Input checks
# ============================================================

for (path in c(INPUT_FILE, TREE_FILE)) {
    if (!file.exists(path)) {
        stop("Required file not found: ", path)
    }
}

for (path in c(OUT_MODEL, OUT_COEFFICIENTS, OUT_SUMMARY)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
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

missing_columns <- setdiff(required_columns, names(dat))

if (length(missing_columns) > 0) {
    stop(
        "Missing columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

if (nrow(dat) == 0) {
    stop("Climate turnover table is empty.")
}

if (anyDuplicated(dat$tree_name)) {
    stop("Duplicate tree_name values in climate turnover table.")
}


# ============================================================
# Prepare response variable
# ============================================================

dat$turnover_candidate <- as.numeric(dat$turnover_candidate)
dat$mapped <- as.numeric(dat$mapped)

if (anyNA(dat$turnover_candidate) || anyNA(dat$mapped)) {
    stop("Non-numeric turnover_candidate or mapped values.")
}

if (any(dat$mapped <= 0)) {
    stop("All species must have at least one mapped CRE.")
}

if (
    any(dat$turnover_candidate < 0) ||
    any(dat$turnover_candidate > dat$mapped)
) {
    stop("Invalid turnover_candidate counts.")
}

# Empirical-logit correction avoids infinite values when the
# observed turnover proportion is exactly zero or one.
dat$turnover_rate_adjusted <- (
    dat$turnover_candidate + 0.5
) / (
    dat$mapped + 1
)

dat$logit_turnover_rate <- qlogis(dat$turnover_rate_adjusted)


# ============================================================
# Climatic-zone predictor
# ============================================================

CLIMATE_LEVELS <- c("TROP", "ARID", "TEMP", "BORE")

unknown_climates <- setdiff(
    unique(dat$climatic_zone),
    CLIMATE_LEVELS
)

if (length(unknown_climates) > 0) {
    stop(
        "Unknown climatic-zone values: ",
        paste(unknown_climates, collapse = ", ")
    )
}

# TROP is the reference level of the climate model.
dat$climatic_zone <- factor(
    dat$climatic_zone,
    levels = CLIMATE_LEVELS
)


# ============================================================
# Load and prune phylogeny
# ============================================================

tree <- ape::read.tree(TREE_FILE)

missing_from_tree <- setdiff(
    dat$tree_name,
    tree$tip.label
)

if (length(missing_from_tree) > 0) {
    stop(
        "Species missing from phylogenetic tree: ",
        paste(missing_from_tree, collapse = ", ")
    )
}

tree <- ape::drop.tip(
    tree,
    setdiff(tree$tip.label, dat$tree_name)
)


# ============================================================
# Match data order to tree
# ============================================================

dat <- dat[
    match(tree$tip.label, dat$tree_name),
]

if (!all(dat$tree_name == tree$tip.label)) {
    stop("Tree and species-data order do not match.")
}

dat$species_id <- dat$tree_name


# ============================================================
# Pagel-lambda correlation structures
# ============================================================

# Separate correlation objects are used because lambda is
# estimated independently for each GLS model.
cor_null <- ape::corPagel(
    value = 0.5,
    phy = tree,
    fixed = FALSE,
    form = ~ species_id
)

cor_climate <- ape::corPagel(
    value = 0.5,
    phy = tree,
    fixed = FALSE,
    form = ~ species_id
)


# ============================================================
# PGLS models
# ============================================================

# ML is used because the two models differ in fixed effects.
model_null <- nlme::gls(
    logit_turnover_rate ~ 1,
    data = dat,
    correlation = cor_null,
    method = "ML",
    na.action = na.fail
)

model_climate <- nlme::gls(
    logit_turnover_rate ~ climatic_zone,
    data = dat,
    correlation = cor_climate,
    method = "ML",
    na.action = na.fail
)


# ============================================================
# Model comparison
# ============================================================

comparison <- anova(
    model_null,
    model_climate
)

comparison_df <- data.frame(
    model = c("null", "climate"),
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
# Climate-model coefficients
# ============================================================

coef_table <- summary(model_climate)$tTable

coef_df <- data.frame(
    term = rownames(coef_table),
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

lambda_estimate <- as.numeric(
    coef(
        model_climate$modelStruct$corStruct,
        unconstrained = FALSE
    )
)

# ============================================================
# Model-based climate estimates
# ============================================================

newdat <- data.frame(
    climatic_zone = factor(CLIMATE_LEVELS, levels = CLIMATE_LEVELS)
)

X <- model.matrix(~ climatic_zone, data = newdat)
beta <- coef(model_climate)
V <- vcov(model_climate)

eta <- as.numeric(X %*% beta)
se_eta <- sqrt(diag(X %*% V %*% t(X)))

predictions <- data.frame(
    climatic_zone = CLIMATE_LEVELS,
    estimate = plogis(eta),
    ci95_low = plogis(eta - 1.96 * se_eta),
    ci95_high = plogis(eta + 1.96 * se_eta)
)

write.table(
    predictions,
    OUT_PREDICTIONS,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

p_global <- comparison$"p-value"[2]

stats <- data.frame(
    n_species = nrow(dat),
    pagel_lambda = lambda_estimate,
    p_global = p_global
)

write.table(
    stats,
    OUT_STATS,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# Human-readable report
# ============================================================

sink(OUT_SUMMARY)

cat("Phylogenetically controlled CRE turnover analysis\n")
cat("=================================================\n\n")

cat("Response:\n")
cat("logit((turnover_candidate + 0.5) / (mapped + 1))\n\n")

cat("Predictor:\n")
cat("climatic_zone\n\n")

cat("Reference climatic zone:\n")
cat("TROP\n\n")

cat("Number of species: ", nrow(dat), "\n", sep = "")
cat("Estimated Pagel lambda: ", lambda_estimate, "\n\n", sep = "")

cat("MODEL COMPARISON\n")
cat("----------------\n")
print(comparison)

cat("\n\nFULL CLIMATE MODEL\n")
cat("------------------\n")
print(summary(model_climate))

sink()


# ============================================================
# Console summary
# ============================================================

cat("\n============================================================\n")
cat("PGLS CLIMATE ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat("Species: ", nrow(dat), "\n", sep = "")
cat("Estimated Pagel lambda: ", lambda_estimate, "\n\n", sep = "")

print(comparison)

cat(
    "\nWrote:\n",
    OUT_MODEL, "\n",
    OUT_COEFFICIENTS, "\n",
    OUT_SUMMARY, "\n",
    sep = ""
)
