#!/usr/bin/env Rscript

# ============================================================
# Species-level CRE turnover and climatic zone
#
# Panel A:
#   observed turnover_rate_evaluable by climate
#
# Panel B:
#   phylogenetically adjusted PGLS estimates with 95% CI
#
# Output:
#   results/figures/climate_PGLS.png
# ============================================================


# ============================================================
# Packages
# ============================================================

packages <- c(
    "ape",
    "nlme"
)

missing <- packages[
    !vapply(
        packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (length(missing) > 0) {
    stop(
        paste(
            "Missing packages:",
            paste(missing, collapse = ", ")
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

SCRIPT_FILE <- sub(
    "^--file=",
    "",
    file_arg
)

SCRIPT_DIR <- dirname(
    normalizePath(SCRIPT_FILE)
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

INPUT_FILE <- file.path(
    ANALYSIS_DIR,
    "results",
    "climate_turnover_summary.tsv"
)

TREE_FILE <- file.path(
    PROJECT_DIR,
    "cre_classification",
    "phylogeny",
    "results",
    "301Fly_HOG_UCLDtree_40species.nw"
)

FIG_DIR <- file.path(
    ANALYSIS_DIR,
    "results",
    "figures"
)

dir.create(
    FIG_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

OUT_PNG <- file.path(
    FIG_DIR,
    "climate_PGLS.png"
)


# ============================================================
# Settings
# ============================================================

SHOW_POINTS <- TRUE

climate_colors <- c(
    TROP = "#E41A1C",
    ARID = "#FBC02D",
    TEMP = "#009E73",
    BORE = "#0072B2"
)

climate_labels <- c(
    TROP = "Tropical",
    ARID = "Arid",
    TEMP = "Temperate",
    BORE = "Boreal"
)

climate_order <- c(
    "TROP",
    "ARID",
    "TEMP",
    "BORE"
)


# ============================================================
# Load data
# ============================================================

dat <- read.delim(
    INPUT_FILE,
    stringsAsFactors = FALSE
)

dat$climatic_zone <- factor(
    dat$climatic_zone,
    levels = climate_order
)


# ============================================================
# Fit PGLS model
# ============================================================

dat$turnover_rate_adjusted <- (
    dat$turnover_candidate + 0.5
) / (
    dat$mapped + 1
)

dat$logit_turnover_rate <- qlogis(
    dat$turnover_rate_adjusted
)

tree <- read.tree(
    TREE_FILE
)

tree <- drop.tip(
    tree,
    setdiff(
        tree$tip.label,
        dat$tree_name
    )
)

dat <- dat[
    match(
        tree$tip.label,
        dat$tree_name
    ),
]

dat$species_id <- dat$tree_name

cor_structure <- corPagel(
    value = 0.5,
    phy = tree,
    fixed = FALSE,
    form = ~ species_id
)

model <- gls(
    logit_turnover_rate ~ climatic_zone,
    data = dat,
    correlation = cor_structure,
    method = "ML"
)

null_cor <- corPagel(
    value = 0.5,
    phy = tree,
    fixed = FALSE,
    form = ~ species_id
)

null_model <- gls(
    logit_turnover_rate ~ 1,
    data = dat,
    correlation = null_cor,
    method = "ML"
)

comparison <- anova(
    null_model,
    model
)

p_global <- comparison$"p-value"[2]

lambda <- coef(
    model$modelStruct$corStruct,
    unconstrained = FALSE
)


# ============================================================
# Model-based adjusted estimates
# ============================================================

newdat <- data.frame(
    climatic_zone = factor(
        climate_order,
        levels = climate_order
    )
)

X <- model.matrix(
    ~ climatic_zone,
    data = newdat
)

beta <- coef(
    model
)

V <- vcov(
    model
)

eta <- as.numeric(
    X %*% beta
)

se_eta <- sqrt(
    diag(
        X %*% V %*% t(X)
    )
)

lower_eta <- eta - 1.96 * se_eta
upper_eta <- eta + 1.96 * se_eta

pred <- plogis(
    eta
)

lower <- plogis(
    lower_eta
)

upper <- plogis(
    upper_eta
)


# ============================================================
# Helpers
# ============================================================

make_fill <- function(cols, alpha = 0.28) {
    grDevices::adjustcolor(
        cols,
        alpha.f = alpha
    )
}

make_point_fill <- function(cols, alpha = 0.85) {
    grDevices::adjustcolor(
        cols,
        alpha.f = alpha
    )
}


# ============================================================
# Drawing function
# ============================================================

draw_figure <- function() {

    layout(
        matrix(c(1, 2), nrow = 1),
        widths = c(1.18, 1.0)
    )

    # --------------------------------------------------------
    # Panel A
    # --------------------------------------------------------

    par(
        mar = c(5.4, 5.6, 3.4, 1.5)
    )

    y <- dat$turnover_rate_evaluable * 100

    grouped_values <- split(
        y,
        dat$climatic_zone
    )

    bp <- boxplot(
        grouped_values,
        outline = FALSE,
        xaxt = "n",
        col = make_fill(
            climate_colors[climate_order],
            alpha = 0.28
        ),
        border = "#666666",
        ylab = "Turnover candidates among evaluable CRE loci (%)",
        xlab = "Climatic zone",
        las = 1,
        frame.plot = TRUE
    )

    axis(
        1,
        at = 1:4,
        labels = climate_labels[
            climate_order
        ]
    )

    if (SHOW_POINTS) {
        set.seed(1)

        for (i in seq_along(climate_order)) {

            zone <- climate_order[i]

            values <- dat$turnover_rate_evaluable[
                dat$climatic_zone == zone
            ] * 100

            x <- jitter(
                rep(i, length(values)),
                amount = 0.10
            )

            points(
                x,
                values,
                pch = 21,
                bg = make_point_fill(
                    climate_colors[zone],
                    alpha = 0.90
                ),
                col = "#222222",
                cex = 0.85,
                lwd = 0.7
            )
        }
    }

    title(
        main = "Observed species-level turnover",
        line = 1.0,
        font.main = 2
    )

    mtext(
        "A",
        side = 3,
        adj = -0.10,
        line = 2.0,
        font = 2,
        cex = 1.6
    )


    # --------------------------------------------------------
    # Panel B
    # --------------------------------------------------------

    par(
        mar = c(5.4, 5.0, 3.4, 1.5)
    )

    ylim <- range(
        c(lower, upper) * 100
    )

    ylim <- ylim + c(-2.5, 2.5)

    plot(
        1:4,
        pred * 100,
        type = "n",
        xaxt = "n",
        xlab = "Climatic zone",
        ylab = "Phylogenetically adjusted turnover rate (%)",
        ylim = ylim,
        las = 1
    )

    axis(
        1,
        at = 1:4,
        labels = climate_labels[
            climate_order
        ]
    )

    segments(
        x0 = 1:4,
        y0 = lower * 100,
        x1 = 1:4,
        y1 = upper * 100,
        lwd = 1.4,
        col = "#222222"
    )

    segments(
        x0 = (1:4) - 0.07,
        y0 = lower * 100,
        x1 = (1:4) + 0.07,
        y1 = lower * 100,
        lwd = 1.4,
        col = "#222222"
    )

    segments(
        x0 = (1:4) - 0.07,
        y0 = upper * 100,
        x1 = (1:4) + 0.07,
        y1 = upper * 100,
        lwd = 1.4,
        col = "#222222"
    )

    points(
        1:4,
        pred * 100,
        pch = 21,
        bg = climate_colors[
            climate_order
        ],
        col = "#222222",
        cex = 1.4,
        lwd = 0.8
    )

    title(
        main = "Phylogenetically adjusted estimates",
        line = 1.0,
        font.main = 2
    )

    mtext(
        "B",
        side = 3,
        adj = -0.02,
        line = 2.0,
        font = 2,
        cex = 1.6
    )

    legend(
        "topright",
        legend = c(
            sprintf("PGLS: p = %.3f", p_global),
            sprintf("Pagel's lambda = %.3f", lambda)
        ),
        bty = "n",
        cex = 0.92,
        inset = c(0.00, 0.02)
    )
}

# ============================================================
# Save PNG
# ============================================================

png(
    OUT_PNG,
    width = 3360,
    height = 1680,
    res = 300
)

draw_figure()
dev.off()


cat(
    "\nWrote:\n",
    OUT_PNG,
    "\n",
    sep = ""
)
