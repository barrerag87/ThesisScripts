# ============================================================
# IDENTIFIABILITY IN VGAMs
# MODULAR SIMULATION
#
# Author: Rodrigo Barrera
# ============================================================
#
# This is the main script.
#
# Each source of non-identifiability is implemented in a
# separate file:
#
#   1. R/kernel_ambiguity.R
#   2. R/intercept_smooth_confounding.R
#   3. R/concurvity.R
#   4. R/conditional_H_theta.R
#
# Shared functions and common simulated objects are defined in
# R/helpers.R.
#
# The algebraic constructions establish the identifiability
# mechanisms directly. The VGAM fit is included as a reference
# implementation connecting the construction with a fitted
# vector generalized additive model.
# ============================================================


# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

library(VGAM)
library(splines)

setwd("C:/Users/rodri/OneDrive/Escritorio/REV/SIM/R")
# ------------------------------------------------------------
# Source modular files
# ------------------------------------------------------------

source("R/helpers.R")
source("R/kernel_ambiguity.R")
source("R/intercept_smooth_confounding.R")
source("R/concurvity.R")
source("R/conditional_H_theta.R")


# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------

set.seed(1987)


# ------------------------------------------------------------
# Shared simulation objects
# ------------------------------------------------------------
#
# The random covariates, response, smooth coordinates, basis,
# intercept and true predictor are generated only once.
#
# This guarantees that every module works with exactly the same
# simulated sample.

common <- build_common_setup(
  n = 500
)


# ------------------------------------------------------------
# Reference VGAM fit
# ------------------------------------------------------------
#
# This is a single fitted VGAM.
#
# It is not used to construct the alternative representations
# in the identifiability examples. Those representations are
# generated algebraically inside the corresponding modules.

fit_vgam <- fit_reference_vgam(
  common
)


eta_hat <- predict(
  fit_vgam,
  type = "link"
)


fit_summary <- summary(
  fit_vgam
)


fit_constraints <- constraints(
  fit_vgam
)


# ------------------------------------------------------------
# Run each identifiability module
# ------------------------------------------------------------

kernel <- run_kernel_ambiguity(
  common
)


intercept_smooth <- run_intercept_smooth_confounding(
  common
)


concurvity <- run_concurvity(
  common
)


conditional_H <- run_conditional_H_theta(
  common
)


# ------------------------------------------------------------
# Summary of the three exact non-identifiability mechanisms
# ------------------------------------------------------------

results <- data.frame(
  mechanism = c(
    "Kernel ambiguity",
    "Intercept--smooth confounding",
    "Concurvity"
  ),
  max_predictor_difference = c(
    kernel$difference_kernel,
    intercept_smooth$difference_intercept,
    concurvity$difference_concurvity
  ),
  columns = c(
    kernel$diagnostic$columns,
    intercept_smooth$diagnostic_confounded$columns,
    concurvity$diagnostic$columns
  ),
  rank = c(
    kernel$diagnostic$rank,
    intercept_smooth$diagnostic_confounded$rank,
    concurvity$diagnostic$rank
  ),
  nullity = c(
    kernel$diagnostic$nullity,
    intercept_smooth$diagnostic_confounded$nullity,
    concurvity$diagnostic$nullity
  )
)


# Conditional H(theta) summary.

conditional_results <- conditional_H$conditional_results


# ------------------------------------------------------------
# Global numerical verification
# ------------------------------------------------------------

stopifnot(
  kernel$difference_kernel < 1e-10
)

stopifnot(
  intercept_smooth$difference_intercept < 1e-10
)

stopifnot(
  concurvity$difference_concurvity < 1e-10
)

stopifnot(
  intercept_smooth$difference_centering < 1e-10
)

stopifnot(
  kernel$kernel_null_direction < 1e-10
)

stopifnot(
  concurvity$concurvity_null_direction < 1e-10
)


# ------------------------------------------------------------
# Optional plots
# ------------------------------------------------------------
#
# Change this value to TRUE to reproduce the bivariate plots.
#
# The plots are kept inside the corresponding simulation files
# so each module remains self-contained.

make_plots <- FALSE


if (make_plots) {

  plot_kernel_ambiguity(
    common,
    kernel
  )

  plot_intercept_smooth_confounding(
    common,
    intercept_smooth
  )

  plot_concurvity(
    common,
    concurvity
  )


  differences <- c(
    Kernel = kernel$difference_kernel,
    Intercept = intercept_smooth$difference_intercept,
    Concurvity = concurvity$difference_concurvity
  )


  barplot(
    differences,
    ylab = "Maximum absolute difference",
    main = "Predictor invariance"
  )
}


# ------------------------------------------------------------
# Objects to inspect after running main.R
# ------------------------------------------------------------
#
# Main outputs:
#
#   results
#   conditional_results
#   eta_hat
#   fit_summary
#   fit_constraints
#
# Detailed outputs:
#
#   kernel
#   intercept_smooth
#   concurvity
#   conditional_H
#
# No print() or cat() calls are required. In RStudio, inspect
# these objects directly in the Environment or type their names
# in the console.
