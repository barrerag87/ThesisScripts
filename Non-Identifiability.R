# ============================================================
# IDENTIFIABILITY IN VGAMs
# MULTIVARIATE EXAMPLES
#
# Author: Rodrigo Barrera
#
# This script studies the main sources of non-identifiability
# at the level of a vector linear predictor.
#
# All examples are multivariate in the covariates:
#
#          x = (x1, x2)' in R^2.
#
# The vector linear predictor has dimension M = 2:
#
#          eta(x) = (eta_1(x), eta_2(x))'.
#
# The script illustrates:
#
#   1. Kernel ambiguity.
#   2. Intercept--smooth confounding.
#   3. Concurvity between smooth terms.
#   4. Conditional identifiability for H(theta).
#
# Every mechanism is verified in two ways:
#
#   a) Directly, by constructing two different internal
#      representations that generate the same vector predictor.
#
#   b) Through the rank and nullity of a finite-dimensional
#      design matrix based on a centered tensor-product
#      B-spline basis.
#
# The examples distinguish:
#
#   - The observed covariate vector x = (x1, x2).
#   - The latent smooth coordinates f(x).
#   - The vector linear predictor eta(x).
#   - The response generated from eta(x).
#
# The numerical examples are finite-sample analogues of the
# functional identifiability statements.
# ============================================================


# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

library(VGAM)
library(splines)


# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------

set.seed(1987)


# ============================================================
# AUXILIARY FUNCTIONS
# ============================================================


# Maximum absolute difference between two numerical objects.
#
# A value close to zero indicates numerical equality.

max_abs_diff <- function(A, B) {
  max(abs(A - B))
}


# Numerical rank of a matrix.

matrix_rank <- function(X, tol = 1e-10) {
  qr(X, tol = tol)$rank
}


# Center every column of a matrix.
#
# After the transformation, each column has empirical mean zero.

center_columns <- function(X) {
  sweep(
    X,
    MARGIN = 2,
    STATS = colMeans(X),
    FUN = "-"
  )
}


# Add a vector-valued intercept.
#
# S has dimension n x M.
# alpha has dimension M.

add_intercept <- function(S, alpha) {
  sweep(
    S,
    MARGIN = 2,
    STATS = alpha,
    FUN = "+"
  )
}


# Construct a tensor-product B-spline basis for the bivariate
# covariate vector x = (x1, x2).
#
# If B1 has d1 columns and B2 has d2 columns, the resulting
# tensor basis has d1*d2 columns.
#
# The (j,k)-th tensor basis function is
#
#          B1_j(x1) B2_k(x2).

tensor_bs_basis <- function(
    x1,
    x2,
    df1 = 4,
    df2 = 4,
    degree = 3
) {
  
  B1 <- bs(
    x1,
    df = df1,
    degree = degree,
    intercept = FALSE
  )
  
  B2 <- bs(
    x2,
    df = df2,
    degree = degree,
    intercept = FALSE
  )
  
  blocks <- vector(
    mode = "list",
    length = ncol(B1)
  )
  
  for (j in seq_len(ncol(B1))) {
    
    blocks[[j]] <- sweep(
      B2,
      MARGIN = 1,
      STATS = B1[, j],
      FUN = "*"
    )
  }
  
  B <- do.call(
    what = cbind,
    args = blocks
  )
  
  colnames(B) <- paste0(
    "B",
    seq_len(ncol(B))
  )
  
  B
}


# Stack a contribution matrix using component-major ordering.
#
# If A1 and A2 are n x d matrices, the result is 2n x d:
#
#   first n rows  -> contribution to eta_1
#   next n rows   -> contribution to eta_2

stack_components <- function(A1, A2) {
  rbind(A1, A2)
}


# ============================================================
# MULTIVARIATE DATA GENERATION
# ============================================================


# Number of observations.

n <- 500


# ------------------------------------------------------------
# Two covariates
# ------------------------------------------------------------
#
# The examples use a genuinely bivariate covariate vector:
#
#             x_i = (x1_i, x2_i)'.
#
# x1 and x2 are generated independently so that the examples do
# not rely on deterministic collinearity between the covariates.

x1 <- runif(
  n = n,
  min = -1,
  max = 1
)

x2 <- runif(
  n = n,
  min = -1,
  max = 1
)


# ------------------------------------------------------------
# Bivariate smooth coordinates
# ------------------------------------------------------------
#
# Every smooth coordinate depends on both x1 and x2.
#
# This makes the examples multivariate at the covariate level.

g <- (
  0.40 * sin(pi * x1) +
    0.25 * cos(pi * x2) +
    0.15 * x1 * x2
)

u1 <- (
  0.18 * cos(2 * pi * x1) -
    0.12 * sin(pi * x2) +
    0.08 * x1 * x2
)

u2 <- (
  -0.16 * sin(2 * pi * x2) +
    0.10 * cos(pi * x1) -
    0.06 * x1 * x2
)


# ------------------------------------------------------------
# Constraint matrix
# ------------------------------------------------------------
#
# The structure
#
#            [1  1  0]
#       H =  [1  0  1]
#
# generates
#
# eta_1(x1,x2) = alpha_1 + g(x1,x2) + u1(x1,x2)
#
# eta_2(x1,x2) = alpha_2 + g(x1,x2) + u2(x1,x2).
#
# H has rank two and three columns.
#
# Therefore ker(H) is nontrivial.

H_kernel <- matrix(
  c(
    1, 1, 0,
    1, 0, 1
  ),
  nrow = 2,
  byrow = TRUE
)

H_kernel


# Matrix containing the evaluations of the latent smooth vector
#
#       f(x1,x2) = (g(x1,x2), u1(x1,x2), u2(x1,x2))'.

F_kernel <- cbind(
  g = g,
  u1 = u1,
  u2 = u2
)


# Vector-valued intercept.

alpha <- c(
  0.8,
  -0.5
)


# ------------------------------------------------------------
# True vector linear predictor
# ------------------------------------------------------------
#
# F_kernel has dimension n x 3.
# t(H_kernel) has dimension 3 x 2.
#
# Hence
#
#       F_kernel %*% t(H_kernel)
#
# has dimension n x 2.

eta_true <- add_intercept(
  F_kernel %*% t(H_kernel),
  alpha
)


# ============================================================
# RESPONSE GENERATION AND VGAM FIT
# ============================================================
#
# The response remains univariate Gaussian.
#
# The model has two linear predictors:
#
#       eta_1 = mu
#       eta_2 = log(sigma).
#
# The multivariate structure in this script refers to:
#
#   - the bivariate covariate vector (x1, x2),
#   - the two-dimensional vector linear predictor eta,
#   - the multivariate latent smooth representation.
#
# This is enough for the identifiability mechanisms studied here.

mu <- eta_true[, 1]

sigma <- exp(
  eta_true[, 2]
)


y <- rnorm(
  n = n,
  mean = mu,
  sd = sigma
)


datos <- data.frame(
  y = y,
  x1 = x1,
  x2 = x2
)


# Actual VGAM fit using both covariates.

fit_vgam <- vgam(
  y ~ s(x1, df = c(5, 5)) +
    s(x2, df = c(5, 5)),
  family = uninormal(zero = NULL),
  data = datos,
  trace = FALSE
)


# Estimated vector linear predictor.

eta_hat <- predict(
  fit_vgam,
  type = "link"
)


head(eta_hat)

summary(fit_vgam)

constraints(fit_vgam)


# ============================================================
# FINITE BIVARIATE SMOOTH BASIS
# ============================================================
#
# The theoretical mechanisms are functional.
#
# For numerical rank diagnostics, every bivariate smooth is
# represented by a tensor-product B-spline basis.
#
# This basis contains products
#
#        B_j(x1) C_k(x2),
#
# so it can represent smooth functions of both covariates.

B_tensor <- tensor_bs_basis(
  x1 = x1,
  x2 = x2,
  df1 = 4,
  df2 = 4,
  degree = 3
)


# Remove empirical constant components.

B_centered <- center_columns(
  B_tensor
)


# Number of bivariate basis functions.

d <- ncol(
  B_centered
)

d


# Verify empirical centering.

basis_means <- round(
  colMeans(B_centered),
  12
)

basis_means


# ------------------------------------------------------------
# Intercept design
# ------------------------------------------------------------
#
# Component-major stacking:
#
# (
#   eta_1(x_1),
#   ...,
#   eta_1(x_n),
#   eta_2(x_1),
#   ...,
#   eta_2(x_n)
# )'.
#
# x_i denotes the bivariate vector (x1_i, x2_i)'.

X_intercept <- rbind(
  cbind(
    rep(1, n),
    rep(0, n)
  ),
  cbind(
    rep(0, n),
    rep(1, n)
  )
)


# ============================================================
# SOURCE 1
# KERNEL AMBIGUITY
# ============================================================
#
# Kernel ambiguity occurs when a nonzero admissible smooth
# vector h(x1,x2) satisfies
#
#                  H h(x1,x2) = 0
#
# for every point (x1,x2).
#
# The latent smooth representation changes.
#
# The vector predictor remains unchanged.


# ------------------------------------------------------------
# Matrix-level kernel
# ------------------------------------------------------------

rank_H_kernel <- matrix_rank(
  H_kernel
)

rank_H_kernel


# The vector
#
#             v = (1, -1, -1)'
#
# belongs to ker(H).

kernel_direction_check <- H_kernel %*% c(
  1,
  -1,
  -1
)

kernel_direction_check


# ------------------------------------------------------------
# Nonzero bivariate functional perturbation
# ------------------------------------------------------------
#
# q_kernel depends on both covariates.

q_kernel <- (
  0.20 * sin(2 * pi * x1) *
    cos(pi * x2) +
    0.08 * x1 * x2
)


# Original representation:
#
# f(x1,x2)
# =
# (g(x1,x2), u1(x1,x2), u2(x1,x2))'.
#
# Alternative representation:
#
# g*(x1,x2)  = g(x1,x2)  + q(x1,x2)
#
# u1*(x1,x2) = u1(x1,x2) - q(x1,x2)
#
# u2*(x1,x2) = u2(x1,x2) - q(x1,x2)
#
# The difference is
#
#     q(x1,x2) (1,-1,-1)',
#
# which lies pointwise in ker(H).

F_kernel_alt <- cbind(
  g = g + q_kernel,
  u1 = u1 - q_kernel,
  u2 = u2 - q_kernel
)


eta_kernel_alt <- add_intercept(
  F_kernel_alt %*% t(H_kernel),
  alpha
)


# The vector predictors must coincide.

difference_kernel <- max_abs_diff(
  eta_true,
  eta_kernel_alt
)

difference_kernel


# The smooth representation must change.

smooth_difference_kernel <- max_abs_diff(
  F_kernel,
  F_kernel_alt
)

smooth_difference_kernel


# ------------------------------------------------------------
# Finite-dimensional design matrix
# ------------------------------------------------------------
#
# Each latent smooth coordinate uses the same centered
# bivariate tensor-product basis.


# g contributes to both components.

Z_g <- stack_components(
  B_centered,
  B_centered
)


# u1 contributes only to eta_1.

Z_u1 <- stack_components(
  B_centered,
  matrix(
    0,
    nrow = n,
    ncol = d
  )
)


# u2 contributes only to eta_2.

Z_u2 <- stack_components(
  matrix(
    0,
    nrow = n,
    ncol = d
  ),
  B_centered
)


X_kernel <- cbind(
  X_intercept,
  Z_g,
  Z_u1,
  Z_u2
)


columns_kernel <- ncol(
  X_kernel
)

rank_kernel <- matrix_rank(
  X_kernel
)

nullity_kernel <- (
  columns_kernel -
    rank_kernel
)


diagnostic_kernel <- data.frame(
  mechanism = "Kernel ambiguity",
  columns = columns_kernel,
  rank = rank_kernel,
  nullity = nullity_kernel
)

diagnostic_kernel


# ------------------------------------------------------------
# Explicit coefficient-space null direction
# ------------------------------------------------------------
#
# Each bivariate basis coefficient can move in the same
# kernel direction.

gamma <- seq(
  from = -0.3,
  to = 0.3,
  length.out = d
)


# Parameter order:
#
# intercepts | g | u1 | u2

v_kernel <- c(
  0,
  0,
  gamma,
  -gamma,
  -gamma
)


kernel_null_direction <- max(
  abs(
    X_kernel %*% v_kernel
  )
)

kernel_null_direction


# ============================================================
# SOURCE 2
# INTERCEPT--SMOOTH CONFOUNDING
# ============================================================
#
# Intercept--smooth confounding occurs when a smooth term can
# generate a nonzero constant vector contribution.
#
# A constant can then be transferred between the vector-valued
# intercept and the smooth term without changing eta(x1,x2).


# ------------------------------------------------------------
# Bivariate smooth
# ------------------------------------------------------------

s0 <- (
  0.50 * sin(pi * x1) +
    0.30 * cos(pi * x2) +
    0.15 * x1 * x2
)


# The same bivariate smooth contributes to both predictor
# components.

smooth_common <- cbind(
  s0,
  s0
)


eta_intercept <- add_intercept(
  smooth_common,
  alpha
)


# ------------------------------------------------------------
# Constant shift
# ------------------------------------------------------------

c0 <- 0.60


# Shift the smooth by a constant.

s0_alt <- s0 - c0


# Compensate through the intercept.

alpha_alt <- alpha + c(
  c0,
  c0
)


eta_intercept_alt <- add_intercept(
  cbind(
    s0_alt,
    s0_alt
  ),
  alpha_alt
)


difference_intercept <- max_abs_diff(
  eta_intercept,
  eta_intercept_alt
)

difference_intercept


intercept_comparison <- rbind(
  original = alpha,
  alternative = alpha_alt
)

intercept_comparison


smooth_difference_intercept <- max_abs_diff(
  s0,
  s0_alt
)

smooth_difference_intercept


# ------------------------------------------------------------
# Centering
# ------------------------------------------------------------
#
# The empirical restriction
#
#             mean{s0(x1_i,x2_i)} = 0
#
# removes the constant direction from the smooth term.

mean_s0 <- mean(
  s0
)


s0_centered <- (
  s0 -
    mean_s0
)


alpha_centered <- alpha + c(
  mean_s0,
  mean_s0
)


eta_centered <- add_intercept(
  cbind(
    s0_centered,
    s0_centered
  ),
  alpha_centered
)


difference_centering <- max_abs_diff(
  eta_intercept,
  eta_centered
)

difference_centering


mean_s0_centered <- mean(
  s0_centered
)

mean_s0_centered


# A new constant displacement violates the centering condition.

mean_shifted_centered_smooth <- mean(
  s0_centered - c0
)

mean_shifted_centered_smooth


# ------------------------------------------------------------
# Rank deficiency before centering
# ------------------------------------------------------------
#
# Add an explicit constant direction to the bivariate smooth
# basis.

B_uncentered <- cbind(
  constant = rep(1, n),
  B_centered
)


Z_uncentered <- stack_components(
  B_uncentered,
  B_uncentered
)


X_confounded <- cbind(
  X_intercept,
  Z_uncentered
)


columns_confounded <- ncol(
  X_confounded
)

rank_confounded <- matrix_rank(
  X_confounded
)

nullity_confounded <- (
  columns_confounded -
    rank_confounded
)


diagnostic_confounded <- data.frame(
  mechanism = "Intercept--smooth confounding",
  columns = columns_confounded,
  rank = rank_confounded,
  nullity = nullity_confounded
)

diagnostic_confounded


# ------------------------------------------------------------
# Rank after centering
# ------------------------------------------------------------

Z_centered <- stack_components(
  B_centered,
  B_centered
)


X_centered <- cbind(
  X_intercept,
  Z_centered
)


columns_centered <- ncol(
  X_centered
)

rank_centered <- matrix_rank(
  X_centered
)

nullity_centered <- (
  columns_centered -
    rank_centered
)


diagnostic_centered <- data.frame(
  mechanism = "After centering",
  columns = columns_centered,
  rank = rank_centered,
  nullity = nullity_centered
)

diagnostic_centered


# ============================================================
# SOURCE 3
# CONCURVITY
# ============================================================
#
# Concurvity occurs when different smooth terms have
# overlapping contribution spaces.
#
# The example uses two distinct bivariate smooth terms.
#
# Both terms are functions of (x1,x2).
#
# Their common contribution structure permits transfer of a
# nonzero bivariate function from one term to the other.


s1 <- (
  0.35 * sin(pi * x1) +
    0.20 * cos(pi * x2) +
    0.08 * x1 * x2
)


s2 <- (
  0.25 * cos(2 * pi * x1) -
    0.18 * sin(pi * x2) +
    0.05 * x1 * x2
)


eta_concurvity <- add_intercept(
  cbind(
    s1 + s2,
    s1 + s2
  ),
  alpha
)


# ------------------------------------------------------------
# Transfer a bivariate smooth function between terms
# ------------------------------------------------------------

q_concurvity <- (
  0.16 * sin(2 * pi * x1) *
    sin(pi * x2) +
    0.04 * x1 * x2
)


s1_alt <- (
  s1 +
    q_concurvity
)

s2_alt <- (
  s2 -
    q_concurvity
)


eta_concurvity_alt <- add_intercept(
  cbind(
    s1_alt + s2_alt,
    s1_alt + s2_alt
  ),
  alpha
)


difference_concurvity <- max_abs_diff(
  eta_concurvity,
  eta_concurvity_alt
)

difference_concurvity


smooth_difference_s1 <- max_abs_diff(
  s1,
  s1_alt
)

smooth_difference_s2 <- max_abs_diff(
  s2,
  s2_alt
)

smooth_difference_s1

smooth_difference_s2


# ------------------------------------------------------------
# Finite-sample rank diagnostic
# ------------------------------------------------------------
#
# Both smooth terms use the same bivariate contribution space.
#
# Therefore their design matrices coincide exactly.

Z_1 <- stack_components(
  B_centered,
  B_centered
)

Z_2 <- stack_components(
  B_centered,
  B_centered
)


X_concurvity <- cbind(
  X_intercept,
  Z_1,
  Z_2
)


columns_concurvity <- ncol(
  X_concurvity
)

rank_concurvity <- matrix_rank(
  X_concurvity
)

nullity_concurvity <- (
  columns_concurvity -
    rank_concurvity
)


diagnostic_concurvity <- data.frame(
  mechanism = "Concurvity",
  columns = columns_concurvity,
  rank = rank_concurvity,
  nullity = nullity_concurvity
)

diagnostic_concurvity


# ------------------------------------------------------------
# Explicit coefficient-space null direction
# ------------------------------------------------------------

v_concurvity <- c(
  0,
  0,
  gamma,
  -gamma
)


concurvity_null_direction <- max(
  abs(
    X_concurvity %*% v_concurvity
  )
)

concurvity_null_direction


# ============================================================
# SUMMARY OF THE THREE MECHANISMS
# ============================================================

results <- data.frame(
  mechanism = c(
    "Kernel ambiguity",
    "Intercept--smooth confounding",
    "Concurvity"
  ),
  max_predictor_difference = c(
    difference_kernel,
    difference_intercept,
    difference_concurvity
  ),
  columns = c(
    columns_kernel,
    columns_confounded,
    columns_concurvity
  ),
  rank = c(
    rank_kernel,
    rank_confounded,
    rank_concurvity
  ),
  nullity = c(
    nullity_kernel,
    nullity_confounded,
    nullity_concurvity
  )
)

results


# ============================================================
# NUMERICAL VERIFICATION
# ============================================================
#
# The script stops if any transformation fails to preserve the
# vector predictor up to numerical precision.

stopifnot(
  difference_kernel < 1e-10
)

stopifnot(
  difference_intercept < 1e-10
)

stopifnot(
  difference_concurvity < 1e-10
)

stopifnot(
  difference_centering < 1e-10
)

stopifnot(
  kernel_null_direction < 1e-10
)

stopifnot(
  concurvity_null_direction < 1e-10
)


# ============================================================
# CONDITIONAL EXTENSION
# PARAMETRIZED CONSTRAINT MATRIX H(theta)
# ============================================================
#
# This section studies conditional identifiability when the
# smooth vector f* is held fixed and theta changes.
#
# The covariates remain bivariate:
#
#             x = (x1,x2)'.
#
# The latent coordinate vector has dimension K = 2.


H_theta <- function(theta) {
  
  matrix(
    c(
      1, theta,
      0, 1
    ),
    nrow = 2,
    byrow = TRUE
  )
}


theta1 <- 0.5

theta2 <- 1.5


# ============================================================
# CASE A
# f* GENERATES ONLY A PROPER SUBSPACE OF R^2
# ============================================================
#
# z1 depends on both x1 and x2.
#
# The second latent coordinate is identically zero.
#
# Therefore the observed latent vectors explore only the
# direction generated by (1,0)'.

z1 <- (
  sin(pi * x1) +
    0.30 * cos(pi * x2) +
    0.10 * x1 * x2
)


z1 <- (
  z1 -
    mean(z1)
)


Fstar_restricted <- cbind(
  z1,
  rep(0, n)
)


rank_Fstar_restricted <- matrix_rank(
  Fstar_restricted
)

rank_Fstar_restricted


eta_theta1_restricted <- (
  Fstar_restricted %*%
    t(
      H_theta(theta1)
    )
)


eta_theta2_restricted <- (
  Fstar_restricted %*%
    t(
      H_theta(theta2)
    )
)


difference_theta_restricted <- max_abs_diff(
  eta_theta1_restricted,
  eta_theta2_restricted
)

difference_theta_restricted


# ============================================================
# CASE B
# f* GENERATES R^2
# ============================================================
#
# Two linearly independent latent coordinates are used.
#
# Both depend on the full bivariate covariate vector.

z1_full <- (
  sin(pi * x1) +
    0.25 * cos(pi * x2) +
    0.08 * x1 * x2
)

z2_full <- (
  cos(2 * pi * x1) -
    0.20 * sin(pi * x2) +
    0.12 * x1 * x2
)


Fstar_full <- cbind(
  z1_full,
  z2_full
)


Fstar_full <- center_columns(
  Fstar_full
)


rank_Fstar_full <- matrix_rank(
  Fstar_full
)

rank_Fstar_full


eta_theta1_full <- (
  Fstar_full %*%
    t(
      H_theta(theta1)
    )
)


eta_theta2_full <- (
  Fstar_full %*%
    t(
      H_theta(theta2)
    )
)


difference_theta_full <- max_abs_diff(
  eta_theta1_full,
  eta_theta2_full
)

difference_theta_full


conditional_results <- data.frame(
  case = c(
    "Proper subspace",
    "Full span"
  ),
  rank_Fstar = c(
    rank_Fstar_restricted,
    rank_Fstar_full
  ),
  predictor_difference = c(
    difference_theta_restricted,
    difference_theta_full
  )
)

conditional_results


# ============================================================
# ADDITIONAL CHECK
# EMPIRICAL DIMENSION OF THE COVARIATE VECTOR
# ============================================================
#
# This check confirms that the examples do not collapse to a
# single observed covariate direction.

X_covariates <- cbind(
  x1,
  x2
)

rank_X_covariates <- matrix_rank(
  X_covariates
)

rank_X_covariates


stopifnot(
  rank_X_covariates == 2
)

stopifnot(
  rank_Fstar_restricted == 1
)

stopifnot(
  rank_Fstar_full == 2
)

stopifnot(
  difference_theta_restricted < 1e-10
)

stopifnot(
  difference_theta_full > 1e-8
)


# ============================================================
# PLOTS
# ============================================================
#
# Since every smooth depends on two covariates, a one-dimensional
# line plot is no longer sufficient.
#
# The following plots use the observed pairs (x1,x2) and encode
# the smooth value through symbols generated by image().
#
# For a compact visualization, the bivariate domain is evaluated
# on a regular grid.


# ------------------------------------------------------------
# Regular grid
# ------------------------------------------------------------

grid_size <- 60

x1_grid <- seq(
  from = -1,
  to = 1,
  length.out = grid_size
)

x2_grid <- seq(
  from = -1,
  to = 1,
  length.out = grid_size
)


grid <- expand.grid(
  x1 = x1_grid,
  x2 = x2_grid
)


# Original g surface.

g_grid <- with(
  grid,
  0.40 * sin(pi * x1) +
    0.25 * cos(pi * x2) +
    0.15 * x1 * x2
)


# Kernel perturbation surface.

q_kernel_grid <- with(
  grid,
  0.20 * sin(2 * pi * x1) *
    cos(pi * x2) +
    0.08 * x1 * x2
)


# Original common smooth.

s0_grid <- with(
  grid,
  0.50 * sin(pi * x1) +
    0.30 * cos(pi * x2) +
    0.15 * x1 * x2
)


# First concurvity term.

s1_grid <- with(
  grid,
  0.35 * sin(pi * x1) +
    0.20 * cos(pi * x2) +
    0.08 * x1 * x2
)


# Arrange values as matrices for image().

G_matrix <- matrix(
  g_grid,
  nrow = grid_size,
  ncol = grid_size
)

G_alt_matrix <- matrix(
  g_grid + q_kernel_grid,
  nrow = grid_size,
  ncol = grid_size
)

S0_matrix <- matrix(
  s0_grid,
  nrow = grid_size,
  ncol = grid_size
)

S0_alt_matrix <- matrix(
  s0_grid - c0,
  nrow = grid_size,
  ncol = grid_size
)

S1_matrix <- matrix(
  s1_grid,
  nrow = grid_size,
  ncol = grid_size
)


q_concurvity_grid <- with(
  grid,
  0.16 * sin(2 * pi * x1) *
    sin(pi * x2) +
    0.04 * x1 * x2
)


S1_alt_matrix <- matrix(
  s1_grid + q_concurvity_grid,
  nrow = grid_size,
  ncol = grid_size
)


# ------------------------------------------------------------
# Plot 1: original bivariate smooth g
# ------------------------------------------------------------

image(
  x1_grid,
  x2_grid,
  G_matrix,
  xlab = "x1",
  ylab = "x2",
  main = "Kernel ambiguity: g(x1,x2)"
)


# ------------------------------------------------------------
# Plot 2: alternative bivariate smooth g*
# ------------------------------------------------------------

image(
  x1_grid,
  x2_grid,
  G_alt_matrix,
  xlab = "x1",
  ylab = "x2",
  main = "Kernel ambiguity: g*(x1,x2)"
)


# ------------------------------------------------------------
# Plot 3: intercept confounding
# ------------------------------------------------------------

image(
  x1_grid,
  x2_grid,
  S0_matrix,
  xlab = "x1",
  ylab = "x2",
  main = "Intercept confounding: s0(x1,x2)"
)

image(
  x1_grid,
  x2_grid,
  S0_alt_matrix,
  xlab = "x1",
  ylab = "x2",
  main = "Intercept confounding: s0*(x1,x2)"
)


# ------------------------------------------------------------
# Plot 4: concurvity
# ------------------------------------------------------------

image(
  x1_grid,
  x2_grid,
  S1_matrix,
  xlab = "x1",
  ylab = "x2",
  main = "Concurvity: s1(x1,x2)"
)

image(
  x1_grid,
  x2_grid,
  S1_alt_matrix,
  xlab = "x1",
  ylab = "x2",
  main = "Concurvity: s1*(x1,x2)"
)


# ------------------------------------------------------------
# Predictor invariance
# ------------------------------------------------------------

differences <- c(
  Kernel = difference_kernel,
  Intercept = difference_intercept,
  Concurvity = difference_concurvity
)


barplot(
  differences,
  ylab = "Maximum absolute difference",
  main = "Predictor invariance"
)


# ============================================================
# FINAL DIAGNOSTIC OUTPUT
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "MULTIVARIATE VGAM IDENTIFIABILITY EXAMPLES\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Covariate dimension p =",
  rank_X_covariates,
  "\n"
)

cat(
  "Vector predictor dimension M =",
  ncol(eta_true),
  "\n"
)

cat(
  "Bivariate tensor-basis dimension d =",
  d,
  "\n\n"
)

print(
  results
)

cat(
  "\nConditional H(theta) examples:\n"
)

print(
  conditional_results
)

cat(
  "\nAll invariance checks passed.\n"
)

cat(
  "============================================================\n"
)
