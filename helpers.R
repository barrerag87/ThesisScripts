# ============================================================
# AUXILIARY FUNCTIONS
# ============================================================
#
# These functions are shared by all simulation modules.
# The objective is to keep each identifiability mechanism in a
# separate file while avoiding repeated code.
# ============================================================


# Maximum absolute difference between two numerical objects.
#
# A value close to zero indicates numerical equality up to
# floating-point precision.
max_abs_diff <- function(A, B) {
  max(abs(A - B))
}


# Numerical rank of a matrix.
#
# The tolerance is fixed so that rank calculations are
# reproducible across the simulation modules.
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


# Add a vector-valued intercept to an n x M matrix.
#
# S     : n x M matrix of smooth contributions.
# alpha : vector of length M.
add_intercept <- function(S, alpha) {
  sweep(
    S,
    MARGIN = 2,
    STATS = alpha,
    FUN = "+"
  )
}


# Construct a tensor-product B-spline basis for the bivariate
# covariate vector x = (x1, x2)'.
#
# If B1 has d1 columns and B2 has d2 columns, the resulting
# tensor-product basis contains d1*d2 columns.
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


# Stack two contribution matrices using component-major ordering.
#
# If A1 and A2 are n x d matrices, the result is 2n x d:
#
#   first n rows  -> contribution to eta_1
#   next n rows   -> contribution to eta_2
stack_components <- function(A1, A2) {
  rbind(A1, A2)
}


# ------------------------------------------------------------
# Common simulation setup
# ------------------------------------------------------------
#
# This function generates all objects shared by the simulation
# modules:
#
#   - bivariate covariates x1 and x2,
#   - the three smooth coordinates g, u1 and u2,
#   - the constraint matrix H,
#   - the vector-valued intercept,
#   - the true vector linear predictor,
#   - a Gaussian response,
#   - the centered tensor-product spline basis,
#   - the finite-dimensional intercept design.
#
# The random seed is set in main.R before this function is called.
build_common_setup <- function(n = 500) {

  # Two independent covariates.
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


  # Bivariate smooth coordinates.
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


  # Constraint matrix used by the kernel example.
  H_kernel <- matrix(
    c(
      1, 1, 0,
      1, 0, 1
    ),
    nrow = 2,
    byrow = TRUE
  )


  # Latent smooth vector evaluated at the observed points.
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


  # True vector linear predictor.
  eta_true <- add_intercept(
    F_kernel %*% t(H_kernel),
    alpha
  )


  # Gaussian response:
  #
  # eta_1 = mu
  # eta_2 = log(sigma)
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


  # Tensor-product B-spline basis used for the finite-sample
  # rank and nullity diagnostics.
  B_tensor <- tensor_bs_basis(
    x1 = x1,
    x2 = x2,
    df1 = 4,
    df2 = 4,
    degree = 3
  )

  B_centered <- center_columns(
    B_tensor
  )

  d <- ncol(
    B_centered
  )


  # Component-major intercept design:
  #
  # (eta_1(x_1), ..., eta_1(x_n),
  #  eta_2(x_1), ..., eta_2(x_n))'.
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


  # Coefficient vector used to construct explicit null
  # directions in the finite-dimensional examples.
  gamma <- seq(
    from = -0.3,
    to = 0.3,
    length.out = d
  )


  list(
    n = n,
    x1 = x1,
    x2 = x2,
    g = g,
    u1 = u1,
    u2 = u2,
    H_kernel = H_kernel,
    F_kernel = F_kernel,
    alpha = alpha,
    eta_true = eta_true,
    mu = mu,
    sigma = sigma,
    y = y,
    datos = datos,
    B_tensor = B_tensor,
    B_centered = B_centered,
    basis_means = round(colMeans(B_centered), 12),
    d = d,
    X_intercept = X_intercept,
    gamma = gamma
  )
}


# ------------------------------------------------------------
# Reference VGAM fit
# ------------------------------------------------------------
#
# This fit is not used to prove non-identifiability.
#
# Its purpose is to connect the algebraic constructions with an
# actual VGAM fit in which
#
#       eta_1 = mu
#       eta_2 = log(sigma).
fit_reference_vgam <- function(common) {

  vgam(
    y ~ s(x1, df = c(5, 5)) +
      s(x2, df = c(5, 5)),
    family = uninormal(zero = NULL),
    data = common$datos,
    trace = FALSE
  )
}
