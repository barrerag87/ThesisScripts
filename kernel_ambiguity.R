# ============================================================
# SOURCE 1: KERNEL AMBIGUITY
# ============================================================
#
# Kernel ambiguity occurs when a nonzero admissible smooth
# perturbation h(x) satisfies
#
#                  H h(x) = 0
#
# for every observed point.
#
# The internal smooth representation changes, while the vector
# linear predictor remains unchanged.
# ============================================================


run_kernel_ambiguity <- function(common) {

  n <- common$n
  x1 <- common$x1
  x2 <- common$x2
  g <- common$g
  u1 <- common$u1
  u2 <- common$u2
  H_kernel <- common$H_kernel
  F_kernel <- common$F_kernel
  alpha <- common$alpha
  eta_true <- common$eta_true
  B_centered <- common$B_centered
  d <- common$d
  X_intercept <- common$X_intercept
  gamma <- common$gamma


  # ----------------------------------------------------------
  # Matrix-level kernel
  # ----------------------------------------------------------

  rank_H_kernel <- matrix_rank(
    H_kernel
  )


  # v = (1, -1, -1)' belongs to ker(H).
  v <- c(
    1,
    -1,
    -1
  )

  kernel_direction_check <- H_kernel %*% v


  # ----------------------------------------------------------
  # Nonzero functional perturbation
  # ----------------------------------------------------------
  #
  # q_kernel depends on both covariates. The alternative smooth
  # representation is
  #
  # f*(x) = f(x) + q_kernel(x) v.
  #
  # Since H v = 0, this perturbation cannot change H f(x).

  q_kernel <- (
    0.20 * sin(2 * pi * x1) *
      cos(pi * x2) +
      0.08 * x1 * x2
  )


  F_kernel_alt <- cbind(
    g = g + q_kernel,
    u1 = u1 - q_kernel,
    u2 = u2 - q_kernel
  )


  eta_kernel_alt <- add_intercept(
    F_kernel_alt %*% t(H_kernel),
    alpha
  )


  # Maximum difference between the predictor induced by f and
  # the predictor induced by f*.
  difference_kernel <- max_abs_diff(
    eta_true,
    eta_kernel_alt
  )


  # The internal smooth representation must actually change.
  smooth_difference_kernel <- max_abs_diff(
    F_kernel,
    F_kernel_alt
  )


  # ----------------------------------------------------------
  # Finite-dimensional design matrix
  # ----------------------------------------------------------
  #
  # The same centered bivariate tensor-product basis is used for
  # each latent smooth coordinate.

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


  # ----------------------------------------------------------
  # Explicit coefficient-space null direction
  # ----------------------------------------------------------
  #
  # Parameter order:
  #
  # intercepts | g | u1 | u2
  #
  # Each basis coefficient is perturbed according to the same
  # kernel direction (1, -1, -1)'.

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


  # Numerical verification.
  stopifnot(
    difference_kernel < 1e-10
  )

  stopifnot(
    kernel_null_direction < 1e-10
  )


  list(
    rank_H_kernel = rank_H_kernel,
    kernel_direction_check = kernel_direction_check,
    q_kernel = q_kernel,
    F_kernel_alt = F_kernel_alt,
    eta_kernel_alt = eta_kernel_alt,
    difference_kernel = difference_kernel,
    smooth_difference_kernel = smooth_difference_kernel,
    X_kernel = X_kernel,
    diagnostic = diagnostic_kernel,
    v_kernel = v_kernel,
    kernel_null_direction = kernel_null_direction
  )
}


# Optional plot for the kernel example.
#
# The function is not called automatically by main.R unless
# make_plots is set to TRUE.
plot_kernel_ambiguity <- function(common, result, grid_size = 60) {

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

  g_grid <- with(
    grid,
    0.40 * sin(pi * x1) +
      0.25 * cos(pi * x2) +
      0.15 * x1 * x2
  )

  q_grid <- with(
    grid,
    0.20 * sin(2 * pi * x1) *
      cos(pi * x2) +
      0.08 * x1 * x2
  )

  G_matrix <- matrix(
    g_grid,
    nrow = grid_size,
    ncol = grid_size
  )

  G_alt_matrix <- matrix(
    g_grid + q_grid,
    nrow = grid_size,
    ncol = grid_size
  )

  image(
    x1_grid,
    x2_grid,
    G_matrix,
    xlab = "x1",
    ylab = "x2",
    main = "Kernel ambiguity: g(x1,x2)"
  )

  image(
    x1_grid,
    x2_grid,
    G_alt_matrix,
    xlab = "x1",
    ylab = "x2",
    main = "Kernel ambiguity: g*(x1,x2)"
  )

  invisible(result)
}
