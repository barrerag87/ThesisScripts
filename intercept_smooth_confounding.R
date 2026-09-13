# ============================================================
# SOURCE 2: INTERCEPT--SMOOTH CONFOUNDING
# ============================================================
#
# Intercept--smooth confounding occurs when a smooth term can
# generate a nonzero constant vector contribution.
#
# A constant can then be transferred between the intercept and
# the smooth term without changing the vector predictor.
#
# Componentwise centering removes this constant-shift
# ambiguity.
# ============================================================


run_intercept_smooth_confounding <- function(common) {

  n <- common$n
  x1 <- common$x1
  x2 <- common$x2
  alpha <- common$alpha
  B_centered <- common$B_centered
  X_intercept <- common$X_intercept


  # ----------------------------------------------------------
  # Bivariate smooth
  # ----------------------------------------------------------

  s0 <- (
    0.50 * sin(pi * x1) +
      0.30 * cos(pi * x2) +
      0.15 * x1 * x2
  )


  # The same smooth contributes to both predictor components.
  smooth_common <- cbind(
    s0,
    s0
  )


  eta_intercept <- add_intercept(
    smooth_common,
    alpha
  )


  # ----------------------------------------------------------
  # Constant shift
  # ----------------------------------------------------------
  #
  # Replace
  #
  #       s0(x)       by s0(x) - c0
  #
  # and compensate with
  #
  #       alpha       by alpha + (c0, c0)'.
  #
  # The predictor remains exactly the same.

  c0 <- 0.60

  s0_alt <- (
    s0 -
      c0
  )

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


  intercept_comparison <- rbind(
    original = alpha,
    alternative = alpha_alt
  )


  smooth_difference_intercept <- max_abs_diff(
    s0,
    s0_alt
  )


  # ----------------------------------------------------------
  # Centering
  # ----------------------------------------------------------
  #
  # The empirical restriction mean{s0(x_i)} = 0 removes the
  # constant direction from the smooth term.

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


  mean_s0_centered <- mean(
    s0_centered
  )


  # A further constant displacement violates the centering
  # condition.
  mean_shifted_centered_smooth <- mean(
    s0_centered - c0
  )


  # ----------------------------------------------------------
  # Rank deficiency before centering
  # ----------------------------------------------------------
  #
  # An explicit constant column is added to the smooth basis.
  # This creates overlap with the intercept space.

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


  # ----------------------------------------------------------
  # Rank after centering
  # ----------------------------------------------------------

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


  # Numerical verification.
  stopifnot(
    difference_intercept < 1e-10
  )

  stopifnot(
    difference_centering < 1e-10
  )


  list(
    s0 = s0,
    eta_intercept = eta_intercept,
    c0 = c0,
    s0_alt = s0_alt,
    alpha_alt = alpha_alt,
    eta_intercept_alt = eta_intercept_alt,
    difference_intercept = difference_intercept,
    intercept_comparison = intercept_comparison,
    smooth_difference_intercept = smooth_difference_intercept,
    mean_s0 = mean_s0,
    s0_centered = s0_centered,
    alpha_centered = alpha_centered,
    eta_centered = eta_centered,
    difference_centering = difference_centering,
    mean_s0_centered = mean_s0_centered,
    mean_shifted_centered_smooth = mean_shifted_centered_smooth,
    X_confounded = X_confounded,
    X_centered = X_centered,
    diagnostic_confounded = diagnostic_confounded,
    diagnostic_centered = diagnostic_centered
  )
}


# Optional plots for the intercept--smooth example.
plot_intercept_smooth_confounding <- function(
    common,
    result,
    grid_size = 60
) {

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

  s0_grid <- with(
    grid,
    0.50 * sin(pi * x1) +
      0.30 * cos(pi * x2) +
      0.15 * x1 * x2
  )

  S0_matrix <- matrix(
    s0_grid,
    nrow = grid_size,
    ncol = grid_size
  )

  S0_alt_matrix <- matrix(
    s0_grid - result$c0,
    nrow = grid_size,
    ncol = grid_size
  )

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

  invisible(result)
}
