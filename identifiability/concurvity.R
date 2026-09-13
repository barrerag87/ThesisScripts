# ============================================================
# SOURCE 3: CONCURVITY
# ============================================================
#
# Concurvity occurs when different smooth terms have
# overlapping contribution spaces.
#
# In this example, a nonzero bivariate function can be
# transferred from one smooth term to another without changing
# their total contribution to the vector predictor.
# ============================================================


run_concurvity <- function(common) {

  x1 <- common$x1
  x2 <- common$x2
  alpha <- common$alpha
  B_centered <- common$B_centered
  X_intercept <- common$X_intercept
  gamma <- common$gamma


  # ----------------------------------------------------------
  # Two bivariate smooth terms
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # Transfer a smooth function between terms
  # ----------------------------------------------------------
  #
  # Define
  #
  #       s1*(x) = s1(x) + q(x)
  #       s2*(x) = s2(x) - q(x).
  #
  # The individual terms change, while their sum remains fixed.

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


  smooth_difference_s1 <- max_abs_diff(
    s1,
    s1_alt
  )

  smooth_difference_s2 <- max_abs_diff(
    s2,
    s2_alt
  )


  # ----------------------------------------------------------
  # Finite-sample rank diagnostic
  # ----------------------------------------------------------
  #
  # Both terms use exactly the same contribution space, so the
  # two design blocks coincide.

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


  # ----------------------------------------------------------
  # Explicit coefficient-space null direction
  # ----------------------------------------------------------
  #
  # An increment in the first term is exactly compensated by
  # the opposite increment in the second term.

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


  # Numerical verification.
  stopifnot(
    difference_concurvity < 1e-10
  )

  stopifnot(
    concurvity_null_direction < 1e-10
  )


  list(
    s1 = s1,
    s2 = s2,
    eta_concurvity = eta_concurvity,
    q_concurvity = q_concurvity,
    s1_alt = s1_alt,
    s2_alt = s2_alt,
    eta_concurvity_alt = eta_concurvity_alt,
    difference_concurvity = difference_concurvity,
    smooth_difference_s1 = smooth_difference_s1,
    smooth_difference_s2 = smooth_difference_s2,
    X_concurvity = X_concurvity,
    diagnostic = diagnostic_concurvity,
    v_concurvity = v_concurvity,
    concurvity_null_direction = concurvity_null_direction
  )
}


# Optional plots for the concurvity example.
plot_concurvity <- function(common, result, grid_size = 60) {

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

  s1_grid <- with(
    grid,
    0.35 * sin(pi * x1) +
      0.20 * cos(pi * x2) +
      0.08 * x1 * x2
  )

  q_grid <- with(
    grid,
    0.16 * sin(2 * pi * x1) *
      sin(pi * x2) +
      0.04 * x1 * x2
  )

  S1_matrix <- matrix(
    s1_grid,
    nrow = grid_size,
    ncol = grid_size
  )

  S1_alt_matrix <- matrix(
    s1_grid + q_grid,
    nrow = grid_size,
    ncol = grid_size
  )

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

  invisible(result)
}
