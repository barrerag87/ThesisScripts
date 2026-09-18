backfitting <- function(y, x_list, family = gaussian(), tol = 1e-6,
                        max_iter = 1000, span = 0.5) {
  
  n <- length(y)
  p <- length(x_list)
  
  link_fun   <- family$linkfun
  inv_link   <- family$linkinv
  var_fun    <- family$variance
  mu_eta_fun <- family$mu.eta
  
  if (family$family == "binomial") {
    mu <- (y + 0.5) / 2
    mu <- pmin(pmax(mu, 0.01), 0.99)
  } else if (family$family == "poisson") {
    mu <- y + 0.1
  } else if (family$family == "Gamma") {
    mu <- pmax(y, 0.1)
  } else if (family$family == "gaussian" && family$link == "log") {
    mu <- pmax(y, 1e-5)
  } else if (family$family == "gaussian" && family$link == "inverse") {
    mu <- pmax(y, 0.1)
  } else {
    mu <- y
  }
  
  eta <- link_fun(mu)
  
  intercept <- mean(eta)
  
  smooth_functions <- lapply(
    x_list,
    function(x) rep(0, n)
  )
  
  deviance_old <- Inf
  deviance_new <- Inf
  iteration <- 0
  
  repeat {
    
    iteration <- iteration + 1
    
    eta <- intercept + Reduce(`+`, smooth_functions)
    mu <- inv_link(eta)
    
    if (family$family == "binomial") {
      mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)
    }
    
    if (family$family %in% c("poisson", "Gamma")) {
      mu <- pmax(mu, 1e-8)
    }
    
    if (any(!is.finite(mu))) {
      break
    }
    
    mu_eta <- mu_eta_fun(eta)
    
    if (any(!is.finite(mu_eta)) || any(abs(mu_eta) < 1e-12)) {
      break
    }
    
    weights <- mu_eta^2 / var_fun(mu)
    
    working_response <- eta + (y - mu) / mu_eta
    
    intercept <- weighted.mean(
      working_response - Reduce(`+`, smooth_functions),
      weights
    )
    
    for (j in seq_len(p)) {
      
      other_smooths <- if (p == 1) {
        rep(0, n)
      } else {
        Reduce(`+`, smooth_functions[-j])
      }
      
      partial_residual <- working_response -
        intercept -
        other_smooths
      
      dat_j <- data.frame(
        x = x_list[[j]],
        r = partial_residual
      )
      
      loess_fit <- tryCatch(
        loess(
          r ~ x,
          data = dat_j,
          weights = weights,
          span = span
        ),
        error = function(e) NULL
      )
      
      if (!is.null(loess_fit)) {
        
        est <- predict(
          loess_fit,
          newdata = data.frame(x = x_list[[j]])
        )
        
        if (all(is.finite(est))) {
          smooth_functions[[j]] <- est - mean(est)
        } else {
          smooth_functions[[j]] <- rep(0, n)
        }
        
      } else {
        
        smooth_functions[[j]] <- rep(0, n)
        
      }
    }
    
    eta <- intercept + Reduce(`+`, smooth_functions)
    mu <- inv_link(eta)
    
    if (family$family == "binomial") {
      mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)
    }
    
    if (family$family %in% c("poisson", "Gamma")) {
      mu <- pmax(mu, 1e-8)
    }
    
    if (any(!is.finite(mu))) {
      break
    }
    
    deviance_new <- sum(
      family$dev.resids(
        y,
        mu,
        rep(1, n)
      ),
      na.rm = TRUE
    )
    
    if (
      abs(deviance_old - deviance_new) < tol ||
      iteration >= max_iter
    ) {
      break
    }
    
    deviance_old <- deviance_new
  }
  
  result <- list(
    intercept = intercept,
    smooth_functions = smooth_functions,
    iterations = iteration,
    deviance = deviance_new,
    family = family,
    x_list = x_list,
    converged = iteration < max_iter &&
      is.finite(deviance_new) &&
      abs(deviance_old - deviance_new) < tol
  )
  
  class(result) <- "backfitting"
  
  result
}