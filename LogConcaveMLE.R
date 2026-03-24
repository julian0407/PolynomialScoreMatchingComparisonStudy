# ============================================================
# LogConcaveMLE.R
# Classical log-concave MLE baselines for the comparison study
#
# Univariate:
#   - logcondens::logConDens
#
# Multivariate:
#   - LogConcDEAD::mlelcd
#
# API:
#   - fit_logconcave_mle_1d()
#   - fit_logconcave_mle_mv()
#   - predict_density_logconcave_1d()
#   - predict_logdensity_logconcave_1d()
#   - predict_score_logconcave_1d()
#   - predict_density_logconcave_mv()
#   - predict_logdensity_logconcave_mv()
#   - predict_score_logconcave_mv()
# ============================================================


# ------------------------------------------------------------
# (1) 1D Log-concave MLE
# ------------------------------------------------------------

# ------------------------------------------------------------
# (1.1) 1D Log-concave MLE
# ------------------------------------------------------------
fit_logconcave_mle_1d <- function(x,
                                  smoothed = FALSE,
                                  print = FALSE) {
  # Check if x is finite and sample contains at least two distinct samples.
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(unique(x)) < 2L) stop("Sample must contain at least two distinct samples.")
  
  # Apply active set approach for log-concave MLE
  fit <- logcondens::logConDens(
    x = sort(x),
    smoothed = smoothed,
    print = print
  )

  # TODO: Class? 
  # Return list of information relevant for evalution metrics
  structure(
    list(
      fit_MLE = fit,
      smoothed = smoothed
    ),
    class = "logconcave_1d_fit"
  )
}

# ------------------------------------------------------------
# (1.2) Calculate Estimates + Evaluation Metrics in 1D
#       Log-concave MLE
# ------------------------------------------------------------
predict_density_logconcave_1d <- function(newx,
                                          fit) {
  # Get the fitting characteristic and check if newx is finite
  smoothed <- fit$smoothed
  newx <- as.numeric(newx)
  newx <- newx[is.finite(newx)]
  
  fit_MLE <- fit$fit_MLE
  
  if (!isTRUE(smoothed)) {
    # exact unsmoothed MLE:
    # log f is piecewise linear on [x_(1), x_(n)] and f = 0 outside.
    xmin <- min(fit_MLE$x)
    xmax <- max(fit_MLE$x)
    
    # Initilaize density vector with zeros and identify points in newx
    # that can be estimated by the fit_MLE
    out <- numeric(length(newx))
    inside <- (newx >= xmin) & (newx <= xmax)
    
    # apply piecewise linear interpolation to get density estimate 
    # of newx inside the support of fit
    if (any(inside)) {
      phi_eval <- approx(
        x = fit_MLE$x,
        y = fit_MLE$phi,
        xout = newx[inside],
        method = "linear",
        rule = 2
      )$y
      out[inside] <- exp(phi_eval)
    }
    return(out)
  }
  
  # if smoothed=True apply smoothed estimate from package
  # logcondens
  ans <- logcondens::evaluateLogConDens(
    x = newx,
    object = fit_MLE,
    which = "smoothed"
  )
  return(as.numeric(ans$y))
}

predict_logdensity_logconcave_1d <- function(newx,
                                             fit,
                                             eps = 1e-300) {
  log(pmax(predict_density_logconcave_1d(newx, fit), eps))
}

predict_score_logconcave_1d <- function(newx,
                                        fit,
                                        h = 1e-4,
                                        eps = 1e-300) {
  # Get the fitting characteristic and check if newx is finite
  smoothed <- fit$smoothed
  newx <- as.numeric(newx)
  newx <- newx[is.finite(newx)]
  
  # if not smoothed than caluclate exact piecewise-constant
  # derivative of log-density:
  if (!isTRUE(smoothed)) {
    # Get fitting data and calculate slopes from piecewise linear interpolation
    fit_MLE <- fit$fit_MLE
    knots <- as.numeric(fit_MLE$x)
    phi   <- as.numeric(fit_MLE$phi)
    slopes <- diff(phi) / diff(knots)
    # Initilaize output vector with zeros and identify points in
    # newx that are inside the fitting range
    out <- rep(0, length(newx))
    inside <- newx >= min(knots) & newx <= max(knots)
    
    if (any(inside)) {
      idx <- findInterval(newx[inside], knots, rightmost.closed = TRUE)
      # error handling for pints at min/max range
      idx[idx < 1L] <- 1L
      idx[idx >= length(knots)] <- length(knots) - 1L
      # TODO: why negative?
      # score = - d/dx log f
      out[inside] <- -slopes[idx]
    }
    
    # outside support, density is zero and score is not well-defined for log f
    out[!inside] <- NA_real_
    return(as.numeric(out))
  }
  
  # smoothed case: use numerical derivative of log density
  f_log <- function(xx) {
    predict_logdensity_logconcave_1d(xx, fit, eps = eps)
  }
  
  # Apply numerical derivative on every sample point in newx
  out <- vapply(
    newx,
    FUN = function(xx) {
      if (!is.finite(xx)) return(NA_real_)
      # TODO: why negative?
      # score = - d/dx log f
      -num_derivative_1d(f_log, xx, h = h)
    },
    FUN.VALUE = numeric(1)
  )
  as.numeric(out)
}

# ------------------------------------------------------------
# (2) Multivariate log-concave MLE
# ------------------------------------------------------------

# ------------------------------------------------------------
# (2.1) Fitting Multivariate log-concave MLE
# ------------------------------------------------------------

fit_logconcave_mle_mv <- function(x,
                                  smoothed = FALSE) {
  
  # Identify rows that contains no infinite entry
  x <- as_obs_matrix(x)
  keep <- apply(x, 1, function(row) all(is.finite(row)))
  # remain matrix property (drop = false) in case only one row remains
  x <- x[keep, , drop = FALSE]
  # sample size n is at least two
  if (nrow(x) < 2L) stop("x must have at least two rows.")

  d <- ncol(x)
  
  # Need full-dimensional convex hull -> compute rank of
  # centralized matrix (by column means)
  if (qr(x - matrix(colMeans(x), nrow(x), d, byrow = TRUE))$rank < d) {
    stop("Data are not full-dimensional; multivariate log-concave MLE is not well-defined.")
  }
  
  # Use LogConcDead package to calculate multivariate log-concave MLE
  # with tent functions
  fit <- LogConcDEAD::mlelcd(x)
  
  # TODO: Class? 
  # Return list of information relevant for evalution metrics
  structure(
    list(
      fit = fit,
      x = x,
      smoothed = smoothed,
      d = d
    ),
    class = "logconcave_mv_fit"
  )
}

# ------------------------------------------------------------------
# (2.2) Calculate Estimates + Evaluation Metrics in
#       Multivariate log-concave MLE
# ------------------------------------------------------------------

predict_density_logconcave_mv <- function(newx,
                                          fit) {
  # Identify rows that contains no infinite entry
  newx <- as_obs_matrix(newx)
  keep <- apply(newx, 1, function(row) all(is.finite(row)))
  # remain matrix property (drop = false) in case only one row remains
  newx <- newx[keep, , drop = FALSE]
  
  # Check if dimension of data to estimate matches with dimension of fitting data
  if (ncol(newx) != fit$d) stop("Dimension mismatch.")
  
  # Get the fitting characteristic
  smoothed <- fit$smoothed

  if (!isTRUE(smoothed)) {
    # Use package function to derive not smoothed estimate for newx
    return(as.numeric(LogConcDEAD::dlcd(newx, fit$fit, uselog = FALSE)))
  }
  
  # if smoothed=True, than estimate smoothness matrix A_hat
  Ahat <- LogConcDEAD::hatA(fit$fit)
  # Use package function to derive smoothed estimate for newx
  as.numeric(LogConcDEAD::dslcd(newx, fit$fit, A = Ahat))
}

predict_logdensity_logconcave_mv <- function(newx,
                                             fit,
                                             eps = 1e-300) {
  log(pmax(predict_density_logconcave_mv(newx, fit), eps))
}

predict_score_logconcave_mv <- function(newx,
                                        fit,
                                        h = 1e-4,
                                        eps = 1e-300) {

  newx <- as_obs_matrix(newx)
  keep <- apply(newx, 1, function(row) all(is.finite(row)))
  # remain matrix property (drop = false) in case only one row remains
  newx <- newx[keep, , drop = FALSE]
  
  # Check if dimension of data to estimate matches with dimension of fitting data
  if (ncol(newx) != fit$d) stop("Dimension mismatch.")
  
  # Calculate numerical derivative of log density
  f_log <- function(xx) {
    # ensure row is formated as matrix
    xx <- matrix(xx, nrow = 1)
    predict_logdensity_logconcave_mv(xx, fit, eps = eps)
  }
  
  out <- t(apply(newx, 1, function(row) {
    # TODO: why negative?
    # score = - d/dx log f
    -num_gradient_mv(f_log, row, h = h)
  }))
  out
}