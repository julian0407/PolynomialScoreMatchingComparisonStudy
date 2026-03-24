# ============================================================
# Basic KDE estimators for the comparison study
#
# Main practical multivariate baseline:
#   TODO: Check relevance of different h
#   Gaussian KDE with ks::Hpi / Hlscv / Hns
#
# API:
#   - fit_kde_1d()
#   - fit_kde_mv()
#   - predict_density_kde_1d()
#   - predict_logdensity_kde_1d()
#   - predict_score_kde_1d()
#   - predict_density_kde_mv()
#   - predict_logdensity_kde_mv()
#   - predict_score_kde_mv()
# ============================================================

# ------------------------------------------------------------
# (1) 1D KDE
# ------------------------------------------------------------

# ------------------------------------------------------------
# (1.1) Fitting 1D KDE
# ------------------------------------------------------------
fit_kde_1d <- function(x,
                       bw = c("SJ", "nrd0", "ucv", "bcv")) {
  
  # Check if type of bandwith selection is allowed               
  bw <- match.arg(bw)
  
  # Check if x is finite and sample size n is at least two.
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) stop("Sample size must be at least two.")
  
  # Deterine bandwidth based on sample x and selected bandwidth type bw
  bw_value <- switch(
    bw,
    SJ   = stats::bw.SJ(x),
    nrd0 = stats::bw.nrd0(x),
    ucv  = stats::bw.ucv(x),
    bcv  = stats::bw.bcv(x)
  )
  
  # Check if bandwidth parameter calculation is valid
  if (!is.finite(bw_value) || bw_value <= 0) {
    stop("Bandwidth selection failed.")
  }
  
  # Return cleaned data and calculated bandwidth-parameter
  # TODO: Why structure?
  structure(
    list(
      x = x,
      n = length(x),
      bw_method = bw,
      h = bw_value
    ),
    class = "kde_1d_fit"
  )
}

# ------------------------------------------------------------
# (1.2) Calculate Estimates + Evaluation Metrics in 1D KDE
# ------------------------------------------------------------
predict_density_kde_1d <- function(newx, fit) {
  # Clean data to be estimated
  newx <- as.numeric(newx)
  newx <- newx[is.finite(newx)]
  
  # Get training data and parameters from fit
  x <- fit$x
  h <- fit$h
  
  # Create matrix that is used to calculate KDE-estimate for each point in newx
  # Each row "represents" one point-estimate
  z <- outer(newx, x, FUN = "-") / h
  dens <- rowMeans(stats::dnorm(z)) / h
  as.numeric(dens)
}

predict_logdensity_kde_1d <- function(newx, fit, eps = 1e-300) {
  log(pmax(as.numeric(predict_density_kde_1d(newx, fit)), eps))
}

predict_score_kde_1d <- function(newx, fit, eps = 1e-300) {
  # Clean data whose score should be calculated
  newx <- as.numeric(newx)
  newx <- newx[is.finite(newx)]
  
  # Get training data and parameters from fit
  x <- fit$x
  h <- fit$h
  
  # Get K to determine the derivative of Gaussian KDE
  z <- outer(newx, x, FUN = "-") / h
  K <- stats::dnorm(z)
  
  # derivative of Gaussian KDE:
  # f'(x) = mean( -(x - Xi)/h^3 * phi((x-Xi)/h) )
  deriv <- rowMeans((-outer(newx, x, FUN = "-") / h^3) * K)
  
  # Get estimated density from fit
  dens <- pmax(as.numeric(predict_density_kde_1d(newx, fit)), eps)
  
  # negative derivative of log density
  # TODO: Check why the negative is used here?
  score <- -deriv / dens
  as.numeric(score)
}

# ------------------------------------------------------------
# (2.1) Fitting Multivariate KDE
# ------------------------------------------------------------
fit_kde_mv <- function(x,
                       H_method = c("Hpi", "Hlscv", "Hns"),
                       diagonal = FALSE) {
  
  # Check if type of bandwith selection is allowed   
  H_method <- match.arg(H_method)
  
  # Identify rows that contains no infinite entry
  x <- as_obs_matrix(x)
  keep <- apply(x, 1, function(row) all(is.finite(row)))
  # remain matrix property (drop = false) in case only one row remains
  x <- x[keep, , drop = FALSE]
  # sample size n is at least two
  if (nrow(x) < 2L) stop("x must have at least two rows.")
  
  # Determine dimension of data
  d <- ncol(x)
  # Deterine bandwidth matrix H based on sample x and selected H_method
  H <- switch(
    H_method,
    Hpi   = ks::Hpi(x = x),
    Hlscv = ks::Hlscv(x = x),
    Hns   = ks::Hns(x = x)
  )
  # check if dimensions of matrix H matches with data such that
  # Multivariate KDE is well defined
  if (!is.matrix(H) || any(dim(H) != c(d, d))) {
    stop("Bandwidth matrix H has wrong dimension.")
  }
  
  # Optional: Reduce bandwidth matrix H to diagonal structure
  if (diagonal) {
    H <- diag(diag(H), nrow = d, ncol = d)
  }
  
  # TODO: Can be interpreted as Ridge Correction to ensure H to
  # be positive definite
  # Determine eigenvalues and adapt matrix H such that H is positive definite
  eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) <= 0) {
    H <- H + diag(abs(min(eig)) + 1e-8, d)
  }
  
  # TODO: check if structure is needed
  structure(
    list(
      x = x,
      n = nrow(x),
      d = d,
      H = H,
      H_inv = solve(H),
      detH = det(H),
      H_method = H_method,
      diagonal = diagonal
    ),
    class = "kde_mv_fit"
  )
}

# ------------------------------------------------------------------
# (2.2) Calculate Estimates + Evaluation Metrics in Multivariate KDE
# ------------------------------------------------------------------

predict_density_kde_mv <- function(newx, fit) {
  # Clean data whose score should be calculated
  newx <- as_obs_matrix(newx)
  keep <- apply(newx, 1, function(row) all(is.finite(row)))
  newx <- newx[keep, , drop = FALSE]
  # Check if dimensions match with data from fit
  if (ncol(newx) != fit$d) stop("Dimension mismatch.")
  
  # Get training data and parameters from fit
  x <- fit$x
  H_inv <- fit$H_inv
  detH <- fit$detH
  d <- fit$d
  n <- fit$n
  
  # TODO: Check mathematics
  # Calculate normalizing constant
  const <- 1 / (((2 * pi)^(d / 2)) * sqrt(detH) * n)
  # Initialize vector for denisties of each sample point in newx
  out <- numeric(nrow(newx))
  
  # Caluclate for each sample point in newx the KDE estimate
  for (j in seq_len(nrow(newx))) {
    # Matrix containing the diff in each dimension (row) to the trainings data (column)
    diffs <- sweep(x, 2, newx[j, ], FUN = "-")
    # Caluclation of estimate
    quad <- rowSums((diffs %*% H_inv) * diffs)
    out[j] <- const * sum(exp(-0.5 * quad))
  }
  out
}

predict_logdensity_kde_mv <- function(newx, fit, eps = 1e-300) {
  log(pmax(as.numeric(predict_density_kde_mv(newx, fit)), eps))
}

predict_score_kde_mv <- function(newx, fit, eps = 1e-300) {
  # Clean data whose score should be calculated
  newx <- as_obs_matrix(newx)
  keep <- apply(newx, 1, function(row) all(is.finite(row)))
  newx <- newx[keep, , drop = FALSE]
  # Check if dimensions match with data from fit
  if (ncol(newx) != fit$d) stop("Dimension mismatch.")
  
  # Get training data and parameters from fit
  x <- fit$x
  H_inv <- fit$H_inv
  detH <- fit$detH
  d <- fit$d
  n <- fit$n
  
  # TODO: Check mathematics
  # Calculate normalizing constant
  const <- 1 / (((2 * pi)^(d / 2)) * sqrt(detH) * n)
  # Initialize vector for denisties of each sample point in newx
  out <- matrix(NA_real_, nrow = nrow(newx), ncol = d)
  
  # Calculate density and corresponding gradient for each point in newx
  # to determine score of the KDE estimate
  for (j in seq_len(nrow(newx))) {
    diffs <- sweep(x, 2, newx[j, ], FUN = "-")
    quad <- rowSums((diffs %*% H_inv) * diffs)
    w <- exp(-0.5 * quad)
    
    # Get estimated densities
    dens_j <- max(const * sum(w), eps)
    
    # gradient of density
    grad_f <- const * colSums(w * (diffs %*% H_inv))
    
    # TODO: why negative
    # negative gradient of log density
    out[j, ] <- -grad_f / dens_j
  }
  out
}