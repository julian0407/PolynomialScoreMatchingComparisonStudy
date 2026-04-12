library(CVXR)
library(pracma)

# ------------------------------------------------------------
# (1) Precompute matrix K and vector l that are used for the
#     estimated score loss J_hat = (1/2)*y^T*K*y - l^T*y
# ------------------------------------------------------------

# Precompute matrices that are later used for entry-wise 
# matrix operations instead of loop-wise executions
precompute_exponents_univariate <- function(m) {
  # Precompute a mxm matrix containing i+j at the (i,j)-th entry
  # Needed for second devirative
  expoN <- as.vector(outer(0:(m - 1), 0:(m - 1), `+`))
  # Add one to each entry for matrix M
  # Needed for first devirative
  expoM <- expoN + 1
  # (Redundant) Save same values for the denominators in M
  denomM <- expoM
  # Return list of matrices for exponents used in M,N and
  # denominators in M
  list(
    expoN = expoN,
    expoM = expoM,
    denomM = denomM
  )
}

# Use precomputed exponents and denominators to build matrix
# K and vector l and a weighting function h and its derivative
build_K_and_l_univariate <- function(
    z,
    m,
    h = function(z) rep(1, length(z)),
    h_prime = function(z) rep(0, length(z))
    ) {
      n <- length(z)
      z <- as.numeric(z)
      # Build exponents
      exps <- precompute_exponents_univariate(m)
      # For every entry in z apply exponents in expoM (treated as vector)
      # Result is a matrix containing z_i^expoM_j as (i,j)-th element
      # --> Matrix with n rows and mxm columns
      A <- outer(z, exps$expoM, `^`)
      # Do the same with expoN
      B <- outer(z, exps$expoN, `^`)
      # For the first derivative we need to apply on each column j the
      # denominator the debominator denomM[j]
      A <- sweep(A, 2, exps$denomM, "/")
      
      # Get weights by applying function h and h_prime
      hx  <- as.numeric(h(z))
      hpx <- as.numeric(h_prime(z))
      # Only allow positive values for h
      # TODO: Check if this make sense
      hx <- pmax(hx, 0)
      
      # Use column scaling for numerical stability (std for each column)
      # as big m can be very unstable for exponentials
      # TODO: Mention in Report
      scale_vec <- apply(A, 2, function(v) max(stats::sd(v), 1e-12))
      # Apply coulmn scaling
      A_sc <- sweep(A, 2, scale_vec, "/")
      B_sc <- sweep(B, 2, scale_vec, "/")
      
      # Calculate aggregated moments
      Aw <- A_sc * sqrt(hx)
      S  <- crossprod(Aw) / n
      t  <- as.vector(crossprod(A_sc, hx) / n)
      u  <- sum(hx) / n
      r  <- as.vector(crossprod(A_sc, hpx) / n)
      q  <- sum(hpx) / n
      b_bar <- as.vector(crossprod(B_sc, hx) / n)
      
      # Use the moments to derive matrix K and vector l
      K <- rbind(
        cbind(S, matrix(t, ncol = 1)),
        cbind(matrix(t, nrow = 1), matrix(u, 1, 1))
      )
      K <- 0.5 * (K + t(K))
      l <- c(r + b_bar, q)
      
      # Return K, l und the column scaling vector
      list(K = K, l = l, scale_vec = scale_vec)
}


# ------------------------------------------------------------
# (2) Core fit: univariate Score Matching
# ------------------------------------------------------------
# - Use identity function as standard for weighting function h

# TODO: ridge evt. entfernen
fit_score_matching_univariate <- function(
    x,
    m,
    h = function(z) rep(1, length(z)),
    h_prime = function(z) rep(0, length(z)),
    standardize = TRUE,
    ridge = 0,
    solver = "SCS",
    scs_control = list(max_iters = 100000, eps = 1e-5, alpha = 1.8, verbose = FALSE)
) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2) stop("Data vector x must contain at least 
                          two finite observations.")
  if (m < 1) stop("Degree m must be >= 1.")
  # TODO:
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
    stop("ridge must be a nonnegative finite scalar.")
  }
  
  # optional standardization
  if (standardize) {
    sc <- scale_vector_1d(x)
    z <- sc$x_scaled
  } else {
    z <- x
    sc <- list(center = 0, scale = 1)
  }
  
  # Precompute matrix K and vector l which determine the estimated score loss
  score_loss_input <- build_K_and_l_univariate(z, m, h, h_prime)
  
  # Specify optimization variables G, c1 and set PSD constraint to True 
  G  <- Variable(m, m, PSD = TRUE)
  c1 <- Variable(1)
  # Transform these variables to vector representation
  p = m*m
  gvec <- tryCatch(vec(G), error = function(e) reshape(G, c(p, 1)))
  y <- vstack(gvec, c1)
  
  # TODO:
  K_raw <- score_loss_input$K
  K_reg <- K_raw + ridge * diag(nrow(K_raw))
  
  # TODO:
  eig_raw <- tryCatch(
    eigen(K_raw, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, nrow(K_raw))
  )
  
  eig_reg <- tryCatch(
    eigen(K_reg, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, nrow(K_reg))
  )
  
  kappa_raw <- tryCatch(kappa(K_raw, exact = TRUE), error = function(e) NA_real_)
  kappa_reg <- tryCatch(kappa(K_reg, exact = TRUE), error = function(e) NA_real_)
  
  rcond_raw <- tryCatch(1 / kappa(K_raw, exact = TRUE), error = function(e) NA_real_)
  rcond_reg <- tryCatch(1 / kappa(K_reg, exact = TRUE), error = function(e) NA_real_)
  
  # Define Score loss as objective function by using the precomputed
  # matrix K and vector l
  obj <- 0.5 * quad_form(y, Constant(K_reg)) -
    t(Constant(matrix(score_loss_input$l, ncol = 1))) %*% y
  
  # Specify CVXR problem
  prob <- Problem(Minimize(obj))
  
  # Solve CVXR problem using provided solver settings
  sol <- solve(
    prob,
    solver = solver,
    max_iters = scs_control$max_iters,
    eps = scs_control$eps,
    alpha = scs_control$alpha,
    verbose = scs_control$verbose
  )
  
  # Back-transform column scaled solution for g_vec
  g_sc  <- as.numeric(sol$getValue(gvec))
  g_org <- g_sc / score_loss_input$scale_vec
  # Back-transform vector representation of g to Matrix G
  G_org <- matrix(g_org, nrow = m, ncol = m)
  
  # TODO: why using structure and class?
  structure(
    list(
      G = G_org,
      c1 = as.numeric(sol$getValue(c1)),
      status = sol$status,
      solution = sol,
      scaling = list(
        center = sc$center,
        scale = sc$scale,
        standardize = standardize
      ),
      column_scaling = list(
        scale_vec = score_loss_input$scale_vec
      ),
      # TODO:
      ridge = ridge,
      diagnostics = list(
        kappa_raw = kappa_raw,
        kappa_reg = kappa_reg,
        rcond_raw = rcond_raw,
        rcond_reg = rcond_reg,
        eigmin_raw = if (all(is.na(eig_raw))) NA_real_ else min(eig_raw, na.rm = TRUE),
        eigmin_reg = if (all(is.na(eig_reg))) NA_real_ else min(eig_reg, na.rm = TRUE)
      ),
      # END TODO
      z_train = z,
      m = m
    ),
    class = "score_matching_univariate_fit"
  )
}


# ------------------------------------------------------------
# (3) Evaluation metrics
# ------------------------------------------------------------

# ------------------------------------------------------------
# (3.1) Evaluate Score on original scale of x
# ------------------------------------------------------------
predict_score_z_univariate <- function(z, fit) {
  # TODO: Desccribe procedure
  z <- as.numeric(z)
  m <- fit$m
  exps <- precompute_exponents_univariate(m)
  A <- outer(z, exps$expoM, `^`)
  A <- sweep(A, 2, exps$denomM, "/")
  gvec <- as.vector(fit$G)
  as.numeric(A %*% gvec) + fit$c1
}

predict_score_univariate <- function(x, fit) {
  z <- apply_scaling_1d(x, fit$scaling)$z
  score_z <- predict_score_z_univariate(z, fit)
  # TODO: Check this
  # r_z(z) = s * r_x(x)  =>  r_x(x) = r_z(z) / s
  score_z / fit$scaling$scale
}


# ------------------------------------------------------------
# (4) Stop: New Code to check!!!
# ------------------------------------------------------------

# ------------------------------------------------------------
# (3.2) Reconstruct univariate density from the fitted score
# ------------------------------------------------------------
# Theory:
#   r(z) = - d/dz log f_Z(z)
# implies
#   log f_Z(z) = C - integral r(u) du.
# We reconstruct an unnormalized log-density in the standardized
# variable z, estimate the log-normalizing constant numerically,
# and then map back to the original x-scale via
#   f_X(x) = f_Z((x-mu)/s) / s.

# Extract polynomial coefficients a_0, ..., a_{2m-1} of the score
# r(z) = a_0 + a_1 z + ... + a_{2m-1} z^{2m-1}.
get_score_polynomial_coefficients_univariate <- function(fit) {
  m <- fit$m
  G <- fit$G
  coeffs <- numeric(2L * m)
  coeffs[1L] <- fit$c1
  
  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      deg <- i + j - 1L
      coeffs[deg + 1L] <- coeffs[deg + 1L] + G[i, j] / deg
    }
  }
  coeffs
}

# Evaluate the antiderivative A(z) = integral_0^z r(u) du.
# Then log f(z) = -A(z) up to an additive constant.
predict_antiderivative_score_z_univariate <- function(z, fit) {
  z <- as.numeric(z)
  coeffs <- get_score_polynomial_coefficients_univariate(fit)
  degs <- 0:(length(coeffs) - 1L)
  out <- numeric(length(z))
  for (k in seq_along(coeffs)) {
    out <- out + coeffs[k] * z^(degs[k] + 1L) / (degs[k] + 1L)
  }
  out
}

# Unnormalized log-density on standardized z-scale.
predict_logdensity_unnormalized_z_univariate <- function(z, fit) {
  -predict_antiderivative_score_z_univariate(z, fit)
}

# Helper to build a reasonable finite search interval around the
# standardized training sample.
get_default_density_bounds_z_univariate <- function(fit,
                                                    n_sd = 8,
                                                    min_half_width = 8) {
  z_train <- as.numeric(fit$z_train)
  if (length(z_train) == 0L || any(!is.finite(z_train))) {
    return(c(-min_half_width, min_half_width))
  }
  mu <- mean(z_train)
  s  <- stats::sd(z_train)
  if (!is.finite(s) || s <= 0) s <- 1
  lo <- min(min(z_train), mu - n_sd * s, -min_half_width)
  hi <- max(max(z_train), mu + n_sd * s,  min_half_width)
  c(lo, hi)
}

# Approximate mode on z-scale to stabilize numerical integration.
find_mode_logdensity_z_univariate <- function(fit,
                                              interval = NULL,
                                              grid_length = 2001L) {
  if (is.null(interval)) {
    interval <- get_default_density_bounds_z_univariate(fit)
  }
  interval <- as.numeric(interval)
  if (length(interval) != 2L || !all(is.finite(interval)) || interval[1L] >= interval[2L]) {
    stop("interval must be a finite vector (lower, upper) with lower < upper.")
  }
  
  f_log <- function(z) predict_logdensity_unnormalized_z_univariate(z, fit)
  grid <- seq(interval[1L], interval[2L], length.out = grid_length)
  vals <- f_log(grid)
  idx <- which.max(vals)
  z0 <- grid[idx]
  
  # refine locally if possible
  left  <- if (idx <= 1L) interval[1L] else grid[idx - 1L]
  right <- if (idx >= length(grid)) interval[2L] else grid[idx + 1L]
  opt <- tryCatch(
    optimize(function(z) -f_log(z), interval = c(left, right), maximum = FALSE),
    error = function(e) NULL
  )
  if (!is.null(opt)) {
    z0 <- opt$minimum
  }
  list(mode = z0, logvalue = f_log(z0), interval = interval)
}

# Numerically compute log Z_z where
#   Z_z = integral exp(log f_Z(z)) dz.
# A shift by the approximate mode is used for stability.
compute_log_normalizer_z_univariate <- function(fit,
                                                interval = NULL,
                                                subdivisions = 200L,
                                                rel.tol = 1e-8,
                                                abs.tol = 0,
                                                stop_on_failure = FALSE) {
  if (is.null(interval)) {
    interval <- get_default_density_bounds_z_univariate(fit)
  }
  mode_info <- find_mode_logdensity_z_univariate(fit, interval = interval)
  shift <- mode_info$logvalue
  
  integrand <- function(z) {
    val <- predict_logdensity_unnormalized_z_univariate(z, fit) - shift
    exp(val)
  }
  
  # TODO: Intervalgrenzen sin nicht -inf, inf warum vertretbar?
  if (is.null(interval)) {
    interval <- get_default_density_bounds_z_univariate(fit, n_sd = 8, min_half_width = 8)
  }
  interval <- as.numeric(interval)
  
  # TODO: Intervalgrenzen sin nicht -inf, inf warum vertretbar?
  integ <- tryCatch(
    stats::integrate(
      f = integrand,
      lower = interval[1L],
      upper = interval[2L],
      subdivisions = subdivisions,
      rel.tol = rel.tol,
      abs.tol = abs.tol,
      stop.on.error = stop_on_failure
    ),
    error = function(e) e
  )
  
  if (inherits(integ, "error") || !is.list(integ) || !is.finite(integ$value) || integ$value <= 0) {
    msg <- paste(
      "Failed to compute a finite normalizing constant for the reconstructed density.",
      "This can happen if the fitted score does not induce an integrable log-density.")
    if (isTRUE(stop_on_failure)) stop(msg)
    return(list(logZ = NA_real_, value = NA_real_, shift = shift, message = msg,
                mode = mode_info$mode, interval = interval))
  }
  
  list(
    logZ = as.numeric(log(integ$value) + shift),
    value = as.numeric(integ$value * exp(shift)),
    shift = shift,
    abs.error = integ$abs.error,
    subdivisions = subdivisions,
    mode = mode_info$mode,
    interval = interval,
    message = if (!is.null(integ$message)) integ$message else NULL
  )
}

# Normalized log-density on standardized z-scale.
predict_logdensity_z_univariate <- function(z,
                                            fit,
                                            eps = 1e-300,
                                            interval = NULL,
                                            subdivisions = 200L,
                                            rel.tol = 1e-8,
                                            abs.tol = 0,
                                            stop_on_failure = FALSE) {
  z <- as.numeric(z)
  z <- z[is.finite(z)]
  
  log_unnorm <- predict_logdensity_unnormalized_z_univariate(z, fit)
  norm_info <- compute_log_normalizer_z_univariate(
    fit = fit,
    interval = interval,
    subdivisions = subdivisions,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    stop_on_failure = stop_on_failure
  )
  
  if (!is.finite(norm_info$logZ)) {
    out <- rep(NA_real_, length(z))
    return(out)
  }
  
  log_unnorm - norm_info$logZ
}

predict_density_z_univariate <- function(z,
                                         fit,
                                         interval = NULL,
                                         subdivisions = 200L,
                                         rel.tol = 1e-8,
                                         abs.tol = 0,
                                         stop_on_failure = FALSE) {
  logdens <- predict_logdensity_z_univariate(
    z = z,
    fit = fit,
    interval = interval,
    subdivisions = subdivisions,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    stop_on_failure = stop_on_failure
  )
  exp(logdens)
}

# Normalized log-density on original x-scale.
predict_logdensity_univariate <- function(x,
                                          fit,
                                          eps = 1e-300,
                                          interval = NULL,
                                          subdivisions = 200L,
                                          rel.tol = 1e-8,
                                          abs.tol = 0,
                                          stop_on_failure = FALSE) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  z <- apply_scaling_1d(x, fit$scaling)$z
  logdens_z <- predict_logdensity_z_univariate(
    z = z,
    fit = fit,
    eps = eps,
    interval = interval,
    subdivisions = subdivisions,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    stop_on_failure = stop_on_failure
  )
  
  ifelse(is.finite(logdens_z), logdens_z - log(fit$scaling$scale), NA_real_)
}

predict_density_univariate <- function(x,
                                       fit,
                                       interval = NULL,
                                       subdivisions = 200L,
                                       rel.tol = 1e-8,
                                       abs.tol = 0,
                                       stop_on_failure = FALSE) {
  logdens_x <- predict_logdensity_univariate(
    x = x,
    fit = fit,
    interval = interval,
    subdivisions = subdivisions,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    stop_on_failure = stop_on_failure
  )
  exp(logdens_x)
}

# Aliases aligned with the naming convention in the KDE / MLE scripts.
predict_logdensity_sm_1d <- function(newx, fit, ...) {
  predict_logdensity_univariate(newx, fit, ...)
}

predict_density_sm_1d <- function(newx, fit, ...) {
  predict_density_univariate(newx, fit, ...)
}



