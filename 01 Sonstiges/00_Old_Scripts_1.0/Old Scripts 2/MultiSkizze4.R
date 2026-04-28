# ============================================================
# Multivariate KDE baseline for comparison
# Requires: install.packages("ks")
# ============================================================

library(ks)

# ------------------------------------------------------------
# Helper: ensure matrix shape
# x: n x d
# ------------------------------------------------------------
as_matrix_nd <- function(x) {
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("x must be numeric.")
  if (ncol(x) < 1) stop("x must have at least one column.")
  x
}

# ------------------------------------------------------------
# Fit multivariate Gaussian KDE
# ------------------------------------------------------------
fit_kde_mv <- function(x,
                       H_method = c("Hpi", "Hlscv", "Hns"),
                       diagonal = FALSE,
                       unit_interval = FALSE) {
  H_method <- match.arg(H_method)
  x <- as_matrix_nd(x)
  d <- ncol(x)
  
  H <- switch(
    H_method,
    Hpi   = ks::Hpi(x = x),
    Hlscv = ks::Hlscv(x = x),
    Hns   = ks::Hns(x = x)
  )
  
  if (diagonal) {
    H <- diag(diag(H), nrow = d, ncol = d)
  }
  
  eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) <= 0) {
    H <- H + diag(abs(min(eig)) + 1e-8, d)
  }
  
  list(
    x = x,
    H = H,
    H_inv = solve(H),
    detH = det(H),
    d = d,
    n = nrow(x),
    H_method = H_method,
    diagonal = diagonal
  )
}

# ------------------------------------------------------------
# Multivariate Gaussian KDE density evaluation
# ------------------------------------------------------------
predict_kde_density_mv <- function(newx, fit) {
  newx <- as_matrix_nd(newx)
  if (ncol(newx) != fit$d) stop("Dimension mismatch in newx.")
  
  x <- fit$x
  H_inv <- fit$H_inv
  detH <- fit$detH
  d <- fit$d
  n <- fit$n
  
  const <- 1 / (((2 * pi)^(d / 2)) * sqrt(detH) * n)
  out <- numeric(nrow(newx))
  
  for (j in seq_len(nrow(newx))) {
    diffs <- sweep(x, 2, newx[j, ], FUN = "-")
    quad  <- rowSums((diffs %*% H_inv) * diffs)
    out[j] <- const * sum(exp(-0.5 * quad))
  }
  
  out
}

# ------------------------------------------------------------
# Negative gradient of log KDE density
# ------------------------------------------------------------
predict_kde_score_mv <- function(newx, fit, eps_dens = 1e-300) {
  newx <- as_matrix_nd(newx)
  if (ncol(newx) != fit$d) stop("Dimension mismatch in newx.")
  
  x <- fit$x
  H_inv <- fit$H_inv
  detH <- fit$detH
  d <- fit$d
  n <- fit$n
  
  const <- 1 / (((2 * pi)^(d / 2)) * sqrt(detH) * n)
  out <- matrix(NA_real_, nrow = nrow(newx), ncol = d)
  
  for (j in seq_len(nrow(newx))) {
    diffs <- sweep(x, 2, newx[j, ], FUN = "-")
    quad  <- rowSums((diffs %*% H_inv) * diffs)
    w     <- exp(-0.5 * quad)
    
    dens_j <- const * sum(w)
    dens_j <- max(dens_j, eps_dens)
    
    grad_f <- const * colSums((w * (diffs %*% H_inv)))
    out[j, ] <- -grad_f / dens_j
  }
  
  colnames(out) <- paste0("dim", seq_len(d))
  out
}

# ------------------------------------------------------------
# Optional Monte Carlo integrated absolute error
# ------------------------------------------------------------
mc_L1_error_mv <- function(eval_points, p_hat, p_true) {
  eval_points <- as_matrix_nd(eval_points)
  p_ref <- p_true(eval_points)
  if (length(p_ref) != nrow(eval_points)) stop("p_true must return one value per row.")
  mean(abs(p_hat - p_ref))
}

# ------------------------------------------------------------
# Direct multivariate score loss
# ------------------------------------------------------------
score_loss_direct_kde_mv <- function(x_test, fit_kde, score_true,
                                     h = function(x) rep(1, nrow(x))) {
  x_test <- as_matrix_nd(x_test)
  r_hat  <- predict_kde_score_mv(x_test, fit_kde)
  r_true <- score_true(x_test)
  
  if (!is.matrix(r_true)) r_true <- as.matrix(r_true)
  if (!all(dim(r_true) == dim(r_hat))) {
    stop("score_true must return a matrix with same dimension as x_test.")
  }
  
  w <- h(x_test)
  if (length(w) != nrow(x_test)) stop("h must return one weight per row.")
  w <- pmax(as.numeric(w), 0)
  
  sqerr <- rowSums((r_hat - r_true)^2)
  mean(0.5 * w * sqerr)
}

# ------------------------------------------------------------
# Example: multivariate Gaussian truth
# ------------------------------------------------------------
rmvnorm_base <- function(n, mu, Sigma) {
  Z <- matrix(rnorm(n * length(mu)), nrow = n)
  L <- chol(Sigma)
  sweep(Z %*% L, 2, mu, FUN = "+")
}

dmvnorm_base <- function(x, mu, Sigma) {
  x <- as_matrix_nd(x)
  d <- ncol(x)
  xc <- sweep(x, 2, mu, FUN = "-")
  Sigma_inv <- solve(Sigma)
  detS <- det(Sigma)
  const <- 1 / (((2 * pi)^(d / 2)) * sqrt(detS))
  quad <- rowSums((xc %*% Sigma_inv) * xc)
  const * exp(-0.5 * quad)
}

score_mvn_true <- function(x, mu, Sigma) {
  x <- as_matrix_nd(x)
  xc <- sweep(x, 2, mu, FUN = "-")
  t(solve(Sigma, t(xc)))
}

# ============================================================
# Pairwise polynomial score matching basis builder
# ============================================================

build_pairwise_score_matching <- function(
    x, m,
    build_Phi = TRUE,
    drop_constant = TRUE,
    standardize = TRUE,
    center = NULL,
    scale = NULL
) {
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("x must be numeric.")
  n <- nrow(x)
  d <- ncol(x)
  if (m < 1) stop("m must be >= 1.")
  
  if (standardize) {
    if (is.null(center)) center <- colMeans(x)
    if (is.null(scale))  scale  <- apply(x, 2, sd)
    scale[!is.finite(scale) | scale <= 0] <- 1
    
    x_work <- sweep(x, 2, center, FUN = "-")
    x_work <- sweep(x_work, 2, scale, FUN = "/")
  } else {
    center <- rep(0, d)
    scale  <- rep(1, d)
    x_work <- x
  }
  
  uni_deg <- 0:(2 * m - 2)
  pair_deg <- 0:(m - 1)
  
  basis <- vector("list", 0)
  
  for (l in seq_len(d)) {
    for (r in uni_deg) {
      if (drop_constant && r == 0) next
      basis[[length(basis) + 1L]] <- list(
        type = "uni",
        l = l,
        r = r
      )
    }
  }
  
  for (l in seq_len(d - 1L)) {
    for (u in (l + 1L):d) {
      for (i in pair_deg) {
        for (j in pair_deg) {
          if (drop_constant && i == 0 && j == 0) next
          basis[[length(basis) + 1L]] <- list(
            type = "pair",
            l = l,
            u = u,
            i = i,
            j = j
          )
        }
      }
    }
  }
  
  p <- length(basis)
  max_deg <- 2 * m - 2
  Xpow <- vector("list", d)
  
  for (l in seq_len(d)) {
    P <- matrix(1, nrow = n, ncol = max_deg + 1L)
    if (max_deg >= 1) {
      P[, 2] <- x_work[, l]
      if (max_deg >= 2) {
        for (r in 2:max_deg) {
          P[, r + 1L] <- P[, r] * x_work[, l]
        }
      }
    }
    Xpow[[l]] <- P
  }
  
  get_pow <- function(l, deg) {
    if (deg < 0) return(rep(0, n))
    Xpow[[l]][, deg + 1L]
  }
  
  Phi <- if (build_Phi) matrix(0, nrow = n, ncol = p) else NULL
  D <- lapply(seq_len(d), function(k) matrix(0, nrow = n, ncol = p))
  Lap <- matrix(0, nrow = n, ncol = p)
  
  for (col in seq_len(p)) {
    b <- basis[[col]]
    
    if (b$type == "uni") {
      l <- b$l
      r <- b$r
      
      if (build_Phi) {
        Phi[, col] <- get_pow(l, r)
      }
      
      if (r >= 1) {
        D[[l]][, col] <- r * get_pow(l, r - 1L)
      }
      
      if (r >= 2) {
        Lap[, col] <- r * (r - 1) * get_pow(l, r - 2L)
      }
      
    } else {
      l <- b$l
      u <- b$u
      i <- b$i
      j <- b$j
      
      xl_i <- get_pow(l, i)
      xu_j <- get_pow(u, j)
      
      if (build_Phi) {
        Phi[, col] <- xl_i * xu_j
      }
      
      if (i >= 1) {
        D[[l]][, col] <- i * get_pow(l, i - 1L) * xu_j
      }
      
      if (j >= 1) {
        D[[u]][, col] <- j * xl_i * get_pow(u, j - 1L)
      }
      
      term_l <- if (i >= 2) i * (i - 1) * get_pow(l, i - 2L) * xu_j else 0
      term_u <- if (j >= 2) j * (j - 1) * xl_i * get_pow(u, j - 2L) else 0
      Lap[, col] <- term_l + term_u
    }
  }
  
  K <- matrix(0, nrow = p, ncol = p)
  for (k in seq_len(d)) {
    K <- K + crossprod(D[[k]])
  }
  K <- K / n
  
  ell <- colMeans(Lap)
  
  basis_names <- vapply(basis, function(b) {
    if (b$type == "uni") {
      sprintf("z%d^%d", b$l, b$r)
    } else {
      sprintf("z%d^%d*z%d^%d", b$l, b$i, b$u, b$j)
    }
  }, character(1))
  
  colnames(K) <- rownames(K) <- basis_names
  names(ell) <- basis_names
  if (!is.null(Phi)) colnames(Phi) <- basis_names
  
  list(
    x_original = x,
    x_used = x_work,
    center = center,
    scale = scale,
    standardized = standardize,
    Phi = Phi,
    K = K,
    ell = ell,
    basis = basis,
    basis_names = basis_names,
    D = D,
    Lap = Lap,
    m = m
  )
}

transform_with_fit <- function(newx, fit_obj) {
  newx <- as.matrix(newx)
  sweep(sweep(newx, 2, fit_obj$center, "-"), 2, fit_obj$scale, "/")
}

# ============================================================
# Unconstrained multivariate pairwise polynomial score matching
# ============================================================

fit_pairwise_score_matching_unconstrained <- function(
    x, m,
    standardize = TRUE,
    center = NULL,
    scale = NULL,
    drop_constant = TRUE,
    ridge = 0,
    use_ridge = TRUE,
    ridge_exclude_constant = TRUE,
    build_Phi = FALSE,
    solver = c("solve", "qr")
) {
  solver <- match.arg(solver)
  
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) stop("x must be numeric.")
  if (nrow(x) < 2) stop("Need at least 2 observations.")
  if (ncol(x) < 1) stop("Need at least 1 column.")
  if (m < 1) stop("m must be >= 1.")
  if (!is.numeric(ridge) || length(ridge) != 1 || ridge < 0) {
    stop("ridge must be a nonnegative scalar.")
  }
  
  prep <- build_pairwise_score_matching(
    x = x,
    m = m,
    build_Phi = build_Phi,
    drop_constant = drop_constant,
    standardize = standardize,
    center = center,
    scale = scale
  )
  
  K <- prep$K
  ell <- prep$ell
  p <- ncol(K)
  
  pen_diag <- rep(1, p)
  
  if (ridge_exclude_constant) {
    is_const <- vapply(prep$basis, function(b) {
      if (b$type == "uni") {
        b$r == 0
      } else {
        b$i == 0 && b$j == 0
      }
    }, logical(1))
    pen_diag[is_const] <- 0
  }
  
  lambda_eff <- if (use_ridge) ridge else 0
  K_reg <- K + lambda_eff * diag(pen_diag, nrow = p, ncol = p)
  
  theta <- switch(
    solver,
    solve = tryCatch(
      solve(K_reg, ell),
      error = function(e) {
        warning("solve() failed; falling back to qr.solve().")
        qr.solve(K_reg, ell)
      }
    ),
    qr = qr.solve(K_reg, ell)
  )
  
  fitted_score_train <- do.call(
    cbind,
    lapply(seq_along(prep$D), function(k) prep$D[[k]] %*% theta)
  )
  
  if (prep$standardized) {
    fitted_score_train <- sweep(fitted_score_train, 2, prep$scale, "/")
  }
  
  colnames(fitted_score_train) <- paste0("dim", seq_len(ncol(x)))
  
  objective_value <- as.numeric(
    0.5 * crossprod(theta, prep$K %*% theta) - crossprod(prep$ell, theta)
  )
  
  objective_value_reg <- as.numeric(
    objective_value + 0.5 * lambda_eff * sum(pen_diag * theta^2)
  )
  
  out <- list(
    theta = as.numeric(theta),
    K = prep$K,
    ell = prep$ell,
    K_reg = K_reg,
    ridge = ridge,
    use_ridge = use_ridge,
    lambda_eff = lambda_eff,
    ridge_penalty_diag = pen_diag,
    objective_value = objective_value,
    objective_value_reg = objective_value_reg,
    fitted_score_train = fitted_score_train,
    prep = prep,
    m = m,
    d = ncol(x),
    n = nrow(x),
    drop_constant = drop_constant,
    standardize = standardize,
    solver = solver,
    call = match.call()
  )
  
  class(out) <- "pairwise_score_matching_fit"
  out
}

# ------------------------------------------------------------
# Evaluate derivative design matrices on new data
# ------------------------------------------------------------
evaluate_pairwise_derivatives <- function(newx, prep) {
  newx <- as.matrix(newx)
  if (!is.numeric(newx)) stop("newx must be numeric.")
  if (ncol(newx) != ncol(prep$x_original)) {
    stop("Dimension mismatch in newx.")
  }
  
  z <- transform_with_fit(newx, prep)
  n <- nrow(z)
  d <- ncol(z)
  p <- length(prep$basis)
  max_deg <- 2 * prep$m - 2
  
  Xpow <- vector("list", d)
  for (l in seq_len(d)) {
    P <- matrix(1, nrow = n, ncol = max_deg + 1L)
    if (max_deg >= 1) {
      P[, 2] <- z[, l]
      if (max_deg >= 2) {
        for (r in 2:max_deg) {
          P[, r + 1L] <- P[, r] * z[, l]
        }
      }
    }
    Xpow[[l]] <- P
  }
  
  get_pow <- function(l, deg) {
    if (deg < 0) return(rep(0, n))
    Xpow[[l]][, deg + 1L]
  }
  
  D_new <- lapply(seq_len(d), function(k) matrix(0, nrow = n, ncol = p))
  
  for (col in seq_len(p)) {
    b <- prep$basis[[col]]
    
    if (b$type == "uni") {
      l <- b$l
      r <- b$r
      
      if (r >= 1) {
        D_new[[l]][, col] <- r * get_pow(l, r - 1L)
      }
      
    } else {
      l <- b$l
      u <- b$u
      i <- b$i
      j <- b$j
      
      xl_i <- get_pow(l, i)
      xu_j <- get_pow(u, j)
      
      if (i >= 1) {
        D_new[[l]][, col] <- i * get_pow(l, i - 1L) * xu_j
      }
      if (j >= 1) {
        D_new[[u]][, col] <- j * xl_i * get_pow(u, j - 1L)
      }
    }
  }
  
  D_new
}

# ------------------------------------------------------------
# Predict negative gradient of log density
# ------------------------------------------------------------
predict_pairwise_score_matching <- function(newx, fit) {
  if (!inherits(fit, "pairwise_score_matching_fit")) {
    stop("fit must be of class 'pairwise_score_matching_fit'.")
  }
  
  newx <- as.matrix(newx)
  D_new <- evaluate_pairwise_derivatives(newx, fit$prep)
  
  out_z <- do.call(
    cbind,
    lapply(seq_along(D_new), function(k) D_new[[k]] %*% fit$theta)
  )
  
  out_x <- sweep(out_z, 2, fit$prep$scale, "/")
  
  colnames(out_x) <- paste0("dim", seq_len(ncol(out_x)))
  out_x
}

score_loss_direct_pairwise_mv <- function(
    x_test, fit, score_true,
    h = function(x) rep(1, nrow(x))
) {
  x_test <- as.matrix(x_test)
  r_hat  <- predict_pairwise_score_matching(x_test, fit)
  r_true <- score_true(x_test)
  
  if (!is.matrix(r_true)) r_true <- as.matrix(r_true)
  if (!all(dim(r_true) == dim(r_hat))) {
    stop("score_true must return a matrix with same dimension as x_test.")
  }
  
  w <- h(x_test)
  if (length(w) != nrow(x_test)) stop("h must return one weight per row.")
  w <- pmax(as.numeric(w), 0)
  
  sqerr <- rowSums((r_hat - r_true)^2)
  mean(0.5 * w * sqerr)
}

# ============================================================
# DENSITY RECONSTRUCTION FOR PAIRWISE SCORE MATCHING
# via importance sampling in standardized z-space
# ============================================================

# ------------------------------------------------------------
# log-sum-exp helper
# ------------------------------------------------------------
logsumexp <- function(v) {
  vmax <- max(v)
  vmax + log(sum(exp(v - vmax)))
}

# ------------------------------------------------------------
# log density of multivariate Gaussian
# ------------------------------------------------------------
log_dmvnorm_base <- function(x, mu, Sigma) {
  x <- as_matrix_nd(x)
  d <- ncol(x)
  xc <- sweep(x, 2, mu, FUN = "-")
  Sigma_inv <- solve(Sigma)
  logdetS <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  quad <- rowSums((xc %*% Sigma_inv) * xc)
  -0.5 * (d * log(2 * pi) + logdetS + quad)
}

# ------------------------------------------------------------
# Evaluate basis matrix Phi on new standardized points z
# ------------------------------------------------------------
evaluate_pairwise_basis_z <- function(z, prep) {
  z <- as.matrix(z)
  if (!is.numeric(z)) stop("z must be numeric.")
  
  n <- nrow(z)
  d <- ncol(z)
  p <- length(prep$basis)
  max_deg <- 2 * prep$m - 2
  
  Xpow <- vector("list", d)
  for (l in seq_len(d)) {
    P <- matrix(1, nrow = n, ncol = max_deg + 1L)
    if (max_deg >= 1) {
      P[, 2] <- z[, l]
      if (max_deg >= 2) {
        for (r in 2:max_deg) {
          P[, r + 1L] <- P[, r] * z[, l]
        }
      }
    }
    Xpow[[l]] <- P
  }
  
  get_pow <- function(l, deg) {
    if (deg < 0) return(rep(0, n))
    Xpow[[l]][, deg + 1L]
  }
  
  Phi <- matrix(0, nrow = n, ncol = p)
  
  for (col in seq_len(p)) {
    b <- prep$basis[[col]]
    
    if (b$type == "uni") {
      Phi[, col] <- get_pow(b$l, b$r)
    } else {
      Phi[, col] <- get_pow(b$l, b$i) * get_pow(b$u, b$j)
    }
  }
  
  colnames(Phi) <- prep$basis_names
  Phi
}

# ------------------------------------------------------------
# Evaluate energy S(z) = theta^T phi(z) on z-scale
# where score = grad_z S(z)
# ------------------------------------------------------------
predict_pairwise_energy_z <- function(z, fit) {
  z <- as.matrix(z)
  Phi <- evaluate_pairwise_basis_z(z, fit$prep)
  as.vector(Phi %*% fit$theta)
}

# ------------------------------------------------------------
# Evaluate energy on x-scale by first standardizing to z
# ------------------------------------------------------------
predict_pairwise_energy <- function(newx, fit) {
  newx <- as.matrix(newx)
  z <- transform_with_fit(newx, fit$prep)
  predict_pairwise_energy_z(z, fit)
}

# ------------------------------------------------------------
# Estimate log normalizing constant in z-space:
#   Z_z = integral exp(-S(z)) dz
#
# via importance sampling with proposal q(z)
# default proposal: N(0, I)
# ------------------------------------------------------------
estimate_pairwise_log_normalizer <- function(
    fit,
    n_mc = 50000,
    proposal_mu = NULL,
    proposal_Sigma = NULL,
    seed = NULL,
    return_details = TRUE
) {
  if (!inherits(fit, "pairwise_score_matching_fit")) {
    stop("fit must be of class 'pairwise_score_matching_fit'.")
  }
  
  d <- fit$d
  
  if (is.null(proposal_mu)) {
    proposal_mu <- rep(0, d)
  }
  if (is.null(proposal_Sigma)) {
    proposal_Sigma <- diag(d)
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  z_mc <- rmvnorm_base(n_mc, proposal_mu, proposal_Sigma)
  log_q <- log_dmvnorm_base(z_mc, proposal_mu, proposal_Sigma)
  S_mc  <- predict_pairwise_energy_z(z_mc, fit)
  
  log_w <- -S_mc - log_q
  logZ  <- logsumexp(log_w) - log(n_mc)
  
  # normalized importance weights for diagnostics
  lw_shift <- log_w - max(log_w)
  w_raw <- exp(lw_shift)
  w_norm <- w_raw / sum(w_raw)
  ess <- 1 / sum(w_norm^2)
  
  out <- list(
    logZ = as.numeric(logZ),
    n_mc = n_mc,
    proposal_mu = proposal_mu,
    proposal_Sigma = proposal_Sigma,
    ess = ess,
    ess_ratio = ess / n_mc
  )
  
  if (return_details) {
    out$z_mc <- z_mc
    out$log_w <- log_w
  }
  
  out
}

# ------------------------------------------------------------
# Predict normalized density on x-scale
#
# p_X(x) = p_Z(z) / prod(scale)
# log p_X(x) = -S(z) - logZ_z - sum(log(scale))
# ------------------------------------------------------------
predict_pairwise_density <- function(
    newx, fit,
    logZ_info = NULL,
    n_mc = 50000,
    proposal_mu = NULL,
    proposal_Sigma = NULL,
    seed = NULL
) {
  if (!inherits(fit, "pairwise_score_matching_fit")) {
    stop("fit must be of class 'pairwise_score_matching_fit'.")
  }
  
  newx <- as.matrix(newx)
  
  if (is.null(logZ_info)) {
    logZ_info <- estimate_pairwise_log_normalizer(
      fit = fit,
      n_mc = n_mc,
      proposal_mu = proposal_mu,
      proposal_Sigma = proposal_Sigma,
      seed = seed,
      return_details = FALSE
    )
  }
  
  z <- transform_with_fit(newx, fit$prep)
  S_new <- predict_pairwise_energy_z(z, fit)
  
  log_jac <- sum(log(fit$prep$scale))
  log_px <- -S_new - logZ_info$logZ - log_jac
  
  exp(log_px)
}

# ============================================================
# Minimal comparison example
# ============================================================

set.seed(1)

d <- 2
mu <- c(0, 0)
Sigma <- matrix(c(1, 0.4,
                  0.4, 1.5), 2, 2)

n_train <- 1000
n_test  <- 5000

x_train <- rmvnorm_base(n_train, mu, Sigma)
x_test  <- rmvnorm_base(n_test,  mu, Sigma)

# ------------------------------------------------------------
# Pairwise score matching:
# for initial testing use ridge everywhere
# ------------------------------------------------------------
fit_pair_ridge_small <- fit_pairwise_score_matching_unconstrained(
  x = x_train,
  m = 3,
  standardize = TRUE,
  ridge = 1e-3,
  use_ridge = TRUE
)

fit_pair_ridge_large <- fit_pairwise_score_matching_unconstrained(
  x = x_train,
  m = 3,
  standardize = TRUE,
  ridge = 1e-2,
  use_ridge = TRUE
)

# ------------------------------------------------------------
# KDE baselines
# ------------------------------------------------------------
fit_kde_pi   <- fit_kde_mv(x_train, H_method = "Hpi", diagonal = FALSE)
fit_kde_lscv <- fit_kde_mv(x_train, H_method = "Hlscv", diagonal = FALSE)

# ------------------------------------------------------------
# True score for comparison
# ------------------------------------------------------------
score_true_fun <- function(x) score_mvn_true(x, mu, Sigma)

# ------------------------------------------------------------
# Direct score-loss
# ------------------------------------------------------------
sl_pair_ridge_small <- score_loss_direct_pairwise_mv(
  x_test, fit_pair_ridge_small,
  score_true = score_true_fun
)

sl_pair_ridge_large <- score_loss_direct_pairwise_mv(
  x_test, fit_pair_ridge_large,
  score_true = score_true_fun
)

sl_kde_pi <- score_loss_direct_kde_mv(
  x_test, fit_kde_pi,
  score_true = score_true_fun
)

sl_kde_lscv <- score_loss_direct_kde_mv(
  x_test, fit_kde_lscv,
  score_true = score_true_fun
)

cat("Direct score-loss pairwise SM (ridge=1e-3): ", sl_pair_ridge_small, "\n")
cat("Direct score-loss pairwise SM (ridge=1e-2): ", sl_pair_ridge_large, "\n")
cat("Direct score-loss KDE (Hpi):                 ", sl_kde_pi, "\n")
cat("Direct score-loss KDE (Hlscv):               ", sl_kde_lscv, "\n")

# ------------------------------------------------------------
# Density comparison via importance sampling normalization
# ------------------------------------------------------------
logZ_pair_small <- estimate_pairwise_log_normalizer(
  fit = fit_pair_ridge_small,
  n_mc = 50000,
  proposal_mu = rep(0, d),
  proposal_Sigma = diag(d),
  seed = 1,
  return_details = FALSE
)

logZ_pair_large <- estimate_pairwise_log_normalizer(
  fit = fit_pair_ridge_large,
  n_mc = 50000,
  proposal_mu = rep(0, d),
  proposal_Sigma = diag(d),
  seed = 1,
  return_details = FALSE
)

cat("Pairwise SM log-normalizer ESS ratio (ridge=1e-3): ",
    logZ_pair_small$ess_ratio, "\n")
cat("Pairwise SM log-normalizer ESS ratio (ridge=1e-2): ",
    logZ_pair_large$ess_ratio, "\n")

p_pair_small <- predict_pairwise_density(
  x_test,
  fit = fit_pair_ridge_small,
  logZ_info = logZ_pair_small
)

p_pair_large <- predict_pairwise_density(
  x_test,
  fit = fit_pair_ridge_large,
  logZ_info = logZ_pair_large
)

p_kde_pi   <- predict_kde_density_mv(x_test, fit_kde_pi)
p_kde_lscv <- predict_kde_density_mv(x_test, fit_kde_lscv)
p_true     <- dmvnorm_base(x_test, mu, Sigma)

L1_pair_small <- mean(abs(p_pair_small - p_true))
L1_pair_large <- mean(abs(p_pair_large - p_true))
L1_kde_pi     <- mean(abs(p_kde_pi - p_true))
L1_kde_lscv   <- mean(abs(p_kde_lscv - p_true))

cat("L1 error pairwise SM (ridge=1e-3): ", L1_pair_small, "\n")
cat("L1 error pairwise SM (ridge=1e-2): ", L1_pair_large, "\n")
cat("L1 error KDE (Hpi):                ", L1_kde_pi, "\n")
cat("L1 error KDE (Hlscv):              ", L1_kde_lscv, "\n")

# ------------------------------------------------------------
# Inspect first few score predictions
# ------------------------------------------------------------
head(
  cbind(
    x_test[1:6, ],
    predict_pairwise_score_matching(x_test[1:6, ], fit_pair_ridge_small),
    score_true_fun(x_test[1:6, ])
  )
)

# ------------------------------------------------------------
# Inspect first few density predictions
# ------------------------------------------------------------
head(
  cbind(
    x_test[1:6, ],
    p_pair_small = p_pair_small[1:6],
    p_kde_pi = p_kde_pi[1:6],
    p_true = p_true[1:6]
  )
)