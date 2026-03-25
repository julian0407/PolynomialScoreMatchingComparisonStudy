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
#
# H_method:
#   "Hpi"   = plug-in bandwidth (recommended main baseline)
#   "Hlscv" = least-squares cross-validation
#   "Hns"   = normal-scale / Scott-like simple rule
#
# diagonal = TRUE can be useful when components are known independent
#            or when you want a simpler constrained baseline
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
  
  # Optional diagonal restriction
  if (diagonal) {
    H <- diag(diag(H), nrow = d, ncol = d)
  }
  
  # small ridge for numerical stability
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
#
# newx: m x d matrix of evaluation points
# returns vector length m
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
    diffs <- sweep(x, 2, newx[j, ], FUN = "-")      # n x d, Xi - x
    quad  <- rowSums((diffs %*% H_inv) * diffs)     # Mahalanobis^2
    out[j] <- const * sum(exp(-0.5 * quad))
  }
  
  out
}

# ------------------------------------------------------------
# Negative gradient of log KDE density
#
# Returns m x d matrix:
#   r_hat(x) = -∇ log \hat f(x)
#
# This matches your convention where score_true is
# negative derivative / negative gradient of log p.
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
    diffs <- sweep(x, 2, newx[j, ], FUN = "-")      # Xi - x
    quad  <- rowSums((diffs %*% H_inv) * diffs)
    w     <- exp(-0.5 * quad)
    
    dens_j <- const * sum(w)
    dens_j <- max(dens_j, eps_dens)
    
    # ∇ f_hat(x) = const * sum_i w_i * H^{-1}(Xi - x)
    grad_f <- const * colSums((w * (diffs %*% H_inv)))
    
    # negative gradient of log density
    out[j, ] <- -grad_f / dens_j
  }
  
  colnames(out) <- paste0("dim", seq_len(d))
  out
}

# ------------------------------------------------------------
# Optional Monte Carlo integrated absolute error
#
# If p_true can evaluate true density on points:
#   IAE ≈ volume * mean(|fhat - ftrue|) on evaluation cloud
#
# For moderate d, user-supplied evaluation cloud is easier than grids.
# ------------------------------------------------------------
mc_L1_error_mv <- function(eval_points, p_hat, p_true) {
  eval_points <- as_matrix_nd(eval_points)
  p_ref <- p_true(eval_points)
  if (length(p_ref) != nrow(eval_points)) stop("p_true must return one value per row.")
  mean(abs(p_hat - p_ref))
}

# ------------------------------------------------------------
# Direct multivariate score loss
#
# score_true(newx) must return m x d matrix with rows
# equal to -∇ log p_true(x)
#
# h can be scalar weights or a function from m x d matrix -> vector length m
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
  t(solve(Sigma, t(xc)))   # = Sigma^{-1}(x-mu), rowwise
}

# ------------------------------------------------------------
# Minimal comparison example
# ------------------------------------------------------------
set.seed(1)

d <- 2
mu <- c(0, 0)
Sigma <- matrix(c(1, 0.4,
                  0.4, 1.5), 2, 2)

n_train <- 1000
n_test  <- 5000

x_train <- rmvnorm_base(n_train, mu, Sigma)
x_test  <- rmvnorm_base(n_test,  mu, Sigma)

# KDE baselines
fit_kde_pi   <- fit_kde_mv(x_train, H_method = "Hpi", diagonal = FALSE)
fit_kde_lscv <- fit_kde_mv(x_train, H_method = "Hlscv", diagonal = FALSE)

# Score-loss against truth
sl_pi <- score_loss_direct_kde_mv(
  x_test, fit_kde_pi,
  score_true = function(x) score_mvn_true(x, mu, Sigma)
)

sl_lscv <- score_loss_direct_kde_mv(
  x_test, fit_kde_lscv,
  score_true = function(x) score_mvn_true(x, mu, Sigma)
)

cat("Direct score-loss KDE (Hpi):  ", sl_pi, "\n")
cat("Direct score-loss KDE (Hlscv):", sl_lscv, "\n")

# Density values on test points
p_hat_test <- predict_kde_density_mv(x_test, fit_kde_pi)
p_true_test <- dmvnorm_base(x_test, mu, Sigma)

cat("Mean absolute density error on test cloud:",
    mean(abs(p_hat_test - p_true_test)), "\n")