library(CVXR)
library(pracma)
library(extraDistr)
library(sn)

# ============================================================
# Robustness helpers: scaling + stable exp normalization
# ============================================================

scale_x <- function(x, clip = 6) {
  mu <- mean(x)
  s  <- sd(x)
  if (!is.finite(s) || s <= 0) s <- 1
  z <- (x - mu) / s
  if (!is.null(clip)) z <- pmax(pmin(z, clip), -clip)
  list(z = z, mu = mu, s = s, clip = clip)
}

stable_exp_neg <- function(v) {
  shift <- min(v)
  list(values = exp(-(v - shift)), shift = shift)
}

# -----------------------------
# Build M(x) and N(x) matrices
# -----------------------------
make_MN <- function(x, m) {
  i0 <- 0:(m - 1)
  j0 <- 0:(m - 1)
  I  <- matrix(rep(i0, times = m), nrow = m, byrow = FALSE)
  J  <- matrix(rep(j0, each  = m), nrow = m, byrow = FALSE)
  
  N <- x^(I + J)
  M <- x^(I + J + 1) / (I + J + 1)
  
  list(M = M, N = N)
}

# ============================================================
# Robust score matching fit (UNCHANGED)
# ============================================================
fit_score_matching_matrixG_robust <- function(
    x,
    m,
    lambda_trace = 1e-2,
    lambda_frob  = 1e-3,
    eps = 1e-6,
    standardize = TRUE,
    clip = 6,
    col_scale = c("maxabs", "sd"),
    solver = "SCS",
    verbose = TRUE
) {
  col_scale <- match.arg(col_scale)
  
  sc <- list(mu = 0, s = 1, clip = NULL)
  x_used <- x
  if (standardize) {
    sc <- scale_x(x, clip = clip)
    x_used <- sc$z
  }
  
  n <- length(x_used)
  p <- m * m
  
  t_pre <- system.time({
    A <- matrix(0, nrow = n, ncol = p)
    B <- matrix(0, nrow = n, ncol = p)
    
    for (k in seq_len(n)) {
      MN <- make_MN(x_used[k], m)
      A[k, ] <- as.vector(MN$M)
      B[k, ] <- as.vector(MN$N)
    }
    b_bar <- colMeans(B)
  })
  
  scale_vec <- switch(
    col_scale,
    maxabs = apply(abs(A), 2, function(v) max(v, 1e-12)),
    sd     = apply(A, 2, function(v) max(sd(v), 1e-12))
  )
  
  A_scaled    <- sweep(A, 2, scale_vec, "/")
  bbar_scaled <- b_bar / scale_vec
  
  t_build <- system.time({
    G  <- Variable(m, m, PSD = TRUE)
    c1 <- Variable(1)
    
    gvec <- tryCatch(vec(G), error = function(e) reshape(G, c(p, 1)))
    
    A_c  <- Constant(A_scaled)
    ones <- Constant(matrix(1, nrow = n, ncol = 1))
    b_c  <- Constant(matrix(bbar_scaled, ncol = 1))
    
    s1 <- A_c %*% gvec + c1 * ones
    
    obj <- (0.5 / n) * sum_squares(s1) -
      t(b_c) %*% gvec +
      lambda_trace * sum(diag(G))
    
    prob <- Problem(
      Minimize(obj),
      constraints = list(
        G == t(G),
        diag(G) >= eps
      )
    )
  })
  
  t_solve <- system.time({
    sol <- solve(prob, solver = solver)
  })
  
  gvec_sol  <- as.numeric(sol$getValue(gvec))
  gvec_orig <- gvec_sol / scale_vec
  G_orig <- matrix(gvec_orig, nrow = m, ncol = m)
  G_orig <- 0.5 * (G_orig + t(G_orig))
  
  if (verbose) {
    cat("\n--- Timing ---\n")
    cat(sprintf("Precompute (A,B):   %.3f sec\n", t_pre[["elapsed"]]))
    cat(sprintf("Build (CVXR):       %.3f sec\n", t_build[["elapsed"]]))
    cat(sprintf("Solve (solver):     %.3f sec\n", t_solve[["elapsed"]]))
    cat("Status:", sol$status, "\n")
    cat("--------------\n\n")
    
    if (standardize) {
      cat(sprintf("Standardization: mu=%.4f, sd=%.4f, clip=%s\n",
                  sc$mu, sc$s, ifelse(is.null(sc$clip), "NULL", as.character(sc$clip))))
    }
    cat(sprintf("Regularization: lambda_trace=%.2e, lambda_frob=%.2e\n",
                lambda_trace, lambda_frob))
  }
  
  list(
    G  = G_orig,
    c1 = as.numeric(sol$getValue(c1)),
    solution = sol,
    timing = list(pre = t_pre, build = t_build, solve = t_solve),
    scaling = list(standardize = standardize, mu = sc$mu, s = sc$s, clip = sc$clip,
                   col_scale = col_scale, scale_vec = scale_vec)
  )
}

# -----------------------------------
# Reconstruct s(x) from G and c1
# -----------------------------------
s_function <- function(x, G, c1) {
  m <- nrow(G)
  val <- 0
  for (i in 0:(m - 1)) {
    for (j in 0:(m - 1)) {
      val <- val +
        G[i + 1, j + 1] *
        x^(i + j + 2) / ((i + j + 1) * (i + j + 2))
    }
  }
  val + c1 * x
}

# -----------------------------------
# Density based on score matching fit (stable)
# -----------------------------------
density_from_fit <- function(xgrid, G, c1, mu = 0, s = 1) {
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)
  
  zgrid <- (xgrid - mu) / s
  svals <- sapply(zgrid, s_function, G = G, c1 = c1)
  
  tmp <- stable_exp_neg(svals)
  unnorm <- tmp$values
  
  Z <- pracma::trapz(zgrid, unnorm)
  pz <- unnorm / Z
  
  (1 / s) * pz
}

# ============================================================
# ONLY NEW PART: Yeo-Johnson transform (works for all real x)
# ============================================================

yeo_johnson <- function(x, lambda = 0) {
  # hart: attributes + names weg, wirklich plain double
  x <- base::as.double(base::unname(x))
  lambda <- base::as.double(lambda)[1]
  
  out <- base::numeric(base::length(x))
  
  pos <- x >= 0
  if (base::any(pos)) {
    xp <- x[pos]
    out[pos] <- if (base::abs(lambda) > 1e-12) {
      (((xp + 1)^lambda) - 1) / lambda
    } else {
      base::log1p(xp)
    }
  }
  
  neg <- !pos
  if (base::any(neg)) {
    xn <- x[neg]  # xn < 0
    out[neg] <- if (base::abs(lambda - 2) > 1e-12) {
      -((((1 - xn)^(2 - lambda)) - 1) / (2 - lambda))
    } else {
      -base::log1p(-xn)
    }
  }
  
  out
}

yeo_johnson_deriv <- function(x, lambda = 0) {
  x <- base::as.double(base::unname(x))
  lambda <- base::as.double(lambda)[1]
  
  out <- base::numeric(base::length(x))
  pos <- x >= 0
  if (base::any(pos)) {
    xp <- x[pos]
    out[pos] <- if (base::abs(lambda) > 1e-12) (xp + 1)^(lambda - 1) else 1 / (xp + 1)
  }
  neg <- !pos
  if (base::any(neg)) {
    xn <- x[neg]
    out[neg] <- if (base::abs(lambda - 2) > 1e-12) (1 - xn)^(1 - lambda) else 1 / (1 - xn)
  }
  out
}

fit_score_matching_yj <- function(x, m, lambda_yj = 0, ...) {
  x <- base::as.double(base::unname(x))
  y <- yeo_johnson(x, lambda = lambda_yj)
  fit <- fit_score_matching_matrixG_robust(y, m, ...)
  fit$yj_lambda <- lambda_yj
  fit
}

# Wrapper: density back on original x-scale via Jacobian
density_from_fit_yj <- function(xgrid, fit) {
  ygrid <- yeo_johnson(xgrid, lambda = fit$yj_lambda)
  py <- density_from_fit(
    ygrid,
    fit$G,
    fit$c1,
    mu = fit$scaling$mu,
    s  = fit$scaling$s
  )
  px <- py * abs(yeo_johnson_deriv(xgrid, lambda = fit$yj_lambda))
  px
}

# ============================================================
# Example run
# ============================================================
# set.seed(1)
# 
# # try any real-valued skewed data:
# x <- rgamma(10000, shape = 2, rate = 1) 
# x <- rsn(10000, xi = 0, omega = 1, alpha = 5)  # can be negative
# x <- rnorm(10000, mean = 50)
# 
# m <- 10
# 
# # Fit with ONE transformation (Yeo-Johnson)
# fit <- fit_score_matching_yj(
#   x, m,
#   lambda_yj = -0.7,         # change only this if you want (0 is a good default)
#   lambda_trace = 1e-2,
#   lambda_frob  = 1e-3,
#   eps = 1e-6,
#   standardize = TRUE,
#   clip = NULL,
#   col_scale = "maxabs",
#   solver = "SCS",
#   verbose = TRUE
# )
# 
# xgrid <- seq(min(x) - 1, max(x) + 1, length.out = 500)
# p_hat <- density_from_fit_yj(xgrid, fit)
# p_ref <- dgamma(xgrid, shape = 2, rate = 1)
# p_norm    <- dnorm(xgrid, 50, 1)
# 
# # reference for the example above:
# p_ref <- dsn(xgrid, xi = 0, omega = 1, alpha = 5)
# 
# hist(x, breaks = 30, freq = FALSE, col = "grey90", border = "grey70",
#      main = "Score matching with Yeo-Johnson transform", xlab = "x")
# lines(xgrid, p_hat,  lwd = 2, col = "blue")
# lines(xgrid, p_ref,  lwd = 2, col = "red", lty = 2)
# legend("topright",
#        legend = c("Score matching (YJ)", "Reference"),
#        col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n")
