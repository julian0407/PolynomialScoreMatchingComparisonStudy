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
  # returns exp(-(v - min(v))) and the shift, avoids underflow
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
# Robust score matching fit
#   - optional x standardization (+ clipping)
#   - column scaling of A and b_bar (preconditioning)
#   - strong PSD: G - eps I >> 0
#   - trace + Frobenius (ridge) regularization
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
  
  # ---- (1) standardize x to avoid huge powers ----
  sc <- list(mu = 0, s = 1, clip = NULL)
  x_used <- x
  if (standardize) {
    sc <- scale_x(x, clip = clip)
    x_used <- sc$z
  }
  
  n <- length(x_used)
  p <- m * m
  
  # ---------- Precompute A and B ----------
  t_pre <- system.time({
    A <- matrix(0, nrow = n, ncol = p)
    B <- matrix(0, nrow = n, ncol = p)
    
    for (k in seq_len(n)) {
      MN <- make_MN(x_used[k], m)
      A[k, ] <- as.vector(MN$M)  # a_k = vec(Mk)
      B[k, ] <- as.vector(MN$N)  # b_k = vec(Nk)
    }
    b_bar <- colMeans(B)         # (1/n) sum_k b_k
  })
  
  # ---------- (2) Column scaling (preconditioning) ----------
  # scale columns of A and corresponding entries of b_bar to similar magnitudes
  scale_vec <- switch(
    col_scale,
    maxabs = apply(abs(A), 2, function(v) max(v, 1e-12)),
    sd     = apply(A, 2, function(v) max(sd(v), 1e-12))
  )
  
  A_scaled    <- sweep(A, 2, scale_vec, "/")
  bbar_scaled <- b_bar / scale_vec
  
  # ---------- Build CVXR problem ----------
  t_build <- system.time({
    G  <- Variable(m, m, PSD = TRUE)
    c1 <- Variable(1)
    
    # vec(G) robust across CVXR versions
    gvec <- tryCatch(vec(G), error = function(e) reshape(G, c(p, 1)))
    
    A_c  <- Constant(A_scaled)                        # n x p
    ones <- Constant(matrix(1, nrow = n, ncol = 1))   # n x 1
    b_c  <- Constant(matrix(bbar_scaled, ncol = 1))   # p x 1
    
    s1 <- A_c %*% gvec + c1 * ones
    
    obj <- (0.5 / n) * sum_squares(s1) - t(b_c) %*% gvec + lambda_trace * sum(diag(G))
    
    # # weights: center high, tails low
    # w_num <- 1 / (1 + x_used^2)     # sehr stabiler default
    # # optional normalize for scale:
    # w_num <- w_num / mean(w_num)
    # 
    # w <- Constant(matrix(w_num, ncol = 1))
    # s1w <- multiply(w, s1)          # elementwise
    # 
    # obj <- (0.5 / n) * sum_squares(s1w) -
    #   t(b_c) %*% gvec +
    #   lambda_trace * sum(diag(G)) +
    #   lambda_frob  * sum_squares(G)
    
    
    prob <- Problem(
      Minimize(obj),
      constraints = list(
        G == t(G),
        diag(G) >= eps
      )
    )
    
    
  })
  
  # ---------- Solve ----------
  t_solve <- system.time({
    sol <- solve(prob, solver = solver)
  })
  
  # ---------- Back-transform G because we solved with scaled columns ----------
  # We solved with A_scaled * gvec, where A_scaled = A / scale_vec.
  # Equivalent to original A * (gvec / scale_vec). Hence:
  # gvec_original = gvec_solution / scale_vec
  gvec_sol <- as.numeric(sol$getValue(gvec))
  gvec_orig <- gvec_sol / scale_vec
  G_orig <- matrix(gvec_orig, nrow = m, ncol = m)
  
  # symmetrize tiny numerical asymmetry
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
# Density based on score matching fit (stable + optional inverse transform)
# If fit was on z = (x - mu) / s, then
#   p_X(x) = (1/s) * p_Z((x - mu)/s)
# -----------------------------------
density_from_fit <- function(xgrid, G, c1, mu = 0, s = 1) {
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)  # for m=1
  
  zgrid <- (xgrid - mu) / s
  svals <- sapply(zgrid, s_function, G = G, c1 = c1)
  
  tmp <- stable_exp_neg(svals)
  unnorm <- tmp$values
  
  Z <- pracma::trapz(zgrid, unnorm)
  pz <- unnorm / Z
  
  (1 / s) * pz
}

# ============================================================
# Example run
# ============================================================
# set.seed(1)
# 
# x <- rnorm(10000, mean = 50)
# x <- rnorm(10000, mean = 0, sd = 1 / sqrt(2))
# x <- rlogis(1000)
# x <- rlaplace(5000)
# x <- rsn(5000, xi = 0, omega = 1, alpha = 5)
# 
# x2 <- sample(x, size = 1000)
# 
# # positive examples
# # x <- rexp(300, rate = 1)
# x <- rchisq(500, df = 5)
# x <- rgamma(3000, shape = 2, rate = 1)
# # x <- rbeta(300, shape1 = 2, shape2 = 3)
# 
# 
# m <- 10
# 
# fit <- fit_score_matching_matrixG_robust(
#   x, m,
#   lambda_trace = 1e-2,
#   lambda_frob  = 1e-3,
#   eps = 1e-6,
#   standardize = TRUE,
#   clip = 6,
#   col_scale = "maxabs",
#   solver = "SCS",
#   verbose = TRUE
# )
# 
# xgrid <- seq(min(x) - 1, max(x) + 1, length.out = 500)
# p_hat <- density_from_fit(xgrid, fit$G, fit$c1, mu = fit$scaling$mu, s = fit$scaling$s)
# 
# p_norm    <- dnorm(xgrid, 50, 1)
# p_squared <- dnorm(xgrid, 0, 1 / sqrt(2))
# p_logis   <- dlogis(xgrid)
# p_laplace <- dlaplace(xgrid)
# p_skewed_norm <- dsn(xgrid, xi = 0, omega = 1, alpha = 5)
# 
# p_exp   <- dexp(xgrid, rate = 1)
# p_chisq <- dchisq(xgrid, df = 5)
# p_gamma <- dgamma(xgrid, shape = 2, rate = 1)
# 
# hist(x, breaks = 30, freq = FALSE, col = "grey90", border = "grey70",
#      main = "Robust score-matching log-concave density", xlab = "x")
# lines(xgrid, p_hat,  lwd = 2, col = "blue")
# lines(xgrid, p_skewed_norm, lwd = 2, col = "red", lty = 2)
# legend("topright",
#        legend = c("Score matching estimate", "Reference"),
#        col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n")
