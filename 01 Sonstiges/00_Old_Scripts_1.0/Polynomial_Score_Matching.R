library(CVXR)
library(pracma)
library(extraDistr)
library(sn)

# ============================================================
# Build M(x) and N(x) matrices assuming polynomials p(x), q(x)
# have degree m.
# ============================================================
make_MN <- function(x, m) {
  i0 <- 0:(m - 1)
  j0 <- 0:(m - 1)
  # times -> append vector, each -> append entry,
  # nrow -> structure of matrix (otherwise one large vector).
  I  <- matrix(rep(i0, times = m), nrow = m, byrow = FALSE)
  J  <- matrix(rep(j0, each  = m), nrow = m, byrow = FALSE)
  N <- x^(I + J)
  M <- x^(I + J + 1) / (I + J + 1)
  list(M = M, N = N)
}


# ============================================================
# Robust score matching fit
#   - optional x standardization
#   - column scaling of A and b_bar (preconditioning)
#   - strong PSD: G - eps I >> 0
# ============================================================
fit_score_matching_matrixG_robust <- function(
    x,
    m,
    h = function(x) rep(1, length(x)),
    h_prime = function(x) rep(0, length(x)),
    eps = 1e-6,
    standardize = TRUE,
    col_scale = c("maxabs", "sd"),
    solver = "SCS"
) {
  # Check if "maxabs" or "sd" was chosen as scaling operator
  col_scale <- match.arg(col_scale)
  
  # ====== (1) Standardize x to avoid huge powers ============
  sc <- list(mu = 0, s = 1)
  x_used <- x
  if (standardize) {
    sc <- scale_x(x)
    x_used <- sc$z
  }
  # matrices dimensions
  n <- length(x_used)
  p <- m * m
  
  # Testing
  cat(sprintf("Tau:   %.3f \n", log(length(x_used))))
  
  # ====== (2) Precompute A and B with sample data x ==================
  # ========== and generalizing function h ============================
  t_pre <- system.time({
    A <- matrix(0, nrow = n, ncol = p)
    B <- matrix(0, nrow = n, ncol = p)
    
    for (k in seq_len(n)) {
      MN <- make_MN(x_used[k], m)
      A[k, ] <- as.vector(MN$M)
      B[k, ] <- as.vector(MN$N)
    }
    hx <- as.numeric(h(x_used))
    stopifnot(length(hx) == n)
    hpx <- as.numeric(h_prime(x_used))
    stopifnot(length(hpx) == n)
    hp_c <- Constant(matrix(hpx, nrow = n, ncol = 1))
    B_h <- B*hx
    b_bar <- colMeans(B_h)
  })
  
  #  ====== (3) Column scaling (preconditioning) =============
  # scale columns of A and corresponding entries of b_bar,
  # b_bar_prime to similar magnitudes for robustness
  scale_vec <- switch(
    col_scale,
    maxabs = apply(abs(A), 2, function(v) max(v, 1e-12)),
    sd     = apply(A, 2, function(v) max(sd(v), 1e-12))
  )
  A_scaled    <- sweep(A, 2, scale_vec, "/")
  bbar_scaled <- b_bar / scale_vec
  
  #  ====== (4) Build CVXR problem ===========================
  t_build <- system.time({
    # PSD = positive semi definite
    G  <- Variable(m, m, PSD = TRUE)
    c1 <- Variable(1)
    
    # vec(G) robust across CVXR versions
    gvec <- tryCatch(vec(G), error = function(e) reshape(G, c(p, 1)))
    
    A_c  <- Constant(A_scaled)                                    # n x p
    ones <- Constant(matrix(1, nrow = n, ncol = 1))               # n x 1
    w_c  <- Constant(matrix(sqrt(pmax(hx, 1e-12)), ncol = 1))     # n x 1
    b_c  <- Constant(matrix(bbar_scaled, ncol = 1))               # p x 1
    hp_c <- Constant(matrix(hpx, ncol = 1))                       # p x 1
    
    s1 <- A_c %*% gvec + c1 * ones
    
    obj <- (0.5 / n) * sum_squares(multiply(w_c, s1)) -
      (1 / n) * t(hp_c) %*% s1 -
      t(b_c) %*% gvec
    
    prob <- Problem(Minimize(obj), constraints = list(G == t(G), diag(G) >= eps))
    
  })
  
  # ====== (5) Solve ==========================================
  t_solve <- system.time({
    sol <- solve(prob, solver = solver)
  })
  
  # ====== (6) Back-transform G because we solved =============
  # ======     with scaled columns ============================
  gvec_sol <- as.numeric(sol$getValue(gvec))
  gvec_orig <- gvec_sol / scale_vec
  G_orig <- matrix(gvec_orig, nrow = m, ncol = m)
  
  # symmetrize tiny numerical asymmetry
  # G_orig <- 0.5 * (G_orig + t(G_orig))

  # ====== (7) Print computation facts =========================
  cat("\n--- Timing ---\n")
  cat(sprintf("Precompute (A,B):   %.3f sec\n", t_pre[["elapsed"]]))
  cat(sprintf("Build (CVXR):       %.3f sec\n", t_build[["elapsed"]]))
  cat(sprintf("Solve (solver):     %.3f sec\n", t_solve[["elapsed"]]))
  cat("Status:", sol$status, "\n")
  cat("--------------\n\n")
  
  if (standardize) {
    cat(sprintf("Standardization: mu=%.4f, sd=%.4f\n",
                sc$mu, sc$s))
  }
  
  list(
    G  = G_orig,
    c1 = as.numeric(sol$getValue(c1)),
    solution = sol,
    timing = list(pre = t_pre, build = t_build, solve = t_solve),
    scaling = list(standardize = standardize, mu = sc$mu, s = sc$s,
                   col_scale = col_scale, scale_vec = scale_vec)
  )
}


# ============================================================
# Robustness helpers: scaling + stable exp normalization
# ============================================================
scale_x <- function(x) {
  mu <- mean(x)
  s  <- sd(x)
  if (!is.finite(s) || s <= 0) s <- 1
  z <- (x - mu) / s
  # ============================================================
  # Add eventually clip for outliers here
  # ============================================================
  list(z = z, mu = mu, s = s)
}

# ============================================================
# Reconstruct s(x) from G and c1
# ============================================================
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

# ============================================================
# Density based on score matching
# If fit was on z = (x - mu) / s, then
# p_X(x) = (1/s) * p_Z((x - mu)/s)
# ============================================================
density_from_fit <- function(xgrid, G, c1, mu = 0, s = 1) {
  # for m=1 error handling (must be matrux view)
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)
  
  zgrid <- (xgrid - mu) / s
  svals <- sapply(zgrid, s_function, G = G, c1 = c1)
  
  tmp <- stable_exp_neg(svals)
  unnorm <- tmp$values
  
  # Constant (Integral) for probability density
  Z <- pracma::trapz(zgrid, unnorm)
  pz <- unnorm / Z
  
  # Back transformation
  (1 / s) * pz
}

stable_exp_neg <- function(v) {
  # ============================================================
  # Why? -> For computing the density -see later
  # ============================================================ 
  shift <- min(v)
  list(values = exp(-(v - shift)), shift = shift)
}


# ============================================================
# Functions h and h'
# ============================================================
h_tanh_sq <- function(z) {
  tau <- log(length(z))
  u <- abs(z) / tau
  tanh(u)^2
}

h_tanh_sq_prime <- function(z) {
  tau <- log(length(z))
  u <- abs(z) / tau
  sech2 <- 1 / cosh(u)^2
  2 * tanh(u) * sech2 * (sign(z) / tau)
}

# ============================================================
# Example run
# ============================================================
# set.seed(1)

x <- rnorm(10000, mean = 50, sd = 10)
# x <- rnorm(20000, mean = 0, sd = 1 / sqrt(2))
# x <- rlogis(5000)
# x <- rlaplace(5000)

# x <- rsn(15000, xi = 0, omega = 1, alpha = 5)
# x <- rgamma(5000, shape = 2, rate = 1)

m <- 5

fit <- fit_score_matching_matrixG_robust(
  x = x,
  m = m,
  h = function(x) h_tanh_sq(x),
  h_prime = function(x) h_tanh_sq_prime(x),
  eps = 1e-6,
  standardize = TRUE
)

xgrid <- seq(min(x) - 1, max(x) + 1, length.out = 500)
p_hat <- density_from_fit(xgrid, fit$G, fit$c1, mu = fit$scaling$mu, s = fit$scaling$s)

p_ref    <- dnorm(xgrid, 50, 10)
# p_ref <- dnorm(xgrid, 0, 1 / sqrt(2))
# p_ref   <- dlogis(xgrid)
# p_ref <- dlaplace(xgrid)

# p_ref <- dsn(xgrid, xi = 0, omega = 1, alpha = 5)
# p_ref <- dgamma(xgrid, shape = 2, rate = 1)

hist(x, breaks = 30, freq = FALSE, col = "grey90", border = "grey70",
     main = "Score-matching log-concave density", xlab = "x")
lines(xgrid, p_hat,  lwd = 2, col = "green")
lines(xgrid, p_ref, lwd = 2, col = "red", lty = 2)
legend("topright",
       legend = c("Score matching estimate", "Reference"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n")

set.seed(1)

# ---- True model (nur für Vergleich; Fit bleibt allgemein) ----
mu0 <- 50; sd0 <- 10
r_true <- function(n) rnorm(n, mean = mu0, sd = sd0)
d_true <- function(x) dnorm(x, mean = mu0, sd = sd0)

# ---- experiment settings ----
m  <- 5
ns <- c(200, 500, 1000, 2000, 5000, 10000, 20000)
R  <- 20

# Speicher für Fehler
errs <- matrix(NA_real_, nrow = length(ns), ncol = R)

# Speicher für repräsentative Stichprobe + Fit pro n (für Plot)
x_store   <- vector("list", length(ns))
fit_store <- vector("list", length(ns))

# FIXES Evaluation grid (does NOT grow with n)
xg_eval <- seq(mu0 - 10 * sd0, mu0 + 10 * sd0, length.out = 2001)
p_ref_eval <- d_true(xg_eval)

for (i in seq_along(ns)) {
  n <- ns[i]
  
  # repräsentative Stichprobe für Density-Panel
  x_rep <- r_true(n)
  x_store[[i]] <- x_rep
  fit_store[[i]] <- fit_score_matching_matrixG_robust(
    x = x_rep,
    m = m,
    h = h_tanh_sq,
    h_prime = h_tanh_sq_prime,
    eps = 1e-6,
    standardize = TRUE,
    col_scale = "sd",
    solver = "SCS"
  )
  
  # L1 über R Wiederholungen
  for (r in 1:R) {
    x <- r_true(n)
    fit <- fit_score_matching_matrixG_robust(
      x = x,
      m = m,
      h = h_tanh_sq,
      h_prime = h_tanh_sq_prime,
      eps = 1e-6,
      standardize = TRUE,
      col_scale = "sd",
      solver = "SCS"
    )
    
    p_hat <- density_from_fit(xg_eval, fit$G, fit$c1, mu = fit$scaling$mu, s = fit$scaling$s)
    errs[i, r] <- L1_error(xg_eval, p_hat, p_ref_eval)
  }
  
  cat("n=", n,
      " status(rep)=", fit_store[[i]]$solution$status,
      " median L1=", median(errs[i, ], na.rm = TRUE), "\n")
}

# =========================
# Plot 1: L1-Konvergenz
# =========================
plot(ns, apply(errs, 1, median, na.rm = TRUE), log = "x", type = "b",
     pch = 19, col = "blue",
     xlab = "Sample size n",
     ylab = "Median L1 error",
     main = paste0("L1-Konvergenz des Score-Matching Schätzers (m=", m, ")"))
grid()

# =========================
# Plot 2: Density-Panels
# =========================
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))

for (i in seq_along(ns)) {
  n   <- ns[i]
  x   <- x_store[[i]]
  fit <- fit_store[[i]]
  
  # quantilbasiertes Grid (stabiler als min/max bei großem n)
  xgrid <- make_xgrid_quantile(x, ngrid = 600, qlo = 1e-3, qhi = 1 - 1e-3, expand = 0.15)
  
  p_hat <- density_from_fit(xgrid, fit$G, fit$c1, mu = fit$scaling$mu, s = fit$scaling$s)
  p_ref <- d_true(xgrid)
  
  hist(x, breaks = 30, freq = FALSE,
       col = "grey90", border = "grey70",
       main = paste0("n = ", n),
       xlab = "", ylab = "")
  
  lines(xgrid, p_hat, lwd = 2, col = "darkgreen")
  lines(xgrid, p_ref, lwd = 2, col = "red", lty = 2)
  
  legend("topright",
         legend = c("Score Matching", "True"),
         col = c("darkgreen", "red"),
         lwd = 2, lty = c(1, 2), bty = "n", cex = 0.8)
}

par(mfrow = c(1, 1))
