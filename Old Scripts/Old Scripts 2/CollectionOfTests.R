library(CVXR)
library(pracma)
library(sn)

# ============================================================
# (1) Basics: scale + stable exp
# ============================================================
# Scale sample with mean/sd -> serach for motivation why this is imporant
# for robustness in Ridge/Lasso
scale_x <- function(x) {
  mu <- mean(x)
  s  <- sd(x)
  if (!is.finite(s) || s <= 0) s <- 1
  list(z = (x - mu) / s, mu = mu, s = s)
}

stable_exp_neg <- function(v) {
  shift <- min(v)
  list(values = exp(-(v - shift)), shift = shift)
}

# ============================================================
# (2) Fast M/N features:
#       N_ij(z) = z^(i+j)
#       M_ij(z) = z^(i+j+1)/(i+j+1)
# ============================================================
# Precombute exponents/denominator
precompute_exponents <- function(m) {
  expoN <- as.vector(outer(0:(m-1), 0:(m-1), `+`))  # length p= m x m
  expoM <- expoN + 1
  denomM <- expoM                                  # divide by (i+j+1)
  list(expoN = expoN, expoM = expoM, denomM = denomM)
}

build_AB_fast <- function(z, m, exps) {
  # A[k,] = vec(M(z_k)), B[k,] = vec(N(z_k))
  A <- outer(z, exps$expoM, `^`)
  A <- sweep(A, 2, exps$denomM, "/")
  B <- outer(z, exps$expoN, `^`)
  list(A = A, B = B)
}

# ============================================================
# (3) Fast s(z) evaluation (no double for-loops)
#     s(z) = sum_{ij} G_ij * z^(i+j+2)/((i+j+1)(i+j+2)) + c1 z
# ============================================================
svals_fast <- function(zgrid, G, c1) {
  m <- nrow(G)
  expoN <- as.vector(outer(0:(m-1), 0:(m-1), `+`))
  expoS <- expoN + 2
  denomS <- (expoS - 1) * expoS
  coeff <- as.vector(G) / denomS
  
  Phi <- outer(zgrid, expoS, `^`)
  as.numeric(Phi %*% coeff) + c1 * zgrid # %*% -> sum of pairwise products
}

# ============================================================
# (4) Truncation for normalization: cheap + stable
#     - find mode on a grid
#     - expand left/right until s(z)-s(mode) >= target_drop
#  Find robust way to compute normalizing constant for density
#  reconstruction
# ============================================================
find_truncation_fast <- function(G, c1, z_hint, target_drop = 45,
                                 grid_len = 1201, max_expand = 30) {
  # initial window from data: robust and cheap
  lo <- as.numeric(quantile(z_hint, 1e-3))
  hi <- as.numeric(quantile(z_hint, 1 - 1e-3))
  span <- hi - lo
  if (!is.finite(span) || span <= 0) span <- sd(z_hint)
  lo <- lo - 0.5 * span
  hi <- hi + 0.5 * span
  
  z0_grid <- seq(lo, hi, length.out = grid_len)
  s0_grid <- svals_fast(z0_grid, G, c1)
  idx0 <- which.min(s0_grid)
  z0 <- z0_grid[idx0]
  s0 <- s0_grid[idx0]
  
  # expand until drop reached
  L <- lo; R <- hi
  for (k in 1:max_expand) {
    if (svals_fast(L, G, c1) - s0 < target_drop) L <- L * 2 else break
  }
  for (k in 1:max_expand) {
    if (svals_fast(R, G, c1) - s0 < target_drop) R <- R * 2 else break
  }
  list(left = L, right = R, z_mode = z0, s_mode = s0)
}


# ============================================================
# (4.1) Uses data that were used for model fit in case you want
# to plot different data -> use z_hint for truncation
# ============================================================
density_from_fit_fast_trunc <- function(xgrid, fit,
                                        target_drop = 45,
                                        n_int = 2500) {
  G  <- fit$G
  c1 <- fit$c1
  mu <- fit$scaling$mu
  s  <- fit$scaling$s
  
  zgrid <- (xgrid - mu) / s
  
  tr <- find_truncation_fast(G, c1, z_hint = fit$z_hint, target_drop = target_drop)
  zL <- tr$left; zR <- tr$right; s0 <- tr$s_mode
  
  # unnorm on z-scale with shift s0
  z_int <- seq(zL, zR, length.out = n_int)
  s_int <- svals_fast(z_int, G, c1)
  un_int <- exp(-(s_int - s0))
  # Get Z, note that s0 is elimenated in fraction
  Z <- trapz(z_int, un_int)
  
  s_g <- svals_fast(zgrid, G, c1)
  pz <- exp(-(s_g - s0)) / Z
  (1 / s) * pz
}

# ============================================================
# (5) FAST FIT via aggregated moments (statistically correct),
#     but structure + scaling like your old code
# ============================================================
fit_score_matching_matrixG_moments_fast <- function(
    x,
    m,
    h = function(z) rep(1, length(z)),
    h_prime = function(z) rep(0, length(z)),
    standardize = TRUE,
    col_scale = c("sd", "maxabs", "none"),
    eps = 1e-8,            # lower bound for diagonal entries
    delta_lead = 1e-6,     # coercivity (fix tails/center) G[m,m]> delta_lead
    ridge = 1e-8,          # L2 regularization
    solver = "SCS",
    scs_control = list(max_iters = 100000, eps = 1e-5, alpha = 1.8, verbose = FALSE),
    retry_if_inaccurate = TRUE
) {
  col_scale <- match.arg(col_scale)
  
  sc <- list(mu = 0, s = 1)
  z <- x
  if (standardize) {
    sc <- scale_x(x)
    z <- sc$z
  }
  
  n <- length(z)
  p <- m * m
  
  exps <- precompute_exponents(m)
  AB <- build_AB_fast(z, m, exps)
  A <- AB$A
  B <- AB$B
  
  hx  <- as.numeric(h(z))
  hpx <- as.numeric(h_prime(z))
  stopifnot(length(hx) == n, length(hpx) == n)
  hx <- pmax(hx, 0)      # set negative weights to zero
  
  # b_bar = E[h(z) * vec(N(z))]
  b_bar <- as.vector(crossprod(B, hx) / n)   # p-vector
  
  # column scaling like old code (but cheap)
  if (col_scale == "none") {
    scale_vec <- rep(1, p)
  } else if (col_scale == "maxabs") {
    scale_vec <- apply(abs(A), 2, function(v) max(v, 1e-12))
  } else { # "sd"
    scale_vec <- apply(A, 2, function(v) max(sd(v), 1e-12))
  }
  
  A_sc <- sweep(A, 2, scale_vec, "/")
  b_sc <- b_bar / scale_vec
  
  # build moments for objective
  # S = E[h * a a^T], t = E[h * a], u = E[h]
  # r = E[hp * a], q = E[hp]
  # all on scaled a
  Aw <- A_sc * sqrt(hx)              # n x p
  S  <- crossprod(Aw) / n            # p x p (fast BLAS)
  t  <- as.vector(crossprod(A_sc, hx)  / n)  # p
  u  <- sum(hx) / n
  r  <- as.vector(crossprod(A_sc, hpx) / n)  # p
  q  <- sum(hpx) / n
  
  # K and lin for 0.5 y^T K y - lin^T y, y=[g;c]
  K <- rbind(
    cbind(S, matrix(t, ncol = 1)),
    cbind(matrix(t, nrow = 1), matrix(u, 1, 1))
  )
  K <- 0.5 * (K + t(K))  # numerical symmetry
  
  lin <- c(r + b_sc, q)
  
  # CVXR variables (same as old: G PSD + c1)
  G  <- Variable(m, m, PSD = TRUE)
  c1 <- Variable(1)
  gvec <- tryCatch(vec(G), error = function(e) reshape(G, c(p, 1)))
  y <- vstack(gvec, c1)
  
  obj <- 0.5 * quad_form(y, Constant(K)) - t(Constant(matrix(lin, ncol = 1))) %*% y
  
  # Add ridge part
  if (ridge > 0) obj <- obj + (ridge / 2) * sum_squares(y)
  
  constr <- list(
    G == t(G),
    diag(G) >= eps,
    G[m, m] >= delta_lead
  )
  
  prob <- Problem(Minimize(obj), constr)
  
  sol <- solve(
    prob, solver = solver,
    max_iters = scs_control$max_iters,
    eps = scs_control$eps,
    alpha = scs_control$alpha,
    verbose = scs_control$verbose
  )
  
  if (retry_if_inaccurate && identical(sol$status, "optimal_inaccurate") && solver == "SCS") {
    sol <- solve(
      prob, solver = "SCS",
      max_iters = 2 * scs_control$max_iters,
      eps = scs_control$eps / 5,
      alpha = scs_control$alpha,
      verbose = scs_control$verbose
    )
  }
  
  # unscale g back: g_org = g_sc / scale_vec
  g_sc  <- as.numeric(sol$getValue(gvec))
  g_org <- g_sc / scale_vec
  G_org <- matrix(g_org, nrow = m, ncol = m)
  
  list(
    G = G_org,
    c1 = as.numeric(sol$getValue(c1)),
    status = sol$status,
    solution = sol,
    scaling = list(mu = sc$mu, s = sc$s, standardize = standardize,
                   col_scale = col_scale, scale_vec = scale_vec),
    # store a hint for truncation (z-data)
    z_hint = z
  )
}

# ============================================================
# (6) Robust h, h' (your originals)
# ============================================================
h_tanh_sq <- function(z) {
  # tau <- log(length(z))
  tau <- 1
  u <- abs(z) / tau
  tanh(u)^2
}
h_tanh_sq_prime <- function(z) {
  # tau <- log(length(z))
  tau <- 1
  u <- abs(z) / tau
  sech2 <- 1 / cosh(u)^2
  2 * tanh(u) * sech2 * (sign(z) / tau)
}

# ============================================================
# (7) Empirical convergence: L1 + density panels (fast)
# ============================================================
L1_error <- function(xg, p1, p2) trapz(xg, abs(p1 - p2))

set.seed(1)
mu0 <- 50; sd0 <- 10
r_true <- function(n) rnorm(n, mu0, sd0)
d_true <- function(x) dnorm(x, mu0, sd0)

ns <- c(200, 500, 1000, 2000, 5000, 10000, 20000, 100000)
R  <- 20
m  <- 1

errs <- matrix(NA_real_, nrow = length(ns), ncol = R)
x_store   <- vector("list", length(ns))
fit_store <- vector("list", length(ns))

# fixed evaluation grid (does NOT depend on min/max -> avoids tail inflation)
xg_eval <- seq(mu0 - 10 * sd0, mu0 + 10 * sd0, length.out = 2001)
p_ref_eval <- d_true(xg_eval)

for (i in seq_along(ns)) {
  n <- ns[i]
  
  x_rep <- r_true(n)
  x_store[[i]] <- x_rep
  fit_store[[i]] <- fit_score_matching_matrixG_moments_fast(
    x_rep, m = m,
    h = h_tanh_sq,
    h_prime = h_tanh_sq_prime,
    standardize = TRUE,
    col_scale = "sd",
    eps = 1e-8,
    delta_lead = 1e-6,
    ridge = 1e-8,
    scs_control = list(max_iters = 100000, eps = 1e-5, alpha = 1.8, verbose = FALSE)
  )
  
  for (r in 1:R) {
    x <- r_true(n)
    fit <- fit_score_matching_matrixG_moments_fast(
      x, m = m,
      h = h_tanh_sq,
      h_prime = h_tanh_sq_prime,
      standardize = TRUE,
      col_scale = "sd",
      eps = 1e-8,
      delta_lead = 1e-6,
      ridge = 1e-8
    )
    
    p_hat <- density_from_fit_fast_trunc(xg_eval, fit, target_drop = 45, n_int = 2500)
    errs[i, r] <- L1_error(xg_eval, p_hat, p_ref_eval)
  }
  
  cat("n=", n,
      " status(rep)=", fit_store[[i]]$status,
      " median L1=", median(errs[i, ], na.rm = TRUE), "\n")
}

# Plot 1: L1 convergence
plot(ns, apply(errs, 1, median, na.rm = TRUE), log="x", type="b",
     pch = 19, col = "blue",
     xlab = "Sample size n",
     ylab = "Median L1 error",
     main = paste0("L1-Konvergenz (moments-fast, m=", m, ")"))
grid()

# Plot 2: density panels
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
for (i in seq_along(ns)) {
  n   <- ns[i]
  x   <- x_store[[i]]
  fit <- fit_store[[i]]
  
  # panel grid: show a bit wider than eval grid but stable
  xgrid <- seq(mu0 - 6 * sd0, mu0 + 6 * sd0, length.out = 600)
  p_hat <- density_from_fit_fast_trunc(xgrid, fit, target_drop = 45, n_int = 2500)
  p_ref <- d_true(xgrid)
  
  hist(x, breaks = 30, freq = FALSE,
       col = "grey90", border = "grey70",
       main = paste0("n = ", n, "\n", fit$status),
       xlab = "", ylab = "")
  lines(xgrid, p_hat, lwd = 2, col = "darkgreen")
  lines(xgrid, p_ref, lwd = 2, col = "red", lty = 2)
}
par(mfrow = c(1, 1))

# ============================================================
# SINGLE EXPERIMENT (n = 10000) + Density Plot (FAST)
# ============================================================
set.seed(1)

n  <- 10000
m  <- 1
mu0 <- 50; sd0 <- 10

# Sample
x <- rnorm(n, mu0, sd0)

# Fit (your new fast moments estimator)
fit <- fit_score_matching_matrixG_moments_fast(
  x, m = m,
  h = h_tanh_sq,
  h_prime = h_tanh_sq_prime,
  standardize = TRUE,
  col_scale = "sd",
  eps = 1e-8,
  delta_lead = 1e-6,
  ridge = 1e-8,
  scs_control = list(max_iters = 80000, eps = 1e-5, alpha = 1.8, verbose = FALSE)
)

cat("Solver status:", fit$status, "\n")

# Stable plotting grid (DO NOT use min/max!)
xgrid <- seq(mu0 - 6 * sd0, mu0 + 6 * sd0, length.out = 600)

# Faster density (reduced integration cost)
p_hat <- density_from_fit_fast_trunc(
  xgrid, fit,
  target_drop = 35,  # faster, still accurate
  n_int = 800        # huge speedup vs 2500
)

# True density (for validation only)
p_ref <- dnorm(xgrid, mu0, sd0)

# Plot
hist(x, breaks = 40, freq = FALSE,
     col = "grey90", border = "grey70",
     main = paste0("Score Matching (moments-fast), n = ", n, ", m = ", m),
     xlab = "x")

lines(xgrid, p_hat, lwd = 2, col = "darkgreen")
lines(xgrid, p_ref, lwd = 2, col = "red", lty = 2)

legend("topright",
       legend = c("Score Matching", "True Density"),
       col = c("darkgreen", "red"),
       lwd = 2, lty = c(1, 2), bty = "n")






library(CVXR)
library(pracma)

# ============================================================
# (0) Deine Funktionen müssen schon im Workspace sein:
#     - fit_score_matching_matrixG_moments_fast
#     - density_from_fit_fast_trunc
#     - h_tanh_sq, h_tanh_sq_prime
#     - L1_error
# ============================================================

# ============================================================
# (1) Univariates Slice Sampling (Stepping-out + Shrinkage)
#     Neal (2003)-Style, funktioniert super für log-concave
# ============================================================
rslice1 <- function(n, logf, x0 = 0, w = 1, mstep = 50, burn = 200, thin = 1) {
  total <- burn + n * thin
  x <- numeric(total)
  x[1] <- x0
  
  for (t in 2:total) {
    xt <- x[t - 1]
    logy <- logf(xt) + log(runif(1))
    
    # Step out
    u <- runif(1)
    L <- xt - w * u
    R <- L + w
    J <- floor(runif(1, 0, mstep))
    K <- (mstep - 1) - J
    
    while (J > 0 && logf(L) > logy) { L <- L - w; J <- J - 1 }
    while (K > 0 && logf(R) > logy) { R <- R + w; K <- K - 1 }
    
    # Shrinkage
    repeat {
      xprop <- runif(1, L, R)
      if (logf(xprop) >= logy) { x[t] <- xprop; break }
      if (xprop < xt) L <- xprop else R <- xprop
    }
  }
  
  out <- x[(burn + 1):total]
  out <- out[seq(1, length(out), by = thin)]
  out
}

# ============================================================
# (2) Numerisch normalisierte "True"-Dichte auf Grid (stabil)
#     p(x) ∝ exp(logf(x))  => normalisiere über breites Intervall
# ============================================================
true_density_on_grid <- function(xgrid, logf, x_hint,
                                 target_drop = 45, n_int = 4000) {
  # robustes Suchintervall aus Datenhint
  lo <- as.numeric(quantile(x_hint, 1e-3))
  hi <- as.numeric(quantile(x_hint, 1 - 1e-3))
  span <- hi - lo
  if (!is.finite(span) || span <= 0) span <- sd(x_hint)
  lo <- lo - 0.5 * span
  hi <- hi + 0.5 * span
  
  # Mode auf grobem Grid finden
  z0_grid <- seq(lo, hi, length.out = 1201)
  lf0 <- vapply(z0_grid, logf, numeric(1))
  i0 <- which.max(lf0)
  x_mode <- z0_grid[i0]
  lf_mode <- lf0[i0]
  
  # links/rechts expandieren bis drop erreicht
  L <- lo; R <- hi
  for (k in 1:30) if ((lf_mode - logf(L)) < target_drop) L <- L * 2 else break
  for (k in 1:30) if ((lf_mode - logf(R)) < target_drop) R <- R * 2 else break
  
  # Normierungsintegral (Mode-shift stabil)
  x_int <- seq(L, R, length.out = n_int)
  lf_int <- vapply(x_int, logf, numeric(1))
  un_int <- exp(lf_int - lf_mode)
  Z <- trapz(x_int, un_int)
  
  # Dichte auf Grid
  lf_g <- vapply(xgrid, logf, numeric(1))
  exp(lf_g - lf_mode) / Z
}

# ============================================================
# (3) Zwei log-concave Ziel-Dichten, exakt in deiner Klasse
#     A: quartisch -> m=2
#     B: sextisch  -> m=3
# ============================================================
make_target_quartic <- function(a = 0.02, b = 0.5, c = 0.0) {
  stopifnot(a > 0, b >= 0)
  logf <- function(x) -(a * x^4 + b * x^2 + c * x)
  list(name = "Quartic: exp(-a x^4 - b x^2 - c x)", m = 2, logf = logf)
}

make_target_sextic <- function(a = 0.002, b = 0.02, d = 0.3, c = 0.0) {
  stopifnot(a > 0, b >= 0, d >= 0)
  logf <- function(x) -(a * x^6 + b * x^4 + d * x^2 + c * x)
  list(name = "Sextic: exp(-a x^6 - b x^4 - d x^2 - c x)", m = 3, logf = logf)
}

# ============================================================
# (4) 1:1 deine Konsistenz-Studie (L1 + Panels), nur "true" anders
# ============================================================
run_consistency_study <- function(target,
                                  ns = c(200, 500, 1000, 2000, 5000, 10000, 20000),
                                  R = 20,
                                  slice_w = 5, slice_mstep = 80,
                                  seed = 1) {
  set.seed(seed)
  
  m <- target$m
  logf <- target$logf
  
  errs <- matrix(NA_real_, nrow = length(ns), ncol = R)
  x_store   <- vector("list", length(ns))
  fit_store <- vector("list", length(ns))
  
  # Referenz-Grid: stabil; wir nehmen "Hint" aus großem Sample
  x_hint_big <- rslice1(n = 50000, logf = logf, x0 = 0, w = slice_w, mstep = slice_mstep, burn = 2000, thin = 2)
  mu_hint <- mean(x_hint_big); sd_hint <- sd(x_hint_big)
  xg_eval <- seq(mu_hint - 10 * sd_hint, mu_hint + 10 * sd_hint, length.out = 2001)
  p_ref_eval <- true_density_on_grid(xg_eval, logf, x_hint = x_hint_big, target_drop = 55, n_int = 5000)
  
  for (i in seq_along(ns)) {
    n <- ns[i]
    
    # "rep" Sample für Panels (wie bei dir)
    x_rep <- rslice1(n = n, logf = logf, x0 = 0, w = slice_w, mstep = slice_mstep, burn = 800, thin = 1)
    x_store[[i]] <- x_rep
    
    fit_store[[i]] <- fit_score_matching_matrixG_moments_fast(
      x_rep, m = m,
      h = h_tanh_sq,
      h_prime = h_tanh_sq_prime,
      standardize = TRUE,
      col_scale = "sd",
      eps = 1e-8,
      delta_lead = 1e-6,
      ridge = 1e-8,
      scs_control = list(max_iters = 100000, eps = 1e-5, alpha = 1.8, verbose = FALSE)
    )
    
    # R Wiederholungen (L1)
    for (r in 1:R) {
      x <- rslice1(n = n, logf = logf, x0 = 0, w = slice_w, mstep = slice_mstep, burn = 800, thin = 1)
      
      fit <- fit_score_matching_matrixG_moments_fast(
        x, m = m,
        h = h_tanh_sq,
        h_prime = h_tanh_sq_prime,
        standardize = TRUE,
        col_scale = "sd",
        eps = 1e-8,
        delta_lead = 1e-6,
        ridge = 1e-8
      )
      
      p_hat <- density_from_fit_fast_trunc(xg_eval, fit, target_drop = 45, n_int = 2500)
      errs[i, r] <- L1_error(xg_eval, p_hat, p_ref_eval)
    }
    
    cat("n=", n,
        " status(rep)=", fit_store[[i]]$status,
        " median L1=", median(errs[i, ], na.rm = TRUE), "\n")
  }
  
  # ---- Plot 1: L1 convergence ----
  plot(ns, apply(errs, 1, median, na.rm = TRUE), log="x", type="b",
       pch = 19, col = "blue",
       xlab = "Sample size n",
       ylab = "Median L1 error",
       main = paste0("L1-Konvergenz (", target$name, ", m=", m, ")"))
  grid()
  
  # ---- Plot 2: density panels ----
  par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
  for (i in seq_along(ns)) {
    n   <- ns[i]
    x   <- x_store[[i]]
    fit <- fit_store[[i]]
    
    # Panel-Grid stabil um Hint-Momente
    xgrid <- seq(mu_hint - 6 * sd_hint, mu_hint + 6 * sd_hint, length.out = 600)
    
    p_hat <- density_from_fit_fast_trunc(xgrid, fit, target_drop = 45, n_int = 2500)
    p_ref <- true_density_on_grid(xgrid, logf, x_hint = x_hint_big, target_drop = 55, n_int = 5000)
    
    hist(x, breaks = 30, freq = FALSE,
         col = "grey90", border = "grey70",
         main = paste0("n = ", n, "\n", fit$status),
         xlab = "", ylab = "")
    lines(xgrid, p_hat, lwd = 2, col = "darkgreen")
    lines(xgrid, p_ref, lwd = 2, col = "red", lty = 2)
  }
  par(mfrow = c(1, 1))
  
  invisible(list(errs = errs, ns = ns, target = target))
}

# ============================================================
# (5) RUN: Beispiel A (quartisch, m=2) und Beispiel B (sextisch, m=3)
# ============================================================

# A) Quartic (m=2)
targetA <- make_target_quartic(a = 0.02, b = 0.5, c = 0.0)
resA <- run_consistency_study(targetA, R = 20, seed = 1)

# B) Sextic (m=3)
targetB <- make_target_sextic(a = 0.002, b = 0.02, d = 0.3, c = 0.0)
resB <- run_consistency_study(targetB, R = 20, seed = 1)







# ============================================================
# TRUE logistic (location=0, scale=1)
# ============================================================
r_true <- function(n) rlogis(n, location = 0, scale = 1)
d_true <- function(x) dlogis(x, location = 0, scale = 1)

# stabile Grid-Wahl (nicht min/max!)
mu0 <- 0
sd0 <- pi / sqrt(3)   # Var(logistic)=pi^2/3
xg_eval <- seq(mu0 - 10 * sd0, mu0 + 10 * sd0, length.out = 2001)
p_ref_eval <- d_true(xg_eval)

# ============================================================
# Consistency / misspecification study: L1 + density panels
# ============================================================
set.seed(1)

ns <- c(200, 500, 1000, 2000, 5000, 10000, 20000)
R  <- 20

# WICHTIG: m=1 entspricht effektiv nur Gaussian-form -> logistic passt schlecht.
# Nimm m=2 oder m=3 um "flexibler" zu sein.
m <- 3   # <- probier auch 3

errs <- matrix(NA_real_, nrow = length(ns), ncol = R)
x_store   <- vector("list", length(ns))
fit_store <- vector("list", length(ns))

for (i in seq_along(ns)) {
  n <- ns[i]
  
  # ein Replikat für Panels
  x_rep <- r_true(n)
  x_store[[i]] <- x_rep
  
  fit_store[[i]] <- fit_score_matching_matrixG_moments_fast(
    x_rep, m = m,
    h = h_tanh_sq,
    h_prime = h_tanh_sq_prime,
    standardize = TRUE,
    col_scale = "sd",
    eps = 1e-8,
    delta_lead = 1e-6,
    ridge = 1e-8,
    scs_control = list(max_iters = 100000, eps = 1e-5, alpha = 1.8, verbose = FALSE)
  )
  
  # R Wiederholungen für L1
  for (r in 1:R) {
    x <- r_true(n)
    
    fit <- fit_score_matching_matrixG_moments_fast(
      x, m = m,
      h = h_tanh_sq,
      h_prime = h_tanh_sq_prime,
      standardize = TRUE,
      col_scale = "sd",
      eps = 1e-8,
      delta_lead = 1e-6,
      ridge = 1e-8
    )
    
    p_hat <- density_from_fit_fast_trunc(xg_eval, fit, target_drop = 45, n_int = 2500)
    errs[i, r] <- L1_error(xg_eval, p_hat, p_ref_eval)
  }
  
  cat("n=", n,
      " status(rep)=", fit_store[[i]]$status,
      " median L1=", median(errs[i, ], na.rm = TRUE), "\n")
}

# Plot 1: L1 convergence
plot(ns, apply(errs, 1, median, na.rm = TRUE), log="x", type="b",
     pch = 19, col = "blue",
     xlab = "Sample size n",
     ylab = "Median L1 error",
     main = paste0("L1 (Logistic misspecified), moments-fast, m=", m))
grid()

# Plot 2: density panels
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
for (i in seq_along(ns)) {
  n   <- ns[i]
  x   <- x_store[[i]]
  fit <- fit_store[[i]]
  
  xgrid <- seq(mu0 - 6 * sd0, mu0 + 6 * sd0, length.out = 600)
  
  p_hat <- density_from_fit_fast_trunc(xgrid, fit, target_drop = 45, n_int = 2500)
  p_ref <- d_true(xgrid)
  
  hist(x, breaks = 30, freq = FALSE,
       col = "grey90", border = "grey70",
       main = paste0("n = ", n, "\n", fit$status),
       xlab = "", ylab = "")
  lines(xgrid, p_hat, lwd = 2, col = "darkgreen")
  lines(xgrid, p_ref, lwd = 2, col = "red", lty = 2)
}
par(mfrow = c(1, 1))





# ============================================================
# TRUE: Skew-Normal (Azzalini) via sn::rsn/dsn
# ============================================================
xi    <- 0
omega <- 1
alpha <- 5

r_true <- function(n) rsn(n, xi = xi, omega = omega, alpha = alpha)
d_true <- function(x) dsn(x, xi = xi, omega = omega, alpha = alpha)

# Theoretische mean/sd (stabile Grid-Wahl, kein min/max!)
delta <- alpha / sqrt(1 + alpha^2)
mu0   <- xi + omega * delta * sqrt(2 / pi)
sd0   <- omega * sqrt(1 - (2 * delta^2) / pi)

# Eval grid für L1 (wie bei dir: +/- 10 sd)
xg_eval    <- seq(mu0 - 10 * sd0, mu0 + 10 * sd0, length.out = 2001)
p_ref_eval <- d_true(xg_eval)

# ============================================================
# Consistency / misspecification study: L1 + density panels
# ============================================================
set.seed(1)

ns <- c(200, 500, 1000, 2000, 5000, 10000, 20000, 100000)
R  <- 50

# Empfehlung: m >= 2 (m=1 ist praktisch "Gaussian-ish")
m <- 5  # probier auch 3

errs <- matrix(NA_real_, nrow = length(ns), ncol = R)
x_store   <- vector("list", length(ns))
fit_store <- vector("list", length(ns))

for (i in seq_along(ns)) {
  n <- ns[i]
  
  # ein Replikat für Panels
  x_rep <- r_true(n)
  x_store[[i]] <- x_rep
  
  fit_store[[i]] <- fit_score_matching_matrixG_moments_fast(
    x_rep, m = m,
    h = h_tanh_sq,
    h_prime = h_tanh_sq_prime,
    standardize = TRUE,
    col_scale = "sd",
    eps = 1e-8,
    delta_lead = 1e-6,
    ridge = 1e-8,
    scs_control = list(max_iters = 100000, eps = 1e-5, alpha = 1.8, verbose = FALSE)
  )
  
  # R Wiederholungen für L1
  for (r in 1:R) {
    x <- r_true(n)
    
    fit <- fit_score_matching_matrixG_moments_fast(
      x, m = m,
      h = h_tanh_sq,
      h_prime = h_tanh_sq_prime,
      standardize = TRUE,
      col_scale = "sd",
      eps = 1e-8,
      delta_lead = 1e-6,
      ridge = 1e-8
    )
    
    p_hat <- density_from_fit_fast_trunc(xg_eval, fit, target_drop = 45, n_int = 2500)
    errs[i, r] <- L1_error(xg_eval, p_hat, p_ref_eval)
  }
  
  cat("n=", n,
      " status(rep)=", fit_store[[i]]$status,
      " median L1=", median(errs[i, ], na.rm = TRUE), "\n")
}

# Plot 1: L1 convergence
plot(ns, apply(errs, 1, median, na.rm = TRUE), log = "x", type = "b",
     pch = 19, col = "blue",
     xlab = "Sample size n",
     ylab = "Median L1 error",
     main = paste0("L1-Konvergenz (Skew-Normal truth, alpha=", alpha, "), m=", m))
grid()

# Plot 2: density panels
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
for (i in seq_along(ns)) {
  n   <- ns[i]
  x   <- x_store[[i]]
  fit <- fit_store[[i]]
  
  xgrid <- seq(mu0 - 6 * sd0, mu0 + 6 * sd0, length.out = 600)
  
  p_hat <- density_from_fit_fast_trunc(xgrid, fit, target_drop = 45, n_int = 2500)
  p_ref <- d_true(xgrid)
  
  hist(x, breaks = 30, freq = FALSE,
       col = "grey90", border = "grey70",
       main = paste0("n = ", n, "\n", fit$status),
       xlab = "", ylab = "")
  lines(xgrid, p_hat, lwd = 2, col = "darkgreen")
  lines(xgrid, p_ref, lwd = 2, col = "red", lty = 2)
}
par(mfrow = c(1, 1))






# ============================================================
# (8) DIAGNOSTIC TEST: Projection vs Numerik
#     - L1 vs m at large n
#     - "Pseudo-truth" proxy at n_big
#     - sensitivity to normalization + solver accuracy
# ============================================================

library(CVXR)
library(pracma)
library(sn)

# ---- helper: symmetrize + min eigenvalue (PSD diagnostic) ----
min_eig_sym <- function(G) {
  Gs <- 0.5 * (G + t(G))
  # use onlyvalues for speed
  ev <- eigen(Gs, symmetric = TRUE, only.values = TRUE)$values
  min(ev)
}

# ---- helper: run one fit + density + diagnostics ----
run_one_fit_eval <- function(x, m, xg_eval, p_ref_eval,
                             target_drop = 45, n_int = 2500,
                             ridge = 1e-8, delta_lead = 1e-6,
                             scs_control = list(max_iters = 100000, eps = 1e-5, alpha = 1.8, verbose = FALSE)) {
  fit <- fit_score_matching_matrixG_moments_fast(
    x, m = m,
    h = h_tanh_sq,
    h_prime = h_tanh_sq_prime,
    standardize = TRUE,
    col_scale = "sd",
    eps = 1e-8,
    delta_lead = delta_lead,
    ridge = ridge,
    scs_control = scs_control
  )
  
  p_hat <- density_from_fit_fast_trunc(xg_eval, fit, target_drop = target_drop, n_int = n_int)
  l1 <- L1_error(xg_eval, p_hat, p_ref_eval)
  
  list(
    fit = fit,
    p_hat = p_hat,
    l1 = l1,
    status = fit$status,
    min_eig = min_eig_sym(fit$G)
  )
}

# ============================================================
# (8A) Setup: Skew-Normal truth + stable eval grid
# ============================================================
set.seed(1)

xi    <- 0
omega <- 1
alpha <- 5

r_true <- function(n) rsn(n, xi = xi, omega = omega, alpha = alpha)
d_true <- function(x) dsn(x, xi = xi, omega = omega, alpha = alpha)

delta <- alpha / sqrt(1 + alpha^2)
mu0   <- xi + omega * delta * sqrt(2 / pi)
sd0   <- omega * sqrt(1 - (2 * delta^2) / pi)

xg_eval    <- seq(mu0 - 10 * sd0, mu0 + 10 * sd0, length.out = 2001)
p_ref_eval <- d_true(xg_eval)

# ============================================================
# (8B) Test 1: L1 vs m at large n  (projection-bias check)
#     If this does NOT decrease (on average), it's either
#     (i) true projection doesn't help in L1 for this criterion/class
#     or (ii) solver/numerics are dominating.
# ============================================================
ns_big <- 100000          # big n for stable estimation
R_m    <- 10              # repetitions per m (keep moderate)
m_grid <- 2:8

# stricter solver for this diagnostic (reduce "optimal_inaccurate")
scs_strict <- list(max_iters = 300000, eps = 2e-6, alpha = 1.8, verbose = FALSE)

res_m <- matrix(NA_real_, nrow = length(m_grid), ncol = R_m)
stat_m <- matrix(NA_integer_, nrow = length(m_grid), ncol = R_m) # 1=optimal, 0=otherwise
mineig_m <- matrix(NA_real_, nrow = length(m_grid), ncol = R_m)

for (im in seq_along(m_grid)) {
  m <- m_grid[im]
  for (r in 1:R_m) {
    x <- r_true(ns_big)
    
    out <- run_one_fit_eval(
      x = x, m = m,
      xg_eval = xg_eval, p_ref_eval = p_ref_eval,
      target_drop = 60, n_int = 8000,           # robust normalization
      ridge = 1e-8, delta_lead = 1e-6,
      scs_control = scs_strict
    )
    
    res_m[im, r] <- out$l1
    stat_m[im, r] <- as.integer(identical(out$status, "optimal"))
    mineig_m[im, r] <- out$min_eig
    
    cat("m=", m, " rep=", r,
        " status=", out$status,
        " minEig=", signif(out$min_eig, 3),
        " L1=", signif(out$l1, 4), "\n")
  }
}

medL1_m <- apply(res_m, 1, median, na.rm = TRUE)
optRate_m <- rowMeans(stat_m, na.rm = TRUE)
medMinEig_m <- apply(mineig_m, 1, median, na.rm = TRUE)

plot(m_grid, medL1_m, type = "b", pch = 19, col = "blue",
     xlab = "m", ylab = "Median L1 (n=100000)",
     main = "Diagnostic: L1 vs m (Skew-Normal truth)")
grid()
mtext(paste0("Median minEig(G): ",
             paste(sprintf("m=%d:%.2g", m_grid, medMinEig_m), collapse = "  ")),
      side = 3, cex = 0.7, line = 0.2)
mtext(paste0("Optimal-rate: ",
             paste(sprintf("m=%d:%.2f", m_grid, optRate_m), collapse = "  ")),
      side = 3, cex = 0.7, line = 1.1)

# ============================================================
# (8C) Test 2: "Pseudo-truth" plateau for a fixed m
#     Fit once with HUGE n as proxy for p_m^*, then compare:
#     L1(p_hat_n, p_ref) should approach L1(p_m^*, p_ref).
# ============================================================
m_fix  <- 5
n_star <- 1000000    # big proxy sample (adjust if too slow)
n_list <- c(200, 500, 1000, 2000, 5000, 10000, 20000, 100000)

cat("\n--- Building pseudo-truth proxy p_m^* using n_star=", n_star, " m=", m_fix, " ---\n")
x_star <- r_true(n_star)

out_star <- run_one_fit_eval(
  x = x_star, m = m_fix,
  xg_eval = xg_eval, p_ref_eval = p_ref_eval,
  target_drop = 70, n_int = 12000,
  ridge = 1e-8, delta_lead = 1e-6,
  scs_control = scs_strict
)

p_star <- out_star$p_hat
L1_star <- L1_error(xg_eval, p_star, p_ref_eval)

cat("Pseudo-truth proxy: status=", out_star$status,
    " minEig=", signif(out_star$min_eig, 3),
    " L1(p_star, true)=", signif(L1_star, 5), "\n")

# Now check convergence to plateau across n
R_n <- 20
errs_to_true <- matrix(NA_real_, nrow = length(n_list), ncol = R_n)
errs_to_star <- matrix(NA_real_, nrow = length(n_list), ncol = R_n)
status_inacc <- matrix(FALSE, nrow = length(n_list), ncol = R_n)

for (i in seq_along(n_list)) {
  n <- n_list[i]
  for (r in 1:R_n) {
    x <- r_true(n)
    
    out <- run_one_fit_eval(
      x = x, m = m_fix,
      xg_eval = xg_eval, p_ref_eval = p_ref_eval,
      target_drop = 60, n_int = 8000,
      ridge = 1e-8, delta_lead = 1e-6,
      scs_control = scs_strict
    )
    
    errs_to_true[i, r] <- out$l1
    errs_to_star[i, r] <- L1_error(xg_eval, out$p_hat, p_star)
    status_inacc[i, r] <- !identical(out$status, "optimal")
    
    cat("n=", n, " rep=", r,
        " status=", out$status,
        " L1(true)=", signif(errs_to_true[i, r], 4),
        " L1(star)=", signif(errs_to_star[i, r], 4), "\n")
  }
}

med_true <- apply(errs_to_true, 1, median, na.rm = TRUE)
med_star <- apply(errs_to_star, 1, median, na.rm = TRUE)
inacc_rate <- rowMeans(status_inacc)

plot(n_list, med_true, log = "x", type = "b", pch = 19, col = "blue",
     xlab = "n", ylab = "Median L1",
     main = paste0("Plateau diagnostic (m=", m_fix, "): L1 to true vs to p*"))
lines(n_list, med_star, type = "b", pch = 19, col = "darkgreen")
abline(h = L1_star, col = "red", lty = 2, lwd = 2)
legend("topright",
       legend = c("Median L1(p_hat, true)", "Median L1(p_hat, p_star)", "L1(p_star, true)"),
       col = c("blue", "darkgreen", "red"),
       lty = c(1, 1, 2), lwd = c(2, 2, 2), pch = c(19, 19, NA), bty = "n")
grid()
mtext(paste0("inaccurate-rate by n: ",
             paste(sprintf("n=%g:%.2f", n_list, inacc_rate), collapse = "  ")),
      side = 3, cex = 0.7, line = 0.2)

# Interpretation:
# - If Median L1(p_hat, p_star) -> 0 but Median L1(p_hat, true) -> L1_star > 0:
#   => classic projection plateau (model mismatch).
# - If L1(p_hat, p_star) does NOT go to 0 or is noisy:
#   => solver/numerics/truncation issues still present.


# ============================================================
# (9) DEBUG TEST: Is unscaling breaking PSD?
#     Run same experiment with col_scale="none"
# ============================================================

set.seed(1)

xi    <- 0
omega <- 1
alpha <- 5

r_true <- function(n) rsn(n, xi = xi, omega = omega, alpha = alpha)
d_true <- function(x) dsn(x, xi = xi, omega = omega, alpha = alpha)

delta <- alpha / sqrt(1 + alpha^2)
mu0   <- xi + omega * delta * sqrt(2 / pi)
sd0   <- omega * sqrt(1 - (2 * delta^2) / pi)

xg_eval    <- seq(mu0 - 10 * sd0, mu0 + 10 * sd0, length.out = 2001)
p_ref_eval <- d_true(xg_eval)

min_eig_sym <- function(G) {
  Gs <- 0.5 * (G + t(G))
  eigen(Gs, symmetric = TRUE, only.values = TRUE)$values |> min()
}

n_big <- 100000
R <- 5
m <- 5

scs_strict <- list(max_iters = 300000, eps = 2e-6, alpha = 1.8, verbose = FALSE)

for (r in 1:R) {
  x <- r_true(n_big)
  
  fit <- fit_score_matching_matrixG_moments_fast(
    x, m = m,
    h = h_tanh_sq,
    h_prime = h_tanh_sq_prime,
    standardize = TRUE,
    col_scale = "none",     # <<< IMPORTANT: no scaling/unscaling
    eps = 1e-8,
    delta_lead = 1e-6,
    ridge = 1e-8,
    scs_control = scs_strict
  )
  
  p_hat <- density_from_fit_fast_trunc(xg_eval, fit, target_drop = 60, n_int = 8000)
  l1 <- L1_error(xg_eval, p_hat, p_ref_eval)
  
  cat("rep=", r,
      " status=", fit$status,
      " minEig=", signif(min_eig_sym(fit$G), 4),
      " L1=", signif(l1, 5), "\n")
}





# ============================================================
# (10) PSD-PROJECTION PATCH (postprocess G)
#      fixes negative eigenvalues after unscaling
# ============================================================
project_to_psd <- function(G, eps_floor = 0) {
  Gs <- 0.5 * (G + t(G))
  eig <- eigen(Gs, symmetric = TRUE)
  vals <- eig$values
  vals[vals < eps_floor] <- eps_floor
  Gpsd <- eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
  0.5 * (Gpsd + t(Gpsd))
}

min_eig_sym <- function(G) {
  Gs <- 0.5 * (G + t(G))
  eigen(Gs, symmetric = TRUE, only.values = TRUE)$values |> min()
}


# ============================================================
# (11) DIAGNOSTIC TEST: L1 vs m, raw vs PSD-fixed
#      robust to solver errors, keeps col_scale="sd"
# ============================================================

set.seed(1)

# --- Skew-Normal truth setup (wie bei dir) ---
xi    <- 0
omega <- 1
alpha <- 5

r_true <- function(n) rsn(n, xi = xi, omega = omega, alpha = alpha)
d_true <- function(x) dsn(x, xi = xi, omega = omega, alpha = alpha)

delta <- alpha / sqrt(1 + alpha^2)
mu0   <- xi + omega * delta * sqrt(2 / pi)
sd0   <- omega * sqrt(1 - (2 * delta^2) / pi)

xg_eval    <- seq(mu0 - 10 * sd0, mu0 + 10 * sd0, length.out = 2001)
p_ref_eval <- d_true(xg_eval)

# --- density helper that uses given (G,c1,scaling,z_hint) ---
density_from_components_trunc <- function(xgrid, G, c1, scaling, z_hint,
                                          target_drop = 60, n_int = 8000) {
  zgrid <- (xgrid - scaling$mu) / scaling$s
  tr <- find_truncation_fast(G, c1, z_hint = z_hint, target_drop = target_drop)
  zL <- tr$left; zR <- tr$right; s0 <- tr$s_mode
  
  z_int <- seq(zL, zR, length.out = n_int)
  s_int <- svals_fast(z_int, G, c1)
  un_int <- exp(-(s_int - s0))
  Z <- trapz(z_int, un_int)
  
  s_g <- svals_fast(zgrid, G, c1)
  pz <- exp(-(s_g - s0)) / Z
  (1 / scaling$s) * pz
}

# --- run config ---
n_big <- 100000
R_m   <- 10
m_grid <- 2:8

scs_strict <- list(max_iters = 300000, eps = 2e-6, alpha = 1.8, verbose = FALSE)

# store results
L1_raw   <- matrix(NA_real_, nrow = length(m_grid), ncol = R_m)
L1_psd   <- matrix(NA_real_, nrow = length(m_grid), ncol = R_m)
minEig_raw <- matrix(NA_real_, nrow = length(m_grid), ncol = R_m)
minEig_psd <- matrix(NA_real_, nrow = length(m_grid), ncol = R_m)
status_mat <- matrix("", nrow = length(m_grid), ncol = R_m)

for (im in seq_along(m_grid)) {
  m <- m_grid[im]
  for (r in 1:R_m) {
    x <- r_true(n_big)
    
    fit <- tryCatch(
      fit_score_matching_matrixG_moments_fast(
        x, m = m,
        h = h_tanh_sq,
        h_prime = h_tanh_sq_prime,
        standardize = TRUE,
        col_scale = "sd",       # keep your scaling
        eps = 1e-8,
        delta_lead = 1e-6,
        ridge = 1e-8,
        scs_control = scs_strict
      ),
      error = function(e) {
        cat("m=", m, " rep=", r, " FIT ERROR: ", conditionMessage(e), "\n")
        return(NULL)
      }
    )
    if (is.null(fit)) next
    
    status_mat[im, r] <- fit$status
    
    # --- raw ---
    minEig_raw[im, r] <- min_eig_sym(fit$G)
    p_raw <- density_from_fit_fast_trunc(xg_eval, fit, target_drop = 60, n_int = 8000)
    L1_raw[im, r] <- L1_error(xg_eval, p_raw, p_ref_eval)
    
    # --- PSD-fixed ---
    G_fix <- project_to_psd(fit$G, eps_floor = 0)
    minEig_psd[im, r] <- min_eig_sym(G_fix)
    p_fix <- density_from_components_trunc(
      xgrid = xg_eval,
      G = G_fix, c1 = fit$c1,
      scaling = fit$scaling, z_hint = fit$z_hint,
      target_drop = 60, n_int = 8000
    )
    L1_psd[im, r] <- L1_error(xg_eval, p_fix, p_ref_eval)
    
    cat("m=", m, " rep=", r,
        " status=", fit$status,
        " minEig(raw)=", signif(minEig_raw[im, r], 3),
        " L1(raw)=", signif(L1_raw[im, r], 4),
        " minEig(psd)=", signif(minEig_psd[im, r], 3),
        " L1(psd)=", signif(L1_psd[im, r], 4), "\n")
  }
}

# summarize + plot
med_raw <- apply(L1_raw, 1, median, na.rm = TRUE)
med_psd <- apply(L1_psd, 1, median, na.rm = TRUE)

plot(m_grid, med_raw, type = "b", pch = 19, col = "blue",
     xlab = "m", ylab = "Median L1 (n=100000)",
     main = "Diagnostic: raw vs PSD-fixed (Skew-Normal truth)")
lines(m_grid, med_psd, type = "b", pch = 19, col = "darkgreen")
legend("topright",
       legend = c("raw (your output)", "PSD-fixed postprocess"),
       col = c("blue", "darkgreen"),
       lty = 1, lwd = 2, pch = 19, bty = "n")
grid()

cat("\nMedian minEig(raw) by m:\n")
print(setNames(apply(minEig_raw, 1, median, na.rm = TRUE), m_grid))
cat("\nMedian minEig(psd) by m:\n")
print(setNames(apply(minEig_psd, 1, median, na.rm = TRUE), m_grid))





# ---------- Alternative h's ----------
h_const1 <- function(z) rep(1, length(z))
h_const1_prime <- function(z) rep(0, length(z))

# langsam wachsend: log(1+z^2)
h_log1pz2 <- function(z) log1p(z^2)
h_log1pz2_prime <- function(z) 2*z/(1 + z^2)

run_h_weight_test <- function(
    r_true, d_true,
    ns = c(200, 500, 1000, 2000, 5000, 10000, 20000, 100000),
    R = 20,
    m = 5,
    xg_eval = NULL,
    target_drop = 60,
    n_int = 8000,
    scs_control = list(max_iters = 300000, eps = 2e-6, alpha = 1.8, verbose = FALSE),
    seed = 1
){
  set.seed(seed)
  
  # stabile eval grid falls nicht gegeben
  if (is.null(xg_eval)) {
    x_hint <- r_true(50000)
    mu0 <- mean(x_hint); sd0 <- sd(x_hint)
    xg_eval <- seq(mu0 - 10*sd0, mu0 + 10*sd0, length.out = 2001)
  }
  p_ref_eval <- d_true(xg_eval)
  
  h_list <- list(
    tanh_sq = list(h = h_tanh_sq, hp = h_tanh_sq_prime),
    const1  = list(h = h_const1,  hp = h_const1_prime),
    log1pz2 = list(h = h_log1pz2, hp = h_log1pz2_prime)
  )
  
  errs <- array(NA_real_, dim = c(length(ns), R, length(h_list)),
                dimnames = list(paste0("n=", ns), paste0("r=", 1:R), names(h_list)))
  status <- array("", dim = c(length(ns), R, length(h_list)),
                  dimnames = dimnames(errs))
  
  for (i in seq_along(ns)) {
    n <- ns[i]
    cat("\n=== n =", n, "===\n")
    
    for (r in 1:R) {
      x <- r_true(n)
      
      for (hh in names(h_list)) {
        fit <- fit_score_matching_matrixG_moments_fast(
          x, m = m,
          h = h_list[[hh]]$h,
          h_prime = h_list[[hh]]$hp,
          standardize = TRUE,
          col_scale = "sd",
          eps = 1e-8,
          delta_lead = 1e-6,
          ridge = 1e-8,
          scs_control = scs_control
        )
        
        p_hat <- density_from_fit_fast_trunc(
          xg_eval, fit,
          target_drop = target_drop,
          n_int = n_int
        )
        
        errs[i, r, hh] <- L1_error(xg_eval, p_hat, p_ref_eval)
        status[i, r, hh] <- fit$status
      }
    }
    
    cat("Median L1 by h:\n")
    print(apply(errs[i, , , drop = FALSE], 3, median, na.rm = TRUE))
  }
  
  # Plot: median L1 vs n für jedes h
  med <- apply(errs, c(1,3), median, na.rm = TRUE)  # [n, h]
  matplot(ns, med, log = "x", type = "b", pch = 19,
          xlab = "n", ylab = "Median L1",
          main = paste0("h-Weight Test (m=", m, ")"))
  legend("topright", legend = colnames(med), col = 1:ncol(med), lty = 1:ncol(med), pch = 19, bty = "n")
  grid()
  
  invisible(list(errs = errs, med = med, ns = ns, xg_eval = xg_eval))
}


# --- Beispiel: Skew-Normal Wahrheit (wie bei dir) ---
xi <- 0; omega <- 1; alpha <- 5
r_true <- function(n) sn::rsn(n, xi = xi, omega = omega, alpha = alpha)
d_true <- function(x) sn::dsn(x, xi = xi, omega = omega, alpha = alpha)

# stabiles Grid anhand theoretischer Momente (optional)
delta <- alpha / sqrt(1 + alpha^2)
mu0 <- xi + omega * delta * sqrt(2/pi)
sd0 <- omega * sqrt(1 - (2*delta^2)/pi)
xg_eval <- seq(mu0 - 10*sd0, mu0 + 10*sd0, length.out = 2001)

res_h <- run_h_weight_test(
  r_true = r_true, d_true = d_true,
  ns = c(500, 1000, 2000, 5000, 10000, 20000, 100000),
  R = 20,
  m = 5,
  xg_eval = xg_eval,
  target_drop = 60, n_int = 8000
)

# m(n) und Normalisierung an n koppeln
m_of_n <- function(n) {
  if (n <= 1000) 2
  else if (n <= 5000) 3
  else if (n <= 20000) 4
  else if (n <= 100000) 5
  else 6
}

# stärkere Normalisierung bei größerem m/n
norm_of_mn <- function(m, n) {
  list(
    target_drop = 45 + 5 * (m - 2),           # z.B. 45,50,55,...
    n_int = as.integer(2500 + 800 * (m - 2))  # z.B. 2500,3300,4100,...
  )
}

run_sieve_test <- function(
    r_true, d_true,
    ns = c(200, 500, 1000, 2000, 5000, 10000, 20000, 100000),
    R = 20,
    m_fixed = 5,
    m_of_n,
    norm_of_mn,
    xg_eval = NULL,
    scs_control = list(max_iters = 300000, eps = 2e-6, alpha = 1.8, verbose = FALSE),
    seed = 1
){
  set.seed(seed)
  
  if (is.null(xg_eval)) {
    x_hint <- r_true(50000)
    mu0 <- mean(x_hint); sd0 <- sd(x_hint)
    xg_eval <- seq(mu0 - 10*sd0, mu0 + 10*sd0, length.out = 2001)
  }
  p_ref_eval <- d_true(xg_eval)
  
  errs_fixed <- matrix(NA_real_, nrow = length(ns), ncol = R,
                       dimnames = list(paste0("n=", ns), paste0("r=", 1:R)))
  errs_sieve <- matrix(NA_real_, nrow = length(ns), ncol = R,
                       dimnames = dimnames(errs_fixed))
  
  for (i in seq_along(ns)) {
    n <- ns[i]
    m_sieve <- m_of_n(n)
    
    # Normalisierung für fixed m und sieve m separat (hier: abhängig von m)
    norm_fixed <- norm_of_mn(m_fixed, n)
    norm_sieve <- norm_of_mn(m_sieve, n)
    
    cat("\n=== n =", n, " | m_fixed =", m_fixed, " | m_sieve =", m_sieve, "===\n")
    cat("norm_fixed:", norm_fixed$target_drop, norm_fixed$n_int, "\n")
    cat("norm_sieve:", norm_sieve$target_drop, norm_sieve$n_int, "\n")
    
    for (r in 1:R) {
      x <- r_true(n)
      
      # --- fixed m ---
      fit_f <- fit_score_matching_matrixG_moments_fast(
        x, m = m_fixed,
        h = h_tanh_sq,
        h_prime = h_tanh_sq_prime,
        standardize = TRUE,
        col_scale = "sd",
        eps = 1e-8,
        delta_lead = 1e-6,
        ridge = 1e-8,
        scs_control = scs_control
      )
      p_f <- density_from_fit_fast_trunc(xg_eval, fit_f,
                                         target_drop = norm_fixed$target_drop,
                                         n_int = norm_fixed$n_int)
      errs_fixed[i, r] <- L1_error(xg_eval, p_f, p_ref_eval)
      
      # --- sieve m(n) ---
      fit_s <- fit_score_matching_matrixG_moments_fast(
        x, m = m_sieve,
        h = h_tanh_sq,
        h_prime = h_tanh_sq_prime,
        standardize = TRUE,
        col_scale = "sd",
        eps = 1e-8,
        delta_lead = 1e-6,
        ridge = 1e-8,
        scs_control = scs_control
      )
      p_s <- density_from_fit_fast_trunc(xg_eval, fit_s,
                                         target_drop = norm_sieve$target_drop,
                                         n_int = norm_sieve$n_int)
      errs_sieve[i, r] <- L1_error(xg_eval, p_s, p_ref_eval)
    }
    
    cat("Median L1 fixed:", median(errs_fixed[i, ], na.rm = TRUE),
        " | sieve:", median(errs_sieve[i, ], na.rm = TRUE), "\n")
  }
  
  med_fixed <- apply(errs_fixed, 1, median, na.rm = TRUE)
  med_sieve <- apply(errs_sieve, 1, median, na.rm = TRUE)
  
  plot(ns, med_fixed, log = "x", type = "b", pch = 19, col = "blue",
       xlab = "n", ylab = "Median L1",
       main = paste0("Sieve test: fixed m=", m_fixed, " vs m(n)"))
  lines(ns, med_sieve, type = "b", pch = 19, col = "darkgreen")
  legend("topright", legend = c("fixed m", "sieve m(n)"),
         col = c("blue", "darkgreen"), lty = 1, lwd = 2, pch = 19, bty = "n")
  grid()
  
  invisible(list(errs_fixed = errs_fixed, errs_sieve = errs_sieve,
                 med_fixed = med_fixed, med_sieve = med_sieve, ns = ns))
}

res_sieve <- run_sieve_test(
  r_true = r_true, d_true = d_true,
  ns = c(500, 1000, 2000, 5000, 10000, 20000, 100000),
  R = 20,
  m_fixed = 5,
  m_of_n = m_of_n,
  norm_of_mn = norm_of_mn,
  xg_eval = xg_eval
)




sieve_to_zero_report <- function(
    r_true, d_true,
    # n-Werte: so wählen, dass es schnell bleibt, aber Trend sichtbar
    ns = c(2000, 5000, 20000, 100000),
    # Sieve schedule: muss m wachsen lassen (sonst testest du nix)
    m_of_n = function(n) {
      if (n <= 5000) 3 else if (n <= 20000) 4 else if (n <= 100000) 5 else 6
    },
    # wir brauchen p*_m für m bis max(m_of_n(ns)) (und optional +1)
    m_grid = NULL,
    # "Pseudo-Truth" sample size (so klein wie möglich, so groß wie nötig)
    n_star = 300000,
    # repetitions pro n (klein -> schnell)
    R = 3,
    # schnelle, aber brauchbare Normalisierung
    target_drop = 55,
    n_int = 2000,
    # L1-grid (klein -> schnell)
    grid_len = 801,
    grid_sd_mult = 10,
    # solver settings (moderat, sonst dauert's ewig)
    scs_control = list(max_iters = 160000, eps = 1e-5, alpha = 1.8, verbose = FALSE),
    seed = 1
){
  stopifnot(is.function(r_true), is.function(d_true))
  set.seed(seed)
  
  # determine m grid
  if (is.null(m_grid)) {
    mmax <- max(sapply(ns, m_of_n))
    m_grid <- sort(unique(2:mmax)) # start at 2; m=1 is almost gaussian-only
  }
  
  # stable evaluation grid from a hint sample
  x_hint <- r_true(50000)
  mu0 <- mean(x_hint); sd0 <- sd(x_hint)
  xg <- seq(mu0 - grid_sd_mult*sd0, mu0 + grid_sd_mult*sd0, length.out = grid_len)
  p_true <- d_true(xg)
  
  # helper: fit -> density
  fit_and_density <- function(x, m) {
    fit <- fit_score_matching_matrixG_moments_fast(
      x, m = m,
      h = h_tanh_sq,
      h_prime = h_tanh_sq_prime,
      standardize = TRUE,
      col_scale = "sd",
      eps = 1e-8,
      delta_lead = 1e-6,
      ridge = 1e-8,
      scs_control = scs_control
    )
    p_hat <- density_from_fit_fast_trunc(
      xg, fit, target_drop = target_drop, n_int = n_int
    )
    list(fit = fit, p = p_hat)
  }
  
  # ------------------------------------------------------------
  # (1) Build pseudo-truths p*_m using one big sample x_star
  # ------------------------------------------------------------
  x_star <- r_true(n_star)
  
  p_star_list <- vector("list", length(m_grid))
  names(p_star_list) <- paste0("m=", m_grid)
  
  L1_star_to_true <- numeric(length(m_grid))
  status_star <- character(length(m_grid))
  
  for (i in seq_along(m_grid)) {
    m <- m_grid[i]
    out <- fit_and_density(x_star, m)
    p_star_list[[i]] <- out$p
    L1_star_to_true[i] <- L1_error(xg, out$p, p_true)
    status_star[i] <- out$fit$status
  }
  
  # ------------------------------------------------------------
  # (2) Main sieve loop: compare p_hat(n,m(n)) to p*_m(n) and to true
  # ------------------------------------------------------------
  res <- data.frame(
    n = ns,
    m_sieve = sapply(ns, m_of_n),
    med_L1_hat_to_true = NA_real_,
    med_L1_hat_to_star = NA_real_,
    med_status_optimal = NA_real_,  # fraction optimal
    stringsAsFactors = FALSE
  )
  
  # store raw vectors for report
  raw_lines <- character(0)
  
  for (j in seq_along(ns)) {
    n <- ns[j]
    m <- m_of_n(n)
    key <- paste0("m=", m)
    if (!(key %in% names(p_star_list))) {
      stop("m_of_n(ns) produced m=", m, " but m_grid does not include it. Add it to m_grid.")
    }
    p_star_m <- p_star_list[[key]]
    
    L1_hat_true <- numeric(R)
    L1_hat_star <- numeric(R)
    st_opt <- logical(R)
    
    for (r in 1:R) {
      x <- r_true(n)
      out <- fit_and_density(x, m)
      p_hat <- out$p
      
      L1_hat_true[r] <- L1_error(xg, p_hat, p_true)
      L1_hat_star[r] <- L1_error(xg, p_hat, p_star_m)
      st_opt[r] <- identical(out$fit$status, "optimal")
      
      raw_lines <- c(raw_lines,
                     sprintf("n=%g rep=%d m=%d status=%s L1(hat,true)=%.6f L1(hat,star_m)=%.6f",
                             n, r, m, out$fit$status, L1_hat_true[r], L1_hat_star[r])
      )
    }
    
    res$med_L1_hat_to_true[j] <- median(L1_hat_true, na.rm = TRUE)
    res$med_L1_hat_to_star[j] <- median(L1_hat_star, na.rm = TRUE)
    res$med_status_optimal[j] <- mean(st_opt, na.rm = TRUE)
  }
  
  # ------------------------------------------------------------
  # (3) Decide if sieve plausibly -> 0 (heuristics)
  # ------------------------------------------------------------
  # Condition A: estimation error to star should decrease with n (toward 0)
  # Condition B: approximation error (star to true) should decrease with m
  approx_tbl <- data.frame(
    m = m_grid,
    status_star = status_star,
    L1_star_to_true = L1_star_to_true
  )
  
  # monotonic-ish check
  approx_decreasing <- all(diff(approx_tbl$L1_star_to_true) <= 1e-3) # allow small noise
  est_decreasing <- all(diff(res$med_L1_hat_to_star) <= 1e-3)
  
  # ------------------------------------------------------------
  # (4) Print a copy-paste report
  # ------------------------------------------------------------
  cat("\n===SIEVE_TO_ZERO_REPORT===\n")
  cat("seed=", seed, "\n", sep = "")
  cat("n_star=", n_star, "\n", sep = "")
  cat("ns=", paste(ns, collapse = ","), "\n", sep = "")
  cat("m_sieve(ns)=", paste(res$m_sieve, collapse = ","), "\n", sep = "")
  cat("m_grid_for_star=", paste(m_grid, collapse = ","), "\n", sep = "")
  cat("grid_len=", grid_len, " grid_sd_mult=", grid_sd_mult, "\n", sep = "")
  cat("target_drop=", target_drop, " n_int=", n_int, "\n", sep = "")
  cat("solver(max_iters,eps,alpha)=", scs_control$max_iters, ",", scs_control$eps, ",", scs_control$alpha, "\n", sep = "")
  
  cat("\n--- PSEUDO-TRUTH (p*_m) QUALITY: L1(p*_m, true) ---\n")
  for (i in seq_len(nrow(approx_tbl))) {
    cat(sprintf("m=%d status=%s L1(star_m,true)=%.6f\n",
                approx_tbl$m[i], approx_tbl$status_star[i], approx_tbl$L1_star_to_true[i]))
  }
  cat("approx_decreasing(m): ", approx_decreasing, "\n", sep = "")
  
  cat("\n--- SIEVE RESULTS (medians over reps) ---\n")
  for (j in seq_len(nrow(res))) {
    cat(sprintf("n=%g m=%d medL1(hat,true)=%.6f medL1(hat,star_m)=%.6f optimal_rate=%.2f\n",
                res$n[j], res$m_sieve[j], res$med_L1_hat_to_true[j],
                res$med_L1_hat_to_star[j], res$med_status_optimal[j]))
  }
  cat("est_decreasing(n): ", est_decreasing, "\n", sep = "")
  
  cat("\n--- RAW REPS ---\n")
  cat(paste(raw_lines, collapse = "\n"), "\n")
  
  cat("\n--- INTERPRETATION RULE ---\n")
  cat("If medL1(hat,star_m) decreases toward ~0 as n increases AND L1(star_m,true) decreases toward ~0 as m increases,\n")
  cat("then sieve -> 0 is supported (estimation + approximation).\n")
  cat("If only medL1(hat,star_m)->0 but L1(star_m,true) plateaus >0, then you have projection bias (model/metric mismatch).\n")
  cat("If medL1(hat,star_m) does NOT go down, solver/numerics dominate (not a sieve issue).\n")
  
  cat("===END_SIEVE_TO_ZERO_REPORT===\n")
  
  invisible(list(res = res, approx_tbl = approx_tbl))
}

sieve_to_zero_report <- function(
    r_true, d_true,
    # n-Werte: so wählen, dass es schnell bleibt, aber Trend sichtbar
    ns = c(2000, 5000, 20000, 100000),
    # Sieve schedule: muss m wachsen lassen (sonst testest du nix)
    m_of_n = function(n) {
      if (n <= 5000) 3 else if (n <= 20000) 4 else if (n <= 100000) 5 else 6
    },
    # wir brauchen p*_m für m bis max(m_of_n(ns)) (und optional +1)
    m_grid = NULL,
    # "Pseudo-Truth" sample size (so klein wie möglich, so groß wie nötig)
    n_star = 300000,
    # repetitions pro n (klein -> schnell)
    R = 3,
    # schnelle, aber brauchbare Normalisierung
    target_drop = 55,
    n_int = 2000,
    # L1-grid (klein -> schnell)
    grid_len = 801,
    grid_sd_mult = 10,
    # solver settings (moderat, sonst dauert's ewig)
    scs_control = list(max_iters = 160000, eps = 1e-5, alpha = 1.8, verbose = FALSE),
    seed = 1
){
  stopifnot(is.function(r_true), is.function(d_true))
  set.seed(seed)
  
  # determine m grid
  if (is.null(m_grid)) {
    mmax <- max(sapply(ns, m_of_n))
    m_grid <- sort(unique(2:mmax)) # start at 2; m=1 is almost gaussian-only
  }
  
  # stable evaluation grid from a hint sample
  x_hint <- r_true(50000)
  mu0 <- mean(x_hint); sd0 <- sd(x_hint)
  xg <- seq(mu0 - grid_sd_mult*sd0, mu0 + grid_sd_mult*sd0, length.out = grid_len)
  p_true <- d_true(xg)
  
  # helper: fit -> density
  fit_and_density <- function(x, m) {
    fit <- fit_score_matching_matrixG_moments_fast(
      x, m = m,
      h = h_tanh_sq,
      h_prime = h_tanh_sq_prime,
      standardize = TRUE,
      col_scale = "sd",
      eps = 1e-8,
      delta_lead = 1e-6,
      ridge = 1e-8,
      scs_control = scs_control
    )
    p_hat <- density_from_fit_fast_trunc(
      xg, fit, target_drop = target_drop, n_int = n_int
    )
    list(fit = fit, p = p_hat)
  }
  
  # ------------------------------------------------------------
  # (1) Build pseudo-truths p*_m using one big sample x_star
  # ------------------------------------------------------------
  x_star <- r_true(n_star)
  
  p_star_list <- vector("list", length(m_grid))
  names(p_star_list) <- paste0("m=", m_grid)
  
  L1_star_to_true <- numeric(length(m_grid))
  status_star <- character(length(m_grid))
  
  for (i in seq_along(m_grid)) {
    m <- m_grid[i]
    out <- fit_and_density(x_star, m)
    p_star_list[[i]] <- out$p
    L1_star_to_true[i] <- L1_error(xg, out$p, p_true)
    status_star[i] <- out$fit$status
  }
  
  # ------------------------------------------------------------
  # (2) Main sieve loop: compare p_hat(n,m(n)) to p*_m(n) and to true
  # ------------------------------------------------------------
  res <- data.frame(
    n = ns,
    m_sieve = sapply(ns, m_of_n),
    med_L1_hat_to_true = NA_real_,
    med_L1_hat_to_star = NA_real_,
    med_status_optimal = NA_real_,  # fraction optimal
    stringsAsFactors = FALSE
  )
  
  # store raw vectors for report
  raw_lines <- character(0)
  
  for (j in seq_along(ns)) {
    n <- ns[j]
    m <- m_of_n(n)
    key <- paste0("m=", m)
    if (!(key %in% names(p_star_list))) {
      stop("m_of_n(ns) produced m=", m, " but m_grid does not include it. Add it to m_grid.")
    }
    p_star_m <- p_star_list[[key]]
    
    L1_hat_true <- numeric(R)
    L1_hat_star <- numeric(R)
    st_opt <- logical(R)
    
    for (r in 1:R) {
      x <- r_true(n)
      out <- fit_and_density(x, m)
      p_hat <- out$p
      
      L1_hat_true[r] <- L1_error(xg, p_hat, p_true)
      L1_hat_star[r] <- L1_error(xg, p_hat, p_star_m)
      st_opt[r] <- identical(out$fit$status, "optimal")
      
      raw_lines <- c(raw_lines,
                     sprintf("n=%g rep=%d m=%d status=%s L1(hat,true)=%.6f L1(hat,star_m)=%.6f",
                             n, r, m, out$fit$status, L1_hat_true[r], L1_hat_star[r])
      )
    }
    
    res$med_L1_hat_to_true[j] <- median(L1_hat_true, na.rm = TRUE)
    res$med_L1_hat_to_star[j] <- median(L1_hat_star, na.rm = TRUE)
    res$med_status_optimal[j] <- mean(st_opt, na.rm = TRUE)
  }
  
  # ------------------------------------------------------------
  # (3) Decide if sieve plausibly -> 0 (heuristics)
  # ------------------------------------------------------------
  # Condition A: estimation error to star should decrease with n (toward 0)
  # Condition B: approximation error (star to true) should decrease with m
  approx_tbl <- data.frame(
    m = m_grid,
    status_star = status_star,
    L1_star_to_true = L1_star_to_true
  )
  
  # monotonic-ish check
  approx_decreasing <- all(diff(approx_tbl$L1_star_to_true) <= 1e-3) # allow small noise
  est_decreasing <- all(diff(res$med_L1_hat_to_star) <= 1e-3)
  
  # ------------------------------------------------------------
  # (4) Print a copy-paste report
  # ------------------------------------------------------------
  cat("\n===SIEVE_TO_ZERO_REPORT===\n")
  cat("seed=", seed, "\n", sep = "")
  cat("n_star=", n_star, "\n", sep = "")
  cat("ns=", paste(ns, collapse = ","), "\n", sep = "")
  cat("m_sieve(ns)=", paste(res$m_sieve, collapse = ","), "\n", sep = "")
  cat("m_grid_for_star=", paste(m_grid, collapse = ","), "\n", sep = "")
  cat("grid_len=", grid_len, " grid_sd_mult=", grid_sd_mult, "\n", sep = "")
  cat("target_drop=", target_drop, " n_int=", n_int, "\n", sep = "")
  cat("solver(max_iters,eps,alpha)=", scs_control$max_iters, ",", scs_control$eps, ",", scs_control$alpha, "\n", sep = "")
  
  cat("\n--- PSEUDO-TRUTH (p*_m) QUALITY: L1(p*_m, true) ---\n")
  for (i in seq_len(nrow(approx_tbl))) {
    cat(sprintf("m=%d status=%s L1(star_m,true)=%.6f\n",
                approx_tbl$m[i], approx_tbl$status_star[i], approx_tbl$L1_star_to_true[i]))
  }
  cat("approx_decreasing(m): ", approx_decreasing, "\n", sep = "")
  
  cat("\n--- SIEVE RESULTS (medians over reps) ---\n")
  for (j in seq_len(nrow(res))) {
    cat(sprintf("n=%g m=%d medL1(hat,true)=%.6f medL1(hat,star_m)=%.6f optimal_rate=%.2f\n",
                res$n[j], res$m_sieve[j], res$med_L1_hat_to_true[j],
                res$med_L1_hat_to_star[j], res$med_status_optimal[j]))
  }
  cat("est_decreasing(n): ", est_decreasing, "\n", sep = "")
  
  cat("\n--- RAW REPS ---\n")
  cat(paste(raw_lines, collapse = "\n"), "\n")
  
  cat("\n--- INTERPRETATION RULE ---\n")
  cat("If medL1(hat,star_m) decreases toward ~0 as n increases AND L1(star_m,true) decreases toward ~0 as m increases,\n")
  cat("then sieve -> 0 is supported (estimation + approximation).\n")
  cat("If only medL1(hat,star_m)->0 but L1(star_m,true) plateaus >0, then you have projection bias (model/metric mismatch).\n")
  cat("If medL1(hat,star_m) does NOT go down, solver/numerics dominate (not a sieve issue).\n")
  
  cat("===END_SIEVE_TO_ZERO_REPORT===\n")
  
  invisible(list(res = res, approx_tbl = approx_tbl))
}


xi <- 0; omega <- 1; alpha <- 5
r_true <- function(n) rsn(n, xi = xi, omega = omega, alpha = alpha)
d_true <- function(x) dsn(x, xi = xi, omega = omega, alpha = alpha)

sieve_to_zero_report(
  r_true = r_true,
  d_true = d_true,
  ns = c(2000, 5000, 20000, 100000),
  m_of_n = function(n) { if (n <= 5000) 3 else if (n <= 20000) 4 else 5 },
  m_grid = 2:5,        # muss alle m enthalten, die m_of_n erzeugt
  n_star = 300000,
  R = 3,
  target_drop = 55,
  n_int = 2000
)
