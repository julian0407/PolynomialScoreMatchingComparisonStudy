library(CVXR)
library(pracma)

# ============================================================
# (1) Basics: scale + stable exp
# ============================================================
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
# (2) Fast M/N features like your old make_MN, but vectorized
#     We keep exactly your basis:
#       N_ij(z) = z^(i+j)
#       M_ij(z) = z^(i+j+1)/(i+j+1)
# ============================================================
precompute_exponents <- function(m) {
  expoN <- as.vector(outer(0:(m-1), 0:(m-1), `+`))  # length p
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
  as.numeric(Phi %*% coeff) + c1 * zgrid
}

# ============================================================
# (4) Truncation for normalization: cheap + stable
#     - find mode on a grid
#     - expand left/right until s(z)-s(mode) >= target_drop
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
    eps = 1e-8,
    delta_lead = 1e-6,     # coercivity (fix tails/center)
    ridge = 1e-8,
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
  hx <- pmax(hx, 0)
  
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
