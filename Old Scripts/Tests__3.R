library(CVXR)
library(pracma)

# ============================================================
# Basis: schnelle MN via outer + vorcomputierter Nenner
# ============================================================
precompute_denoms <- function(m) {
  i0 <- 0:(m - 1)
  j0 <- 0:(m - 1)
  I  <- matrix(rep(i0, times = m), nrow = m, byrow = FALSE)
  J  <- matrix(rep(j0, each  = m), nrow = m, byrow = FALSE)
  denomM <- 1 / (I + J + 1)  # for M = x * N * denomM
  list(denomM = denomM)
}

make_MN_fast <- function(x, m, denoms) {
  v <- x^(0:(m - 1))
  N <- tcrossprod(v, v)                 # v v^T
  M <- x * (N * denoms$denomM)          # elementwise denom
  list(M = M, N = N)
}

scale_x <- function(x) {
  mu <- mean(x)
  s  <- sd(x)
  if (!is.finite(s) || s <= 0) s <- 1
  list(z = (x - mu) / s, mu = mu, s = s)
}

# ============================================================
# RICHTIGER "FAST" FIT: aggregiert 1. und 2. Momente
# (kein Gaussian-/m=1-Spezialfall, identisches Objective wie groß)
# ============================================================
fit_score_matching_moments <- function(
    x,
    m,
    h       = function(z) rep(1, length(z)),
    h_prime = function(z) rep(0, length(z)),
    standardize = TRUE,
    col_scale = c("sd", "maxabs", "none"),
    delta_lead = 1e-8,   # coercivity / proper tails
    ridge = 1e-8,
    solver = "SCS",
    scs_control = list(max_iters = 200000, eps = 1e-5, alpha = 1.8, verbose = FALSE),
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
  denoms <- precompute_denoms(m)
  
  hx  <- as.numeric(h(z))
  hpx <- as.numeric(h_prime(z))
  stopifnot(length(hx) == n, length(hpx) == n)
  hx <- pmax(hx, 0)
  
  # Aggregate moments (weighted)
  sum_S  <- matrix(0, p, p)  # sum h * a a^T
  sum_t  <- rep(0, p)        # sum h * a
  sum_u  <- 0                # sum h
  sum_r  <- rep(0, p)        # sum hp * a
  sum_q  <- 0                # sum hp
  sum_b  <- rep(0, p)        # sum h * b_k  where b_k=vec(N)
  
  # also accumulate column scaling stats for a (unweighted is fine; weighted also ok)
  sum_a  <- rep(0, p)
  sum_a2 <- rep(0, p)
  
  for (k in seq_len(n)) {
    MN <- make_MN_fast(z[k], m, denoms)
    a  <- as.vector(MN$M)           # p-vector
    b  <- as.vector(MN$N)           # p-vector (for b_bar)
    
    hk  <- hx[k]
    hpk <- hpx[k]
    
    # weighted moments for objective
    if (hk != 0) {
      sum_u <- sum_u + hk
      sum_t <- sum_t + hk * a
      sum_S <- sum_S + hk * tcrossprod(a, a)
      sum_b <- sum_b + hk * b
    }
    if (hpk != 0) {
      sum_q <- sum_q + hpk
      sum_r <- sum_r + hpk * a
    }
    
    # scaling stats
    sum_a  <- sum_a  + a
    sum_a2 <- sum_a2 + a * a
  }
  
  # empirical expectations
  S <- sum_S / n
  t <- sum_t / n
  u <- sum_u / n
  r <- sum_r / n
  q <- sum_q / n
  bbar <- sum_b / n
  
  # column scaling for a (like your robust version)
  if (col_scale == "none") {
    scale_vec <- rep(1, p)
  } else if (col_scale == "maxabs") {
    # approximate maxabs via mean+sd proxy isn't safe, so just use sd if unsure
    # but to keep option: use sqrt(E[a^2]) as stable proxy
    scale_vec <- sqrt(pmax((sum_a2 / n), 1e-12))
  } else { # "sd"
    mean_a <- sum_a / n
    var_a  <- pmax((sum_a2 / n) - mean_a^2, 1e-12)
    scale_vec <- sqrt(var_a)
  }
  
  # Apply scaling analytically: A_sc = A D^{-1}
  # => S_sc = D^{-1} S D^{-1}, t_sc = D^{-1} t, r_sc = D^{-1} r, b_sc = D^{-1} b
  Dinv <- diag(1 / scale_vec, p, p)
  S_sc <- Dinv %*% S %*% Dinv
  t_sc <- as.vector(Dinv %*% t)
  r_sc <- as.vector(Dinv %*% r)
  b_sc <- as.vector(Dinv %*% bbar)
  
  # Symmetrize for numerical stability
  S_sc <- 0.5 * (S_sc + t(S_sc))
  
  # Build small QP in y = [g; c]
  # Objective:
  # 0.5 * (g,c)^T [S_sc  t_sc; t_sc^T  u] (g,c) - (r_sc + b_sc)^T g - q*c + ridge/2*(||g||^2 + c^2)
  K <- rbind(
    cbind(S_sc, matrix(t_sc, ncol = 1)),
    cbind(matrix(t_sc, nrow = 1), matrix(u, 1, 1))
  )
  K <- 0.5 * (K + t(K))
  
  lin_g <- (r_sc + b_sc)
  lin_c <- q
  lin   <- c(lin_g, lin_c)
  
  # CVXR variables
  G  <- Variable(m, m, PSD = TRUE)
  c1 <- Variable(1)
  gvec <- vec(G)
  y <- vstack(gvec, c1)
  
  K_c   <- Constant(K)
  lin_c <- Constant(matrix(lin, ncol = 1))
  
  obj <- 0.5 * quad_form(y, K_c) - t(lin_c) %*% y
  if (ridge > 0) {
    obj <- obj + (ridge / 2) * sum_squares(y)
  }
  
  constraints <- list(
    G == t(G),
    G[m, m] >= delta_lead
  )
  
  prob <- Problem(Minimize(obj), constraints)
  
  solve_one <- function(eps_override = NULL, it_override = NULL) {
    ctl <- scs_control
    if (!is.null(eps_override)) ctl$eps <- eps_override
    if (!is.null(it_override))  ctl$max_iters <- it_override
    solve(prob, solver = solver,
          max_iters = ctl$max_iters, eps = ctl$eps,
          alpha = ctl$alpha, verbose = ctl$verbose)
  }
  
  sol <- solve_one()
  if (retry_if_inaccurate && identical(sol$status, "optimal_inaccurate") && solver == "SCS") {
    sol <- solve_one(eps_override = 1e-6, it_override = 400000)
  }
  
  # Unscale g back: g_org = D^{-1} g_sc
  g_sc  <- as.numeric(sol$getValue(gvec))
  g_org <- g_sc / scale_vec
  G_org <- matrix(g_org, nrow = m, ncol = m)
  
  list(
    G = G_org,
    c1 = as.numeric(sol$getValue(c1)),
    status = sol$status,
    scaling = list(mu = sc$mu, s = sc$s, standardize = standardize),
    hyper = list(delta_lead = delta_lead, ridge = ridge, col_scale = col_scale),
    moments = list(u = u)
  )
}

# ============================================================
# s(z) + stabile Dichte (wirklich adaptive Trunkierung)
# ============================================================
s_function <- function(z, G, c1) {
  m <- nrow(G)
  val <- 0
  for (i in 0:(m - 1)) {
    for (j in 0:(m - 1)) {
      val <- val +
        G[i + 1, j + 1] *
        z^(i + j + 2) / ((i + j + 1) * (i + j + 2))
    }
  }
  val + c1 * z
}

find_mode_z <- function(G, c1, grow = 2, max_grow = 10) {
  f <- function(z) s_function(z, G, c1)
  lo <- -grow; hi <- grow
  best <- optimize(f, c(lo, hi))
  for (k in 1:max_grow) {
    lo2 <- lo * 2; hi2 <- hi * 2
    best2 <- optimize(f, c(lo2, hi2))
    if (best2$objective < best$objective - 1e-10) {
      lo <- lo2; hi <- hi2; best <- best2
    } else break
  }
  list(z_mode = best$minimum, s_mode = best$objective)
}

find_truncation <- function(G, c1, target = 45, start = 2, max_expand = 20) {
  md <- find_mode_z(G, c1, grow = start)
  z0 <- md$z_mode; s0 <- md$s_mode
  f <- function(z) s_function(z, G, c1)
  
  L <- start; R <- start
  for (k in 1:max_expand) {
    if (f(z0 - L) - s0 < target) L <- L * 2 else break
  }
  for (k in 1:max_expand) {
    if (f(z0 + R) - s0 < target) R <- R * 2 else break
  }
  list(left = z0 - L, right = z0 + R, z_mode = z0, s_mode = s0)
}

density_from_fit <- function(xgrid, fit, target = 45, n_int = 8000) {
  G  <- fit$G
  c1 <- fit$c1
  mu <- fit$scaling$mu
  s  <- fit$scaling$s
  
  zgrid <- (xgrid - mu) / s
  
  tr <- find_truncation(G, c1, target = target)
  zL <- tr$left; zR <- tr$right
  z0 <- tr$z_mode; s0 <- tr$s_mode
  
  unnorm <- function(z) exp(-(s_function(z, G, c1) - s0))
  
  z_int <- seq(zL, zR, length.out = n_int)
  Z <- trapz(z_int, unnorm(z_int))
  
  pz <- unnorm(zgrid) / Z
  (1 / s) * pz
}

L1_error <- function(xg, p1, p2) trapz(xg, abs(p1 - p2))

# ============================================================
# EXPERIMENT: n=200..20000, L1 + Dichtepanels (wie bei dir)
# ============================================================
set.seed(1)
mu0 <- 50; sd0 <- 10
ns <- c(200, 500, 1000, 2000, 5000, 10000, 20000)
R  <- 30
m  <- 1  # allgemein: m=1 ist Teilmodell; kein closed-form, kein Gaussian-fit

errs <- matrix(NA_real_, nrow = length(ns), ncol = R)
x_store   <- vector("list", length(ns))
fit_store <- vector("list", length(ns))

for (i in seq_along(ns)) {
  n <- ns[i]
  
  x_rep <- rnorm(n, mu0, sd0)
  x_store[[i]] <- x_rep
  fit_store[[i]] <- fit_score_matching_moments(
    x_rep, m = m,
    standardize = TRUE,
    col_scale = "sd",
    delta_lead = 1e-8,
    ridge = 1e-8,
    solver = "SCS"
  )
  
  for (r in 1:R) {
    x <- rnorm(n, mu0, sd0)
    
    fit <- fit_score_matching_moments(
      x, m = m,
      standardize = TRUE,
      col_scale = "sd",
      delta_lead = 1e-8,
      ridge = 1e-8,
      solver = "SCS"
    )
    
    # Evaluationsgrid (fair + stabil)
    xg <- seq(mu0 - 12 * sd0, mu0 + 12 * sd0, length.out = 4001)
    p_hat <- density_from_fit(xg, fit, target = 45, n_int = 8000)
    p_ref <- dnorm(xg, mu0, sd0)
    
    errs[i, r] <- L1_error(xg, p_hat, p_ref)
  }
  
  cat("n=", n, " median L1=", median(errs[i, ]), " status(rep)=", fit_store[[i]]$status, "\n")
}

# Plot 1: L1-Konvergenz
plot(ns, apply(errs, 1, median), log="x", type="b",
     pch=19, col="blue",
     xlab="Sample size n",
     ylab="Median L1 error",
     main=paste0("L1-Konvergenz (allgemein, m=", m, ")"))
grid()

# Plot 2: Dichte-Panels
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
for (i in seq_along(ns)) {
  n   <- ns[i]
  x   <- x_store[[i]]
  fit <- fit_store[[i]]
  
  sd_ref <- sd0; mu_ref <- mu0
  xgrid <- seq(min(x) - 4*sd_ref, max(x) + 4*sd_ref, length.out = 600)
  
  p_hat <- density_from_fit(xgrid, fit, target = 45, n_int = 8000)
  p_ref <- dnorm(xgrid, mu_ref, sd_ref)
  
  hist(x, breaks=30, freq=FALSE,
       col="grey90", border="grey70",
       main=paste0("n = ", n),
       xlab="", ylab="")
  lines(xgrid, p_hat, lwd=2, col="darkgreen")
  lines(xgrid, p_ref, lwd=2, col="red", lty=2)
  
  legend("topright",
         legend=c("Score Matching (allg.)", "True N(50,10)"),
         col=c("darkgreen", "red"),
         lwd=2, lty=c(1,2), bty="n", cex=0.8)
}
par(mfrow = c(1, 1))