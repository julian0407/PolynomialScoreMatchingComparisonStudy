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
# (4.2) Empirical score matching score from fitted model
#     J_n(fit) = mean( 0.5 h(z) s'(z)^2 + h'(z) s'(z) - h(z) s''(z) )
# where z is on the standardized scale used in the fit
# ============================================================
score_matching_score_from_fit <- function(
    x, fit,
    h = function(z) rep(1, length(z)),
    h_prime = function(z) rep(0, length(z))
) {
  # same scale as in fitting
  mu <- fit$scaling$mu
  s  <- fit$scaling$s
  z  <- (x - mu) / s
  
  G  <- fit$G
  c1 <- fit$c1
  m  <- nrow(G)
  
  exps <- precompute_exponents(m)
  AB <- build_AB_fast(z, m, exps)
  
  gvec <- as.vector(G)
  s1 <- as.numeric(AB$A %*% gvec) + c1   # s'(z)
  s2 <- as.numeric(AB$B %*% gvec)        # s''(z)
  
  hz  <- as.numeric(h(z))
  hpz <- as.numeric(h_prime(z))
  hz  <- pmax(hz, 0)
  
  mean(0.5 * hz * s1^2 + hpz * s1 - hz * s2)
}

# ============================================================
# (4.3) Empirical score matching score of the true model
#     using USER-SUPPLIED true derivatives on x-scale
#
# Required:
#   score_true(x)       = - d/dx log p(x)
#   score_true_prime(x) = - d^2/dx^2 log p(x)
#
# The fit lives on z-scale: z = (x - mu)/s
# Hence
#   ds_true/dz   = s * score_true(x)
#   d2s_true/dz2 = s^2 * score_true_prime(x)
# ============================================================
score_matching_score_true <- function(
    x, fit,
    score_true,
    score_true_prime,
    h = function(z) rep(1, length(z)),
    h_prime = function(z) rep(0, length(z))
) {
  stopifnot(is.function(score_true), is.function(score_true_prime))
  
  mu <- fit$scaling$mu
  s  <- fit$scaling$s
  z  <- (x - mu) / s
  
  hz  <- as.numeric(h(z))
  hpz <- as.numeric(h_prime(z))
  hz  <- pmax(hz, 0)
  
  # supplied on x-scale, converted to z-scale
  s1_true <- s^1 * as.numeric(score_true(x))
  s2_true <- s^2 * as.numeric(score_true_prime(x))
  
  mean(0.5 * hz * s1_true^2 + hpz * s1_true - hz * s2_true)
}

# ============================================================
# (4.4) Empirical score loss
# ============================================================
score_loss_from_fit <- function(
    x, fit,
    score_true,
    score_true_prime,
    h = function(z) rep(1, length(z)),
    h_prime = function(z) rep(0, length(z))
) {
  sm_fit <- score_matching_score_from_fit(
    x = x, fit = fit, h = h, h_prime = h_prime
  )
  
  sm_true <- score_matching_score_true(
    x = x, fit = fit,
    score_true = score_true,
    score_true_prime = score_true_prime,
    h = h, h_prime = h_prime
  )
  
  sm_fit - sm_true
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
# (7) Test function
# ============================================================
# ============================================================
# (7) Empirical convergence: L1 + density panels (fast)
# ============================================================
L1_error <- function(xg, p1, p2) trapz(xg, abs(p1 - p2))

library(CVXR)
library(pracma)
library(sn)   # nur wenn du Skew-Normal nutzt

# ============================================================
# Wrapper: Density Panels + letztes Panel = L1-Plot
#         + Timing: Median fit/dens pro n
# ============================================================
run_L1_and_density_panels <- function(
    name,
    ns,
    R,
    m,
    r_true = NULL,               # function(n) -> sample
    d_true = NULL,               # function(x) -> true density on x
    logf  = NULL,                # function(x) -> log unnormalized density (if d_true missing)
    score_true = NULL,           # REQUIRED for score-loss: - d/dx log p(x)
    score_true_prime = NULL,     # REQUIRED for score-loss: - d^2/dx^2 log p(x)
    seed  = 1,
    
    # Grid-Definition stabil über "Hint"-Sample:
    hint_n = 50000,
    eval_sd_mult  = 10,
    panel_sd_mult = 6,
    eval_grid_len  = 2001,
    panel_grid_len = 600,
    
    # Fit-Parameter (Defaults wie in deinen Tests)
    fit_control = list(
      standardize = TRUE,
      col_scale   = "sd",
      eps         = 1e-8,
      delta_lead  = 1e-6,
      ridge       = 0,
      scs_control = list(max_iters = 100000, eps = 1e-5, alpha = 1.8, verbose = FALSE)
    ),
    
    # Dichte-Rekonstruktion aus Fit
    dens_control = list(target_drop = 45, n_int = 2500),
    
    # Nur relevant, falls logf benutzt wird (Sampling + Referenz-Normierung)
    slice_control = list(x0 = 0, w = 5, mstep = 80, burn = 800, thin = 1),
    true_control  = list(target_drop = 55, n_int = 5000),
    
    print_timing = TRUE
) {
  stopifnot(length(ns) >= 1, R >= 1, m >= 1)
  stopifnot(!is.null(r_true) || !is.null(logf))
  stopifnot(!is.null(d_true) || !is.null(logf))
  stopifnot(is.function(score_true), is.function(score_true_prime))
  
  # ----------------------------
  # Helpers (nur falls logf)
  # ----------------------------
  rslice1 <- function(n, logf, x0 = 0, w = 1, mstep = 50, burn = 200, thin = 1) {
    total <- burn + n * thin
    x <- numeric(total)
    x[1] <- x0
    
    for (t in 2:total) {
      xt <- x[t - 1]
      logy <- logf(xt) + log(runif(1))
      
      u <- runif(1)
      L <- xt - w * u
      R <- L + w
      J <- floor(runif(1, 0, mstep))
      K <- (mstep - 1) - J
      
      while (J > 0 && logf(L) > logy) { L <- L - w; J <- J - 1 }
      while (K > 0 && logf(R) > logy) { R <- R + w; K <- K - 1 }
      
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
  
  true_density_on_grid <- function(xgrid, logf, x_hint, target_drop = 45, n_int = 4000) {
    lo <- as.numeric(quantile(x_hint, 1e-3))
    hi <- as.numeric(quantile(x_hint, 1 - 1e-3))
    span <- hi - lo
    if (!is.finite(span) || span <= 0) span <- sd(x_hint)
    lo <- lo - 0.5 * span
    hi <- hi + 0.5 * span
    
    x0_grid <- seq(lo, hi, length.out = 1201)
    lf0 <- vapply(x0_grid, logf, numeric(1))
    i0 <- which.max(lf0)
    lf_mode <- lf0[i0]
    
    L <- lo; R <- hi
    for (k in 1:30) if ((lf_mode - logf(L)) < target_drop) L <- L * 2 else break
    for (k in 1:30) if ((lf_mode - logf(R)) < target_drop) R <- R * 2 else break
    
    x_int <- seq(L, R, length.out = n_int)
    lf_int <- vapply(x_int, logf, numeric(1))
    un_int <- exp(lf_int - lf_mode)
    Z <- pracma::trapz(x_int, un_int)
    
    lf_g <- vapply(xgrid, logf, numeric(1))
    exp(lf_g - lf_mode) / Z
  }
  
  # ----------------------------
  # Sampling-Wrapper
  # ----------------------------
  draw_sample <- function(n) {
    if (!is.null(r_true)) return(r_true(n))
    rslice1(
      n = n, logf = logf,
      x0 = slice_control$x0, w = slice_control$w, mstep = slice_control$mstep,
      burn = slice_control$burn, thin = slice_control$thin
    )
  }
  
  # ----------------------------
  # Setup: Hint-sample -> stabile Grids
  # ----------------------------
  set.seed(seed)
  x_hint <- draw_sample(hint_n)
  mu0 <- mean(x_hint)
  sd0 <- sd(x_hint)
  if (!is.finite(sd0) || sd0 <= 0) sd0 <- 1
  
  xg_eval  <- seq(mu0 - eval_sd_mult * sd0,  mu0 + eval_sd_mult * sd0,  length.out = eval_grid_len)
  xg_panel <- seq(mu0 - panel_sd_mult * sd0, mu0 + panel_sd_mult * sd0, length.out = panel_grid_len)
  
  p_ref_eval <- if (!is.null(d_true)) {
    d_true(xg_eval)
  } else {
    true_density_on_grid(
      xg_eval, logf, x_hint,
      target_drop = true_control$target_drop,
      n_int       = true_control$n_int
    )
  }
  
  # ----------------------------
  # Storage
  # ----------------------------
  errs <- matrix(NA_real_, nrow = length(ns), ncol = R)
  score_losses <- matrix(NA_real_, nrow = length(ns), ncol = R)
  x_store   <- vector("list", length(ns))
  fit_store <- vector("list", length(ns))
  
  # Timing (Sekunden)
  t_fit  <- matrix(NA_real_, nrow = length(ns), ncol = R)
  t_dens <- matrix(NA_real_, nrow = length(ns), ncol = R)
  
  for (i in seq_along(ns)) {
    n <- ns[i]
    
    # 1 Replicate für Panels
    x_rep <- draw_sample(n)
    x_store[[i]] <- x_rep
    fit_store[[i]] <- do.call(
      fit_score_matching_matrixG_moments_fast,
      c(list(x = x_rep, m = m, h = h_tanh_sq, h_prime = h_tanh_sq_prime), fit_control)
    )
    
    # R Wiederholungen für L1 + Score-Loss + Timing
    for (r in 1:R) {
      x <- draw_sample(n)
      
      tf <- system.time({
        fit <- do.call(
          fit_score_matching_matrixG_moments_fast,
          c(list(x = x, m = m, h = h_tanh_sq, h_prime = h_tanh_sq_prime), fit_control)
        )
      })[["elapsed"]]
      
      td <- system.time({
        p_hat <- density_from_fit_fast_trunc(
          xg_eval, fit,
          target_drop = dens_control$target_drop,
          n_int       = dens_control$n_int
        )
      })[["elapsed"]]
      
      errs[i, r] <- L1_error(xg_eval, p_hat, p_ref_eval)
      
      score_losses[i, r] <- score_loss_from_fit(
        x = x,
        fit = fit,
        score_true = score_true,
        score_true_prime = score_true_prime,
        h = h_tanh_sq,
        h_prime = h_tanh_sq_prime
      )
      
      t_fit[i, r]  <- tf
      t_dens[i, r] <- td
    }
    
    if (print_timing) {
      cat(
        "n=", n, " m=", m,
        " | status(rep)=", fit_store[[i]]$status,
        " | median L1=", median(errs[i, ], na.rm = TRUE),
        " | median ScoreLoss=", median(score_losses[i, ], na.rm = TRUE),
        " | median times [fit/dens] = ",
        sprintf("%.4f", median(t_fit[i, ],  na.rm = TRUE)), "/",
        sprintf("%.4f", median(t_dens[i, ], na.rm = TRUE)),
        " sec\n",
        sep = ""
      )
    }
  }
  
  timing_summary <- data.frame(
    n = ns,
    m = rep(m, length(ns)),
    median_fit   = apply(t_fit,  1, median, na.rm = TRUE),
    median_dens  = apply(t_dens, 1, median, na.rm = TRUE),
    median_total = apply(t_fit + t_dens, 1, median, na.rm = TRUE)
  )
  
  if (print_timing) {
    cat("\nTiming summary (median seconds per model fit):\n")
    print(timing_summary, row.names = FALSE)
  }
  
  # ----------------------------
  # Plot: Density panels + L1 + Score-Loss
  # ----------------------------
  k <- length(ns)
  total_panels <- k + 2
  ncolp <- ceiling(sqrt(total_panels))
  nrowp <- ceiling(total_panels / ncolp)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(nrowp, ncolp), mar = c(3, 3, 2, 1))
  
  # Density panels
  for (i in seq_along(ns)) {
    x   <- x_store[[i]]
    fit <- fit_store[[i]]
    
    p_hat <- density_from_fit_fast_trunc(
      xg_panel, fit,
      target_drop = dens_control$target_drop,
      n_int       = dens_control$n_int
    )
    
    p_ref <- if (!is.null(d_true)) {
      d_true(xg_panel)
    } else {
      true_density_on_grid(
        xg_panel, logf, x_hint,
        target_drop = true_control$target_drop,
        n_int       = true_control$n_int
      )
    }
    
    hist(x, breaks = 30, freq = FALSE,
         col = "grey90", border = "grey70",
         main = paste0("n = ", ns[i], "\n", fit$status),
         xlab = "", ylab = "")
    lines(xg_panel, p_hat, lwd = 2, col = "darkgreen")
    lines(xg_panel, p_ref, lwd = 2, col = "red", lty = 2)
  }
  
  # Summary panel 1: L1
  l1_med <- apply(errs, 1, median, na.rm = TRUE)
  plot(ns, l1_med,
       log = "x", type = "b", pch = 19, col = "blue",
       xlab = "n (log)", ylab = "Median L1",
       main = "L1 (median over R)")
  grid()
  
  # Summary panel 2: Score-Loss
  score_loss_med <- apply(score_losses, 1, median, na.rm = TRUE)
  plot(ns, score_loss_med,
       log = "x", type = "b", pch = 19, col = "purple",
       xlab = "n (log)", ylab = "Median Score-Loss",
       main = "Score-Loss (median over R)")
  grid()
  
  invisible(list(
    name = name, ns = ns, R = R, m = m,
    mu0 = mu0, sd0 = sd0,
    xg_eval = xg_eval, p_ref_eval = p_ref_eval,
    errs = errs,
    score_losses = score_losses,
    t_fit = t_fit, t_dens = t_dens,
    timing_summary = timing_summary,
    x_store = x_store,
    fit_store = fit_store
  ))
}

# ============================================================
# TESTS (kurz) – kommentiert nach Verteilung
# ============================================================

ns_default <- c(200, 500, 1000, 2000, 5000, 10000, 20000)

# ------------------------------------------------------------
# (A) Gaussian N(mu, sd^2)
# ------------------------------------------------------------
mu_g <- 50
sd_g <- 10

score_gauss <- function(x) (x - mu_g) / sd_g^2
score_gauss_prime <- function(x) rep(1 / sd_g^2, length(x))

res_gauss <- run_L1_and_density_panels(
  name = "Gaussian N(50,10^2)",
  ns = c(200, 500, 1000, 2000, 5000, 10000, 20000, 100000),
  R  = 20,
  m  = 1,
  r_true = function(n) rnorm(n, mu_g, sd_g),
  d_true = function(x) dnorm(x, mu_g, sd_g),
  score_true = score_gauss,
  score_true_prime = score_gauss_prime,
  seed = 1
)

# ------------------------------------------------------------
# (B) Logistic(0,1)
# ------------------------------------------------------------
score_logistic <- function(x) {
  2 / (1 + exp(-x)) - 1
}

score_logistic_prime <- function(x) {
  p <- 1 / (1 + exp(-x))
  2 * p * (1 - p)
}

res_logistic <- run_L1_and_density_panels(
  name = "Logistic(0,1)",
  ns = ns_default,
  R  = 20,
  m  = 3,
  r_true = function(n) rlogis(n, location = 0, scale = 1),
  d_true = function(x) dlogis(x, location = 0, scale = 1),
  score_true = score_logistic,
  score_true_prime = score_logistic_prime,
  seed = 1
)

# ------------------------------------------------------------
# (C) Skew-Normal (sn::rsn/dsn)
# ------------------------------------------------------------
score_skewnorm <- function(x, xi, omega, alpha) {
  z <- (x - xi) / omega
  za <- alpha * z
  
  ratio <- dnorm(za) / pnorm(za)
  
  (x - xi) / omega^2 - (alpha / omega) * ratio
}

score_skewnorm_prime <- function(x, xi, omega, alpha) {
  z <- (x - xi) / omega
  za <- alpha * z
  
  ratio <- dnorm(za) / pnorm(za)
  
  (1 / omega^2) +
    (alpha^2 / omega^2) *
    (za * ratio + ratio^2)
}

xi <- 0
omega <- 1
alpha <- 5

score_true <- function(x) score_skewnorm(x, xi, omega, alpha)
score_true_prime <- function(x) score_skewnorm_prime(x, xi, omega, alpha)

res_skewnorm <- run_L1_and_density_panels(
  name = paste0("Skew-Normal (alpha=", alpha, ")"),
  ns = c(200, 500, 1000, 2000, 5000, 10000, 20000, 100000),
  R  = 20,
  m  = 10,
  r_true = function(n) rsn(n, xi = xi, omega = omega, alpha = alpha),
  d_true = function(x) dsn(x, xi = xi, omega = omega, alpha = alpha),
  score_true = score_true,
  score_true_prime = score_true_prime,
  seed = 1
)

# ------------------------------------------------------------
# (D) Quartic: p(x) ∝ exp(-a x^4 - b x^2 - c x)  (log target)
# ------------------------------------------------------------
a <- 0.02
b <- 0.5
c0 <- 0.0

logf_quartic <- function(x) -(a * x^4 + b * x^2 + c0 * x)

score_quartic <- function(x) 4 * a * x^3 + 2 * b * x + c0
score_quartic_prime <- function(x) 12 * a * x^2 + 2 * b

res_quartic <- run_L1_and_density_panels(
  name = "Quartic: exp(-a x^4 - b x^2 - c x)",
  ns = ns_default,
  R  = 20,
  m  = 2,
  logf = logf_quartic,
  score_true = score_quartic,
  score_true_prime = score_quartic_prime,
  seed = 1
)

# ------------------------------------------------------------
# (E) Sextic: p(x) ∝ exp(-a x^6 - b x^4 - d x^2 - c x) (log target)
# ------------------------------------------------------------
a <- 0.002
b <- 0.02
d <- 0.3
c0 <- 0.0

logf_sextic <- function(x) -(a * x^6 + b * x^4 + d * x^2 + c0 * x)

score_sextic <- function(x) 6 * a * x^5 + 4 * b * x^3 + 2 * d * x + c0
score_sextic_prime <- function(x) 30 * a * x^4 + 12 * b * x^2 + 2 * d

res_sextic <- run_L1_and_density_panels(
  name = "Sextic: exp(-a x^6 - b x^4 - d x^2 - c x)",
  ns = ns_default,
  R  = 20,
  m  = 3,
  logf = logf_sextic,
  score_true = score_sextic,
  score_true_prime = score_sextic_prime,
  seed = 1
)
