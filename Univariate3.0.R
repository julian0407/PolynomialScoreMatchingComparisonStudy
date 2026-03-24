library(CVXR)
library(pracma)
library(sn)

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
# (2) Fast M/N features
#       N_ij(z) = z^(i+j)
#       M_ij(z) = z^(i+j+1)/(i+j+1)
# ============================================================
precompute_exponents <- function(m) {
  expoN <- as.vector(outer(0:(m-1), 0:(m-1), `+`))
  expoM <- expoN + 1
  denomM <- expoM
  list(expoN = expoN, expoM = expoM, denomM = denomM)
}

build_AB_fast <- function(z, m, exps) {
  A <- outer(z, exps$expoM, `^`)
  A <- sweep(A, 2, exps$denomM, "/")
  B <- outer(z, exps$expoN, `^`)
  list(A = A, B = B)
}

# ============================================================
# (3) Fast s(z) evaluation
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
# (4) Truncation for normalization
# ============================================================
find_truncation_fast <- function(G, c1, z_hint, target_drop = 45,
                                 grid_len = 1201, max_expand = 30) {
  lo <- as.numeric(quantile(z_hint, 1e-3))
  hi <- as.numeric(quantile(z_hint, 1 - 1e-3))
  span <- hi - lo
  if (!is.finite(span) || span <= 0) span <- sd(z_hint)
  if (!is.finite(span) || span <= 0) span <- 1
  
  lo <- lo - 0.5 * span
  hi <- hi + 0.5 * span
  
  z0_grid <- seq(lo, hi, length.out = grid_len)
  s0_grid <- svals_fast(z0_grid, G, c1)
  idx0 <- which.min(s0_grid)
  z0 <- z0_grid[idx0]
  s0 <- s0_grid[idx0]
  
  L <- lo
  R <- hi
  for (k in 1:max_expand) {
    if (svals_fast(L, G, c1) - s0 < target_drop) L <- L * 2 else break
  }
  for (k in 1:max_expand) {
    if (svals_fast(R, G, c1) - s0 < target_drop) R <- R * 2 else break
  }
  
  list(left = L, right = R, z_mode = z0, s_mode = s0)
}

# ============================================================
# (4.1) Density reconstruction from fit
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
  zL <- tr$left
  zR <- tr$right
  s0 <- tr$s_mode
  
  z_int <- seq(zL, zR, length.out = n_int)
  s_int <- svals_fast(z_int, G, c1)
  un_int <- exp(-(s_int - s0))
  Z <- trapz(z_int, un_int)
  
  s_g <- svals_fast(zgrid, G, c1)
  pz <- exp(-(s_g - s0)) / Z
  (1 / s) * pz
}

# ============================================================
# (4.2) Direkter Score-Loss (echte Definition)
#     1/2 * E[h(z) * (r_hat(z) - r_true(z))^2]
#
# score_true(x) muss - d/dx log p(x) auf x-Skala liefern
# Der Fit lebt auf z-Skala: z = (x - mu)/s
# Daher: r_true(z) = s * score_true(x)
# ============================================================
score_loss_direct_from_fit <- function(
    x,
    fit,
    score_true,
    h = function(z) rep(1, length(z))
) {
  stopifnot(is.function(score_true))
  
  mu <- fit$scaling$mu
  s  <- fit$scaling$s
  z  <- (x - mu) / s
  
  G  <- fit$G
  c1 <- fit$c1
  m  <- nrow(G)
  
  exps <- precompute_exponents(m)
  AB <- build_AB_fast(z, m, exps)
  
  gvec <- as.vector(G)
  
  # geschätzter Score auf z-Skala
  r_hat <- as.numeric(AB$A %*% gvec) + c1
  
  # wahrer Score auf x-Skala -> transformiert auf z-Skala
  r_true_z <- s * as.numeric(score_true(x))
  
  hz <- as.numeric(h(z))
  stopifnot(length(hz) == length(z))
  hz <- pmax(hz, 0)
  
  mean(0.5 * hz * (r_hat - r_true_z)^2)
}

# ============================================================
# (5) Fast fit via aggregated moments
# ============================================================
fit_score_matching_matrixG_moments_fast <- function(
    x,
    m,
    h = function(z) rep(1, length(z)),
    h_prime = function(z) rep(0, length(z)),
    standardize = TRUE,
    col_scale = c("sd", "maxabs", "none"),
    eps = 1e-8,
    delta_lead = 1e-6,
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
  
  b_bar <- as.vector(crossprod(B, hx) / n)
  
  if (col_scale == "none") {
    scale_vec <- rep(1, p)
  } else if (col_scale == "maxabs") {
    scale_vec <- apply(abs(A), 2, function(v) max(v, 1e-12))
  } else {
    scale_vec <- apply(A, 2, function(v) max(sd(v), 1e-12))
  }
  
  A_sc <- sweep(A, 2, scale_vec, "/")
  b_sc <- b_bar / scale_vec
  
  Aw <- A_sc * sqrt(hx)
  S  <- crossprod(Aw) / n
  t  <- as.vector(crossprod(A_sc, hx) / n)
  u  <- sum(hx) / n
  r  <- as.vector(crossprod(A_sc, hpx) / n)
  q  <- sum(hpx) / n
  
  K <- rbind(
    cbind(S, matrix(t, ncol = 1)),
    cbind(matrix(t, nrow = 1), matrix(u, 1, 1))
  )
  K <- 0.5 * (K + t(K))
  
  lin <- c(r + b_sc, q)
  
  G  <- Variable(m, m, PSD = TRUE)
  c1 <- Variable(1)
  gvec <- tryCatch(vec(G), error = function(e) reshape(G, c(p, 1)))
  y <- vstack(gvec, c1)
  
  lin_vec <- matrix(as.numeric(lin), ncol = 1)
  
  obj <- 0.5 * quad_form(y, Constant(K)) -
    sum_entries(multiply(Constant(lin_vec), y))
  
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
  
  g_sc  <- as.numeric(sol$getValue(gvec))
  g_org <- g_sc / scale_vec
  G_org <- matrix(g_org, nrow = m, ncol = m)
  
  list(
    G = G_org,
    c1 = as.numeric(sol$getValue(c1)),
    status = sol$status,
    solution = sol,
    scaling = list(
      mu = sc$mu,
      s = sc$s,
      standardize = standardize,
      col_scale = col_scale,
      scale_vec = scale_vec
    ),
    z_hint = z
  )
}

# ============================================================
# (6) Robust h, h'
# ============================================================
h_tanh_sq <- function(z) {
  tau <- 1
  u <- abs(z) / tau
  tanh(u)^2
}

h_tanh_sq_prime <- function(z) {
  tau <- 1
  u <- abs(z) / tau
  sech2 <- 1 / cosh(u)^2
  2 * tanh(u) * sech2 * (sign(z) / tau)
}

# ============================================================
# (7) L1 helper
# ============================================================
L1_error <- function(xg, p1, p2) {
  trapz(xg, abs(p1 - p2))
}

# ============================================================
# (8) Wrapper:
#     Density panels + letztes Panel = L1-Plot + Score-Loss-Plot
#     NEU:
#     - show_plots = TRUE/FALSE
#     - wenn FALSE: keine Density-Panels, keine L1/Score-Loss-Panels
# ============================================================
run_L1_and_density_panels <- function(
    name,
    ns,
    R,
    m,
    r_true = NULL,
    d_true = NULL,
    logf  = NULL,
    score_true = NULL,
    seed  = 1,
    
    # unabhängiges Testsample für direkten Score-Loss
    n_test_score = 10000,
    
    hint_n = 50000,
    eval_sd_mult  = 10,
    panel_sd_mult = 6,
    eval_grid_len  = 2001,
    panel_grid_len = 600,
    
    fit_control = list(
      standardize = TRUE,
      col_scale   = "sd",
      eps         = 1e-8,
      delta_lead  = 1e-6,
      ridge       = 0,
      scs_control = list(max_iters = 100000, eps = 1e-5, alpha = 1.8, verbose = FALSE)
    ),
    
    dens_control = list(target_drop = 45, n_int = 2500),
    
    slice_control = list(x0 = 0, w = 5, mstep = 80, burn = 800, thin = 1),
    true_control  = list(target_drop = 55, n_int = 5000),
    
    print_timing = TRUE,
    show_plots   = TRUE
) {
  stopifnot(length(ns) >= 1, R >= 1, m >= 1)
  stopifnot(!is.null(r_true) || !is.null(logf))
  stopifnot(!is.null(d_true) || !is.null(logf))
  stopifnot(is.function(score_true))
  
  # ----------------------------
  # Helpers nur falls logf
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
    if (!is.finite(span) || span <= 0) span <- 1
    
    lo <- lo - 0.5 * span
    hi <- hi + 0.5 * span
    
    x0_grid <- seq(lo, hi, length.out = 1201)
    lf0 <- vapply(x0_grid, logf, numeric(1))
    i0 <- which.max(lf0)
    lf_mode <- lf0[i0]
    
    L <- lo
    R <- hi
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
  # Sampling wrapper
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
  # Hint-sample -> stabile Grids
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
  
  t_fit  <- matrix(NA_real_, nrow = length(ns), ncol = R)
  t_dens <- matrix(NA_real_, nrow = length(ns), ncol = R)
  
  for (i in seq_along(ns)) {
    n <- ns[i]
    
    # ein Replikat nur für optionale Panels
    x_rep <- draw_sample(n)
    x_store[[i]] <- x_rep
    fit_store[[i]] <- do.call(
      fit_score_matching_matrixG_moments_fast,
      c(list(x = x_rep, m = m, h = h_tanh_sq, h_prime = h_tanh_sq_prime), fit_control)
    )
    
    # R Wiederholungen für L1 + direkten Score-Loss + Timing
    for (r in 1:R) {
      x_train <- draw_sample(n)
      x_test  <- draw_sample(n_test_score)
      
      tf <- system.time({
        fit <- do.call(
          fit_score_matching_matrixG_moments_fast,
          c(list(x = x_train, m = m, h = h_tanh_sq, h_prime = h_tanh_sq_prime), fit_control)
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
      
      score_losses[i, r] <- score_loss_direct_from_fit(
        x = x_test,
        fit = fit,
        score_true = score_true,
        h = h_tanh_sq
      )
      
      t_fit[i, r]  <- tf
      t_dens[i, r] <- td
    }
    
    if (print_timing) {
      cat(
        "n=", n, " m=", m,
        " | status(rep)=", fit_store[[i]]$status,
        " | median L1=", median(errs[i, ], na.rm = TRUE),
        " | median DirectScoreLoss=", median(score_losses[i, ], na.rm = TRUE),
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
  # Plot nur falls show_plots = TRUE
  # ----------------------------
  if (show_plots) {
    k <- length(ns)
    total_panels <- k + 2
    ncolp <- ceiling(sqrt(total_panels))
    nrowp <- ceiling(total_panels / ncolp)
    
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mfrow = c(nrowp, ncolp), mar = c(3, 3, 2, 1))
    
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
    
    l1_med <- apply(errs, 1, median, na.rm = TRUE)
    plot(ns, l1_med,
         log = "x", type = "b", pch = 19, col = "blue",
         xlab = "n (log)", ylab = "Median L1",
         main = "L1 (median over R)")
    grid()
    
    score_loss_med <- apply(score_losses, 1, median, na.rm = TRUE)
    plot(ns, score_loss_med,
         log = "x", type = "b", pch = 19, col = "purple",
         xlab = "n (log)", ylab = "Median Direct Score-Loss",
         main = "Direct Score-Loss (median over R)")
    grid()
  }
  
  invisible(list(
    name = name,
    ns = ns,
    R = R,
    m = m,
    mu0 = mu0,
    sd0 = sd0,
    xg_eval = xg_eval,
    p_ref_eval = p_ref_eval,
    errs = errs,
    score_losses = score_losses,
    t_fit = t_fit,
    t_dens = t_dens,
    timing_summary = timing_summary,
    x_store = x_store,
    fit_store = fit_store
  ))
}

# ============================================================
# (9) Helper für Vergleich über mehrere m:
#     Für jedes n: Plot gegen m
# ============================================================
plot_error_vs_m_by_n <- function(res_list, ms,
                                 metric = c("L1", "score"),
                                 use_median = TRUE) {
  metric <- match.arg(metric)
  
  stopifnot(length(res_list) == length(ms))
  ns_ref <- res_list[[1]]$ns
  stopifnot(all(sapply(res_list, function(res) identical(res$ns, ns_ref))))
  
  err_mat <- sapply(res_list, function(res) {
    if (metric == "L1") {
      apply(res$errs, 1, median, na.rm = TRUE)
    } else {
      apply(res$score_losses, 1, median, na.rm = TRUE)
    }
  })
  
  if (is.vector(err_mat)) {
    err_mat <- matrix(err_mat, nrow = length(ns_ref), ncol = length(ms))
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  k <- length(ns_ref)
  ncolp <- ceiling(sqrt(k))
  nrowp <- ceiling(k / ncolp)
  par(mfrow = c(nrowp, ncolp), mar = c(4, 4, 3, 1))
  
  ylab_txt <- if (metric == "L1") "Median L1" else "Median Direct Score-Loss"
  main_prefix <- if (metric == "L1") "L1 vs m" else "Score-Loss vs m"
  
  for (i in seq_along(ns_ref)) {
    plot(ms, err_mat[i, ],
         type = "b", pch = 19,
         xlab = "m", ylab = ylab_txt,
         main = paste0(main_prefix, "\n n = ", ns_ref[i]))
    grid()
  }
  
  invisible(err_mat)
}

# ============================================================
# TESTS
# ============================================================

ns_default <- c(200, 500, 1000, 2000, 5000, 10000, 20000)

# ------------------------------------------------------------
# (A) Gaussian N(mu, sd^2)
# ------------------------------------------------------------
mu_g <- 50
sd_g <- 10

score_gauss <- function(x) {
  (x - mu_g) / sd_g^2
}

res_gauss <- run_L1_and_density_panels(
  name = "Gaussian N(50,10^2)",
  ns = ns_default,
  R  = 20,
  m  = 1,
  r_true = function(n) rnorm(n, mu_g, sd_g),
  d_true = function(x) dnorm(x, mu_g, sd_g),
  score_true = score_gauss,
  seed = 1,
  n_test_score = 10000,
  show_plots = TRUE
)

# ------------------------------------------------------------
# (B) Logistic(0,1)
# ------------------------------------------------------------
score_logistic <- function(x) {
  2 / (1 + exp(-x)) - 1
}

res_logistic <- run_L1_and_density_panels(
  name = "Logistic(0,1)",
  ns = ns_default,
  R  = 20,
  m  = 3,
  r_true = function(n) rlogis(n, location = 0, scale = 1),
  d_true = function(x) dlogis(x, location = 0, scale = 1),
  score_true = score_logistic,
  seed = 1,
  n_test_score = 10000,
  show_plots = TRUE
)

# ------------------------------------------------------------
# (C) Skew-Normal
# ------------------------------------------------------------
score_skewnorm <- function(x, xi, omega, alpha) {
  z <- (x - xi) / omega
  za <- alpha * z
  ratio <- dnorm(za) / pnorm(za)
  (x - xi) / omega^2 - (alpha / omega) * ratio
}

xi <- 0
omega <- 1
alpha <- 5

score_true_skew <- function(x) score_skewnorm(x, xi, omega, alpha)

res_skewnorm <- run_L1_and_density_panels(
  name = paste0("Skew-Normal (alpha=", alpha, ")"),
  ns = ns_default,
  R  = 20,
  m  = 3,
  r_true = function(n) rsn(n, xi = xi, omega = omega, alpha = alpha),
  d_true = function(x) dsn(x, xi = xi, omega = omega, alpha = alpha),
  score_true = score_true_skew,
  seed = 1,
  n_test_score = 10000,
  show_plots = TRUE
)

# ------------------------------------------------------------
# (D) Quartic: p(x) proportional zu exp(-a x^4 - b x^2 - c x)
# ------------------------------------------------------------
a <- 0.02
b <- 0.5
c0 <- 0.0

logf_quartic <- function(x) {
  -(a * x^4 + b * x^2 + c0 * x)
}

score_quartic <- function(x) {
  4 * a * x^3 + 2 * b * x + c0
}

res_quartic <- run_L1_and_density_panels(
  name = "Quartic: exp(-a x^4 - b x^2 - c x)",
  ns = ns_default,
  R  = 20,
  m  = 2,
  logf = logf_quartic,
  score_true = score_quartic,
  seed = 1,
  n_test_score = 10000,
  show_plots = TRUE
)

# ------------------------------------------------------------
# (E) Sextic: p(x) proportional zu exp(-a x^6 - b x^4 - d x^2 - c x)
# ------------------------------------------------------------
a <- 0.002
b <- 0.02
d <- 0.3
c0 <- 0.0

logf_sextic <- function(x) {
  -(a * x^6 + b * x^4 + d * x^2 + c0 * x)
}

score_sextic <- function(x) {
  6 * a * x^5 + 4 * b * x^3 + 2 * d * x + c0
}

res_sextic <- run_L1_and_density_panels(
  name = "Sextic: exp(-a x^6 - b x^4 - d x^2 - c x)",
  ns = ns_default,
  R  = 20,
  m  = 3,
  logf = logf_sextic,
  score_true = score_sextic,
  seed = 1,
  n_test_score = 10000,
  show_plots = TRUE
)

# ------------------------------------------------------------
# (F) Gumbel(0,1): log-concave, misspecified
#     1) normaler Test ohne Panels
#     2) Vergleich über mehrere m und n
# ------------------------------------------------------------
rgumbel <- function(n, mu = 0, beta = 1) {
  mu - beta * log(-log(runif(n)))
}

dgumbel <- function(x, mu = 0, beta = 1) {
  z <- (x - mu) / beta
  (1 / beta) * exp(-(z + exp(-z)))
}

score_gumbel <- function(x) {
  1 - exp(-x)
}

gumbel_ns <- ns_default

res_gumbel <- run_L1_and_density_panels(
  name = "Gumbel(0,1) -- log-concave, misspecified",
  ns = gumbel_ns,
  R  = 20,
  m  = 3,
  r_true = function(n) rgumbel(n),
  d_true = function(x) dgumbel(x),
  score_true = score_gumbel,
  seed = 1,
  n_test_score = 10000,
  show_plots = FALSE
)

# ------------------------------------------------------------
# Vergleich über m für Gumbel:
# Für jedes n ein Plot mit x-Achse = m und y-Achse = Fehler
# ------------------------------------------------------------
ms <- c(2, 3, 4, 5, 6)
gumbel_ns_compare <- c(5000, 10000, 20000, 100000, 1000000)

res_gumbel_by_m <- setNames(vector("list", length(ms)), paste0("m_", ms))

for (j in seq_along(ms)) {
  mj <- ms[j]
  
  cat("\n============================\n")
  cat("Running Gumbel comparison for m =", mj, "\n")
  cat("============================\n")
  
  res_gumbel_by_m[[j]] <- run_L1_and_density_panels(
    name = paste0("Gumbel(0,1) -- log-concave, misspecified, m=", mj),
    ns = gumbel_ns_compare,
    R  = 20,
    m  = mj,
    r_true = function(n) rgumbel(n),
    d_true = function(x) dgumbel(x),
    score_true = score_gumbel,
    seed = 1,
    n_test_score = 10000,
    show_plots = FALSE
  )
}

gumbel_summary <- do.call(rbind, lapply(seq_along(ms), function(j) {
  res <- res_gumbel_by_m[[j]]
  
  data.frame(
    m = ms[j],
    n = res$ns,
    median_L1 = apply(res$errs, 1, median, na.rm = TRUE),
    median_score_loss = apply(res$score_losses, 1, median, na.rm = TRUE)
  )
}))

print(gumbel_summary)

# ------------------------------------------------------------
# Separate Plot-Figur 1: L1 vs m, für jedes n ein eigenes Panel
# ------------------------------------------------------------
plot_error_vs_m_by_n(
  res_list = res_gumbel_by_m,
  ms = ms,
  metric = "L1"
)

# ------------------------------------------------------------
# Separate Plot-Figur 2: Score-Loss vs m, für jedes n ein eigenes Panel
# ------------------------------------------------------------
plot_error_vs_m_by_n(
  res_list = res_gumbel_by_m,
  ms = ms,
  metric = "score"
)