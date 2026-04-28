# ============================================================
# Test_LogConcave_Density_3D.R
#
# Separates Testskript für dein modulares multivariates
# Polynomial-Score-Matching-Setup.
#
# Es simuliert 2D-Gaussian-Daten, fitttet
#   (1) Polynomial SM ohne Log-Konkavität
#   (2) Polynomial SM mit Log-Konkavität
# und vergleicht die geschätzten Dichten in 3D gegen die
# wahre Gaussian-Dichte.
#
# Erwartete Dateien im gleichen Ordner:
#   - helper_functions.R
#   - Multivariate_Polynomial_Score_Matching_LogConcave.R
# ============================================================

# ------------------------------------------------------------
# (0) Dateien sourcen
# ------------------------------------------------------------
source("helper_functions.R")
source("Multivariate_Polynomial_Score_Matching_LogConcave_2.R")

# ------------------------------------------------------------
# (1) Kleine Hilfsfunktionen
# ------------------------------------------------------------

# Bivariate Gaussian-Simulation ohne weitere Paketabhängigkeit
rmvnorm_chol <- function(n, mu, Sigma) {
  mu <- as.numeric(mu)
  d <- length(mu)
  if (!all(dim(Sigma) == c(d, d))) stop("Sigma has wrong dimension.")

  R <- chol(Sigma)
  Z <- matrix(rnorm(n * d), nrow = n, ncol = d)
  X <- Z %*% R
  X <- sweep(X, 2, mu, FUN = "+")
  X
}

# Wahre bivariate Gaussian-Dichte
true_dmvnorm_2d <- function(x, mu, Sigma) {
  x <- as.matrix(x)
  mu <- as.numeric(mu)
  d <- ncol(x)
  if (d != 2L) stop("This helper is intended for d = 2.")

  Sigma_inv <- solve(Sigma)
  det_Sigma <- det(Sigma)
  xc <- sweep(x, 2, mu, FUN = "-")
  quad <- rowSums((xc %*% Sigma_inv) * xc)

  (2 * pi)^(-d / 2) * det_Sigma^(-1 / 2) * exp(-0.5 * quad)
}

# Polynom-Basis auswerten: liefert Potential V(z) = sum_j theta_j b_j(z)
eval_poly_potential_mv_basic <- function(z, fit) {
  z <- as_obs_matrix(z)
  if (ncol(z) != fit$d) stop("Dimension mismatch in eval_poly_potential_mv_basic().")

  n <- nrow(z)
  p <- length(fit$basis)
  max_deg <- fit$m
  P <- build_power_cache_mv_basic(z, max_deg = max_deg)

  B <- matrix(0, nrow = n, ncol = p)

  for (col in seq_len(p)) {
    b <- fit$basis[[col]]

    if (b$type == "uni") {
      B[, col] <- get_power_mv_basic(P, b$j, b$r)
    } else if (b$type == "pair") {
      B[, col] <- get_power_mv_basic(P, b$j, b$r) * get_power_mv_basic(P, b$k, b$s)
    } else {
      stop("Unknown basis type.")
    }
  }

  as.numeric(B %*% fit$theta)
}

# Auswertung der geschätzten Dichte auf x-Skala
# score = grad(V), log f = const - V
predict_density_mv_basic_grid <- function(newx, fit) {
  newx <- as_obs_matrix(newx)
  if (ncol(newx) != fit$d) stop("Dimension mismatch in predict_density_mv_basic_grid().")

  if (isTRUE(fit$scaling$standardize)) {
    z <- apply_scaling_matrix(newx, fit$scaling)$z
    jacobian_const <- 1 / prod(fit$scaling$scale)
  } else {
    z <- newx
    jacobian_const <- 1
  }

  Vz <- eval_poly_potential_mv_basic(z, fit)

  # numerisch stabil
  dens_unnorm <- exp(-(Vz - min(Vz))) * jacobian_const
  dens_unnorm
}

# numerische Grid-Normalisierung auf 2D-Rechteck
normalize_density_on_grid <- function(dens_vec, x1_seq, x2_seq) {
  nx <- length(x1_seq)
  ny <- length(x2_seq)
  if (length(dens_vec) != nx * ny) stop("dens_vec has wrong length.")

  dx <- mean(diff(x1_seq))
  dy <- mean(diff(x2_seq))

  zmat <- matrix(dens_vec, nrow = nx, ncol = ny)
  mass <- sum(zmat) * dx * dy
  if (!is.finite(mass) || mass <= 0) stop("Normalization failed: non-positive mass.")

  zmat / mass
}

# 3D-Plot-Helfer
plot_density_surface <- function(zmat, x1_seq, x2_seq, main = "") {
  persp(
    x = x1_seq,
    y = x2_seq,
    z = zmat,
    theta = 35,
    phi = 25,
    expand = 0.7,
    ticktype = "detailed",
    shade = 0.45,
    border = NA,
    xlab = "x1",
    ylab = "x2",
    zlab = "density",
    main = main
  )
}

# ------------------------------------------------------------
# (2) Hauptfunktion für den Test
# ------------------------------------------------------------
run_logconcave_density_3d_test <- function(
    n = 1000,
    mu_true = c(0, 0),
    Sigma_true = matrix(c(1, 0.5,
                          0.5, 1.5), 2, 2),
    m = 3,
    include_interactions = TRUE,
    standardize = TRUE,
    ridge = 1e-8,
    solver = "solve",
    log_concave_lc_method = "grid",
    lc_grid_size = 5,
    lc_max_points = 300,
    lc_box_expand = 0.25,
    lc_penalty = 1e4,
    grid_n = 80,
    grid_sd_expand = 3.5,
    seed = 123
) {
  set.seed(seed)

  mu_true <- as.numeric(mu_true)
  if (length(mu_true) != 2L) stop("This test script is intended for d = 2.")
  if (!all(dim(Sigma_true) == c(2, 2))) stop("Sigma_true must be 2x2.")

  # ----------------------------------------------------------
  # Daten simulieren
  # ----------------------------------------------------------
  x <- rmvnorm_chol(n = n, mu = mu_true, Sigma = Sigma_true)

  # ----------------------------------------------------------
  # Modelle fitten
  # ----------------------------------------------------------
  fit_plain <- fit_score_matching_mv_basic(
    x,
    m = m,
    include_interactions = include_interactions,
    standardize = standardize,
    ridge = ridge,
    solver = solver,
    log_concave = FALSE
  )

  fit_lc <- fit_score_matching_mv_basic(
    x,
    m = m,
    include_interactions = include_interactions,
    standardize = standardize,
    ridge = ridge,
    solver = solver,
    log_concave = TRUE,
    lc_method = log_concave_lc_method,
    lc_grid_size = lc_grid_size,
    lc_max_points = lc_max_points,
    lc_box_expand = lc_box_expand,
    lc_penalty = lc_penalty
  )

  # ----------------------------------------------------------
  # Auswertegitter definieren
  # ----------------------------------------------------------
  sds_true <- sqrt(diag(Sigma_true))
  x1_seq <- seq(mu_true[1] - grid_sd_expand * sds_true[1],
                mu_true[1] + grid_sd_expand * sds_true[1],
                length.out = grid_n)
  x2_seq <- seq(mu_true[2] - grid_sd_expand * sds_true[2],
                mu_true[2] + grid_sd_expand * sds_true[2],
                length.out = grid_n)

  grid <- expand.grid(x1_seq, x2_seq)
  grid_mat <- as.matrix(grid)
  colnames(grid_mat) <- c("x1", "x2")

  # ----------------------------------------------------------
  # Wahre und geschätzte Dichten
  # ----------------------------------------------------------
  dens_true_vec  <- true_dmvnorm_2d(grid_mat, mu = mu_true, Sigma = Sigma_true)
  dens_plain_vec <- predict_density_mv_basic_grid(grid_mat, fit_plain)
  dens_lc_vec    <- predict_density_mv_basic_grid(grid_mat, fit_lc)

  dens_true  <- normalize_density_on_grid(dens_true_vec,  x1_seq, x2_seq)
  dens_plain <- normalize_density_on_grid(dens_plain_vec, x1_seq, x2_seq)
  dens_lc    <- normalize_density_on_grid(dens_lc_vec,    x1_seq, x2_seq)

  # ----------------------------------------------------------
  # 3D-Plots
  # ----------------------------------------------------------
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mfrow = c(1, 3), mar = c(2.5, 2.5, 3, 1))

  plot_density_surface(
    dens_true,
    x1_seq,
    x2_seq,
    main = "Wahre Gaussian-Dichte"
  )

  plot_density_surface(
    dens_plain,
    x1_seq,
    x2_seq,
    main = sprintf("SM ohne Log-Konkavitaet\n(m = %d)", m)
  )

  plot_density_surface(
    dens_lc,
    x1_seq,
    x2_seq,
    main = sprintf("SM mit Log-Konkavitaet\n(m = %d, penalty = %g)", m, lc_penalty)
  )

  # ----------------------------------------------------------
  # Rückgabe
  # ----------------------------------------------------------
  invisible(list(
    x = x,
    fit_plain = fit_plain,
    fit_lc = fit_lc,
    grid = grid_mat,
    x1_seq = x1_seq,
    x2_seq = x2_seq,
    density_true = dens_true,
    density_plain = dens_plain,
    density_lc = dens_lc
  ))
}

# ------------------------------------------------------------
# (3) Beispielaufruf
# ------------------------------------------------------------
# Einfach dieses Skript sourcen; der Test läuft dann direkt los.
res_3d_test <- run_logconcave_density_3d_test(
  n = 500,
  mu_true = c(0, 0),
  Sigma_true = matrix(c(1, 0.6,
                        0.6, 1.4), 2, 2),
  m = 3,
  include_interactions = TRUE,
  standardize = TRUE,
  ridge = 1e-8,
  solver = "solve",
  log_concave_lc_method = "grid",
  lc_grid_size = 5,
  lc_max_points = 300,
  lc_box_expand = 0.25,
  lc_penalty = 1e4,
  grid_n = 80,
  grid_sd_expand = 3.5,
  seed = 123
)
