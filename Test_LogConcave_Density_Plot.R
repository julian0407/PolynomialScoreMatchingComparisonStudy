# ============================================================
# Test_LogConcave_Density_Plot.R
#
# Separate test script for visualizing the fitted multivariate
# density from the polynomial score matching model and comparing
# it against the true 2D Gaussian density.
#
# Intended use:
#   - source helper_functions.R
#   - source Multivariate_Polynomial_Score_Matching_LogConcave.R
#   - simulate or provide 2D Gaussian data
#   - fit unconstrained and log-concave versions
#   - plot true density vs estimated density
#
# IMPORTANT:
#   This script is written for d = 2, because a density surface
#   plot is otherwise not directly visualizable.
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# (1) Load code
# ------------------------------------------------------------
source("helper_functions.R")
source("Multivariate_Polynomial_Score_Matching_LogConcave_2.R")

# ------------------------------------------------------------
# (2) Helper: simulate 2D Gaussian data without extra packages
# ------------------------------------------------------------
rmvnorm_basic_2d <- function(n, mu, Sigma, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)
  if (length(mu) != 2L) stop("mu must have length 2.")
  if (!all(dim(Sigma) == c(2L, 2L))) stop("Sigma must be 2x2.")

  L <- chol(Sigma)
  z <- matrix(stats::rnorm(n * 2L), ncol = 2L)
  x <- z %*% L
  x <- sweep(x, 2, mu, FUN = "+")
  colnames(x) <- c("x1", "x2")
  x
}

# ------------------------------------------------------------
# (3) Helper: evaluate basis matrix itself (not derivative)
# ------------------------------------------------------------
evaluate_basis_mv_basic <- function(z, basis_obj) {
  z <- as_obs_matrix(z)
  if (ncol(z) != basis_obj$d) stop("Dimension mismatch in evaluate_basis_mv_basic().")

  n <- nrow(z)
  p <- basis_obj$p
  P <- build_power_cache_mv_basic(z, max_deg = basis_obj$m)
  B <- matrix(0, nrow = n, ncol = p)

  for (col in seq_len(p)) {
    b <- basis_obj$basis[[col]]

    if (b$type == "uni") {
      B[, col] <- get_power_mv_basic(P, b$j, b$r)
    } else if (b$type == "pair") {
      B[, col] <- get_power_mv_basic(P, b$j, b$r) * get_power_mv_basic(P, b$k, b$s)
    } else {
      stop("Unknown basis type.")
    }
  }

  colnames(B) <- basis_obj$basis_names
  B
}

# ------------------------------------------------------------
# (4) Helper: unnormalized log density under fitted model
#     The fitted score is r(x) = - grad log f(x) = grad psi(x)
#     with potential psi(x) = sum_j theta_j b_j(x).
#     Hence: log f(x) = const - psi(x).
# ------------------------------------------------------------
predict_log_unnormalized_density_mv_basic <- function(newx, fit) {
  if (!inherits(fit, "score_matching_mv_basic_fit")) {
    stop("fit must be of class 'score_matching_mv_basic_fit'.")
  }

  newx <- as_obs_matrix(newx)
  if (ncol(newx) != fit$d) stop("Dimension mismatch.")

  if (isTRUE(fit$scaling$standardize)) {
    z <- apply_scaling_matrix(newx, fit$scaling)$z
  } else {
    z <- newx
  }

  basis_obj <- list(
    basis = fit$basis,
    basis_names = fit$basis_names,
    p = length(fit$basis),
    d = fit$d,
    m = fit$m,
    include_interactions = fit$include_interactions
  )

  B <- evaluate_basis_mv_basic(z, basis_obj)
  psi <- as.numeric(B %*% fit$theta)

  # Jacobian term from standardization is constant in x and therefore
  # irrelevant for the shape. It is ignored here and absorbed into the
  # later numerical normalization on the plotting grid.
  -psi
}

# ------------------------------------------------------------
# (5) Helper: true 2D Gaussian density
# ------------------------------------------------------------
dmvnorm_basic_2d <- function(x, mu, Sigma, log = FALSE) {
  x <- as_obs_matrix(x)
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)

  if (ncol(x) != 2L) stop("x must be n x 2.")
  if (length(mu) != 2L) stop("mu must have length 2.")
  if (!all(dim(Sigma) == c(2L, 2L))) stop("Sigma must be 2x2.")

  Sinv <- solve(Sigma)
  logdet <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)

  xc <- sweep(x, 2, mu, FUN = "-")
  quad <- rowSums((xc %*% Sinv) * xc)

  logdens <- -0.5 * (2L * log(2 * pi) + logdet + quad)
  if (isTRUE(log)) logdens else exp(logdens)
}

# ------------------------------------------------------------
# (6) Helper: create 2D grid and normalize density numerically
# ------------------------------------------------------------
make_grid_2d <- function(x, n_grid = 120L, expand = 0.15) {
  x <- as_obs_matrix(x)
  if (ncol(x) != 2L) stop("This plotting script is restricted to 2D data.")

  mins <- apply(x, 2, min)
  maxs <- apply(x, 2, max)
  spans <- pmax(maxs - mins, 1e-6)

  mins <- mins - expand * spans
  maxs <- maxs + expand * spans

  x1 <- seq(mins[1], maxs[1], length.out = n_grid)
  x2 <- seq(mins[2], maxs[2], length.out = n_grid)
  g <- expand.grid(x1 = x1, x2 = x2)

  list(
    x1 = x1,
    x2 = x2,
    grid = as.matrix(g),
    dx = x1[2] - x1[1],
    dy = x2[2] - x2[1]
  )
}

normalize_from_log_grid <- function(logu, dx, dy) {
  m <- max(logu)
  mass <- sum(exp(logu - m)) * dx * dy
  logZ <- m + log(mass)
  exp(logu - logZ)
}

matrix_from_grid_values <- function(vals, n1, n2) {
  matrix(vals, nrow = n2, ncol = n1, byrow = FALSE)
}

# ------------------------------------------------------------
# (7) Build density surfaces for true / unconstrained / log-concave
# ------------------------------------------------------------
build_density_surfaces <- function(x,
                                   fit_unconstrained,
                                   fit_logconcave,
                                   mu_true,
                                   Sigma_true,
                                   n_grid = 120L,
                                   expand = 0.15) {
  x <- as_obs_matrix(x)
  if (ncol(x) != 2L) stop("This plotting script requires d = 2.")

  gd <- make_grid_2d(x, n_grid = n_grid, expand = expand)
  grid <- gd$grid

  true_vals <- dmvnorm_basic_2d(grid, mu = mu_true, Sigma = Sigma_true, log = FALSE)

  logu_uncon <- predict_log_unnormalized_density_mv_basic(grid, fit_unconstrained)
  dens_uncon <- normalize_from_log_grid(logu_uncon, dx = gd$dx, dy = gd$dy)

  logu_lc <- predict_log_unnormalized_density_mv_basic(grid, fit_logconcave)
  dens_lc <- normalize_from_log_grid(logu_lc, dx = gd$dx, dy = gd$dy)

  list(
    x1 = gd$x1,
    x2 = gd$x2,
    true = matrix_from_grid_values(true_vals, length(gd$x1), length(gd$x2)),
    unconstrained = matrix_from_grid_values(dens_uncon, length(gd$x1), length(gd$x2)),
    logconcave = matrix_from_grid_values(dens_lc, length(gd$x1), length(gd$x2)),
    grid = grid,
    dx = gd$dx,
    dy = gd$dy
  )
}

# ------------------------------------------------------------
# (8) Plot helper
# ------------------------------------------------------------
plot_density_comparison_2d <- function(surfaces,
                                       x,
                                       main_prefix = "2D density comparison",
                                       draw_points = TRUE) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  x <- as_obs_matrix(x)
  zlim <- range(c(surfaces$true, surfaces$unconstrained, surfaces$logconcave))

  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

  image(
    x = surfaces$x1,
    y = surfaces$x2,
    z = surfaces$true,
    xlab = "x1",
    ylab = "x2",
    main = paste(main_prefix, "\nTrue Gaussian"),
    zlim = zlim
  )
  contour(surfaces$x1, surfaces$x2, surfaces$true, add = TRUE)
  if (isTRUE(draw_points)) points(x[, 1], x[, 2], pch = 16, cex = 0.35)

  image(
    x = surfaces$x1,
    y = surfaces$x2,
    z = surfaces$unconstrained,
    xlab = "x1",
    ylab = "x2",
    main = paste(main_prefix, "\nPolynomial SM"),
    zlim = zlim
  )
  contour(surfaces$x1, surfaces$x2, surfaces$unconstrained, add = TRUE)
  if (isTRUE(draw_points)) points(x[, 1], x[, 2], pch = 16, cex = 0.35)

  image(
    x = surfaces$x1,
    y = surfaces$x2,
    z = surfaces$logconcave,
    xlab = "x1",
    ylab = "x2",
    main = paste(main_prefix, "\nPolynomial SM + log-concave"),
    zlim = zlim
  )
  contour(surfaces$x1, surfaces$x2, surfaces$logconcave, add = TRUE)
  if (isTRUE(draw_points)) points(x[, 1], x[, 2], pch = 16, cex = 0.35)
}

# ------------------------------------------------------------
# (9) Integrated test function
# ------------------------------------------------------------
run_logconcave_density_test <- function(n = 1000L,
                                        mu_true = c(0, 0),
                                        Sigma_true = matrix(c(1.0, 0.5,
                                                              0.5, 1.5),
                                                            nrow = 2,
                                                            byrow = TRUE),
                                        m = 3,
                                        include_interactions = TRUE,
                                        ridge = 1e-8,
                                        standardize = TRUE,
                                        lc_method = "grid",
                                        lc_grid_size = 6L,
                                        lc_max_points = 400L,
                                        lc_penalty = 1e4,
                                        seed = 123,
                                        n_grid = 120L,
                                        expand = 0.15,
                                        save_pdf = FALSE,
                                        pdf_file = "logconcave_density_comparison.pdf") {
  x <- rmvnorm_basic_2d(n = n, mu = mu_true, Sigma = Sigma_true, seed = seed)

  fit_unconstrained <- fit_score_matching_mv_basic(
    x,
    m = m,
    include_interactions = include_interactions,
    standardize = standardize,
    ridge = ridge,
    log_concave = FALSE
  )

  fit_logconcave <- fit_score_matching_mv_basic(
    x,
    m = m,
    include_interactions = include_interactions,
    standardize = standardize,
    ridge = ridge,
    log_concave = TRUE,
    lc_method = lc_method,
    lc_grid_size = lc_grid_size,
    lc_max_points = lc_max_points,
    lc_penalty = lc_penalty,
    lc_seed = seed
  )

  surfaces <- build_density_surfaces(
    x = x,
    fit_unconstrained = fit_unconstrained,
    fit_logconcave = fit_logconcave,
    mu_true = mu_true,
    Sigma_true = Sigma_true,
    n_grid = n_grid,
    expand = expand
  )

  if (isTRUE(save_pdf)) {
    grDevices::pdf(pdf_file, width = 14, height = 4.8)
    plot_density_comparison_2d(surfaces, x)
    grDevices::dev.off()
  }

  plot_density_comparison_2d(surfaces, x)

  cat("\n--- Log-concavity diagnostics ---\n")
  if (!is.null(fit_logconcave$lc_diagnostics)) {
    print(fit_logconcave$lc_diagnostics)
  } else {
    cat("No log-concavity diagnostics available.\n")
  }

  invisible(list(
    x = x,
    fit_unconstrained = fit_unconstrained,
    fit_logconcave = fit_logconcave,
    surfaces = surfaces
  ))
}

# ------------------------------------------------------------
# (10) Example run
# ------------------------------------------------------------
# Adjust these settings as needed. This block can stay at the end,
# so the script works directly as a standalone test file.

res <- run_logconcave_density_test(
  n = 1000,
  mu_true = c(0.5, -0.25),
  Sigma_true = matrix(c(1.0, 0.6,
                        0.6, 1.8),
                      nrow = 2,
                      byrow = TRUE),
  m = 3,
  include_interactions = TRUE,
  ridge = 1e-8,
  standardize = TRUE,
  lc_method = "grid",
  lc_grid_size = 6,
  lc_max_points = 400,
  lc_penalty = 1e4,
  seed = 123,
  n_grid = 130,
  expand = 0.20,
  save_pdf = FALSE
)
