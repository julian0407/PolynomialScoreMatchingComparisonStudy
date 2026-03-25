############################################################
# Test_LogConcave_Density_3D_Full.R
#
# Separates Testskript für dein multivariates Polynomial-
# Score-Matching mit optionaler Log-Konkavität.
#
# Was das Skript macht:
# 1) simuliert 2D-Daten aus einer wählbaren wahren Dichte
#    - gaussian
#    - indep_logistic
#    - indep_gumbel
# 2) fitted dein Modell
#    - ohne log_concave
#    - mit log_concave = TRUE
# 3) rekonstruiert aus dem geschätzten Score eine 2D-Dichte
# 4) berechnet L2-Fehler gegen die wahre Dichte
# 5) macht
#    - separate 3D-Plots
#    - einen überlagerten 3D-Plot
# 6) liefert Diagnostik, wie log-konkav der unconstrained Fit
#    "von selbst" schon ist (ohne die Log-Concavity-Option)
#
# Erwartete Dateien im gleichen Ordner:
#   - helper_functions.R
#   - Multivariate_Polynomial_Score_Matching_LogConcave.R
#
# Benötigte Pakete:
#   - MASS
#   - mvtnorm
#   - plotly
############################################################

## ==========================================================
## (0) Konfiguration
## ==========================================================

SCENARIO <- "indep_logistic"      # "gaussian", "indep_logistic", "indep_gumbel"
SEED <- 124
N_TRAIN <- 10000
GRID_SIZE <- 60
M_POLY <- 3
STANDARDIZE <- TRUE
INCLUDE_INTERACTIONS <- TRUE
RIDGE <- 1e-8

# Log-concavity settings
LC_METHOD <- "grid"         # "grid" oder "data"
LC_GRID_SIZE <- 6
LC_BOX_EXPAND <- 0.15
LC_MAX_POINTS <- 400
LC_PENALTY <- 1e4
LC_TOL <- 1e-8
LC_OPTIM_METHOD <- "BFGS"

# Plot settings
MAKE_SEPARATE_PLOTS <- TRUE
MAKE_OVERLAY_PLOT <- TRUE
SEPARATE_PLOTS_PERSPECTIVE <- list(theta = 35, phi = 28, expand = 0.7)
OVERLAY_OPACITY <- 0.55

# Diagnose der "natürlichen" Log-Konkavität des unconstrained Fits
N_DIAG_POINTS <- 1500
DIAG_TOL <- 1e-8

## ==========================================================
## (1) Pakete laden
## ==========================================================

required_pkgs <- c("MASS", "mvtnorm")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Bitte installiere zuerst diese Pakete: ", paste(missing_pkgs, collapse = ", "))
}

has_plotly <- requireNamespace("plotly", quietly = TRUE)
if (!has_plotly && isTRUE(MAKE_OVERLAY_PLOT)) {
  message("Hinweis: Paket 'plotly' nicht gefunden. Der Overlay-Plot wird übersprungen.")
}

## ==========================================================
## (2) Hilfsfunktion zum robusten Sourcen
## ==========================================================

source_local_or_data <- function(filename) {
  candidates <- c(
    file.path(getwd(), filename),
    file.path("/mnt/data", filename)
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("Datei nicht gefunden: ", filename)
  source(hit, local = FALSE)
  invisible(hit)
}

source_local_or_data("helper_functions.R")
source_local_or_data("Multivariate_Pairwise_Polynomial_Score_Matching.R")

## ==========================================================
## (3) Wahre Verteilungen: Sampling, Dichte, Score
##     Score-Konvention: score(x) = - grad log f(x)
## ==========================================================

r_indep_logistic_2d <- function(n, location = c(0, 0), scale = c(1, 1)) {
  cbind(
    stats::rlogis(n, location = location[1], scale = scale[1]),
    stats::rlogis(n, location = location[2], scale = scale[2])
  )
}

r_indep_gumbel_2d <- function(n, location = c(0, 0), scale = c(1, 1)) {
  # Inverse CDF sampling: X = mu - beta * log(-log(U))
  U1 <- stats::runif(n)
  U2 <- stats::runif(n)
  cbind(
    location[1] - scale[1] * log(-log(U1)),
    location[2] - scale[2] * log(-log(U2))
  )
}

true_log_density_gaussian <- function(x, mu, Sigma) {
  mvtnorm::dmvnorm(x, mean = mu, sigma = Sigma, log = TRUE)
}

true_density_gaussian <- function(x, mu, Sigma) {
  mvtnorm::dmvnorm(x, mean = mu, sigma = Sigma, log = FALSE)
}

true_score_gaussian <- function(x, mu, Sigma) {
  x <- as.matrix(x)
  centered <- sweep(x, 2, mu, FUN = "-")
  centered %*% solve(Sigma)
}

true_log_density_indep_logistic <- function(x, location = c(0, 0), scale = c(1, 1)) {
  x <- as.matrix(x)
  sumlog <- numeric(nrow(x))
  for (j in seq_len(ncol(x))) {
    sumlog <- sumlog + stats::dlogis(x[, j], location = location[j], scale = scale[j], log = TRUE)
  }
  sumlog
}

true_density_indep_logistic <- function(x, location = c(0, 0), scale = c(1, 1)) {
  exp(true_log_density_indep_logistic(x, location = location, scale = scale))
}

true_score_indep_logistic <- function(x, location = c(0, 0), scale = c(1, 1)) {
  x <- as.matrix(x)
  out <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  for (j in seq_len(ncol(x))) {
    z <- (x[, j] - location[j]) / scale[j]
    out[, j] <- tanh(z / 2) / scale[j]
  }
  out
}

true_log_density_indep_gumbel <- function(x, location = c(0, 0), scale = c(1, 1)) {
  x <- as.matrix(x)
  out <- numeric(nrow(x))
  for (j in seq_len(ncol(x))) {
    z <- (x[, j] - location[j]) / scale[j]
    out <- out + (-log(scale[j]) - z - exp(-z))
  }
  out
}

true_density_indep_gumbel <- function(x, location = c(0, 0), scale = c(1, 1)) {
  exp(true_log_density_indep_gumbel(x, location = location, scale = scale))
}

true_score_indep_gumbel <- function(x, location = c(0, 0), scale = c(1, 1)) {
  x <- as.matrix(x)
  out <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  for (j in seq_len(ncol(x))) {
    z <- (x[, j] - location[j]) / scale[j]
    out[, j] <- (1 - exp(-z)) / scale[j]
  }
  out
}

make_scenario <- function(scenario) {
  scenario <- match.arg(scenario, c("gaussian", "indep_logistic", "indep_gumbel"))

  if (scenario == "gaussian") {
    mu <- c(0, 0)
    Sigma <- matrix(c(1.0, 0.65,
                      0.65, 1.4), 2, 2)
    return(list(
      name = "Gaussian",
      r_sample = function(n) MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma),
      true_density = function(x) true_density_gaussian(x, mu = mu, Sigma = Sigma),
      true_log_density = function(x) true_log_density_gaussian(x, mu = mu, Sigma = Sigma),
      true_score = function(x) true_score_gaussian(x, mu = mu, Sigma = Sigma)
    ))
  }

  if (scenario == "indep_logistic") {
    location <- c(0, 0)
    scale <- c(1.0, 0.7)
    return(list(
      name = "Independent logistic marginals",
      r_sample = function(n) r_indep_logistic_2d(n, location = location, scale = scale),
      true_density = function(x) true_density_indep_logistic(x, location = location, scale = scale),
      true_log_density = function(x) true_log_density_indep_logistic(x, location = location, scale = scale),
      true_score = function(x) true_score_indep_logistic(x, location = location, scale = scale)
    ))
  }

  location <- c(0, 0)
  scale <- c(0.9, 1.1)
  list(
    name = "Independent Gumbel marginals",
    r_sample = function(n) r_indep_gumbel_2d(n, location = location, scale = scale),
    true_density = function(x) true_density_indep_gumbel(x, location = location, scale = scale),
    true_log_density = function(x) true_log_density_indep_gumbel(x, location = location, scale = scale),
    true_score = function(x) true_score_indep_gumbel(x, location = location, scale = scale)
  )
}

## ==========================================================
## (4) Grid und Rekonstruktion der Dichte aus dem Score
## ==========================================================

make_2d_grid <- function(x, grid_size = 60, expand_factor = 0.20) {
  x <- as.matrix(x)
  mins <- apply(x, 2, min)
  maxs <- apply(x, 2, max)
  spans <- pmax(maxs - mins, 1e-6)
  lower <- mins - expand_factor * spans
  upper <- maxs + expand_factor * spans

  x1_seq <- seq(lower[1], upper[1], length.out = grid_size)
  x2_seq <- seq(lower[2], upper[2], length.out = grid_size)

  grid <- expand.grid(x1_seq, x2_seq)
  colnames(grid) <- c("x1", "x2")

  list(
    x1 = x1_seq,
    x2 = x2_seq,
    grid = as.matrix(grid),
    dx = x1_seq[2] - x1_seq[1],
    dy = x2_seq[2] - x2_seq[1],
    grid_size = grid_size
  )
}

score_matrix_on_grid <- function(fit, grid_obj) {
  s <- predict_score_mv_basic(grid_obj$grid, fit)
  list(
    s1 = matrix(s[, 1], nrow = grid_obj$grid_size, ncol = grid_obj$grid_size),
    s2 = matrix(s[, 2], nrow = grid_obj$grid_size, ncol = grid_obj$grid_size)
  )
}

# Rekonstruktion des Potentials V mit score = grad V.
# Pfad: zuerst entlang x1, dann entlang x2.
reconstruct_potential_from_score <- function(score_grid, grid_obj) {
  nx <- grid_obj$grid_size
  ny <- grid_obj$grid_size
  dx <- grid_obj$dx
  dy <- grid_obj$dy
  s1 <- score_grid$s1
  s2 <- score_grid$s2

  V <- matrix(0, nrow = nx, ncol = ny)

  if (nx >= 2L) {
    for (i in 2:nx) {
      V[i, 1] <- V[i - 1, 1] + 0.5 * (s1[i - 1, 1] + s1[i, 1]) * dx
    }
  }

  if (ny >= 2L) {
    for (j in 2:ny) {
      for (i in 1:nx) {
        V[i, j] <- V[i, j - 1] + 0.5 * (s2[i, j - 1] + s2[i, j]) * dy
      }
    }
  }

  V
}

normalize_log_density_on_grid <- function(logf, dx, dy) {
  f_shift <- exp(logf - max(logf))
  Z <- sum(f_shift) * dx * dy
  if (!is.finite(Z) || Z <= 0) stop("Numerische Normalisierung fehlgeschlagen.")
  f_shift / Z
}

estimate_density_from_fit <- function(fit, grid_obj) {
  score_grid <- score_matrix_on_grid(fit, grid_obj)
  V_hat <- reconstruct_potential_from_score(score_grid, grid_obj)
  logf_hat <- -V_hat
  f_hat <- normalize_log_density_on_grid(logf_hat, dx = grid_obj$dx, dy = grid_obj$dy)

  list(
    score_grid = score_grid,
    potential = V_hat,
    log_density = logf_hat,
    density = f_hat
  )
}

true_density_on_grid <- function(true_density_fun, grid_obj) {
  f <- true_density_fun(grid_obj$grid)
  matrix(f, nrow = grid_obj$grid_size, ncol = grid_obj$grid_size)
}

l2_error_grid <- function(f_hat, f_true, dx, dy) {
  sqrt(sum((f_hat - f_true)^2) * dx * dy)
}

## ==========================================================
## (5) Log-Konkavitäts-Diagnostik ohne Constraint
## ==========================================================

sample_diagnostic_points <- function(x, n_points = 1000, box_expand = 0.15, seed = NULL) {
  x <- as.matrix(x)
  if (!is.null(seed)) set.seed(seed)
  mins <- apply(x, 2, min)
  maxs <- apply(x, 2, max)
  spans <- pmax(maxs - mins, 1e-6)
  lower <- mins - box_expand * spans
  upper <- maxs + box_expand * spans

  cbind(
    stats::runif(n_points, lower[1], upper[1]),
    stats::runif(n_points, lower[2], upper[2])
  )
}

# score = grad V, daher braucht man Hess(V) >= 0 für log-concavity von f.
logconcavity_diagnostic_from_fit <- function(fit, points, tol = 1e-8) {
  points <- as.matrix(points)
  z <- if (isTRUE(fit$scaling$standardize)) {
    apply_scaling_matrix(points, fit$scaling)$z
  } else {
    points
  }

  basis_obj <- list(
    basis = fit$basis,
    basis_names = fit$basis_names,
    p = length(fit$basis),
    d = fit$d,
    m = fit$m,
    include_interactions = fit$include_interactions
  )

  hess_design <- build_hessian_design_mv_basic(z, basis_obj)
  H_arr <- hessian_array_from_theta_mv_basic(fit$theta, hess_design, d = fit$d)
  
  get_min_eigs_from_H_array <- function(H_arr) {
    dm <- dim(H_arr)
    
    if (length(dm) != 3L) {
      stop("H_arr must be a 3D array.")
    }
    
    if (dm[1] == dm[2]) {
      min_eigs <- vapply(seq_len(dm[3]), function(i) {
        H <- H_arr[, , i]
        min(eigen(H, symmetric = TRUE, only.values = TRUE)$values)
      }, numeric(1))
      return(min_eigs)
    }
    
    if (dm[2] == dm[3]) {
      H_arr <- aperm(H_arr, c(2, 3, 1))
      dm <- dim(H_arr)
      
      min_eigs <- vapply(seq_len(dm[3]), function(i) {
        H <- H_arr[, , i]
        min(eigen(H, symmetric = TRUE, only.values = TRUE)$values)
      }, numeric(1))
      return(min_eigs)
    }
    
    stop("H_arr has incompatible dimensions.")
  }
  min_eigs <- get_min_eigs_from_H_array(H_arr)
  
  violations <- pmax(0, -min_eigs)

  list(
    n_points = ncol(H_arr),
    min_eig_min = min(min_eigs),
    min_eig_median = stats::median(min_eigs),
    min_eig_mean = mean(min_eigs),
    prop_violating = mean(min_eigs < -tol),
    max_violation = max(violations),
    mean_violation = mean(violations),
    tol = tol
  )
}

## ==========================================================
## (6) Plot-Funktionen
## ==========================================================

plot_surface_base <- function(x1, x2, z, main = "", col = "lightblue", zlim = NULL) {
  do.call(
    graphics::persp,
    c(
      list(
        x = x1,
        y = x2,
        z = z,
        main = main,
        xlab = "x1",
        ylab = "x2",
        zlab = "density",
        col = col,
        border = NA,
        shade = 0.35,
        zlim = zlim,
        ticktype = "detailed"
      ),
      SEPARATE_PLOTS_PERSPECTIVE
    )
  )
}

plot_overlay_plotly <- function(x1, x2, z_true, z_plain, z_lc, title = "Overlay 3D") {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    message("plotly nicht installiert - Overlay wird übersprungen.")
    return(invisible(NULL))
  }

  p <- plotly::plot_ly()
  p <- plotly::add_surface(
    p,
    x = x1, y = x2, z = z_true,
    opacity = OVERLAY_OPACITY,
    showscale = FALSE,
    colorscale = list(c(0, 1), c("blue", "blue")),
    name = "True"
  )
  p <- plotly::add_surface(
    p,
    x = x1, y = x2, z = z_plain,
    opacity = OVERLAY_OPACITY,
    showscale = FALSE,
    colorscale = list(c(0, 1), c("red", "red")),
    name = "SM"
  )
  p <- plotly::add_surface(
    p,
    x = x1, y = x2, z = z_lc,
    opacity = OVERLAY_OPACITY,
    showscale = FALSE,
    colorscale = list(c(0, 1), c("green", "green")),
    name = "SM + log-concave"
  )
  p <- plotly::layout(
    p,
    title = title,
    scene = list(
      xaxis = list(title = "x1"),
      yaxis = list(title = "x2"),
      zaxis = list(title = "density")
    )
  )
  p
}

## ==========================================================
## (7) Hauptlauf
## ==========================================================

run_density_test <- function(scenario = SCENARIO) {
  set.seed(SEED)

  sc <- make_scenario(scenario)
  cat("\n====================================================\n")
  cat("Szenario:", sc$name, "\n")
  cat("====================================================\n")

  x_train <- sc$r_sample(N_TRAIN)
  x_train <- as.matrix(x_train)

  fit_plain <- fit_score_matching_mv_basic(
    x_train,
    m = M_POLY,
    include_interactions = INCLUDE_INTERACTIONS,
    standardize = STANDARDIZE,
    ridge = RIDGE,
    log_concave = FALSE
  )

  fit_lc <- fit_score_matching_mv_basic(
    x_train,
    m = M_POLY,
    include_interactions = INCLUDE_INTERACTIONS,
    standardize = STANDARDIZE,
    ridge = RIDGE,
    log_concave = TRUE,
    lc_method = LC_METHOD,
    lc_grid_size = LC_GRID_SIZE,
    lc_box_expand = LC_BOX_EXPAND,
    lc_max_points = LC_MAX_POINTS,
    lc_penalty = LC_PENALTY,
    lc_tol = LC_TOL,
    lc_optim_method = LC_OPTIM_METHOD
  )

  grid_obj <- make_2d_grid(x_train, grid_size = GRID_SIZE, expand_factor = 0.20)

  est_plain <- estimate_density_from_fit(fit_plain, grid_obj)
  est_lc <- estimate_density_from_fit(fit_lc, grid_obj)
  f_true <- true_density_on_grid(sc$true_density, grid_obj)

  l2_plain <- l2_error_grid(est_plain$density, f_true, grid_obj$dx, grid_obj$dy)
  l2_lc <- l2_error_grid(est_lc$density, f_true, grid_obj$dx, grid_obj$dy)

  # Post-hoc Diagnose: wie log-konkav ist der unkonstruierte Fit bereits?
  diag_points <- sample_diagnostic_points(
    x_train,
    n_points = N_DIAG_POINTS,
    box_expand = 0.20,
    seed = SEED + 1L
  )
  lc_diag_plain <- logconcavity_diagnostic_from_fit(fit_plain, diag_points, tol = DIAG_TOL)
  lc_diag_lc <- logconcavity_diagnostic_from_fit(fit_lc, diag_points, tol = DIAG_TOL)

  cat("\n--- L2-Fehler gegen wahre Dichte ---\n")
  cat(sprintf("SM ohne log-concave     : %.6f\n", l2_plain))
  cat(sprintf("SM mit log-concave      : %.6f\n", l2_lc))

  cat("\n--- Natürliche Log-Konkavitäts-Diagnostik (ohne Grid-Constraint im Test) ---\n")
  cat("Unconstrained Fit:\n")
  cat(sprintf("  min(min_eig Hess(V))) : %.6f\n", lc_diag_plain$min_eig_min))
  cat(sprintf("  prop violating        : %.4f\n", lc_diag_plain$prop_violating))
  cat(sprintf("  max violation         : %.6f\n", lc_diag_plain$max_violation))

  cat("Constrained Fit:\n")
  cat(sprintf("  min(min_eig Hess(V))) : %.6f\n", lc_diag_lc$min_eig_min))
  cat(sprintf("  prop violating        : %.4f\n", lc_diag_lc$prop_violating))
  cat(sprintf("  max violation         : %.6f\n", lc_diag_lc$max_violation))

  if (!is.null(fit_lc$lc_diagnostics)) {
    cat("\n--- In-fit Log-Concavity-Diagnostik des constrained Fits ---\n")
    print(fit_lc$lc_diagnostics)
  }

  zlim_all <- range(c(f_true, est_plain$density, est_lc$density), finite = TRUE)

  if (isTRUE(MAKE_SEPARATE_PLOTS)) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(1, 3))

    plot_surface_base(
      grid_obj$x1, grid_obj$x2, f_true,
      main = paste0("True density\n", sc$name),
      col = "lightblue",
      zlim = zlim_all
    )
    plot_surface_base(
      grid_obj$x1, grid_obj$x2, est_plain$density,
      main = paste0("Estimated density\nSM, L2 = ", format(round(l2_plain, 5), nsmall = 5)),
      col = "mistyrose",
      zlim = zlim_all
    )
    plot_surface_base(
      grid_obj$x1, grid_obj$x2, est_lc$density,
      main = paste0("Estimated density\nSM + LC, L2 = ", format(round(l2_lc, 5), nsmall = 5)),
      col = "palegreen3",
      zlim = zlim_all
    )
  }

  overlay_plot <- NULL
  if (isTRUE(MAKE_OVERLAY_PLOT) && requireNamespace("plotly", quietly = TRUE)) {
    overlay_plot <- plot_overlay_plotly(
      grid_obj$x1,
      grid_obj$x2,
      f_true,
      est_plain$density,
      est_lc$density,
      title = paste0("Overlay 3D: ", sc$name,
                     " | L2 plain = ", round(l2_plain, 5),
                     " | L2 LC = ", round(l2_lc, 5))
    )
    print(overlay_plot)
  }

  invisible(list(
    scenario = sc$name,
    x_train = x_train,
    fit_plain = fit_plain,
    fit_lc = fit_lc,
    grid = grid_obj,
    true_density = f_true,
    est_plain = est_plain,
    est_lc = est_lc,
    l2_plain = l2_plain,
    l2_lc = l2_lc,
    lc_diag_plain = lc_diag_plain,
    lc_diag_lc = lc_diag_lc,
    overlay_plot = overlay_plot
  ))
}

## ==========================================================
## (8) Optional: mehrere Szenarien nacheinander laufen lassen
## ==========================================================

run_all_default_scenarios <- function() {
  scenarios <- c("gaussian", "indep_logistic", "indep_gumbel")
  out <- vector("list", length(scenarios))
  names(out) <- scenarios
  for (s in scenarios) {
    out[[s]] <- run_density_test(s)
  }
  invisible(out)
}

## ==========================================================
## (9) Direkter Start
## ==========================================================

result <- run_density_test(SCENARIO)

# Für weitere Tests z.B.:
# result_logistic <- run_density_test("indep_logistic")
# result_gumbel   <- run_density_test("indep_gumbel")
# all_results     <- run_all_default_scenarios()

