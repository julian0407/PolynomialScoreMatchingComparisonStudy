# ============================================================
# conditional_mean_regression_logconcave_test.R
#
# A simple downstream regression test for multivariate density
# estimators in the log-concave setting.
#
# Main idea:
#   We fit a JOINT density estimator for (X, Y) and then derive
#   the conditional mean m(x) = E[Y | X = x] from the fitted
#   conditional density.
#
# Why this is a sensible test:
#   - It is a classical regression functional derived from the
#     conditional density p(y | x).
#   - It is simple and interpretable: we compare how well KDE,
#     multivariate log-concave MLE, and score matching recover
#     the regression function.
#   - For score matching, we only need an UNNORMALIZED conditional
#     density in y for each fixed x. The x-dependent normalizing
#     constant cancels in the conditional mean ratio.
#
# Why this is especially appropriate here:
#   - Hyndman, Bashtannyk and Grunwald (1996) explicitly connect
#     conditional density estimation and regression, noting that
#     the conditional mean implied by a conditional density estimator
#     is itself a regression smoother.
#   - Sugiyama et al. (2010) motivate conditional density estimation
#     as a regression generalization when one wants more than just
#     the conditional mean.
#   - In our benchmark we deliberately go in the opposite direction:
#     we use a classical regression target (the conditional mean)
#     as a simple downstream quantity that ALL methods can be judged on.
#
# Design choice for the data-generating process:
#   We use a bivariate Gaussian model. This is intentionally simple:
#   - the JOINT density is log-concave,
#   - the conditional density Y|X=x is Gaussian and therefore log-concave,
#   - the true conditional mean is known in closed form.
#
# Interpretation:
#   This script is best viewed as a SANITY CHECK / baseline test in
#   a clean log-concave model, not as the final difficult benchmark.
# ============================================================

# ------------------------------------------------------------
# (0) Load your existing code
# ------------------------------------------------------------
source("helper_functions.R")
source("KDE.R")
source("LogConcaveMLE.R")
source("Multivariate_Pairwise_Polynomial_Score_Matching.R")


# ------------------------------------------------------------
# (1) Small numerical helper functions
# ------------------------------------------------------------

# Trapezoidal integration on a finite grid.
# We use it repeatedly for numerical expectations over a y-grid.
trapz_1d_local <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

# Cumulative trapezoidal integration.
# If dy/dt = g(t), then cumtrapz gives an approximation to
# integral_{t0}^t g(u) du along the supplied grid.
#
# We need this for score matching because SM gives us the score
# component r_y(x, y) = - d/dy log f(x, y). Integrating -r_y over y
# reconstructs log f(x, y) up to an additive constant depending on x.
cumtrapz_1d_local <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  n <- length(x)
  if (n != length(y)) stop("x and y must have the same length.")
  if (n < 2L) return(rep(0, n))

  out <- numeric(n)
  for (i in 2:n) {
    out[i] <- out[i - 1L] + 0.5 * (x[i] - x[i - 1L]) * (y[i] + y[i - 1L])
  }
  out
}

# Stable normalization helper.
# Given log-weights on a grid, subtract the maximum first so that
# exponentiation is numerically stable.
normalize_from_log_grid <- function(grid, logw) {
  grid <- as.numeric(grid)
  logw <- as.numeric(logw)
  keep <- is.finite(grid) & is.finite(logw)
  grid <- grid[keep]
  logw <- logw[keep]
  if (length(grid) < 2L) {
    return(list(grid = grid, dens = rep(NA_real_, length(grid)), mass = NA_real_))
  }

  shift <- max(logw)
  w <- exp(logw - shift)
  mass <- trapz_1d_local(grid, w)
  if (!is.finite(mass) || mass <= 0) {
    return(list(grid = grid, dens = rep(NA_real_, length(grid)), mass = NA_real_))
  }

  list(grid = grid, dens = w / mass, mass = mass * exp(shift), shift = shift)
}


# ------------------------------------------------------------
# (2) Data-generating process: simple log-concave regression model
# ------------------------------------------------------------

# We simulate a bivariate Gaussian vector (X, Y).
# This gives a JOINT log-concave density and the classical conditional
# mean formula
#
#   E[Y | X = x] = mu_y + rho * (sigma_y / sigma_x) * (x - mu_x).
#
# Here we choose mean zero and sigma_x = sigma_y = 1 for simplicity,
# so the truth becomes
#
#   E[Y | X = x] = rho * x.
#
# This is deliberately simple because it isolates the density estimation
# step from unnecessary modeling complexity.
r_joint_bvn <- function(n,
                        rho = 0.7,
                        mu = c(0, 0),
                        sigma_x = 1,
                        sigma_y = 1) {
  if (abs(rho) >= 1) stop("rho must satisfy |rho| < 1.")
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  x <- mu[1] + sigma_x * z1
  y <- mu[2] + sigma_y * (rho * z1 + sqrt(1 - rho^2) * z2)
  cbind(x, y)
}

true_cond_mean_bvn <- function(x,
                               rho = 0.7,
                               mu = c(0, 0),
                               sigma_x = 1,
                               sigma_y = 1) {
  mu[2] + rho * (sigma_y / sigma_x) * (x - mu[1])
}


# ------------------------------------------------------------
# (3) Recover E[Y|X=x] from a fitted model
# ------------------------------------------------------------

# Generic conditional mean from a JOINT density evaluator.
#
# For a fixed x0 we approximate
#
#   E[Y | X=x0] = int y f(x0, y) dy / int f(x0, y) dy.
#
# IMPORTANT:
# The x-specific factor p(x0) is irrelevant because it cancels in the ratio.
# Therefore an unnormalized conditional density in y is enough.
conditional_mean_from_logjoint <- function(x0,
                                           y_grid,
                                           predict_logjoint_fun) {
  eval_grid <- cbind(rep(x0, length(y_grid)), y_grid)
  logjoint <- as.numeric(predict_logjoint_fun(eval_grid))

  norm_obj <- normalize_from_log_grid(y_grid, logjoint)
  dens <- norm_obj$dens
  if (all(!is.finite(dens))) return(NA_real_)

  num <- trapz_1d_local(norm_obj$grid, norm_obj$grid * dens)
  den <- trapz_1d_local(norm_obj$grid, dens)
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

# Conditional mean for multivariate score matching.
#
# Your score functions return
#   r(x, y) = - grad log f(x, y).
#
# Hence the second component satisfies
#   r_y(x, y) = - d/dy log f(x, y).
#
# Therefore, for fixed x,
#   log f(x, y) = C(x) - integral r_y(x, t) dt.
#
# The unknown C(x) depends only on x and cancels in the ratio for E[Y|X=x].
# So we only need to reconstruct the y-shape of the conditional density.
conditional_mean_from_sm_score <- function(x0,
                                           y_grid,
                                           fit_sm) {
  eval_grid <- cbind(rep(x0, length(y_grid)), y_grid)
  score_mat <- predict_score_mv_basic(eval_grid, fit_sm)

  # second component corresponds to derivative with respect to y
  r_y <- as.numeric(score_mat[, 2])
  if (length(r_y) != length(y_grid)) return(NA_real_)
  if (any(!is.finite(r_y))) return(NA_real_)

  # log f(x, y) up to an additive constant in x
  log_unnorm <- -cumtrapz_1d_local(y_grid, r_y)

  norm_obj <- normalize_from_log_grid(y_grid, log_unnorm)
  dens <- norm_obj$dens
  if (all(!is.finite(dens))) return(NA_real_)

  num <- trapz_1d_local(norm_obj$grid, norm_obj$grid * dens)
  den <- trapz_1d_local(norm_obj$grid, dens)
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(NA_real_)
  num / den
}


# ------------------------------------------------------------
# (4) Fit the three competing estimators on joint data (X, Y)
# ------------------------------------------------------------

fit_all_models <- function(train_xy) {
  list(
    KDE = fit_kde_mv(train_xy, H_method = "Hpi", diagonal = FALSE),
    MLE = fit_logconcave_mle_mv(train_xy, smoothed = FALSE),
    SM  = fit_score_matching_mv_basic(
      train_xy,
      m = 2,
      include_interactions = TRUE,
      standardize = TRUE,
      ridge = 1e-6,
      log_concave = TRUE,
      lc_method = "m2"
    )
  )
}


# ------------------------------------------------------------
# (5) Evaluate conditional mean regression on a test set
# ------------------------------------------------------------

# We use two error notions:
#
# (A) mean_function_mse:
#     squared error against the TRUE conditional mean E[Y|X=x].
#     This is the cleanest estimator-comparison metric in simulation,
#     because it measures only how well the method recovers the regression
#     functional implied by the density.
#
# (B) prediction_mse:
#     squared error against the realized Y.
#     This looks more like classical regression evaluation, but it includes
#     irreducible noise and is therefore less targeted.
#
# Both are reported. If you only keep one number in the thesis, I would
# use mean_function_mse as the main metric for this synthetic experiment.
evaluate_conditional_mean_regression <- function(fits,
                                                 x_test,
                                                 y_test,
                                                 true_mean_test,
                                                 y_grid) {
  pred_kde <- vapply(
    x_test,
    FUN = function(x0) {
      conditional_mean_from_logjoint(
        x0 = x0,
        y_grid = y_grid,
        predict_logjoint_fun = function(newxy) predict_logdensity_kde_mv(newxy, fits$KDE)
      )
    },
    FUN.VALUE = numeric(1)
  )

  pred_mle <- vapply(
    x_test,
    FUN = function(x0) {
      conditional_mean_from_logjoint(
        x0 = x0,
        y_grid = y_grid,
        predict_logjoint_fun = function(newxy) predict_logdensity_logconcave_mv(newxy, fits$MLE)
      )
    },
    FUN.VALUE = numeric(1)
  )

  pred_sm <- vapply(
    x_test,
    FUN = function(x0) conditional_mean_from_sm_score(x0, y_grid, fits$SM),
    FUN.VALUE = numeric(1)
  )

  out <- data.frame(
    x = x_test,
    y = y_test,
    true_mean = true_mean_test,
    pred_KDE = pred_kde,
    pred_MLE = pred_mle,
    pred_SM  = pred_sm
  )

  metric_row <- function(pred) {
    keep <- is.finite(pred) & is.finite(true_mean_test) & is.finite(y_test)
    if (!any(keep)) {
      return(c(mean_function_mse = NA_real_, prediction_mse = NA_real_, n_valid = 0))
    }
    c(
      mean_function_mse = mean((pred[keep] - true_mean_test[keep])^2),
      prediction_mse    = mean((pred[keep] - y_test[keep])^2),
      n_valid           = sum(keep)
    )
  }

  metrics <- rbind(
    KDE = metric_row(pred_kde),
    MLE = metric_row(pred_mle),
    SM  = metric_row(pred_sm)
  )

  list(predictions = out, metrics = as.data.frame(metrics))
}


# ------------------------------------------------------------
# (6) One complete experiment
# ------------------------------------------------------------

run_one_regression_experiment <- function(n_train = 400,
                                          n_test = 300,
                                          rho = 0.7,
                                          seed = 123) {
  set.seed(seed)

  # Training data for fitting the JOINT density of (X, Y)
  train_xy <- r_joint_bvn(n_train, rho = rho)

  # Independent test data
  test_xy <- r_joint_bvn(n_test, rho = rho)
  x_test <- test_xy[, 1]
  y_test <- test_xy[, 2]
  true_mean_test <- true_cond_mean_bvn(x_test, rho = rho)

  # We need a finite y-grid for the numerical conditional expectation.
  # We choose a wide grid around the observed y-range.
  # This is sufficient here because the Gaussian tails decay quickly.
  y_all <- c(train_xy[, 2], test_xy[, 2])
  y_sd <- stats::sd(y_all)
  if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1
  y_grid <- seq(min(y_all) - 4 * y_sd,
                max(y_all) + 4 * y_sd,
                length.out = 601)

  fits <- fit_all_models(train_xy)
  eval <- evaluate_conditional_mean_regression(
    fits = fits,
    x_test = x_test,
    y_test = y_test,
    true_mean_test = true_mean_test,
    y_grid = y_grid
  )

  list(
    train_xy = train_xy,
    test_xy = test_xy,
    y_grid = y_grid,
    fits = fits,
    predictions = eval$predictions,
    metrics = eval$metrics
  )
}


# ------------------------------------------------------------
# (7) Replicated benchmark
# ------------------------------------------------------------

# Repeating the experiment several times gives a more stable picture.
# This is useful because KDE and MLE can vary across random samples,
# and the score-based reconstruction for SM also has numerical variation.
run_replicated_regression_benchmark <- function(n_rep = 20,
                                                n_train = 400,
                                                n_test = 300,
                                                rho = 0.7,
                                                seed = 123) {
  all_rows <- vector("list", n_rep)

  for (r in seq_len(n_rep)) {
    one <- run_one_regression_experiment(
      n_train = n_train,
      n_test = n_test,
      rho = rho,
      seed = seed + r - 1L
    )

    tab <- one$metrics
    tab$method <- rownames(tab)
    rownames(tab) <- NULL
    tab$rep <- r
    all_rows[[r]] <- tab[, c("rep", "method", "mean_function_mse", "prediction_mse", "n_valid")]
  }

  res <- do.call(rbind, all_rows)

  aggregate_table <- aggregate(
    cbind(mean_function_mse, prediction_mse) ~ method,
    data = res,
    FUN = function(z) c(mean = mean(z, na.rm = TRUE), sd = stats::sd(z, na.rm = TRUE))
  )

  list(raw = res, summary = aggregate_table)
}


# ------------------------------------------------------------
# (8) Minimal example run
# ------------------------------------------------------------

# A single run (useful for inspecting predictions visually)
example_run <- run_one_regression_experiment(
  n_train = 400,
  n_test = 300,
  rho = 0.7,
  seed = 1
)

print("Single-run metrics:")
print(example_run$metrics)

# A replicated benchmark (useful for thesis tables)
benchmark_res <- run_replicated_regression_benchmark(
  n_rep = 20,
  n_train = 400,
  n_test = 300,
  rho = 0.7,
  seed = 100
)

print("Replicated benchmark: raw results")
print(head(benchmark_res$raw))

print("Replicated benchmark: summary")
print(benchmark_res$summary)


# ------------------------------------------------------------
# (9) Optional plots
# ------------------------------------------------------------

# Simple plot of predicted conditional means against the true mean curve.
# This helps you understand WHETHER the method gets the regression shape right,
# not only whether one summary number is small.
plot_conditional_mean_curves <- function(example_run) {
  pred_df <- example_run$predictions
  ord <- order(pred_df$x)
  pred_df <- pred_df[ord, ]

  plot(pred_df$x, pred_df$true_mean,
       type = "l", lwd = 3,
       xlab = "x", ylab = "conditional mean",
       main = "Conditional mean recovered from fitted joint density")
  lines(pred_df$x, pred_df$pred_KDE, lwd = 2, lty = 2)
  lines(pred_df$x, pred_df$pred_MLE, lwd = 2, lty = 3)
  lines(pred_df$x, pred_df$pred_SM,  lwd = 2, lty = 4)
  legend("topleft",
         legend = c("Truth", "KDE", "Log-concave MLE", "Score Matching"),
         lwd = c(3, 2, 2, 2),
         lty = c(1, 2, 3, 4),
         bty = "n")
}

# Uncomment if you want the plot immediately:
# plot_conditional_mean_curves(example_run)


# ------------------------------------------------------------
# (10) Final interpretation notes
# ------------------------------------------------------------
#
# What this script tells you:
#   - It compares the methods on a classical downstream target from
#     conditional density estimation: the conditional mean.
#   - It avoids the unfairness of comparing only score loss, while still
#     allowing score matching to participate.
#   - It uses only information that is naturally available from each method:
#       * KDE / MLE: explicit log joint density
#       * SM: score field, integrated only in the response direction y
#
# What this script does NOT tell you:
#   - It is not a full conditional density benchmark.
#   - It does not prove that one method is globally better as a density estimator.
#   - In the Gaussian case, SM with quadratic basis is close to a well-specified model,
#     so this experiment should be viewed as a clean baseline rather than a hard stress test.
#
# If you later want a harder but still plausible test, the next step would be:
#   - keep the same regression-functional idea,
#   - but replace the Gaussian joint law by another log-concave joint distribution
#     whose conditional mean is still known or numerically tractable.
# ============================================================
