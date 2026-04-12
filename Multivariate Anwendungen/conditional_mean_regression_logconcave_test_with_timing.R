# ============================================================
# conditional_mean_regression_logconcave_test_with_timing.R
#
# Same regression benchmark as before, but now with explicit timing
# separation for
#   (i) training / fitting time
#   (ii) evaluation / prediction time
#
# The script reports:
#   - single-run regression metrics
#   - single-run timing table
#   - replicated raw results including times
#   - replicated timing summary including medians
#   - optional timing plots
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

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::median(x)
}


# ------------------------------------------------------------
# (2) Data-generating process: simple log-concave regression model
# ------------------------------------------------------------
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

conditional_mean_from_sm_score <- function(x0,
                                           y_grid,
                                           fit_sm) {
  eval_grid <- cbind(rep(x0, length(y_grid)), y_grid)
  score_mat <- predict_score_mv_basic(eval_grid, fit_sm)

  r_y <- as.numeric(score_mat[, 2])
  if (length(r_y) != length(y_grid)) return(NA_real_)
  if (any(!is.finite(r_y))) return(NA_real_)

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
#     with separate training-time measurements
# ------------------------------------------------------------
fit_all_models_timed <- function(train_xy) {
  fits <- list()
  timing_rows <- list()

  t_kde <- system.time({
    fits$KDE <- fit_kde_mv(train_xy, H_method = "Hpi", diagonal = FALSE)
  })
  timing_rows[[1L]] <- data.frame(
    method = "KDE",
    train_time_sec = as.numeric(t_kde["elapsed"]),
    stringsAsFactors = FALSE
  )

  t_mle <- system.time({
    fits$MLE <- fit_logconcave_mle_mv(train_xy, smoothed = FALSE)
  })
  timing_rows[[2L]] <- data.frame(
    method = "MLE",
    train_time_sec = as.numeric(t_mle["elapsed"]),
    stringsAsFactors = FALSE
  )

  t_sm <- system.time({
    fits$SM <- fit_score_matching_mv_basic(
      train_xy,
      m = 2,
      include_interactions = TRUE,
      standardize = TRUE,
      ridge = 1e-6,
      log_concave = TRUE,
      lc_method = "m2"
    )
  })
  timing_rows[[3L]] <- data.frame(
    method = "SM",
    train_time_sec = as.numeric(t_sm["elapsed"]),
    stringsAsFactors = FALSE
  )

  list(
    fits = fits,
    timing = do.call(rbind, timing_rows)
  )
}


# ------------------------------------------------------------
# (5) Evaluate conditional mean regression on a test set
#     with separate evaluation-time measurements
# ------------------------------------------------------------
evaluate_conditional_mean_regression_timed <- function(fits,
                                                       x_test,
                                                       y_test,
                                                       true_mean_test,
                                                       y_grid) {
  pred_kde <- pred_mle <- pred_sm <- NULL

  t_kde <- system.time({
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
  })

  t_mle <- system.time({
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
  })

  t_sm <- system.time({
    pred_sm <- vapply(
      x_test,
      FUN = function(x0) conditional_mean_from_sm_score(x0, y_grid, fits$SM),
      FUN.VALUE = numeric(1)
    )
  })

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

  timing <- data.frame(
    method = c("KDE", "MLE", "SM"),
    eval_time_sec = c(
      as.numeric(t_kde["elapsed"]),
      as.numeric(t_mle["elapsed"]),
      as.numeric(t_sm["elapsed"])
    ),
    stringsAsFactors = FALSE
  )

  list(
    predictions = out,
    metrics = as.data.frame(metrics),
    timing = timing
  )
}


# ------------------------------------------------------------
# (6) One complete experiment
# ------------------------------------------------------------
run_one_regression_experiment <- function(n_train = 400,
                                          n_test = 300,
                                          rho = 0.7,
                                          seed = 123) {
  set.seed(seed)

  train_xy <- r_joint_bvn(n_train, rho = rho)

  test_xy <- r_joint_bvn(n_test, rho = rho)
  x_test <- test_xy[, 1]
  y_test <- test_xy[, 2]
  true_mean_test <- true_cond_mean_bvn(x_test, rho = rho)

  y_all <- c(train_xy[, 2], test_xy[, 2])
  y_sd <- stats::sd(y_all)
  if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1
  y_grid <- seq(min(y_all) - 4 * y_sd,
                max(y_all) + 4 * y_sd,
                length.out = 601)

  fit_obj <- fit_all_models_timed(train_xy)
  eval_obj <- evaluate_conditional_mean_regression_timed(
    fits = fit_obj$fits,
    x_test = x_test,
    y_test = y_test,
    true_mean_test = true_mean_test,
    y_grid = y_grid
  )

  timing <- merge(
    fit_obj$timing,
    eval_obj$timing,
    by = "method",
    all = TRUE,
    sort = FALSE
  )
  timing$total_time_sec <- timing$train_time_sec + timing$eval_time_sec

  list(
    train_xy = train_xy,
    test_xy = test_xy,
    y_grid = y_grid,
    fits = fit_obj$fits,
    predictions = eval_obj$predictions,
    metrics = eval_obj$metrics,
    timing = timing[match(c("KDE", "MLE", "SM"), timing$method), ]
  )
}


# ------------------------------------------------------------
# (7) Replicated benchmark
# ------------------------------------------------------------
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

    metrics_tab <- one$metrics
    metrics_tab$method <- rownames(metrics_tab)
    rownames(metrics_tab) <- NULL

    timing_tab <- one$timing

    tab <- merge(metrics_tab, timing_tab, by = "method", all.x = TRUE, sort = FALSE)
    tab$rep <- r

    all_rows[[r]] <- tab[, c(
      "rep", "method", "mean_function_mse", "prediction_mse", "n_valid",
      "train_time_sec", "eval_time_sec", "total_time_sec"
    )]
  }

  res <- do.call(rbind, all_rows)

  summary_metrics <- do.call(rbind, lapply(split(res, res$method), function(dd) {
    data.frame(
      method = dd$method[1],
      mean_function_mse.mean = mean(dd$mean_function_mse, na.rm = TRUE),
      mean_function_mse.sd = stats::sd(dd$mean_function_mse, na.rm = TRUE),
      prediction_mse.mean = mean(dd$prediction_mse, na.rm = TRUE),
      prediction_mse.sd = stats::sd(dd$prediction_mse, na.rm = TRUE),
      train_time_sec.median = safe_median(dd$train_time_sec),
      train_time_sec.mean = mean(dd$train_time_sec, na.rm = TRUE),
      train_time_sec.sd = stats::sd(dd$train_time_sec, na.rm = TRUE),
      eval_time_sec.median = safe_median(dd$eval_time_sec),
      eval_time_sec.mean = mean(dd$eval_time_sec, na.rm = TRUE),
      eval_time_sec.sd = stats::sd(dd$eval_time_sec, na.rm = TRUE),
      total_time_sec.median = safe_median(dd$total_time_sec),
      total_time_sec.mean = mean(dd$total_time_sec, na.rm = TRUE),
      total_time_sec.sd = stats::sd(dd$total_time_sec, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary_metrics) <- NULL

  list(raw = res, summary = summary_metrics)
}


# ------------------------------------------------------------
# (8) Timing extraction helpers
# ------------------------------------------------------------
get_single_run_timing_table <- function(example_run) {
  example_run$timing[, c("method", "train_time_sec", "eval_time_sec", "total_time_sec")]
}

get_benchmark_time_medians <- function(benchmark_res) {
  benchmark_res$summary[, c(
    "method",
    "train_time_sec.median",
    "eval_time_sec.median",
    "total_time_sec.median"
  )]
}


# ------------------------------------------------------------
# (9) Minimal example run
# ------------------------------------------------------------
example_run <- run_one_regression_experiment(
  n_train = 400,
  n_test = 300,
  rho = 0.7,
  seed = 1
)

print("Single-run metrics:")
print(example_run$metrics)

print("Single-run timings (seconds):")
print(get_single_run_timing_table(example_run))

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

print("Replicated benchmark: timing medians (seconds)")
print(get_benchmark_time_medians(benchmark_res))


# ------------------------------------------------------------
# (10) Optional plots
# ------------------------------------------------------------
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

plot_single_run_timings <- function(example_run) {
  timing <- get_single_run_timing_table(example_run)
  mat <- t(as.matrix(timing[, c("train_time_sec", "eval_time_sec")]))
  colnames(mat) <- timing$method

  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  barplot(
    mat,
    beside = FALSE,
    legend.text = c("Training", "Evaluation"),
    args.legend = list(x = "topright", bty = "n"),
    ylab = "seconds",
    main = "Single-run timing by method"
  )
}

plot_benchmark_time_medians <- function(benchmark_res) {
  med <- get_benchmark_time_medians(benchmark_res)
  mat <- t(as.matrix(med[, c(
    "train_time_sec.median",
    "eval_time_sec.median"
  )]))
  colnames(mat) <- med$method

  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  barplot(
    mat,
    beside = FALSE,
    legend.text = c("Median training", "Median evaluation"),
    args.legend = list(x = "topright", bty = "n"),
    ylab = "seconds",
    main = "Median timing across replicated benchmark"
  )
}

plot_benchmark_time_boxplots <- function(benchmark_res) {
  raw <- benchmark_res$raw
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(1, 2))
  boxplot(train_time_sec ~ method,
          data = raw,
          ylab = "seconds",
          main = "Training times")
  boxplot(eval_time_sec ~ method,
          data = raw,
          ylab = "seconds",
          main = "Evaluation times")
}

# Uncomment if you want the plots immediately:
# plot_conditional_mean_curves(example_run)
# plot_single_run_timings(example_run)
# plot_benchmark_time_medians(benchmark_res)
# plot_benchmark_time_boxplots(benchmark_res)

