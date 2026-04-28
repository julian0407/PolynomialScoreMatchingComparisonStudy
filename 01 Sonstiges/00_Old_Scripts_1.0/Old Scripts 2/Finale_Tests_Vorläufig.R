# ============================================================
# final_tests.R
# Standardisierte finale Tests für Dichte-, Score-, Form- und
# Numerik-Metriken
# ============================================================

source("helper_functions.R")
source("KDE.R")
source("LogConcaveMLE.R")
source("Univariate_Polynomial_Score_Matching_1.0.R")
source("Multivariate_Pairwise_Polynomial_Score_Matching.R")
source("evaluation_metrics.R")

# ------------------------------------------------------------
# (1) Eine Methode spezifizieren
# ------------------------------------------------------------
# Beispiel:
# list(
#   label = "SM_m2",
#   method = "SM",
#   smoothed = FALSE,
#   fit_args = list(m = 2, include_interactions = TRUE, log_concave = TRUE, lc_method = "m2"),
#   predict_args = list()
# )

# ------------------------------------------------------------
# (2) Ein einzelner Lauf
# ------------------------------------------------------------

run_one_final_experiment <- function(n,
                                     family = c("univariate", "multivariate"),
                                     estimator_spec,
                                     r_sample,
                                     metrics,
                                     n_test = 2000,
                                     true_density = NULL,
                                     true_logdensity = NULL,
                                     true_score = NULL,
                                     grid_1d = NULL) {
  family <- match.arg(family)
  
  x_train <- r_sample(n)
  if (family == "univariate") x_train <- as.numeric(x_train) else x_train <- as.matrix(x_train)
  
  timing <- system.time({
    fit <- tryCatch(
      do.call(
        fit_estimator_generic,
        c(
          list(
            x = x_train,
            family = family,
            method = estimator_spec$method,
            smoothed = estimator_spec$smoothed
          ),
          estimator_spec$fit_args
        )
      ),
      error = function(e) structure(list(error_message = conditionMessage(e)), class = "fit_error")
    )
  })
  
  if (inherits(fit, "fit_error")) {
    row <- data.frame(
      n = n,
      method_label = estimator_spec$label,
      method = estimator_spec$method,
      runtime_sec = as.numeric(timing["elapsed"]),
      success = FALSE,
      status = fit$error_message,
      iterations = NA_real_,
      objective_value = NA_real_,
      stringsAsFactors = FALSE
    )
    
    for (nm in c("negloglik", "kl", "hellinger2", "ise", "score_loss", "fisher", "score_rmse",
                 "min_hessian_eigenvalue", "share_violated", "mean_violation", "max_violation")) {
      row[[nm]] <- NA_real_
    }
    return(row)
  }
  
  x_test <- r_sample(n_test)
  if (family == "univariate") x_test <- as.numeric(x_test) else x_test <- as.matrix(x_test)
  
  diags <- extract_fit_diagnostics(fit)
  
  metric_values <- tryCatch(
    evaluate_requested_metrics(
      metrics = metrics,
      x_test = x_test,
      fit = fit,
      family = family,
      method = estimator_spec$method,
      true_density = true_density,
      true_logdensity = true_logdensity,
      true_score = true_score,
      grid_1d = grid_1d,
      predict_args = estimator_spec$predict_args
    ),
    error = function(e) {
      warning(sprintf("Metric evaluation failed for %s at n=%s: %s",
                      estimator_spec$label, n, conditionMessage(e)))
      list()
    }
  )
  
  row <- data.frame(
    n = n,
    method_label = estimator_spec$label,
    method = estimator_spec$method,
    runtime_sec = as.numeric(timing["elapsed"]),
    success = diags$success,
    status = diags$status,
    iterations = diags$iterations,
    objective_value = diags$objective_value,
    stringsAsFactors = FALSE
  )
  
  for (nm in names(metric_values)) {
    row[[nm]] <- metric_values[[nm]]
  }
  
  row
}

# ------------------------------------------------------------
# (3) Gesamtbenchmark
# ------------------------------------------------------------

run_final_benchmark <- function(sample_sizes,
                                family = c("univariate", "multivariate"),
                                estimator_specs,
                                r_sample,
                                metrics,
                                n_rep = 20,
                                n_test = 2000,
                                true_density = NULL,
                                true_logdensity = NULL,
                                true_score = NULL,
                                grid_1d = NULL,
                                seed = NULL,
                                verbose = TRUE) {
  family <- match.arg(family)
  
  if (!is.null(seed)) set.seed(seed)
  
  out <- list()
  counter <- 1L
  
  for (spec in estimator_specs) {
    if (verbose) message("Method: ", spec$label)
    
    for (n in sample_sizes) {
      if (verbose) message("  n = ", n)
      
      for (rep in seq_len(n_rep)) {
        if (verbose) message("    repetition ", rep, "/", n_rep)
        
        ans <- run_one_final_experiment(
          n = n,
          family = family,
          estimator_spec = spec,
          r_sample = r_sample,
          metrics = metrics,
          n_test = n_test,
          true_density = true_density,
          true_logdensity = true_logdensity,
          true_score = true_score,
          grid_1d = grid_1d
        )
        ans$repetition <- rep
        out[[counter]] <- ans
        counter <- counter + 1L
      }
    }
  }
  
  raw <- do.call(rbind, out)
  rownames(raw) <- NULL
  
  structure(
    list(
      raw = raw,
      settings = list(
        sample_sizes = sample_sizes,
        family = family,
        metrics = metrics,
        n_rep = n_rep,
        n_test = n_test
      )
    ),
    class = "final_benchmark"
  )
}

# ------------------------------------------------------------
# (4) Aggregation
# ------------------------------------------------------------

aggregate_final_benchmark <- function(obj, metric) {
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")
  df <- obj$raw
  
  if (!metric %in% names(df)) {
    stop(sprintf("Metric '%s' not found in benchmark output.", metric))
  }
  
  split_key <- interaction(df$method_label, df$n, drop = TRUE)
  
  agg_list <- lapply(split(df, split_key), function(dd) {
    data.frame(
      method_label = dd$method_label[1],
      n = dd$n[1],
      mean = mean(dd[[metric]], na.rm = TRUE),
      median = median(dd[[metric]], na.rm = TRUE),
      sd = stats::sd(dd[[metric]], na.rm = TRUE),
      q25 = stats::quantile(dd[[metric]], probs = 0.25, na.rm = TRUE, names = FALSE),
      q75 = stats::quantile(dd[[metric]], probs = 0.75, na.rm = TRUE, names = FALSE),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, agg_list)
}

# ------------------------------------------------------------
# (5) Plot
# ------------------------------------------------------------

plot_metric_comparison <- function(obj,
                                   metric,
                                   log_x = TRUE,
                                   log_y = FALSE,
                                   use_median = FALSE,
                                   main = NULL,
                                   xlab = "sample size n",
                                   ylab = NULL) {
  agg <- aggregate_final_benchmark(obj, metric)
  methods <- unique(agg$method_label)
  
  if (is.null(main)) main <- paste("Metric:", metric)
  if (is.null(ylab)) ylab <- metric
  
  ycol <- if (use_median) "median" else "mean"
  
  xlim <- range(agg$n, na.rm = TRUE)
  ylim <- range(agg[[ycol]], na.rm = TRUE)
  
  plot(
    NA, NA,
    xlim = xlim,
    ylim = ylim,
    log = paste0(if (log_x) "x" else "", if (log_y) "y" else ""),
    xlab = xlab,
    ylab = ylab,
    main = main
  )
  
  for (i in seq_along(methods)) {
    dd <- agg[agg$method_label == methods[i], ]
    dd <- dd[order(dd$n), ]
    lines(dd$n, dd[[ycol]], type = "b", lwd = 2, pch = 19 + i - 1, col = i)
  }
  
  legend("topright", legend = methods, col = seq_along(methods),
         lty = 1, pch = 19 + seq_along(methods) - 1, bty = "n")
}

plot_metric_panel <- function(obj,
                              metrics,
                              log_x = TRUE,
                              log_y = FALSE,
                              use_median = FALSE) {
  k <- length(metrics)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(1, k))
  for (met in metrics) {
    plot_metric_comparison(
      obj = obj,
      metric = met,
      log_x = log_x,
      log_y = log_y,
      use_median = use_median
    )
  }
}

# ------------------------------------------------------------
# (6) Beispiel 1: univariat
# ------------------------------------------------------------

r_sample_norm_1d <- function(n) rnorm(n, mean = 0, sd = 1)
true_density_norm_1d <- function(x) dnorm(x, mean = 0, sd = 1)
true_logdensity_norm_1d <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)
true_score_norm_1d <- function(x) matrix(x, ncol = 1)

grid_1d <- seq(-5, 5, length.out = 2001)

estimators_1d <- list(
  list(
    label = "KDE_SJ",
    method = "KDE",
    smoothed = FALSE,
    fit_args = list(bw = "SJ"),
    predict_args = list()
  ),
  # list(
  #   label = "MLE_unsmoothed",
  #   method = "MLE",
  #   smoothed = FALSE,
  #   fit_args = list(),
  #   predict_args = list()
  # ),
  list(
    label = "SM_m2",
    method = "SM",
    smoothed = FALSE,
    fit_args = list(
      m = 2,
      standardize = TRUE
    ),
    predict_args = list(
      subdivisions = 200L,
      rel.tol = 1e-8,
      stop_on_failure = FALSE
    )
  )
)

res_1d <- run_final_benchmark(
  sample_sizes = c(50, 100, 200, 500, 1000),
  family = "univariate",
  estimator_specs = estimators_1d,
  r_sample = r_sample_norm_1d,
  metrics = c("negloglik", "ise", "score_loss", "score_rmse", "logconcavity"),
  n_rep = 20,
  n_test = 2000,
  true_density = true_density_norm_1d,
  true_logdensity = true_logdensity_norm_1d,
  true_score = true_score_norm_1d,
  grid_1d = grid_1d,
  seed = 123
)

plot_metric_panel(res_1d, metrics = c("negloglik", "score_loss"))
plot_metric_panel(res_1d, metrics = c("ise", "score_rmse"))

# Numerik zusammenfassen:
# by(res_1d$raw, res_1d$raw$method_label, summarize_numerics)

# ------------------------------------------------------------
# (7) Beispiel 2: multivariat
# ------------------------------------------------------------

r_sample_norm_2d <- function(n) {
  cbind(rnorm(n), rnorm(n))
}

true_score_norm_2d <- function(x) {
  x <- as.matrix(x)
  x
}

estimators_2d <- list(
  list(
    label = "KDE_Hpi",
    method = "KDE",
    smoothed = FALSE,
    fit_args = list(H_method = "Hpi"),
    predict_args = list()
  ),
  list(
    label = "MLE_mv",
    method = "MLE",
    smoothed = FALSE,
    fit_args = list(),
    predict_args = list()
  ),
  list(
    label = "SM_basic",
    method = "SM",
    smoothed = FALSE,
    fit_args = list(
      m = 2,
      include_interactions = TRUE,
      standardize = TRUE,
      ridge = 1e-6,
      log_concave = FALSE
    ),
    predict_args = list()
  ),
  list(
    label = "SM_grid",
    method = "SM",
    smoothed = FALSE,
    fit_args = list(
      m = 2,
      include_interactions = TRUE,
      standardize = TRUE,
      ridge = 1e-6,
      log_concave = TRUE,
      lc_method = "grid",
      lc_grid_size = 5L,
      lc_penalty = 1e4
    ),
    predict_args = list()
  )
)

# Achtung:
# Für multivariates SM werden hier bewusst keine Dichte-Metriken gerechnet.
res_2d <- run_final_benchmark(
  sample_sizes = c(100, 200, 500, 1000),
  family = "multivariate",
  estimator_specs = estimators_2d,
  r_sample = r_sample_norm_2d,
  metrics = c("score_loss", "score_rmse", "logconcavity"),
  n_rep = 20,
  n_test = 2000,
  true_score = true_score_norm_2d,
  seed = 123
)

plot_metric_panel(res_2d, metrics = c("score_loss", "score_rmse"))

# Log-Konkavitätsmetriken ansehen:
# aggregate_final_benchmark(res_2d, "min_hessian_eigenvalue")
# aggregate_final_benchmark(res_2d, "share_violated")