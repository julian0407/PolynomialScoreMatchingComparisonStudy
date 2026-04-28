# ============================================================
# final_tests_robust.R
# Robuste finale Tests für Dichte-, Score-, Form- und
# Numerik-Metriken
#
# Wichtige Änderungen gegenüber der Vorversion:
# 1) Robustes Zusammenführen der Läufe auch bei unterschiedlichen
#    Rückgabe-Spalten (kein rbind-Fehler mehr).
# 2) Einheitliches Spaltenschema pro Benchmark, damit alle Methoden
#    dieselben Metrikspalten besitzen.
# 3) Schnelle Metrik-Auswertung bleibt erhalten:
#    zuerst werden alle angefragten Metriken gemeinsam berechnet;
#    nur falls das fehlschlägt, wird auf eine robuste Fallback-Logik
#    pro Metrik gewechselt.
# 4) Plot/Aggregation brechen nicht mehr hart ab, wenn eine Metrik für
#    eine Methode nicht verfügbar ist; stattdessen werden NAs geführt
#    und Plots für vollständig fehlende Metriken sauber übersprungen.
# ============================================================

source("helper_functions.R")
source("KDE.R")
source("LogConcaveMLE.R")
source("Univariate_Polynomial_Score_Matching_1.0.R")
source("Multivariate_Pairwise_Polynomial_Score_Matching.R")
source("Draft_Evaluation_Metrics.R")
# source("BiasVariance_Score.R")

# ------------------------------------------------------------
# Hilfsfunktionen
# ------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::median(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x)
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

# Welche Ausgabespalten können aus einer angeforderten Metrik entstehen?
# "logconcavity" ist eine Gruppen-Metrik und erzeugt mehrere Spalten.
metric_to_output_columns <- function(metrics) {
  metrics <- unique(metrics)
  out <- character(0)
  
  for (m in metrics) {
    if (identical(m, "logconcavity")) {
      out <- c(
        out,
        "min_hessian_eigenvalue",
        "share_violated",
        "mean_violation",
        "max_violation"
      )
    } else {
      out <- c(out, m)
      
      if (m %in% c("negloglik", "kl")) {
        out <- c(
          out,
          paste0(m, "_na_share"),
          paste0(m, "_tail_share"),
          paste0(m, "_outlier_dominated")
        )
      }
    }
  }
  
  unique(out)
}

make_empty_metric_list <- function(metric_columns) {
  out <- as.list(rep(NA_real_, length(metric_columns)))
  names(out) <- metric_columns
  out
}

coerce_row_to_schema <- function(row, schema_names) {
  missing_cols <- setdiff(schema_names, names(row))
  if (length(missing_cols) > 0L) {
    for (nm in missing_cols) row[[nm]] <- NA
  }
  row <- row[, schema_names, drop = FALSE]
  row
}

bind_rows_fill_base <- function(rows) {
  if (length(rows) == 0L) return(data.frame())
  
  schema <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows2 <- lapply(rows, coerce_row_to_schema, schema_names = schema)
  out <- do.call(rbind, rows2)
  rownames(out) <- NULL
  out
}

# Robuste, aber möglichst schnelle Metrikauswertung:
# - zuerst alles auf einmal
# - nur bei Fehler Fallback pro Metrik / Metrik-Gruppe
safe_evaluate_requested_metrics <- function(metrics,
                                            x_test,
                                            fit,
                                            family,
                                            method,
                                            true_density = NULL,
                                            true_logdensity = NULL,
                                            true_score = NULL,
                                            grid_1d = NULL,
                                            density_predict_args = list(),
                                            score_predict_args = list(),
                                            density_metric_args = list(),
                                            score_metric_args = list(),
                                            estimator_label = method,
                                            n = NA_integer_,
                                            verbose = TRUE) {
  metric_columns <- metric_to_output_columns(metrics)
  empty_out <- make_empty_metric_list(metric_columns)
  
  # schneller Pfad: alles in einem Aufruf
  full_res <- tryCatch(
    evaluate_requested_metrics(
      metrics = metrics,
      x_test = x_test,
      fit = fit,
      family = family,
      method = method,
      true_density = true_density,
      true_logdensity = true_logdensity,
      true_score = true_score,
      grid_1d = grid_1d,
      density_predict_args = density_predict_args,
      score_predict_args = score_predict_args,
      density_metric_args = density_metric_args,
      score_metric_args = score_metric_args
    ),
    error = function(e) e
  )
  
  if (!inherits(full_res, "error")) {
    for (nm in names(full_res)) empty_out[[nm]] <- full_res[[nm]]
    return(empty_out)
  }
  
  if (isTRUE(verbose)) {
    warning(sprintf(
      "Joint metric evaluation failed for %s at n=%s. Falling back to robust metric-wise evaluation. Reason: %s",
      estimator_label, n, conditionMessage(full_res)
    ))
  }
  
  # Fallback: pro angefragter Metrik / Metrik-Gruppe
  for (m in metrics) {
    one_res <- tryCatch(
      evaluate_requested_metrics(
        metrics = m,
        x_test = x_test,
        fit = fit,
        family = family,
        method = method,
        true_density = true_density,
        true_logdensity = true_logdensity,
        true_score = true_score,
        grid_1d = grid_1d,
        density_predict_args = density_predict_args,
        score_predict_args = score_predict_args,
        density_metric_args = density_metric_args,
        score_metric_args = score_metric_args
      ),
      error = function(e) e
    )
    
    if (inherits(one_res, "error")) {
      if (isTRUE(verbose)) {
        warning(sprintf(
          "Metric '%s' failed for %s at n=%s: %s",
          m, estimator_label, n, conditionMessage(one_res)
        ))
      }
      next
    }
    
    for (nm in names(one_res)) empty_out[[nm]] <- one_res[[nm]]
  }
  
  empty_out
}

# ------------------------------------------------------------
# (1) Eine Methode spezifizieren
# ------------------------------------------------------------
# Beispiel:
# list(
#   label = "SM_m2",
#   method = "SM",
#   smoothed = FALSE,
#   fit_args = list(m = 2, include_interactions = TRUE, log_concave = TRUE, lc_method = "m2"),
#   density_predict_args = list(),
#   score_predict_args = list(),
#   score_metric_args = list()
# )

# ------------------------------------------------------------
# (2) Ein einzelner Lauf
# ------------------------------------------------------------

run_one_final_experiment <- function(n,
                                     family = c("univariate", "multivariate"),
                                     estimator_spec,
                                     r_sample,
                                     metrics,
                                     metric_columns = NULL,
                                     n_test = 2000,
                                     true_density = NULL,
                                     true_logdensity = NULL,
                                     true_score = NULL,
                                     grid_1d = NULL,
                                     verbose = TRUE) {
  family <- match.arg(family)
  metric_columns <- metric_columns %||% metric_to_output_columns(metrics)
  
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
  
  base_row <- data.frame(
    n = n,
    method_label = estimator_spec$label,
    method = estimator_spec$method,
    runtime_sec = as.numeric(timing["elapsed"]),
    success = FALSE,
    status = NA_character_,
    iterations = NA_real_,
    objective_value = NA_real_,
    kappa_raw = NA_real_,
    kappa_reg = NA_real_,
    rcond_raw = NA_real_,
    rcond_reg = NA_real_,
    eigmin_raw = NA_real_,
    eigmin_reg = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (nm in metric_columns) base_row[[nm]] <- NA_real_
  
  if (inherits(fit, "fit_error")) {
    base_row$status <- fit$error_message
    return(base_row)
  }
  
  x_test <- r_sample(n_test)
  if (family == "univariate") x_test <- as.numeric(x_test) else x_test <- as.matrix(x_test)
  
  diags <- tryCatch(
    extract_fit_diagnostics(fit),
    error = function(e) {
      list(
        success = TRUE,
        status = paste("diagnostics_error:", conditionMessage(e)),
        iterations = NA_real_,
        objective_value = NA_real_
      )
    }
  )
  
  metric_values <- safe_evaluate_requested_metrics(
    metrics = metrics,
    x_test = x_test,
    fit = fit,
    family = family,
    method = estimator_spec$method,
    true_density = true_density,
    true_logdensity = true_logdensity,
    true_score = true_score,
    grid_1d = grid_1d,
    density_predict_args = estimator_spec$density_predict_args %||% list(),
    score_predict_args = estimator_spec$score_predict_args %||% list(),
    density_metric_args = estimator_spec$density_metric_args %||% list(),
    score_metric_args = estimator_spec$score_metric_args %||% list(),
    estimator_label = estimator_spec$label,
    n = n,
    verbose = verbose
  )
  
  base_row$success <- diags$success %||% TRUE
  base_row$status <- diags$status %||% NA_character_
  base_row$iterations <- diags$iterations %||% NA_real_
  base_row$objective_value <- diags$objective_value %||% NA_real_
  base_row$kappa_raw <- diags$kappa_raw %||% NA_real_
  base_row$kappa_reg <- diags$kappa_reg %||% NA_real_
  base_row$rcond_raw <- diags$rcond_raw %||% NA_real_
  base_row$rcond_reg <- diags$rcond_reg %||% NA_real_
  base_row$eigmin_raw <- diags$eigmin_raw %||% NA_real_
  base_row$eigmin_reg <- diags$eigmin_reg %||% NA_real_
  
  for (nm in names(metric_values)) {
    base_row[[nm]] <- metric_values[[nm]]
  }
  
  base_row
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
  metric_columns <- metric_to_output_columns(metrics)
  
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
          metric_columns = metric_columns,
          n_test = n_test,
          true_density = true_density,
          true_logdensity = true_logdensity,
          true_score = true_score,
          grid_1d = grid_1d,
          verbose = verbose
        )
        ans$repetition <- rep
        out[[counter]] <- ans
        counter <- counter + 1L
      }
    }
  }
  
  raw <- bind_rows_fill_base(out)
  
  preferred_order <- c(
    "n", "repetition", "method_label", "method", "runtime_sec",
    "success", "status", "iterations", "objective_value",
    metric_columns
  )
  other_cols <- setdiff(names(raw), preferred_order)
  raw <- raw[, c(intersect(preferred_order, names(raw)), other_cols), drop = FALSE]
  rownames(raw) <- NULL
  
  structure(
    list(
      raw = raw,
      settings = list(
        sample_sizes = sample_sizes,
        family = family,
        metrics = metrics,
        metric_columns = metric_columns,
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

aggregate_final_benchmark <- function(obj, metric, drop_all_na = FALSE) {
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")
  df <- obj$raw
  
  if (!metric %in% names(df)) {
    available <- intersect(
      c("negloglik", "kl", "hellinger2", "ise", "score_loss", "fisher", "score_rmse",
        "min_hessian_eigenvalue", "share_violated", "mean_violation", "max_violation"),
      names(df)
    )
    stop(sprintf(
      "Metric '%s' not found in benchmark output. Available metrics: %s",
      metric,
      paste(available, collapse = ", ")
    ))
  }
  
  split_key <- interaction(df$method_label, df$n, drop = TRUE)
  
  na_share_col <- paste0(metric, "_na_share")
  tail_share_col <- paste0(metric, "_tail_share")
  outlier_dom_col <- paste0(metric, "_outlier_dominated")
  
  has_na_diag <- na_share_col %in% names(df)
  has_tail_diag <- tail_share_col %in% names(df)
  has_outlier_diag <- outlier_dom_col %in% names(df)
  
  agg_list <- lapply(split(df, split_key), function(dd) {
    x <- dd[[metric]]
    
    out <- data.frame(
      method_label = dd$method_label[1],
      n = dd$n[1],
      n_non_na = sum(is.finite(x)),
      mean = safe_mean(x),
      median = safe_median(x),
      sd = safe_sd(x),
      q25 = safe_quantile(x, 0.25),
      q75 = safe_quantile(x, 0.75),
      stringsAsFactors = FALSE
    )
    
    if (has_na_diag) {
      na_sh <- dd[[na_share_col]]
      out$run_with_any_na_share <- safe_mean(na_sh > 0)
      out$mean_na_share_within_run <- safe_mean(na_sh)
    }
    
    if (has_tail_diag) {
      out$mean_tail_share_within_run <- safe_mean(dd[[tail_share_col]])
    }
    
    if (has_outlier_diag) {
      out$run_outlier_dominated_share <- safe_mean(dd[[outlier_dom_col]])
    }
    
    out
  })
  
  agg <- do.call(rbind, agg_list)
  rownames(agg) <- NULL
  
  if (isTRUE(drop_all_na)) {
    agg <- agg[agg$n_non_na > 0L, , drop = FALSE]
  }
  
  agg
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
                                   ylab = NULL,
                                   warn_if_missing = TRUE) {
  agg <- aggregate_final_benchmark(obj, metric)
  ycol <- if (use_median) "median" else "mean"
  
  keep <- is.finite(agg[[ycol]]) & is.finite(agg$n)
  agg2 <- agg[keep, , drop = FALSE]
  
  if (nrow(agg2) == 0L) {
    if (isTRUE(warn_if_missing)) {
      warning(sprintf("Metric '%s' has no finite values to plot. Plot skipped.", metric))
    }
    plot.new()
    title(main = paste("Metric:", metric))
    text(0.5, 0.5, labels = sprintf("No finite values for '%s'", metric))
    return(invisible(NULL))
  }
  
  methods <- unique(agg2$method_label)
  
  if (is.null(main)) main <- paste("Metric:", metric)
  if (is.null(ylab)) ylab <- metric
  
  xlim <- range(agg2$n, na.rm = TRUE)
  ylim <- range(agg2[[ycol]], na.rm = TRUE)
  
  # Falls log_y angefragt ist, dürfen nur positive Werte geplottet werden.
  if (isTRUE(log_y)) {
    pos <- agg2[[ycol]] > 0 & is.finite(agg2[[ycol]])
    agg2 <- agg2[pos, , drop = FALSE]
    if (nrow(agg2) == 0L) {
      warning(sprintf("Metric '%s' has no positive finite values for log_y plot. Plot skipped.", metric))
      plot.new()
      title(main = paste("Metric:", metric))
      text(0.5, 0.5, labels = sprintf("No positive values for log_y: '%s'", metric))
      return(invisible(NULL))
    }
    xlim <- range(agg2$n, na.rm = TRUE)
    ylim <- range(agg2[[ycol]], na.rm = TRUE)
    methods <- unique(agg2$method_label)
  }
  
  plot(
    NA, NA,
    xlim = xlim,
    ylim = ylim,
    log = paste0(if (log_x) "x" else "", if (log_y) "y" else ""),
    xlab = xlab,
    ylab = ylab,
    main = main
  )
  
  drawn_methods <- character(0)
  
  for (i in seq_along(methods)) {
    dd <- agg2[agg2$method_label == methods[i], , drop = FALSE]
    dd <- dd[order(dd$n), , drop = FALSE]
    if (nrow(dd) == 0L) next
    
    lines(dd$n, dd[[ycol]], type = "b", lwd = 2, pch = 19 + i - 1, col = i)
    drawn_methods <- c(drawn_methods, methods[i])
  }
  
  if (length(drawn_methods) > 0L) {
    idx <- match(drawn_methods, methods)
    legend(
      "topright",
      legend = drawn_methods,
      col = idx,
      lty = 1,
      pch = 19 + idx - 1,
      bty = "n"
    )
  }
  
  invisible(agg2)
}

plot_metric_panel <- function(obj,
                              metrics,
                              log_x = TRUE,
                              log_y = FALSE,
                              use_median = FALSE) {
  k <- length(metrics)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  if (k>1) {
    par(mfrow = c(2, k/2))
  } else {
    par(mfrow = c(1, k))
  }
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

# r_sample_norm_1d <- function(n) rnorm(n, mean = 0, sd = 1)
# true_density_norm_1d <- function(x) dnorm(x, mean = 0, sd = 1)
# true_logdensity_norm_1d <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)
# true_score_norm_1d <- function(x) matrix(x, ncol = 1)
# 
# grid_1d <- seq(-5, 5, length.out = 2001)
# 
# estimators_1d <- list(
#   list(
#     label = "KDE_SJ",
#     method = "KDE",
#     smoothed = FALSE,
#     fit_args = list(bw = "SJ"),
#     predict_args = list()
#   ),
#   list(
#     label = "MLE_unsmoothed",
#     method = "MLE",
#     smoothed = FALSE,
#     fit_args = list(),
#     predict_args = list()
#   ),
#   list(
#     label = "SM_m2",
#     method = "SM",
#     smoothed = FALSE,
#     fit_args = list(
#       m = 2,
#       standardize = TRUE
#     ),
#     predict_args = list(
#       subdivisions = 200L,
#       rel.tol = 1e-8,
#       stop_on_failure = FALSE
#     )
#   )
# )
# 
# res_1d <- run_final_benchmark(
#   sample_sizes = c(50, 100),
#   family = "univariate",
#   estimator_specs = estimators_1d,
#   r_sample = r_sample_norm_1d,
#   metrics = c("negloglik", "ise", "score_loss", "score_rmse", "logconcavity"),
#   n_rep = 20,
#   n_test = 2000,
#   true_density = true_density_norm_1d,
#   true_logdensity = true_logdensity_norm_1d,
#   true_score = true_score_norm_1d,
#   grid_1d = grid_1d,
#   seed = 123
# )
# 
# plot_metric_panel(res_1d, metrics = c("negloglik", "score_loss"))
# plot_metric_panel(res_1d, metrics = c("ise", "score_rmse"))
# 
# # Numerik zusammenfassen:
# # by(res_1d$raw, res_1d$raw$method_label, summarize_numerics)
# 
# # ------------------------------------------------------------
# # (7) Beispiel 2: multivariat
# # ------------------------------------------------------------
# 
# r_sample_norm_2d <- function(n) {
#   cbind(rnorm(n), rnorm(n))
# }
# 
# true_score_norm_2d <- function(x) {
#   x <- as.matrix(x)
#   x
# }
# 
# estimators_2d <- list(
#   list(
#     label = "KDE_Hpi",
#     method = "KDE",
#     smoothed = FALSE,
#     fit_args = list(H_method = "Hpi"),
#     predict_args = list()
#   ),
#   list(
#     label = "MLE_mv",
#     method = "MLE",
#     smoothed = FALSE,
#     fit_args = list(),
#     predict_args = list()
#   ),
#   list(
#     label = "SM_basic",
#     method = "SM",
#     smoothed = FALSE,
#     fit_args = list(
#       m = 2,
#       include_interactions = TRUE,
#       standardize = TRUE,
#       ridge = 1e-6,
#       log_concave = FALSE
#     ),
#     predict_args = list()
#   ),
#   list(
#     label = "SM_grid",
#     method = "SM",
#     smoothed = FALSE,
#     fit_args = list(
#       m = 2,
#       include_interactions = TRUE,
#       standardize = TRUE,
#       ridge = 1e-6,
#       log_concave = TRUE,
#       lc_method = "grid",
#       lc_grid_size = 5L,
#       lc_penalty = 1e4
#     ),
#     predict_args = list()
#   )
# )

# Achtung:
# Für multivariates SM werden hier bewusst keine Dichte-Metriken gerechnet.
# res_2d <- run_final_benchmark(
#   sample_sizes = c(100, 200, 500, 1000),
#   family = "multivariate",
#   estimator_specs = estimators_2d,
#   r_sample = r_sample_norm_2d,
#   metrics = c("score_loss", "score_rmse", "logconcavity"),
#   n_rep = 20,
#   n_test = 2000,
#   true_score = true_score_norm_2d,
#   seed = 123
# )
# 
# plot_metric_panel(res_2d, metrics = c("score_loss", "score_rmse"))

# Log-Konkavitätsmetriken ansehen:
# aggregate_final_benchmark(res_2d, "min_hessian_eigenvalue")
# aggregate_final_benchmark(res_2d, "share_violated")






# ============================================================
# Iterative 1D smoke tests:
#   (A) Gaussian
#   (B) Laplace
# Kleine n, wenige Wiederholungen
# ============================================================
# 
# # ----------------------------
# # Estimatoren wie bisher
# # ----------------------------
# estimators_1d <- list(
#   # list(
#   #   label = "KDE_SJ",
#   #   method = "KDE",
#   #   smoothed = FALSE,
#   #   fit_args = list(bw = "SJ"),
#   #   # score_metric_args = list(robust = "none")
#   #   score_metric_args = list(robust = "trim", trim_alpha = 0.05)
#   # ),
#   # list(
#   #   label = "MLE_smoothed",
#   #   method = "MLE",
#   #   smoothed = TRUE,
#   #   fit_args = list(),
#   #   # score_metric_args = list(robust = "none")
#   #   score_metric_args = list(robust = "trim", trim_alpha = 0.05)
#   # ),
#   list(
#     label = "SM",
#     method = "SM",
#     smoothed = FALSE,
#     fit_args = list(
#       m = 6,
#       standardize = TRUE,
#       ridge = 1e-2
#       # h = function(z) h_tanh_sq(z, tau = 1),
#       # h_prime = function(z) h_tanh_sq_prime(z, tau = 1)
#     ),
#     density_predict_args = list(
#       subdivisions = 200L,
#       rel.tol = 1e-8,
#       stop_on_failure = TRUE
#     ),
#     density_metric_args = list(
#       robust = "median",
#       trim_alpha = 0.05,
#       outlier_dom_threshold = 0.25
#     ),
#     score_predict_args = list(),
#     # score_metric_args = list(robust = "none")
#     score_metric_args = list(robust = "trim", trim_alpha = 0.05)
#   )
# )
# 
# # Für die ersten Checks eher etwas breiterer Grid
# grid_1d <- seq(-6, 6, length.out = 1501)
# 
# # ============================================================
# # (A) 1D Gaussian
# # Score-Konvention in deinem Code: r(x) = - d/dx log f(x)
# # Für N(0,1): r(x) = x
# # ============================================================
# 
# r_sample_norm_1d <- function(n) {
#   rnorm(n, mean = -50, sd = 4)
# }
# 
# true_density_norm_1d <- function(x) {
#   dnorm(x, mean = -50, sd = 4)
# }
# 
# true_logdensity_norm_1d <- function(x) {
#   dnorm(x, mean = -50, sd = 4, log = TRUE)
# }
# 
# true_score_norm_1d <- function(x) {
#   x <- as.numeric(x)
#   matrix((x + 50) / 16, ncol = 1)
# }
# 
# # metrics_all_1d <- c("negloglik", "ise", "kl", "hellinger2", "score_loss", "score_rmse", "logconcavity")
# metrics_all_1d <- c("negloglik", "ise", "kl", "hellinger2", "score_loss", "score_rmse")
# 
# res_1d_norm_smoke <- run_final_benchmark(
#   sample_sizes = c(50, 200, 500, 1000),
#   family = "univariate",
#   estimator_specs = estimators_1d,
#   r_sample = r_sample_norm_1d,
#   metrics = metrics_all_1d,
#   n_rep = 20,
#   n_test = 1000,
#   true_density = true_density_norm_1d,
#   true_logdensity = true_logdensity_norm_1d,
#   true_score = true_score_norm_1d,
#   grid_1d = grid_1d,
#   seed = 123,
#   verbose = TRUE
# )
# 
# plot_metric_panel(res_1d_norm_smoke, metrics = c("negloglik"))
# plot_metric_panel(res_1d_norm_smoke, metrics = c("ise"))
# plot_metric_panel(res_1d_norm_smoke, metrics = c("kl"))
# plot_metric_panel(res_1d_norm_smoke, metrics = c("hellinger2"))
# plot_metric_panel(res_1d_norm_smoke, metrics = c("score_loss"))
# plot_metric_panel(res_1d_norm_smoke, metrics = c("score_rmse"))
# # plot_metric_panel(res_1d_norm_smoke, metrics = c("ise", "score_rmse"))
# plot_metric_panel(res_1d_norm_smoke, metrics = c("negloglik", "ise", "kl", "hellinger2", "score_loss", "score_rmse"))
# 
# # optional
# # aggregate_final_benchmark(res_1d_norm_smoke, "negloglik")
# # aggregate_final_benchmark(res_1d_norm_smoke, "score_loss")
# 
# # TODO:
# bv_1d_norm <- run_bias_variance_score_benchmark(
#   sample_sizes = c(20, 50, 100, 200, 500, 1000, 5000, 10000),
#   family = "univariate",
#   estimator_specs = estimators_1d,
#   r_sample = r_sample_norm_1d,
#   true_score = true_score_norm_1d,
#   n_rep = 100,
#   seed = 123,
#   verbose = TRUE,
#   n_probe_grid = 5000,
#   grid_size_1d = 401
# )
# 
# bv_1d_norm$summary
# plot_metrics = c("integrated_bias2", "integrated_variance", "integrated_mse")
# plot_metrics = c("integrated_bias2_trimmed_mean", "integrated_variance", "integrated_mse")
# plot_metrics = c("integrated_bias2_median", "integrated_variance", "integrated_mse")
# plot_bias_variance_panel(
#   bv_1d_norm,
#   metrics = plot_metrics
# )
# 
# 
# 
# # ============================================================
# # (B) 1D Laplace(0, b=1)
# # Dichte: f(x)=0.5 exp(-|x|)
# # Logdichte: log(0.5)-|x|
# # Score: r(x) = - d/dx log f(x) = sign(x), bei x=0 nicht differenzierbar
# # Für numerische Robustheit setzen wir dort 0
# # ============================================================
# 
# r_sample_laplace_1d <- function(n) {
#   u <- runif(n, min = -0.5, max = 0.5)
#   -sign(u) * log(1 - 2 * abs(u))
# }
# 
# true_density_laplace_1d <- function(x) {
#   0.5 * exp(-abs(x))
# }
# 
# true_logdensity_laplace_1d <- function(x) {
#   log(0.5) - abs(x)
# }
# 
# true_score_laplace_1d <- function(x) {
#   x <- as.numeric(x)
#   s <- sign(x)
#   s[x == 0] <- 0
#   matrix(s, ncol = 1)
# }
# 
# res_1d_laplace_smoke <- run_final_benchmark(
#   sample_sizes = c(50, 200, 500, 1000, 5000),
#   family = "univariate",
#   estimator_specs = estimators_1d,
#   r_sample = r_sample_laplace_1d,
#   metrics = metrics_all_1d,
#   n_rep = 20,
#   n_test = 2000,
#   true_density = true_density_laplace_1d,
#   true_logdensity = true_logdensity_laplace_1d,
#   true_score = true_score_laplace_1d,
#   grid_1d = grid_1d,
#   seed = 124,
#   verbose = TRUE
# )
# 
# plot_metric_panel(res_1d_laplace_smoke, metrics = c("score_loss"))
# # plot_metric_panel(res_1d_laplace_smoke, metrics = c("ise", "score_rmse"))
# plot_metric_panel(res_1d_laplace_smoke, metrics = c("negloglik", "ise", "kl", "hellinger2", "score_loss", "score_rmse"))
# 
# # optional
# # aggregate_final_benchmark(res_1d_laplace_smoke, "negloglik")
# # aggregate_final_benchmark(res_1d_laplace_smoke, "score_loss")
# 
# # ============================================================
# # (C) 1D Gumbel(mu, beta)  [Maximum-Typ]
# # CDF: F(x) = exp(-exp(-(x-mu)/beta))
# # Dichte: f(x) = (1/beta) * exp(-z - exp(-z)), z = (x-mu)/beta
# # Logdichte: -log(beta) - z - exp(-z)
# # Score: r(x) = - d/dx log f(x) = (1 - exp(-z)) / beta
# # ============================================================
# 
# mu_gumbel <- 0
# beta_gumbel <- 1
# 
# r_sample_gumbel_1d <- function(n, mu = mu_gumbel, beta = beta_gumbel) {
#   u <- runif(n)
#   mu - beta * log(-log(u))
# }
# 
# true_density_gumbel_1d <- function(x, mu = mu_gumbel, beta = beta_gumbel) {
#   x <- as.numeric(x)
#   z <- (x - mu) / beta
#   (1 / beta) * exp(-z - exp(-z))
# }
# 
# true_logdensity_gumbel_1d <- function(x, mu = mu_gumbel, beta = beta_gumbel) {
#   x <- as.numeric(x)
#   z <- (x - mu) / beta
#   -log(beta) - z - exp(-z)
# }
# 
# true_score_gumbel_1d <- function(x, mu = mu_gumbel, beta = beta_gumbel) {
#   x <- as.numeric(x)
#   z <- (x - mu) / beta
#   matrix((1 - exp(-z)) / beta, ncol = 1)
# }
# 
# res_1d_gumbel_smoke <- run_final_benchmark(
#   sample_sizes = c(50, 200, 500, 1000, 5000),
#   family = "univariate",
#   estimator_specs = estimators_1d,
#   r_sample = r_sample_gumbel_1d,
#   metrics = metrics_all_1d,
#   n_rep = 20,
#   n_test = 2000,
#   true_density = true_density_gumbel_1d,
#   true_logdensity = true_logdensity_gumbel_1d,
#   true_score = true_score_gumbel_1d,
#   grid_1d = grid_1d,
#   seed = 125,
#   verbose = TRUE
# )
# 
# plot_metric_panel(res_1d_gumbel_smoke, metrics = c("negloglik", "ise", "kl", "hellinger2", "score_loss", "score_rmse"))
# plot_metric_panel(res_1d_gumbel_smoke, metrics = c("ise"))
# plot_metric_panel(res_1d_gumbel_smoke, metrics = c("kl"))
# plot_metric_panel(res_1d_gumbel_smoke, metrics = c("hellinger2"))
# plot_metric_panel(res_1d_gumbel_smoke, metrics = c("score_loss"))
# plot_metric_panel(res_1d_gumbel_smoke, metrics = c("score_rmse"))
# 
# 
# 
# # ============================================================
# # (C) Logstic
# # ============================================================
# mu_logistic <- 0
# s_logistic <- 1
# 
# r_sample_logistic_1d <- function(n, mu = mu_logistic, s = s_logistic) {
#   rlogis(n, location = mu, scale = s)
# }
# 
# true_density_logistic_1d <- function(x, mu = mu_logistic, s = s_logistic) {
#   dlogis(x, location = mu, scale = s)
# }
# 
# true_logdensity_logistic_1d <- function(x, mu = mu_logistic, s = s_logistic) {
#   dlogis(x, location = mu, scale = s, log = TRUE)
# }
# 
# true_score_logistic_1d <- function(x, mu = mu_logistic, s = s_logistic) {
#   x <- as.numeric(x)
#   z <- (x - mu) / s
#   matrix((1 - 2 / (1 + exp(z))) / s, ncol = 1)
# }
# 
# res_1d_logistic_smoke <- run_final_benchmark(
#   sample_sizes = c(50, 200, 500, 1000, 5000),
#   family = "univariate",
#   estimator_specs = estimators_1d,
#   r_sample = r_sample_logistic_1d,
#   metrics = metrics_all_1d,
#   n_rep = 10,
#   n_test = 1000,
#   true_density = true_density_logistic_1d,
#   true_logdensity = true_logdensity_logistic_1d,
#   true_score = true_score_logistic_1d,
#   grid_1d = grid_1d,
#   # seed = 125,
#   seed = 177,
#   verbose = TRUE
# )
# 
# plot_metric_panel(res_1d_logistic_smoke, metrics = c("negloglik", "ise", "kl", "hellinger2"))
# 
# # ============================================================
# # DEBUG: punktweise KL-Beiträge eines einzelnen problematischen Runs
# # ============================================================
# 
# debug_one_kl_run <- function(n,
#                              repetition,
#                              estimator_spec,
#                              r_sample,
#                              true_logdensity,
#                              family = "univariate",
#                              n_test = 1000,
#                              seed = NULL,
#                              top_k = 10,
#                              verbose = TRUE) {
#   if (!is.null(seed)) set.seed(seed)
#   
#   # Reproduziere exakt die Schleifenlogik aus run_final_benchmark():
#   # für alle vorherigen Kombinationen Samples "verbrauchen"
#   for (rep in seq_len(repetition)) {
#     x_train <- r_sample(n)
#     if (family == "univariate") x_train <- as.numeric(x_train) else x_train <- as.matrix(x_train)
#     
#     fit <- tryCatch(
#       do.call(
#         fit_estimator_generic,
#         c(
#           list(
#             x = x_train,
#             family = family,
#             method = estimator_spec$method,
#             smoothed = estimator_spec$smoothed
#           ),
#           estimator_spec$fit_args
#         )
#       ),
#       error = function(e) structure(list(error_message = conditionMessage(e)), class = "fit_error")
#     )
#     
#     x_test <- r_sample(n_test)
#     if (family == "univariate") x_test <- as.numeric(x_test) else x_test <- as.matrix(x_test)
#   }
#   
#   if (inherits(fit, "fit_error")) {
#     stop(sprintf("Fit failed in debug run: %s", fit$error_message))
#   }
#   
#   log_hat <- do.call(
#     predict_logdensity_estimator_generic,
#     c(
#       list(
#         newx = x_test,
#         fit = fit,
#         family = family,
#         method = estimator_spec$method
#       ),
#       estimator_spec$density_predict_args %||% list()
#     )
#   )
#   
#   log_true <- as.numeric(true_logdensity(x_test))
#   kl_point <- log_true - log_hat
#   
#   ord_desc <- order(kl_point, decreasing = TRUE, na.last = TRUE)
#   top_idx <- head(ord_desc, top_k)
#   
#   finite_kl <- kl_point[is.finite(kl_point)]
#   
#   trim_alpha <- estimator_spec$density_metric_args$trim_alpha %||% 0.05
#   cutoff <- if (length(finite_kl) > 0L) {
#     as.numeric(stats::quantile(finite_kl, probs = 1 - trim_alpha, na.rm = TRUE, names = FALSE))
#   } else {
#     NA_real_
#   }
#   
#   summary_list <- list(
#     n_test = length(kl_point),
#     n_finite = sum(is.finite(kl_point)),
#     n_bad = sum(!is.finite(kl_point)),
#     bad_share = mean(!is.finite(kl_point)),
#     kl_mean_full = if (length(finite_kl) > 0L) mean(finite_kl) else NA_real_,
#     kl_median = if (length(finite_kl) > 0L) stats::median(finite_kl) else NA_real_,
#     kl_trimmed = if (length(finite_kl) > 0L) mean(finite_kl[finite_kl <= cutoff]) else NA_real_,
#     trim_cutoff = cutoff,
#     top1_share_of_sum = if (length(finite_kl) > 0L && sum(abs(finite_kl)) > 0) {
#       max(abs(finite_kl)) / sum(abs(finite_kl))
#     } else {
#       NA_real_
#     },
#     top5_share_of_sum = if (length(finite_kl) >= 5L && sum(abs(finite_kl)) > 0) {
#       sum(sort(abs(finite_kl), decreasing = TRUE)[1:5]) / sum(abs(finite_kl))
#     } else {
#       NA_real_
#     }
#   )
#   
#   top_df <- data.frame(
#     idx = top_idx,
#     x_test = x_test[top_idx],
#     log_true = log_true[top_idx],
#     log_hat = log_hat[top_idx],
#     kl_point = kl_point[top_idx],
#     is_finite = is.finite(kl_point[top_idx]),
#     stringsAsFactors = FALSE
#   )
#   
#   if (isTRUE(verbose)) {
#     cat("\n===== DEBUG KL RUN =====\n")
#     print(summary_list)
#     cat("\nTop punktweise KL-Beiträge:\n")
#     print(top_df)
#   }
#   
#   invisible(list(
#     summary = summary_list,
#     top = top_df,
#     x_test = x_test,
#     log_true = log_true,
#     log_hat = log_hat,
#     kl_point = kl_point,
#     fit = fit
#   ))
# }
# 
# dbg_kl <- debug_one_kl_run(
#   n = 1000,
#   repetition = 8,
#   estimator_spec = estimators_1d[[1]],
#   r_sample = r_sample_logistic_1d,
#   true_logdensity = true_logdensity_logistic_1d,
#   family = "univariate",
#   n_test = 1000,
#   seed = 177,
#   top_k = 10,
#   verbose = TRUE
# )
# 
# hist(
#   dbg_kl$kl_point[is.finite(dbg_kl$kl_point)],
#   breaks = 80,
#   main = "Punktweise KL-Beiträge",
#   xlab = "log_true(x) - log_hat(x)"
# )
# abline(v = dbg_kl$summary$kl_median, lwd = 2, lty = 2)
# abline(v = dbg_kl$summary$kl_trimmed, lwd = 2, lty = 3)
# abline(v = dbg_kl$summary$kl_mean_full, lwd = 2)
# 
# # ============================================================
# # Iterative 2D smoke tests:
# #   (A) abhängige Gaussian
# #   (B) abhängige multivariate t
# # ============================================================
# 
# # ------------------------------------------------------------
# # (7) Beispiel 2: multivariat  --- angepasst an neue Implementierung
# # ------------------------------------------------------------
# 
# r_sample_norm_2d <- function(n) {
#   cbind(rnorm(n), rnorm(n))
# }
# 
# true_score_norm_2d <- function(x) {
#   x <- as.matrix(x)
#   x
# }
# 
# estimators_2d <- list(
#   list(
#     label = "KDE_Hpi",
#     method = "KDE",
#     smoothed = FALSE,
#     fit_args = list(H_method = "Hpi"),
#     density_predict_args = list(),
#     score_predict_args = list(),
#     # score_metric_args = list(robust = "none")
#     # alternativ:
#     score_metric_args = list(robust = "trim", trim_alpha = 0.05)
#   ),
#   list(
#     label = "MLE_mv",
#     method = "MLE",
#     smoothed = FALSE,
#     fit_args = list(),
#     density_predict_args = list(),
#     score_predict_args = list(),
#     # score_metric_args = list(robust = "none")
#     # alternativ:
#     score_metric_args = list(robust = "trim", trim_alpha = 0.05)
#   ),
#   list(
#     label = "SM_basic",
#     method = "SM",
#     smoothed = FALSE,
#     fit_args = list(
#       m = 2,
#       include_interactions = TRUE,
#       standardize = TRUE,
#       ridge = 1e-6,
#       log_concave = FALSE
#     ),
#     density_predict_args = list(),   # bleibt leer; Dichte-Metriken für mv-SM derzeit nicht unterstützt
#     score_predict_args = list(),
#     # score_metric_args = list(robust = "none")
#     # alternativ:
#     score_metric_args = list(robust = "trim", trim_alpha = 0.05)
#   ),
#   list(
#     label = "SM_grid",
#     method = "SM",
#     smoothed = FALSE,
#     fit_args = list(
#       m = 2,
#       include_interactions = TRUE,
#       standardize = TRUE,
#       ridge = 1e-6,
#       log_concave = TRUE,
#       lc_method = "grid",
#       lc_grid_size = 5L,
#       lc_penalty = 1e4
#     ),
#     density_predict_args = list(),   # bleibt leer; Dichte-Metriken für mv-SM derzeit nicht unterstützt
#     score_predict_args = list(),
#     # score_metric_args = list(robust = "none")
#     # alternativ:
#     score_metric_args = list(robust = "trim", trim_alpha = 0.05)
#   )
# )
# 
# # Achtung:
# # Für multivariates SM werden hier bewusst keine Dichte-Metriken gerechnet.
# res_2d <- run_final_benchmark(
#   sample_sizes = c(100, 200, 500, 1000),
#   family = "multivariate",
#   estimator_specs = estimators_2d,
#   r_sample = r_sample_norm_2d,
#   metrics = c("score_loss", "score_rmse", "logconcavity"),
#   n_rep = 5,
#   n_test = 500,
#   true_score = true_score_norm_2d,
#   seed = 123,
#   verbose = TRUE
# )
# 
# plot_metric_panel(res_2d, metrics = c("score_loss", "score_rmse"))
# plot_metric_panel(res_2d, metrics = c("score_loss"))
# plot_metric_panel(res_2d, metrics = c("score_rmse"))
# 
# # Log-Konkavitätsmetriken ansehen:
# aggregate_final_benchmark(res_2d, "min_hessian_eigenvalue")
# aggregate_final_benchmark(res_2d, "share_violated")
# 
# 
# 
# # 
# # 
# # # ============================================================
# # # (B) 2D abhängige multivariate t mit df = 4
# # # Konstruktion:
# # #   X = Z / sqrt(W/nu),  Z ~ N(0,Sigma), W ~ chi^2_nu
# # #
# # # Score:
# # #   r(x) = ((nu + d) / (nu + x^T Sigma^{-1} x)) * Sigma^{-1} x
# # # mit d = 2
# # # ============================================================
# # 
# # nu_t <- 4
# # d_t <- 2
# # 
# # r_sample_t_2d_dep <- function(n) {
# #   z <- matrix(rnorm(2 * n), ncol = 2) %*% chol_Sigma_2d
# #   w <- rchisq(n, df = nu_t)
# #   z / sqrt(w / nu_t)
# # }
# # 
# # true_score_t_2d_dep <- function(x) {
# #   x <- as.matrix(x)
# #   qf <- rowSums((x %*% Sigma_inv_2d) * x)
# #   fac <- (nu_t + d_t) / (nu_t + qf)
# #   (x %*% Sigma_inv_2d) * fac
# # }
# # 
# # res_2d_t_smoke <- run_final_benchmark(
# #   sample_sizes = c(40, 80),
# #   family = "multivariate",
# #   estimator_specs = estimators_2d,
# #   r_sample = r_sample_t_2d_dep,
# #   metrics = c("score_loss", "score_rmse", "logconcavity"),
# #   n_rep = 5,
# #   n_test = 500,
# #   true_score = true_score_t_2d_dep,
# #   seed = 126,
# #   verbose = TRUE
# # )
# # 
# # plot_metric_panel(res_2d_t_smoke, metrics = c("score_loss", "score_rmse"))
# # 
# # aggregate_final_benchmark(res_2d_t_smoke, "score_loss")
# # aggregate_final_benchmark(res_2d_t_smoke, "share_violated")
# 
