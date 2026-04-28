# ============================================================
# Tests_Clean.R
# Modulare Benchmark-Funktionen für univariate / multivariate Tests
# Fokus auf die aktuell verwendeten Metriken:
#   - negloglik
#   - kl
#   - score_loss
# Zusätzlich:
#   - einheitliche Fit- und Inferenzzeiten
#   - einheitliche numerische Diagnostik
# ============================================================

source("helper_functions.R")
source("KDE.R")
source("LogConcaveMLE.R")
source("Univariate_Polynomial_Score_Matching_1.0.R")
source("Multivariate_Pairwise_Polynomial_Score_Matching.R")
source("Evaluation_Metrics_Clean.R")

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

metric_to_output_columns <- function(metrics) {
  unique(metrics)
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
  row[, schema_names, drop = FALSE]
}

bind_rows_fill_base <- function(rows) {
  if (length(rows) == 0L) return(data.frame())
  schema <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows2 <- lapply(rows, coerce_row_to_schema, schema_names = schema)
  out <- do.call(rbind, rows2)
  rownames(out) <- NULL
  out
}

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
                                            score_metric_args = list(),
                                            estimator_label = method,
                                            n = NA_integer_,
                                            verbose = TRUE) {
  metric_columns <- metric_to_output_columns(metrics)
  empty_out <- make_empty_metric_list(metric_columns)

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
      "Joint metric evaluation failed for %s at n=%s. Falling back to metric-wise evaluation. Reason: %s",
      estimator_label, n, conditionMessage(full_res)
    ))
  }

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

  fit_time <- system.time({
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
    fit_time_sec = as.numeric(fit_time["elapsed"]),
    density_inference_time_sec = NA_real_,
    score_inference_time_sec = NA_real_,
    total_inference_time_sec = NA_real_,
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
    error = function(e) list(
      success = TRUE,
      status = paste("diagnostics_error:", conditionMessage(e)),
      iterations = NA_real_,
      objective_value = NA_real_,
      kappa_raw = NA_real_,
      kappa_reg = NA_real_,
      rcond_raw = NA_real_,
      rcond_reg = NA_real_,
      eigmin_raw = NA_real_,
      eigmin_reg = NA_real_
    )
  )

  density_metrics_requested <- any(metrics %in% c("negloglik", "kl"))
  score_metrics_requested <- any(metrics %in% c("score_loss"))

  density_time <- if (density_metrics_requested) {
    system.time({
      invisible(do.call(
        predict_logdensity_estimator_generic,
        c(
          list(newx = x_test, fit = fit, family = family, method = estimator_spec$method),
          estimator_spec$density_predict_args %||% list()
        )
      ))
    })
  } else {
    c(elapsed = NA_real_)
  }

  score_time <- if (score_metrics_requested) {
    system.time({
      invisible(do.call(
        predict_score_estimator_generic,
        c(
          list(newx = x_test, fit = fit, family = family, method = estimator_spec$method),
          estimator_spec$score_predict_args %||% list()
        )
      ))
    })
  } else {
    c(elapsed = NA_real_)
  }

  base_row$density_inference_time_sec <- as.numeric(density_time["elapsed"])
  base_row$score_inference_time_sec <- as.numeric(score_time["elapsed"])
  base_row$total_inference_time_sec <- sum(
    c(base_row$density_inference_time_sec, base_row$score_inference_time_sec),
    na.rm = TRUE
  )
  if (!is.finite(base_row$total_inference_time_sec)) base_row$total_inference_time_sec <- NA_real_

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

  for (nm in names(metric_values)) base_row[[nm]] <- metric_values[[nm]]

  base_row
}

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
    "n", "repetition", "method_label", "method",
    "fit_time_sec", "density_inference_time_sec", "score_inference_time_sec", "total_inference_time_sec",
    "success", "status", "iterations", "objective_value",
    "kappa_raw", "kappa_reg", "rcond_raw", "rcond_reg", "eigmin_raw", "eigmin_reg",
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

aggregate_final_benchmark <- function(obj, metric, drop_all_na = FALSE, keep_n = NULL, drop_n = NULL) {
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")
  df <- filter_benchmark_df_by_n(obj$raw, keep_n = keep_n, drop_n = drop_n)

  if (!metric %in% names(df)) {
    available <- intersect(
      c(
        "negloglik", "kl", "score_loss",
        "fit_time_sec", "density_inference_time_sec", "score_inference_time_sec", "total_inference_time_sec",
        "kappa_raw", "kappa_reg", "rcond_raw", "rcond_reg", "eigmin_raw", "eigmin_reg"
      ),
      names(df)
    )
    stop(sprintf(
      "Metric '%s' not found in benchmark output. Available metrics: %s",
      metric,
      paste(available, collapse = ", ")
    ))
  }

  split_key <- interaction(df$method_label, df$n, drop = TRUE)
  agg_list <- lapply(split(df, split_key), function(dd) {
    x <- dd[[metric]]
    data.frame(
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
  })

  agg <- do.call(rbind, agg_list)
  rownames(agg) <- NULL
  if (isTRUE(drop_all_na)) agg <- agg[agg$n_non_na > 0L, , drop = FALSE]
  agg
}


normalize_n_filter <- function(keep_n = NULL, drop_n = NULL) {
  if (!is.null(keep_n) && !is.null(drop_n)) {
    stop("Use either keep_n or drop_n, not both.")
  }
  list(
    keep_n = if (is.null(keep_n)) NULL else unique(as.numeric(keep_n)),
    drop_n = if (is.null(drop_n)) NULL else unique(as.numeric(drop_n))
  )
}

filter_benchmark_df_by_n <- function(df, keep_n = NULL, drop_n = NULL) {
  filt <- normalize_n_filter(keep_n = keep_n, drop_n = drop_n)
  out <- df
  if (!is.null(filt$keep_n)) out <- out[out$n %in% filt$keep_n, , drop = FALSE]
  if (!is.null(filt$drop_n)) out <- out[!out$n %in% filt$drop_n, , drop = FALSE]
  rownames(out) <- NULL
  out
}

subset_final_benchmark <- function(obj, keep_n = NULL, drop_n = NULL, update_settings = TRUE) {
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")
  out <- obj
  out$raw <- filter_benchmark_df_by_n(out$raw, keep_n = keep_n, drop_n = drop_n)

  if (isTRUE(update_settings)) {
    out$settings$sample_sizes <- sort(unique(out$raw$n))
  }
  out
}


plot_final_benchmark <- function(obj,
                                 metric,
                                 center = c("mean", "median"),
                                 interval = c("sd", "iqr", "none"),
                                 keep_n = NULL,
                                 drop_n = NULL,
                                 drop_all_na = TRUE,
                                 log_x = FALSE,
                                 log_y = FALSE) {
  center <- match.arg(center)
  interval <- match.arg(interval)

  agg <- aggregate_final_benchmark(
    obj = obj,
    metric = metric,
    drop_all_na = drop_all_na,
    keep_n = keep_n,
    drop_n = drop_n
  )

  if (nrow(agg) == 0L) {
    stop("No rows left after applying the sample-size filter.")
  }

  y_col <- if (center == "mean") "mean" else "median"
  agg$y <- agg[[y_col]]

  if (interval == "sd") {
    agg$ymin <- agg$y - agg$sd
    agg$ymax <- agg$y + agg$sd
  } else if (interval == "iqr") {
    agg$ymin <- agg$q25
    agg$ymax <- agg$q75
  } else {
    agg$ymin <- NA_real_
    agg$ymax <- NA_real_
  }

  p <- ggplot2::ggplot(
    agg,
    ggplot2::aes(x = n, y = y, color = method_label, group = method_label)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Sample size n",
      y = metric,
      color = "Method",
      title = sprintf("%s by sample size", metric),
      subtitle = if (!is.null(keep_n)) {
        sprintf("Shown n: %s", paste(sort(unique(keep_n)), collapse = ", "))
      } else if (!is.null(drop_n)) {
        sprintf("Dropped n: %s", paste(sort(unique(drop_n)), collapse = ", "))
      } else {
        NULL
      }
    ) +
    ggplot2::theme_minimal()

  if (interval != "none") {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ymin, ymax = ymax, fill = method_label),
      alpha = 0.15,
      colour = NA,
      show.legend = FALSE
    )
  }

  if (isTRUE(log_x)) p <- p + ggplot2::scale_x_log10()
  if (isTRUE(log_y)) p <- p + ggplot2::scale_y_log10()

  p
}
