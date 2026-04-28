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

metric_to_output_columns <- function(metrics,
                                     density_metric_args = list(),
                                     score_metric_args = list()) {
  metrics <- unique(metrics)
  out <- character(0)

  for (m in metrics) {
    cfgs <- if (m %in% c("negloglik", "kl")) {
      normalize_metric_configurations(m, density_metric_args)
    } else if (m == "score_loss") {
      normalize_metric_configurations(m, score_metric_args)
    } else {
      list(default = list())
    }

    n_cfg <- length(cfgs)
    metric_names <- vapply(
      names(cfgs),
      function(cfg_name) metric_variant_name(m, cfg_name, n_cfg),
      character(1)
    )

    out <- c(out, metric_names)
    out <- c(
      out,
      paste0(metric_names, "_na_share"),
      paste0(metric_names, "_tail_share"),
      paste0(metric_names, "_outlier_dominated")
    )
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
                                            density_metric_args = list(),
                                            score_metric_args = list(),
                                            estimator_label = method,
                                            n = NA_integer_,
                                            verbose = TRUE) {
  metric_columns <- metric_to_output_columns(
    metrics,
    density_metric_args = density_metric_args,
    score_metric_args = score_metric_args
  )
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
                                     run_seed = NULL,
                                     verbose = TRUE) {
  family <- match.arg(family)
  metric_columns <- metric_columns %||% metric_to_output_columns(
    metrics,
    density_metric_args = estimator_spec$density_metric_args %||% list(),
    score_metric_args = estimator_spec$score_metric_args %||% list()
  )

  if (!is.null(run_seed)) set.seed(run_seed)

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
    run_seed = if (is.null(run_seed)) NA_integer_ else as.integer(run_seed),
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
  metric_columns <- unique(unlist(lapply(estimator_specs, function(spec) {
    metric_to_output_columns(
      metrics,
      density_metric_args = spec$density_metric_args %||% list(),
      score_metric_args = spec$score_metric_args %||% list()
    )
  }), use.names = FALSE))

  for (spec in estimator_specs) {
    if (verbose) message("Method: ", spec$label)

    for (n in sample_sizes) {
      if (verbose) message("  n = ", n)

      for (rep in seq_len(n_rep)) {
        if (verbose) message("    repetition ", rep, "/", n_rep)

        run_seed <- sample.int(.Machine$integer.max, size = 1L)

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
          run_seed = run_seed,
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
    "n", "repetition", "run_seed", "method_label", "method",
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
        n_test = n_test,
        benchmark_seed = seed
      ),
      estimator_specs = estimator_specs,
      benchmark_inputs = list(
        r_sample = r_sample,
        true_density = true_density,
        true_logdensity = true_logdensity,
        true_score = true_score,
        grid_1d = grid_1d
      )
    ),
    class = "final_benchmark"
  )
}

normalize_method_filter <- function(keep_method_labels = NULL,
                                   drop_method_labels = NULL,
                                   keep_methods = NULL,
                                   drop_methods = NULL) {
  if (!is.null(keep_method_labels) && !is.null(drop_method_labels)) {
    stop("Use either keep_method_labels or drop_method_labels, not both.")
  }
  if (!is.null(keep_methods) && !is.null(drop_methods)) {
    stop("Use either keep_methods or drop_methods, not both.")
  }
  list(
    keep_method_labels = if (is.null(keep_method_labels)) NULL else unique(as.character(keep_method_labels)),
    drop_method_labels = if (is.null(drop_method_labels)) NULL else unique(as.character(drop_method_labels)),
    keep_methods = if (is.null(keep_methods)) NULL else unique(as.character(keep_methods)),
    drop_methods = if (is.null(drop_methods)) NULL else unique(as.character(drop_methods))
  )
}

filter_benchmark_df_by_method <- function(df,
                                          keep_method_labels = NULL,
                                          drop_method_labels = NULL,
                                          keep_methods = NULL,
                                          drop_methods = NULL) {
  filt <- normalize_method_filter(
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )
  out <- df
  if (!is.null(filt$keep_method_labels)) out <- out[out$method_label %in% filt$keep_method_labels, , drop = FALSE]
  if (!is.null(filt$drop_method_labels)) out <- out[!out$method_label %in% filt$drop_method_labels, , drop = FALSE]
  if (!is.null(filt$keep_methods)) out <- out[out$method %in% filt$keep_methods, , drop = FALSE]
  if (!is.null(filt$drop_methods)) out <- out[!out$method %in% filt$drop_methods, , drop = FALSE]
  rownames(out) <- NULL
  out
}

robust_across_runs_center <- function(x,
                                      center = c("mean", "median", "trim"),
                                      trim_alpha = 0.05) {
  center <- match.arg(center)
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) return(NA_real_)
  
  if (center == "mean") {
    return(mean(x))
  }
  
  if (center == "median") {
    return(stats::median(x))
  }
  
  mean(x, trim = trim_alpha)
}

aggregate_final_benchmark <- function(obj,
                                      metric,
                                      drop_all_na = FALSE,
                                      keep_n = NULL,
                                      drop_n = NULL,
                                      keep_method_labels = NULL,
                                      drop_method_labels = NULL,
                                      keep_methods = NULL,
                                      drop_methods = NULL,
                                      across_runs_center = c("mean", "median", "trim"),
                                      trim_alpha = 0.05) {
  across_runs_center <- match.arg(across_runs_center)
  
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")
  df <- filter_benchmark_df_by_n(obj$raw, keep_n = keep_n, drop_n = drop_n)
  df <- filter_benchmark_df_by_method(
    df,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )
  
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
      trim = robust_across_runs_center(x, center = "trim", trim_alpha = trim_alpha),
      selected = robust_across_runs_center(x, center = across_runs_center, trim_alpha = trim_alpha),
      sd = safe_sd(x),
      q25 = safe_quantile(x, 0.25),
      q75 = safe_quantile(x, 0.75),
      across_runs_center = across_runs_center,
      trim_alpha = trim_alpha,
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

subset_final_benchmark <- function(obj,
                                   keep_n = NULL,
                                   drop_n = NULL,
                                   keep_method_labels = NULL,
                                   drop_method_labels = NULL,
                                   keep_methods = NULL,
                                   drop_methods = NULL,
                                   update_settings = TRUE) {
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")
  out <- obj
  out$raw <- filter_benchmark_df_by_n(out$raw, keep_n = keep_n, drop_n = drop_n)
  out$raw <- filter_benchmark_df_by_method(
    out$raw,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )

  if (isTRUE(update_settings)) {
    out$settings$sample_sizes <- sort(unique(out$raw$n))
  }
  out
}


debug_benchmark_outliers <- function(obj,
                                     metric_pattern = "^(negloglik|kl|score_loss)",
                                     top_n = 10,
                                     group_cols = c("method_label", "n"),
                                     min_group_size = 3L) {
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")

  df <- obj$raw
  metric_cols <- grep(metric_pattern, names(df), value = TRUE)
  metric_cols <- metric_cols[!grepl("_(na_share|tail_share|outlier_dominated)$", metric_cols)]

  if (length(metric_cols) == 0L) {
    stop("No metric columns matched metric_pattern.")
  }

  rows <- list()
  counter <- 1L

  for (metric in metric_cols) {
    split_key <- interaction(df[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)

    for (dd in split(df, split_key)) {
      x <- as.numeric(dd[[metric]])
      keep <- is.finite(x)
      if (sum(keep) < min_group_size) next

      center <- stats::median(x[keep])
      mad_val <- stats::mad(x[keep], center = center, constant = 1, na.rm = TRUE)
      scale_val <- max(mad_val, 1e-12)
      abs_dev <- abs(x - center)
      robust_z <- abs_dev / scale_val
      rel_dev <- abs_dev / pmax(abs(center), 1e-12)

      ord <- order(robust_z, decreasing = TRUE, na.last = NA)
      take <- head(ord, top_n)

      for (idx in take) {
        rows[[counter]] <- data.frame(
          metric = metric,
          method_label = dd$method_label[idx],
          method = dd$method[idx],
          n = dd$n[idx],
          repetition = dd$repetition[idx],
          run_seed = dd$run_seed[idx],
          value = x[idx],
          group_center = center,
          abs_dev = abs_dev[idx],
          rel_dev = rel_dev[idx],
          robust_z = robust_z[idx],
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }
  }

  out <- bind_rows_fill_base(rows)
  if (nrow(out) == 0L) return(out)
  out <- out[order(out$robust_z, out$abs_dev, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

replay_benchmark_run <- function(obj,
                                 method_label,
                                 n,
                                 repetition = NULL,
                                 run_seed = NULL,
                                 metrics = NULL,
                                 n_test = NULL,
                                 verbose = TRUE) {
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")

  df <- obj$raw
  sel <- df$method_label == method_label & df$n == n
  if (!is.null(repetition)) sel <- sel & df$repetition == repetition
  if (!is.null(run_seed)) sel <- sel & df$run_seed == run_seed

  hits <- df[sel, , drop = FALSE]
  if (nrow(hits) == 0L) stop("No matching run found.")
  if (nrow(hits) > 1L) stop("Selection is ambiguous. Please specify repetition or run_seed.")

  row <- hits[1, , drop = FALSE]
  spec_idx <- which(vapply(obj$estimator_specs, function(spec) identical(spec$label, row$method_label), logical(1)))
  if (length(spec_idx) != 1L) stop("Could not uniquely resolve estimator_spec from method_label.")
  spec <- obj$estimator_specs[[spec_idx]]

  metrics <- metrics %||% obj$settings$metrics
  n_test <- n_test %||% obj$settings$n_test
  family <- obj$settings$family

  set.seed(row$run_seed)
  x_train <- obj$benchmark_inputs$r_sample(row$n)
  if (family == "univariate") x_train <- as.numeric(x_train) else x_train <- as.matrix(x_train)

  fit <- do.call(
    fit_estimator_generic,
    c(
      list(
        x = x_train,
        family = family,
        method = spec$method,
        smoothed = spec$smoothed
      ),
      spec$fit_args
    )
  )

  x_test <- obj$benchmark_inputs$r_sample(n_test)
  if (family == "univariate") x_test <- as.numeric(x_test) else x_test <- as.matrix(x_test)

  aggregate_metrics <- evaluate_requested_metrics(
    metrics = metrics,
    x_test = x_test,
    fit = fit,
    family = family,
    method = spec$method,
    true_density = obj$benchmark_inputs$true_density,
    true_logdensity = obj$benchmark_inputs$true_logdensity,
    true_score = obj$benchmark_inputs$true_score,
    grid_1d = obj$benchmark_inputs$grid_1d,
    density_predict_args = spec$density_predict_args %||% list(),
    score_predict_args = spec$score_predict_args %||% list(),
    density_metric_args = spec$density_metric_args %||% list(),
    score_metric_args = spec$score_metric_args %||% list()
  )

  pointwise <- data.frame(stringsAsFactors = FALSE)
  if (family == "univariate") {
    pointwise$x_test <- as.numeric(x_test)
  } else {
    x_mat <- as.matrix(x_test)
    pointwise <- as.data.frame(x_mat)
    names(pointwise) <- paste0("x", seq_len(ncol(x_mat)))
  }

  if (any(metrics %in% c("negloglik", "kl"))) {
    negloglik_pt <- compute_pointwise_negloglik(
      x_test = x_test,
      fit = fit,
      family = family,
      method = spec$method,
      predict_args = spec$density_predict_args %||% list()
    )
    pointwise$negloglik_point <- negloglik_pt$losses
    pointwise$log_hat <- negloglik_pt$log_hat

    if (!is.null(obj$benchmark_inputs$true_logdensity)) {
      kl_pt <- compute_pointwise_kl(
        x_test = x_test,
        fit = fit,
        family = family,
        method = spec$method,
        true_logdensity = obj$benchmark_inputs$true_logdensity,
        predict_args = spec$density_predict_args %||% list()
      )
      pointwise$kl_point <- kl_pt$losses
      pointwise$log_true <- kl_pt$log_true
    }
  }

  if (any(metrics %in% c("score_loss")) && !is.null(obj$benchmark_inputs$true_score)) {
    score_pt <- compute_pointwise_score_loss(
      x_test = x_test,
      fit = fit,
      family = family,
      method = spec$method,
      true_score = obj$benchmark_inputs$true_score,
      predict_args = spec$score_predict_args %||% list()
    )
    score_df <- data.frame(score_loss_point = NA_real_)
    score_df$score_loss_point <- NA_real_
    score_df <- score_df[rep(1, nrow(pointwise)), , drop = FALSE]
    score_df$score_loss_point[] <- NA_real_
    score_df$score_loss_point[score_pt$keep] <- score_pt$losses
    pointwise$score_loss_point <- score_df$score_loss_point

    score_hat <- score_pt$score_hat
    score_true <- score_pt$score_true
    if (ncol(score_hat) == 1L) {
      pointwise$score_hat <- NA_real_
      pointwise$score_true <- NA_real_
      pointwise$score_hat[score_pt$keep] <- score_hat[, 1]
      pointwise$score_true[score_pt$keep] <- score_true[, 1]
    } else {
      for (j in seq_len(ncol(score_hat))) {
        pointwise[[paste0("score_hat_", j)]] <- NA_real_
        pointwise[[paste0("score_true_", j)]] <- NA_real_
        pointwise[[paste0("score_hat_", j)]][score_pt$keep] <- score_hat[, j]
        pointwise[[paste0("score_true_", j)]][score_pt$keep] <- score_true[, j]
      }
    }
  }

  list(
    selected_run = row,
    estimator_spec = spec,
    fit = fit,
    x_train = x_train,
    x_test = x_test,
    aggregate_metrics = aggregate_metrics,
    pointwise = pointwise
  )
}

plot_final_benchmark <- function(obj,
                                 metric,
                                 center = c("mean", "median"),
                                 interval = c("sd", "iqr", "none"),
                                 keep_n = NULL,
                                 drop_n = NULL,
                                 keep_method_labels = NULL,
                                 drop_method_labels = NULL,
                                 keep_methods = NULL,
                                 drop_methods = NULL,
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
    drop_n = drop_n,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
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
      subtitle = paste(na.omit(c(
        if (!is.null(keep_n)) sprintf("Shown n: %s", paste(sort(unique(keep_n)), collapse = ", ")) else NULL,
        if (!is.null(drop_n)) sprintf("Dropped n: %s", paste(sort(unique(drop_n)), collapse = ", ")) else NULL,
        if (!is.null(keep_method_labels)) sprintf("Shown labels: %s", paste(unique(keep_method_labels), collapse = ", ")) else NULL,
        if (!is.null(drop_method_labels)) sprintf("Dropped labels: %s", paste(unique(drop_method_labels), collapse = ", ")) else NULL,
        if (!is.null(keep_methods)) sprintf("Shown methods: %s", paste(unique(keep_methods), collapse = ", ")) else NULL,
        if (!is.null(drop_methods)) sprintf("Dropped methods: %s", paste(unique(drop_methods), collapse = ", ")) else NULL
      )), collapse = " | ")
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
