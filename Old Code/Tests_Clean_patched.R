# ============================================================
# Tests_Clean_patched.R
# Modular benchmark functions for univariate / multivariate tests
# Focus metrics:
#   - kl
#   - score_loss
# Optional metric variants are generated automatically from
# central_trim and robust_trim settings.
# ============================================================

source("helper_functions.R")
source("KDE.R")
source("LogConcaveMLE.R")
source("Univariate_Polynomial_Score_Matching_1.0_Kopie.R")
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
  unique(build_metric_names(
    metrics = metrics,
    density_metric_args = density_metric_args,
    score_metric_args = score_metric_args
  ))
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
    condition_number = NA_real_,
    normalization_ok = NA,
    normalization_suspect = NA,
    normalization_suspect_0_5 = NA,
    normalization_suspect_1_0 = NA,
    normalization_loghat_finite_share = NA_real_,
    normalization_median_kl_shift = NA_real_,
    stringsAsFactors = FALSE
  )
  for (nm in metric_columns) base_row[[nm]] <- NA_real_

  if (inherits(fit, "fit_error")) {
    base_row$status <- fit$error_message
    return(base_row)
  }

  x_test <- r_sample(n_test)
  if (family == "univariate") x_test <- as.numeric(x_test) else x_test <- as.matrix(x_test)
  
  if (family == "univariate") {
    train_min <- min(x_train, na.rm = TRUE)
    train_max <- max(x_train, na.rm = TRUE)
    test_min  <- min(x_test, na.rm = TRUE)
    test_max  <- max(x_test, na.rm = TRUE)
    
    base_row$left_gap  <- max(0, train_min - test_min)
    base_row$right_gap <- max(0, test_max - train_max)
  }

  diags <- tryCatch(
    extract_fit_diagnostics(fit),
    error = function(e) list(
      success = TRUE,
      status = paste("diagnostics_error:", conditionMessage(e)),
      iterations = NA_real_,
      objective_value = NA_real_,
      condition_number = NA_real_
    )
  )

  density_metrics_requested <- any(metrics %in% c("kl"))
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
  } else c(elapsed = NA_real_)

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
  } else c(elapsed = NA_real_)

  base_row$density_inference_time_sec <- as.numeric(density_time["elapsed"])
  base_row$score_inference_time_sec <- as.numeric(score_time["elapsed"])
  base_row$total_inference_time_sec <- sum(
    c(base_row$density_inference_time_sec, base_row$score_inference_time_sec),
    na.rm = TRUE
  )
  if (!is.finite(base_row$total_inference_time_sec)) base_row$total_inference_time_sec <- NA_real_

  if (density_metrics_requested) {
    dens_diag <- tryCatch(
      extract_density_diagnostic(
        fit = fit,
        family = family,
        method = estimator_spec$method,
        predict_args = estimator_spec$density_predict_args %||% list()
      ),
      error = function(e) list(normalization_ok = FALSE, normalization_message = conditionMessage(e))
    )
    
    base_row$normalization_ok <- dens_diag$normalization_ok
    
    # Neue, separate Diagnose:
    base_row$normalization_suspect <- NA
    base_row$normalization_loghat_finite_share <- NA_real_
    base_row$normalization_median_kl_shift <- NA_real_
    
    if (identical(family, "univariate") &&
        identical(estimator_spec$method, "SM") &&
        !is.null(true_logdensity)) {
      
      log_hat <- tryCatch(
        do.call(
          predict_logdensity_estimator_generic,
          c(
            list(newx = x_test, fit = fit, family = family, method = estimator_spec$method),
            estimator_spec$density_predict_args %||% list()
          )
        ),
        error = function(e) rep(NA_real_, length(x_test))
      )
      
      log_hat <- as.numeric(log_hat)
      log_true <- as.numeric(true_logdensity(x_test))
      
      finite_share <- mean(is.finite(log_hat))
      kl_point <- log_true - log_hat
      median_kl_shift <- if (any(is.finite(kl_point))) {
        stats::median(kl_point[is.finite(kl_point)])
      } else {
        Inf
      }
      
      base_row$normalization_loghat_finite_share <- finite_share
      base_row$normalization_median_kl_shift <- median_kl_shift
      
      kl_point <- kl_point[is.finite(kl_point)]
      mean_kl <- mean(kl_point)
      sd_kl   <- sd(kl_point)
      
      constant_shift_flag_0_5 <-
        length(kl_point) > 10 &&
        is.finite(mean_kl) &&
        is.finite(sd_kl) &&
        abs(mean_kl) > 0.1 &&
        (sd_kl / abs(mean_kl)) < 0.1
      
      constant_shift_flag_1_0 <-
        length(kl_point) > 10 &&
        is.finite(mean_kl) &&
        is.finite(sd_kl) &&
        abs(mean_kl) > 1.0 &&
        (sd_kl / abs(mean_kl)) < 0.1
      
      base_row$normalization_suspect_0_5 <-
        isFALSE(dens_diag$normalization_ok) ||
        (!is.na(finite_share) && finite_share < 0.99) ||
        constant_shift_flag_0_5
      
      base_row$normalization_suspect_1_0 <-
        isFALSE(dens_diag$normalization_ok) ||
        (!is.na(finite_share) && finite_share < 0.99) ||
        constant_shift_flag_1_0
      
      # Standard bleibt 0.5
      base_row$normalization_suspect <- base_row$normalization_suspect_0_5
    }
  }

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
  base_row$condition_number <- diags$condition_number %||% NA_real_
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
                                truth_name = NULL,
                                seed = NULL,
                                verbose = TRUE,
                                save = FALSE,
                                save_dir = ".",
                                save_name = NULL) {
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
    "success", "status", "iterations", "objective_value", "condition_number",
    "left_gap", "right_gap",
    "normalization_ok", "normalization_suspect",
    "normalization_suspect_0_1", "normalization_suspect_1_0",
    "normalization_loghat_finite_share", "normalization_median_kl_shift",
    metric_columns
  )
  other_cols <- setdiff(names(raw), preferred_order)
  raw <- raw[, c(intersect(preferred_order, names(raw)), other_cols), drop = FALSE]
  rownames(raw) <- NULL

  obj <- structure(
    list(
      raw = raw,
      settings = list(
        sample_sizes = sample_sizes,
        family = family,
        truth_name = truth_name,
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

  if (isTRUE(save)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    save_name <- save_name %||% make_benchmark_filename(
      kind = "final",
      truth_name = truth_name,
      estimator_specs = estimator_specs,
      family = family
    )
    save_path <- file.path(save_dir, save_name)
    saveRDS(obj, save_path)
    obj$settings$save_path <- save_path
  }

  obj
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

normalize_n_filter <- function(keep_n = NULL, drop_n = NULL) {
  if (!is.null(keep_n) && !is.null(drop_n)) stop("Use either keep_n or drop_n, not both.")
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

aggregate_final_benchmark <- function(obj,
                                      metric,
                                      drop_all_na = FALSE,
                                      keep_n = NULL,
                                      drop_n = NULL,
                                      keep_method_labels = NULL,
                                      drop_method_labels = NULL,
                                      keep_methods = NULL,
                                      drop_methods = NULL,
                                      across_runs_center = c("mean", "median", "sd"),
                                      exclude_normalization_suspect = FALSE) {
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
  
  if (isTRUE(exclude_normalization_suspect)) {
    if (!"normalization_suspect" %in% names(df)) {
      stop("Column 'normalization_suspect' not found in benchmark output.")
    }
    df <- df[!(df$method %in% "SM" & df$normalization_suspect %in% TRUE), , drop = FALSE]
  }
  
  if (!metric %in% names(df)) {
    available <- intersect(c(
      obj$settings$metric_columns,
      "fit_time_sec", "density_inference_time_sec", "score_inference_time_sec", "total_inference_time_sec",
      "condition_number"
    ), names(df))
    stop(sprintf(
      "Metric '%s' not found in benchmark output. Available metrics: %s",
      metric,
      paste(available, collapse = ", ")
    ))
  }
  
  if (nrow(df) == 0L) {
    return(data.frame())
  }
  
  split_key <- interaction(df$method_label, df$n, drop = TRUE)
  
  agg_list <- lapply(split(df, split_key), function(dd) {
    x <- as.numeric(dd[[metric]])
    x_finite <- x[is.finite(x)]
    
    q1 <- safe_quantile(x, 0.25)
    q3 <- safe_quantile(x, 0.75)
    iqr <- if (is.finite(q1) && is.finite(q3)) q3 - q1 else NA_real_
    
    outlier_share <- if (length(x_finite) == 0L || !is.finite(iqr)) {
      NA_real_
    } else if (iqr == 0) {
      mean(x_finite != stats::median(x_finite))
    } else {
      lo <- q1 - 1.5 * iqr
      hi <- q3 + 1.5 * iqr
      mean(x_finite < lo | x_finite > hi)
    }
    
    failure_rate <- mean(!(dd$success) | !is.finite(x), na.rm = TRUE)
    
    normalization_failure_rate <- if ("normalization_ok" %in% names(dd)) {
      mean(dd$normalization_ok %in% FALSE, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    normalization_suspect_rate <- if ("normalization_suspect" %in% names(dd)) {
      mean(dd$normalization_suspect %in% TRUE, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    data.frame(
      method_label = dd$method_label[1],
      method = dd$method[1],
      n = dd$n[1],
      n_non_na = sum(is.finite(x)),
      mean = safe_mean(x),
      median = safe_median(x),
      selected = switch(
        across_runs_center,
        mean = safe_mean(x),
        median = safe_median(x),
        sd = safe_sd(x)
      ),
      q25 = q1,
      q75 = q3,
      iqr = iqr,
      sd = safe_sd(x),
      across_runs_center = across_runs_center,
      outlier_run_share = outlier_share,
      failure_rate = failure_rate,
      normalization_failure_rate = normalization_failure_rate,
      normalization_suspect_rate = normalization_suspect_rate,
      stringsAsFactors = FALSE
    )
  })
  
  agg <- do.call(rbind, agg_list)
  rownames(agg) <- NULL
  
  if (isTRUE(drop_all_na)) {
    agg <- agg[agg$n_non_na > 0L, , drop = FALSE]
  }
  
  agg
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
  if (isTRUE(update_settings)) out$settings$sample_sizes <- sort(unique(out$raw$n))
  out
}

debug_benchmark_outliers <- function(obj,
                                     metric_pattern = "^(kl|score_loss)",
                                     top_n = 10,
                                     group_cols = c("method_label", "n"),
                                     min_group_size = 3L) {
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")

  df <- obj$raw
  metric_cols <- grep(metric_pattern, names(df), value = TRUE)
  if (length(metric_cols) == 0L) stop("No metric columns matched metric_pattern.")

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
          condition_number = dd$condition_number[idx],
          normalization_suspect = dd$normalization_suspect[idx],
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
      list(x = x_train, family = family, method = spec$method, smoothed = spec$smoothed),
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

  if (family == "univariate") {
    pointwise <- data.frame(
      x_test = as.numeric(x_test),
      stringsAsFactors = FALSE
    )
  } else {
    x_mat <- as.matrix(x_test)
    pointwise <- as.data.frame(x_mat, stringsAsFactors = FALSE)
    names(pointwise) <- paste0("x", seq_len(ncol(x_mat)))
  }

  if ("kl" %in% metrics) {
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
    pointwise$log_hat <- kl_pt$log_hat
    pointwise$central_keep_kl <- make_central_mask(
      x_test,
      family = family,
      central_trim = spec$density_metric_args$central_trim %||% NULL
    )
  }

  if ("score_loss" %in% metrics && !is.null(obj$benchmark_inputs$true_score)) {
    score_pt <- compute_pointwise_score_loss(
      x_test = x_test,
      fit = fit,
      family = family,
      method = spec$method,
      true_score = obj$benchmark_inputs$true_score,
      predict_args = spec$score_predict_args %||% list()
    )
    tmp <- rep(NA_real_, if (family == "univariate") length(x_test) else nrow(as.matrix(x_test)))
    tmp[score_pt$keep] <- score_pt$losses
    pointwise$score_loss_point <- tmp
    tmp_keep <- rep(FALSE, length(tmp))
    tmp_keep[score_pt$keep] <- make_central_mask(
      if (family == "univariate") x_test[score_pt$keep] else as.matrix(x_test)[score_pt$keep, , drop = FALSE],
      family = family,
      central_trim = spec$score_metric_args$central_trim %||% NULL
    )
    pointwise$central_keep_score <- tmp_keep
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
                                 center = c("mean", "median", "sd"),
                                 interval = c("iqr", "none"),
                                 keep_n = NULL,
                                 drop_n = NULL,
                                 keep_method_labels = NULL,
                                 drop_method_labels = NULL,
                                 keep_methods = NULL,
                                 drop_methods = NULL,
                                 drop_all_na = TRUE,
                                 log_x = FALSE,
                                 log_y = FALSE,
                                 exclude_normalization_suspect = FALSE) {
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
    drop_methods = drop_methods,
    across_runs_center = center,
    exclude_normalization_suspect = exclude_normalization_suspect
  )

  if (nrow(agg) == 0L) stop("No rows left after applying the sample-size filter.")
  y_col <- switch(
    center,
    mean = "mean",
    median = "median",
    sd = "sd"
  )
  agg$y <- agg[[y_col]]
  if (interval == "iqr" && center != "sd") {
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
      title = sprintf("%s (%s) by sample size", metric, center),
      subtitle = paste(na.omit(c(
        sprintf("Center: %s", center),
        sprintf("Exclude normalization suspect: %s", exclude_normalization_suspect),
        
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

  p <- apply_optional_log_scale(p, agg$n, axis = "x", requested = log_x)
  p <- apply_optional_log_scale(p, agg$y, axis = "y", requested = log_y)
  p
}
