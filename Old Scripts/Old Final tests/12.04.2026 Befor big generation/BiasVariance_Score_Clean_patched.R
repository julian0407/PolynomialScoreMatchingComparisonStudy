# ============================================================
# BiasVariance_Score_Clean.R
# Standardisierte Bias-Varianz-Analyse für Score-Schätzer
# Schlank gehalten für das aktuelle Studiendesign
# ============================================================

# Erwartet:
# source("helper_functions.R")
# source("Evaluation_Metrics_Clean.R")
# plus die Estimator-Dateien wie bisher

`%||%` <- function(x, y) if (is.null(x)) y else x

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_var <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::var(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::median(x)
}

trapz_vec <- function(x, y) {
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

make_eval_grid <- function(r_sample,
                           family = c("univariate", "multivariate"),
                           n_probe = 5000,
                           grid_size_1d = 401,
                           probs_1d = c(0.005, 0.995),
                           seed = NULL) {
  family <- match.arg(family)
  if (!is.null(seed)) set.seed(seed)

  x_probe <- r_sample(n_probe)

  if (family == "univariate") {
    x_probe <- as.numeric(x_probe)
    x_probe <- x_probe[is.finite(x_probe)]
    qs <- stats::quantile(x_probe, probs = probs_1d, na.rm = TRUE, names = FALSE)
    return(seq(qs[1], qs[2], length.out = grid_size_1d))
  }

  x_probe <- as.matrix(x_probe)
  keep <- apply(x_probe, 1, function(row) all(is.finite(row)))
  x_probe[keep, , drop = FALSE]
}

run_bias_variance_score_one <- function(n,
                                        family = c("univariate", "multivariate"),
                                        estimator_spec,
                                        r_sample,
                                        true_score,
                                        eval_grid,
                                        n_rep = 50,
                                        seed = NULL,
                                        verbose = TRUE) {
  family <- match.arg(family)
  if (!is.null(seed)) set.seed(seed)

  eval_grid_use <- if (family == "univariate") as.numeric(eval_grid) else as.matrix(eval_grid)
  true_mat <- as_score_matrix(true_score(eval_grid_use))

  B <- n_rep
  n_eval <- nrow(true_mat)
  d <- ncol(true_mat)

  score_array <- array(NA_real_, dim = c(B, n_eval, d))
  diag_rows <- vector("list", B)

  for (b in seq_len(B)) {
    if (isTRUE(verbose)) message(sprintf("    repetition %d/%d", b, B))

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

    base_diag <- data.frame(
      n = n,
      repetition = b,
      method_label = estimator_spec$label,
      method = estimator_spec$method,
      fit_time_sec = as.numeric(fit_time["elapsed"]),
      score_inference_time_sec = NA_real_,
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

    if (inherits(fit, "fit_error")) {
      base_diag$status <- fit$error_message
      diag_rows[[b]] <- base_diag
      next
    }

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

    score_time <- system.time({
      score_hat <- tryCatch(
        do.call(
          predict_score_estimator_generic,
          c(
            list(
              newx = eval_grid_use,
              fit = fit,
              family = family,
              method = estimator_spec$method
            ),
            estimator_spec$score_predict_args %||% list()
          )
        ),
        error = function(e) NULL
      )
    })

    base_diag$score_inference_time_sec <- as.numeric(score_time["elapsed"])

    if (!is.null(score_hat)) {
      score_hat <- tryCatch(as_score_matrix(score_hat, n_expected = n_eval), error = function(e) NULL)
    }
    if (!is.null(score_hat)) score_array[b, , ] <- score_hat

    base_diag$success <- diags$success %||% TRUE
    base_diag$status <- diags$status %||% NA_character_
    base_diag$iterations <- diags$iterations %||% NA_real_
    base_diag$objective_value <- diags$objective_value %||% NA_real_
    base_diag$kappa_raw <- diags$kappa_raw %||% NA_real_
    base_diag$kappa_reg <- diags$kappa_reg %||% NA_real_
    base_diag$rcond_raw <- diags$rcond_raw %||% NA_real_
    base_diag$rcond_reg <- diags$rcond_reg %||% NA_real_
    base_diag$eigmin_raw <- diags$eigmin_raw %||% NA_real_
    base_diag$eigmin_reg <- diags$eigmin_reg %||% NA_real_

    diag_rows[[b]] <- base_diag
  }

  diag_df <- do.call(rbind, diag_rows)
  rownames(diag_df) <- NULL

  structure(
    list(
      n = n,
      family = family,
      method_label = estimator_spec$label,
      method = estimator_spec$method,
      eval_grid = eval_grid_use,
      true_score = true_mat,
      score_array = score_array,
      diagnostics = diag_df
    ),
    class = "score_bias_variance_run"
  )
}

summarize_bias_variance_score_run <- function(obj,
                                              integrate = TRUE,
                                              bias_center = c("mean", "trimmed_mean", "median", "all"),
                                              trim_alpha = 0.05) {
  if (!inherits(obj, "score_bias_variance_run")) {
    stop("obj must be of class 'score_bias_variance_run'.")
  }
  
  bias_center <- match.arg(bias_center)
  
  if (!is.numeric(trim_alpha) || length(trim_alpha) != 1L ||
      !is.finite(trim_alpha) || trim_alpha < 0 || trim_alpha >= 0.5) {
    stop("trim_alpha must be a scalar in [0, 0.5).")
  }
  
  arr <- obj$score_array
  true_mat <- obj$true_score
  
  B <- dim(arr)[1]
  n_eval <- dim(arr)[2]
  d <- dim(arr)[3]
  
  safe_mean_local <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    mean(x)
  }
  
  safe_var_local <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(NA_real_)
    stats::var(x)
  }
  
  safe_median_local <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    stats::median(x)
  }
  
  safe_trimmed_mean_local <- function(x, trim = 0.05) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    mean(x, trim = trim)
  }
  
  mean_hat <- matrix(NA_real_, nrow = n_eval, ncol = d)
  trimmed_mean_hat <- matrix(NA_real_, nrow = n_eval, ncol = d)
  median_hat <- matrix(NA_real_, nrow = n_eval, ncol = d)
  
  var_hat <- matrix(NA_real_, nrow = n_eval, ncol = d)
  mse <- matrix(NA_real_, nrow = n_eval, ncol = d)
  
  bias2_mean <- matrix(NA_real_, nrow = n_eval, ncol = d)
  bias2_trimmed_mean <- matrix(NA_real_, nrow = n_eval, ncol = d)
  bias2_median <- matrix(NA_real_, nrow = n_eval, ncol = d)
  
  n_finite_mat <- matrix(0L, nrow = n_eval, ncol = d)
  
  for (j in seq_len(n_eval)) {
    for (k in seq_len(d)) {
      vals <- arr[, j, k]
      vals <- vals[is.finite(vals)]
      
      n_finite_mat[j, k] <- length(vals)
      
      if (length(vals) == 0L) next
      
      mean_hat[j, k] <- safe_mean_local(vals)
      trimmed_mean_hat[j, k] <- safe_trimmed_mean_local(vals, trim = trim_alpha)
      median_hat[j, k] <- safe_median_local(vals)
      
      var_hat[j, k] <- safe_var_local(vals)
      mse[j, k] <- mean((vals - true_mat[j, k])^2)
      
      bias2_mean[j, k] <- (mean_hat[j, k] - true_mat[j, k])^2
      bias2_trimmed_mean[j, k] <- (trimmed_mean_hat[j, k] - true_mat[j, k])^2
      bias2_median[j, k] <- (median_hat[j, k] - true_mat[j, k])^2
    }
  }
  
  bias2_mean_point <- rowSums(bias2_mean, na.rm = TRUE)
  bias2_trimmed_mean_point <- rowSums(bias2_trimmed_mean, na.rm = TRUE)
  bias2_median_point <- rowSums(bias2_median, na.rm = TRUE)
  var_point <- rowSums(var_hat, na.rm = TRUE)
  mse_point <- rowSums(mse, na.rm = TRUE)
  
  integrate_vec <- function(y) {
    if (obj$family == "univariate" && isTRUE(integrate)) {
      trapz_vec(as.numeric(obj$eval_grid), y)
    } else {
      safe_mean_local(y)
    }
  }
  
  integrated_bias2_mean <- integrate_vec(bias2_mean_point)
  integrated_bias2_trimmed_mean <- integrate_vec(bias2_trimmed_mean_point)
  integrated_bias2_median <- integrate_vec(bias2_median_point)
  integrated_variance <- integrate_vec(var_point)
  integrated_mse <- integrate_vec(mse_point)
  
  summary_df <- data.frame(
    n = obj$n,
    method_label = obj$method_label,
    method = obj$method,
    integrated_bias2_mean = integrated_bias2_mean,
    integrated_bias2_trimmed_mean = integrated_bias2_trimmed_mean,
    integrated_bias2_median = integrated_bias2_median,
    integrated_variance = integrated_variance,
    integrated_mse = integrated_mse,
    finite_rep_mean = safe_mean_local(as.numeric(n_finite_mat) / B),
    trim_alpha = trim_alpha,
    stringsAsFactors = FALSE
  )
  
  if (bias_center == "mean") {
    summary_df$integrated_bias2 <- summary_df$integrated_bias2_mean
  } else if (bias_center == "trimmed_mean") {
    summary_df$integrated_bias2 <- summary_df$integrated_bias2_trimmed_mean
  } else if (bias_center == "median") {
    summary_df$integrated_bias2 <- summary_df$integrated_bias2_median
  } else {
    summary_df$integrated_bias2 <- summary_df$integrated_bias2_mean
  }
  
  list(
    pointwise = data.frame(
      eval_id = seq_len(n_eval),
      bias2_mean = bias2_mean_point,
      bias2_trimmed_mean = bias2_trimmed_mean_point,
      bias2_median = bias2_median_point,
      variance = var_point,
      mse = mse_point,
      stringsAsFactors = FALSE
    ),
    summary = summary_df,
    mean_hat = mean_hat,
    trimmed_mean_hat = trimmed_mean_hat,
    median_hat = median_hat,
    var_hat = var_hat,
    mse = mse,
    bias2_mean = bias2_mean,
    bias2_trimmed_mean = bias2_trimmed_mean,
    bias2_median = bias2_median,
    n_finite = n_finite_mat
  )
}

run_bias_variance_score_benchmark <- function(sample_sizes,
                                              family = c("univariate", "multivariate"),
                                              estimator_specs,
                                              r_sample,
                                              true_score,
                                              eval_grid = NULL,
                                              n_rep = 50,
                                              seed = NULL,
                                              verbose = TRUE,
                                              n_probe_grid = 5000,
                                              grid_size_1d = 401) {
  family <- match.arg(family)
  if (!is.null(seed)) set.seed(seed)

  if (is.null(eval_grid)) {
    eval_grid <- make_eval_grid(
      r_sample = r_sample,
      family = family,
      n_probe = n_probe_grid,
      grid_size_1d = grid_size_1d,
      seed = seed
    )
  }

  runs <- list()
  summaries <- list()
  diagnostics <- list()
  counter <- 1L

  for (spec in estimator_specs) {
    if (isTRUE(verbose)) message("Method: ", spec$label)

    for (n in sample_sizes) {
      if (isTRUE(verbose)) message("  n = ", n)

      obj <- run_bias_variance_score_one(
        n = n,
        family = family,
        estimator_spec = spec,
        r_sample = r_sample,
        true_score = true_score,
        eval_grid = eval_grid,
        n_rep = n_rep,
        verbose = verbose
      )

      summ <- summarize_bias_variance_score_run(
        obj,
        bias_center = "all",
        trim_alpha = 0.05
      )

      runs[[counter]] <- obj
      summaries[[counter]] <- summ$summary
      dd <- obj$diagnostics
      dd$bv_target <- "score"
      diagnostics[[counter]] <- dd
      counter <- counter + 1L
    }
  }

  summary_df <- do.call(rbind, summaries)
  rownames(summary_df) <- NULL

  diagnostics_df <- do.call(rbind, diagnostics)
  rownames(diagnostics_df) <- NULL

  structure(
    list(
      summary = summary_df,
      diagnostics = diagnostics_df,
      runs = runs,
      settings = list(
        sample_sizes = sample_sizes,
        family = family,
        n_rep = n_rep
      )
    ),
    class = "score_bias_variance_benchmark"
  )
}


subset_score_bias_variance_benchmark <- function(obj,
                                                 keep_n = NULL,
                                                 drop_n = NULL,
                                                 keep_method_labels = NULL,
                                                 drop_method_labels = NULL,
                                                 keep_methods = NULL,
                                                 drop_methods = NULL,
                                                 update_settings = TRUE) {
  if (!inherits(obj, "score_bias_variance_benchmark")) {
    stop("obj must be of class 'score_bias_variance_benchmark'.")
  }

  filt_n <- normalize_n_filter(keep_n = keep_n, drop_n = drop_n)
  filt_m <- normalize_method_filter(
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )
  keep_fun <- function(run) {
    out <- TRUE
    if (!is.null(filt_n$keep_n)) out <- out && (run$n %in% filt_n$keep_n)
    if (!is.null(filt_n$drop_n)) out <- out && !(run$n %in% filt_n$drop_n)
    if (!is.null(filt_m$keep_method_labels)) out <- out && (run$method_label %in% filt_m$keep_method_labels)
    if (!is.null(filt_m$drop_method_labels)) out <- out && !(run$method_label %in% filt_m$drop_method_labels)
    if (!is.null(filt_m$keep_methods)) out <- out && (run$method %in% filt_m$keep_methods)
    if (!is.null(filt_m$drop_methods)) out <- out && !(run$method %in% filt_m$drop_methods)
    out
  }

  out <- obj
  out$summary <- filter_benchmark_df_by_n(out$summary, keep_n = keep_n, drop_n = drop_n)
  out$summary <- filter_benchmark_df_by_method(
    out$summary,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )
  out$diagnostics <- filter_benchmark_df_by_n(out$diagnostics, keep_n = keep_n, drop_n = drop_n)
  out$diagnostics <- filter_benchmark_df_by_method(
    out$diagnostics,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )
  out$runs <- Filter(keep_fun, out$runs)

  if (isTRUE(update_settings)) {
    out$settings$sample_sizes <- sort(unique(out$summary$n))
  }

  out
}

aggregate_score_bias_variance_benchmark <- function(obj,
                                                    value_cols = c("integrated_bias2", "integrated_variance", "integrated_mse"),
                                                    keep_n = NULL,
                                                    drop_n = NULL,
                                                    keep_method_labels = NULL,
                                                    drop_method_labels = NULL,
                                                    keep_methods = NULL,
                                                    drop_methods = NULL) {
  if (!inherits(obj, "score_bias_variance_benchmark")) {
    stop("obj must be of class 'score_bias_variance_benchmark'.")
  }

  df <- filter_benchmark_df_by_n(obj$summary, keep_n = keep_n, drop_n = drop_n)
  df <- filter_benchmark_df_by_method(
    df,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )
  miss <- setdiff(value_cols, names(df))
  if (length(miss) > 0L) {
    stop(sprintf("Requested columns not found: %s", paste(miss, collapse = ", ")))
  }
  df[, c("method_label", "method", "n", value_cols), drop = FALSE]
}

plot_score_bias_variance <- function(obj,
                                     metric = c("integrated_bias2", "integrated_variance", "integrated_mse"),
                                     keep_n = NULL,
                                     drop_n = NULL,
                                     keep_method_labels = NULL,
                                     drop_method_labels = NULL,
                                     keep_methods = NULL,
                                     drop_methods = NULL,
                                     log_x = FALSE,
                                     log_y = FALSE) {
  metric <- match.arg(metric)
  df <- aggregate_score_bias_variance_benchmark(
    obj = obj,
    value_cols = metric,
    keep_n = keep_n,
    drop_n = drop_n,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )

  if (nrow(df) == 0L) {
    stop("No rows left after applying the sample-size filter.")
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = n, y = .data[[metric]], color = method_label, group = method_label)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Sample size n",
      y = metric,
      color = "Method",
      title = sprintf("Score bias-variance: %s", metric),
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

  if (isTRUE(log_x)) p <- p + ggplot2::scale_x_log10()
  if (isTRUE(log_y)) p <- p + ggplot2::scale_y_log10()

  p
}
