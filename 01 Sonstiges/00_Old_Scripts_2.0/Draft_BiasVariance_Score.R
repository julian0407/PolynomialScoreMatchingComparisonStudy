# ============================================================
# BiasVariance_Score.R
# Standardisierte Bias-Varianz-Analyse für Score-Schätzer
# ============================================================

# Erwartet:
# source("helper_functions.R")
# source("Draft_Evaluation_Metrics.R")
# plus die Estimator-Dateien wie bisher

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
    grid <- seq(qs[1], qs[2], length.out = grid_size_1d)
    return(grid)
  }
  
  x_probe <- as.matrix(x_probe)
  keep <- apply(x_probe, 1, function(row) all(is.finite(row)))
  x_probe <- x_probe[keep, , drop = FALSE]
  x_probe
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
  
  if (family == "univariate") {
    eval_grid_use <- as.numeric(eval_grid)
    true_mat <- as_score_matrix(true_score(eval_grid_use))
  } else {
    eval_grid_use <- as.matrix(eval_grid)
    true_mat <- as_score_matrix(true_score(eval_grid_use))
  }
  
  B <- n_rep
  n_eval <- nrow(true_mat)
  d <- ncol(true_mat)
  
  score_array <- array(NA_real_, dim = c(B, n_eval, d))
  diag_rows <- vector("list", B)
  
  for (b in seq_len(B)) {
    if (isTRUE(verbose)) {
      message(sprintf("    repetition %d/%d", b, B))
    }
    
    x_train <- r_sample(n)
    if (family == "univariate") x_train <- as.numeric(x_train) else x_train <- as.matrix(x_train)
    
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
    
    base_diag <- data.frame(
      n = n,
      repetition = b,
      method_label = estimator_spec$label,
      method = estimator_spec$method,
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
    
    if (!is.null(score_hat)) {
      score_hat <- tryCatch(as_score_matrix(score_hat, n_expected = n_eval), error = function(e) NULL)
    }
    
    if (!is.null(score_hat)) {
      score_array[b, , ] <- score_hat
    }
    
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
  
  trapz_vec_local <- function(x, y) {
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
      trapz_vec_local(as.numeric(obj$eval_grid), y)
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
    point_bias2 <- bias2_mean_point
  } else if (bias_center == "trimmed_mean") {
    summary_df$integrated_bias2 <- summary_df$integrated_bias2_trimmed_mean
    point_bias2 <- bias2_trimmed_mean_point
  } else if (bias_center == "median") {
    summary_df$integrated_bias2 <- summary_df$integrated_bias2_median
    point_bias2 <- bias2_median_point
  } else {
    summary_df$integrated_bias2 <- summary_df$integrated_bias2_mean
    point_bias2 <- bias2_mean_point
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
  counter_sum <- 1L
  
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
      summaries[[counter_sum]] <- summ$summary
      
      dd <- obj$diagnostics
      dd$bv_target <- "score"
      diagnostics[[counter_sum]] <- dd
      
      counter <- counter + 1L
      counter_sum <- counter_sum + 1L
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

plot_bias_variance_comparison <- function(obj,
                                          metric = c("integrated_bias2",
                                                     "integrated_variance",
                                                     "integrated_bias2_trimmed_mean",
                                                     "integrated_bias2_median",
                                                     "integrated_mse"),
                                          log_x = TRUE,
                                          log_y = TRUE,
                                          main = NULL,
                                          xlab = "sample size n",
                                          ylab = NULL) {
  if (!inherits(obj, "score_bias_variance_benchmark")) {
    stop("obj must be 'score_bias_variance_benchmark'.")
  }
  
  metric <- match.arg(metric)
  df <- obj$summary
  keep <- is.finite(df[[metric]]) & is.finite(df$n)
  df <- df[keep, , drop = FALSE]
  
  if (nrow(df) == 0L) {
    plot.new()
    title(main = metric)
    text(0.5, 0.5, labels = sprintf("No finite values for '%s'", metric))
    return(invisible(NULL))
  }
  
  methods <- unique(df$method_label)
  
  if (is.null(main)) main <- metric
  if (is.null(ylab)) ylab <- metric
  
  if (isTRUE(log_y)) {
    df <- df[df[[metric]] > 0, , drop = FALSE]
  }
  
  if (nrow(df) == 0L) {
    plot.new()
    title(main = metric)
    text(0.5, 0.5, labels = sprintf("No positive values for '%s'", metric))
    return(invisible(NULL))
  }
  
  plot(
    NA, NA,
    xlim = range(df$n, na.rm = TRUE),
    ylim = range(df[[metric]], na.rm = TRUE),
    log = paste0(if (log_x) "x" else "", if (log_y) "y" else ""),
    xlab = xlab,
    ylab = ylab,
    main = main
  )
  
  drawn <- character(0)
  
  for (i in seq_along(methods)) {
    dd <- df[df$method_label == methods[i], , drop = FALSE]
    dd <- dd[order(dd$n), , drop = FALSE]
    if (nrow(dd) == 0L) next
    lines(dd$n, dd[[metric]], type = "b", lwd = 2, pch = 19 + i - 1, col = i)
    drawn <- c(drawn, methods[i])
  }
  
  if (length(drawn) > 0L) {
    idx <- match(drawn, methods)
    legend(
      "topright",
      legend = drawn,
      col = idx,
      lty = 1,
      pch = 19 + idx - 1,
      bty = "n"
    )
  }
  
  invisible(df)
}

plot_bias_variance_panel <- function(obj,
                                     metrics = c("integrated_bias2",
                                                 "integrated_variance",
                                                 "integrated_mse"),
                                     log_x = TRUE,
                                     log_y = TRUE) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  k <- length(metrics)
  par(mfrow = c(1, k))
  
  for (met in metrics) {
    plot_bias_variance_comparison(
      obj = obj,
      metric = met,
      log_x = log_x,
      log_y = log_y,
      main = met,
      ylab = met
    )
  }
  
  invisible(NULL)
}