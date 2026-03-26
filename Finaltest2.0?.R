# ============================================================
# Finale_Tests_Robust.R
# Robustes finales Testskript
# - mehrere Methoden vergleichbar
# - mehrere Metriken auswählbar
# - optional 2er-Panelplot
# - toleriert fehlende diagnostics
# ============================================================

source("helper_functions.R")
source("KDE.R")
source("LogConcaveMLE.R")
source("Univariate_Polynomial_Score_Matching_1.0.R")
source("Multivariate_Pairwise_Polynomial_Score_Matching.R")

# ============================================================
# (1) Generische Wrapper
# ============================================================

fit_estimator <- function(x,
                          family = c("univariate", "multivariate"),
                          method = c("MLE", "KDE", "SM"),
                          smoothed = FALSE,
                          ...) {
  family <- match.arg(family)
  method <- match.arg(method)
  
  if (family == "univariate") {
    if (method == "KDE") return(fit_kde_1d(x, ...))
    if (method == "MLE") return(fit_logconcave_mle_1d(x, smoothed = smoothed, ...))
    if (method == "SM")  return(fit_score_matching_univariate(x, ...))
  }
  
  if (family == "multivariate") {
    if (method == "KDE") return(fit_kde_mv(x, ...))
    if (method == "MLE") return(fit_logconcave_mle_mv(x, smoothed = smoothed, ...))
    if (method == "SM")  return(fit_score_matching_mv_basic(x, ...))
  }
  
  stop("Unsupported family / method combination.")
}

predict_density_estimator <- function(newx,
                                      fit,
                                      family = c("univariate", "multivariate"),
                                      method = c("MLE", "KDE", "SM"),
                                      ...) {
  family <- match.arg(family)
  method <- match.arg(method)
  
  if (family == "univariate") {
    if (method == "KDE") return(as.numeric(predict_density_kde_1d(newx, fit, ...)))
    if (method == "MLE") return(as.numeric(predict_density_logconcave_1d(newx, fit, ...)))
    if (method == "SM")  return(as.numeric(predict_density_sm_1d(newx, fit, ...)))
  }
  
  if (family == "multivariate") {
    if (method == "KDE") return(as.numeric(predict_density_kde_mv(newx, fit, ...)))
    if (method == "MLE") return(as.numeric(predict_density_logconcave_mv(newx, fit, ...)))
    if (method == "SM")  stop("Density metrics are currently disabled for multivariate SM.")
  }
  
  stop("Unsupported family / method combination.")
}

predict_logdensity_estimator <- function(newx,
                                         fit,
                                         family = c("univariate", "multivariate"),
                                         method = c("MLE", "KDE", "SM"),
                                         ...) {
  family <- match.arg(family)
  method <- match.arg(method)
  
  if (family == "univariate") {
    if (method == "KDE") return(as.numeric(predict_logdensity_kde_1d(newx, fit, ...)))
    if (method == "MLE") return(as.numeric(predict_logdensity_logconcave_1d(newx, fit, ...)))
    if (method == "SM")  return(as.numeric(predict_logdensity_sm_1d(newx, fit, ...)))
  }
  
  if (family == "multivariate") {
    if (method == "KDE") return(as.numeric(predict_logdensity_kde_mv(newx, fit, ...)))
    if (method == "MLE") return(as.numeric(predict_logdensity_logconcave_mv(newx, fit, ...)))
    if (method == "SM")  stop("Log-density metrics are currently disabled for multivariate SM.")
  }
  
  stop("Unsupported family / method combination.")
}

predict_score_estimator <- function(newx,
                                    fit,
                                    family = c("univariate", "multivariate"),
                                    method = c("MLE", "KDE", "SM"),
                                    ...) {
  family <- match.arg(family)
  method <- match.arg(method)
  
  if (family == "univariate") {
    if (method == "KDE") return(as.numeric(predict_score_kde_1d(newx, fit, ...)))
    if (method == "MLE") return(as.numeric(predict_score_logconcave_1d(newx, fit, ...)))
    if (method == "SM")  return(as.numeric(predict_score_univariate(newx, fit, ...)))
  }
  
  if (family == "multivariate") {
    if (method == "KDE") return(as.matrix(predict_score_kde_mv(newx, fit, ...)))
    if (method == "MLE") return(as.matrix(predict_score_logconcave_mv(newx, fit, ...)))
    if (method == "SM")  return(as.matrix(predict_score_mv_basic(newx, fit, ...)))
  }
  
  stop("Unsupported family / method combination.")
}

# ============================================================
# (2) Hilfsfunktionen
# ============================================================

as_score_matrix <- function(s, n_expected = NULL) {
  if (is.null(dim(s))) {
    s <- matrix(as.numeric(s), ncol = 1)
  } else {
    s <- as.matrix(s)
  }
  
  if (!is.null(n_expected) && nrow(s) != n_expected) {
    stop("Score output has unexpected number of rows.")
  }
  
  s
}

clean_complete_cases_pair <- function(a, b) {
  a <- as_score_matrix(a)
  b <- as_score_matrix(b)
  
  if (nrow(a) != nrow(b) || ncol(a) != ncol(b)) {
    stop("Score matrices have incompatible dimensions.")
  }
  
  keep <- apply(a, 1, function(row) all(is.finite(row))) &
    apply(b, 1, function(row) all(is.finite(row)))
  
  list(
    a = a[keep, , drop = FALSE],
    b = b[keep, , drop = FALSE],
    keep = keep
  )
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_rmse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  sqrt(mean(x^2))
}

trapz_1d <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2L) return(NA_real_)
  idx <- order(x)
  x <- x[idx]
  y <- y[idx]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

num_second_derivative_1d <- function(fun, x, h = 1e-4) {
  (fun(x + h) - 2 * fun(x) + fun(x - h)) / (h^2)
}

# ============================================================
# (3) Robuste Diagnostics-Auslese
# ============================================================

extract_fit_diagnostics_safe <- function(fit) {
  out <- list(
    success = NA,
    status = NA_character_,
    iterations = NA_real_,
    objective_value = NA_real_
  )
  
  if (is.null(fit) || inherits(fit, "fit_error")) {
    out$success <- FALSE
    out$status <- "fit_error"
    return(out)
  }
  
  # Falls fit$diagnostics schon existiert
  if (!is.null(fit$diagnostics) && is.list(fit$diagnostics)) {
    d <- fit$diagnostics
    if (!is.null(d$success))         out$success <- d$success
    if (!is.null(d$status))          out$status <- as.character(d$status)
    if (!is.null(d$iterations))      out$iterations <- suppressWarnings(as.numeric(d$iterations))
    if (!is.null(d$objective_value)) out$objective_value <- suppressWarnings(as.numeric(d$objective_value))
    return(out)
  }
  
  # Univariates SM: status + solution
  if (!is.null(fit$status)) {
    out$status <- as.character(fit$status)
    out$success <- grepl("optimal", tolower(out$status))
  }
  
  if (!is.null(fit$solution)) {
    val <- tryCatch(as.numeric(fit$solution$value), error = function(e) NA_real_)
    if (is.finite(val)) out$objective_value <- val
    
    iter_candidates <- c(
      tryCatch(as.numeric(fit$solution$num_iters), error = function(e) NA_real_),
      tryCatch(as.numeric(fit$solution$solver_stats$num_iters), error = function(e) NA_real_),
      tryCatch(as.numeric(fit$solution$solver_stats$iter), error = function(e) NA_real_)
    )
    iter_candidates <- iter_candidates[is.finite(iter_candidates)]
    if (length(iter_candidates) > 0L) out$iterations <- iter_candidates[1L]
  }
  
  # Multivariates SM: optimizer
  if (!is.null(fit$optimizer) && is.list(fit$optimizer)) {
    opt <- fit$optimizer
    if (!is.null(opt$convergence)) out$success <- identical(opt$convergence, 0)
    if (!is.null(opt$message)) out$status <- as.character(opt$message)
    if (!is.null(opt$value)) {
      val <- suppressWarnings(as.numeric(opt$value))
      if (is.finite(val)) out$objective_value <- val
    }
    if (!is.null(opt$counts)) {
      it <- suppressWarnings(sum(as.numeric(unlist(opt$counts)), na.rm = TRUE))
      if (is.finite(it)) out$iterations <- it
    }
  }
  
  out
}

# ============================================================
# (4) Metriken
# ============================================================

metric_negloglik <- function(x_test, fit, family, method, predict_args = list()) {
  log_hat <- do.call(
    predict_logdensity_estimator,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  -safe_mean(log_hat)
}

metric_kl <- function(x_test, fit, family, method, true_logdensity, predict_args = list()) {
  log_hat <- do.call(
    predict_logdensity_estimator,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  log_true <- as.numeric(true_logdensity(x_test))
  keep <- is.finite(log_hat) & is.finite(log_true)
  safe_mean(log_true[keep] - log_hat[keep])
}

metric_hellinger2 <- function(x_test, fit, family, method, true_density, predict_args = list()) {
  f_hat <- do.call(
    predict_density_estimator,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  f_true <- as.numeric(true_density(x_test))
  keep <- is.finite(f_hat) & is.finite(f_true) & (f_hat >= 0) & (f_true >= 0)
  safe_mean((sqrt(f_hat[keep]) - sqrt(f_true[keep]))^2)
}

metric_ise_1d <- function(grid_1d, fit, family, method, true_density, predict_args = list()) {
  if (family != "univariate") return(NA_real_)
  f_hat <- do.call(
    predict_density_estimator,
    c(list(newx = grid_1d, fit = fit, family = family, method = method), predict_args)
  )
  f_true <- as.numeric(true_density(grid_1d))
  trapz_1d(grid_1d, (f_hat - f_true)^2)
}

metric_score_loss <- function(x_test, fit, family, method, true_score, predict_args = list()) {
  score_hat <- do.call(
    predict_score_estimator,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  score_true <- true_score(x_test)
  
  if (family == "univariate") {
    score_hat <- as_score_matrix(score_hat, n_expected = length(x_test))
    score_true <- as_score_matrix(score_true, n_expected = length(x_test))
  } else {
    score_hat <- as_score_matrix(score_hat, n_expected = nrow(x_test))
    score_true <- as_score_matrix(score_true, n_expected = nrow(x_test))
  }
  
  tmp <- clean_complete_cases_pair(score_hat, score_true)
  err <- tmp$a - tmp$b
  safe_mean(rowSums(err^2))
}

metric_fisher <- function(x_test, fit, family, method, true_score, predict_args = list()) {
  metric_score_loss(x_test, fit, family, method, true_score, predict_args)
}

metric_score_rmse <- function(x_test, fit, family, method, true_score, predict_args = list()) {
  score_hat <- do.call(
    predict_score_estimator,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  score_true <- true_score(x_test)
  
  if (family == "univariate") {
    score_hat <- as_score_matrix(score_hat, n_expected = length(x_test))
    score_true <- as_score_matrix(score_true, n_expected = length(x_test))
  } else {
    score_hat <- as_score_matrix(score_hat, n_expected = nrow(x_test))
    score_true <- as_score_matrix(score_true, n_expected = nrow(x_test))
  }
  
  tmp <- clean_complete_cases_pair(score_hat, score_true)
  err <- tmp$a - tmp$b
  sqrt(safe_mean(rowSums(err^2)))
}

metric_logconcavity <- function(fit, family, method, grid_1d = NULL, predict_args = list(), tol = 1e-8) {
  out <- list(
    min_hessian_eigenvalue = NA_real_,
    share_violated = NA_real_,
    mean_violation = NA_real_,
    max_violation = NA_real_
  )
  
  # Multivariates SM: vorhandene lc_diagnostics benutzen
  if (family == "multivariate" && method == "SM" && !is.null(fit$lc_diagnostics)) {
    diag <- fit$lc_diagnostics
    if (!is.null(diag$min_eigenvalues)) {
      out$min_hessian_eigenvalue <- min(diag$min_eigenvalues, na.rm = TRUE)
      out$share_violated <- if (!is.null(diag$n_violated) && length(diag$min_eigenvalues) > 0L) {
        diag$n_violated / length(diag$min_eigenvalues)
      } else {
        NA_real_
      }
    }
    if (!is.null(diag$mean_violation)) out$mean_violation <- diag$mean_violation
    if (!is.null(diag$max_violation)) out$max_violation <- diag$max_violation
    return(out)
  }
  
  # Univariat: numerisch über log-density auf Grid
  if (family == "univariate" && !is.null(grid_1d)) {
    f_log <- function(xx) {
      do.call(
        predict_logdensity_estimator,
        c(list(newx = xx, fit = fit, family = family, method = method), predict_args)
      )
    }
    
    d2 <- vapply(grid_1d, function(xx) num_second_derivative_1d(f_log, xx, h = 1e-4), numeric(1))
    violations <- pmax(d2 - tol, 0)
    
    out$min_hessian_eigenvalue <- suppressWarnings(min(-d2, na.rm = TRUE))
    out$share_violated <- suppressWarnings(mean(violations > 0, na.rm = TRUE))
    out$mean_violation <- suppressWarnings(mean(violations, na.rm = TRUE))
    out$max_violation <- suppressWarnings(max(violations, na.rm = TRUE))
  }
  
  out
}

# ============================================================
# (5) Spalten-Handling, damit rbind immer klappt
# ============================================================

get_metric_columns <- function(metrics) {
  out <- character(0)
  
  for (met in metrics) {
    if (met %in% c("negloglik", "kl", "hellinger2", "ise",
                   "score_loss", "fisher", "score_rmse")) {
      out <- c(out, met)
    } else if (met == "logconcavity") {
      out <- c(out,
               "min_hessian_eigenvalue",
               "share_violated",
               "mean_violation",
               "max_violation")
    } else {
      stop(sprintf("Unknown metric: %s", met))
    }
  }
  
  unique(out)
}

align_result_row <- function(row, target_names) {
  miss <- setdiff(target_names, names(row))
  if (length(miss) > 0L) {
    for (nm in miss) row[[nm]] <- NA
  }
  
  row <- row[, target_names, drop = FALSE]
  row
}

# ============================================================
# (6) Eine Wiederholung
# ============================================================

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
  
  metric_cols <- get_metric_columns(metrics)
  target_cols <- c(
    "n",
    "method_label",
    "method",
    "runtime_sec",
    "success",
    "status",
    "iterations",
    "objective_value",
    metric_cols
  )
  
  row <- data.frame(
    n = n,
    method_label = estimator_spec$label,
    method = estimator_spec$method,
    runtime_sec = NA_real_,
    success = NA,
    status = NA_character_,
    iterations = NA_real_,
    objective_value = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (nm in metric_cols) row[[nm]] <- NA_real_
  
  x_train <- r_sample(n)
  if (family == "univariate") x_train <- as.numeric(x_train) else x_train <- as.matrix(x_train)
  
  timing <- system.time({
    fit <- tryCatch(
      do.call(
        fit_estimator,
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
  
  row$runtime_sec <- as.numeric(timing["elapsed"])
  
  if (inherits(fit, "fit_error")) {
    row$success <- FALSE
    row$status <- fit$error_message
    return(row[, target_cols, drop = FALSE])
  }
  
  diags <- extract_fit_diagnostics_safe(fit)
  row$success <- diags$success
  row$status <- diags$status
  row$iterations <- diags$iterations
  row$objective_value <- diags$objective_value
  
  x_test <- r_sample(n_test)
  if (family == "univariate") x_test <- as.numeric(x_test) else x_test <- as.matrix(x_test)
  
  for (met in metrics) {
    val <- tryCatch({
      if (met == "negloglik") {
        metric_negloglik(x_test, fit, family, estimator_spec$method, estimator_spec$predict_args)
      } else if (met == "kl") {
        metric_kl(x_test, fit, family, estimator_spec$method, true_logdensity, estimator_spec$predict_args)
      } else if (met == "hellinger2") {
        metric_hellinger2(x_test, fit, family, estimator_spec$method, true_density, estimator_spec$predict_args)
      } else if (met == "ise") {
        metric_ise_1d(grid_1d, fit, family, estimator_spec$method, true_density, estimator_spec$predict_args)
      } else if (met == "score_loss") {
        metric_score_loss(x_test, fit, family, estimator_spec$method, true_score, estimator_spec$predict_args)
      } else if (met == "fisher") {
        metric_fisher(x_test, fit, family, estimator_spec$method, true_score, estimator_spec$predict_args)
      } else if (met == "score_rmse") {
        metric_score_rmse(x_test, fit, family, estimator_spec$method, true_score, estimator_spec$predict_args)
      } else if (met == "logconcavity") {
        metric_logconcavity(fit, family, estimator_spec$method, grid_1d, estimator_spec$predict_args)
      } else {
        NA_real_
      }
    }, error = function(e) {
      if (met == "logconcavity") {
        list(
          min_hessian_eigenvalue = NA_real_,
          share_violated = NA_real_,
          mean_violation = NA_real_,
          max_violation = NA_real_
        )
      } else {
        NA_real_
      }
    })
    
    if (met == "logconcavity") {
      row$min_hessian_eigenvalue <- val$min_hessian_eigenvalue
      row$share_violated <- val$share_violated
      row$mean_violation <- val$mean_violation
      row$max_violation <- val$max_violation
    } else {
      row[[met]] <- val
    }
  }
  
  row[, target_cols, drop = FALSE]
}

# ============================================================
# (7) Gesamter Benchmark
# ============================================================

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
  
  sample_sizes <- as.integer(sample_sizes)
  sample_sizes <- sample_sizes[is.finite(sample_sizes) & sample_sizes >= 2L]
  if (length(sample_sizes) == 0L) stop("sample_sizes must contain at least one integer >= 2.")
  
  metric_cols <- get_metric_columns(metrics)
  target_cols <- c(
    "n", "repetition", "method_label", "method",
    "runtime_sec", "success", "status", "iterations", "objective_value",
    metric_cols
  )
  
  out <- vector("list", length(estimator_specs) * length(sample_sizes) * n_rep)
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
        ans <- align_result_row(ans, target_cols)
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

# ============================================================
# (8) Aggregation
# ============================================================

aggregate_final_benchmark <- function(obj, metric) {
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")
  df <- obj$raw
  
  if (!metric %in% names(df)) {
    stop(sprintf("Metric '%s' not found in benchmark output.", metric))
  }
  
  split_key <- interaction(df$method_label, df$n, drop = TRUE)
  
  out <- lapply(split(df, split_key), function(dd) {
    vals <- dd[[metric]]
    vals <- vals[is.finite(vals)]
    
    data.frame(
      method_label = dd$method_label[1],
      n = dd$n[1],
      mean = if (length(vals) > 0L) mean(vals) else NA_real_,
      median = if (length(vals) > 0L) median(vals) else NA_real_,
      sd = if (length(vals) > 1L) stats::sd(vals) else NA_real_,
      q25 = if (length(vals) > 0L) stats::quantile(vals, 0.25, names = FALSE) else NA_real_,
      q75 = if (length(vals) > 0L) stats::quantile(vals, 0.75, names = FALSE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}

# ============================================================
# (9) Plot
# ============================================================

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
  
  ycol <- if (use_median) "median" else "mean"
  
  if (is.null(main)) main <- paste("Metric:", metric)
  if (is.null(ylab)) ylab <- metric
  
  xlim <- range(agg$n, na.rm = TRUE)
  ylim <- range(agg[[ycol]], na.rm = TRUE, finite = TRUE)
  
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
    dd <- agg[agg$method_label == methods[i], , drop = FALSE]
    dd <- dd[order(dd$n), , drop = FALSE]
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

# ============================================================
# (10) Beispiel univariat
# ============================================================

r_sample_norm_1d <- function(n) rnorm(n, mean = 0, sd = 1)
true_density_norm_1d <- function(x) dnorm(x, mean = 0, sd = 1)
true_logdensity_norm_1d <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)
true_score_norm_1d <- function(x) matrix(x, ncol = 1)

grid_1d <- seq(-5, 5, length.out = 2001)

estimators_1d <- list(
  # list(
  #   label = "KDE_SJ",
  #   method = "KDE",
  #   smoothed = FALSE,
  #   fit_args = list(bw = "SJ"),
  #   predict_args = list()
  # ),
  # list(
  #   label = "MLE_unsmoothed",
  #   method = "MLE",
  #   smoothed = TRUE,
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

# Beispiel:
res_1d <- run_final_benchmark(
  sample_sizes = c(50, 100, 1000),
  family = "univariate",
  estimator_specs = estimators_1d,
  r_sample = r_sample_norm_1d,
  metrics = c("score_loss", "score_rmse"),
  n_rep = 5,
  n_test = 500,
  true_density = true_density_norm_1d,
  true_logdensity = true_logdensity_norm_1d,
  true_score = true_score_norm_1d,
  grid_1d = grid_1d,
  seed = 123
)

plot_metric_panel(res_1d, metrics = c("score_loss"))
plot_metric_panel(res_1d, metrics = c("ise", "score_rmse"))
plot_metric_panel(res_1d, metrics = c("logconcavity"))

# ============================================================
# (11) Beispiel multivariat
# ============================================================

r_sample_norm_2d <- function(n) cbind(rnorm(n), rnorm(n))
true_score_norm_2d <- function(x) as.matrix(x)

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
  )
)

# Beispiel:
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






x_train <- rnorm(50)
x_test  <- rnorm(2000)

fit_kde <- fit_estimator(
  x_train, family = "univariate", method = "KDE",
  bw = "SJ"
)

fit_mle <- fit_estimator(
  x_train, family = "univariate", method = "MLE",
  smoothed = TRUE
)

lk_kde <- predict_logdensity_estimator(x_test, fit_kde, family = "univariate", method = "KDE")
lk_mle <- predict_logdensity_estimator(x_test, fit_mle, family = "univariate", method = "MLE")

summary(lk_kde)
summary(lk_mle)

mean(-lk_kde[is.finite(lk_kde)])
mean(-lk_mle[is.finite(lk_mle)])

sum(!is.finite(lk_kde))
sum(!is.finite(lk_mle))
