# ============================================================
# evaluation_metrics.R
# Generische Evaluierungsmetriken für 1D / multivariat
# ============================================================

# Erwartet:
# source("helper_functions.R")
# source("KDE.R")
# source("LogConcaveMLE.R")
# source("Univariate_Polynomial_Score_Matching_1.0.R")
# source("Multivariate_Pairwise_Polynomial_Score_Matching.R")

# ------------------------------------------------------------
# (1) Generische Estimator-Wrapper
# ------------------------------------------------------------

fit_estimator_generic <- function(x,
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

predict_density_estimator_generic <- function(newx,
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
    if (method == "SM") {
      stop("Density metrics are currently not supported for multivariate SM.")
    }
  }
  
  stop("Unsupported family / method combination.")
}

predict_logdensity_estimator_generic <- function(newx,
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
    if (method == "SM") {
      stop("Density metrics are currently not supported for multivariate SM.")
    }
  }
  
  stop("Unsupported family / method combination.")
}

predict_score_estimator_generic <- function(newx,
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

# ------------------------------------------------------------
# (2) Hilfsfunktionen
# ------------------------------------------------------------

as_score_matrix <- function(s, n_expected = NULL) {
  if (is.null(dim(s))) {
    s <- matrix(as.numeric(s), ncol = 1)
  } else {
    s <- as.matrix(s)
  }
  
  if (!is.null(n_expected) && nrow(s) != n_expected) {
    stop("Unexpected number of score rows.")
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

robust_pointwise_loss_summary <- function(x,
                                          robust = c("none", "trim", "winsor", "median"),
                                          trim_alpha = 0.01,
                                          outlier_dom_threshold = 0.25) {
  robust <- match.arg(robust)
  
  x_raw <- as.numeric(x)
  is_bad <- !is.finite(x_raw)
  x <- x_raw[!is_bad]
  
  na_share <- mean(is_bad)
  
  if (length(x) == 0L) {
    return(list(
      value = NA_real_,
      na_share = na_share,
      tail_share = NA_real_,
      outlier_dominated = NA_real_
    ))
  }
  
  full_mean <- mean(x)
  
  cutoff <- stats::quantile(
    x,
    probs = 1 - trim_alpha,
    na.rm = TRUE,
    names = FALSE
  )
  
  tail_idx <- x > cutoff
  
  tail_share <- if (sum(x, na.rm = TRUE) > 0) {
    sum(x[tail_idx], na.rm = TRUE) / sum(x, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  value <- switch(
    robust,
    none = full_mean,
    trim = {
      x2 <- x[!tail_idx]
      if (length(x2) == 0L) NA_real_ else mean(x2)
    },
    winsor = mean(pmin(x, cutoff)),
    median = stats::median(x)
  )
  
  rel_gap <- if (is.finite(full_mean) && is.finite(value)) {
    abs(full_mean - value) / max(abs(full_mean), 1e-12)
  } else {
    NA_real_
  }
  
  outlier_dominated <- if (is.finite(rel_gap)) rel_gap > outlier_dom_threshold else NA_real_
  
  list(
    value = value,
    na_share = na_share,
    tail_share = tail_share,
    outlier_dominated = as.numeric(outlier_dominated)
  )
}

attach_metric_diagnostics <- function(value, diag) {
  attr(value, "na_share") <- diag$na_share
  attr(value, "tail_share") <- diag$tail_share
  attr(value, "outlier_dominated") <- diag$outlier_dominated
  value
}

strip_metric_value <- function(x) {
  out <- as.numeric(x)[1]
  attributes(out) <- NULL
  out
}

# ------------------------------------------------------------
# (3) Dichte-Metriken
# ------------------------------------------------------------

metric_neg_loglik <- function(x_test, fit, family, method,
                              predict_args = list(),
                              robust = c("none", "trim", "winsor", "median"),
                              trim_alpha = 0.01,
                              outlier_dom_threshold = 0.25) {
  robust <- match.arg(robust)
  
  logdens <- do.call(
    predict_logdensity_estimator_generic,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  
  losses <- -as.numeric(logdens)
  
  diag <- robust_pointwise_loss_summary(
    losses,
    robust = robust,
    trim_alpha = trim_alpha,
    outlier_dom_threshold = outlier_dom_threshold
  )
  
  attach_metric_diagnostics(diag$value, diag)
}

metric_kl_mc <- function(x_test, fit, family, method, true_logdensity,
                         predict_args = list(),
                         robust = c("none", "trim", "winsor", "median"),
                         trim_alpha = 0.01,
                         outlier_dom_threshold = 0.25) {
  robust <- match.arg(robust)
  
  log_hat <- do.call(
    predict_logdensity_estimator_generic,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  log_true <- as.numeric(true_logdensity(x_test))
  
  losses <- log_true - log_hat
  
  diag <- robust_pointwise_loss_summary(
    losses,
    robust = robust,
    trim_alpha = trim_alpha,
    outlier_dom_threshold = outlier_dom_threshold
  )
  
  attach_metric_diagnostics(diag$value, diag)
}

metric_hellinger2_mc <- function(x_test, fit, family, method, true_density,
                                 predict_args = list()) {
  f_hat  <- do.call(
    predict_density_estimator_generic,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  f_true <- as.numeric(true_density(x_test))
  keep <- is.finite(f_hat) & is.finite(f_true) & (f_hat >= 0) & (f_true >= 0)
  0.5 * safe_mean((sqrt(f_hat[keep]) - sqrt(f_true[keep]))^2 / pmax(f_true[keep], 1e-300) * f_true[keep])
}

metric_ise_1d <- function(grid, fit, family, method, true_density,
                          predict_args = list()) {
  if (family != "univariate") stop("ISE currently implemented only for 1D grids.")
  grid <- as.numeric(grid)
  
  f_hat <- do.call(
    predict_density_estimator_generic,
    c(list(newx = grid, fit = fit, family = family, method = method), predict_args)
  )
  f_true <- as.numeric(true_density(grid))
  
  trapz_1d(grid, (f_hat - f_true)^2)
}

# ------------------------------------------------------------
# (4) Score-Metriken
# ------------------------------------------------------------

metric_score_loss <- function(x_test, fit, family, method, true_score,
                              predict_args = list(),
                              robust = c("none", "trim", "winsor"),
                              trim_alpha = 0.01) {
  robust <- match.arg(robust)
  
  score_hat <- do.call(
    predict_score_estimator_generic,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  score_true <- true_score(x_test)
  
  score_hat <- as_score_matrix(score_hat)
  score_true <- as_score_matrix(score_true)
  
  tmp <- clean_complete_cases_pair(score_hat, score_true)
  err2 <- rowSums((tmp$a - tmp$b)^2)
  err2 <- err2[is.finite(err2)]
  
  if (length(err2) == 0L) return(NA_real_)
  
  if (robust == "trim") {
    cutoff <- stats::quantile(err2, probs = 1 - trim_alpha, na.rm = TRUE, names = FALSE)
    err2 <- err2[err2 <= cutoff]
  } else if (robust == "winsor") {
    cutoff <- stats::quantile(err2, probs = 1 - trim_alpha, na.rm = TRUE, names = FALSE)
    err2 <- pmin(err2, cutoff)
  }
  
  mean(err2)
}

metric_fisher_divergence <- function(x_test, fit, family, method, true_score,
                                     predict_args = list(), half = FALSE) {
  out <- metric_score_loss(
    x_test = x_test,
    fit = fit,
    family = family,
    method = method,
    true_score = true_score,
    predict_args = predict_args
  )
  if (half) 0.5 * out else out
}

metric_score_rmse <- function(x_test, fit, family, method, true_score,
                              predict_args = list(),
                              robust = c("none", "trim", "winsor"),
                              trim_alpha = 0.01) {
  sqrt(metric_score_loss(
    x_test = x_test,
    fit = fit,
    family = family,
    method = method,
    true_score = true_score,
    predict_args = predict_args,
    robust = robust,
    trim_alpha = trim_alpha
  ))
}

# ------------------------------------------------------------
# (5) Log-Konkavitäts- / Form-Metriken
# ------------------------------------------------------------

metric_logconcavity_mv_from_fit <- function(fit) {
  if (!inherits(fit, "score_matching_mv_fit")) {
    return(list(
      min_hessian_eigenvalue = NA_real_,
      share_violated = NA_real_,
      mean_violation = NA_real_,
      max_violation = NA_real_
    ))
  }
  
  diag <- fit$lc_diagnostics
  if (is.null(diag)) {
    return(list(
      min_hessian_eigenvalue = NA_real_,
      share_violated = NA_real_,
      mean_violation = NA_real_,
      max_violation = NA_real_
    ))
  }
  
  n_pts <- length(diag$min_eigenvalues)
  list(
    min_hessian_eigenvalue = min(diag$min_eigenvalues, na.rm = TRUE),
    share_violated = if (n_pts > 0L) diag$n_violated / n_pts else NA_real_,
    mean_violation = diag$mean_violation,
    max_violation = diag$max_violation
  )
}

num_second_derivative_1d <- function(fun, x, h = 1e-4) {
  (fun(x + h) - 2 * fun(x) + fun(x - h)) / (h^2)
}

metric_logconcavity_1d_grid <- function(grid, fit, family, method,
                                        predict_args = list(),
                                        tol = 1e-8,
                                        h = 1e-4) {
  if (family != "univariate") {
    stop("metric_logconcavity_1d_grid is only for 1D.")
  }
  
  f_log <- function(xx) {
    do.call(
      predict_logdensity_estimator_generic,
      c(list(newx = xx, fit = fit, family = family, method = method), predict_args)
    )
  }
  
  d2 <- vapply(grid, function(xx) num_second_derivative_1d(f_log, xx, h = h), numeric(1))
  violations <- pmax(d2 - tol, 0)
  
  list(
    min_hessian_eigenvalue = min(-d2, na.rm = TRUE),  # in 1D: -logf'' als "Formstärke"
    share_violated = mean(violations > 0, na.rm = TRUE),
    mean_violation = mean(violations, na.rm = TRUE),
    max_violation = max(violations, na.rm = TRUE)
  )
}

# ------------------------------------------------------------
# (6) Numerik / Solver-Diagnostik
# ------------------------------------------------------------

extract_fit_diagnostics <- function(fit) {
  status <- NA_character_
  iterations <- NA_real_
  objective_value <- NA_real_
  converged <- NA
  
  # univariates SM
  if (!is.null(fit$status)) {
    status <- as.character(fit$status)
  }
  if (!is.null(fit$solution$value)) {
    objective_value <- suppressWarnings(as.numeric(fit$solution$value))
  }
  
  # CVXR/SCS zählt Iterationen nicht immer gleich zugänglich;
  # wir lesen mehrere mögliche Felder robust aus
  it_candidates <- c(
    fit$solution$num_iters,
    fit$solution$solver_stats$num_iters,
    fit$solution$solver_stats$iter
  )
  it_candidates <- unlist(it_candidates)
  it_candidates <- it_candidates[is.finite(it_candidates)]
  if (length(it_candidates) > 0L) iterations <- it_candidates[1L]
  
  # multivariates SM
  if (!is.null(fit$optimizer)) {
    if (!is.null(fit$optimizer$convergence)) {
      converged <- identical(fit$optimizer$convergence, 0)
    }
    if (!is.null(fit$optimizer$value)) {
      objective_value <- suppressWarnings(as.numeric(fit$optimizer$value))
    }
    if (!is.null(fit$optimizer$counts)) {
      iterations <- sum(unlist(fit$optimizer$counts), na.rm = TRUE)
    }
    if (!is.null(fit$optimizer$message) && is.na(status)) {
      status <- as.character(fit$optimizer$message)
    }
  }
  
  if (is.na(converged)) {
    if (!is.na(status)) {
      converged <- !grepl("fail|error|infeasible", tolower(status))
    }
  }
  
  list(
    success = ifelse(is.na(converged), NA, converged),
    status = status,
    iterations = iterations,
    objective_value = objective_value,
    kappa_raw = fit$diagnostics$kappa_raw %||% NA_real_,
    kappa_reg = fit$diagnostics$kappa_reg %||% NA_real_,
    rcond_raw = fit$diagnostics$rcond_raw %||% NA_real_,
    rcond_reg = fit$diagnostics$rcond_reg %||% NA_real_,
    eigmin_raw = fit$diagnostics$eigmin_raw %||% NA_real_,
    eigmin_reg = fit$diagnostics$eigmin_reg %||% NA_real_
  )
}

summarize_numerics <- function(df) {
  data.frame(
    runtime_median = stats::median(df$runtime_sec, na.rm = TRUE),
    runtime_IQR = stats::IQR(df$runtime_sec, na.rm = TRUE),
    success_rate = mean(df$success, na.rm = TRUE),
    failure_rate = mean(!df$success, na.rm = TRUE),
    iterations_median = stats::median(df$iterations, na.rm = TRUE),
    objective_median = stats::median(df$objective_value, na.rm = TRUE)
  )
}

# ------------------------------------------------------------
# (7) Metrik-Dispatcher
# ------------------------------------------------------------

evaluate_requested_metrics <- function(metrics,
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
                                       score_metric_args = list()) {
  out <- list()
  
  for (met in metrics) {
    if (met == "negloglik") {
      val <- do.call(
        metric_neg_loglik,
        c(
          list(
            x_test = x_test,
            fit = fit,
            family = family,
            method = method,
            predict_args = density_predict_args
          ),
          density_metric_args
        )
      )
      out[[met]] <- strip_metric_value(val)
      out[[paste0(met, "_na_share")]] <- attr(val, "na_share")
      out[[paste0(met, "_tail_share")]] <- attr(val, "tail_share")
      out[[paste0(met, "_outlier_dominated")]] <- attr(val, "outlier_dominated")
      
    } else if (met == "kl") {
      val <- do.call(
        metric_kl_mc,
        c(
          list(
            x_test = x_test,
            fit = fit,
            family = family,
            method = method,
            true_logdensity = true_logdensity,
            predict_args = density_predict_args
          ),
          density_metric_args
        )
      )
      out[[met]] <- strip_metric_value(val)
      out[[paste0(met, "_na_share")]] <- attr(val, "na_share")
      out[[paste0(met, "_tail_share")]] <- attr(val, "tail_share")
      out[[paste0(met, "_outlier_dominated")]] <- attr(val, "outlier_dominated")
    } else if (met == "hellinger2") {
      out[[met]] <- metric_hellinger2_mc(
        x_test, fit, family, method, true_density, density_predict_args
      )
    } else if (met == "ise") {
      out[[met]] <- metric_ise_1d(
        grid_1d, fit, family, method, true_density, density_predict_args
      )
    } else if (met == "score_loss") {
      out[[met]] <- do.call(
        metric_score_loss,
        c(
          list(
            x_test = x_test,
            fit = fit,
            family = family,
            method = method,
            true_score = true_score,
            predict_args = score_predict_args
          ),
          score_metric_args
        )
      )
    } else if (met == "fisher") {
      out[[met]] <- metric_fisher_divergence(
        x_test, fit, family, method, true_score, score_predict_args
      )
    } else if (met == "score_rmse") {
      out[[met]] <- do.call(
        metric_score_rmse,
        c(
          list(
            x_test = x_test,
            fit = fit,
            family = family,
            method = method,
            true_score = true_score,
            predict_args = score_predict_args
          ),
          score_metric_args
        )
      )
    } else if (met == "logconcavity") {
      if (family == "univariate") {
        lc <- metric_logconcavity_1d_grid(
          grid = grid_1d,
          fit = fit,
          family = family,
          method = method,
          predict_args = density_predict_args
        )
      } else {
        lc <- metric_logconcavity_mv_from_fit(fit)
      }
      
      out[["min_hessian_eigenvalue"]] <- lc$min_hessian_eigenvalue
      out[["share_violated"]] <- lc$share_violated
      out[["mean_violation"]] <- lc$mean_violation
      out[["max_violation"]] <- lc$max_violation
    } else {
      stop(sprintf("Unknown metric: %s", met))
    }
  }
  
  out
}