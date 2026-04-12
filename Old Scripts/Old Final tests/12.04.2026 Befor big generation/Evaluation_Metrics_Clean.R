# ============================================================
# Evaluation_Metrics_Clean.R
# Schlanke Evaluierungsmetriken für 1D / multivariat
# Enthält nur die im aktuellen Studiendesign verwendeten Metriken:
#   - negloglik
#   - kl
#   - score_loss
# sowie einheitliche Fit-/Solver-Diagnostik
# ============================================================

# Erwartet:
# source("helper_functions.R")
# source("KDE.R")
# source("LogConcaveMLE.R")
# source("Univariate_Polynomial_Score_Matching_1.0.R")
# source("Multivariate_Pairwise_Polynomial_Score_Matching.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

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

  pointwise <- compute_pointwise_negloglik(
    x_test = x_test,
    fit = fit,
    family = family,
    method = method,
    predict_args = predict_args
  )

  diag <- robust_pointwise_loss_summary(
    pointwise$losses,
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

  pointwise <- compute_pointwise_kl(
    x_test = x_test,
    fit = fit,
    family = family,
    method = method,
    true_logdensity = true_logdensity,
    predict_args = predict_args
  )

  diag <- robust_pointwise_loss_summary(
    pointwise$losses,
    robust = robust,
    trim_alpha = trim_alpha,
    outlier_dom_threshold = outlier_dom_threshold
  )

  attach_metric_diagnostics(diag$value, diag)
}


# ------------------------------------------------------------
# (4) Score-Metrik
# ------------------------------------------------------------

compute_pointwise_negloglik <- function(x_test,
                                     fit,
                                     family,
                                     method,
                                     predict_args = list()) {
  logdens <- do.call(
    predict_logdensity_estimator_generic,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )

  list(
    losses = -as.numeric(logdens),
    log_hat = as.numeric(logdens)
  )
}

compute_pointwise_kl <- function(x_test,
                                 fit,
                                 family,
                                 method,
                                 true_logdensity,
                                 predict_args = list()) {
  if (is.null(true_logdensity)) {
    stop("true_logdensity must be supplied for KL.")
  }

  log_hat <- do.call(
    predict_logdensity_estimator_generic,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  log_true <- as.numeric(true_logdensity(x_test))

  list(
    losses = as.numeric(log_true - log_hat),
    log_true = log_true,
    log_hat = as.numeric(log_hat)
  )
}

compute_pointwise_score_loss <- function(x_test,
                                         fit,
                                         family,
                                         method,
                                         true_score,
                                         predict_args = list()) {
  if (is.null(true_score)) {
    stop("true_score must be supplied for score_loss.")
  }

  score_hat <- do.call(
    predict_score_estimator_generic,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  score_true <- true_score(x_test)

  score_hat <- as_score_matrix(score_hat)
  score_true <- as_score_matrix(score_true)

  tmp <- clean_complete_cases_pair(score_hat, score_true)
  err2 <- rowSums((tmp$a - tmp$b)^2)

  list(
    losses = as.numeric(err2),
    score_hat = tmp$a,
    score_true = tmp$b,
    keep = tmp$keep
  )
}

metric_score_loss <- function(x_test,
                              fit,
                              family,
                              method,
                              true_score,
                              predict_args = list(),
                              robust = c("none", "trim", "winsor", "median"),
                              trim_alpha = 0.01,
                              center = NULL,
                              outlier_dom_threshold = 0.25) {
  robust <- match.arg(robust)

  if (!is.null(center)) {
    warning("Argument 'center' is deprecated for metric_score_loss() and will be ignored.", call. = FALSE)
  }

  pointwise <- compute_pointwise_score_loss(
    x_test = x_test,
    fit = fit,
    family = family,
    method = method,
    true_score = true_score,
    predict_args = predict_args
  )

  diag <- robust_pointwise_loss_summary(
    pointwise$losses,
    robust = robust,
    trim_alpha = trim_alpha,
    outlier_dom_threshold = outlier_dom_threshold
  )

  attach_metric_diagnostics(diag$value, diag)
}

normalize_metric_configurations <- function(metric_name,
                                            metric_args = list(),
                                            default_name = "default") {
  metric_args <- metric_args %||% list()

  if (!is.null(metric_args$configs)) {
    cfgs <- metric_args$configs
    if (!is.list(cfgs) || length(cfgs) == 0L) {
      stop(sprintf("%s metric_args$configs must be a non-empty list.", metric_name))
    }
    base_args <- metric_args
    base_args$configs <- NULL
    if (is.null(names(cfgs)) || any(names(cfgs) == "")) {
      names(cfgs) <- paste0(default_name, seq_along(cfgs))
    }
    return(lapply(cfgs, function(cfg) c(base_args, cfg %||% list())))
  }

  list(default = metric_args)
}

metric_variant_name <- function(metric, config_name, n_configs) {
  if (n_configs == 1L && identical(config_name, "default")) {
    return(metric)
  }
  paste0(metric, "__", config_name)
}

# ------------------------------------------------------------
# (5) Numerik / Solver-Diagnostik
# ------------------------------------------------------------

extract_fit_diagnostics <- function(fit) {
  status <- NA_character_
  iterations <- NA_real_
  objective_value <- NA_real_
  converged <- NA

  if (!is.null(fit$status)) {
    status <- as.character(fit$status)
  }
  if (!is.null(fit$solution$value)) {
    objective_value <- suppressWarnings(as.numeric(fit$solution$value))
  }

  it_candidates <- c(
    fit$solution$num_iters,
    fit$solution$solver_stats$num_iters,
    fit$solution$solver_stats$iter
  )
  it_candidates <- unlist(it_candidates)
  it_candidates <- it_candidates[is.finite(it_candidates)]
  if (length(it_candidates) > 0L) iterations <- it_candidates[1L]

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
    fit_time_median = stats::median(df$fit_time_sec, na.rm = TRUE),
    inference_time_median = stats::median(df$total_inference_time_sec, na.rm = TRUE),
    success_rate = mean(df$success, na.rm = TRUE),
    failure_rate = mean(!df$success, na.rm = TRUE),
    iterations_median = stats::median(df$iterations, na.rm = TRUE),
    objective_median = stats::median(df$objective_value, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# (6) Metrik-Dispatcher
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

  for (met in unique(metrics)) {
    if (met %in% c("negloglik", "kl")) {
      cfgs <- normalize_metric_configurations(met, density_metric_args)
      n_cfg <- length(cfgs)

      for (cfg_name in names(cfgs)) {
        val <- if (met == "negloglik") {
          do.call(
            metric_neg_loglik,
            c(
              list(
                x_test = x_test,
                fit = fit,
                family = family,
                method = method,
                predict_args = density_predict_args
              ),
              cfgs[[cfg_name]]
            )
          )
        } else {
          do.call(
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
              cfgs[[cfg_name]]
            )
          )
        }

        metric_name <- metric_variant_name(met, cfg_name, n_cfg)
        out[[metric_name]] <- strip_metric_value(val)
        out[[paste0(metric_name, "_na_share")]] <- attr(val, "na_share")
        out[[paste0(metric_name, "_tail_share")]] <- attr(val, "tail_share")
        out[[paste0(metric_name, "_outlier_dominated")]] <- attr(val, "outlier_dominated")
      }

    } else if (met == "score_loss") {
      cfgs <- normalize_metric_configurations(met, score_metric_args)
      n_cfg <- length(cfgs)

      for (cfg_name in names(cfgs)) {
        val <- do.call(
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
            cfgs[[cfg_name]]
          )
        )

        metric_name <- metric_variant_name(met, cfg_name, n_cfg)
        out[[metric_name]] <- strip_metric_value(val)
        out[[paste0(metric_name, "_na_share")]] <- attr(val, "na_share")
        out[[paste0(metric_name, "_tail_share")]] <- attr(val, "tail_share")
        out[[paste0(metric_name, "_outlier_dominated")]] <- attr(val, "outlier_dominated")
      }
    } else {
      stop(sprintf("Unknown metric: %s", met))
    }
  }

  out
}
