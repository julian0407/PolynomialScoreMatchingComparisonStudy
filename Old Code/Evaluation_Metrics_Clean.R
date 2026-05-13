# ============================================================
# Evaluation_Metrics_Clean.R
# Lean evaluation metrics for 1D / multivariate benchmarks
# Used metrics:
#   - kl
#   - score_loss
# Optional metric variants:
#   - central_* : metric computed only on the empirical bulk region
#   - *_trim    : score metric after trimming the largest pointwise losses
# ============================================================

# Expected:
# source("helper_functions.R")
# source("KDE.R")
# source("LogConcaveMLE.R")
# source("Univariate_Polynomial_Score_Matching_1.0.R")
# source("Multivariate_Pairwise_Polynomial_Score_Matching.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------
# (1) Generic estimator wrappers
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
    if (method == "SM") stop("Density metrics are currently not supported for multivariate SM.")
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
    if (method == "SM") stop("Density metrics are currently not supported for multivariate SM.")
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
# (2) Helpers
# ------------------------------------------------------------
as_score_matrix <- function(s, n_expected = NULL) {
  if (is.null(dim(s))) {
    s <- matrix(as.numeric(s), ncol = 1)
  } else {
    s <- as.matrix(s)
  }
  if (!is.null(n_expected) && nrow(s) != n_expected) stop("Unexpected number of score rows.")
  s
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

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

validate_trim_scalar <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x >= 0.5) {
    stop(sprintf("%s must be NULL or a scalar in [0, 0.5).", name))
  }
  as.numeric(x)
}

# central_trim is interpreted as evaluation on the empirical bulk region of the test sample.
# In 1D: keep points between empirical alpha and 1-alpha quantiles.
# In multivariate settings: keep rows that lie within these marginal quantiles in every coordinate.
make_central_mask <- function(x, family = c("univariate", "multivariate"), central_trim = NULL) {
  family <- match.arg(family)
  central_trim <- validate_trim_scalar(central_trim, "central_trim")
  if (is.null(central_trim) || central_trim <= 0) {
    n <- if (family == "univariate") length(as.numeric(x)) else nrow(as.matrix(x))
    return(rep(TRUE, n))
  }

  if (family == "univariate") {
    xx <- as.numeric(x)
    lo <- safe_quantile(xx, central_trim)
    hi <- safe_quantile(xx, 1 - central_trim)
    return(is.finite(xx) & xx >= lo & xx <= hi)
  }

  xx <- as.matrix(x)
  keep <- apply(xx, 1, function(row) all(is.finite(row)))
  xx2 <- xx[keep, , drop = FALSE]
  if (nrow(xx2) == 0L) return(rep(FALSE, nrow(xx)))
  lo <- apply(xx2, 2, safe_quantile, p = central_trim)
  hi <- apply(xx2, 2, safe_quantile, p = 1 - central_trim)
  out <- rep(FALSE, nrow(xx))
  out[keep] <- apply(xx2, 1, function(row) all(row >= lo & row <= hi))
  out
}

# robust_trim is interpreted as post-hoc trimming of the largest pointwise score losses.
# It changes the target away from the full Fisher-type average and should therefore be read
# as a robust score diagnostic, not as the untrimmed score loss itself.
trim_top_pointwise_losses <- function(losses, robust_trim = NULL) {
  robust_trim <- validate_trim_scalar(robust_trim, "robust_trim")
  losses <- as.numeric(losses)
  keep <- is.finite(losses)
  if (sum(keep) == 0L) return(losses)
  if (is.null(robust_trim) || robust_trim <= 0) return(losses[keep])
  cutoff <- stats::quantile(losses[keep], probs = 1 - robust_trim, na.rm = TRUE, names = FALSE)
  losses[keep & losses <= cutoff]
}

compute_pointwise_kl <- function(x_test,
                                 fit,
                                 family,
                                 method,
                                 true_logdensity,
                                 predict_args = list()) {
  if (is.null(true_logdensity)) stop("true_logdensity must be supplied for KL.")
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
  if (is.null(true_score)) stop("true_score must be supplied for score_loss.")
  score_hat <- do.call(
    predict_score_estimator_generic,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  score_true <- true_score(x_test)
  score_hat <- as_score_matrix(score_hat)
  score_true <- as_score_matrix(score_true)
  if (nrow(score_hat) != nrow(score_true) || ncol(score_hat) != ncol(score_true)) {
    stop("Score matrices have incompatible dimensions.")
  }
  keep <- apply(score_hat, 1, function(row) all(is.finite(row))) &
    apply(score_true, 1, function(row) all(is.finite(row)))
  err2 <- rowSums((score_hat[keep, , drop = FALSE] - score_true[keep, , drop = FALSE])^2)
  list(
    losses = as.numeric(err2),
    keep = keep,
    score_hat = score_hat[keep, , drop = FALSE],
    score_true = score_true[keep, , drop = FALSE]
  )
}

extract_fit_diagnostics <- function(fit) {
  status <- NA_character_
  iterations <- NA_real_
  objective_value <- NA_real_
  converged <- NA

  if (!is.null(fit$status)) status <- as.character(fit$status)
  if (!is.null(fit$solution$value)) objective_value <- suppressWarnings(as.numeric(fit$solution$value))

  it_candidates <- c(
    fit$solution$num_iters,
    fit$solution$solver_stats$num_iters,
    fit$solution$solver_stats$iter
  )
  it_candidates <- unlist(it_candidates)
  it_candidates <- it_candidates[is.finite(it_candidates)]
  if (length(it_candidates) > 0L) iterations <- it_candidates[1L]

  if (!is.null(fit$optimizer)) {
    if (!is.null(fit$optimizer$convergence)) converged <- identical(fit$optimizer$convergence, 0)
    if (!is.null(fit$optimizer$value)) objective_value <- suppressWarnings(as.numeric(fit$optimizer$value))
    if (!is.null(fit$optimizer$counts)) iterations <- sum(unlist(fit$optimizer$counts), na.rm = TRUE)
    if (!is.null(fit$optimizer$message) && is.na(status)) status <- as.character(fit$optimizer$message)
  }

  if (is.na(converged) && !is.na(status)) converged <- !grepl("fail|error|infeasible", tolower(status))

  list(
    success = ifelse(is.na(converged), NA, converged),
    status = status,
    iterations = iterations,
    objective_value = objective_value,
    condition_number = fit$diagnostics$kappa_reg %||% fit$diagnostics$kappa_raw %||% NA_real_
  )
}

extract_density_diagnostic <- function(fit,
                                       family,
                                       method,
                                       predict_args = list()) {
  if (!(identical(family, "univariate") && identical(method, "SM"))) {
    return(list(normalization_ok = NA, normalization_message = NA_character_))
  }
  if (!exists("compute_log_normalizer_z_univariate", mode = "function")) {
    return(list(normalization_ok = NA, normalization_message = NA_character_))
  }
  args <- list(fit = fit)
  for (nm in c("interval", "subdivisions", "rel.tol", "abs.tol", "stop_on_failure")) {
    if (!is.null(predict_args[[nm]])) args[[nm]] <- predict_args[[nm]]
  }
  info <- tryCatch(do.call(compute_log_normalizer_z_univariate, args), error = function(e) e)
  if (inherits(info, "error")) {
    return(list(normalization_ok = FALSE, normalization_message = conditionMessage(info)))
  }
  list(
    normalization_ok = is.finite(info$logZ),
    normalization_message = info$message %||% NA_character_
  )
}

build_metric_names <- function(metrics,
                               density_metric_args = list(),
                               score_metric_args = list()) {
  metrics <- unique(metrics)
  out <- character(0)
  d_trim <- validate_trim_scalar(density_metric_args$central_trim %||% NULL, "central_trim")
  s_central <- validate_trim_scalar(score_metric_args$central_trim %||% NULL, "central_trim")
  s_robust <- validate_trim_scalar(score_metric_args$robust_trim %||% NULL, "robust_trim")

  if ("kl" %in% metrics) {
    out <- c(out, "kl")
    if (!is.null(d_trim) && d_trim > 0) out <- c(out, "kl_central")
  }
  if ("score_loss" %in% metrics) {
    out <- c(out, "score_loss")
    if (!is.null(s_central) && s_central > 0) out <- c(out, "score_loss_central")
    if (!is.null(s_robust) && s_robust > 0) out <- c(out, "score_loss_trim")
    if ((!is.null(s_central) && s_central > 0) && (!is.null(s_robust) && s_robust > 0)) {
      out <- c(out, "score_loss_central_trim")
    }
  }
  unique(out)
}

metric_kl_mc <- function(x_test, fit, family, method, true_logdensity,
                         predict_args = list(), central_trim = NULL) {
  pointwise <- compute_pointwise_kl(
    x_test = x_test,
    fit = fit,
    family = family,
    method = method,
    true_logdensity = true_logdensity,
    predict_args = predict_args
  )
  mask <- make_central_mask(x_test, family = family, central_trim = central_trim)
  losses <- pointwise$losses[mask]
  safe_mean(losses)
}

metric_score_loss <- function(x_test, fit, family, method, true_score,
                              predict_args = list(), central_trim = NULL, robust_trim = NULL) {
  pointwise <- compute_pointwise_score_loss(
    x_test = x_test,
    fit = fit,
    family = family,
    method = method,
    true_score = true_score,
    predict_args = predict_args
  )
  base_mask <- pointwise$keep
  full_losses <- pointwise$losses
  central_mask_all <- make_central_mask(x_test, family = family, central_trim = central_trim)
  central_mask <- central_mask_all[base_mask]

  out <- list(score_loss = safe_mean(full_losses))
  if (!is.null(validate_trim_scalar(central_trim, "central_trim")) && central_trim > 0) {
    out$score_loss_central <- safe_mean(full_losses[central_mask])
  }
  if (!is.null(validate_trim_scalar(robust_trim, "robust_trim")) && robust_trim > 0) {
    out$score_loss_trim <- safe_mean(trim_top_pointwise_losses(full_losses, robust_trim = robust_trim))
    if (!is.null(validate_trim_scalar(central_trim, "central_trim")) && central_trim > 0) {
      out$score_loss_central_trim <- safe_mean(trim_top_pointwise_losses(full_losses[central_mask], robust_trim = robust_trim))
    }
  }
  out
}

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
  if ("kl" %in% metrics) {
    out$kl <- metric_kl_mc(
      x_test = x_test,
      fit = fit,
      family = family,
      method = method,
      true_logdensity = true_logdensity,
      predict_args = density_predict_args,
      central_trim = NULL
    )
    d_trim <- validate_trim_scalar(density_metric_args$central_trim %||% NULL, "central_trim")
    if (!is.null(d_trim) && d_trim > 0) {
      out$kl_central <- metric_kl_mc(
        x_test = x_test,
        fit = fit,
        family = family,
        method = method,
        true_logdensity = true_logdensity,
        predict_args = density_predict_args,
        central_trim = d_trim
      )
    }
  }
  if ("score_loss" %in% metrics) {
    out <- c(out, metric_score_loss(
      x_test = x_test,
      fit = fit,
      family = family,
      method = method,
      true_score = true_score,
      predict_args = score_predict_args,
      central_trim = score_metric_args$central_trim %||% NULL,
      robust_trim = score_metric_args$robust_trim %||% NULL
    ))
  }
  out
}
