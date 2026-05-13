# ============================================================
# Evaluation_Metrics_Clean.R
# Evaluation metrics for 1D / multivariate benchmarks
# Used metrics:
#   - kl
#   - score_loss
# Optional metric variants:
#   - central_* : metric computed only on the empirical bulk region (Be careful in higher dimensions
#                 cut of componentwise -> ~ (1-2c)^d *100% of data remains)
#   - *_trim    : score metric after trimming the largest pointwise losses
# ============================================================

# Expected:
# source("01_Rscripts/02_Helper_and_Metrics/helper_functions.R")
# source("01_Rscripts/01_Estimator/KDE.R")
# source("01_Rscripts/01_Estimator/LogConcaveMLE.R")
# source("01_Rscripts/01_Estimator/Univariate_Polynomial_Score_Matching.R")
# source("01_Rscripts/01_Estimator/Multivariate_Pairwise_Polynomial_Score_Matching.R")

# ------------------------------------------------------------
# (1) Generic estimator wrappers
# ------------------------------------------------------------

# Wrapper for fit an estimator based on method and family (uni or multivariate) tag
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
    if (method == "SM")  return(fit_score_matching_mv(x, ...))
  }

  stop("Unsupported family / method combination.")
}

# Wrapper to predict density given a fit object and test data newx 
# - based on method and family (uni or multivariate) tag
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

# Wrapper to predict logdensity given a fit object and test data newx 
# - based on method and family (uni or multivariate) tag
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

# Wrapper to predict score given a fit object and test data newx 
# - based on method and family (uni or multivariate) tag
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
    if (method == "SM")  return(as.matrix(predict_score_mv(newx, fit, ...)))
  }

  stop("Unsupported family / method combination.")
}

# ------------------------------------------------------------
# (2) Helpers to derive the score and KL loss
# ------------------------------------------------------------

# Helper function to ensure data is in matrix format
as_score_matrix <- function(s, n_expected = NULL) {
  if (is.null(dim(s))) {
    s <- matrix(as.numeric(s), ncol = 1)
  } else {
    s <- as.matrix(s)
  }
  if (!is.null(n_expected) && nrow(s) != n_expected) stop("Unexpected number of score rows.")
  s
}

# Error handling for trim parameter (finite, in [0, 0.5), must be scalar, etc.)
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
# Be aware!!!: In higher dimensions d this may cut off too many samples - (1-2central_trim)^d *100% samples remain
make_central_mask <- function(x, family = c("univariate", "multivariate"), central_trim = NULL) {
  family <- match.arg(family)
  central_trim <- validate_trim_scalar(central_trim, "central_trim")
  # Set everything to true (means remain point when evaluate central_loss) if central_trim is NULL or <=0
  if (is.null(central_trim) || central_trim <= 0) {
    n <- if (family == "univariate") length(as.numeric(x)) else nrow(as.matrix(x))
    return(rep(TRUE, n))
  }

  # Univariate: Save false values for points that should be neglected in central_loss   
  if (family == "univariate") {
    xx <- as.numeric(x)
    lo <- safe_quantile(xx, central_trim)
    hi <- safe_quantile(xx, 1 - central_trim)
    return(is.finite(xx) & xx >= lo & xx <= hi)
  }

  # Multivariate: Save false values for points that should be neglected in central_loss 
  xx <- as.matrix(x)
  keep <- apply(xx, 1, function(row) all(is.finite(row)))
  xx2 <- xx[keep, , drop = FALSE]
  if (nrow(xx2) == 0L) return(rep(FALSE, nrow(xx)))
  lo <- apply(xx2, 2, safe_quantile, p = central_trim)
  hi <- apply(xx2, 2, safe_quantile, p = 1 - central_trim)
  # Create mask and set every value to false that corresponds to a sample point 
  # that violates at least one component (not in quantile range for this component)
  out <- rep(FALSE, nrow(xx))
  out[keep] <- apply(xx2, 1, function(row) all(row >= lo & row <= hi))
  out
}

# robust_trim is interpreted as post-hoc trimming of the largest pointwise score losses.
# It changes the target away from the full Fisher-type average and should therefore be read
# as a robust score diagnostic, not as the untrimmed score loss itself.
trim_top_pointwise_losses <- function(losses, robust_trim = NULL) {
  # Error handling and removing of infinite values
  robust_trim <- validate_trim_scalar(robust_trim, "robust_trim")
  losses <- as.numeric(losses)
  keep <- is.finite(losses)
  if (sum(keep) == 0L) return(losses)
  if (is.null(robust_trim) || robust_trim <= 0) return(losses[keep])
  # if after preprocessing losses still exists apply trimming in quantile range
  cutoff <- stats::quantile(losses[keep], probs = 1 - robust_trim, na.rm = TRUE, names = FALSE)
  losses[keep & losses <= cutoff]
}

# ------------------------------------------------------------
# (3) Pointwise score and KL loss
# ------------------------------------------------------------

# Compute pointwise Kl loss given 
#   - a fit object
#   - family and method tag
#   - rue log_density
#   - optional preditction arguments
compute_pointwise_kl <- function(x_test,
                                 fit,
                                 family,
                                 method,
                                 true_logdensity,
                                 predict_args = list()) {
  if (is.null(true_logdensity)) stop("true_logdensity must be supplied for KL.")
  # Call generic AI for predict log density
  log_hat <- do.call(
    predict_logdensity_estimator_generic,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  # true log density
  log_true <- as.numeric(true_logdensity(x_test))
  # Error handling
  if (length(log_hat) != length(log_true)) {
    stop("Logdensity predictions and true logdensities have incompatible lengths.")
  }
  # Return prediction, true logden and pointwise difference
  list(
    losses = as.numeric(log_true - log_hat),
    log_true = log_true,
    log_hat = as.numeric(log_hat)
  )
}

# Compute pointwise score loss given 
#   - a fit object
#   - family and method tag
#   - rue log_density
#   - optional preditction arguments
compute_pointwise_score_loss <- function(x_test,
                                         fit,
                                         family,
                                         method,
                                         true_score,
                                         predict_args = list()) {
  if (is.null(true_score)) stop("true_score must be supplied for score_loss.")
  # Call generic AI for predict pointwise score
  score_hat <- do.call(
    predict_score_estimator_generic,
    c(list(newx = x_test, fit = fit, family = family, method = method), predict_args)
  )
  # Compute true pointwise score for each component
  score_true <- true_score(x_test)
  n_expected <- if (family == "univariate") {
    length(as.numeric(x_test))
  } else {
    nrow(as.matrix(x_test))
  }
  score_hat <- as_score_matrix(score_hat, n_expected = n_expected)
  score_true <- as_score_matrix(score_true, n_expected = n_expected)
  if (nrow(score_hat) != nrow(score_true) || ncol(score_hat) != ncol(score_true)) {
    stop("Score matrices have incompatible dimensions.")
  }
  # Check if score of a sample is finite in each component
  keep <- apply(score_hat, 1, function(row) all(is.finite(row))) &
    apply(score_true, 1, function(row) all(is.finite(row)))
  # sum up Squared Difference of true and estimated scores in each component for each sample
  err2 <- rowSums((score_hat[keep, , drop = FALSE] - score_true[keep, , drop = FALSE])^2)
  # return pointwise score loss and matrices of pointwise scores
  list(
    losses = as.numeric(err2),
    keep = keep,
    score_hat = score_hat[keep, , drop = FALSE],
    score_true = score_true[keep, , drop = FALSE]
  )
}

# ------------------------------------------------------------
# (4) Save diagnostic data
# ------------------------------------------------------------

# Diagnostics of solver used during fitting the estimator
extract_fit_diagnostics <- function(fit) {
  status <- NA_character_
  iterations <- NA_real_
  objective_value <- NA_real_
  converged <- NA
  
  # Save solver status and solution
  if (!is.null(fit$status)) status <- as.character(fit$status)
  if (!is.null(fit$solution$value)) objective_value <- suppressWarnings(as.numeric(fit$solution$value))

  # Get number of iterations for solving problem
  it_candidates <- c(
    fit$solution$num_iters,
    fit$solution$solver_stats$num_iters,
    fit$solution$solver_stats$iter
  )
  it_candidates <- unlist(it_candidates)
  it_candidates <- it_candidates[is.finite(it_candidates)]
  if (length(it_candidates) > 0L) iterations <- it_candidates[1L]

  # Read optimizer object and use its values if this object exists
  if (!is.null(fit$optimizer)) {
    if (!is.null(fit$optimizer$convergence)) converged <- identical(fit$optimizer$convergence, 0)
    if (!is.null(fit$optimizer$value)) objective_value <- suppressWarnings(as.numeric(fit$optimizer$value))
    if (!is.null(fit$optimizer$counts)) iterations <- sum(unlist(fit$optimizer$counts), na.rm = TRUE)
    if (!is.null(fit$optimizer$message) && is.na(status)) status <- as.character(fit$optimizer$message)
  }

  if (is.na(converged) && !is.na(status)) converged <- !grepl("fail|error|infeasible", tolower(status))
  
  # Save diagnostics in list.
  # Save the conditon number of the regularized matrix if existent else the unregularized 
  list(
    success = ifelse(is.na(converged), NA, converged),
    status = status,
    iterations = iterations,
    objective_value = objective_value,
    condition_number = fit$diagnostics$kappa_reg %||% fit$diagnostics$kappa_raw %||% NA_real_
  )
}

# Only relevant for univariate SM, check if normalization function was successfull
extract_density_diagnostic <- function(fit,
                                       family,
                                       method,
                                       predict_args = list()) {
  if (!(identical(family, "univariate") && identical(method, "SM"))) {
    return(list(normalization_ok = NA, normalization_message = NA_character_))
  }
  # Check if normalizing function exists
  if (!exists("compute_log_normalizer_z_univariate", mode = "function")) {
    return(list(normalization_ok = NA, normalization_message = NA_character_))
  }
  # Create args for normalization based on input
  args <- list(fit = fit)
  for (nm in c("interval", "subdivisions", "rel.tol", "abs.tol", "stop_on_failure")) {
    if (!is.null(predict_args[[nm]])) args[[nm]] <- predict_args[[nm]]
  }
  # catch error if normalization was not succesfull
  info <- tryCatch(do.call(compute_log_normalizer_z_univariate, args), error = function(e) e)
  # Return error if function failed or logZ infinite
  if (inherits(info, "error")) {
    return(list(normalization_ok = FALSE, normalization_message = conditionMessage(info)))
  }
  list(
    normalization_ok = is.finite(info$logZ),
    normalization_message = info$message %||% NA_character_
  )
}

# ------------------------------------------------------------
# (5) Build string for metric names that are used to tag a points
#     in final testing objects
# ------------------------------------------------------------

# Create all combinations of possible tags that might be used in a test
# given the args for density and score estimation
build_metric_names <- function(metrics,
                               density_metric_args = list(),
                               score_metric_args = list()) {
  metrics <- unique(metrics)
  out <- character(0)
  # Either use inout as names or the standard naming for trim objects
  d_trim <- validate_trim_scalar(density_metric_args$central_trim %||% NULL, "central_trim")
  s_central <- validate_trim_scalar(score_metric_args$central_trim %||% NULL, "central_trim")
  s_robust <- validate_trim_scalar(score_metric_args$robust_trim %||% NULL, "robust_trim")

  # Create all combinations out of given input or standard names
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

# ------------------------------------------------------------
# (5) Compute kl and score loss
# ------------------------------------------------------------

# Compute Monte Carlo Kl loss
metric_kl_mc <- function(x_test, fit, family, method, true_logdensity,
                         predict_args = list(), central_trim = NULL) {
  # Use pointwise helper function for KL loss
  pointwise <- compute_pointwise_kl(
    x_test = x_test,
    fit = fit,
    family = family,
    method = method,
    true_logdensity = true_logdensity,
    predict_args = predict_args
  )
  # get mask if valid central parameter was provided and filter pointwise losses
  mask <- make_central_mask(x_test, family = family, central_trim = central_trim)
  losses <- pointwise$losses[mask]
  # Final loss
  safe_mean(losses)
}

# Compute score loss
metric_score_loss <- function(x_test, fit, family, method, true_score,
                              predict_args = list(), central_trim = NULL, robust_trim = NULL) {
  # Use pointwise helper function for score loss
  pointwise <- compute_pointwise_score_loss(
    x_test = x_test,
    fit = fit,
    family = family,
    method = method,
    true_score = true_score,
    predict_args = predict_args
  )
  # Mask of points that are valid (finite in every component)
  base_mask <- pointwise$keep
  # all valid losses
  full_losses <- pointwise$losses
  # Mask of points that are valid according to central trim
  central_mask_all <- make_central_mask(x_test, family = family, central_trim = central_trim)
  central_mask <- central_mask_all[base_mask]

  # Final score loss without central
  out <- list(score_loss = safe_mean(full_losses))
  if (!is.null(validate_trim_scalar(central_trim, "central_trim")) && central_trim > 0) {
    # Final score loss with central
    out$score_loss_central <- safe_mean(full_losses[central_mask])
  }
  if (!is.null(validate_trim_scalar(robust_trim, "robust_trim")) && robust_trim > 0) {
    # Final score loss with trim
    out$score_loss_trim <- safe_mean(trim_top_pointwise_losses(full_losses, robust_trim = robust_trim))
    if (!is.null(validate_trim_scalar(central_trim, "central_trim")) && central_trim > 0) {
      # Final score loss with central trim
      out$score_loss_central_trim <- safe_mean(trim_top_pointwise_losses(full_losses[central_mask], robust_trim = robust_trim))
    }
  }
  out
}

# ------------------------------------------------------------
# (6) Final Function that evaluates requested metrics given
#     - fit, x_test
#     - metric, family and method tags
#     - true density + score function

# ------------------------------------------------------------
evaluate_requested_metrics <- function(metrics,
                                       x_test,
                                       fit,
                                       family,
                                       method,
                                       true_logdensity = NULL,
                                       true_score = NULL,
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
