# ============================================================
# ScoreLossTest.R
# Standardized score-loss benchmark for univariate / multivariate
# estimators across sample sizes.
#
# Expected sourced files before use:
source("helper_functions.R")
source("KDE.R")
source("LogConcaveMLE.R")
source("Univariate_Polynomial_Score_Matching_1.0.R")   # only if method = "SM"
#
# Score convention:
#   This script assumes ALL estimator score functions return
#   score(x) = - grad log f(x)
#   and that true_score(x) uses the SAME convention.
# ============================================================

# ------------------------------------------------------------
# (1) Generic estimator wrappers
# ------------------------------------------------------------

fit_estimator <- function(x,
                          family = c("univariate", "multivariate"),
                          method = c("MLE", "KDE", "SM"),
                          smoothed = FALSE,
                          ...) {
  family <- match.arg(family)
  method <- match.arg(method)
  
  if (family == "univariate") {
    if (method == "KDE") {
      return(fit_kde_1d(x, ...))
    }
    if (method == "MLE") {
      return(fit_logconcave_mle_1d(x, smoothed = smoothed, ...))
    }
    if (method == "SM") {
      return(fit_score_matching_univariate(x, ...))
    }
  }
  
  if (family == "multivariate") {
    if (method == "KDE") {
      return(fit_kde_mv(x, ...))
    }
    if (method == "MLE") {
      return(fit_logconcave_mle_mv(x, smoothed = smoothed, ...))
    }
    if (method == "SM") {
      return(fit_score_matching_mv_basic(x, ...))
    }
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
    if (method == "KDE") {
      return(as.numeric(predict_score_kde_1d(newx, fit, ...)))
    }
    if (method == "MLE") {
      return(as.numeric(predict_score_logconcave_1d(newx, fit, ...)))
    }
    if (method == "SM") {
      return(as.numeric(predict_score_univariate(newx, fit)))
    }
  }
  
  if (family == "multivariate") {
    if (method == "KDE") {
      return(as.matrix(predict_score_kde_mv(newx, fit, ...)))
    }
    if (method == "MLE") {
      return(as.matrix(predict_score_logconcave_mv(newx, fit, ...)))
    }
    if (method == "SM") {
      return(as.matrix(predict_score_mv_basic(newx, fit)))
    }
  }
  
  stop("Unsupported family / method combination.")
}

# ------------------------------------------------------------
# (2) Helpers for shape / format consistency
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# (3) Score loss
# ------------------------------------------------------------
# Monte Carlo approximation of
#   E || \hat s(X) - s(X) ||^2
# based on an independent test sample X_1,...,X_m.
# For d = 1 this reduces to mean((shat - strue)^2).

score_loss_mc <- function(score_hat, score_true) {
  tmp <- clean_complete_cases_pair(score_hat, score_true)
  err <- tmp$a - tmp$b
  mean(rowSums(err^2))
}

# Optional variant with 1/2 factor as often used in score matching style losses
score_loss_mc_half <- function(score_hat, score_true) {
  0.5 * score_loss_mc(score_hat, score_true)
}

# ------------------------------------------------------------
# (4) One benchmark run for one n
# ------------------------------------------------------------

run_one_score_loss_experiment <- function(n,
                                          family = c("univariate", "multivariate"),
                                          method = c("MLE", "KDE", "SM"),
                                          r_sample,
                                          true_score,
                                          n_test = 2000,
                                          smoothed = FALSE,
                                          loss_half = FALSE,
                                          fit_args = list(),
                                          predict_args = list()) {
  family <- match.arg(family)
  method <- match.arg(method)
  
  # Generate training sample
  x_train <- r_sample(n)
  
  if (family == "univariate") {
    x_train <- as.numeric(x_train)
  } else {
    x_train <- as.matrix(x_train)
  }
  
  # Fit estimator
  fit <- do.call(
    fit_estimator,
    c(
      list(
        x = x_train,
        family = family,
        method = method,
        smoothed = smoothed
      ),
      fit_args
    )
  )
  
  # Generate independent test sample
  x_test <- r_sample(n_test)
  
  if (family == "univariate") {
    x_test <- as.numeric(x_test)
    
    score_hat <- do.call(
      predict_score_estimator,
      c(
        list(
          newx = x_test,
          fit = fit,
          family = family,
          method = method
        ),
        predict_args
      )
    )
    
    score_true <- true_score(x_test)
    
    score_hat <- as_score_matrix(score_hat, n_expected = length(x_test))
    score_true <- as_score_matrix(score_true, n_expected = length(x_test))
    
  } else {
    x_test <- as.matrix(x_test)
    
    score_hat <- do.call(
      predict_score_estimator,
      c(
        list(
          newx = x_test,
          fit = fit,
          family = family,
          method = method
        ),
        predict_args
      )
    )
    
    score_true <- true_score(x_test)
    
    score_hat <- as_score_matrix(score_hat, n_expected = nrow(x_test))
    score_true <- as_score_matrix(score_true, n_expected = nrow(x_test))
  }
  
  loss_value <- if (loss_half) {
    score_loss_mc_half(score_hat, score_true)
  } else {
    score_loss_mc(score_hat, score_true)
  }
  
  list(
    n = n,
    family = family,
    method = method,
    smoothed = smoothed,
    fit = fit,
    score_loss = loss_value,
    n_test = n_test
  )
}

# ------------------------------------------------------------
# (5) Benchmark over many sample sizes and repetitions
# ------------------------------------------------------------

run_score_loss_benchmark <- function(sample_sizes,
                                     family = c("univariate", "multivariate"),
                                     method = c("MLE", "KDE", "SM"),
                                     r_sample,
                                     true_score,
                                     n_rep = 50,
                                     n_test = 2000,
                                     smoothed = FALSE,
                                     loss_half = FALSE,
                                     seed = NULL,
                                     fit_args = list(),
                                     predict_args = list(),
                                     verbose = TRUE) {
  family <- match.arg(family)
  method <- match.arg(method)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  sample_sizes <- as.integer(sample_sizes)
  sample_sizes <- sample_sizes[is.finite(sample_sizes) & sample_sizes >= 2L]
  
  if (length(sample_sizes) == 0L) {
    stop("sample_sizes must contain at least one integer >= 2.")
  }
  
  if (family == "multivariate" && method == "SM") {
    stop("Multivariate SM is not implemented yet.")
  }
  
  out_list <- vector("list", length(sample_sizes) * n_rep)
  counter <- 1L
  
  for (n in sample_sizes) {
    if (verbose) {
      message("Running n = ", n, " ...")
    }
    
    for (rep in seq_len(n_rep)) {
      if (verbose) {
        message("  repetition ", rep, "/", n_rep)
      }
      
      ans <- run_one_score_loss_experiment(
        n = n,
        family = family,
        method = method,
        r_sample = r_sample,
        true_score = true_score,
        n_test = n_test,
        smoothed = smoothed,
        loss_half = loss_half,
        fit_args = fit_args,
        predict_args = predict_args
      )
      
      out_list[[counter]] <- data.frame(
        n = ans$n,
        repetition = rep,
        family = ans$family,
        method = ans$method,
        smoothed = ans$smoothed,
        score_loss = ans$score_loss,
        n_test = ans$n_test,
        stringsAsFactors = FALSE
      )
      
      counter <- counter + 1L
    }
  }
  
  results <- do.call(rbind, out_list)
  rownames(results) <- NULL
  
  agg_mean <- aggregate(score_loss ~ n, data = results, FUN = mean)
  names(agg_mean)[2] <- "mean_score_loss"
  
  agg_sd <- aggregate(score_loss ~ n, data = results, FUN = stats::sd)
  names(agg_sd)[2] <- "sd_score_loss"
  
  agg_median <- aggregate(score_loss ~ n, data = results, FUN = stats::median)
  names(agg_median)[2] <- "median_score_loss"
  
  summary_df <- Reduce(
    function(a, b) merge(a, b, by = "n", all = TRUE),
    list(agg_mean, agg_sd, agg_median)
  )
  
  structure(
    list(
      raw = results,
      summary = summary_df,
      settings = list(
        sample_sizes = sample_sizes,
        family = family,
        method = method,
        n_rep = n_rep,
        n_test = n_test,
        smoothed = smoothed,
        loss_half = loss_half,
        fit_args = fit_args,
        predict_args = predict_args
      )
    ),
    class = "score_loss_benchmark"
  )
}

# ------------------------------------------------------------
# (6) Plotting
# ------------------------------------------------------------

plot_score_loss_benchmark <- function(obj,
                                      log_x = TRUE,
                                      log_y = TRUE,
                                      type = "b",
                                      pch = 19,
                                      lwd = 2,
                                      col = "black",
                                      xlab = "sample size n",
                                      ylab = "mean score loss",
                                      main = NULL,
                                      ylim = NULL) {
  if (!inherits(obj, "score_loss_benchmark")) {
    stop("obj must be of class 'score_loss_benchmark'.")
  }
  
  df <- obj$summary
  
  if (is.null(main)) {
    main <- paste0(
      obj$settings$family, " - ",
      obj$settings$method,
      if (isTRUE(obj$settings$method == "MLE")) {
        paste0(" (smoothed = ", obj$settings$smoothed, ")")
      } else {
        ""
      }
    )
  }
  
  x <- df$n
  y <- df$mean_score_loss
  
  if (is.null(ylim)) {
    ylim <- range(y, na.rm = TRUE)
  }
  
  plot(
    x, y,
    log = paste0(if (log_x) "x" else "", if (log_y) "y" else ""),
    type = type,
    pch = pch,
    lwd = lwd,
    col = col,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ylim = ylim
  )
  
  invisible(df)
}

# ------------------------------------------------------------
# (7) Optional multi-panel plot for several benchmark objects
# ------------------------------------------------------------

plot_score_loss_panel <- function(benchmark_list,
                                  names_list = NULL,
                                  nrow = NULL,
                                  ncol = NULL,
                                  log_x = TRUE,
                                  log_y = TRUE) {
  if (!is.list(benchmark_list) || length(benchmark_list) == 0L) {
    stop("benchmark_list must be a non-empty list.")
  }
  
  k <- length(benchmark_list)
  
  if (is.null(names_list)) {
    names_list <- names(benchmark_list)
  }
  if (is.null(names_list) || any(names_list == "")) {
    names_list <- paste("Benchmark", seq_len(k))
  }
  
  if (is.null(nrow) || is.null(ncol)) {
    nrow <- ceiling(sqrt(k))
    ncol <- ceiling(k / nrow)
  }
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(nrow, ncol))
  
  for (i in seq_len(k)) {
    plot_score_loss_benchmark(
      benchmark_list[[i]],
      log_x = log_x,
      log_y = log_y,
      main = names_list[i]
    )
  }
  
  invisible(NULL)
}

# ------------------------------------------------------------
# (8) Example usage
# ------------------------------------------------------------
# Example 1: univariate Gaussian, KDE
#
r_sample_norm_1d <- function(n) rnorm(n, mean = 0, sd = 1)
true_score_norm_1d <- function(x) matrix(x, ncol = 1)   # because -d/dx log f = x

res_kde_1d <- run_score_loss_benchmark(
  sample_sizes = c(50, 100, 200, 500, 1000, 10000),
  family = "univariate",
  method = "KDE",
  r_sample = r_sample_norm_1d,
  true_score = true_score_norm_1d,
  n_rep = 25,
  n_test = 2000,
  seed = 123
)

plot_score_loss_benchmark(res_kde_1d)
#
#
# Example 2: univariate Gaussian, log-concave MLE
#
res_mle_1d <- run_score_loss_benchmark(
  sample_sizes = c(50, 100, 200, 500, 1000),
  family = "univariate",
  method = "MLE",
  r_sample = r_sample_norm_1d,
  true_score = true_score_norm_1d,
  n_rep = 25,
  n_test = 2000,
  smoothed = FALSE,
  seed = 123
)

plot_score_loss_benchmark(res_mle_1d)
#
#
# Example 3: univariate SM
#
res_sm_1d <- run_score_loss_benchmark(
  sample_sizes = c(50, 100, 200, 500, 1000, 10000),
  family = "univariate",
  method = "SM",
  r_sample = r_sample_norm_1d,
  true_score = true_score_norm_1d,
  n_rep = 25,
  n_test = 2000,
  seed = 123,
  fit_args = list(m = 3)
)

plot_score_loss_benchmark(res_sm_1d)
#
#
# Example 4: multivariate Gaussian, KDE
#
r_sample_norm_2d <- function(n) {
  matrix(rnorm(2 * n), ncol = 2)
}

true_score_norm_2d <- function(x) {
  x <- as.matrix(x)
  x   # because -grad log f(x) = x for N(0, I)
}

res_kde_2d <- run_score_loss_benchmark(
  sample_sizes = c(100, 200, 500, 1000),
  family = "multivariate",
  method = "KDE",
  r_sample = r_sample_norm_2d,
  true_score = true_score_norm_2d,
  n_rep = 20,
  n_test = 2000,
  seed = 123
)

plot_score_loss_benchmark(res_kde_2d)

# Example 4: multivariate Gaussian, MLE

res_mle_2d <- run_score_loss_benchmark(
  sample_sizes = c(100, 200, 500, 1000),
  family = "multivariate",
  method = "MLE",
  r_sample = r_sample_norm_2d,
  true_score = true_score_norm_2d,
  n_rep = 20,
  n_test = 2000,
  smoothed = FALSE,
  seed = 123
)

plot_score_loss_benchmark(res_mle_2d)
#
#
# Example 6: compare several methods in one panel
#
# plot_score_loss_panel(
#   list(
#     "KDE 1D" = res_kde_1d,
#     "MLE 1D" = res_mle_1d,
#     "SM 1D"  = res_sm_1d
#   )
# )