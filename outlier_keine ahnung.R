# ============================================================
# outlier_detection_multivariate_compare_sm_scores.R
#
# Multivariate outlier-detection benchmark for:
#   - KDE
#   - multivariate log-concave MLE
#   - multivariate Score Matching with two anomaly scores:
#       (1) squared score norm
#       (2) line-integral score
#
# Main idea:
#   - Fit all methods on INLIER training data only
#   - Evaluate on mixed test data = inliers + outliers
#   - Compare methods mainly by ranking metrics:
#       * ROC-AUC
#       * Average Precision
#       * Best balanced accuracy
#
# SM anomaly score variants:
#   (A) Score norm:
#       S(x) = || r_hat(x) ||^2
#
#   (B) Line integral:
#       S(x) = \int_0^1 r_hat(gamma(t))^T gamma'(t) dt
#       where gamma(t) is the straight line from an anchor point c
#       to x:
#           gamma(t) = c + t (x - c)
#
#       Since r(x) = - grad log f(x), this approximates
#           -log f(x) + constant
#       up to numerical error and path approximation.
#
# ============================================================

# ------------------------------------------------------------
# (0) Load project files
# ------------------------------------------------------------
source("helper_functions.R")
source("KDE.R")
source("LogConcaveMLE.R")
source("Multivariate_Pairwise_Polynomial_Score_Matching.R")

# ------------------------------------------------------------
# (1) Small helpers
# ------------------------------------------------------------

as_matrix_safe <- function(x) {
  x <- as.matrix(x)
  if (!is.numeric(x)) storage.mode(x) <- "double"
  x
}

clean_xy <- function(x, y = NULL) {
  x <- as_matrix_safe(x)
  keep <- apply(x, 1, function(row) all(is.finite(row)))
  x <- x[keep, , drop = FALSE]
  if (is.null(y)) return(list(x = x, keep = keep))
  list(x = x, y = y[keep], keep = keep)
}

trapz_1d_local <- function(x, y) {
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

# ------------------------------------------------------------
# (2) Ranking metrics for anomaly detection
# ------------------------------------------------------------

# labels: 0 = inlier, 1 = outlier
roc_auc_rank <- function(labels, scores) {
  labels <- as.integer(labels)
  scores <- as.numeric(scores)
  keep <- is.finite(scores) & !is.na(labels)
  labels <- labels[keep]
  scores <- scores[keep]
  
  n_pos <- sum(labels == 1L)
  n_neg <- sum(labels == 0L)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  
  r <- rank(scores, ties.method = "average")
  sum_r_pos <- sum(r[labels == 1L])
  (sum_r_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

average_precision <- function(labels, scores) {
  labels <- as.integer(labels)
  scores <- as.numeric(scores)
  keep <- is.finite(scores) & !is.na(labels)
  labels <- labels[keep]
  scores <- scores[keep]
  
  n_pos <- sum(labels == 1L)
  if (n_pos == 0L) return(NA_real_)
  
  ord <- order(scores, decreasing = TRUE)
  y <- labels[ord]
  
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  precision <- tp / (tp + fp)
  
  mean(precision[y == 1L])
}

best_balanced_accuracy <- function(labels, scores) {
  labels <- as.integer(labels)
  scores <- as.numeric(scores)
  keep <- is.finite(scores) & !is.na(labels)
  labels <- labels[keep]
  scores <- scores[keep]
  
  pos <- sum(labels == 1L)
  neg <- sum(labels == 0L)
  if (pos == 0L || neg == 0L) return(NA_real_)
  
  thr <- sort(unique(scores))
  vals <- vapply(thr, function(t) {
    pred <- as.integer(scores >= t)
    tpr <- sum(pred == 1L & labels == 1L) / pos
    tnr <- sum(pred == 0L & labels == 0L) / neg
    0.5 * (tpr + tnr)
  }, numeric(1))
  
  max(vals)
}

# ------------------------------------------------------------
# (3) Data-generating mechanism
# ------------------------------------------------------------

r_inlier_gaussian_2d <- function(n, Sigma) {
  z <- matrix(rnorm(2 * n), ncol = 2)
  z %*% chol(Sigma)
}

r_outlier_shifted_2d <- function(n, mu, Sigma) {
  z <- matrix(rnorm(2 * n), ncol = 2)
  sweep(z %*% chol(Sigma), 2, mu, FUN = "+")
}

# ------------------------------------------------------------
# (4) Method-specific anomaly scores
# ------------------------------------------------------------

# KDE / MLE use negative log-density
anomaly_score_kde_mv <- function(newx, fit, eps = 1e-300) {
  -predict_logdensity_kde_mv(newx, fit, eps = eps)
}

anomaly_score_mle_mv <- function(newx, fit, eps = 1e-300) {
  -predict_logdensity_logconcave_mv(newx, fit, eps = eps)
}

# ------------------------------------------------------------
# (4.1) SM variant A: score norm
# ------------------------------------------------------------

anomaly_score_sm_norm_mv <- function(newx, fit) {
  s_hat <- predict_score_mv_basic(newx, fit)
  rowSums(as.matrix(s_hat)^2)
}

# ------------------------------------------------------------
# (4.2) SM variant B: line-integral score
# ------------------------------------------------------------

# Choose an anchor point c.
# Preferred choice: training center stored in the fitted object.
# Fallback: zero vector.
get_sm_anchor <- function(fit) {
  if (!is.null(fit$scaling) && !is.null(fit$scaling$center)) {
    c0 <- as.numeric(fit$scaling$center)
    if (all(is.finite(c0))) return(c0)
  }
  
  # fallback
  d <- if (!is.null(fit$d)) fit$d else NA_integer_
  if (is.finite(d)) return(rep(0, d))
  
  stop("Could not determine anchor point for line-integral score.")
}

# For one single point x, approximate
#   integral_0^1 r(c + t(x-c))^T (x-c) dt
# via trapezoidal integration on a regular t-grid.
line_integral_score_one <- function(x, fit, anchor = NULL, n_steps = 50L) {
  x <- as.numeric(x)
  
  if (is.null(anchor)) {
    anchor <- get_sm_anchor(fit)
  }
  anchor <- as.numeric(anchor)
  
  if (length(anchor) != length(x)) {
    stop("anchor and x must have the same dimension.")
  }
  if (n_steps < 2L) stop("n_steps must be at least 2.")
  
  direction <- x - anchor
  
  # If x equals anchor, the path length is zero.
  if (sum(direction^2) == 0) return(0)
  
  t_grid <- seq(0, 1, length.out = n_steps)
  
  integrand <- vapply(t_grid, function(t) {
    xt <- matrix(anchor + t * direction, nrow = 1)
    r_t <- as.numeric(predict_score_mv_basic(xt, fit))
    sum(r_t * direction)
  }, numeric(1))
  
  trapz_1d_local(t_grid, integrand)
}

anomaly_score_sm_line_mv <- function(newx,
                                     fit,
                                     anchor = NULL,
                                     n_steps = 50L) {
  newx <- as_matrix_safe(newx)
  
  if (is.null(anchor)) {
    anchor <- get_sm_anchor(fit)
  }
  
  out <- apply(newx, 1, function(row) {
    line_integral_score_one(
      x = row,
      fit = fit,
      anchor = anchor,
      n_steps = n_steps
    )
  })
  
  as.numeric(out)
}

# ------------------------------------------------------------
# (4.3) Unified wrappers
# ------------------------------------------------------------

predict_anomaly_score <- function(newx,
                                  fit,
                                  method,
                                  sm_line_n_steps = 50L,
                                  sm_line_anchor = NULL) {
  method <- match.arg(method, c("KDE", "MLE", "SM_norm", "SM_line"))
  
  if (method == "KDE") {
    return(as.numeric(anomaly_score_kde_mv(newx, fit)))
  }
  if (method == "MLE") {
    return(as.numeric(anomaly_score_mle_mv(newx, fit)))
  }
  if (method == "SM_norm") {
    return(as.numeric(anomaly_score_sm_norm_mv(newx, fit)))
  }
  if (method == "SM_line") {
    return(as.numeric(
      anomaly_score_sm_line_mv(
        newx = newx,
        fit = fit,
        anchor = sm_line_anchor,
        n_steps = sm_line_n_steps
      )
    ))
  }
  
  stop("Unknown method.")
}

# ------------------------------------------------------------
# (5) One single experiment replicate: shared engine
# ------------------------------------------------------------

run_one_outlier_experiment_compare <- function(
    n_train = 200,
    n_test_in = 500,
    n_test_out = 100,
    Sigma_in = matrix(c(1, 0.6, 0.6, 1), 2, 2),
    mu_out = c(3, 3),
    Sigma_out = matrix(c(1, 0.2, 0.2, 1), 2, 2),
    kde_fit_args = list(H_method = "Hpi"),
    mle_fit_args = list(smoothed = FALSE),
    sm_fit_args = list(
      m = 2,
      include_interactions = TRUE,
      standardize = TRUE,
      ridge = 1e-6,
      log_concave = FALSE
    ),
    sm_line_n_steps = 50L,
    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # ----------------------------
  # Step 1: training sample = only inliers
  # ----------------------------
  x_train <- r_inlier_gaussian_2d(n_train, Sigma = Sigma_in)
  
  # ----------------------------
  # Step 2: fit all methods
  # ----------------------------
  fit_kde <- do.call(fit_kde_mv, c(list(x = x_train), kde_fit_args))
  fit_mle <- do.call(fit_logconcave_mle_mv, c(list(x = x_train), mle_fit_args))
  fit_sm  <- do.call(fit_score_matching_mv_basic, c(list(x = x_train), sm_fit_args))
  
  # ----------------------------
  # Step 3: mixed test sample
  # labels: 0 = inlier, 1 = outlier
  # ----------------------------
  x_test_in  <- r_inlier_gaussian_2d(n_test_in, Sigma = Sigma_in)
  x_test_out <- r_outlier_shifted_2d(n_test_out, mu = mu_out, Sigma = Sigma_out)
  
  x_test <- rbind(x_test_in, x_test_out)
  y_test <- c(rep(0L, n_test_in), rep(1L, n_test_out))
  
  tmp <- clean_xy(x_test, y_test)
  x_test <- tmp$x
  y_test <- tmp$y
  
  # ----------------------------
  # Step 4: anomaly scores
  # ----------------------------
  scores_kde <- predict_anomaly_score(x_test, fit_kde, method = "KDE")
  scores_mle <- predict_anomaly_score(x_test, fit_mle, method = "MLE")
  scores_sm_norm <- predict_anomaly_score(x_test, fit_sm, method = "SM_norm")
  
  sm_anchor <- get_sm_anchor(fit_sm)
  scores_sm_line <- predict_anomaly_score(
    x_test,
    fit_sm,
    method = "SM_line",
    sm_line_n_steps = sm_line_n_steps,
    sm_line_anchor = sm_anchor
  )
  
  # ----------------------------
  # Step 5: ranking metrics
  # ----------------------------
  out <- rbind(
    data.frame(
      method = "KDE",
      roc_auc = roc_auc_rank(y_test, scores_kde),
      average_precision = average_precision(y_test, scores_kde),
      best_bal_acc = best_balanced_accuracy(y_test, scores_kde)
    ),
    data.frame(
      method = "MLE",
      roc_auc = roc_auc_rank(y_test, scores_mle),
      average_precision = average_precision(y_test, scores_mle),
      best_bal_acc = best_balanced_accuracy(y_test, scores_mle)
    ),
    data.frame(
      method = "SM_norm",
      roc_auc = roc_auc_rank(y_test, scores_sm_norm),
      average_precision = average_precision(y_test, scores_sm_norm),
      best_bal_acc = best_balanced_accuracy(y_test, scores_sm_norm)
    ),
    data.frame(
      method = "SM_line",
      roc_auc = roc_auc_rank(y_test, scores_sm_line),
      average_precision = average_precision(y_test, scores_sm_line),
      best_bal_acc = best_balanced_accuracy(y_test, scores_sm_line)
    )
  )
  
  list(
    results = out,
    fits = list(KDE = fit_kde, MLE = fit_mle, SM = fit_sm),
    test_data = list(x_test = x_test, y_test = y_test),
    test_scores = list(
      KDE = scores_kde,
      MLE = scores_mle,
      SM_norm = scores_sm_norm,
      SM_line = scores_sm_line
    ),
    sm_line_info = list(
      anchor = sm_anchor,
      n_steps = sm_line_n_steps
    )
  )
}

# ------------------------------------------------------------
# (6) Repeated benchmark: shared engine
# ------------------------------------------------------------

run_outlier_benchmark_compare <- function(
    n_rep = 20,
    n_train = 200,
    n_test_in = 500,
    n_test_out = 100,
    Sigma_in = matrix(c(1, 0.6, 0.6, 1), 2, 2),
    mu_out = c(3, 3),
    Sigma_out = matrix(c(1, 0.2, 0.2, 1), 2, 2),
    kde_fit_args = list(H_method = "Hpi"),
    mle_fit_args = list(smoothed = FALSE),
    sm_fit_args = list(
      m = 2,
      include_interactions = TRUE,
      standardize = TRUE,
      ridge = 1e-6,
      log_concave = FALSE
    ),
    sm_line_n_steps = 50L,
    seed = 123) {
  
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n_rep)
  
  all_results <- vector("list", n_rep)
  last_object <- NULL
  
  for (b in seq_len(n_rep)) {
    ans <- run_one_outlier_experiment_compare(
      n_train = n_train,
      n_test_in = n_test_in,
      n_test_out = n_test_out,
      Sigma_in = Sigma_in,
      mu_out = mu_out,
      Sigma_out = Sigma_out,
      kde_fit_args = kde_fit_args,
      mle_fit_args = mle_fit_args,
      sm_fit_args = sm_fit_args,
      sm_line_n_steps = sm_line_n_steps,
      seed = seeds[b]
    )
    
    tmp <- ans$results
    tmp$rep <- b
    all_results[[b]] <- tmp
    last_object <- ans
  }
  
  raw <- do.call(rbind, all_results)
  
  summary <- aggregate(
    raw[, c("roc_auc", "average_precision", "best_bal_acc")],
    by = list(method = raw$method),
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE))
  )
  
  summary_clean <- data.frame(
    method = summary$method,
    roc_auc_mean = summary$roc_auc[, "mean"],
    roc_auc_sd = summary$roc_auc[, "sd"],
    average_precision_mean = summary$average_precision[, "mean"],
    average_precision_sd = summary$average_precision[, "sd"],
    best_bal_acc_mean = summary$best_bal_acc[, "mean"],
    best_bal_acc_sd = summary$best_bal_acc[, "sd"]
  )
  
  list(
    raw = raw,
    summary = summary_clean,
    last_run = last_object
  )
}

# ------------------------------------------------------------
# (7) Convenience wrappers: same test, separate SM variants
# ------------------------------------------------------------

# Benchmark with original SM score norm only
run_outlier_benchmark_sm_norm <- function(
    n_rep = 20,
    n_train = 200,
    n_test_in = 500,
    n_test_out = 100,
    Sigma_in = matrix(c(1, 0.6, 0.6, 1), 2, 2),
    mu_out = c(3, 3),
    Sigma_out = matrix(c(1, 0.2, 0.2, 1), 2, 2),
    kde_fit_args = list(H_method = "Hpi"),
    mle_fit_args = list(smoothed = FALSE),
    sm_fit_args = list(
      m = 2,
      include_interactions = TRUE,
      standardize = TRUE,
      ridge = 1e-6,
      log_concave = FALSE
    ),
    seed = 123) {
  
  full <- run_outlier_benchmark_compare(
    n_rep = n_rep,
    n_train = n_train,
    n_test_in = n_test_in,
    n_test_out = n_test_out,
    Sigma_in = Sigma_in,
    mu_out = mu_out,
    Sigma_out = Sigma_out,
    kde_fit_args = kde_fit_args,
    mle_fit_args = mle_fit_args,
    sm_fit_args = sm_fit_args,
    sm_line_n_steps = 50L,
    seed = seed
  )
  
  keep_methods <- c("KDE", "MLE", "SM_norm")
  raw_sub <- full$raw[full$raw$method %in% keep_methods, , drop = FALSE]
  summary_sub <- full$summary[full$summary$method %in% keep_methods, , drop = FALSE]
  
  list(
    raw = raw_sub,
    summary = summary_sub,
    last_run = full$last_run
  )
}

# Benchmark with line-integral SM only
run_outlier_benchmark_sm_line <- function(
    n_rep = 20,
    n_train = 200,
    n_test_in = 500,
    n_test_out = 100,
    Sigma_in = matrix(c(1, 0.6, 0.6, 1), 2, 2),
    mu_out = c(3, 3),
    Sigma_out = matrix(c(1, 0.2, 0.2, 1), 2, 2),
    kde_fit_args = list(H_method = "Hpi"),
    mle_fit_args = list(smoothed = FALSE),
    sm_fit_args = list(
      m = 2,
      include_interactions = TRUE,
      standardize = TRUE,
      ridge = 1e-6,
      log_concave = FALSE
    ),
    sm_line_n_steps = 50L,
    seed = 123) {
  
  full <- run_outlier_benchmark_compare(
    n_rep = n_rep,
    n_train = n_train,
    n_test_in = n_test_in,
    n_test_out = n_test_out,
    Sigma_in = Sigma_in,
    mu_out = mu_out,
    Sigma_out = Sigma_out,
    kde_fit_args = kde_fit_args,
    mle_fit_args = mle_fit_args,
    sm_fit_args = sm_fit_args,
    sm_line_n_steps = sm_line_n_steps,
    seed = seed
  )
  
  keep_methods <- c("KDE", "MLE", "SM_line")
  raw_sub <- full$raw[full$raw$method %in% keep_methods, , drop = FALSE]
  summary_sub <- full$summary[full$summary$method %in% keep_methods, , drop = FALSE]
  
  list(
    raw = raw_sub,
    summary = summary_sub,
    last_run = full$last_run
  )
}

# ------------------------------------------------------------
# (8) Plot helpers
# ------------------------------------------------------------

plot_benchmark_boxplots <- function(benchmark_result) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(1, 3))
  boxplot(roc_auc ~ method, data = benchmark_result$raw,
          main = "ROC-AUC", ylab = "value")
  boxplot(average_precision ~ method, data = benchmark_result$raw,
          main = "Average Precision", ylab = "value")
  boxplot(best_bal_acc ~ method, data = benchmark_result$raw,
          main = "Best balanced accuracy", ylab = "value")
}

plot_last_scores_scatter <- function(last_run,
                                     method = c("KDE", "MLE", "SM_norm", "SM_line")) {
  method <- match.arg(method)
  
  x <- last_run$test_data$x_test
  y <- last_run$test_data$y_test
  s <- last_run$test_scores[[method]]
  
  ord <- order(s)
  
  plot(x[ord, 1], x[ord, 2],
       pch = 19,
       cex = 0.8,
       col = ifelse(y[ord] == 1L, "tomato", "grey40"),
       xlab = "x1",
       ylab = "x2",
       main = paste("Test sample colored by label -", method))
  
  s_scaled <- (s - min(s, na.rm = TRUE)) /
    (max(s, na.rm = TRUE) - min(s, na.rm = TRUE) + 1e-12)
  
  points(x[ord, 1], x[ord, 2], pch = 1, cex = 0.7 + 1.3 * s_scaled[ord])
}

# ------------------------------------------------------------
# (9) Example runs
# ------------------------------------------------------------

# --------------------------------
# (9.1) Original benchmark:
# KDE vs MLE vs SM with ||score||^2
# --------------------------------
benchmark_sm_norm <- run_outlier_benchmark_sm_norm(
  n_rep = 20,
  n_train = 200,
  n_test_in = 500,
  n_test_out = 100,
  Sigma_in = matrix(c(1, 0.6, 0.6, 1), 2, 2),
  mu_out = c(3, 3),
  Sigma_out = matrix(c(1, 0.2, 0.2, 1), 2, 2),
  kde_fit_args = list(H_method = "Hpi"),
  mle_fit_args = list(smoothed = FALSE),
  sm_fit_args = list(
    m = 2,
    include_interactions = TRUE,
    standardize = TRUE,
    ridge = 1e-6,
    log_concave = FALSE
  ),
  seed = 123
)

cat("\n==============================\n")
cat("Benchmark with SM score norm\n")
cat("==============================\n")
print(benchmark_sm_norm$summary)

# Optional plots
plot_benchmark_boxplots(benchmark_sm_norm)
plot_last_scores_scatter(benchmark_sm_norm$last_run, method = "KDE")
plot_last_scores_scatter(benchmark_sm_norm$last_run, method = "MLE")
plot_last_scores_scatter(benchmark_sm_norm$last_run, method = "SM_norm")

# --------------------------------
# (9.2) New benchmark:
# KDE vs MLE vs SM with line integral
# --------------------------------
benchmark_sm_line <- run_outlier_benchmark_sm_line(
  n_rep = 20,
  n_train = 200,
  n_test_in = 500,
  n_test_out = 100,
  Sigma_in = matrix(c(1, 0.6, 0.6, 1), 2, 2),
  mu_out = c(3, 3),
  Sigma_out = matrix(c(1, 0.2, 0.2, 1), 2, 2),
  kde_fit_args = list(H_method = "Hpi"),
  mle_fit_args = list(smoothed = FALSE),
  sm_fit_args = list(
    m = 2,
    include_interactions = TRUE,
    standardize = TRUE,
    ridge = 1e-6,
    log_concave = FALSE
  ),
  sm_line_n_steps = 50L,
  seed = 123
)

cat("\n==============================\n")
cat("Benchmark with SM line integral\n")
cat("==============================\n")
print(benchmark_sm_line$summary)

# Optional plots
plot_benchmark_boxplots(benchmark_sm_line)
plot_last_scores_scatter(benchmark_sm_line$last_run, method = "KDE")
plot_last_scores_scatter(benchmark_sm_line$last_run, method = "MLE")
plot_last_scores_scatter(benchmark_sm_line$last_run, method = "SM_line")

# --------------------------------
# (9.3) Direct comparison:
# KDE vs MLE vs both SM variants together
# --------------------------------
benchmark_compare_both <- run_outlier_benchmark_compare(
  n_rep = 20,
  n_train = 200,
  n_test_in = 500,
  n_test_out = 100,
  Sigma_in = matrix(c(1, 0.6, 0.6, 1), 2, 2),
  mu_out = c(3, 3),
  Sigma_out = matrix(c(1, 0.2, 0.2, 1), 2, 2),
  kde_fit_args = list(H_method = "Hpi"),
  mle_fit_args = list(smoothed = FALSE),
  sm_fit_args = list(
    m = 2,
    include_interactions = TRUE,
    standardize = TRUE,
    ridge = 1e-6,
    log_concave = FALSE
  ),
  sm_line_n_steps = 50L,
  seed = 123
)

cat("\n==============================\n")
cat("Benchmark with both SM variants\n")
cat("==============================\n")
print(benchmark_compare_both$summary)

# Optional plots
plot_benchmark_boxplots(benchmark_compare_both)
plot_last_scores_scatter(benchmark_compare_both$last_run, method = "SM_norm")
plot_last_scores_scatter(benchmark_compare_both$last_run, method = "SM_line")

# ------------------------------------------------------------
# (10) Optional extension:
# log-concave SM with both anomaly-score variants
# ------------------------------------------------------------
# Example:
#
# benchmark_compare_both_lc <- run_outlier_benchmark_compare(
#   n_rep = 20,
#   n_train = 200,
#   n_test_in = 500,
#   n_test_out = 100,
#   sm_fit_args = list(
#     m = 2,
#     include_interactions = TRUE,
#     standardize = TRUE,
#     ridge = 1e-6,
#     log_concave = TRUE,
#     lc_method = "grid",
#     lc_grid_size = 5L,
#     lc_penalty = 1e4
#   ),
#   sm_line_n_steps = 50L,
#   seed = 123
# )
#
# print(benchmark_compare_both_lc$summary)
# plot_benchmark_boxplots(benchmark_compare_both_lc)
# ============================================================