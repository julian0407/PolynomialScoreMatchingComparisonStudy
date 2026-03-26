# ============================================================
# Simple multivariate outlier-detection benchmark
# for KDE, multivariate log-concave MLE and multivariate
# Score Matching, based on the existing project files.
#
# Main idea:
#   - Fit each method only on INLIER training data.
#   - Evaluate on a test set containing inliers + outliers.
#   - For KDE and MLE, use NEGATIVE LOG-DENSITY as anomaly score.
#       Reason: lower estimated density means "more atypical".
#   - For multivariate Score Matching, a normalized density is not
#     available in the current code base, so we use the SQUARED
#     NORM OF THE ESTIMATED SCORE FIELD as anomaly score.
#       Reason: for unimodal log-concave distributions, the score
#       norm is typically small near the mode and larger in the tails,
#       so it can act as a simple outlier score.
#
# Why this test is useful:
#   - It is simple and application-driven.
#   - It avoids forcing density metrics on multivariate SM.
#   - ROC-AUC and Average Precision only require a ranking of points,
#     not scores on the same numerical scale.
#
# IMPORTANT:
#   The anomaly scores are method-specific and live on different scales.
#   Therefore, we compare methods mainly via ranking metrics:
#     * ROC-AUC
#     * Average Precision (AP)
#   and NOT via raw score magnitudes.
# ============================================================

# ------------------------------------------------------------
# (0) Load project files
# ------------------------------------------------------------
# Adjust the paths if needed.
source("helper_functions.R")
source("KDE.R")
source("LogConcaveMLE.R")
source("Multivariate_Pairwise_Polynomial_Score_Matching.R")

# ------------------------------------------------------------
# (1) Small helpers
# ------------------------------------------------------------

# Ensure input is stored as observation matrix.
as_matrix_safe <- function(x) {
  x <- as.matrix(x)
  if (!is.numeric(x)) storage.mode(x) <- "double"
  x
}

# Remove rows with non-finite values and keep labels aligned.
clean_xy <- function(x, y = NULL) {
  x <- as_matrix_safe(x)
  keep <- apply(x, 1, function(row) all(is.finite(row)))
  x <- x[keep, , drop = FALSE]
  if (is.null(y)) return(list(x = x, keep = keep))
  list(x = x, y = y[keep], keep = keep)
}

# ------------------------------------------------------------
# (2) Metrics for anomaly detection
# ------------------------------------------------------------

# ROC-AUC via rank formula.
# labels must be 0/1 with 1 = outlier, 0 = inlier.
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

# Average Precision (area under precision-recall curve in step form).
# This is often informative when outliers are rare.
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

  # AP = mean precision at ranks where a positive appears.
  mean(precision[y == 1L])
}

# Optional threshold-based metric:
# Best balanced accuracy over all score thresholds.
# This is NOT the main metric, but can help interpretation.
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
# We use a very simple and plausible anomaly-detection setup:
#
#   Inliers:  X ~ N(0, Sigma)
#   Outliers: X ~ N(mu_out, Sigma_out)
#
# Why this is a good first test:
#   - simple and transparent
#   - genuinely multivariate
#   - not artificially tailored to one method
#   - a shifted Gaussian contamination model is a standard anomaly setup
#
# We choose a 2D correlated Gaussian so that the task is not purely trivial,
# but still easy to visualize and debug.

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
# Convention used here:
#   larger score  => more anomalous
#
# KDE / MLE:
#   anomaly_score(x) = - log f_hat(x)
#   because low density means atypical observation.
#
# SM:
#   anomaly_score(x) = || r_hat(x) ||^2
#   where r_hat is the estimated score in your code convention
#   r(x) = - grad log f(x).
#   This is only a proxy, but it is simple, available, and natural for SM.

anomaly_score_kde_mv <- function(newx, fit, eps = 1e-300) {
  -predict_logdensity_kde_mv(newx, fit, eps = eps)
}

anomaly_score_mle_mv <- function(newx, fit, eps = 1e-300) {
  -predict_logdensity_logconcave_mv(newx, fit, eps = eps)
}

anomaly_score_sm_mv <- function(newx, fit) {
  s_hat <- predict_score_mv_basic(newx, fit)
  rowSums(as.matrix(s_hat)^2)
}

# Unified wrapper so the benchmark loop stays simple.
predict_anomaly_score <- function(newx, fit, method) {
  method <- match.arg(method, c("KDE", "MLE", "SM"))

  if (method == "KDE") return(as.numeric(anomaly_score_kde_mv(newx, fit)))
  if (method == "MLE") return(as.numeric(anomaly_score_mle_mv(newx, fit)))
  if (method == "SM")  return(as.numeric(anomaly_score_sm_mv(newx, fit)))

  stop("Unknown method.")
}

# ------------------------------------------------------------
# (5) One single experiment replicate
# ------------------------------------------------------------
# This function:
#   1) generates training inliers
#   2) fits all three methods
#   3) generates test inliers + test outliers
#   4) computes anomaly scores
#   5) evaluates ROC-AUC / AP / best balanced accuracy
#
# Returning one row per method makes later aggregation very easy.

run_one_outlier_experiment <- function(
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
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # ----------------------------
  # Step 1: training sample = only nominal observations
  # ----------------------------
  x_train <- r_inlier_gaussian_2d(n_train, Sigma = Sigma_in)

  # ----------------------------
  # Step 2: fit all methods on the same training sample
  # ----------------------------
  fit_kde <- do.call(fit_kde_mv, c(list(x = x_train), kde_fit_args))
  fit_mle <- do.call(fit_logconcave_mle_mv, c(list(x = x_train), mle_fit_args))
  fit_sm  <- do.call(fit_score_matching_mv_basic, c(list(x = x_train), sm_fit_args))

  # ----------------------------
  # Step 3: build mixed test sample
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
  # Step 4: compute method-specific anomaly scores
  # ----------------------------
  scores_kde <- predict_anomaly_score(x_test, fit_kde, method = "KDE")
  scores_mle <- predict_anomaly_score(x_test, fit_mle, method = "MLE")
  scores_sm  <- predict_anomaly_score(x_test, fit_sm,  method = "SM")

  # ----------------------------
  # Step 5: ranking-based evaluation metrics
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
      method = "SM",
      roc_auc = roc_auc_rank(y_test, scores_sm),
      average_precision = average_precision(y_test, scores_sm),
      best_bal_acc = best_balanced_accuracy(y_test, scores_sm)
    )
  )

  # We also return the last fitted objects and the last test scores,
  # because this is useful for debugging and plotting.
  list(
    results = out,
    fits = list(KDE = fit_kde, MLE = fit_mle, SM = fit_sm),
    test_data = list(x_test = x_test, y_test = y_test),
    test_scores = list(KDE = scores_kde, MLE = scores_mle, SM = scores_sm)
  )
}

# ------------------------------------------------------------
# (6) Repeated benchmark
# ------------------------------------------------------------
# A single repetition can be noisy.
# Therefore, we repeat the experiment several times and average the
# resulting metrics. This is the same logic as in your simulation tests.

run_outlier_benchmark <- function(
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

  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n_rep)

  all_results <- vector("list", n_rep)
  last_object <- NULL

  for (b in seq_len(n_rep)) {
    ans <- run_one_outlier_experiment(
      n_train = n_train,
      n_test_in = n_test_in,
      n_test_out = n_test_out,
      Sigma_in = Sigma_in,
      mu_out = mu_out,
      Sigma_out = Sigma_out,
      kde_fit_args = kde_fit_args,
      mle_fit_args = mle_fit_args,
      sm_fit_args = sm_fit_args,
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

  # Convert aggregate output into a cleaner data frame.
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
# (7) Simple plotting helpers
# ------------------------------------------------------------
# These are intentionally basic, so the script stays lightweight.

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

plot_last_scores_scatter <- function(last_run, method = c("KDE", "MLE", "SM")) {
  method <- match.arg(method)

  x <- last_run$test_data$x_test
  y <- last_run$test_data$y_test
  s <- last_run$test_scores[[method]]

  # Order points so that low-score points are plotted first and
  # high-score points remain visible on top.
  ord <- order(s)

  plot(x[ord, 1], x[ord, 2],
       pch = 19,
       cex = 0.8,
       col = ifelse(y[ord] == 1L, "tomato", "grey40"),
       xlab = "x1",
       ylab = "x2",
       main = paste("Test sample colored by label -", method))

  # Add a second layer that highlights very anomalous points by size.
  # This is only for quick visual debugging.
  s_scaled <- (s - min(s, na.rm = TRUE)) / (max(s, na.rm = TRUE) - min(s, na.rm = TRUE) + 1e-12)
  points(x[ord, 1], x[ord, 2], pch = 1, cex = 0.7 + 1.3 * s_scaled[ord])
}

# ------------------------------------------------------------
# (8) Example run
# ------------------------------------------------------------
# You can run this block directly as a smoke test.
# Expected behavior in this simple setup:
#   - all methods should usually achieve decent ROC-AUC
#   - KDE / MLE use -log f_hat(x)
#   - SM uses ||score_hat(x)||^2
#
# If you want a slightly harder problem, reduce mu_out from c(3,3)
# to c(2,2). If you want an easier problem, increase it.

benchmark_res <- run_outlier_benchmark(
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

# Print aggregated results.
print(benchmark_res$summary)

# Optional quick plots.
plot_benchmark_boxplots(benchmark_res)
plot_last_scores_scatter(benchmark_res$last_run, method = "KDE")
plot_last_scores_scatter(benchmark_res$last_run, method = "MLE")
plot_last_scores_scatter(benchmark_res$last_run, method = "SM")

# ------------------------------------------------------------
# (9) Optional extension: log-concave score matching variant
# ------------------------------------------------------------
# If you also want to compare your penalized/log-concave SM variant,
# you can simply run a second benchmark with log_concave = TRUE.
# Example:
#
# benchmark_res_sm_lc <- run_outlier_benchmark(
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
#   )
# )
#
# This keeps the test design fixed and only changes the SM fit.
# That is a good way to study whether the log-concavity constraint
# improves anomaly ranking in this simple application.
# ============================================================
