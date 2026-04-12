# ============================================================
# Final_Univariate_Test_Template_patched.R
# Template for univariate benchmark runs under the simplified
# evaluation design.
#
# Core metrics:
#   - kl
#   - score_loss
# Optional variants from metric args:
#   - *_central      : metric on empirical bulk region only
#   - score_loss_trim: score metric after trimming the largest
#                      pointwise score losses
# ============================================================

source("helper_functions.R")
source("KDE.R")
source("LogConcaveMLE.R")
source("Univariate_Polynomial_Score_Matching_1.0.R")
source("Multivariate_Pairwise_Polynomial_Score_Matching.R")
source("Evaluation_Metrics_Clean.R")
source("BiasVariance_Score_Clean_patched.R")
source("Tests_Clean_patched.R")

# ------------------------------------------------------------
# A) Global settings
# ------------------------------------------------------------
sample_sizes_main <- c(50, 100, 200, 500, 1000)
sample_sizes_bias_variance <- c(50, 100, 200, 500, 1000)
metrics_main <- c("kl", "score_loss")
seed_main <- 123

# central_trim evaluates the metric on the empirical bulk of the test set.
# robust_trim only applies to score_loss and trims the largest pointwise
# score losses after the score errors have been computed.
score_metric_args_default <- list(
  central_trim = NULL,
  robust_trim = NULL
)

density_metric_args_default <- list(
  central_trim = NULL
)

# ------------------------------------------------------------
# B) True univariate distributions
# ------------------------------------------------------------
make_truth_gaussian <- function(mean = 0, sd = 1) {
  list(
    name = "gaussian",
    family = "univariate",
    r_sample = function(n) stats::rnorm(n, mean = mean, sd = sd),
    true_density = function(x) stats::dnorm(x, mean = mean, sd = sd),
    true_logdensity = function(x) stats::dnorm(x, mean = mean, sd = sd, log = TRUE),
    true_score = function(x) matrix((as.numeric(x) - mean) / (sd^2), ncol = 1)
  )
}

make_truth_logistic <- function(location = 0, scale = 1) {
  list(
    name = "logistic",
    family = "univariate",
    r_sample = function(n) stats::rlogis(n, location = location, scale = scale),
    true_density = function(x) stats::dlogis(x, location = location, scale = scale),
    true_logdensity = function(x) stats::dlogis(x, location = location, scale = scale, log = TRUE),
    true_score = function(x) {
      z <- (as.numeric(x) - location) / scale
      matrix(tanh(z / 2) / scale, ncol = 1)
    }
  )
}

make_truth_gumbel <- function(location = 0, scale = 1) {
  rgumbel <- function(n) location - scale * log(-log(stats::runif(n)))
  dgumbel <- function(x, log = FALSE) {
    z <- (as.numeric(x) - location) / scale
    ld <- -log(scale) - z - exp(-z)
    if (log) ld else exp(ld)
  }
  list(
    name = "gumbel",
    family = "univariate",
    r_sample = rgumbel,
    true_density = function(x) dgumbel(x, log = FALSE),
    true_logdensity = function(x) dgumbel(x, log = TRUE),
    true_score = function(x) {
      z <- (as.numeric(x) - location) / scale
      matrix((1 - exp(-z)) / scale, ncol = 1)
    }
  )
}

make_truth_student <- function(df = 3) {
  list(
    name = sprintf("student_t_df%d", df),
    family = "univariate",
    r_sample = function(n) stats::rt(n, df = df),
    true_density = function(x) stats::dt(x, df = df),
    true_logdensity = function(x) stats::dt(x, df = df, log = TRUE),
    true_score = function(x) {
      xx <- as.numeric(x)
      matrix(((df + 1) * xx) / (df + xx^2), ncol = 1)
    }
  )
}

make_truth_laplace <- function(location = 0, scale = 1) {
  rlaplace <- function(n) {
    u <- stats::runif(n, min = -0.5, max = 0.5)
    location - scale * sign(u) * log(1 - 2 * abs(u))
  }
  dlaplace <- function(x, log = FALSE) {
    z <- abs(as.numeric(x) - location) / scale
    ld <- -log(2 * scale) - z
    if (log) ld else exp(ld)
  }
  list(
    name = "laplace",
    family = "univariate",
    r_sample = rlaplace,
    true_density = function(x) dlaplace(x, log = FALSE),
    true_logdensity = function(x) dlaplace(x, log = TRUE),
    true_score = function(x) {
      xx <- as.numeric(x)
      out <- sign(xx - location) / scale
      out[abs(xx - location) < .Machine$double.eps^0.5] <- 0
      matrix(out, ncol = 1)
    }
  )
}

truth_gaussian <- make_truth_gaussian(-50, 4)
truth_logistic <- make_truth_logistic()
truth_gumbel <- make_truth_gumbel()
truth_student <- make_truth_student(df = 3)
truth_laplace <- make_truth_laplace()

all_truths <- list(
  gaussian = truth_gaussian,
  logistic = truth_logistic,
  gumbel = truth_gumbel,
  student = truth_student,
  laplace = truth_laplace
)

# ------------------------------------------------------------
# C) Candidate lists per method family
# ------------------------------------------------------------
make_kde_specs_1d <- function(score_metric_args = score_metric_args_default,
                              density_metric_args = density_metric_args_default) {
  lapply(c("SJ", "nrd0", "ucv", "bcv"), function(bw) {
    list(
      label = paste0("KDE_", bw),
      method = "KDE",
      smoothed = FALSE,
      fit_args = list(bw = bw),
      density_predict_args = list(),
      score_predict_args = list(),
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args
    )
  })
}

make_mle_specs_1d <- function(score_metric_args = score_metric_args_default,
                              density_metric_args = density_metric_args_default) {
  list(
    list(
      label = "MLE_unsmoothed",
      method = "MLE",
      smoothed = FALSE,
      fit_args = list(),
      density_predict_args = list(),
      score_predict_args = list(),
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args
    ),
    list(
      label = "MLE_smoothed",
      method = "MLE",
      smoothed = TRUE,
      fit_args = list(),
      density_predict_args = list(),
      score_predict_args = list(),
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args
    )
  )
}

make_sm_specs_1d <- function(m_values = c(1, 2, 3, 4, 5),
                             score_metric_args = score_metric_args_default,
                             density_metric_args = density_metric_args_default) {
  lapply(m_values, function(m) {
    list(
      label = paste0("SM_m", m),
      method = "SM",
      smoothed = FALSE,
      fit_args = list(m = m, standardize = TRUE, ridge = 1e-2),
      density_predict_args = list(subdivisions = 200L, rel.tol = 1e-8, stop_on_failure = FALSE),
      score_predict_args = list(),
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args
    )
  })
}

# ------------------------------------------------------------
# D) Benchmark wrappers
# ------------------------------------------------------------
run_family_selection_benchmark <- function(truth,
                                           estimator_specs,
                                           sample_sizes = sample_sizes_main,
                                           metrics = metrics_main,
                                           n_rep = 50,
                                           n_test = 2000,
                                           seed = seed_main,
                                           verbose = TRUE,
                                           save = FALSE,
                                           save_dir = ".",
                                           save_name = NULL) {
  run_final_benchmark(
    sample_sizes = sample_sizes,
    family = truth$family,
    estimator_specs = estimator_specs,
    r_sample = truth$r_sample,
    metrics = metrics,
    n_rep = n_rep,
    n_test = n_test,
    true_density = truth$true_density,
    true_logdensity = truth$true_logdensity,
    true_score = truth$true_score,
    truth_name = truth$name,
    seed = seed,
    verbose = verbose,
    save = save,
    save_dir = save_dir,
    save_name = save_name
  )
}

run_family_bias_variance <- function(truth,
                                     estimator_specs,
                                     sample_sizes = sample_sizes_bias_variance,
                                     n_rep = 50,
                                     seed = seed_main,
                                     verbose = TRUE,
                                     save = FALSE,
                                     save_dir = ".",
                                     save_name = NULL) {
  eval_grid <- make_eval_grid(r_sample = truth$r_sample, family = truth$family, seed = seed)
  run_bias_variance_score_benchmark(
    sample_sizes = sample_sizes,
    family = truth$family,
    estimator_specs = estimator_specs,
    r_sample = truth$r_sample,
    true_score = truth$true_score,
    eval_grid = eval_grid,
    truth_name = truth$name,
    n_rep = n_rep,
    seed = seed,
    verbose = verbose,
    save = save,
    save_dir = save_dir,
    save_name = save_name
  )
}

# ------------------------------------------------------------
# E) Example runs
# ------------------------------------------------------------
# Example: compare SM candidates on Gaussian data and save the object.
# sm_gaussian_candidates <- run_family_selection_benchmark(
#   truth = truth_gaussian,
#   estimator_specs = make_sm_specs_1d(
#     m_values = c(1,2,3,4, 5),
#     score_metric_args = list(central_trim = 0.05, robust_trim = 0.01),
#     density_metric_args = list(central_trim = 0.05)
#   ),
#   n_rep = 10,
#   n_test = 1000,
#   save = TRUE,
#   save_dir = "results"
# )
# #
# aggregate_final_benchmark(sm_gaussian_candidates, metric = "kl", across_runs_center = "median")
# aggregate_final_benchmark(sm_gaussian_candidates, metric = "score_loss_trim", across_runs_center = "median")
# plot_final_benchmark(sm_gaussian_candidates, metric = "kl", center = "mean", interval = "none", log_y = TRUE)
# plot_final_benchmark(sm_gaussian_candidates, metric = "score_loss", center = "median", interval = "none", log_y = TRUE)
#
# Bias-variance example:
sm_gaussian_bv <- run_family_bias_variance(
  truth = truth_gaussian,
  estimator_specs = make_sm_specs_1d(m_values = c(1, 2, 3, 4, 5)),
  n_rep = 75,
  save = TRUE,
  save_dir = "results"
)
plot_score_bias_variance(sm_gaussian_bv, metric = "integrated_variance", log_y = TRUE)

# sm_logistic_candidates <- run_family_selection_benchmark(
#   truth = truth_logistic,
#   estimator_specs = make_sm_specs_1d(
#     m_values = c(1, 2, 3, 4, 5, 6),
#     score_metric_args = list(central_trim = 0.05, robust_trim = 0.01),
#     density_metric_args = list(central_trim = 0.05)
#   ),
#   n_rep = 10,
#   n_test = 1000,
#   save = TRUE,
#   save_dir = "results"
# )
# 
# aggregate_final_benchmark(sm_logistic_candidates, metric = "kl", across_runs_center = "median")
# aggregate_final_benchmark(sm_logistic_candidates, metric = "score_loss_trim", across_runs_center = "median")
# plot_final_benchmark(sm_logistic_candidates, metric = "kl", center = "mean", interval = "none", log_y = TRUE)
# plot_final_benchmark(sm_logistic_candidates, metric = "score_loss", center = "median", interval = "none", log_y = TRUE)

# Bias-variance example:
sm_logistic_bv <- run_family_bias_variance(
  truth = truth_logistic,
  estimator_specs = make_sm_specs_1d(m_values = c(1, 2, 3, 4, 5, 6)),
  n_rep = 75,
  save = TRUE,
  save_dir = "results"
)
plot_score_bias_variance(sm_logistic_bv, metric = "integrated_variance", log_y = TRUE)

