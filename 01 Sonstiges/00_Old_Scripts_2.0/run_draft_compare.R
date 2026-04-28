# ============================================================
# run_draft_compare.R
# Führt Vergleichsläufe mit den Draft-Skripten aus
# und speichert die Ergebnisse als RDS
# ============================================================

rm(list = ls())

source("helper_functions.R")
source("KDE.R")
source("LogConcaveMLE.R")
source("Univariate_Polynomial_Score_Matching_1.0.R")
source("Multivariate_Pairwise_Polynomial_Score_Matching.R")
source("Draft_Evaluation_Metrics.R")
# source("BiasVariance_Score.R")
source("Draft_Tests.R")   # optional, falls du dort bestehende Helfer nutzen willst

compare_seed <- 123
compare_sample_sizes <- c(50, 100, 200, 500, 1000)
compare_metrics <- c("kl", "score_loss", "negloglik")
compare_n_rep <- 10
compare_n_test <- 1000

compare_score_metric_args_draft <- list(
  robust = "none",
  trim_alpha = 0.05
)

compare_density_metric_args <- list(
  robust = "none",
  trim_alpha = 0.05,
  outlier_dom_threshold = 0.25
)

make_truth_gaussian_compare <- function() {
  list(
    name = "gaussian",
    family = "univariate",
    r_sample = function(n) stats::rnorm(n, mean = -50, sd = 4),
    true_density = function(x) stats::dnorm(x, mean = -50, sd = 4),
    true_logdensity = function(x) stats::dnorm(x, mean = -50, sd = 4, log = TRUE),
    true_score = function(x) matrix((as.numeric(x) + 50) / 16, ncol = 1)
  )
}

make_truth_logistic_compare <- function() {
  list(
    name = "logistic",
    family = "univariate",
    r_sample = function(n) stats::rlogis(n, location = 0, scale = 1),
    true_density = function(x) stats::dlogis(x, location = 0, scale = 1),
    true_logdensity = function(x) stats::dlogis(x, location = 0, scale = 1, log = TRUE),
    true_score = function(x) {
      z <- as.numeric(x)
      matrix(tanh(z / 2), ncol = 1)
    }
  )
}

make_sm_specs <- function(m_values) {
  lapply(m_values, function(m) {
    list(
      label = paste0("SM_m", m),
      method = "SM",
      smoothed = FALSE,
      fit_args = list(
        m = m,
        standardize = TRUE,
        ridge = 1e-2
      ),
      density_predict_args = list(
        subdivisions = 200L,
        rel.tol = 1e-8,
        stop_on_failure = FALSE
      ),
      score_predict_args = list(),
      score_metric_args = compare_score_metric_args,
      density_metric_args = compare_density_metric_args
    )
  })
}

run_case <- function(truth, estimator_specs) {
  run_final_benchmark(
    sample_sizes = compare_sample_sizes,
    family = truth$family,
    estimator_specs = estimator_specs,
    r_sample = truth$r_sample,
    metrics = compare_metrics,
    n_rep = compare_n_rep,
    n_test = compare_n_test,
    true_density = truth$true_density,
    true_logdensity = truth$true_logdensity,
    true_score = truth$true_score,
    seed = compare_seed,
    verbose = FALSE
  )
}

draft_results <- list(
  gaussian = run_case(make_truth_gaussian_compare(), make_sm_specs(1:5)),
  logistic = run_case(make_truth_logistic_compare(), make_sm_specs(1:6))
)

saveRDS(draft_results, file = "draft_compare_results_seed123.rds")
cat("Draft results saved to draft_compare_results_seed123.rds\n")