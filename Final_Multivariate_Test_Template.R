# ============================================================
# Final_Multivariate_Benchmark_Runs_simple.R
# Simple final multivariate benchmark template
#
# Design:
#   For each Gaussian multivariate truth object, call run_final_benchmark once
#   and save one benchmark object for that density configuration.
#
# Truths:
#   - Gaussian with independent components, d = 2, 3, 4
#   - Gaussian with dependent components,   d = 2, 3, 4
#
# Estimators per density configuration, exactly 11 total:
#   - 10 SM estimators, one for each m = 1, 2, 3, 4, 5.
#     Each SM estimator is pairwise polynomial SM with ridge and
#     (no) grid log-concavity penalty.
#   - 1 KDE estimators: Hns.

# ============================================================

source("01_Rscripts/03_Test_Framework/Unified_Testing_Framework.R")

# ------------------------------------------------------------
# (1) Global settings
# ------------------------------------------------------------

sample_sizes_mv_main <- c(50, 100, 200, 500, 1000, 5000)
metrics_mv_main <- c("score_loss")
seed_mv_main <- 123

score_metric_args_mv_default <- list(
  central_trim = NULL,
  robust_trim  = NULL
)

density_metric_args_mv_default <- list(
  central_trim = NULL
)

# ------------------------------------------------------------
# (2) Base-R multivariate Gaussian helpers
# ------------------------------------------------------------

rmvnorm_base <- function(n, mean, Sigma) {
  mean <- as.numeric(mean)
  Sigma <- as.matrix(Sigma)
  d <- length(mean)
  if (!all(dim(Sigma) == c(d, d))) stop("Sigma has incompatible dimension.")

  R <- chol(Sigma)
  z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  sweep(z %*% R, 2, mean, FUN = "+")
}

dmvnorm_base <- function(x, mean, Sigma, log = FALSE) {
  x <- as.matrix(x)
  mean <- as.numeric(mean)
  Sigma <- as.matrix(Sigma)
  d <- length(mean)
  if (ncol(x) != d) stop("x has incompatible dimension.")

  R <- chol(Sigma)
  xc <- sweep(x, 2, mean, FUN = "-")
  sol <- backsolve(R, t(xc), transpose = TRUE)
  quad <- colSums(sol^2)
  logdet <- 2 * sum(log(diag(R)))
  ld <- -0.5 * (d * log(2 * pi) + logdet + quad)

  if (isTRUE(log)) ld else exp(ld)
}

# In this code base the score is the negative gradient of log-density.
# For N(mu, Sigma), this is Sigma^{-1}(x - mu).
make_truth_gaussian_mv <- function(d = 2,
                                   scenario = c("independent", "dependent"),
                                   mean = rep(0, d),
                                   sd = NULL,
                                   rho = 0.6) {
  scenario <- match.arg(scenario)
  mean <- as.numeric(mean)
  if (length(mean) != d) stop("mean must have length d.")

  if (is.null(sd)) sd <- seq(1, 1 + 0.25 * (d - 1), length.out = d)
  sd <- as.numeric(sd)
  if (length(sd) != d || any(sd <= 0)) stop("sd must be positive and have length d.")

  if (scenario == "independent") {
    Sigma <- diag(sd^2, nrow = d, ncol = d)
  } else {
    Corr <- outer(seq_len(d), seq_len(d), function(i, j) rho^abs(i - j))
    Sigma <- diag(sd, d) %*% Corr %*% diag(sd, d)
  }

  Sigma_inv <- solve(Sigma)

  list(
    name = sprintf("mv_gaussian_%s_d%d", scenario, d),
    family = "multivariate",
    d = d,
    scenario = scenario,
    mean = mean,
    Sigma = Sigma,
    Sigma_inv = Sigma_inv,
    r_sample = function(n) rmvnorm_base(n = n, mean = mean, Sigma = Sigma),
    true_density = function(x) dmvnorm_base(x, mean = mean, Sigma = Sigma, log = FALSE),
    true_logdensity = function(x) dmvnorm_base(x, mean = mean, Sigma = Sigma, log = TRUE),
    true_score = function(x) {
      x <- as.matrix(x)
      if (ncol(x) != d) stop("x has incompatible dimension.")
      sweep(x, 2, mean, FUN = "-") %*% Sigma_inv
    }
  )
}

# ------------------------------------------------------------
# (3) Truth objects
# ------------------------------------------------------------

make_mv_gaussian_truths <- function(d_values = c(2, 3, 4), rho = 0.6) {
  truths <- list()
  for (d in d_values) {
    truths[[paste0("independent_d", d)]] <- make_truth_gaussian_mv(
      d = d,
      scenario = "independent",
      rho = rho
    )
    truths[[paste0("dependent_d", d)]] <- make_truth_gaussian_mv(
      d = d,
      scenario = "dependent",
      rho = rho
    )
  }
  truths
}

all_mv_truths <- make_mv_gaussian_truths()

# ------------------------------------------------------------
# (4) Estimator specs
# ------------------------------------------------------------

make_kde_specs_mv <- function(H_methods = c("Hpi", "Hns"),
                              diagonal = FALSE,
                              score_metric_args = score_metric_args_mv_default,
                              density_metric_args = density_metric_args_mv_default) {
  lapply(H_methods, function(H_method) {
    list(
      label = paste0("KDE_", H_method, if (isTRUE(diagonal)) "_diag" else "_full"),
      method = "KDE",
      smoothed = FALSE,
      fit_args = list(
        H_method = H_method,
        diagonal = diagonal
      ),
      density_predict_args = list(),
      score_predict_args = list(),
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args
    )
  })
}

make_mle_specs_mv <- function(score_metric_args = score_metric_args_mv_default,
                              density_metric_args = density_metric_args_mv_default) {
  list(
    list(
      label = "MLE_smoothed",
      method = "MLE",
      smoothed = TRUE,
      fit_args = list(),
      density_predict_args = list(),
      score_predict_args = list(h = 1e-4),
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args
    )
  )
}

make_sm_specs_mv_lc_compare <- function(m_values = c(1, 2, 3, 4, 5),
                                        ridge = 1e-2,
                                        lc_grid_size = 5L,
                                        lc_penalty = 1e2,
                                        optim_grad = TRUE,
                                        score_metric_args = score_metric_args_mv_trimmed,
                                        density_metric_args = density_metric_args_mv_trimmed) {
  specs <- list()
  
  for (m in m_values) {
    specs <- c(specs, list(
      list(
        label = paste0("SM_no_logconcave_m", m),
        method = "SM",
        smoothed = FALSE,
        fit_args = list(
          m = m,
          include_interactions = TRUE,
          standardize = TRUE,
          ridge = ridge,
          log_concave = FALSE,
          optim_grad = TRUE
        ),
        density_predict_args = list(),
        score_predict_args = list(),
        score_metric_args = score_metric_args,
        density_metric_args = density_metric_args
      ),
      list(
        label = paste0("SM_grid_logconcave_m", m),
        method = "SM",
        smoothed = FALSE,
        fit_args = list(
          m = m,
          include_interactions = TRUE,
          standardize = TRUE,
          ridge = ridge,
          log_concave = TRUE,
          lc_method = "grid",
          lc_grid_size = lc_grid_size,
          lc_penalty = lc_penalty,
          optim_grad = TRUE
        ),
        density_predict_args = list(),
        score_predict_args = list(),
        score_metric_args = score_metric_args,
        density_metric_args = density_metric_args
      )
    ))
  }
  
  specs
}


make_estimator_specs_for_mv_truth <- function(kde_H_methods = c("Hns"),
                                              kde_diagonal = FALSE,
                                              sm_m_values = c(1, 2, 3, 4, 5),
                                              sm_ridge = 1e-2,
                                              sm_lc_grid_size = 5L,
                                              sm_lc_penalty = 1e2,
                                              optim_grad = TRUE,
                                              score_metric_args = score_metric_args_mv_default,
                                              density_metric_args = density_metric_args_mv_default) {
  specs <- list()

  specs <- c(
    make_sm_specs_mv_lc_compare(
      m_values = sm_m_values,
      ridge = sm_ridge,
      lc_grid_size = sm_lc_grid_size,
      lc_penalty = sm_lc_penalty,
      optim_grad = optim_grad,
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args
    ),
    make_kde_specs_mv(
      H_methods = kde_H_methods,
      diagonal = kde_diagonal,
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args
    )
  )

  specs
}

# ------------------------------------------------------------
# (5) Benchmark wrapper
# ------------------------------------------------------------

run_mv_family_selection_benchmark <- function(truth,
                                              estimator_specs,
                                              sample_sizes = sample_sizes_mv_main,
                                              metrics = metrics_mv_main,
                                              n_rep = 20,
                                              n_test = 3000,
                                              seed = seed_mv_main,
                                              verbose = TRUE,
                                              save = FALSE,
                                              save_dir = "resultsMulti",
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

# ------------------------------------------------------------
# (6) Final benchmark settings
# ------------------------------------------------------------
# Pool of exactly the 11 requested multivariate estimators per density:
#   10 SM = pairwise + ridge + (no) grid log-concavity for m = 1, 2, 3, 4, 5
#   1 KDE = Hns

score_metric_args_mv_trimmed <- list(
  central_trim = 0.05,
  robust_trim  = 0.01
)

density_metric_args_mv_trimmed <- list(
  central_trim = 0.05
)

estimator_specs <- make_estimator_specs_for_mv_truth(
                                              kde_H_methods = c("Hns"),
                                              sm_m_values = c(1, 5),
                                              sm_ridge = 1e-2,
                                              sm_lc_grid_size = 10L,
                                              sm_lc_penalty = 1e2,
                                              score_metric_args = score_metric_args_mv_trimmed,
                                              density_metric_args = density_metric_args_mv_trimmed) 


# ------------------------------------------------------------
# (7) Final benchmark runs, one object per density configuration
# ------------------------------------------------------------

# d = 2: 3 combined SM + KDE Hpi/Hns + smoothed MLE.
res_compare_mv_gaussian_independent_d2 <- run_mv_family_selection_benchmark(
  truth = all_mv_truths$independent_d2,
  estimator_specs = estimator_specs,
  n_rep = 20,
  n_test = 3000,
  save = TRUE,
  save_dir = "resultsMulti",
  save_name = "res_compare_mv_gaussian_independent_d2.rds"
)

res_compare_mv_gaussian_dependent_d2 <- run_mv_family_selection_benchmark(
  truth = all_mv_truths$dependent_d2,
  estimator_specs = estimator_specs,
  n_rep = 20,
  n_test = 3000,
  save = TRUE,
  save_dir = "resultsMulti",
  save_name = "res_compare_mv_gaussian_dependent_d2.rds"
)

# d = 3: 3 combined SM + KDE Hpi/Hns + smoothed MLE.
res_compare_mv_gaussian_independent_d3 <- run_mv_family_selection_benchmark(
  truth = all_mv_truths$independent_d3,
  estimator_specs = estimator_specs,
  n_rep = 20,
  n_test = 3000,
  save = TRUE,
  save_dir = "resultsMulti",
  save_name = "res_compare_mv_gaussian_independent_d3.rds"
)

res_compare_mv_gaussian_dependent_d3 <- run_mv_family_selection_benchmark(
  truth = all_mv_truths$dependent_d3,
  estimator_specs = estimator_specs,
  n_rep = 20,
  n_test = 3000,
  save = TRUE,
  save_dir = "resultsMulti",
  save_name = "res_compare_mv_gaussian_dependent_d3.rds"
)

# d = 4: 3 combined SM + KDE Hpi/Hns + smoothed MLE.
res_compare_mv_gaussian_independent_d4 <- run_mv_family_selection_benchmark(
  truth = all_mv_truths$independent_d4,
  estimator_specs = estimator_specs,
  n_rep = 20,
  n_test = 3000,
  save = TRUE,
  save_dir = "resultsMulti",
  save_name = "res_compare_mv_gaussian_independent_d4.rds"
)

res_compare_mv_gaussian_dependent_d4 <- run_mv_family_selection_benchmark(
  truth = all_mv_truths$dependent_d4,
  estimator_specs = estimator_specs,
  n_rep = 20,
  n_test = 3000,
  save = TRUE,
  save_dir = "resultsMulti",
  save_name = "res_compare_mv_gaussian_dependent_d4.rds"
)









truth_d2 <- make_truth_gaussian_mv(
  d = 2,
  scenario = "independent",
  mean = c(0, 0),
  rho = 0.6
)

# ------------------------------------------------------------
# Estimator: smoothed log-concave MLE only
# ------------------------------------------------------------

mle_spec_only <- list(
  list(
    label = "MLE_smoothed",
    method = "MLE",
    smoothed = TRUE,
    fit_args = list(),
    density_predict_args = list(),
    score_predict_args = list(h = 1e-4),
    score_metric_args = list(
      central_trim = 0.05,
      robust_trim  = 0.01
    ),
    density_metric_args = list(
      central_trim = 0.05
    )
  )
)

# ------------------------------------------------------------
# Run benchmark
# ------------------------------------------------------------

# res_test_mv_gaussian_d2_mle <- run_final_benchmark(
#   sample_sizes = c(50, 200, 500, 1000),
#   family = truth_d2$family,
#   estimator_specs = mle_spec_only,
#   r_sample = truth_d2$r_sample,
#   metrics = c("kl"),
#   n_rep = 5,
#   n_test = 50,
#   true_density = truth_d2$true_density,
#   true_logdensity = truth_d2$true_logdensity,
#   true_score = truth_d2$true_score,
#   truth_name = truth_d2$name,
#   seed = 123,
#   verbose = TRUE,
#   save = TRUE,
#   save_dir = "resultsMulti",
#   save_name = "test_mv_gaussian_d2_mle_nrep5_n1000.rds"
# )




# ------------------------------------------------------------
# Debug test: m = 4, 5 without log-concavity vs weak grid penalty
# ------------------------------------------------------------

make_sm_specs_mv_debug_m45 <- function(m_values = c(4, 5),
                                       ridge = 1e-8,
                                       lc_grid_size = 5L,
                                       weak_lc_penalty = 1e2,
                                       optim_grad = FALSE,
                                       score_metric_args = score_metric_args_mv_trimmed,
                                       density_metric_args = density_metric_args_mv_trimmed) {
  specs <- list()
  
  for (m in m_values) {
    specs <- c(specs, list(
      list(
        label = paste0("SM_pairwise_no_logconcave_m", m),
        method = "SM",
        smoothed = FALSE,
        fit_args = list(
          m = m,
          include_interactions = TRUE,
          standardize = TRUE,
          ridge = ridge,
          log_concave = FALSE
        ),
        density_predict_args = list(),
        score_predict_args = list(),
        score_metric_args = score_metric_args,
        density_metric_args = density_metric_args
      ),
      list(
        label = paste0("SM_pairwise_weak_grid_logconcave_m", m),
        method = "SM",
        smoothed = FALSE,
        fit_args = list(
          m = m,
          include_interactions = TRUE,
          standardize = TRUE,
          ridge = ridge,
          log_concave = TRUE,
          lc_method = "grid",
          lc_grid_size = lc_grid_size,
          lc_penalty = weak_lc_penalty,
          optim_grad = optim_grad
        ),
        density_predict_args = list(),
        score_predict_args = list(),
        score_metric_args = score_metric_args,
        density_metric_args = density_metric_args
      )
    ))
  }
  
  specs
}


estimator_specs_debug_m45 <- make_sm_specs_mv_debug_m45(
  m_values = c(4, 5),
  ridge = 1e-2,
  lc_grid_size = 10L,
  weak_lc_penalty = 1e2,
  optim_grad = TRUE
)

# 
# res_debug_mv_gaussian_dependent_d2_m45 <- run_mv_family_selection_benchmark(
#   truth = all_mv_truths$dependent_d2,
#   estimator_specs = estimator_specs_debug_m45,
#   sample_sizes = c(50, 100, 200, 500, 1000, 5000),
#   metrics = c("score_loss"),
#   n_rep = 5,
#   n_test = 1000,
#   seed = 123,
#   verbose = TRUE,
#   save = TRUE,
#   save_dir = "resultsMulti",
#   save_name = "debug_mv_gaussian_dependent_d2_m45_noLC_vs_weakGridLC.rds"
# )
# 
# 
# plot_final_benchmark(res_debug_mv_gaussian_dependent_d2_m45, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
#                      exclude_normalization_suspect = FALSE)
# 
# plot_final_benchmark(res_debug_mv_gaussian_dependent_d2_m45, metric = "fit_time_sec", center = "mean", interval = "none",
#                      exclude_normalization_suspect = FALSE, log_y = TRUE)
# 
# plot_final_benchmark(res_debug_mv_gaussian_dependent_d2_m45, metric = "lc_n_violated", center = "mean", interval = "none",
#                      exclude_normalization_suspect = FALSE)
# 
# aggregate_final_benchmark(res_debug_mv_gaussian_dependent_d2_m45, metric = "lc_min_eigenvalue", across_runs_center = "mean",
#                           exclude_normalization_suspect = FALSE)
# 
# 




