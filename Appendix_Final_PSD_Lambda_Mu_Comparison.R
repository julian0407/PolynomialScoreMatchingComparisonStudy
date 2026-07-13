
# ---------------------------------------------------------------------------
# Test Framework for coefficent Comparison of PSD and Lambda-Mu optimization
# ---------------------------------------------------------------------------
source("01_Rscripts/03_Test_Framework/Unified_Testing_Framework.R")

# ---------------------------------------------------------------------------
# (1) Plot Function
# ---------------------------------------------------------------------------

# Takes sm_parameterization_comparison object, comparison metric and option if
# logscale should be applied
plot_sm_coef_l2_summary <- function(obj,
                                    metric = "coef_l2_psd_vs_lambda_mu",
                                    log_y = TRUE) {
  # Checks if sm_parameterization_comparison object
  if (!inherits(obj, "sm_parameterization_comparison")) {
    stop("obj must be a sm_parameterization_comparison object.")
  }
  
  # gets data and check if data contains values for this metric
  df <- obj$summary
  if (!metric %in% names(df)) {
    stop(sprintf(
      "Metric '%s' not found. Available columns are: %s",
      metric,
      paste(names(df), collapse = ", ")
    ))
  }
  
  # Generates plot for this metric
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = n,
      y = .data[[metric]],
      color = factor(m),
      group = factor(m)
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Sample Size n",
      y = "Average squared L2 coefficient distance",
      color = "Degree m",
      title = "PSD vs Lambda-Mu: Coefficient L2 Distance"
    ) +
    ggplot2::theme_minimal()
  
  # Optional: Apply log-scale
  if (isTRUE(log_y)) {
    p <- p + ggplot2::scale_y_log10()
  }
  
  p
}


# ---------------------------------------------------------------------------
# (2) Generate Test Objects
# ---------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# (2.1) Generate true univariate distributions: gaussian, Logistic, Gumbel, Laplace
# --------------------------------------------------------------------------------
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

# Create density objects
truth_gaussian <- make_truth_gaussian(-50, 4)
truth_logistic <- make_truth_logistic()
truth_gumbel <- make_truth_gumbel()
truth_laplace <- make_truth_laplace()

# --------------------------------------------------------------------------------
# (2.2) Call the comparison run function to generate test objects for each distribution
# --------------------------------------------------------------------------------

# Comparison Coefficients
sm_gaussian_solution_comparison_ridge <- run_sm_parameterization_comparison(
  sample_sizes = c(50, 100, 200, 500, 1000, 5000),
  m_values = c(1, 2, 3, 4, 5, 6),
  ridge = 0,
  n_rep = 100,
  n_test = 3000,
  seed = 123,
  r_sample = truth_gaussian$r_sample,
  true_score = truth_gaussian$true_score,
  save = TRUE,
  save_dir = "02_Results",
  save_name = "solution_comparison_gaussian_sm_parameterizations_ridge.rds",
  lambda_mu_n_starts = 1
)

sm_logistic_solution_comparison_ridge <- run_sm_parameterization_comparison(
  sample_sizes = c(50, 100, 200, 500, 1000, 5000),
  m_values = c(1, 2, 3, 4, 5, 6),
  ridge = 0,
  n_rep = 100,
  n_test = 3000,
  seed = 123,
  r_sample = truth_logistic$r_sample,
  true_score = truth_logistic$true_score,
  save = TRUE,
  save_dir = "02_Results",
  save_name = "solution_comparison_logistic_sm_parameterizations_ridge.rds",
  lambda_mu_n_starts = 1
)

sm_gumbel_solution_comparison_ridge <- run_sm_parameterization_comparison(
  sample_sizes = c(50, 100, 200, 500, 1000, 5000),
  m_values = c(1, 2, 3, 4, 5, 6),
  ridge = 0,
  n_rep = 100,
  n_test = 3000,
  seed = 123,
  r_sample = truth_gumbel$r_sample,
  true_score = truth_gumbel$true_score,
  save = TRUE,
  save_dir = "02_Results",
  save_name = "solution_comparison_gumbel_sm_parameterizations_ridge.rds",
  lambda_mu_n_starts = 1
)

sm_laplace_solution_comparison_ridge <- run_sm_parameterization_comparison(
  sample_sizes = c(50, 100, 200, 500, 1000, 5000),
  m_values = c(1, 2, 3, 4, 5, 6),
  ridge = 0,
  n_rep = 100,
  n_test = 3000,
  seed = 123,
  r_sample = truth_laplace$r_sample,
  true_score = truth_laplace$true_score,
  save = TRUE,
  save_dir = "02_Results",
  save_name = "solution_comparison_laplace_sm_parameterizations_ridge.rds",
  lambda_mu_n_starts = 1
)


# --------------------------------------------------------------------------------
# (2.3) Plot Results
# --------------------------------------------------------------------------------
gaussian_parameter_compare <- readRDS("02_Results/solution_comparison_gaussian_sm_parameterizations_ridge.rds")
logistic_parameter_compare <- readRDS("02_Results/solution_comparison_logistic_sm_parameterizations_ridge.rds")
gumbel_parameter_compare <- readRDS("02_Results/solution_comparison_gumbel_sm_parameterizations_ridge.rds")
laplace_parameter_compare <- readRDS("02_Results/solution_comparison_laplace_sm_parameterizations_ridge.rds")

plot_sm_coef_l2_summary(gaussian_parameter_compare)
plot_sm_coef_l2_summary(logistic_parameter_compare)
plot_sm_coef_l2_summary(gumbel_parameter_compare)
plot_sm_coef_l2_summary(laplace_parameter_compare)

