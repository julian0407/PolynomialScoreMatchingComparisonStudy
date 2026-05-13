# ============================================================
# Final_Univariate_Benchmark_Runs.R
# Template for univariate benchmark runs 
#
# Core metrics:
#   - kl
#   - score_loss
# Optional variants from metric args:
#   - *_central      : metric on empirical bulk region only
#   - score_loss_trim: score metric after trimming the largest
#                      pointwise score losses
# ============================================================

source("01_Rscripts/Unified_Testing_Framework.R")

# ------------------------------------------------------------
# (1) Global settings
# ------------------------------------------------------------
sample_sizes_main <- c(50, 100)
metrics_main <- c("kl", "score_loss")
seed_main <- 123

# central_trim evaluates the metric on the empirical bulk of the test set.
# robust_trim only applies to score_loss and trims the largest pointwise
# score losses after the score errors have been computed.
# args are set in specific estimators
score_metric_args_default <- list(
  central_trim = NULL,
  robust_trim = NULL
)

density_metric_args_default <- list(
  central_trim = NULL
)

# ------------------------------------------------------------
# (2) True univariate distributions
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

# ------------------------------------------------------------
# (3) Create true density objects
# ------------------------------------------------------------

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
# (4) Create estimator candidate lists per method family
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

make_sm_specs_1d_noridge <- function(m_values = c(1, 2, 3, 4, 5, 6),
                                     score_metric_args = score_metric_args_default,
                                     density_metric_args = density_metric_args_default) {
  lapply(m_values, function(m) {
    list(
      label = paste0("SM_m", m, "_noridge_std"),
      method = "SM",
      smoothed = FALSE,
      fit_args = list(
        m = m,
        standardize = TRUE,
        ridge = 0
      ),
      density_predict_args = list(
        subdivisions = 200L,
        rel.tol = 1e-8,
        stop_on_failure = FALSE
      ),
      score_predict_args = list(),
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args
    )
  })
}

make_sm_specs_1d_ridge <- function(m_values = c(1, 2, 3, 4, 5, 6),
                                   ridge = 1e-2,
                                   score_metric_args = score_metric_args_default,
                                   density_metric_args = density_metric_args_default) {
  lapply(m_values, function(m) {
    list(
      label = paste0("SM_m", m, "_ridge1e-02_std"),
      method = "SM",
      smoothed = FALSE,
      fit_args = list(
        m = m,
        standardize = TRUE,
        ridge = ridge
      ),
      density_predict_args = list(
        subdivisions = 200L,
        rel.tol = 1e-8,
        stop_on_failure = FALSE
      ),
      score_predict_args = list(),
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args
    )
  })
}

# ------------------------------------------------------------
# (5) Benchmark wrapper 
# ------------------------------------------------------------
run_family_selection_benchmark <- function(truth,
                                           estimator_specs,
                                           sample_sizes = sample_sizes_main,
                                           metrics = metrics_main,
                                           n_rep = 20,
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

# ------------------------------------------------------------
# (6)) Example runs
# ------------------------------------------------------------

# ------------------------------------------------------------
# (6.1) Set central and trim values
# ------------------------------------------------------------

score_metric_args_trimmed <- list(
  central_trim = 0.05,
  robust_trim  = 0.01
)

density_metric_args_trimmed <- list(
  central_trim = 0.05
)

# ------------------------------------------------------------
# (6.2) Tests within estimator families
# ------------------------------------------------------------

# --- 1) Gaussian -----------------------------------------------------------

kde_gaussian_candidates <- run_family_selection_benchmark(
  truth = truth_gaussian,
  estimator_specs = make_kde_specs_1d(
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

mle_gaussian_candidates <- run_family_selection_benchmark(
  truth = truth_gaussian,
  estimator_specs = make_mle_specs_1d(
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

sm_gaussian_candidates2_noridge <- run_family_selection_benchmark(
  truth = truth_gaussian,
  estimator_specs = make_sm_specs_1d_noridge(
    m_values = c(1, 2, 3, 4, 5, 6),
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

sm_gaussian_candidates2_ridge <- run_family_selection_benchmark(
  truth = truth_gaussian,
  estimator_specs = make_sm_specs_1d_ridge(
    m_values = c(1, 2, 3, 4, 5, 6),
    ridge = 1e-2,
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

# --- 2) Logistic -----------------------------------------------------------

kde_logistic_candidates <- run_family_selection_benchmark(
  truth = truth_logistic,
  estimator_specs = make_kde_specs_1d(
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

mle_logistic_candidates <- run_family_selection_benchmark(
  truth = truth_logistic,
  estimator_specs = make_mle_specs_1d(
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

sm_logistic_candidates_noridge <- run_family_selection_benchmark(
  truth = truth_logistic,
  estimator_specs = make_sm_specs_1d_noridge(
    m_values = c(1, 2, 3, 4, 5, 6),
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

sm_logistic_candidates_ridge <- run_family_selection_benchmark(
  truth = truth_logistic,
  estimator_specs = make_sm_specs_1d_ridge(
    m_values = c(1, 2, 3, 4, 5, 6),
    ridge = 1e-2,
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

# --- 3) Gumbel -------------------------------------------------------------

kde_gumbel_candidates <- run_family_selection_benchmark(
  truth = truth_gumbel,
  estimator_specs = make_kde_specs_1d(
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

mle_gumbel_candidates <- run_family_selection_benchmark(
  truth = truth_gumbel,
  estimator_specs = make_mle_specs_1d(
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

sm_gumbel_candidates_noridge <- run_family_selection_benchmark(
  truth = truth_gumbel,
  estimator_specs = make_sm_specs_1d_noridge(
    m_values = c(1, 2, 3, 4, 5, 6),
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

sm_gumbel_candidates_ridge <- run_family_selection_benchmark(
  truth = truth_gumbel,
  estimator_specs = make_sm_specs_1d_ridge(
    m_values = c(1, 2, 3, 4, 5, 6),
    ridge = 1e-2,
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

# --- 4) Laplace ------------------------------------------------------------

kde_laplace_candidates <- run_family_selection_benchmark(
  truth = truth_laplace,
  estimator_specs = make_kde_specs_1d(
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

mle_laplace_candidates <- run_family_selection_benchmark(
  truth = truth_laplace,
  estimator_specs = make_mle_specs_1d(
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

sm_laplace_candidates_noridge <- run_family_selection_benchmark(
  truth = truth_laplace,
  estimator_specs = make_sm_specs_1d_noridge(
    m_values = c(1, 2, 3, 4, 5, 6),
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

sm_laplace_candidates_ridge <- run_family_selection_benchmark(
  truth = truth_laplace,
  estimator_specs = make_sm_specs_1d_ridge(
    m_values = c(1, 2, 3, 4, 5, 6),
    ridge = 1e-2,
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

# --- 5) Student-t ----------------------------------------------------------

kde_student_candidates <- run_family_selection_benchmark(
  truth = truth_student,
  estimator_specs = make_kde_specs_1d(
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

mle_student_candidates <- run_family_selection_benchmark(
  truth = truth_student,
  estimator_specs = make_mle_specs_1d(
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

sm_student_candidates_noridge <- run_family_selection_benchmark(
  truth = truth_student,
  estimator_specs = make_sm_specs_1d_noridge(
    m_values = c(1, 2, 3, 4, 5, 6),
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

sm_student_candidates_ridge <- run_family_selection_benchmark(
  truth = truth_student,
  estimator_specs = make_sm_specs_1d_ridge(
    m_values = c(1, 2, 3, 4, 5, 6),
    ridge = 1e-2,
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

# ------------------------------------------------------------
# (6.3) Helper for Comparison across estimators
# ------------------------------------------------------------

# Create tags of methods that should be applied for each density
initial_best_guess <- list(
  gaussian = list(
    kde = "KDE_SJ",
    mle = "MLE_smoothed",
    sm1 = "SM_m1_noridge_std",
    sm2 = "SM_m2_noridge_std"
  ),
  logistic = list(
    kde = "KDE_ucv",
    mle = "MLE_smoothed",
    sm1 = "SM_m3_ridge1e-02_std",
    sm2 = "SM_m4_ridge1e-02_std"
  ),
  gumbel = list(
    kde = "KDE_ucv",
    mle = "MLE_smoothed",
    sm1 = "SM_m3_ridge1e-02_std",
    sm2 = "SM_m4_ridge1e-02_std"
  ),
  laplace = list(
    kde = "KDE_SJ",
    mle = "MLE_unsmoothed",
    sm1 = "SM_m2_ridge1e-02_std",
    sm2 = "SM_m3_ridge1e-02_std"
  ),
  student = list(
    kde = "KDE_SJ",
    mle = "MLE_smoothed",
    sm1 = "SM_m2_ridge1e-02_std",
    sm2 = "SM_m3_ridge1e-02_std"
  )
)

# Create pool of possible estimators
manual_spec_pool <- c(
  make_kde_specs_1d(
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  make_mle_specs_1d(
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  make_sm_specs_1d_noridge(
    m_values = c(1, 2, 3, 4, 5, 6),
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  ),
  make_sm_specs_1d_ridge(
    m_values = c(1, 2, 3, 4, 5, 6),
    ridge = 1e-2,
    score_metric_args = score_metric_args_trimmed,
    density_metric_args = density_metric_args_trimmed
  )
)

# Get estimators by labels
get_specs_by_labels <- function(labels, spec_pool = manual_spec_pool) {
  keep <- vapply(spec_pool, function(sp) sp$label %in% labels, logical(1))
  out <- spec_pool[keep]
  
  found_labels <- vapply(out, function(sp) sp$label, character(1))
  missing_labels <- setdiff(labels, found_labels)
  if (length(missing_labels) > 0L) {
    stop("Missing estimator specs for labels: ", paste(missing_labels, collapse = ", "))
  }
  
  out[match(labels, found_labels)]
}

# Get final estimator list per density based on initial list above
make_manual_compare_specs <- function(truth_name) {
  guess <- initial_best_guess[[truth_name]]
  if (is.null(guess)) stop("Unknown truth_name: ", truth_name)
  
  get_specs_by_labels(c(guess$kde, guess$mle, guess$sm1, guess$sm2))
}

# ------------------------------------------------------------
#  (6.4) Final Benchmark Runs for Comparison across estimator families
# ------------------------------------------------------------
manual_compare_gaussian <- make_manual_compare_specs("gaussian")
manual_compare_logistic <- make_manual_compare_specs("logistic")
manual_compare_gumbel  <- make_manual_compare_specs("gumbel")
manual_compare_laplace <- make_manual_compare_specs("laplace")
manual_compare_student <- make_manual_compare_specs("student")

res_compare_gaussian <- run_family_selection_benchmark(
  truth = truth_gaussian,
  estimator_specs = manual_compare_gaussian,
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

res_compare_logistic <- run_family_selection_benchmark(
  truth = truth_logistic,
  estimator_specs = manual_compare_logistic,
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

res_compare_gumbel <- run_family_selection_benchmark(
  truth = truth_gumbel,
  estimator_specs = manual_compare_gumbel,
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results",
  seed = 187
)

res_compare_laplace <- run_family_selection_benchmark(
  truth = truth_laplace,
  estimator_specs = manual_compare_laplace,
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)

res_compare_student <- run_family_selection_benchmark(
  truth = truth_student,
  estimator_specs = manual_compare_student,
  n_rep = 2,
  n_test = 100,
  save = TRUE,
  save_dir = "results"
)


