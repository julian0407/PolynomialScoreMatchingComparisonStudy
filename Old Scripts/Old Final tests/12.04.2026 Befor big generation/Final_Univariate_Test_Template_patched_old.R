# ============================================================
# Final_Univariate_Test_Template.R
# Vorlage für die univariaten Tests gemäß Methodik
#
# Szenarien:
#   1) log-concave + im Modell: Gaussian
#   2) log-concave, nicht direkt im Modell: Logistic, Gumbel, Laplace
#   3) nicht log-concave: Student-t (heavy tails)
#
# Die Vorlage ist bewusst in Abschnitte unterteilt:
#   A) Setup und Wahrheiten
#   B) estimatorinterne Kandidatentests
#   C) manuelle schätzerübergreifende Vergleiche
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
# Abschnitt A: globale Testeinstellungen
# ------------------------------------------------------------

sample_sizes_main <- c(50, 100, 200, 500, 1000)
sample_sizes_bias_variance <- c(50, 100, 200, 500, 1000)
metrics_main <- c("negloglik", "kl", "score_loss")
seed_main <- 123

# score_loss-Defaults für alle Spezifikationen
# Eine einzelne Konfiguration: robust = "none" | "trim" | "winsor" | "median"
# score_metric_args_default <- list(
#   robust = "median",
#   trim_alpha = 0.05
# )

# Beispiel für mehrere Konfigurationen pro Run:
score_metric_args_default <- list(
  configs = list(
    mean = list(robust = "none"),
    median = list(robust = "median"),
    trim = list(robust = "trim", trim_alpha = 0.05)
  )
)

# density_metric_args_default <- list(
#   robust = "median",
#   trim_alpha = 0.05,
#   outlier_dom_threshold = 0.25
# )

# Beispiel für mehrere Dichte-Konfigurationen pro Run:
density_metric_args_default <- list(
  outlier_dom_threshold = 0.25,
  configs = list(
    mean = list(robust = "none"),
    median = list(robust = "median"),
    trim = list(robust = "trim", trim_alpha = 0.05)
  )
)

# ------------------------------------------------------------
# Abschnitt B: Wahre univariate Verteilungen
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
# Abschnitt C: Kandidatenlisten pro Methodenfamilie
# ------------------------------------------------------------

make_kde_specs_1d <- function(score_metric_args = score_metric_args_default) {
  lapply(c("SJ", "nrd0", "ucv", "bcv"), function(bw) {
    list(
      label = paste0("KDE_", bw),
      method = "KDE",
      smoothed = FALSE,
      fit_args = list(bw = bw),
      density_predict_args = list(),
      score_predict_args = list(),
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args_default
    )
  })
}

make_mle_specs_1d <- function(score_metric_args = score_metric_args_default) {
  list(
    list(
      label = "MLE_unsmoothed",
      method = "MLE",
      smoothed = FALSE,
      fit_args = list(),
      density_predict_args = list(),
      score_predict_args = list(),
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args_default
    ),
    list(
      label = "MLE_smoothed",
      method = "MLE",
      smoothed = TRUE,
      fit_args = list(),
      density_predict_args = list(),
      score_predict_args = list(),
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args_default
    )
  )
}

make_sm_specs_1d <- function(m_values = c(1, 2, 3, 4, 5),
                             score_metric_args = score_metric_args_default) {
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
      score_metric_args = score_metric_args,
      density_metric_args = density_metric_args_default
    )
  })
}

# ------------------------------------------------------------
# Abschnitt D: estimatorinterne Tests
# Für jede Dichte separat innerhalb einer Methodenfamilie
# ------------------------------------------------------------

run_family_selection_benchmark <- function(truth,
                                           estimator_specs,
                                           sample_sizes = sample_sizes_main,
                                           metrics = metrics_main,
                                           n_rep = 50,
                                           n_test = 2000,
                                           seed = seed_main,
                                           verbose = TRUE) {
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
    seed = seed,
    verbose = verbose
  )
}

run_family_bias_variance <- function(truth,
                                     estimator_specs,
                                     sample_sizes = sample_sizes_bias_variance,
                                     n_rep = 50,
                                     seed = seed_main,
                                     verbose = TRUE) {
  eval_grid <- make_eval_grid(
    r_sample = truth$r_sample,
    family = truth$family,
    seed = seed
  )

  run_bias_variance_score_benchmark(
    sample_sizes = sample_sizes,
    family = truth$family,
    estimator_specs = estimator_specs,
    r_sample = truth$r_sample,
    true_score = truth$true_score,
    eval_grid = eval_grid,
    n_rep = n_rep,
    seed = seed,
    verbose = verbose
  )
}

# ------------------------------------------------------------
# Abschnitt E: Vorlagen für estimatorinterne Runs
# ------------------------------------------------------------

# --- 1) Gaussian: log-concave + im Modell ---------------------------------

kde_gaussian_candidates <- run_family_selection_benchmark(
  truth = truth_gaussian,
  estimator_specs = make_kde_specs_1d(),
  n_rep = 50,
  n_test = 3000
)

# plot_metric_comparison(kde_gaussian_candidates, metric = "negloglik", use_median = TRUE)
# plot_metric_comparison(kde_gaussian_candidates, metric = "kl", use_median = TRUE)
# plot_metric_comparison(kde_gaussian_candidates, metric = "score_loss", use_median = TRUE)
# plot_metric_comparison(kde_gaussian_candidates, metric = "fit_time_sec", use_median = TRUE)
# plot_metric_comparison(kde_gaussian_candidates, metric = "total_inference_time_sec", use_median = TRUE)
# 
# aggregate_final_benchmark(kde_gaussian_candidates, metric = "score_loss")



mle_gaussian_candidates <- run_family_selection_benchmark(
  truth = truth_gaussian,
  estimator_specs = make_mle_specs_1d(),
  n_rep = 50,
  n_test = 3000
)

sm_gaussian_candidates <- run_family_selection_benchmark(
  truth = truth_gaussian,
  estimator_specs = make_sm_specs_1d(m_values = c(1, 2, 3, 4, 5)),
  n_rep = 10,
  n_test = 1000
)

plot_final_benchmark(sm_gaussian_candidates, metric = "negloglik__mean", center = "median", interval = "none")
plot_final_benchmark(sm_gaussian_candidates, metric = "negloglik__mean", interval = "none")

plot_final_benchmark(sm_gaussian_candidates, metric = "kl__mean", center = "median", interval = "none")
plot_final_benchmark(sm_gaussian_candidates, metric = "kl__mean", interval = "none")


plot_final_benchmark(sm_gaussian_candidates, metric = "score_loss__trim", center = "median", interval = "none")
plot_final_benchmark(sm_gaussian_candidates, metric = "score_loss__trim", interval = "none")

sm_gaussian_candidates_filtered <- subset_final_benchmark(
  sm_gaussian_candidates,
  drop_method_labels = c("SM_m6"),
  drop_n = c(50, 100)
)

plot_final_benchmark(sm_gaussian_candidates_filtered, metric = "negloglik__mean", center = "median", interval = "none")
plot_final_benchmark(sm_gaussian_candidates_filtered, metric = "negloglik__mean", interval = "none")

plot_final_benchmark(sm_gaussian_candidates_filtered, metric = "kl__mean", center = "median", interval = "none")
plot_final_benchmark(sm_gaussian_candidates_filtered, metric = "kl__mean", interval = "none")


plot_final_benchmark(sm_gaussian_candidates_filtered, metric = "score_loss__trim", center = "median", interval = "none")
plot_final_benchmark(sm_gaussian_candidates_filtered, metric = "score_loss__trim", interval = "none")

# plot_metric_comparison(sm_gaussian_candidates, metric = "negloglik", use_median = TRUE)
# plot_metric_comparison(sm_gaussian_candidates, metric = "kl", use_median = TRUE)
# plot_metric_comparison(sm_gaussian_candidates, metric = "score_loss", use_median = TRUE)
# plot_metric_comparison(sm_gaussian_candidates, metric = "fit_time_sec", use_median = TRUE)
# plot_metric_comparison(sm_gaussian_candidates, metric = "total_inference_time_sec", use_median = TRUE)
# 
# 
# plot_metric_comparison(sm_gaussian_candidates, metric = "kappa_raw", log_y = TRUE)
# plot_metric_comparison(sm_gaussian_candidates, metric = "eigmin_raw")
# plot_metric_comparison(sm_gaussian_candidates, metric = "rcond_raw", log_y = TRUE)

#
sm_gaussian_bv <- run_family_bias_variance(
  truth = truth_gaussian,
  estimator_specs = make_sm_specs_1d(m_values = c(1, 2, 3, 4, 5)),
  n_rep = 50
)

# --- 2) Logistic: log-concave, nicht direkt im Modell ---------------------

kde_logistic_candidates <- run_family_selection_benchmark(
  truth = truth_logistic,
  estimator_specs = make_kde_specs_1d(),
  n_rep = 50,
  n_test = 3000
)

mle_logistic_candidates <- run_family_selection_benchmark(
  truth = truth_logistic,
  estimator_specs = make_mle_specs_1d(),
  n_rep = 50,
  n_test = 3000
)

sm_logistic_candidates <- run_family_selection_benchmark(
  truth = truth_logistic,
  estimator_specs = make_sm_specs_1d(m_values = c(1,2,3,4,5, 6)),
  n_rep = 10,
  n_test = 1000
)

plot_final_benchmark(sm_logistic_candidates, metric = "negloglik", center = "median", interval = "none")
plot_final_benchmark(sm_logistic_candidates, metric = "kl__mean", center = "median", interval = "none")
plot_final_benchmark(sm_logistic_candidates, metric = "score_loss")

sm_logistic_candidates_filtered <- subset_final_benchmark(
  sm_logistic_candidates,
  drop_method_labels = c("SM_m6", "SM_m5")
  # drop_n = c(50, 100)
)

plot_final_benchmark(sm_logistic_candidates_filtered, metric = "negloglik", center = "median", interval = "none")
plot_final_benchmark(sm_logistic_candidates_filtered, metric = "kl__mean", center = "median", interval = "none")
plot_final_benchmark(sm_logistic_candidates_filtered, metric = "score_loss")

# --- 3) Gumbel: log-concave, nicht direkt im Modell -----------------------

kde_gumbel_candidates <- run_family_selection_benchmark(
  truth = truth_gumbel,
  estimator_specs = make_kde_specs_1d(),
  n_rep = 50,
  n_test = 3000
)

mle_gumbel_candidates <- run_family_selection_benchmark(
  truth = truth_gumbel,
  estimator_specs = make_mle_specs_1d(),
  n_rep = 50,
  n_test = 3000
)

sm_gumbel_candidates <- run_family_selection_benchmark(
  truth = truth_gumbel,
  estimator_specs = make_sm_specs_1d(m_values = c(1, 2, 3, 4, 5, 6)),
  n_rep = 50,
  n_test = 3000
)

# --- 4) Laplace: log-concave, eckiger Spezialfall -------------------------

kde_laplace_candidates <- run_family_selection_benchmark(
  truth = truth_laplace,
  estimator_specs = make_kde_specs_1d(),
  n_rep = 50,
  n_test = 3000
)

mle_laplace_candidates <- run_family_selection_benchmark(
  truth = truth_laplace,
  estimator_specs = make_mle_specs_1d(),
  n_rep = 50,
  n_test = 3000
)

sm_laplace_candidates <- run_family_selection_benchmark(
  truth = truth_laplace,
  estimator_specs = make_sm_specs_1d(m_values = c(1, 2, 3, 4, 5, 6)),
  n_rep = 50,
  n_test = 3000
)

# --- 5) Student-t: nicht log-concave / heavy tails ------------------------

kde_student_candidates <- run_family_selection_benchmark(
  truth = truth_student,
  estimator_specs = make_kde_specs_1d(),
  n_rep = 50,
  n_test = 3000
)

mle_student_candidates <- run_family_selection_benchmark(
  truth = truth_student,
  estimator_specs = make_mle_specs_1d(),
  n_rep = 50,
  n_test = 3000
)

sm_student_candidates <- run_family_selection_benchmark(
  truth = truth_student,
  estimator_specs = make_sm_specs_1d(m_values = c(1, 2, 3, 4, 5, 6)),
  n_rep = 50,
  n_test = 3000
)

# ------------------------------------------------------------
# Abschnitt F: manuelle schätzerübergreifende Vergleiche
# Hier trägst du pro Dichte die intern selektierten Konfigurationen ein.
# ------------------------------------------------------------

# Beispiel: Falls du nach den Kandidatentests entscheidest,
#   - Gaussian: KDE_SJ, MLE_unsmoothed, SM_m2
#   - Logistic: KDE_ucv, MLE_smoothed, SM_m4
# dann definierst du dir die Listen manuell und vergleichst sie direkt.

manual_compare_gaussian <- list(
  list(
    label = "KDE_manual",
    method = "KDE",
    smoothed = FALSE,
    fit_args = list(bw = "SJ"),
    density_predict_args = list(),
    score_predict_args = list(),
    score_metric_args = score_metric_args_default
  ),
  list(
    label = "MLE_manual",
    method = "MLE",
    smoothed = FALSE,
    fit_args = list(),
    density_predict_args = list(),
    score_predict_args = list(),
    score_metric_args = score_metric_args_default
  ),
  list(
    label = "SM_manual",
    method = "SM",
    smoothed = FALSE,
    fit_args = list(m = 2, standardize = TRUE, ridge = 0),
    density_predict_args = list(),
    score_predict_args = list(),
    score_metric_args = score_metric_args_default
  )
)

manual_compare_logistic <- manual_compare_gaussian
manual_compare_gumbel <- manual_compare_gaussian
manual_compare_student <- manual_compare_gaussian
manual_compare_laplace <- manual_compare_gaussian

# Beispielaufrufe:
res_compare_gaussian <- run_family_selection_benchmark(
  truth = truth_gaussian,
  estimator_specs = manual_compare_gaussian,
  n_rep = 100,
  n_test = 5000
)

res_compare_logistic <- run_family_selection_benchmark(
  truth = truth_logistic,
  estimator_specs = manual_compare_logistic,
  n_rep = 100,
  n_test = 5000
)

res_compare_gumbel <- run_family_selection_benchmark(
  truth = truth_gumbel,
  estimator_specs = manual_compare_gumbel,
  n_rep = 100,
  n_test = 5000
)

res_compare_student <- run_family_selection_benchmark(
  truth = truth_student,
  estimator_specs = manual_compare_student,
  n_rep = 100,
  n_test = 5000
)

res_compare_laplace <- run_family_selection_benchmark(
  truth = truth_laplace,
  estimator_specs = manual_compare_laplace,
  n_rep = 100,
  n_test = 5000
)

# ------------------------------------------------------------
# Abschnitt G: typische Auswertung nach einem Run
# ------------------------------------------------------------

aggregate_final_benchmark(res_compare_gaussian, "negloglik")
aggregate_final_benchmark(res_compare_gaussian, "kl")
aggregate_final_benchmark(res_compare_gaussian, "score_loss")
aggregate_final_benchmark(res_compare_gaussian, "fit_time_sec")
aggregate_final_benchmark(res_compare_gaussian, "score_inference_time_sec")
aggregate_final_benchmark(res_compare_gaussian, "kappa_raw")
aggregate_final_benchmark(res_compare_gaussian, "eigmin_raw")

# Beispiel: kleine Sample Sizes nur für die Plot-Auswertung ausblenden
# plot_final_benchmark(res_compare_gaussian, "score_loss", drop_n = c(50, 100))
# plot_final_benchmark(res_compare_gaussian, "kl", keep_n = c(200, 500, 1000))
#
# Falls du mit Bias-Variance-Runs arbeitest:
# sm_gaussian_bv_filtered <- subset_score_bias_variance_benchmark(
#   sm_gaussian_bv,
#   drop_n = c(50, 100)
# )
# plot_score_bias_variance(sm_gaussian_bv, "integrated_mse", drop_n = c(50, 100))
# plot_score_bias_variance(sm_gaussian_bv, "integrated_variance", keep_n = c(200, 500, 1000))



# Beispiele fuer Debugging / Replay einzelner problematischer Runs
# outliers <- debug_benchmark_outliers(res_compare_gaussian, metric_pattern = "^kl")
# replay <- replay_benchmark_run(
#   res_compare_gaussian,
#   method_label = outliers$method_label[1],
#   n = outliers$n[1],
#   repetition = outliers$repetition[1]
# )
# head(replay$pointwise)

# Beispiele fuer nachtraegliches Filtern
# p_kde_only <- plot_final_benchmark(res_compare_gaussian, "score_loss", keep_method_labels = c("KDE_SJ", "KDE_ucv"))
# p_drop_small_and_m1 <- plot_score_bias_variance(sm_gaussian_bv, "integrated_mse", drop_n = c(50, 100), drop_method_labels = "SM_m1")
# p_mle_only <- plot_final_benchmark(res_compare_gaussian, "kl", keep_method_labels = c("MLE_smoothed", "MLE_unsmoothed"))
