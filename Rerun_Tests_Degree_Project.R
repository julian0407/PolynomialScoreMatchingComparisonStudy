# ============================================================
# Rerun_Tests_Degree_Project.R
# Rerun Tests based on saved data objects in project repository
# To rerun complete fit -> Source Final_Univariate_Test_Template
# ============================================================

source("01_Rscripts/03_Test_Framework/Unified_Testing_Framework.R")

# ------------------------------------------------------------
# (1) Univariate Tests on Score Loss to Gaussian density
# ------------------------------------------------------------

# Load Final Benchmark objects for Gaussian (Ridge and no ridge)
sm_gaussian_candidates_noridge <- readRDS("02_Results/final_gaussian_sm_noridge.rds")
sm_gaussian_candidates_ridge <- readRDS("02_Results/final_gaussian_sm_ridge.rds")

# Score Loss, Mean, No rdige
renaming <- c(
  "SM_m1_noridge_std" = "SM1",
  "SM_m2_noridge_std" = "SM2",
  "SM_m3_noridge_std" = "SM3",
  "SM_m4_noridge_std" = "SM4",
  "SM_m5_noridge_std" = "SM5",
  "SM_m6_noridge_std" = "SM6"
)
plot_final_benchmark(sm_gaussian_candidates_noridge, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, method_label_map = renaming)
aggregate_final_benchmark(sm_gaussian_candidates_noridge, metric = "score_loss", across_runs_center = "mean",
                          exclude_normalization_suspect = FALSE)
# Average condition number per estimator, No ridge
average_condition_number_by_estimator(sm_gaussian_candidates_noridge)

# Score Loss No ridge, Mean / Median Ratio
score_loss_mean_median_factor(
  sm_gaussian_candidates_noridge,
  estimator_pattern = "^SM_m(5|6)_",
  metric = "score_loss",
  exclude_normalization_suspect = FALSE
)

# Score Loss No ridge, proportion of extreme runs per n
tail_probability_table(
  sm_gaussian_candidates_noridge,
  estimator_pattern = "^SM_m(5|6)_",
  thresholds = c(1, 0.1, 0.05)
)

# Filter top3 outlier runs
outliers_gaussian_sm_noridge <- debug_benchmark_outliers(sm_gaussian_candidates_noridge, metric_pattern = "^score_loss$", top_n = 1)
top3_outliers_gaussian_sm_noridge <- head(outliers_gaussian_sm_noridge, 3)
top3_outliers_gaussian_sm_noridge

# Get gap summary for the different benchmark configurations
gap_summary_gaussian_no_ridge  <- summarise_gap_score_compact(benchmark_obj = sm_gaussian_candidates_noridge, group_vars = c("method_label", "method", "n"), 
                                                              gap_vars = c("left_gap", "right_gap"), metric = "score_loss"
)

# Get topn1000 and method SM6
gap_summary_gaussian_no_ridge_n1000_SM6 <- gap_summary_gaussian_no_ridge[gap_summary_gaussian_no_ridge$n == 1000 &
                                                     gap_summary_gaussian_no_ridge$method_label == "SM_m6_noridge_std",]
gap_summary_gaussian_no_ridge_n1000_SM6
# Check if the maximal gap corresponds also to the largest outlier run from above (run seed= 54768478) and identify size of max_gap
gap_summary_gaussian_no_ridge_n1000_SM6[
  , c("method_label", "n", "max_gap", "run_seed")
]
# replay the run that generates the max training test gap
Gap_run_top1 <- replay_benchmark_run(sm_gaussian_candidates_noridge, "SM_m6_noridge_std", 1000, run_seed = 54768478)
# get index that corresponds to the maximum pointwise score loss
which.max(Gap_run_top1[["pointwise"]][["score_loss_point"]])
# Get largest test sample
max(Gap_run_top1[["x_test"]])
# Check if this test sample coincide with sample that generates largest pointwise score loss
Gap_run_top1[["x_test"]][2017]
# Check: Does this sample generate the largest training test gap (= 8.4)
Gap_run_top1[["x_test"]][2017]-max(Gap_run_top1[["x_train"]])
# Get the value of the largest and second largest pointwise score loss
Gap_run_top1[["pointwise"]][["score_loss_point"]][2017]
max(Gap_run_top1[["pointwise"]][["score_loss_point"]][-2017])

# Check second and third laregst gaps
max(Gap_run_top1[["x_test"]][-2017])-max(Gap_run_top1[["x_train"]][-2017])
min(Gap_run_top1[["x_train"]][-2017])-min(Gap_run_top1[["x_test"]][-2017])
which.max(Gap_run_top1[["pointwise"]][["score_loss_point"]][-2017])
min(Gap_run_top1[["x_train"]])-Gap_run_top1[["x_test"]][143]

# Score Loss No ridge, mean rations of test sample within (0) and out of (1) training sample
summarise_gap_score_compact(benchmark_obj = sm_gaussian_candidates_noridge, group_vars = c("method_label", "method", "n"), 
                                                              gap_vars = c("left_gap", "right_gap"), metric = "score_loss"
)

# Plot central score loss
plot_final_benchmark(sm_gaussian_candidates_noridge, metric = "score_loss_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, method_label_map = renaming)

# Average condition number per estimator, ridge
average_condition_number_by_estimator(sm_gaussian_candidates_ridge)

# plot for score loss with ridge
renaming_Ridge <- c(
  "SM_m1_ridge1e-02_std" = "SM1",
  "SM_m2_ridge1e-02_std" = "SM2",
  "SM_m3_ridge1e-02_std" = "SM3",
  "SM_m4_ridge1e-02_std" = "SM4",
  "SM_m5_ridge1e-02_std" = "SM5",
  "SM_m6_ridge1e-02_std" = "SM6"
)
plot_final_benchmark(sm_gaussian_candidates_ridge, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, method_label_map = renaming_Ridge)
plot_final_benchmark(sm_gaussian_candidates_ridge, metric = "score_loss_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, method_label_map = renaming_Ridge)

compare_factor(
  obj_left  = sm_gaussian_candidates_ridge,
  obj_right = sm_gaussian_candidates_noridge,
  metric_left  = "score_loss",
  metric_right = "score_loss",
  exclude_normalization_suspect = TRUE
)

# Bias Variance Tradeoff
# drop sm6
only_show_sm_methods <- c(
  "SM_m6_ridge1e-02_std"
)
# --- 2) Logistic -----------------------------------------------------------
sm_logistic_candidates_ridge <- readRDS("02_Results/final_logistic_sm_ridge.rds")
plot_final_benchmark(sm_logistic_candidates_ridge, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE,  method_label_map = renaming_Ridge , drop_method_labels = only_show_sm_methods)
# --- 3) Gumbel -------------------------------------------------------------
sm_gumbel_candidates_ridge <- readRDS("02_Results/final_gumbel_sm_ridge.rds")
plot_final_benchmark(sm_gumbel_candidates_ridge, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE,  method_label_map = renaming_Ridge, drop_method_labels = only_show_sm_methods )



# ------------------------------------------------------------
# (2) Univariate Tests on KL loss to Gaussian density
# ------------------------------------------------------------

# KL Loss analysis
plot_final_benchmark(sm_gaussian_candidates_noridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, method_label_map = renaming)

# proportion of outlier runs
tail_probability_table(
  sm_gaussian_candidates_noridge,
  estimator_pattern = "^SM_m(5|6)_",
  metric = "kl",
  thresholds = c(1, 0.5, 0.1),
  exclude_normalization_suspect = TRUE
) 


# ---------------------------------------------------------------
# (3) Benchmarking SM in KL loss against log-concave MLE and KDE
# ---------------------------------------------------------------
renaming_Mixed <- c(
  "SM_m1_ridge1e-02_std" = "SM1",
  "SM_m2_ridge1e-02_std" = "SM2",
  "SM_m3_ridge1e-02_std" = "SM3",
  "SM_m4_ridge1e-02_std" = "SM4",
  "SM_m5_ridge1e-02_std" = "SM5",
  "SM_m6_ridge1e-02_std" = "SM6",
  "KDE_bcv" = "KDE BCV",
  "KDE_nrd0" = "KDE NRD0",
  "KDE_SJ" = "KDE SJ",
  "KDE_ucv" = "KDE UCV",
  "MLE_smoothed" = "MLE",
  "MLE_unsmoothed" = "MLE"
)
# --------------------------------------------------------------------
# (3.1) Preselection of degree parameter in SM and bandwidth method KDE
# --------------------------------------------------------------------
# --- 1) Gaussian -----------------------------------------------------------
sm_gaussian_candidates_ridge <- readRDS("02_Results/final_gaussian_sm_ridge.rds")
kde_gaussian_candidates <- readRDS("02_Results/final_gaussian_kde.rds")
plot_final_benchmark(sm_gaussian_candidates_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, method_label_map = renaming_Mixed)
plot_final_benchmark(kde_gaussian_candidates, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, method_label_map = renaming_Mixed)
# --- 2) Logistic -----------------------------------------------------------
sm_logistic_candidates_ridge <- readRDS("02_Results/final_logistic_sm_ridge.rds")
kde_logistic_candidates <- readRDS("02_Results/final_logistic_kde.rds")
plot_final_benchmark(sm_logistic_candidates_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE,  method_label_map = renaming_Mixed)
plot_final_benchmark(kde_logistic_candidates, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, method_label_map = renaming_Mixed)
# --- 3) Gumbel -------------------------------------------------------------
sm_gumbel_candidates_ridge <- readRDS("02_Results/final_gumbel_sm_ridge.rds")
kde_gumbel_candidates <- readRDS("02_Results/final_gumbel_kde.rds")
plot_final_benchmark(sm_gumbel_candidates_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE,  method_label_map = renaming_Mixed)
plot_final_benchmark(kde_gumbel_candidates, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, method_label_map = renaming_Mixed)
# --- 4) Laplace -------------------------------------------------------------
sm_laplace_candidates_ridge <- readRDS("02_Results/final_laplace_sm_ridge.rds")
kde_laplace_candidates <- readRDS("02_Results/final_laplace_kde.rds")
plot_final_benchmark(sm_laplace_candidates_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE,  method_label_map = renaming_Mixed)
plot_final_benchmark(kde_laplace_candidates, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, method_label_map = renaming_Mixed)

# --------------------------------------------------------------------
# (3.1) Benchmarking
# --------------------------------------------------------------------

# --- 1) Gaussian -----------------------------------------------------------
gaussian_mixed <- readRDS("02_Results/final_gaussian_mixed.rds")
plot_final_benchmark(gaussian_mixed, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE, method_label_map = renaming_Mixed)

# --- 2) Logistic -----------------------------------------------------------

logistic_mixed <- readRDS("02_Results/final_logistic_mixed.rds")
plot_final_benchmark(logistic_mixed, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE, method_label_map = renaming_Mixed)

# --- 3) Gumbel -------------------------------------------------------------
gumbel_mixed <- readRDS("02_Results/final_gumbel_mixed.rds")
plot_final_benchmark(gumbel_mixed, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, method_label_map = renaming_Mixed)

# --- 4) Laplace -------------------------------------------------------------
laplace_mixed <- readRDS("02_Results/final_laplace_mixed.rds")
plot_final_benchmark(laplace_mixed, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, method_label_map = renaming_Mixed)




# ------------------------------------------------------------
# (4) Multivaraiate Tests on Gaussian density
# ------------------------------------------------------------

sm_gaussian_cmultivariate_d2 <- readRDS("02_Results/Multivariate/res_compare_mv_gaussian_independent_d2.rds")
sm_gaussian_cmultivariate_d3 <- readRDS("02_Results/Multivariate/res_compare_mv_gaussian_independent_d3.rds")
sm_gaussian_cmultivariate_d4 <- readRDS("02_Results/Multivariate/res_compare_mv_gaussian_independent_d4.rds")
sm_gaussian_cmultivariate_d2_dependent <- readRDS("02_Results/Multivariate/res_compare_mv_gaussian_dependent_d2.rds")
sm_gaussian_cmultivariate_d3_dependent <- readRDS("02_Results/Multivariate/res_compare_mv_gaussian_dependent_d3.rds")
sm_gaussian_cmultivariate_d4_dependent <- readRDS("02_Results/Multivariate/res_compare_mv_gaussian_dependent_d4.rds")

mle_gaussian_d2_independent <- readRDS("02_Results/Multivariate/test_mv_gaussian_d2_mle_nrep5_n1000.rds")
mle_gaussian_d3_independent <- readRDS("02_Results/Multivariate/test_mv_gaussian_d3_mle_nrep5_n1000.rds")
mle_gaussian_d4_independent <- readRDS("02_Results/Multivariate/test_mv_gaussian_d4_mle_nrep5_n1000.rds")

only_show_sm_methods <- c(
  "SM_grid_logconcave_m1",
  "SM_no_logconcave_m1",
  "SM_grid_logconcave_m2",
  "SM_no_logconcave_m2",
  "SM_grid_logconcave_m5",
  "SM_no_logconcave_m5"
)

rename_sm_methods_multi <- c(
  "SM_grid_logconcave_m1" = "SM1 Grid",
  "SM_no_logconcave_m1" = "SM1 No Grid" ,
  "SM_grid_logconcave_m2" = "SM2 Grid",
  "SM_no_logconcave_m2" = "SM2 No Grid",
  "SM_grid_logconcave_m5" = "SM5 Grid",
  "SM_no_logconcave_m5" = "SM5 No Grid"
)

# Convergence
plot_final_benchmark(sm_gaussian_cmultivariate_d2, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, keep_method_labels = only_show_sm_methods,  method_label_map = rename_sm_methods_multi)
plot_final_benchmark(sm_gaussian_cmultivariate_d3, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE,  keep_method_labels = only_show_sm_methods, method_label_map = rename_sm_methods_multi)
plot_final_benchmark(sm_gaussian_cmultivariate_d2_dependent, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE, keep_method_labels = only_show_sm_methods, method_label_map = rename_sm_methods_multi)
plot_final_benchmark(sm_gaussian_cmultivariate_d3_dependent, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE,  keep_method_labels = only_show_sm_methods, method_label_map = rename_sm_methods_multi)

# Log-concavity
aggregate_final_benchmark(sm_gaussian_cmultivariate_d3_dependent, metric = "lc_n_violated", across_runs_center = "mean",
                          exclude_normalization_suspect = FALSE)
aggregate_final_benchmark(sm_gaussian_cmultivariate_d3_dependent, metric = "lc_min_eigenvalue", across_runs_center = "mean",
                          exclude_normalization_suspect = FALSE)

# Runtime
aggregate_final_benchmark(sm_gaussian_cmultivariate_d2_dependent, metric = "fit_time_sec", across_runs_center = "mean",
                          exclude_normalization_suspect = FALSE)
aggregate_final_benchmark(sm_gaussian_cmultivariate_d3_dependent, metric = "fit_time_sec", across_runs_center = "mean",
                          exclude_normalization_suspect = FALSE)
aggregate_final_benchmark(sm_gaussian_cmultivariate_d4_dependent, metric = "fit_time_sec", across_runs_center = "mean",
                          exclude_normalization_suspect = FALSE)

aggregate_final_benchmark(mle_gaussian_d2_independent, metric = "fit_time_sec", across_runs_center = "mean",
                          exclude_normalization_suspect = FALSE)
aggregate_final_benchmark(mle_gaussian_d3_independent, metric = "fit_time_sec", across_runs_center = "mean",
                          exclude_normalization_suspect = FALSE)
aggregate_final_benchmark(mle_gaussian_d4_independent, metric = "fit_time_sec", across_runs_center = "mean",
                          exclude_normalization_suspect = FALSE)


