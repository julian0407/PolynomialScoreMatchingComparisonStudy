source("BiasVariance_Score_Clean_patched.R")
source("Tests_Clean_patched.R")

# mle_gaussian_candidates <- readRDS("results/final_gaussian_mle_20260413-002050.rds")

# Gaussian

# Score Loss, Mean, No rdige
plot_final_benchmark(sm_gaussian_candidates_noridge, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)
aggregate_final_benchmark(sm_gaussian_candidates2_noridge, metric = "score_loss", across_runs_center = "mean",
                          exclude_normalization_suspect = FALSE)
# Average condition number per estimator, No ridge
average_condition_number_by_estimator(sm_gaussian_candidates2_noridge)

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

# Analysis of Training vs Test gap in the tails
outliers_gaussian_sm_noridge <- debug_benchmark_outliers(sm_gaussian_candidates_noridge, metric_pattern = "^score_loss$", top_n = 1)
outliers_gaussian_sm_noridge <- outliers_gaussian_sm_noridge[order(-abs(outliers_gaussian_sm_noridge$rel_dev)), ]
top3_outliers_gaussian_sm_noridge <- head(outliers_gaussian_sm_noridge, 3)
top3_outliers_gaussian_sm_noridge

gap_summary_gaussian_no_ridge  <- summarise_gap_score_compact(benchmark_obj = sm_gaussian_candidates_noridge, group_vars = c("method_label", "method", "n"), 
                                                              gap_vars = c("left_gap", "right_gap"), metric = "score_loss"
                                                              )

# Overalll
gap_summary_gaussian_no_ridge <- gap_summary_gaussian_no_ridge[
  order(-gap_summary_gaussian_no_ridge$max_minus_median_overall),
]
top3_training_test_gaps_gaussian_sm_noridge <- head(gap_summary_gaussian_no_ridge, 3)
top3_training_test_gaps_gaussian_sm_noridge

# n1000
gap_summary_gaussian_no_ridge_n1000 <- gap_summary_gaussian_no_ridge[gap_summary_gaussian_no_ridge$n == 1000,]
gap_summary_gaussian_no_ridge_n1000 <- gap_summary_gaussian_no_ridge_n1000[
  order(-gap_summary_gaussian_no_ridge_n1000$max_minus_median_overall),
]
top3_training_test_gaps_gaussian_sm_noridge_n1000 <- head(gap_summary_gaussian_no_ridge_n1000, 3)
top3_training_test_gaps_gaussian_sm_noridge_n1000

Gap_run_top1 <- replay_benchmark_run(sm_gaussian_candidates_noridge, "SM_m6_noridge_std", 1000, run_seed = 54768478)
which.max(Gap_run_top1[["pointwise"]][["score_loss_point"]])
max(Gap_run_top1[["x_test"]])
Gap_run_top1[["x_test"]][2017]
max(Gap_run_top1[["pointwise"]][["score_loss_point"]])
max(Gap_run_top1[["pointwise"]][["score_loss_point"]][-2017])
Gap_run_top1[["pointwise"]][["score_loss_point"]][2017]

# Central Score Loss, Mean, No Ridge
plot_final_benchmark(sm_gaussian_candidates_noridge, metric = "score_loss_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)

# Score Loss, Mean, Ridge
plot_final_benchmark(sm_gaussian_candidates_ridge, metric = "score_loss_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)
aggregate_final_benchmark(sm_gaussian_candidates_noridge, metric = "score_loss", across_runs_center = "mean",
                          exclude_normalization_suspect = FALSE)

average_condition_number_by_estimator(sm_gaussian_candidates_ridge)

compare_factor(
  obj_left  = sm_gaussian_candidates_ridge,
  obj_right = sm_gaussian_candidates_noridge,
  metric_left  = "score_loss",
  metric_right = "score_loss",
  exclude_normalization_suspect = TRUE
)

# Central Score Loss, Mean, Ridge
compare_factor(
  obj_left  = sm_gaussian_candidates_ridge,
  obj_right = sm_gaussian_candidates_noridge,
  metric_left  = "score_loss_central",
  metric_right = "score_loss_central",
  exclude_normalization_suspect = TRUE
)

# KL Loss, Mean, No ridge
plot_final_benchmark(sm_gaussian_candidates_noridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)
plot_final_benchmark(sm_gaussian_candidates_noridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_gaussian_candidates_noridge, metric = "kl_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)

plot_final_benchmark(sm_gaussian_candidates_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)


compare_factor(
  obj_left  = sm_gaussian_candidates_noridge,
  obj_right = sm_gaussian_candidates_ridge,
  metric_left  = "kl",
  metric_right = "kl",
  exclude_normalization_suspect = TRUE
)

# kl vs kl_central innerhalb desselben Objekts
compare_factor(
  obj_left  = sm_gaussian_candidates_ridge,
  obj_right = sm_gaussian_candidates_ridge,
  metric_left  = "kl",
  metric_right = "kl_central",
  exclude_normalization_suspect = TRUE
)

gap_summary_gaussian_no_ridge_kl  <- summarise_gap_score_compact(benchmark_obj = sm_gaussian_candidates_noridge, group_vars = c("method_label", "method", "n"), 
                                                              gap_vars = c("left_gap", "right_gap"), metric = "kl", exclude_normalization_suspect = TRUE)

tail_probability_table(
    sm_gaussian_candidates_noridge,
    estimator_pattern = "^SM_m(5|6)_",
    metric = "kl",
    thresholds = c(1, 0.5, 0.1),
    exclude_normalization_suspect = TRUE
) 








plot_final_benchmark(compare_gaussian, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)
plot_final_benchmark(compare_gaussian, metric = "score_loss_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)

plot_final_benchmark(compare_gaussian, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)

plot_final_benchmark(compare_logistic, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)
plot_final_benchmark(compare_logistic, metric = "score_loss_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)

plot_final_benchmark(compare_logistic, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)

plot_final_benchmark(compare_laplace, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)
plot_final_benchmark(compare_laplace, metric = "score_loss_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)

plot_final_benchmark(compare_laplace, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)

plot_final_benchmark(compare_gumbel, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)
plot_final_benchmark(compare_gumbel, metric = "score_loss_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)

plot_final_benchmark(compare_gumbel, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)

outliers_gaussian_compare <- debug_benchmark_outliers(compare_gaussian, metric_pattern = "^score_loss$", top_n = 3)
outliers_gaussian_compare <- outliers_gaussian_compare[order(-abs(outliers_gaussian_compare$rel_dev)), ]
top3_outliers_gaussian_compare <- head(outliers_gaussian_compare, 3)
top3_outliers_gaussian_compare

outlierrun <- replay_benchmark_run(compare_gaussian, "MLE_smoothed", 1000, run_seed = 22000300)
which.max(outlierrun[["pointwise"]][["score_loss_point"]])
min(outlierrun[["x_test"]])
outlierrun[["x_test"]][4420]
max(outlierrun[["pointwise"]][["score_loss_point"]])
max(outlierrun[["pointwise"]][["score_loss_point"]][-4420])
outlierrun[["pointwise"]][["score_loss_point"]][4420]




