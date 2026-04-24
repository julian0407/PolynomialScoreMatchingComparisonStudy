source("BiasVariance_Score_Clean_patched.R")
source("Tests_Clean_patched.R")

# mle_gaussian_candidates <- readRDS("results/final_gaussian_mle_20260413-002050.rds")


# Gaussian

plot_final_benchmark(sm_gaussian_2_candidates_noridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE)
plot_final_benchmark(sm_gaussian_2_candidates_noridge, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)

plot_final_benchmark(sm_gaussian_2_candidates_noridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)

plot_final_benchmark(sm_gaussian_2_candidates_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE)
plot_final_benchmark(sm_gaussian_2_candidates_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_gaussian_2_candidates_ridge, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE)

aggregate_final_benchmark(sm_gaussian_2_candidates_noridge, metric = "kl", across_runs_center = "mean",
                          exclude_normalization_suspect = TRUE)
aggregate_final_benchmark(sm_gaussian_2_candidates_ridge, metric = "kl", across_runs_center = "mean")

aggregate_final_benchmark(sm_gaussian_2_candidates_noridge, metric = "score_loss_central", across_runs_center = "sd")
aggregate_final_benchmark(sm_gaussian_2_candidates_ridge, metric = "score_loss", across_runs_center = "mean")


agg_n <- aggregate_final_benchmark(
  sm_gaussian_candidates_noridge,
  metric = "condition_number",
  across_runs_center = "mean"
)



outliers <- debug_benchmark_outliers(sm_gaussian_2_candidates_noridge, metric_pattern = "^kl$", top_n = 20)
outliers <- outliers[order(-abs(outliers$robust_z)), ]
outliers <- outliers[outliers$normalization_suspect == FALSE, ]
outliers

test <- replay_benchmark_run(sm_gaussian_2_candidates_noridge, "SM_m6_noridge_std", 5000, run_seed = 1342057551)
test2 <- replay_benchmark_run(sm_gaussian_2_candidates_noridge, "SM_m6_noridge_std", 5000, run_seed = 1055318384)

breaks2 = c(13,16,18, 200, 220)
table(cut(test2[["pointwise"]][["kl_point"]], breaks = breaks2))

aggregate_final_benchmark(sm_gaussian_2_candidates_noridge, metric = "kl", across_runs_center = "mean", exclude_normalization_suspect = TRUE)

test3 <- replay_benchmark_run(sm_gaussian_2_candidates_ridge, "SM_m6_ridge1e-02_std", 1000, run_seed = 54768478)

# plot_final_benchmark(res_compare_gaussian, metric = "kl", center = "median", interval = "none", log_y = TRUE)
# plot_final_benchmark(res_compare_gaussian, metric = "score_loss_central_trim", center = "median", interval = "none", log_y = TRUE)

# plot_final_benchmark(sm_gaussian_candidates_noridge, metric = "kl_central", center = "median", interval = "none", log_y = TRUE,
#                      exclude_normalization_suspect = TRUE)
# plot_final_benchmark(sm_gaussian_candidates_noridge, metric = "score_loss_central", center = "mean", interval = "none", log_y = TRUE,
#                      exclude_normalization_suspect = FALSE)
# plot_final_benchmark(sm_gaussian_candidates_noridge, metric = "score_loss", center = "sd", interval = "none", log_y = TRUE)

# plot_final_benchmark(
#   compare_gaussian_candidates,
#   metric = "fit_time_sec",
#   center = "median",
#   interval = "none",
#   log_y = TRUE
# )
# 
# plot_final_benchmark(
#   compare_gaussian_candidates,
#   metric = "density_inference_time_sec",
#   center = "median",
#   interval = "none",
#   log_y = TRUE
# )


# Logistic

plot_final_benchmark(sm_logistic_candidates_no_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE)
plot_final_benchmark(sm_logistic_candidates_no_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_logistic_candidates_no_ridge, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)

plot_final_benchmark(sm_logistic_candidates_ridge, metric = "kl", center = "median", interval = "none", log_y = TRUE)
plot_final_benchmark(sm_logistic_candidates_ridge, metric = "kl", center = "median", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_logistic_candidates_ridge, metric = "score_loss", center = "median", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = FALSE)

aggregate_final_benchmark(
  sm_logistic_candidates_ridge,
  metric = "condition_number",
  across_runs_center = "mean"
)


outliers <- debug_benchmark_outliers(sm_logistic_candidates_ridge, metric_pattern = "^score_loss$", top_n = 20)
outliers <- outliers[order(-abs(outliers$robust_z)), ]
outliers <- outliers[outliers$normalization_suspect == FALSE, ]
outliers

test3 <- replay_benchmark_run(sm_gaussian_2_candidates_noridge, "SM_m6_noridge_std", 1000, run_seed = 54768478)
test4 <- replay_benchmark_run(sm_gaussian_2_candidates_noridge, "SM_m6_noridge_std", 500, run_seed = 1359536412)

which.max(test3[["pointwise"]][["kl_point"]])
test3[["x_test"]][2017]
max(test3[["pointwise"]][["kl_point"]])
test3[["pointwise"]][["kl_point"]][2017]
max(test3[["pointwise"]][["kl_point"]][-2017])

which.max(test3[["pointwise"]][["score_loss_point"]])
test3[["x_test"]][2017]
max(test3[["pointwise"]][["score_loss_point"]])
test3[["pointwise"]][["score_loss_point"]][2017]
max(test3[["pointwise"]][["score_loss_point"]][-2017])

max(test3[["x_test"]][-2017])
max(test3[["x_test"]])
min(test3[["x_train"]])
max(test3[["x_train"]])
max(test3[["pointwise"]][["kl_point"]])
max(test3[["pointwise"]][["kl_point"]][-1047])
max(test3[["x_train"]])

which.max(test4[["pointwise"]][["kl_point"]])
test4[["x_test"]][215]
max(test4[["pointwise"]][["kl_point"]])
test4[["pointwise"]][["kl_point"]][215]

which.max(test4[["pointwise"]][["score_loss_point"]])
test4[["x_test"]][215]
max(test4[["pointwise"]][["score_loss_point"]])
test4[["pointwise"]][["score_loss_point"]][215]
max(test4[["pointwise"]][["score_loss_point"]][-215])


# Gaussian
plot_final_benchmark(sm_gaussian_2_candidates_noridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_gaussian_2_candidates_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_gaussian_2_candidates_noridge, metric = "kl_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_gaussian_2_candidates_ridge, metric = "kl_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_gaussian_2_candidates_noridge, metric = "score_loss", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_gaussian_2_candidates_noridge, metric = "score_loss_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_gaussian_2_candidates_ridge, metric = "score_loss_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)

# more
# plot_final_benchmark(sm_logistic_candidates_no_ridge, metric = "kl", center = "median", interval = "none", log_y = TRUE)

plot_final_benchmark(sm_logistic_candidates_no_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_logistic_candidates_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)

plot_final_benchmark(sm_logistic_candidates_no_ridge, metric = "kl", center = "median", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_logistic_candidates_ridge, metric = "kl", center = "median", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)

plot_final_benchmark(sm_logistic_candidates_no_ridge, metric = "kl_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_logistic_candidates_ridge, metric = "kl_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)

# plot_final_benchmark(sm_logistic_candidates_no_ridge, metric = "score_loss", center = "median", interval = "none", log_y = TRUE,
#                      exclude_normalization_suspect = FALSE)

plot_final_benchmark(sm_gumble_candidates_noridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_gumble_candidates_ridge, metric = "kl", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)

plot_final_benchmark(sm_gumble_candidates_noridge, metric = "kl", center = "median", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_gumble_candidates_ridge, metric = "kl", center = "median", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)

plot_final_benchmark(sm_gumble_candidates_noridge, metric = "kl_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
plot_final_benchmark(sm_gumble_candidates_ridge, metric = "kl_central", center = "mean", interval = "none", log_y = TRUE,
                     exclude_normalization_suspect = TRUE)
