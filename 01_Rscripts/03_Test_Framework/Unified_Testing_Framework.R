# ============================================================
# Unified_Testing_Framework.R
# Modular benchmark functions for univariate / multivariate tests
# Focus metrics:
#   - kl
#   - score_loss
# Optional metric variants are generated automatically from
# central_trim and robust_trim settings.
# ============================================================

source("01_Rscripts/02_Helper_and_Metrics/helper_functions.R")
source("01_Rscripts/01_Estimator/KDE.R")
source("01_Rscripts/01_Estimator/LogConcaveMLE.R")
source("01_Rscripts/01_Estimator/Univariate_Polynomial_Score_Matching.R")
source("01_Rscripts/01_Estimator/Multivariate_Pairwise_Polynomial_Score_Matching.R")
source("01_Rscripts/02_Helper_and_Metrics/Evaluation_Metrics.R")


# ---------------------------------------------------------------------------------
# (1) Helper Functions to initialize empty data frames for metric evaluation
# ---------------------------------------------------------------------------------

# Get unique metric names given metric tags and metric_args (e.g central and trim)
metric_to_output_columns <- function(metrics,
                                     density_metric_args = list(),
                                     score_metric_args = list()) {
  unique(build_metric_names(
    metrics = metrics,
    density_metric_args = density_metric_args,
    score_metric_args = score_metric_args
  ))
}

# Create named list of NA_real_ given metric names from metric_to_output_columns
make_empty_metric_list <- function(metric_columns) {
  out <- as.list(rep(NA_real_, length(metric_columns)))
  names(out) <- metric_columns
  out
}

# For multivariate Tests if "KL" loss is selected fallback for SM
resolve_metrics_for_estimator <- function(metrics, family, estimator_spec) {
  if (identical(family, "multivariate") &&
      identical(estimator_spec$method, "SM") &&
      "kl" %in% metrics) {
    return(unique(c(setdiff(metrics, "kl"), "score_loss")))
  }
  metrics
}

# Fill missing columns with NA
coerce_row_to_schema <- function(row, schema_names) {
  missing_cols <- setdiff(schema_names, names(row))
  if (length(missing_cols) > 0L) {
    for (nm in missing_cols) row[[nm]] <- NA
  }
  row[, schema_names, drop = FALSE]
}

# Create dataframe with common scheme based on rows with different named columns or missing columns
bind_rows_fill_base <- function(rows) {
  if (length(rows) == 0L) return(data.frame())
  # create common scheme
  schema <- unique(unlist(lapply(rows, names), use.names = FALSE))
  # Fill missing columns according common scheme
  rows2 <- lapply(rows, coerce_row_to_schema, schema_names = schema)
  # bind these rows 
  out <- do.call(rbind, rows2)
  # reset row names
  rownames(out) <- NULL
  out
}


# ---------------------------------------------------------------------------------
# (2) Safety function to evaluate all requested metrics and safe them into named df
# ---------------------------------------------------------------------------------

# This function tries to evaluate all requested metric "together" such that 
# runtime does not increase substantially
# Important arguments
#   - metric, family and method tags
#   - test data
#   - true density, log density and score function
#   - predict_args for score based and kl based metrics (central and trim)
#   - verbose to control console outputs
safe_evaluate_requested_metrics <- function(metrics,
                                            x_test,
                                            fit,
                                            family,
                                            method,
                                            true_density = NULL,
                                            true_logdensity = NULL,
                                            true_score = NULL,
                                            density_predict_args = list(),
                                            score_predict_args = list(),
                                            density_metric_args = list(),
                                            score_metric_args = list(),
                                            estimator_label = method,
                                            n = NA_integer_,
                                            verbose = TRUE) {
  # initialize empty named list given metric tags
  metric_columns <- metric_to_output_columns(
    metrics,
    density_metric_args = density_metric_args,
    score_metric_args = score_metric_args
  )
  empty_out <- make_empty_metric_list(metric_columns)
  
  # Compute requested metrics at once
  full_res <- tryCatch(
    evaluate_requested_metrics(
      metrics = metrics,
      x_test = x_test,
      fit = fit,
      family = family,
      method = method,
      true_logdensity = true_logdensity,
      true_score = true_score,
      density_predict_args = density_predict_args,
      score_predict_args = score_predict_args,
      density_metric_args = density_metric_args,
      score_metric_args = score_metric_args
    ),
    error = function(e) e
  )
  
  # No error -> Take evaluated metrics and save them in empty out -> return
  if (!inherits(full_res, "error")) {
    for (nm in names(full_res)) empty_out[[nm]] <- full_res[[nm]]
    return(empty_out)
  }

  # Warning message that fallback to evaluate metrics one after another
  if (isTRUE(verbose)) {
    warning(sprintf(
      "Joint metric evaluation failed for %s at n=%s. Falling back to metric-wise evaluation. Reason: %s",
      estimator_label, n, conditionMessage(full_res)
    ))
  }
  
  
  for (m in metrics) {
    # Evaluation of a single metric
    one_res <- tryCatch(
      evaluate_requested_metrics(
        metrics = m,
        x_test = x_test,
        fit = fit,
        family = family,
        method = method,
        true_logdensity = true_logdensity,
        true_score = true_score,
        density_predict_args = density_predict_args,
        score_predict_args = score_predict_args,
        density_metric_args = density_metric_args,
        score_metric_args = score_metric_args
      ),
      error = function(e) e
    )
    # Warning message if this metric fails, jump to next evaluation
    if (inherits(one_res, "error")) {
      if (isTRUE(verbose)) {
        warning(sprintf(
          "Metric '%s' failed for %s at n=%s: %s",
          m, estimator_label, n, conditionMessage(one_res)
        ))
      }
      next
    }
    # if succesfull save these values zo empty_out
    for (nm in names(one_res)) empty_out[[nm]] <- one_res[[nm]]
  }
  # return
  empty_out
}


# ---------------------------------------------------------------------------------
# (3) One single Benchmark run
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# (3.1) Seed helper for local lambda-mu starts
# ---------------------------------------------------------------------------------

# Generates new run_seed with offset, offset=1 -> run_seed does not change
make_lambda_mu_seed_from_run_seed <- function(run_seed, offset = 1L) {
  # check if run_seed is valid input argument
  if (is.null(run_seed) || length(run_seed) != 1L || !is.finite(run_seed)) {
    return(NULL)
  }
  # Upper bound for valid run_seed
  modulus <- as.double(.Machine$integer.max - 1L)
  # generates run seed wirh offset
  as.integer(((as.double(run_seed) + as.double(offset) - 1) %% modulus) + 1L)
}

# Helper function to generate run seeds such that every estimator uses same seeds in training and testing
make_paired_run_seed_grid <- function(sample_sizes, n_rep, seed = NULL) {
  # set global seed
  if (!is.null(seed)) set.seed(seed)
  
  # empty seed list
  run_seed_grid <- list()
  # generates for each sample size n_rep random seeds
  for (n in sample_sizes) {
    run_seed_grid[[as.character(n)]] <- sample.int(.Machine$integer.max, size = n_rep)
  }
  # return list where each sample size "n" -> c(seed_1, ..., seed_n_rep)
  run_seed_grid
}

# This function performs one single benchmark run
# Important arguments
#   - metric and family tags
#   - sample function of true data (r_sample)
#   - training size n
#   - test size n_test
#   - true density, log density and score function
#   - Estimator object with
#         - predict_args for score based and kl based metrics (central and trim)
#         - method tags
#         - smoothed for log concave MLE (set this to false for other estimators)
#   - run_seed to make run repeatable for outlier analysis
#   - verbose to control console outputs
run_one_final_experiment <- function(n,
                                     family = c("univariate", "multivariate"),
                                     estimator_spec,
                                     r_sample,
                                     metrics,
                                     metric_columns = NULL,
                                     n_test = 2000,
                                     true_density = NULL,
                                     true_logdensity = NULL,
                                     true_score = NULL,
                                     run_seed = NULL,
                                     verbose = TRUE) {
  family <- match.arg(family)
  # fallback option for multivariate tests
  metrics <- resolve_metrics_for_estimator(
    metrics = metrics,
    family = family,
    estimator_spec = estimator_spec
  )
  
  # initalize named list of requested metricy
  metric_columns <- metric_columns %||% metric_to_output_columns(
    metrics,
    density_metric_args = estimator_spec$density_metric_args %||% list(),
    score_metric_args = estimator_spec$score_metric_args %||% list()
  )
  
  # Set run seed, preprocessing of data into correct format
  if (!is.null(run_seed)) set.seed(run_seed)
  x_train <- r_sample(n)
  if (family == "univariate") x_train <- as.numeric(x_train) else x_train <- as.matrix(x_train)
  
  # Optional: derive a local lambda_mu_seed from run_seed.
  # This is only active for new estimator_specs that explicitly set lambda_mu_seed_from_run_seed = TRUE.
  
  # Only for SM estimator
  if (identical(family, "univariate") &&
      identical(estimator_spec$method, "SM")) {
    
    # Get args or initialize them with standard values
    estimator_spec$fit_args <- estimator_spec$fit_args %||% list()
    parameterization <- estimator_spec$fit_args$parameterization %||% "psd"
    # Only generate seed if direct lambda and mu optimization is used, lambda_mu_seed_from_run_seed==True
    # and lambda_mu_seed not NULL
    if (identical(parameterization, "lambda_mu") &&
        isTRUE(estimator_spec$lambda_mu_seed_from_run_seed) &&
        is.null(estimator_spec$fit_args$lambda_mu_seed)) {
      
      # Get offset or set offset to one
      lambda_mu_seed_offset <- estimator_spec$lambda_mu_seed_offset %||% 1L
      # Generate lambda_mu_seed with offset based on run_seed
      estimator_spec$fit_args$lambda_mu_seed <-
        make_lambda_mu_seed_from_run_seed(
          run_seed = run_seed,
          offset = lambda_mu_seed_offset
        )
    }
  }
  
  # Save runtime for fit
  fit_time <- system.time({
    # Call generic fit function
    fit <- tryCatch(
      do.call(
        fit_estimator_generic,
        c(
          list(
            x = x_train,
            family = family,
            method = estimator_spec$method,
            smoothed = estimator_spec$smoothed
          ),
          # optional arguments
          estimator_spec$fit_args
        )
      ),
      error = function(e) structure(list(error_message = conditionMessage(e)), class = "fit_error")
    )
  })
  
  # initialize basic result row
  base_row <- data.frame(
    n = n,
    run_seed = if (is.null(run_seed)) NA_integer_ else as.integer(run_seed),
    method_label = estimator_spec$label,
    method = estimator_spec$method,
    fit_time_sec = as.numeric(fit_time["elapsed"]),
    density_inference_time_sec = NA_real_,
    score_inference_time_sec = NA_real_,
    total_inference_time_sec = NA_real_,
    success = FALSE,
    status = NA_character_,
    iterations = NA_real_,
    objective_value = NA_real_,
    condition_number = NA_real_,
    normalization_ok = NA,
    normalization_suspect = NA,
    normalization_suspect_0_1 = NA,
    normalization_suspect_1_0 = NA,
    normalization_loghat_finite_share = NA_real_,
    normalization_median_kl_shift = NA_real_,
    stringsAsFactors = FALSE,
    lc_min_eigenvalue = NA_real_,
    lc_min_eigenvalue_raw = NA_real_,
    lc_max_min_eigenvalue = NA_real_,
    lc_mean_min_eigenvalue = NA_real_,
    lc_n_grid_points = NA_real_,
    lc_n_violated = NA_real_,
    lc_max_violation = NA_real_,
    lc_mean_violation = NA_real_,
    lc_is_log_concave_tol0 = NA,
    lc_min_eigenvalues = I(list(numeric(0)))
  )
  for (nm in metric_columns) base_row[[nm]] <- NA_real_

  # If error return this base row
  if (inherits(fit, "fit_error")) {
    base_row$status <- fit$error_message
    return(base_row)
  }
  
  # else we perform testing
  # generate test sample 
  x_test <- r_sample(n_test)
  if (family == "univariate") x_test <- as.numeric(x_test) else x_test <- as.matrix(x_test)
  
  # get some diagnostics if data univariate
  if (family == "univariate") {
    train_min <- min(x_train, na.rm = TRUE)
    train_max <- max(x_train, na.rm = TRUE)
    test_min  <- min(x_test, na.rm = TRUE)
    test_max  <- max(x_test, na.rm = TRUE)
    
    base_row$left_gap  <- max(0, train_min - test_min)
    base_row$right_gap <- max(0, test_max - train_max)
  }
  
  # create diagnostic object from fit
  diags <- tryCatch(
    extract_fit_diagnostics(
      fit,
      x_diag = if (identical(family, "multivariate")) x_test else NULL
    ),
    error = function(e) list(
      success = TRUE,
      status = paste("diagnostics_error:", conditionMessage(e)),
      iterations = NA_real_,
      objective_value = NA_real_,
      condition_number = NA_real_
    )
  )

  # which metrics are needed?
  density_metrics_requested <- any(metrics %in% c("kl"))
  score_metrics_requested <- any(metrics %in% c("score_loss"))

  # save runtime for determine logdensity of this run
  density_time <- if (density_metrics_requested) {
    system.time({
      invisible(do.call(
        predict_logdensity_estimator_generic,
        c(
          list(newx = x_test, fit = fit, family = family, method = estimator_spec$method),
          estimator_spec$density_predict_args %||% list()
        )
      ))
    })
  } else c(elapsed = NA_real_)

  # save runtime for determine score of this run
  score_time <- if (score_metrics_requested) {
    system.time({
      invisible(do.call(
        predict_score_estimator_generic,
        c(
          list(newx = x_test, fit = fit, family = family, method = estimator_spec$method),
          estimator_spec$score_predict_args %||% list()
        )
      ))
    })
  } else c(elapsed = NA_real_)

  # save runtimes and computeb total runtime
  base_row$density_inference_time_sec <- as.numeric(density_time["elapsed"])
  base_row$score_inference_time_sec <- as.numeric(score_time["elapsed"])
  base_row$total_inference_time_sec <- sum(
    c(base_row$density_inference_time_sec, base_row$score_inference_time_sec),
    na.rm = TRUE
  )
  if (!is.finite(base_row$total_inference_time_sec)) base_row$total_inference_time_sec <- NA_real_

  # Check if normalization was succesfull
  if (density_metrics_requested) {
    dens_diag <- tryCatch(
      extract_density_diagnostic(
        fit = fit,
        family = family,
        method = estimator_spec$method,
        predict_args = estimator_spec$density_predict_args %||% list()
      ),
      error = function(e) list(normalization_ok = FALSE, normalization_message = conditionMessage(e))
    )
    
    base_row$normalization_ok <- dens_diag$normalization_ok
    
    # More important: Check if normalization seems reasonable for SM
    base_row$normalization_suspect <- NA
    base_row$normalization_loghat_finite_share <- NA_real_
    base_row$normalization_median_kl_shift <- NA_real_
    
    if (identical(family, "univariate") &&
        identical(estimator_spec$method, "SM") &&
        !is.null(true_logdensity)) {
      
      # Get log density
      log_hat <- tryCatch(
        do.call(
          predict_logdensity_estimator_generic,
          c(
            list(newx = x_test, fit = fit, family = family, method = estimator_spec$method),
            estimator_spec$density_predict_args %||% list()
          )
        ),
        error = function(e) rep(NA_real_, length(x_test))
      )
      
      log_hat <- as.numeric(log_hat)
      log_true <- as.numeric(true_logdensity(x_test))
      
      # get finite share of estimates
      finite_share <- mean(is.finite(log_hat))
      kl_point <- log_true - log_hat
      # get median kl_point
      median_kl_shift <- if (any(is.finite(kl_point))) {
        stats::median(kl_point[is.finite(kl_point)])
      } else {
        Inf
      }
      
      base_row$normalization_loghat_finite_share <- finite_share
      base_row$normalization_median_kl_shift <- median_kl_shift
      
      # get mean kl and sd
      kl_point <- kl_point[is.finite(kl_point)]
      mean_kl <- mean(kl_point)
      sd_kl   <- sd(kl_point)
      
      # Check if constant shift of KL "almost everywhere" with mean kl > 0.1
      constant_shift_flag_0_1 <-
        length(kl_point) > 10 &&
        is.finite(mean_kl) &&
        is.finite(sd_kl) &&
        abs(mean_kl) > 0.1 &&
        (sd_kl / abs(mean_kl)) < 0.1
      
      # Check if constant shift of KL "almost everywhere"  with mean kl > 1.0
      constant_shift_flag_1_0 <-
        length(kl_point) > 10 &&
        is.finite(mean_kl) &&
        is.finite(sd_kl) &&
        abs(mean_kl) > 1.0 &&
        (sd_kl / abs(mean_kl)) < 0.1
      
      # normalization is suspect if at least one of this diagnostics is true
      base_row$normalization_suspect_0_1 <-
        isFALSE(dens_diag$normalization_ok) ||
        (!is.na(finite_share) && finite_share < 0.99) ||
        constant_shift_flag_0_1
      
      base_row$normalization_suspect_1_0 <-
        isFALSE(dens_diag$normalization_ok) ||
        (!is.na(finite_share) && finite_share < 0.99) ||
        constant_shift_flag_1_0
      
      # Standard diagnostic is with mean kl > 0.1
      base_row$normalization_suspect <- base_row$normalization_suspect_0_1
    }
  }
  
  # Now evaluate finally the requested metrics
  metric_values <- safe_evaluate_requested_metrics(
    metrics = metrics,
    x_test = x_test,
    fit = fit,
    family = family,
    method = estimator_spec$method,
    true_logdensity = true_logdensity,
    true_score = true_score,
    density_predict_args = estimator_spec$density_predict_args %||% list(),
    score_predict_args = estimator_spec$score_predict_args %||% list(),
    density_metric_args = estimator_spec$density_metric_args %||% list(),
    score_metric_args = estimator_spec$score_metric_args %||% list(),
    estimator_label = estimator_spec$label,
    n = n,
    verbose = verbose
  )
  
  # save diagnostics from fit
  base_row$success <- diags$success %||% TRUE
  base_row$status <- diags$status %||% NA_character_
  base_row$iterations <- diags$iterations %||% NA_real_
  base_row$objective_value <- diags$objective_value %||% NA_real_
  base_row$condition_number <- diags$condition_number %||% NA_real_
  base_row$lc_min_eigenvalue <- diags$lc_min_eigenvalue %||% NA_real_
  base_row$lc_min_eigenvalue_raw <- diags$lc_min_eigenvalue_raw %||% NA_real_
  base_row$lc_max_min_eigenvalue <- diags$lc_max_min_eigenvalue %||% NA_real_
  base_row$lc_mean_min_eigenvalue <- diags$lc_mean_min_eigenvalue %||% NA_real_
  base_row$lc_n_grid_points <- diags$lc_n_grid_points %||% NA_real_
  base_row$lc_n_violated <- diags$lc_n_violated %||% NA_real_
  base_row$lc_max_violation <- diags$lc_max_violation %||% NA_real_
  base_row$lc_mean_violation <- diags$lc_mean_violation %||% NA_real_
  
  base_row$lc_is_log_concave_tol0 <- if (is.finite(base_row$lc_min_eigenvalue)) {
    base_row$lc_min_eigenvalue >= 0
  } else {
    NA
  }
  
  base_row$lc_min_eigenvalues <- I(list(diags$lc_min_eigenvalues %||% numeric(0)))
  # save evaluated metrics
  for (nm in names(metric_values)) base_row[[nm]] <- metric_values[[nm]]
  base_row
}


# -------------------------------------------------------------------------------------------------
# (4) Run a family of estimators to create a data object that can be easily used to create plots
#     and interpret results
# -------------------------------------------------------------------------------------------------

# This function performs multiple benchmarks for multiple estimators run
# Important arguments
#   - most input as for run_one_final_experiment is needed (see above)
#   - estimator_specs is now a list of dfferent estimators
#   - sample sizes as list for training
#   - n_rep is number of repetiotions (fit + test) to evaluate the average score+Kl loss over multiple runs
#   - save, save_dir, and save_name specify if the data should be saved in the directory

run_final_benchmark <- function(sample_sizes,
                                family = c("univariate", "multivariate"),
                                estimator_specs,
                                r_sample,
                                metrics,
                                n_rep = 20,
                                n_test = 2000,
                                true_density = NULL,
                                true_logdensity = NULL,
                                true_score = NULL,
                                truth_name = NULL,
                                seed = NULL,
                                verbose = TRUE,
                                save = FALSE,
                                save_dir = ".",
                                save_name = NULL) {
  family <- match.arg(family)
  # set seed if desired
  if (!is.null(seed)) set.seed(seed)

  
  out <- list()
  counter <- 1L
  # get desired metrics
  metric_columns <- unique(unlist(lapply(estimator_specs, function(spec) {
    metrics_eff <- resolve_metrics_for_estimator(
      metrics = metrics,
      family = family,
      estimator_spec = spec
    )
    
    metric_to_output_columns(
      metrics_eff,
      density_metric_args = spec$density_metric_args %||% list(),
      score_metric_args = spec$score_metric_args %||% list()
    )
  }), use.names = FALSE))
  
  # iterate over estimators
  for (spec in estimator_specs) {
    if (verbose) message("Method: ", spec$label)
    # iterate over requested sample sizes for training
    for (n in sample_sizes) {
      if (verbose) message("  n = ", n)
      # perform n_rep benchmark runs for each of these configurations
      for (rep in seq_len(n_rep)) {
        if (verbose) message("    repetition ", rep, "/", n_rep)
        # per repetition create own seed for outlier analysis
        run_seed <- sample.int(.Machine$integer.max, size = 1L)
        # compute final benchmark
        ans <- run_one_final_experiment(
          n = n,
          family = family,
          estimator_spec = spec,
          r_sample = r_sample,
          metrics = metrics,
          metric_columns = metric_columns,
          n_test = n_test,
          true_density = true_density,
          true_logdensity = true_logdensity,
          true_score = true_score,
          run_seed = run_seed,
          verbose = verbose
        )
        # flag run
        ans$repetition <- rep
        # save in output
        out[[counter]] <- ans
        counter <- counter + 1L
      }
    }
  }
  
  # create table
  raw <- bind_rows_fill_base(out)
  # reorder
  preferred_order <- c(
    "n", "repetition", "run_seed", "method_label", "method",
    "fit_time_sec", "density_inference_time_sec", "score_inference_time_sec", "total_inference_time_sec",
    "success", "status", "iterations", "objective_value", "condition_number",
    "left_gap", "right_gap",
    "normalization_ok", "normalization_suspect",
    "normalization_suspect_0_1", "normalization_suspect_1_0",
    "normalization_loghat_finite_share", "normalization_median_kl_shift",
    metric_columns
  )
  other_cols <- setdiff(names(raw), preferred_order)
  raw <- raw[, c(intersect(preferred_order, names(raw)), other_cols), drop = FALSE]
  rownames(raw) <- NULL
  
  # create final benchmark object
  obj <- structure(
    list(
      raw = raw,
      settings = list(
        sample_sizes = sample_sizes,
        family = family,
        truth_name = truth_name,
        metrics = metrics,
        metric_columns = metric_columns,
        n_rep = n_rep,
        n_test = n_test,
        benchmark_seed = seed
      ),
      estimator_specs = estimator_specs,
      benchmark_inputs = list(
        r_sample = r_sample,
        true_density = true_density,
        true_logdensity = true_logdensity,
        true_score = true_score
      )
    ),
    class = "final_benchmark"
  )
  
  # optinally save object as .rds 
  if (isTRUE(save)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    save_name <- save_name %||% make_benchmark_filename(
      kind = "final",
      truth_name = truth_name,
      estimator_specs = estimator_specs,
      family = family
    )
    save_path <- file.path(save_dir, save_name)
    saveRDS(obj, save_path)
    obj$settings$save_path <- save_path
  }

  obj
}

# ---------------------------------------------------------------------------------
# (5) Paired comparison of univariate SM parameterizations
# ---------------------------------------------------------------------------------
# This code is to compare the SM estiator using the PSD constraint with the direct 
# lambda - mu optimization

# ------------------------------------------------------------
# (5.1) Diagnostics and Helper for comparing SM fits
# ------------------------------------------------------------

# Get coefficients from G that correspond to second derivative of polynomial s
second_derivative_coefficients_from_G <- function(G) {
  # q(z) = phi(z)^T G phi(z) with max degree 2m-2
  m <- nrow(G)
  # Create empty coefficient vector
  coeff <- numeric(2L * m - 1L)
  
  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      # degree = i+j-2
      deg <- (i - 1L) + (j - 1L)
      # sum up coefficients of G that belongs to basis z^deg
      coeff[deg + 1L] <- coeff[deg + 1L] + G[i, j]
    }
  }
  coeff
}

# Derive score coefficients from G and c1
score_coefficients_from_G_c1 <- function(G, c1) {
  # Get coefficients from secon derivative
  q <- second_derivative_coefficients_from_G(G)
  # Initialize empty vector
  out <- numeric(length(q) + 1L)
  # first entry is the constant term c1
  out[1L] <- c1
  
  # Derive coefficents after integration
  for (r in seq_along(q)) {
    deg <- r - 1L
    out[deg + 2L] <- q[r] / (deg + 1L)
  }
  out
}

# Compute l2 difference of coefficients of score polynomials obtained by two fit objects
coef_l2_diff <- function(fit_a, fit_b) {
  ca <- score_coefficients_from_G_c1(fit_a$G, fit_a$c1)
  cb <- score_coefficients_from_G_c1(fit_b$G, fit_b$c1)
  
  # Fill missing coefficoents with 0
  L <- max(length(ca), length(cb))
  ca <- c(ca, rep(0, L - length(ca)))
  cb <- c(cb, rep(0, L - length(cb)))
  # l2 difference of coefficients
  sum((ca - cb)^2)
}

# l2 difference of matrices G obtained by two fit objects
gram_l2_diff <- function(fit_a, fit_b) {
  sum((as.vector(fit_a$G) - as.vector(fit_b$G))^2)
}

# Get objective value after fitting
empirical_objective_value <- function(fit,
                                      z_train = fit$z_train) {
  # Get degree parameter m
  m <- fit$m
  # This function is only used for not simple SM
  h = function(z) rep(1, length(z))
  h_prime = function(z) rep(0, length(z))
  
  # Build input of objecctive function
  inp <- build_K_and_l_univariate(z_train, m, h, h_prime)
  K <- inp$K + fit$ridge * diag(nrow(inp$K))
  l <- inp$l
  
  # Use scaling to backtransform to internal scaling
  g_sc <- as.vector(fit$G) * inp$scale_vec
  y <- c(g_sc, fit$c1)
  # Calculate objective value
  as.numeric(0.5 * crossprod(y, K %*% y) - sum(l * y))
}

# Compute score loss based on fit and test data
score_loss_univariate_fit <- function(fit, x_test, true_score) {
  # Use metric evaluation from other script
  res <- evaluate_requested_metrics(
    metrics = "score_loss",
    x_test = x_test,
    fit = fit,
    family = "univariate",
    method = "SM",
    true_score = true_score
  )
  # only return score loss
  as.numeric(res$score_loss)
}

# Helper function to measure fitting time and determine score loss
fit_score_time_sm_univariate <- function(parameterization,
                                         x_train,
                                         x_test,
                                         m,
                                         ridge,
                                         true_score,
                                         standardize = TRUE,
                                         optim_control = list(maxit = 1000, reltol = 1e-8),
                                         lambda_mu_seed = NULL,
                                         lambda_mu_n_starts = 10,
                                         lambda_mu_init_sd = 0.1) {
  # Call fitting and save time needed for fit
  fit_time <- system.time({
    fit <- fit_score_matching_univariate(
      x = x_train,
      m = m,
      standardize = standardize,
      ridge = ridge,
      parameterization = parameterization,
      optim_control = optim_control,
      lambda_mu_seed = lambda_mu_seed,
      lambda_mu_n_starts = lambda_mu_n_starts,
      lambda_mu_init_sd = lambda_mu_init_sd
    )
  })
  # Call score_loss calculation and save time needed for it
  score_time <- system.time({
    score_loss <- score_loss_univariate_fit(
      fit = fit,
      x_test = x_test,
      true_score = true_score
    )
  })
  
  #return fit, times, score loss, objective value, and solver status
  list(
    fit = fit,
    fit_time = as.numeric(fit_time["elapsed"]),
    score_time = as.numeric(score_time["elapsed"]),
    score_loss = score_loss,
    objective = empirical_objective_value(fit),
    status = fit$status
  )
}

# ------------------------------------------------------------
# (5.2) Run Comparison experiment
# ------------------------------------------------------------

run_one_sm_parameterization_comparison <- function(n,
                                                   m,
                                                   ridge = 0,
                                                   n_test = 3000,
                                                   run_seed,
                                                   r_sample,
                                                   true_score,
                                                   standardize = TRUE,
                                                   optim_control = list(maxit = 1000, reltol = 1e-8),
                                                   lambda_mu_n_starts = 10,
                                                   lambda_mu_init_sd = 0.1) {
  # set global seed and sample training and test data
  set.seed(run_seed)
  x_train <- as.numeric(r_sample(n))
  x_test  <- as.numeric(r_sample(n_test))
  # do not change lambda mu seed
  lambda_mu_seed <- make_lambda_mu_seed_from_run_seed(
    run_seed = run_seed,
    offset = 1L
  )
  # Fitting + measuring for psd
  psd <- fit_score_time_sm_univariate(
    parameterization = "psd",
    x_train = x_train,
    x_test = x_test,
    m = m,
    ridge = ridge,
    true_score = true_score,
    standardize = standardize,
    optim_control = optim_control,
  )
  # Fitting + measuring for lambda mu
  lambda_mu <- fit_score_time_sm_univariate(
    parameterization = "lambda_mu",
    x_train = x_train,
    x_test = x_test,
    m = m,
    ridge = ridge,
    true_score = true_score,
    standardize = standardize,
    optim_control = optim_control,
    lambda_mu_seed = lambda_mu_seed,
    lambda_mu_n_starts = lambda_mu_n_starts,
    lambda_mu_init_sd = lambda_mu_init_sd
  )
  
  # Save data into dataframe
  data.frame(
    n = n,
    m = m,
    ridge = ridge,
    run_seed = run_seed,
    lambda_mu_seed = lambda_mu_seed,
    
    fit_time_psd = psd$fit_time,
    fit_time_lambda_mu = lambda_mu$fit_time,
    
    score_inference_time_psd = psd$score_time,
    score_inference_time_lambda_mu = lambda_mu$score_time,

    score_loss_psd = psd$score_loss,
    score_loss_lambda_mu = lambda_mu$score_loss,
    
    objective_psd = psd$objective,
    objective_lambda_mu = lambda_mu$objective,
    
    coef_l2_psd_vs_lambda_mu =
      coef_l2_diff(psd$fit, lambda_mu$fit),
    
    gram_l2_psd_vs_lambda_mu =
      gram_l2_diff(psd$fit, lambda_mu$fit),

    
    status_psd = psd$status,
    status_lambda_mu = lambda_mu$status,
    
    stringsAsFactors = FALSE
  )
}

# Use many single runs to compare the average differences in coefficients and other metrics
run_sm_parameterization_comparison <- function(sample_sizes = c(200, 500, 1000),
                                               m_values = c(5, 6),
                                               ridge = 0,
                                               n_rep = 10,
                                               n_test = 3000,
                                               seed = 123,
                                               r_sample,
                                               true_score,
                                               standardize = TRUE,
                                               optim_control = list(maxit = 1000, reltol = 1e-8),
                                               lambda_mu_n_starts = 10,
                                               lambda_mu_init_sd = 0.1,
                                               save = FALSE,
                                               save_dir = ".",
                                               save_name = NULL) {
  # Generate seed grid to use the same training and test data for both estimators (psd and lambda mu)
  run_seed_grid <- make_paired_run_seed_grid(
    sample_sizes = sample_sizes,
    n_rep = n_rep,
    seed = seed
  )
  
  # Initialize empty result list and counter
  rows <- list()
  counter <- 1L
  
  # Iterate over all configurations
  for (m in m_values) {
    for (n in sample_sizes) {
      for (rep in seq_len(n_rep)) {
        # Get corresponding run seed
        run_seed <- run_seed_grid[[as.character(n)]][rep]
        
        message(
          "SM parameterization comparison: m=", m,
          ", n=", n,
          ", rep=", rep, "/", n_rep
        )
        # run one comparison_bemchmark_run
        rows[[counter]] <- tryCatch(
          run_one_sm_parameterization_comparison(
            n = n,
            m = m,
            ridge = ridge,
            n_test = n_test,
            run_seed = run_seed,
            r_sample = r_sample,
            true_score = true_score,
            standardize = standardize,
            optim_control = optim_control,
            lambda_mu_n_starts = lambda_mu_n_starts,
            lambda_mu_init_sd = lambda_mu_init_sd
          ),
          error = function(e) {
            data.frame(
              n = n,
              m = m,
              ridge = ridge,
              run_seed = run_seed,
              error = conditionMessage(e),
              stringsAsFactors = FALSE
            )
          }
        )
        # iterate counter
        rows[[counter]]$repetition <- rep
        counter <- counter + 1L
      }
    }
  }
  # Summarize values in one dataframe
  raw <- bind_rows_fill_base(rows)
  # Initlize list of names that will be filled with average values
  summary_vars <- c(
    "fit_time_psd",
    "fit_time_lambda_mu",
    "score_inference_time_psd",
    "score_inference_time_lambda_mu",
    "score_loss_psd",
    "score_loss_lambda_mu",
    "objective_psd",
    "objective_lambda_mu",
    "coef_l2_psd_vs_lambda_mu",
    "gram_l2_psd_vs_lambda_mu"
  )
  # Only take variables that exist in raw
  summary_vars <- intersect(summary_vars, names(raw))
  # aggregate and group by n, m and ridge/noridge
  summary <- aggregate(
    raw[, summary_vars, drop = FALSE],
    by = list(n = raw$n, m = raw$m, ridge = raw$ridge),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  # Create final structure
  out <- structure(
    list(
      raw = raw,
      summary = summary,
      settings = list(
        sample_sizes = sample_sizes,
        m_values = m_values,
        ridge = ridge,
        n_rep = n_rep,
        n_test = n_test,
        seed = seed,
        lambda_mu_n_starts = lambda_mu_n_starts,
        lambda_mu_init_sd = lambda_mu_init_sd
      )
    ),
    class = "sm_parameterization_comparison"
  )
  
  # Optional: Save object 
  if (isTRUE(save)) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    }
    save_name <- save_name %||% paste0(
      "sm_parameterization_comparison_",
      ".rds"
    )
    save_path <- file.path(save_dir, save_name)
    saveRDS(out, save_path)
  }
  out
}

# --------------------------------------------------------------------------
# (6) Helper functions to filter result object to create specific plots
# --------------------------------------------------------------------------

# Helper function to check if input arguments are valid to filter (either state methods 
# to keep or to drop)
normalize_method_filter <- function(keep_method_labels = NULL,
                                    drop_method_labels = NULL,
                                    keep_methods = NULL,
                                    drop_methods = NULL) {
  if (!is.null(keep_method_labels) && !is.null(drop_method_labels)) {
    stop("Use either keep_method_labels or drop_method_labels, not both.")
  }
  if (!is.null(keep_methods) && !is.null(drop_methods)) {
    stop("Use either keep_methods or drop_methods, not both.")
  }
  list(
    keep_method_labels = if (is.null(keep_method_labels)) NULL else unique(as.character(keep_method_labels)),
    drop_method_labels = if (is.null(drop_method_labels)) NULL else unique(as.character(drop_method_labels)),
    keep_methods = if (is.null(keep_methods)) NULL else unique(as.character(keep_methods)),
    drop_methods = if (is.null(drop_methods)) NULL else unique(as.character(drop_methods))
  )
}

# Provide dataframe and filter for/out method labels/names
filter_benchmark_df_by_method <- function(df,
                                          keep_method_labels = NULL,
                                          drop_method_labels = NULL,
                                          keep_methods = NULL,
                                          drop_methods = NULL) {
  filt <- normalize_method_filter(
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )
  out <- df
  # apply filtering if argument is not null
  if (!is.null(filt$keep_method_labels)) out <- out[out$method_label %in% filt$keep_method_labels, , drop = FALSE]
  if (!is.null(filt$drop_method_labels)) out <- out[!out$method_label %in% filt$drop_method_labels, , drop = FALSE]
  if (!is.null(filt$keep_methods)) out <- out[out$method %in% filt$keep_methods, , drop = FALSE]
  if (!is.null(filt$drop_methods)) out <- out[!out$method %in% filt$drop_methods, , drop = FALSE]
  rownames(out) <- NULL
  # return df without rownames
  out
}

# analogous to normalize_method_filter but for the sample sizes for training
normalize_n_filter <- function(keep_n = NULL, drop_n = NULL) {
  if (!is.null(keep_n) && !is.null(drop_n)) stop("Use either keep_n or drop_n, not both.")
  list(
    keep_n = if (is.null(keep_n)) NULL else unique(as.numeric(keep_n)),
    drop_n = if (is.null(drop_n)) NULL else unique(as.numeric(drop_n))
  )
}

# analogous to filter_benchmark_df_by_method but for the sample sizes for training
filter_benchmark_df_by_n <- function(df, keep_n = NULL, drop_n = NULL) {
  filt <- normalize_n_filter(keep_n = keep_n, drop_n = drop_n)
  out <- df
  if (!is.null(filt$keep_n)) out <- out[out$n %in% filt$keep_n, , drop = FALSE]
  if (!is.null(filt$drop_n)) out <- out[!out$n %in% filt$drop_n, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# --------------------------------------------------------------------------
# (7) Use the raw data in the result object to create aggregated benchmarks 
#     (average score/Kl loss)
# --------------------------------------------------------------------------
# Important arguments
#   - obj -> benchmark object
#   - across_runs_center -> which aggregated metric should be used
#   - metric -> which metric should be analyzed (kl, score_loss, central etc)
#   - exclude_normalization_suspect = FALSE (important fopr density metrics in univariate SM)

aggregate_final_benchmark <- function(obj,
                                      metric,
                                      drop_all_na = FALSE,
                                      keep_n = NULL,
                                      drop_n = NULL,
                                      keep_method_labels = NULL,
                                      drop_method_labels = NULL,
                                      keep_methods = NULL,
                                      drop_methods = NULL,
                                      across_runs_center = c("mean", "median", "sd"),
                                      exclude_normalization_suspect = FALSE,
                                      conf_level = 0.95) {
  across_runs_center <- match.arg(across_runs_center)
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")
  
  # Filter benchmark object
  df <- filter_benchmark_df_by_n(obj$raw, keep_n = keep_n, drop_n = drop_n)
  df <- filter_benchmark_df_by_method(
    df,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )
  
  # optionally exclude suspect runs
  if (isTRUE(exclude_normalization_suspect)) {
    if (!"normalization_suspect" %in% names(df)) {
      stop("Column 'normalization_suspect' not found in benchmark output.")
    }
    df <- df[!(df$method %in% "SM" & df$normalization_suspect %in% TRUE), , drop = FALSE]
  }
  
  # check if desired metric is in benchmark object
  if (!metric %in% names(df)) {
    available <- intersect(c(
      obj$settings$metric_columns,
      "fit_time_sec", "density_inference_time_sec", "score_inference_time_sec", "total_inference_time_sec",
      "condition_number"
    ), names(df))
    stop(sprintf(
      "Metric '%s' not found in benchmark output. Available metrics: %s",
      metric,
      paste(available, collapse = ", ")
    ))
  }
  
  if (nrow(df) == 0L) {
    return(data.frame())
  }
  
  # grouping key for aggregates
  split_key <- interaction(df$method_label, df$n, drop = TRUE)
  
  # perform aggregations according to split key and filtered data df
  agg_list <- lapply(split(df, split_key), function(dd) {
    x <- as.numeric(dd[[metric]])
    x_finite <- x[is.finite(x)]
    
    # outlier share
    q1 <- safe_quantile(x, 0.25)
    q3 <- safe_quantile(x, 0.75)
    iqr <- if (is.finite(q1) && is.finite(q3)) q3 - q1 else NA_real_
    
    outlier_share <- if (length(x_finite) == 0L || !is.finite(iqr)) {
      NA_real_
    } else if (iqr == 0) {
      mean(x_finite != stats::median(x_finite))
    } else {
      lo <- q1 - 1.5 * iqr
      hi <- q3 + 1.5 * iqr
      mean(x_finite < lo | x_finite > hi)
    }
    
    # fit failure and normalization failure rates
    failure_rate <- mean((dd$success %in% FALSE) | is.na(dd$success) | !is.finite(x))
    
    normalization_failure_rate <- if ("normalization_ok" %in% names(dd)) {
      mean(dd$normalization_ok %in% FALSE, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    # Normalization suspect rate for univariate SM
    normalization_suspect_rate <- if ("normalization_suspect" %in% names(dd)) {
      mean(dd$normalization_suspect %in% TRUE, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    # Get data length (number of benchmark runs inside a group)
    R <- length(x_finite)
    
    # Get mean, sd and se
    mean_x <- if (R > 0L) mean(x_finite) else NA_real_
    sd_x   <- if (R > 1L) stats::sd(x_finite) else NA_real_
    se_x   <- if (R > 1L) sd_x / sqrt(R) else NA_real_
    
    # Calculate t-quantile for given conf_level and number of benchmark_runs - 1
    alpha <- 1 - conf_level
    tcrit <- if (R > 1L) stats::qt(1 - alpha / 2, df = R - 1L) else NA_real_
    
    # Calculate CI bounds
    ci_low  <- mean_x - tcrit * se_x
    ci_high <- mean_x + tcrit * se_x
    
    # several standard aggregates and information
    data.frame(
      method_label = dd$method_label[1],
      method = dd$method[1],
      n = dd$n[1],
      n_non_na = sum(is.finite(x)),
      mean = safe_mean(x),
      median = safe_median(x),
      selected = switch(
        across_runs_center,
        mean = safe_mean(x),
        median = safe_median(x),
        sd = safe_sd_finite(x)
      ),
      q25 = q1,
      q75 = q3,
      iqr = iqr,
      sd = safe_sd_finite(x),
      across_runs_center = across_runs_center,
      outlier_run_share = outlier_share,
      failure_rate = failure_rate,
      normalization_failure_rate = normalization_failure_rate,
      normalization_suspect_rate = normalization_suspect_rate,
      stringsAsFactors = FALSE,
      n_rep_used = R,
      se = se_x,
      ci_low = ci_low,
      ci_high = ci_high,
      conf_level = conf_level
    )
  })
  
  # combine rows in df
  agg <- do.call(rbind, agg_list)
  rownames(agg) <- NULL
  
  if (isTRUE(drop_all_na)) {
    agg <- agg[agg$n_non_na > 0L, , drop = FALSE]
  }
  
  agg
}

# Filtering not only raw data, but provide an object and return this object with filtered
# data according to method label and sample size filter
subset_final_benchmark <- function(obj,
                                   keep_n = NULL,
                                   drop_n = NULL,
                                   keep_method_labels = NULL,
                                   drop_method_labels = NULL,
                                   keep_methods = NULL,
                                   drop_methods = NULL,
                                   update_settings = TRUE) {
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")
  out <- obj
  out$raw <- filter_benchmark_df_by_n(out$raw, keep_n = keep_n, drop_n = drop_n)
  out$raw <- filter_benchmark_df_by_method(
    out$raw,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )
  if (isTRUE(update_settings)) out$settings$sample_sizes <- sort(unique(out$raw$n))
  out
}

# --------------------------------------------------------------------------
# (9) Debugging outlier runs
# --------------------------------------------------------------------------

# Function to derive the outlier runs that contributes comparatively the most to the median value of a metric
# - provide benchmark object and metrics to identify outliers
# - topn is number of topn largest runs according to these metrics per group (methodlabel, samplesize)
debug_benchmark_outliers <- function(obj,
                                     metric_pattern = "^(kl|score_loss)",
                                     top_n = 10) {
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")
  # grouping according methid label and sample size + minimum runs that mus exists
  group_cols = c("method_label", "n")
  min_group_size = 3L
  # Get raw data
  df <- obj$raw
  # get all metrics og Kl or score loss
  metric_cols <- grep(metric_pattern, names(df), value = TRUE)
  if (length(metric_cols) == 0L) stop("No metric columns matched metric_pattern.")

  rows <- list()
  counter <- 1L
  for (metric in metric_cols) {
    # iterate over data according to groups (method label, sample size)
    split_key <- interaction(df[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
    for (dd in split(df, split_key)) {
      x <- as.numeric(dd[[metric]])
      keep <- is.finite(x)
      if (sum(keep) < min_group_size) next
      center <- stats::median(x[keep])
      # mad_val = median(abs(x - median(x))) "interpret as robust variance"
      mad_val <- stats::mad(x[keep], center = center, constant = 1, na.rm = TRUE)
      scale_val <- max(mad_val, 1e-12)
      abs_dev <- abs(x - center)
      # "z score with median metrics"
      robust_z <- abs_dev / scale_val
      # More relevant: (Size of difference in terms of size of median -> "relative discrepancy")
      rel_dev <- abs_dev / pmax(abs(center), 1e-12)
      # Order according relative discrepancy
      ord <- order(rel_dev, decreasing = TRUE, na.last = NA)
      # Take only topn values
      take <- head(ord, top_n)
      # Save all diagnsoticy for topn ids
      for (idx in take) {
        rows[[counter]] <- data.frame(
          metric = metric,
          method_label = dd$method_label[idx],
          method = dd$method[idx],
          n = dd$n[idx],
          repetition = dd$repetition[idx],
          run_seed = dd$run_seed[idx],
          condition_number = dd$condition_number[idx],
          normalization_suspect = dd$normalization_suspect[idx],
          value = x[idx],
          group_center = center,
          abs_dev = abs_dev[idx],
          rel_dev = rel_dev[idx],
          robust_z = robust_z[idx],
          stringsAsFactors = FALSE
        )
        counter <- counter + 1L
      }
    }
  }

  out <- bind_rows_fill_base(rows)
  if (nrow(out) == 0L) return(out)
  # return sorted outliers with relative biggest outlier at the Top
  out <- out[order(out$rel_dev, out$abs_dev, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Function to replay a specific Benchmark run to analyze a specific run in detail (only for univariate data)
# Replay one benchmark run and receive pointiwse metrics given a specific seed, sample size, method label, and benchmark object
replay_benchmark_run <- function(obj,
                                     run_seed,
                                     method_label = NULL,
                                     n = NULL) {
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")
  
  df <- obj$raw
  
  # Check if such a combination exists (seed, method label, sample size)
  sel <- df$run_seed == run_seed
  if (!is.null(method_label)) sel <- sel & df$method_label == method_label
  if (!is.null(n)) sel <- sel & df$n == n
  hits <- df[sel, , drop = FALSE]
  if (nrow(hits) == 0L) stop("No matching run found.")
  if (nrow(hits) > 1L) stop("Selection is ambiguous. Add method_label or n.")
  # Take the row of this run given the raw data in the benchmark object
  row <- hits[1, , drop = FALSE]
  
  # Lookup estimator specification
  spec_idx <- which(vapply(
    obj$estimator_specs,
    function(spec) identical(spec$label, row$method_label),
    logical(1)
  ))
  if (length(spec_idx) != 1L) stop("Could not uniquely resolve estimator_spec.")
  spec <- obj$estimator_specs[[spec_idx]]
  
  # Get settings from benchgmark object
  metrics <- obj$settings$metrics
  n_test <- obj$settings$n_test
  family <- obj$settings$family
  # set seed
  set.seed(row$run_seed)
  # Create training sample with this seed
  x_train <- obj$benchmark_inputs$r_sample(row$n)
  if (family == "univariate") {
    x_train <- as.numeric(x_train)
  } else {
    x_train <- as.matrix(x_train)
  }
  # Call generic fit function
  fit <- do.call(
    fit_estimator_generic,
    c(
      list(
        x = x_train,
        family = family,
        method = spec$method,
        smoothed = spec$smoothed
      ),
      spec$fit_args
    )
  )
  # Create test sample (seed is still set)
  x_test <- obj$benchmark_inputs$r_sample(n_test)
  if (family == "univariate") {
    x_test <- as.numeric(x_test)
  } else {
    x_test <- as.matrix(x_test)
  }
  # Get aggregated metrics
  aggregate_metrics <- evaluate_requested_metrics(
    metrics = metrics,
    x_test = x_test,
    fit = fit,
    family = family,
    method = spec$method,
    true_logdensity = obj$benchmark_inputs$true_logdensity,
    true_score = obj$benchmark_inputs$true_score,
    density_predict_args = spec$density_predict_args %||% list(),
    score_predict_args = spec$score_predict_args %||% list(),
    density_metric_args = spec$density_metric_args %||% list(),
    score_metric_args = spec$score_metric_args %||% list()
  )
  # Get pointwise estimates and metrics
  if (family == "univariate") {
    pointwise <- data.frame(x_test = as.numeric(x_test))
  } else {
    x_mat <- as.matrix(x_test)
    pointwise <- as.data.frame(x_mat)
    names(pointwise) <- paste0("x", seq_len(ncol(x_mat)))
  }
  if ("kl" %in% metrics) {
    kl_pt <- compute_pointwise_kl(
      x_test = x_test,
      fit = fit,
      family = family,
      method = spec$method,
      true_logdensity = obj$benchmark_inputs$true_logdensity,
      predict_args = spec$density_predict_args %||% list()
    )
    
    pointwise$kl_point <- kl_pt$losses
    pointwise$log_true <- kl_pt$log_true
    pointwise$log_hat <- kl_pt$log_hat
  }
  if ("score_loss" %in% metrics && !is.null(obj$benchmark_inputs$true_score)) {
    score_pt <- compute_pointwise_score_loss(
      x_test = x_test,
      fit = fit,
      family = family,
      method = spec$method,
      true_score = obj$benchmark_inputs$true_score,
      predict_args = spec$score_predict_args %||% list()
    )
    
    score_loss <- rep(NA_real_, nrow(pointwise))
    score_loss[score_pt$keep] <- score_pt$losses
    pointwise$score_loss_point <- score_loss
  }
  # return pointwise and aggregated metrics
  list(
    selected_run = row,
    estimator_spec = spec,
    fit = fit,
    x_train = x_train,
    x_test = x_test,
    aggregate_metrics = aggregate_metrics,
    pointwise = pointwise
  )
}


# --------------------------------------------------------------------------
# (10) Crete plots based on raw data by applying aggregate function
# --------------------------------------------------------------------------
plot_final_benchmark <- function(obj,
                                 metric,
                                 center = c("mean", "median", "sd"),
                                 interval = c("iqr", "sd", "ci", "none"),
                                 interval_geom = c("ribbon", "errorbar", "linerange"),
                                 conf_level = 0.95,
                                 keep_n = NULL,
                                 drop_n = NULL,
                                 keep_method_labels = NULL,
                                 drop_method_labels = NULL,
                                 keep_methods = NULL,
                                 drop_methods = NULL,
                                 method_label_map = NULL,
                                 drop_all_na = TRUE,
                                 facet_methods = FALSE,
                                 log_x = FALSE,
                                 log_y = FALSE,
                                 exclude_normalization_suspect = FALSE) {
  center <- match.arg(center)
  interval <- match.arg(interval)
  interval_geom <- match.arg(interval_geom)
  # get aggregated benchmark object
  agg <- aggregate_final_benchmark(
    obj = obj,
    metric = metric,
    drop_all_na = drop_all_na,
    keep_n = keep_n,
    drop_n = drop_n,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods,
    across_runs_center = center,
    exclude_normalization_suspect = exclude_normalization_suspect,
    conf_level = conf_level
  )
  
  if (nrow(agg) == 0L) stop("No rows left after applying the sample-size filter.")
  
  # Rename method labels according dictionary
  if (!is.null(method_label_map)) {
    hit <- agg$method_label %in% names(method_label_map)
    agg$method_label[hit] <- unname(method_label_map[agg$method_label[hit]])
  }
  # Build title from data in readable name
  metric_label <- tools::toTitleCase(gsub("_", " ", metric))
  center_label <- if (center == "mean") {
    "Average"
  } else if (center == "median") {
    "Median"
  } else {
    "SD"
  }
  # Replace Score Loss with Score MSE
  metric_label <- gsub("Score Loss", "Score MSE", metric_label)
  # Replace KL with EKL
  metric_label <- gsub("Kl", "EKL", metric_label)
  axis_label <- paste(center_label, metric_label)
  plot_title <- paste(axis_label, "by Sample Size")
  # get data
  agg$y <- agg[[center]]
  # optional iqr band
  if (interval == "iqr") {
    agg$ymin <- agg$q25
    agg$ymax <- agg$q75
  } else if (interval == "sd") {
    agg$ymin <- agg$mean - agg$sd
    agg$ymax <- agg$mean + agg$sd
  } else if (interval == "ci") {
    agg$ymin <- agg$ci_low
    agg$ymax <- agg$ci_high
  } else {
    agg$ymin <- NA_real_
    agg$ymax <- NA_real_
  }
  
  # Apply logscale if True
  if (isTRUE(log_y) && interval != "none") {
    vals <- c(agg$y, agg$ymin, agg$ymax)
    vals_pos <- vals[is.finite(vals) & vals > 0]
    
    if (length(vals_pos) > 0L) {
      eps <- min(vals_pos, na.rm = TRUE) / 10
      agg$ymin[!is.finite(agg$ymin) | agg$ymin <= 0] <- eps
      agg$ymax[!is.finite(agg$ymax) | agg$ymax <= 0] <- NA_real_
    }
  }
  # subtitle if KL
  subtitle <- NULL
  if (grepl("kl", metric, ignore.case = TRUE) & exclude_normalization_suspect) {
    subtitle <- sprintf(
      "Exclude normalization suspect: %s",
      exclude_normalization_suspect
    )
  }
  # build plot
  p <- ggplot2::ggplot(
    agg,
    ggplot2::aes(x = n, y = y, color = method_label, group = method_label)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Sample Size n",
      y = axis_label,
      color = "Method",
      title = plot_title,
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
  # Optionally Add iqr, ci, sd band
  if (interval != "none") {
    if (interval_geom == "ribbon") {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ymin, ymax = ymax, fill = method_label),
        alpha = 0.15,
        colour = NA,
        show.legend = FALSE
      )
    } else if (interval_geom == "errorbar") {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = ymin, ymax = ymax),
        width = 0.04,
        alpha = 0.75,
        show.legend = FALSE
      )
    } else if (interval_geom == "linerange") {
      p <- p + ggplot2::geom_linerange(
        ggplot2::aes(ymin = ymin, ymax = ymax),
        alpha = 0.75,
        show.legend = FALSE
      )
    }
  }
  # optional logscale
  p <- apply_optional_log_scale(p, agg$n, axis = "x", requested = log_x)
  p <- apply_optional_log_scale(
    p,
    c(agg$y, agg$ymin, agg$ymax),
    axis = "y",
    requested = log_y
  )
  # Optional: Split plots
  if (isTRUE(facet_methods)) {
    p <- p +
      ggplot2::scale_x_continuous(
        breaks = c(0, max(agg$n, na.rm = TRUE) / 2, max(agg$n, na.rm = TRUE))
      ) +
      ggplot2::facet_wrap(~ method_label)
  }
  p
}


# --------------------------------------------------------------------------
# (10) Additional analysis tools that were needed during our analysis
# --------------------------------------------------------------------------

# Condition number
# Get average condition number across all runs (independent of sample size) per estimator
average_condition_number_by_estimator <- function(obj) {
  df <- obj$raw
  
  split_list <- split(df, df$method_label, drop = TRUE)
  
  out <- lapply(split_list, function(dd) {
    x <- as.numeric(dd$condition_number)
    x <- x[is.finite(x)]
    
    data.frame(
      method_label = dd$method_label[1],
      average_condition_number = mean(x),
      n_runs_used = length(x)
    )
  })
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$average_condition_number), ]
}


# Mean / Median Ratio
# - given estimator labels, metric , benchmark object
score_loss_mean_median_factor <- function(
    obj,
    estimator_pattern = "^SM_m(5|6)_",
    metric = "score_loss",
    exclude_normalization_suspect = FALSE,
    finite_only = TRUE
) {
  if (!inherits(obj, "final_benchmark")) {
    stop("obj must be 'final_benchmark'.")
  }
  
  df <- obj$raw
  
  if (!metric %in% names(df)) {
    stop(sprintf("Metric '%s' not found in obj$raw.", metric))
  }
  
  # get the data that contains a method label according to the provided pattern
  df <- df[grepl(estimator_pattern, df$method_label), , drop = FALSE]
  
  # if normalizingsuspect = true in provided input filter only for rows that contain a false value here
  if (isTRUE(exclude_normalization_suspect)) {
    if (!"normalization_suspect" %in% names(df)) {
      stop("Column 'normalization_suspect' not found in benchmark output.")
    }
    df <- df[df$normalization_suspect %in% FALSE, , drop = FALSE]
  }
  
  if (nrow(df) == 0L) {
    return(data.frame())
  }
  # group by method label and sample size
  split_key <- interaction(df$method_label, df$n, drop = TRUE)
  # Compute per group the mean, median, and the ratio (some other metrics are also included but not used at the moment)
  out <- lapply(split(df, split_key), function(dd) {
    x <- as.numeric(dd[[metric]])
    
    if (isTRUE(finite_only)) {
      x <- x[is.finite(x)]
    }
    
    med <- if (length(x) == 0L) NA_real_ else stats::median(x)
    mn  <- if (length(x) == 0L) NA_real_ else mean(x)
    
    data.frame(
      method_label = dd$method_label[1],
      method = if ("method" %in% names(dd)) dd$method[1] else NA_character_,
      n = dd$n[1],
      n_used = length(x),
      median = med,
      mean = mn,
      mean_div_median = if (is.finite(mn) && is.finite(med) && med != 0) mn / med else NA_real_,
      median_div_mean = if (is.finite(mn) && is.finite(med) && mn != 0) med / mn else NA_real_,
      log2_mean_div_median = if (is.finite(mn) && is.finite(med) && mn > 0 && med > 0) log2(mn / med) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out <- out[order(out$method_label, out$n), , drop = FALSE]
  out
}


# get proportion of extreme runs per method label and sample size group
tail_probability_table <- function(
    obj,
    estimator_pattern = "^SM_m(5|6)_",
    metric = "score_loss",
    thresholds = c(1, 0.1, 0.05),
    exclude_normalization_suspect = FALSE
) {
  # Filter for rows with method label according to provided pattern
  df <- obj$raw
  df <- df[grepl(estimator_pattern, df$method_label), , drop = FALSE]
  # Optionally Filter out suspicious runs (normalizing cosntant is suspect)
  if (isTRUE(exclude_normalization_suspect)) {
    df <- df[df$normalization_suspect %in% FALSE, , drop = FALSE]
  }
  # group per method label and sample size and compute the proportion that is greater than thresholds
  out <- do.call(rbind, lapply(thresholds, function(eps) {
    tmp <- aggregate(df[[metric]] > eps,
                     by = list(method_label = df$method_label, n = df$n),
                     FUN = function(z) mean(z, na.rm = TRUE))
    names(tmp)[3] <- "prob"
    tmp$threshold <- eps
    tmp
  }))
  
  out[order(out$method_label, out$threshold, out$n), ]
}

# Needed for the score loss analysis considering dependence of perfromance on training data
# This function computes losses by aggregating runs that uses test samples that are within the training interval
# Additionally the max right and left gaps to the training data is reported and set to zero if largest/smallest test 
# data is within training ingterval
summarise_gap_score_compact <- function(benchmark_obj,
                                        group_vars = c("method_label", "method", "n"),
                                        gap_vars = c("left_gap", "right_gap"),
                                        metric = "score_loss",
                                        eps = 1e-12) {
  df <- benchmark_obj$raw
  # Check if all needed columns are availbale in benchmark object
  needed_cols <- unique(c(group_vars, gap_vars, metric, "run_seed"))
  missing_cols <- setdiff(needed_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  # get raw data according to grouping
  groups <- split(
    df,
    interaction(df[, group_vars, drop = FALSE], drop = TRUE, lex.order = TRUE)
  )
  # compute for each group the gap anylsis
  out <- lapply(groups, function(d) {
    res <- d[1, group_vars, drop = FALSE]
    
    # initialize variables for max gap
    best_gap_type <- NA_character_
    best_max_gap <- NA_real_
    best_max_minus_median <- NA_real_
    best_run_seed <- NA
    
    for (gap_var in gap_vars) {
      z_raw <- d[[gap_var]]
      valid_idx <- which(is.finite(z_raw))
      
      if (length(valid_idx) > 0L) {
        # Compute gap data
        z <- pmax(z_raw[valid_idx], 0)
        z_median <- median(z)
        local_max_idx <- which.max(z)
        row_idx <- valid_idx[local_max_idx]
        z_max <- z[local_max_idx]
        # Check if max gap is greater than the best_max_gap so far
        # and overwrite in this case (or if best_max_gap is not intliazed)
        if (is.na(best_max_gap) || z_max > best_max_gap) {
          best_gap_type <- gap_var
          best_max_gap <- z_max
          best_max_minus_median <- z_max - z_median
          best_run_seed <- d$run_seed[row_idx]
        }
      }
    }
    # Save largest identified gap for each group + related data
    res$gap_type <- best_gap_type
    res$max_gap <- best_max_gap
    res$max_minus_median_overall <- best_max_minus_median
    res$run_seed <- best_run_seed
    # Filter runs of this group whose test data lies within training interval
    both_zero <- is.finite(d$left_gap) & abs(d$left_gap) <= eps &
      is.finite(d$right_gap) & abs(d$right_gap) <= eps
    # Save the metric values of runs satisfy  and not satisfy both_zero in dofferent vectors
    x0 <- as.numeric(d[[metric]][both_zero])
    x1 <- as.numeric(d[[metric]][!both_zero])
    
    x0 <- x0[is.finite(x0)]
    x1 <- x1[is.finite(x1)]
    
    # get number and prooportion of runs that have test sample within training interval
    res$n_both_zero <- sum(both_zero, na.rm = TRUE)
    res$prop_both_zero <- mean(both_zero, na.rm = TRUE)
    
    # calculate mean/median of this vectors x0, x1
    res$mean_both_zero <- if (length(x0) == 0L) NA_real_ else mean(x0)
    res$median_both_zero <- if (length(x0) == 0L) NA_real_ else median(x0)
    
    res$mean_at_least_one_not_zero <- if (length(x1) == 0L) NA_real_ else mean(x1)
    res$median_at_least_one_not_zero <- if (length(x1) == 0L) NA_real_ else median(x1)
    
    # Compute ratios of mean and median to compare difference across estimators
    res$ratio_m0_m1 <- if (is.na(res$mean_both_zero) || is.na(res$mean_at_least_one_not_zero) ||
                           abs(res$mean_at_least_one_not_zero) <= eps) {
      NA_real_
    } else {
      res$mean_both_zero / res$mean_at_least_one_not_zero
    }
    
    res$ratio_m1_m0 <- if (is.na(res$mean_both_zero) || is.na(res$mean_at_least_one_not_zero) ||
                           abs(res$mean_both_zero) <= eps) {
      NA_real_
    } else {
      res$mean_at_least_one_not_zero / res$mean_both_zero
    }
    
    res
  })
  # return final table
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

# Calculate Factor of Decrease / Increase of ridge
# provide two benchmark objects and compute the ratios of specific metrics of them
compare_factor <- function(obj_left,
                           obj_right,
                           metric_left = "score_loss",
                           metric_right = metric_left,
                           exclude_normalization_suspect = FALSE) {
  
  # helper function: Get only degree parameter from method label
  extract_m <- function(x) {
    suppressWarnings(as.integer(sub(".*_m([0-9]+).*", "\\1", x)))
  }
  # get aggregated benchmark objects with desired metric
  agg_left <- aggregate_final_benchmark(
    obj_left,
    metric = metric_left,
    exclude_normalization_suspect = exclude_normalization_suspect
  )
  agg_right <- aggregate_final_benchmark(
    obj_right,
    metric = metric_right,
    exclude_normalization_suspect = exclude_normalization_suspect
  )
  # Save degree parameter
  agg_left$m  <- extract_m(agg_left$method_label)
  agg_right$m <- extract_m(agg_right$method_label)
  # drop unnecessary columns
  left  <- agg_left[,  c("m", "n", "mean")]
  right <- agg_right[, c("m", "n", "method_label", "mean")]
  # name columns
  names(left)  <- c("m", "n", "mean_left")
  names(right) <- c("m", "n", "method_label_right", "mean_right")
  # merge seperate df by group (degree, samplesize)
  out <- merge(right, left, by = c("m", "n"))
  
  # Compute ratio
  factor_numeric <- out$mean_left / out$mean_right
  out$Enhance_Performance <- factor_numeric > 1
  # Update fnumber ormat 
  out$mean_right <- formatC(out$mean_right, format = "e", digits = 3)
  out$mean_left  <- formatC(out$mean_left,  format = "e", digits = 3)
  out$factor     <- formatC(factor_numeric, format = "e", digits = 3)
  # sort columns
  out <- out[, c(
    "m", "n", "method_label_right",
    "mean_right", "mean_left",
    "factor", "Enhance_Performance"
  )]
  
  out <- out[order(out$m, out$n), ]
  rownames(out) <- NULL
  out
}

# Generate Boxplot for a metric in a final benchmark object
plot_final_benchmark_boxplot <- function(obj,
                                         metric,
                                         keep_n = NULL,
                                         drop_n = NULL,
                                         keep_method_labels = NULL,
                                         drop_method_labels = NULL,
                                         keep_methods = NULL,
                                         drop_methods = NULL,
                                         method_label_map = NULL,
                                         log_y = FALSE,
                                         exclude_normalization_suspect = FALSE) {
  # Check if final benchmark object
  if (!inherits(obj, "final_benchmark")) stop("obj must be 'final_benchmark'.")
  
  # Get filtered raw data
  df <- filter_benchmark_df_by_n(obj$raw, keep_n = keep_n, drop_n = drop_n)
  df <- filter_benchmark_df_by_method(
    df,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )
  
  # optionally exclude runs with suspicious normlaizing constant
  if (isTRUE(exclude_normalization_suspect)) {
    if (!"normalization_suspect" %in% names(df)) {
      stop("Column 'normalization_suspect' not found in benchmark output.")
    }
    df <- df[!(df$method %in% "SM" & df$normalization_suspect %in% TRUE), , drop = FALSE]
  }
  
  # Only metrics that are contained in df can be plotted
  if (!metric %in% names(df)) stop(sprintf("Metric '%s' not found.", metric))
  
  # Get desired finite metric values
  df$value <- as.numeric(df[[metric]])
  df <- df[is.finite(df$value), , drop = FALSE]
  
  # Apply renaming for method labels in the plots
  if (!is.null(method_label_map)) {
    hit <- df$method_label %in% names(method_label_map)
    df$method_label[hit] <- unname(method_label_map[df$method_label[hit]])
  }
  
  # x labels should be categorical values
  df$n <- factor(df$n, levels = sort(unique(df$n)))
  
  # Generate axis names
  metric_label <- tools::toTitleCase(gsub("_", " ", metric))
  metric_label <- gsub("Score Loss", "Score MSE", metric_label)
  metric_label <- gsub("Kl", "EKL", metric_label)
  
  # generate boxplot (0.25, 0.75 quartile, 1.5 x IQR = Whisker)
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = n, y = value, fill = method_label)
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0.35) +
    ggplot2::labs(
      x = "Sample Size n",
      y = metric_label,
      fill = "Method",
      title = paste(metric_label, "distribution across repetitions")
    ) +
    ggplot2::theme_minimal()
  
  # Optional: Apply log scale
  p <- apply_optional_log_scale(p, df$value, axis = "y", requested = log_y)
  # return plot
  p
}


# --------------------------------------------------------------------------
# (11) Additional analysis for Confidence Intervals
# --------------------------------------------------------------------------

# Create a summary table that only reports the CI related values and not the full aggregated data
ci_summary_table <- function(obj,
                             metric = "score_loss",
                             conf_level = 0.95,
                             method_label_map = NULL,
                             ...) {
  # Call aggregate_final_benchmark() using mean
  tab <- aggregate_final_benchmark(
    obj = obj,
    metric = metric,
    across_runs_center = "mean",
    conf_level = conf_level,
    ...
  )
  
  # Apply optional renaming for plot descriptions in later functions
  if (!is.null(method_label_map)) {
    hit <- tab$method_label %in% names(method_label_map)
    tab$method_label[hit] <- unname(method_label_map[tab$method_label[hit]])
  }
  
  # Check if lower interval is negative and if both bounds are positive
  tab$ci_low_negative <- tab$ci_low < 0
  tab$log_scale_ci_ok <- tab$ci_low > 0 & tab$ci_high > 0
  
  # Filtere aggregated data only for relevant CI data
  tab <- tab[, c(
    "method_label",
    "method",
    "n",
    "n_rep_used",
    "mean",
    "ci_low",
    "ci_high",
    "ci_low_negative",
    "log_scale_ci_ok"
  )]
  
  # Return sorted table group by (method_label, training size)
  tab[order(tab$method_label, tab$n), ]
}


# Plots the relative halfwidth of a confidence interval
# Is used to interpret CIs with negative lower bound
# needs CI summary table as input and apply optional logscale
# Plots descriotion is only for Score MSE. This function is not for other metrics
plot_relative_ci_halfwidth <- function(ci_tab, log_y = TRUE) {
  # Computes relative halfwidth
  ci_tab$relative_ci_halfwidth <- (ci_tab$ci_high - ci_tab$mean) / ci_tab$mean
  
  # Plots values (recall: method renaming was already generated in CI summary table)
  p <- ggplot2::ggplot(
    ci_tab,
    ggplot2::aes(
      x = n,
      y = relative_ci_halfwidth,
      color = method_label,
      group = method_label
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Sample Size n",
      y = "Relative CI half-width",
      color = "Method",
      title = "Relative Monte Carlo Uncertainty of Mean Score MSE"
    ) +
    ggplot2::theme_minimal()
  
  # Apply optional logscale
  if (isTRUE(log_y)) {
    p <- p + ggplot2::scale_y_log10()
  }
  
  p
}


# Alternative Transformed CI Intervals with delta method. Hardcoded for Score MSE. For other metrics this method must be adapted
plot_score_mse_log_ci <- function(obj,
                                  conf_level = 0.95,
                                  method_label_map = NULL) {
  # Takes and checks final benchmark object
  if (!inherits(obj, "final_benchmark")) {
    stop("obj must be a 'final_benchmark' object.")
  }
  
  # Splits data according group key (method, training size)
  groups <- split(
    obj$raw,
    interaction(obj$raw$method_label, obj$raw$n, drop = TRUE)
  )
  
  # Apply this function for every grouped data (Computes the CI with delta method)
  ci <- do.call(rbind, lapply(groups, function(d) {
    # Get finite score loss
    x <- as.numeric(d$score_loss)
    x <- x[is.finite(x)]
    # get data vector length
    R <- length(x)
    
    # Must at least two Score MSE values
    if (R < 2L) return(NULL)
    
    # Compute mean
    xbar <- mean(x)
    
    # t-quantile for (1+ conf_level)/2 and df=R-1
    h <- stats::qt((1 + conf_level) / 2, df = R - 1L) *
      # adapted sd part because of delta method with log
      stats::sd(x) / (xbar * sqrt(R))
    
    # Backtransform CI and save final data
    data.frame(
      method_label = d$method_label[1L],
      n = d$n[1L],
      mean = xbar,
      ci_low = xbar * exp(-h),
      ci_high = xbar * exp(h)
    )
  }))
  
  # Apply renaming of methods
  if (!is.null(method_label_map)) {
    hit <- ci$method_label %in% names(method_label_map)
    ci$method_label[hit] <-
      unname(method_label_map[ci$method_label[hit]])
  }
  
  # Reorder data according (method, training size)
  ci <- ci[order(ci$method_label, ci$n), ]
  
  # Generate Plot with CI ribbon
  ggplot2::ggplot(
    ci,
    ggplot2::aes(
      x = n,
      y = mean,
      color = method_label,
      group = method_label
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = ci_low,
        ymax = ci_high,
        fill = method_label
      ),
      alpha = 0.15,
      colour = NA,
      show.legend = FALSE
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    # apply logscale
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Sample Size n",
      y = "Average Score MSE",
      color = "Method",
      title = "Average Score MSE with Back-Transformed Log CI"
    ) +
    ggplot2::theme_minimal()
}
