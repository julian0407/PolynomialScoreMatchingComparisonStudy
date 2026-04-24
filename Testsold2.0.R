

aggregate_kl_without_suspect <- function(obj,
                                         metric = "kl",
                                         keep_n = NULL,
                                         drop_n = NULL,
                                         keep_method_labels = NULL,
                                         drop_method_labels = NULL,
                                         keep_methods = NULL,
                                         drop_methods = NULL,
                                         drop_all_na = FALSE) {
  if (!inherits(obj, "final_benchmark")) {
    stop("obj must be 'final_benchmark'.")
  }
  
  df <- filter_benchmark_df_by_n(obj$raw, keep_n = keep_n, drop_n = drop_n)
  df <- filter_benchmark_df_by_method(
    df,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )
  
  if (!"normalization_suspect" %in% names(df)) {
    stop("Column 'normalization_suspect' not found in obj$raw.")
  }
  
  if (!metric %in% names(df)) {
    stop(sprintf("Metric '%s' not found in obj$raw.", metric))
  }
  
  # Nur Runs behalten, die NICHT verdächtig sind
  df <- df[df$normalization_suspect %in% FALSE, , drop = FALSE]

  split_key <- interaction(df$method_label, df$n, drop = TRUE)

  agg_list <- lapply(split(df, split_key), function(dd) {
    x <- as.numeric(dd[[metric]])

    data.frame(
      method_label = dd$method_label[1],
      method = dd$method[1],
      n = dd$n[1],
      n_non_na = sum(is.finite(x)),
      mean = safe_mean(x),
      median = safe_median(x),
      sd = safe_sd(x),
      stringsAsFactors = FALSE
    )
  })

  if (length(agg_list) == 0L) return(data.frame())

  agg <- do.call(rbind, agg_list)
  rownames(agg) <- NULL

  if (isTRUE(drop_all_na)) {
    agg <- agg[agg$n_non_na > 0L, , drop = FALSE]
  }

  agg
}


aggregate_kl_without_suspect(sm_logistic_candidates_ridge, metric = "kl")
aggregate_kl_without_suspect(res_compare_gaussian, metric = "kl_central")




# m extrahieren
agg_n$m <- ifelse(
  grepl("m[0-9]+", agg_n$method_label),
  as.numeric(sub(".*m([0-9]+).*", "\\1", agg_n$method_label)),
  NA_real_
)

# jetzt über n aggregieren
agg_m <- do.call(rbind, lapply(split(agg_n, agg_n$m), function(dd) {
  x <- dd$selected
  
  data.frame(
    m = dd$m[1],
    mean_over_n = mean(x, na.rm = TRUE),
    median_over_n = median(x, na.rm = TRUE),
    sd_over_n = sd(x, na.rm = TRUE),
    n_n_values = length(unique(dd$n)),
    stringsAsFactors = FALSE
  )
}))

agg_m



# Test normalization constant

integrate_density_from_fit <- function(fit, spec, lower = -Inf, upper = Inf) {
  f <- function(x) {
    do.call(
      predict_density_estimator_generic,
      c(
        list(
          newx = x,
          fit = fit,
          family = "univariate",
          method = spec$method
        ),
        spec$density_predict_args %||% list()
      )
    )
  }
  
  integrate(
    f,
    lower = lower,
    upper = upper,
    subdivisions = 1000L,
    rel.tol = 1e-8,
    abs.tol = 1e-8,
    stop.on.error = FALSE
  )
}

mass_check <- integrate_density_from_fit(test2$fit, test2$spec)
mass_check

# nicht (-ing, inf) sondern endliches intervall
xr <- range(test$x_train, finite = TRUE)
pad <- 5 * stats::sd(test$x_train)

mass_check_finite <- integrate_density_from_fit(
  test$fit,
  test$spec,
  lower = xr[1] - pad,
  upper = xr[2] + pad
)

mass_check_finite



# repair density
repair_kl_with_optimal_constant <- function(test_obj) {
  stopifnot(!is.null(test_obj$pointwise))
  
  log_true <- test_obj$pointwise$log_true
  log_hat  <- test_obj$pointwise$log_hat
  
  keep <- is.finite(log_true) & is.finite(log_hat)
  
  log_true <- log_true[keep]
  log_hat  <- log_hat[keep]
  
  # pointwise KL
  kl_point <- log_true - log_hat
  
  # optimale Konstante
  C_opt <- mean(kl_point)
  
  # korrigierte log-density
  log_hat_corrected <- log_hat + C_opt
  
  # korrigierter KL
  kl_point_corrected <- log_true - log_hat_corrected
  
  list(
    original = list(
      mean_kl = mean(kl_point),
      median_kl = median(kl_point),
      sd_kl = sd(kl_point)
    ),
    corrected = list(
      mean_kl = mean(kl_point_corrected),
      median_kl = median(kl_point_corrected),
      sd_kl = sd(kl_point_corrected)
    ),
    C_opt = C_opt,
    improvement = mean(kl_point) - mean(kl_point_corrected)
  )
}

diag <- repair_kl_with_optimal_constant(test2)
diag



# Left right gap
summarize_gap_extremes <- function(benchmark_obj,
                                   group_vars = c("method_label", "method", "n"),
                                   gap_vars = c("left_gap", "right_gap")) {
  if (is.null(benchmark_obj$raw)) {
    stop("benchmark_obj$raw fehlt.")
  }
  
  df <- benchmark_obj$raw
  
  needed_cols <- c(group_vars, gap_vars, "run_seed", "repetition")
  missing_cols <- setdiff(needed_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Fehlende Spalten: ", paste(missing_cols, collapse = ", "))
  }
  
  out_all <- list()
  k <- 1L
  
  for (gap_var in gap_vars) {
    x <- df[, c(group_vars, gap_var, "run_seed", "repetition"), drop = FALSE]
    names(x)[names(x) == gap_var] <- "gap"
    
    grp <- split(x, interaction(x[, group_vars, drop = FALSE], drop = TRUE, lex.order = TRUE))
    
    for (d in grp) {
      valid_idx <- which(is.finite(d$gap))
      
      out <- d[1, group_vars, drop = FALSE]
      out$gap_type <- gap_var
      
      if (length(valid_idx) == 0) {
        out$mean_gap <- NA_real_
        out$max_gap <- NA_real_
        out$max_minus_mean <- NA_real_
        out$max_div_mean <- NA_real_
        out$max_gap_seed <- NA_real_
        out$max_gap_repetition <- NA_real_
      } else {
        z <- d$gap[valid_idx]
        mean_gap <- mean(z)
        local_max_idx <- which.max(z)
        row_max_idx <- valid_idx[local_max_idx]
        
        max_gap <- d$gap[row_max_idx]
        max_seed <- d$run_seed[row_max_idx]
        max_rep <- d$repetition[row_max_idx]
        
        out$mean_gap <- mean_gap
        out$max_gap <- max_gap
        out$max_minus_mean <- max_gap - mean_gap
        out$max_div_mean <- if (is.finite(mean_gap) && mean_gap != 0) max_gap / mean_gap else NA_real_
        out$max_gap_seed <- max_seed
        out$max_gap_repetition <- max_rep
      }
      
      out_all[[k]] <- out
      k <- k + 1L
    }
  }
  
  out <- do.call(rbind, out_all)
  rownames(out) <- NULL
  out
}
