# Additional Tests

# Condition number
average_condition_number_by_estimator <- function(
    obj,
    keep_n = NULL,
    drop_n = NULL,
    keep_method_labels = NULL,
    drop_method_labels = NULL,
    keep_methods = NULL,
    drop_methods = NULL,
    finite_only = TRUE,
    only_success = FALSE,
    exclude_normalization_suspect = FALSE,
    decreasing = FALSE
) {
  if (!inherits(obj, "final_benchmark")) {
    stop("obj must be 'final_benchmark'.")
  }
  
  df <- obj$raw
  
  # optional dieselben Filter wie im restlichen Benchmark-Code
  df <- filter_benchmark_df_by_n(df, keep_n = keep_n, drop_n = drop_n)
  df <- filter_benchmark_df_by_method(
    df,
    keep_method_labels = keep_method_labels,
    drop_method_labels = drop_method_labels,
    keep_methods = keep_methods,
    drop_methods = drop_methods
  )
  
  if (isTRUE(exclude_normalization_suspect)) {
    if (!"normalization_suspect" %in% names(df)) {
      stop("Column 'normalization_suspect' not found in benchmark output.")
    }
    df <- df[df$normalization_suspect %in% FALSE, , drop = FALSE]
  }
  
  if (!"condition_number" %in% names(df)) {
    stop("Column 'condition_number' not found in obj$raw.")
  }
  
  if (isTRUE(only_success) && "success" %in% names(df)) {
    df <- df[df$success %in% TRUE, , drop = FALSE]
  }
  
  if (nrow(df) == 0L) {
    return(data.frame())
  }
  
  split_list <- split(df, df$method_label, drop = TRUE)
  
  out <- lapply(split_list, function(dd) {
    x <- as.numeric(dd$condition_number)
    
    if (isTRUE(finite_only)) {
      x <- x[is.finite(x)]
    }
    
    data.frame(
      method_label = dd$method_label[1],
      method = if ("method" %in% names(dd)) dd$method[1] else NA_character_,
      average_condition_number = if (length(x) == 0L) NA_real_ else mean(x),
      n_runs_used = length(x),
      n_rows_total = nrow(dd),
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out <- out[order(out$average_condition_number, decreasing = decreasing), , drop = FALSE]
  out
}


# Mean / Median Ratio
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
  
  # nur SM5 und SM6
  df <- df[grepl(estimator_pattern, df$method_label), , drop = FALSE]
  
  if (isTRUE(exclude_normalization_suspect)) {
    if (!"normalization_suspect" %in% names(df)) {
      stop("Column 'normalization_suspect' not found in benchmark output.")
    }
    df <- df[df$normalization_suspect %in% FALSE, , drop = FALSE]
  }
  
  if (nrow(df) == 0L) {
    return(data.frame())
  }
  
  split_key <- interaction(df$method_label, df$n, drop = TRUE)
  
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

tail_probability_table <- function(
    obj,
    estimator_pattern = "^SM_m(5|6)_",
    metric = "score_loss",
    thresholds = c(1, 10, 100, 1000, 1e4),
    exclude_normalization_suspect = FALSE
) {
  df <- obj$raw
  df <- df[grepl(estimator_pattern, df$method_label), , drop = FALSE]
  
  if (isTRUE(exclude_normalization_suspect)) {
    df <- df[df$normalization_suspect %in% FALSE, , drop = FALSE]
  }
  
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

# Left right gap
summarise_gap_score_compact <- function(benchmark_obj,
                                        group_vars = c("method_label", "method", "n"),
                                        gap_vars = c("left_gap", "right_gap"),
                                        metric = "score_loss",
                                        eps = 1e-12,
                                        only_success = FALSE,
                                        exclude_normalization_suspect = FALSE) {
  if (is.null(benchmark_obj$raw)) {
    stop("benchmark_obj$raw fehlt.")
  }
  
  df <- benchmark_obj$raw
  
  needed_cols <- unique(c(group_vars, gap_vars, metric, "run_seed"))
  missing_cols <- setdiff(needed_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Fehlende Spalten: ", paste(missing_cols, collapse = ", "))
  }
  
  if (only_success) {
    if (!"success" %in% names(df)) {
      stop("Spalte 'success' fehlt, aber only_success = TRUE wurde gesetzt.")
    }
    df <- df[df$success %in% TRUE, , drop = FALSE]
  }
  
  if (exclude_normalization_suspect) {
    if (!"normalization_suspect" %in% names(df)) {
      stop("Spalte 'normalization_suspect' fehlt, aber exclude_normalization_suspect = TRUE wurde gesetzt.")
    }
    df <- df[!isTRUE(df$normalization_suspect) & !df$normalization_suspect %in% TRUE, , drop = FALSE]
  }
  
  grp <- split(
    df,
    interaction(df[, group_vars, drop = FALSE], drop = TRUE, lex.order = TRUE)
  )
  
  out <- lapply(grp, function(d) {
    res <- d[1, group_vars, drop = FALSE]
    
    gap_info <- vector("list", length(gap_vars))
    names(gap_info) <- gap_vars
    
    for (gap_var in gap_vars) {
      z_raw <- d[[gap_var]]
      valid_idx <- which(is.finite(z_raw))
      
      if (length(valid_idx) == 0L) {
        gap_info[[gap_var]] <- list(
          gap_type = gap_var,
          max_gap = NA_real_,
          max_minus_median = NA_real_,
          run_seed = NA
        )
      } else {
        z <- pmax(z_raw[valid_idx], 0)
        z_median <- stats::median(z)
        
        local_max_idx <- which.max(z)
        row_idx <- valid_idx[local_max_idx]
        z_max <- z[local_max_idx]
        
        gap_info[[gap_var]] <- list(
          gap_type = gap_var,
          max_gap = z_max,
          max_minus_median = z_max - z_median,
          run_seed = d$run_seed[row_idx]
        )
      }
    }
    
    max_gaps <- vapply(gap_info, function(x) x$max_gap, numeric(1))
    
    if (all(is.na(max_gaps))) {
      best_gap <- list(
        gap_type = NA_character_,
        max_gap = NA_real_,
        max_minus_median = NA_real_,
        run_seed = NA
      )
    } else {
      best_name <- names(which.max(replace(max_gaps, is.na(max_gaps), -Inf)))[1]
      best_gap <- gap_info[[best_name]]
    }
    
    res$gap_type <- best_gap$gap_type
    res$max_gap <- best_gap$max_gap
    res$max_minus_median_overall <- best_gap$max_minus_median
    res$run_seed <- best_gap$run_seed
    
    is_zero_left  <- is.finite(d$left_gap)  & abs(d$left_gap)  <= eps
    is_zero_right <- is.finite(d$right_gap) & abs(d$right_gap) <= eps
    
    both_zero <- is_zero_left & is_zero_right
    at_least_one_not_zero <- !both_zero
    
    x_both <- as.numeric(d[[metric]][both_zero])
    x_both <- x_both[is.finite(x_both)]
    
    x_other <- as.numeric(d[[metric]][at_least_one_not_zero])
    x_other <- x_other[is.finite(x_other)]
    
    res$n_both_zero <- sum(both_zero, na.rm = TRUE)
    res$mean_both_zero <- if (length(x_both) == 0L) NA_real_ else mean(x_both)
    res$median_both_zero <- if (length(x_both) == 0L) NA_real_ else stats::median(x_both)
    
    res$mean_at_least_one_not_zero <- if (length(x_other) == 0L) NA_real_ else mean(x_other)
    res$median_at_least_one_not_zero <- if (length(x_other) == 0L) NA_real_ else stats::median(x_other)
    
    res
  })
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

# Calculate Factor od Decrease / Increase of ridge
compare_factor <- function(obj_left,
                                  obj_right,
                                  metric_left = "score_loss",
                                  metric_right = metric_left,
                                  exclude_normalization_suspect = FALSE) {
  
  extract_m <- function(x) {
    suppressWarnings(as.integer(sub(".*_m([0-9]+).*", "\\1", x)))
  }
  
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
  
  agg_left$m  <- extract_m(agg_left$method_label)
  agg_right$m <- extract_m(agg_right$method_label)
  
  left  <- agg_left[,  c("m", "n", "mean")]
  right <- agg_right[, c("m", "n", "method_label", "mean")]
  
  names(left)  <- c("m", "n", "mean_left")
  names(right) <- c("m", "n", "method_label_right", "mean_right")
  
  out <- merge(right, left, by = c("m", "n"))
  
  factor_numeric <- out$mean_left / out$mean_right
  out$Enhance_Performance <- factor_numeric > 1
  
  out$mean_right <- formatC(out$mean_right, format = "e", digits = 3)
  out$mean_left  <- formatC(out$mean_left,  format = "e", digits = 3)
  out$factor     <- formatC(factor_numeric, format = "e", digits = 3)
  
  out <- out[, c(
    "m", "n", "method_label_right",
    "mean_right", "mean_left",
    "factor", "Enhance_Performance"
  )]
  
  out <- out[order(out$m, out$n), ]
  rownames(out) <- NULL
  out
}
