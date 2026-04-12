# ============================================================
# compare_draft_vs_patched.R
# Vergleicht Draft- und patched-Resultate
# ============================================================

compare_one <- function(draft_obj, patched_obj, label) {
  # final_benchmark -> raw Data Frame ziehen
  draft_df <- if (is.data.frame(draft_obj)) draft_obj else draft_obj$raw
  patched_df <- if (is.data.frame(patched_obj)) patched_obj else patched_obj$raw
  
  if (!is.data.frame(draft_df)) {
    stop(sprintf("Draft object for '%s' does not contain a data frame in $raw.", label))
  }
  if (!is.data.frame(patched_df)) {
    stop(sprintf("Patched object for '%s' does not contain a data frame in $raw.", label))
  }
  
  join_keys <- intersect(
    c("n", "repetition", "method_label", "method"),
    intersect(names(draft_df), names(patched_df))
  )
  
  metric_keys <- intersect(
    c("kl", "score_loss"),
    intersect(names(draft_df), names(patched_df))
  )
  
  if (length(join_keys) == 0L) {
    stop(sprintf("No common join keys found for '%s'.", label))
  }
  if (length(metric_keys) == 0L) {
    stop(sprintf("No common metric columns found for '%s'.", label))
  }
  
  d <- draft_df[, c(join_keys, metric_keys), drop = FALSE]
  p <- patched_df[, c(join_keys, metric_keys), drop = FALSE]
  
  cmp <- merge(
    d, p,
    by = join_keys,
    suffixes = c("_draft", "_patched"),
    all = TRUE,
    sort = TRUE
  )
  
  if ("kl_draft" %in% names(cmp) && "kl_patched" %in% names(cmp)) {
    cmp$kl_absdiff <- abs(cmp$kl_draft - cmp$kl_patched)
    identical_kl <- isTRUE(all.equal(cmp$kl_draft, cmp$kl_patched, tolerance = 0))
    max_absdiff_kl <- if (all(is.na(cmp$kl_absdiff))) NA_real_ else max(cmp$kl_absdiff, na.rm = TRUE)
  } else {
    cmp$kl_absdiff <- NA_real_
    identical_kl <- NA
    max_absdiff_kl <- NA_real_
  }
  
  if ("score_loss_draft" %in% names(cmp) && "score_loss_patched" %in% names(cmp)) {
    cmp$score_loss_absdiff <- abs(cmp$score_loss_draft - cmp$score_loss_patched)
    identical_score_loss <- isTRUE(all.equal(cmp$score_loss_draft, cmp$score_loss_patched, tolerance = 0))
    max_absdiff_score_loss <- if (all(is.na(cmp$score_loss_absdiff))) NA_real_ else max(cmp$score_loss_absdiff, na.rm = TRUE)
  } else {
    cmp$score_loss_absdiff <- NA_real_
    identical_score_loss <- NA
    max_absdiff_score_loss <- NA_real_
  }
  
  summary <- data.frame(
    case = label,
    identical_kl = identical_kl,
    identical_score_loss = identical_score_loss,
    max_absdiff_kl = max_absdiff_kl,
    max_absdiff_score_loss = max_absdiff_score_loss
  )
  
  list(summary = summary, comparison = cmp)
}

draft_results <- readRDS("draft_compare_results_seed123.rds")
patched_results <- readRDS("patched_compare_results_seed123.rds")

gaussian_cmp <- compare_one(draft_results$gaussian, patched_results$gaussian, "gaussian")
logistic_cmp <- compare_one(draft_results$logistic, patched_results$logistic, "logistic")

summary_out <- rbind(gaussian_cmp$summary, logistic_cmp$summary)
print(summary_out)

cat("\n--- Gaussian differences ---\n")
print(subset(gaussian_cmp$comparison, kl_absdiff > 0 | score_loss_absdiff > 0))

cat("\n--- Logistic differences ---\n")
print(subset(logistic_cmp$comparison, kl_absdiff > 0 | score_loss_absdiff > 0))


if (!all(summary_out$identical_kl) || !all(summary_out$identical_score_loss)) {
  stop("Draft and patched results are NOT identical.")
} else {
  cat("\nAll compared KL and score_loss values are identical.\n")
}