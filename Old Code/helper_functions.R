# ============================================================
# helper_functions.R
# Commonly used functions for univariate and multivariate
# score matching and benchmarking
# ============================================================

# ------------------------------------------------------------
# (1) Scaling helpers
# ------------------------------------------------------------

# Helper functions for scaling if data is constant
# Very unlikely edge-case: In general this behavior will not except
# sample size equals one
safe_sd <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) s <- 1
  s
}

# ------------------------------------------------------------
# (1.1) Shift one dimensional data to zero and scale by standard deviation
# ------------------------------------------------------------

# Initilaize and apply scaling on some data
scale_vector_1d <- function(x) {
  # Get numeric vector from input x
  # TODO: Check why necessary 
  x <- as.numeric(x)
  # Count number and proportion of infinite entries
  c_inf <- sum(is.infinite(x))
  prop_inf <- mean(is.infinite(x))
  # Remove infinte entries
  x <- x[is.finite(x)]
  
  # Get mean and standard deviation
  mu <- mean(x)
  s  <- safe_sd(x)
  # Apply shifting and scaling to x and
  # return transformed data and parameters + info about infinite entries
  list(
    x_scaled = (x - mu) / s,
    center = mu,
    scale = s,
    c_inf,
    prop_inf
  )
}

# Apply scaling to other data
apply_scaling_1d <- function(x, scaling) {
  x <- as.numeric(x)
  # Count number and proportion of infinite entries
  c_inf <- sum(is.infinite(x))
  prop_inf <- mean(is.infinite(x))
  # Remove infinte entries
  x <- x[is.finite(x)]
  # Apply shifting and scaling to x and
  # return transformed data + info about infinite entries
  z <- (x - scaling$center) / scaling$scale
  list(
    z = z,
    c_inf = c_inf,
    prop_inf = prop_inf
  )
}

# Back transform data to initial representation
reverse_scaling_1d <- function(z, scaling) {
  z <- as.numeric(z)
  z <- z[is.finite(z)]
  # Apply re-shifting and re-scaling to z and return back-transformed data
  x <- z * scaling$scale + scaling$center
}

# ------------------------------------------------------------
# (1.2) Shift multi-dimensional data to zero and scale by standard deviation
# ------------------------------------------------------------

# Initilaize and apply scaling on some data
scale_matrix_cols <- function(x) {
  x <- as.matrix(x)
  # Identify rows that contains no infinite entry
  keep <- apply(x, 1, function(row) all(is.finite(row)))
  # Count number and proportion of infinite entries
  c_inf <- sum(!keep)
  prop_inf <- mean(!keep)
  # Remove rows with at least one infinite entry and
  # remain matrix property (drop = false) in case only one row remains
  x <- x[keep, , drop = FALSE]

  # Get mean and scale for every dimension
  centers <- colMeans(x)
  scales  <- apply(x, 2, safe_sd)
  
  # Apply shifting and scaling on each data column
  x_scaled <- sweep(x, 2, centers, FUN = "-")
  x_scaled <- sweep(x_scaled, 2, scales, FUN = "/")
  
  # Return scaled data matrix, transforming parameters and infinite
  # data count + proportion
  list(
    x_scaled = x_scaled,
    center = centers,
    scale = scales,
    c_inf,
    prop_inf
  )
}

# Apply scaling to other data
apply_scaling_matrix <- function(x, scaling) {
  x <- as.matrix(x)
  # Identify rows that contains no infinite entry
  keep <- apply(x, 1, function(row) all(is.finite(row)))
  # Count number and proportion of infinite entries
  c_inf <- sum(!keep)
  prop_inf <- mean(!keep)
  # Remove rows with at least one infinite entry and
  # remain matrix property (drop = false) in case only one row remains
  x <- x[keep, , drop = FALSE]
  
  # Apply shifting and scaling to x and
  # return transformed data + info about infinite entries
  z <- sweep(x, 2, scaling$center, FUN = "-")
  z <- sweep(z, 2, scaling$scale, FUN = "/")
  list(
    z = z,
    c_inf = c_inf,
    prop_inf = prop_inf
  )
}

# Back transform data to initial representation
reverse_scaling_matrix <- function(z, scaling) {
  z <- as.matrix(z)
  # Identify rows that contains no infinite entry
  keep <- apply(z, 1, function(row) all(is.finite(row)))
  # remain matrix property (drop = false) in case only one row remains
  z <- z[keep, , drop = FALSE]
  
  # # Apply re-shifting and re-scaling to z and return back-transformed data
  x <- sweep(z, 2, scaling$scale, FUN = "*")
  x <- sweep(x, 2, scaling$center, FUN = "+")
  x
}


# ------------------------------------------------------------
# (2) Weighting function h helpers
# ------------------------------------------------------------
h_tanh_sq <- function(z, tau = 1) {
  u <- abs(z) / tau
  out <- tanh(u)^2
  out
}

h_tanh_sq_prime <- function(z, tau = 1) {
  u <- abs(z) / tau
  sech2 <- 1 / cosh(u)^2
  out <- 2 * tanh(u) * sech2 * sign(z) / tau
  # TODO: Ceck this and motivate this weighting function
  # at z=0 derivative is not smooth because of |z| ->
  # choose zero in this case
  out[z == 0] <- 0
  
  out
}

# ------------------------------------------------------------
# (3) Approximate derivative/gradient to estimate score in
#     log-concave MLE
# ------------------------------------------------------------
# central finite-difference derivative of scalar function
num_derivative_1d <- function(fun, x, h = 1e-4) {
  (fun(x + h) - fun(x - h)) / (2 * h)
}

# central finite-difference gradient of scalar-valued multivariate function
num_gradient_mv <- function(fun, x, h = 1e-4) {
  x <- as.numeric(x)
  d <- length(x)
  out <- numeric(d)
  
  for (j in seq_len(d)) {
    e <- rep(0, d)
    e[j] <- h
    out[j] <- (fun(x + e) - fun(x - e)) / (2 * h)
  }
  out
}

# ------------------------------------------------------------
# (4) Helpers to ensure correct format for input matrices
# ------------------------------------------------------------

# Helper function for the unlikely case that a single point
# estimate is desired
as_obs_matrix <- function(x) {
  if (is.null(dim(x))) {
    matrix(x, nrow = 1)
  } else {
    as.matrix(x)
  }
}




# ------------------------------------------------------------
# (5) To be checked !!!!!
# ------------------------------------------------------------
sanitize_filename_component <- function(x) {
  x <- tolower(as.character(x %||% "unknown"))
  x <- gsub("[^a-z0-9]+", "-", x)
  x <- gsub("(^-+|-+$)", "", x)
  if (!nzchar(x)) x <- "unknown"
  x
}

make_method_tag <- function(estimator_specs = NULL) {
  if (is.null(estimator_specs) || length(estimator_specs) == 0L) return("unknown")
  methods <- unique(vapply(estimator_specs, function(s) s$method %||% "unknown", character(1)))
  if (length(methods) == 1L) return(sanitize_filename_component(methods))
  "mixed"
}

make_benchmark_filename <- function(kind = c("final", "biasvar"),
                                    truth_name = NULL,
                                    estimator_specs = NULL,
                                    family = NULL,
                                    ext = "rds") {
  kind <- match.arg(kind)
  parts <- c(
    kind,
    sanitize_filename_component(truth_name %||% family %||% "benchmark"),
    make_method_tag(estimator_specs),
    format(Sys.time(), "%Y%m%d-%H%M%S")
  )
  paste0(paste(parts, collapse = "_"), ".", ext)
}

apply_optional_log_scale <- function(p, values, axis = c("x", "y"), requested = FALSE) {
  axis <- match.arg(axis)
  if (!isTRUE(requested)) return(p)
  vals <- as.numeric(values)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L || any(vals <= 0)) {
    warning(sprintf("Log scale on axis '%s' skipped because non-positive values are present.", axis), call. = FALSE)
    return(p)
  }
  if (axis == "x") p + ggplot2::scale_x_log10() else p + ggplot2::scale_y_log10()
}
