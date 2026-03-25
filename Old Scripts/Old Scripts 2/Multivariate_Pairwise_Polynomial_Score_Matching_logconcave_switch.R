# ============================================================
# Multivariate_Pairwise_Polynomial_Score_Matching.R
# Minimal multivariate polynomial score matching
# with optional log-concavity switch:
#   - log_concave = "off"  : original model
#   - log_concave = "grid" : original basis + grid penalty on Hessian
#   - log_concave = "m2"   : fit quadratic model (m = 2) and project
#                             Hessian to be negative semidefinite
#
# API:
#   - fit_score_matching_mv_basic()
#   - predict_score_mv_basic()
# ============================================================

# ------------------------------------------------------------
# (1) Basis construction
# ------------------------------------------------------------

build_poly_basis_mv_basic <- function(d, m, include_interactions = TRUE) {
  if (!is.numeric(d) || length(d) != 1L || d < 1L) stop("d must be >= 1.")
  if (!is.numeric(m) || length(m) != 1L || m < 1L) stop("m must be >= 1.")
  
  basis <- vector("list", 0L)
  
  # univariate monomials z_j^r, r = 1,...,m
  for (j in seq_len(d)) {
    for (r in seq_len(m)) {
      basis[[length(basis) + 1L]] <- list(
        type = "uni",
        j = j,
        r = r
      )
    }
  }
  
  # pairwise monomials z_j^r z_k^s, r,s = 1,...,m-1
  if (isTRUE(include_interactions) && d >= 2L && m >= 2L) {
    for (j in seq_len(d - 1L)) {
      for (k in (j + 1L):d) {
        for (r in seq_len(m - 1L)) {
          for (s in seq_len(m - 1L)) {
            basis[[length(basis) + 1L]] <- list(
              type = "pair",
              j = j,
              k = k,
              r = r,
              s = s
            )
          }
        }
      }
    }
  }
  
  basis_names <- vapply(basis, function(b) {
    if (b$type == "uni") {
      sprintf("z%d^%d", b$j, b$r)
    } else {
      sprintf("z%d^%d*z%d^%d", b$j, b$r, b$k, b$s)
    }
  }, character(1))
  
  list(
    basis = basis,
    basis_names = basis_names,
    p = length(basis),
    d = d,
    m = m,
    include_interactions = include_interactions
  )
}

# ------------------------------------------------------------
# (2) Power cache
# ------------------------------------------------------------

build_power_cache_mv_basic <- function(z, max_deg) {
  z <- as_obs_matrix(z)
  n <- nrow(z)
  d <- ncol(z)
  
  out <- vector("list", d)
  for (j in seq_len(d)) {
    P <- matrix(1, nrow = n, ncol = max_deg + 1L)
    if (max_deg >= 1L) {
      P[, 2L] <- z[, j]
      if (max_deg >= 2L) {
        for (r in 2:max_deg) {
          P[, r + 1L] <- P[, r] * z[, j]
        }
      }
    }
    out[[j]] <- P
  }
  out
}

get_power_mv_basic <- function(P, j, deg) {
  P[[j]][, deg + 1L]
}

# ------------------------------------------------------------
# (3) Build derivative design + Laplacian
# ------------------------------------------------------------

build_sm_design_mv_basic <- function(z, basis_obj) {
  z <- as_obs_matrix(z)
  if (ncol(z) != basis_obj$d) stop("Dimension mismatch.")
  
  n <- nrow(z)
  d <- basis_obj$d
  p <- basis_obj$p
  
  max_deg <- basis_obj$m
  P <- build_power_cache_mv_basic(z, max_deg = max_deg)
  
  D <- lapply(seq_len(d), function(k) matrix(0, nrow = n, ncol = p))
  Lap <- matrix(0, nrow = n, ncol = p)
  
  for (col in seq_len(p)) {
    b <- basis_obj$basis[[col]]
    
    if (b$type == "uni") {
      if (b$r >= 1L) {
        D[[b$j]][, col] <- b$r * get_power_mv_basic(P, b$j, b$r - 1L)
      }
      if (b$r >= 2L) {
        Lap[, col] <- b$r * (b$r - 1L) * get_power_mv_basic(P, b$j, b$r - 2L)
      }
      
    } else {
      zj_r <- get_power_mv_basic(P, b$j, b$r)
      zk_s <- get_power_mv_basic(P, b$k, b$s)
      
      D[[b$j]][, col] <- b$r * get_power_mv_basic(P, b$j, b$r - 1L) * zk_s
      D[[b$k]][, col] <- b$s * zj_r * get_power_mv_basic(P, b$k, b$s - 1L)
      
      term_j <- if (b$r >= 2L) {
        b$r * (b$r - 1L) * get_power_mv_basic(P, b$j, b$r - 2L) * zk_s
      } else {
        0
      }
      
      term_k <- if (b$s >= 2L) {
        b$s * (b$s - 1L) * zj_r * get_power_mv_basic(P, b$k, b$s - 2L)
      } else {
        0
      }
      
      Lap[, col] <- term_j + term_k
    }
  }
  
  for (k in seq_len(d)) colnames(D[[k]]) <- basis_obj$basis_names
  colnames(Lap) <- basis_obj$basis_names
  
  list(D = D, Lap = Lap)
}

# ------------------------------------------------------------
# (4) Hessian helpers for log-concavity
# ------------------------------------------------------------

build_hessian_design_mv_basic <- function(z, basis_obj) {
  z <- as_obs_matrix(z)
  if (ncol(z) != basis_obj$d) stop("Dimension mismatch.")
  
  n <- nrow(z)
  d <- basis_obj$d
  p <- basis_obj$p
  P <- build_power_cache_mv_basic(z, max_deg = basis_obj$m)
  
  H <- vector("list", d)
  for (a in seq_len(d)) {
    H[[a]] <- vector("list", d)
    for (b in seq_len(d)) {
      H[[a]][[b]] <- matrix(0, nrow = n, ncol = p)
    }
  }
  
  for (col in seq_len(p)) {
    b <- basis_obj$basis[[col]]
    
    if (b$type == "uni") {
      if (b$r >= 2L) {
        H[[b$j]][[b$j]][, col] <- b$r * (b$r - 1L) * get_power_mv_basic(P, b$j, b$r - 2L)
      }
      
    } else {
      zj_r <- get_power_mv_basic(P, b$j, b$r)
      zk_s <- get_power_mv_basic(P, b$k, b$s)
      
      if (b$r >= 2L) {
        H[[b$j]][[b$j]][, col] <- b$r * (b$r - 1L) * get_power_mv_basic(P, b$j, b$r - 2L) * zk_s
      }
      if (b$s >= 2L) {
        H[[b$k]][[b$k]][, col] <- b$s * (b$s - 1L) * zj_r * get_power_mv_basic(P, b$k, b$s - 2L)
      }
      if (b$r >= 1L && b$s >= 1L) {
        cross_term <- b$r * b$s * get_power_mv_basic(P, b$j, b$r - 1L) * get_power_mv_basic(P, b$k, b$s - 1L)
        H[[b$j]][[b$k]][, col] <- cross_term
        H[[b$k]][[b$j]][, col] <- cross_term
      }
    }
  }
  
  H
}

make_logconcavity_grid_mv_basic <- function(z,
                                            grid_size = 5L,
                                            box_expand = 0,
                                            max_points = 500L,
                                            seed = NULL) {
  z <- as_obs_matrix(z)
  d <- ncol(z)
  
  if (!is.numeric(grid_size) || length(grid_size) != 1L || grid_size < 2L) {
    stop("lc_grid_size must be an integer >= 2.")
  }
  if (!is.numeric(box_expand) || length(box_expand) != 1L || box_expand < 0) {
    stop("lc_box_expand must be a nonnegative scalar.")
  }
  if (!is.numeric(max_points) || length(max_points) != 1L || max_points < 1L) {
    stop("lc_max_points must be a positive integer.")
  }
  
  mins <- apply(z, 2, min)
  maxs <- apply(z, 2, max)
  spans <- pmax(maxs - mins, 1e-8)
  mins <- mins - box_expand * spans
  maxs <- maxs + box_expand * spans
  
  grid_list <- lapply(seq_len(d), function(j) seq(mins[j], maxs[j], length.out = grid_size))
  grid_df <- do.call(expand.grid, grid_list)
  grid <- as.matrix(grid_df)
  colnames(grid) <- paste0("z", seq_len(d))
  
  if (nrow(grid) > max_points) {
    if (!is.null(seed)) set.seed(seed)
    idx <- sample.int(nrow(grid), size = max_points)
    grid <- grid[idx, , drop = FALSE]
  }
  
  grid
}

hessian_array_from_theta_mv_basic <- function(theta, hess_design, d) {
  n <- nrow(hess_design[[1L]][[1L]])
  out <- array(0, dim = c(n, d, d))
  for (a in seq_len(d)) {
    for (b in seq_len(d)) {
      out[, a, b] <- as.numeric(hess_design[[a]][[b]] %*% theta)
    }
  }
  out
}

logconcavity_violation_mv_basic <- function(theta, hess_design, d, tol = 1e-8) {
  H_arr <- hessian_array_from_theta_mv_basic(theta, hess_design, d = d)
  n <- dim(H_arr)[1L]
  max_eigs <- numeric(n)
  
  for (i in seq_len(n)) {
    Hi <- H_arr[i, , , drop = TRUE]
    Hi <- 0.5 * (Hi + t(Hi))
    max_eigs[i] <- max(eigen(Hi, symmetric = TRUE, only.values = TRUE)$values)
  }
  
  violations <- pmax(max_eigs - tol, 0)
  list(
    max_eigenvalues = max_eigs,
    violations = violations,
    max_violation = max(violations),
    mean_violation = mean(violations),
    n_violated = sum(violations > 0)
  )
}

penalized_sm_objective_mv_basic <- function(theta,
                                            K_reg,
                                            ell,
                                            hess_design = NULL,
                                            d = NULL,
                                            penalty = 0,
                                            tol = 1e-8) {
  base_obj <- 0.5 * drop(crossprod(theta, K_reg %*% theta)) - sum(theta * ell)
  
  if (is.null(hess_design) || penalty <= 0) {
    return(base_obj)
  }
  
  lc_diag <- logconcavity_violation_mv_basic(theta, hess_design, d = d, tol = tol)
  base_obj + penalty * mean(lc_diag$violations^2)
}

project_theta_logconcave_m2_mv_basic <- function(theta, basis, d, eps = 1e-8) {
  beta <- numeric(d)
  H <- matrix(0, nrow = d, ncol = d)
  
  for (col in seq_along(basis)) {
    b <- basis[[col]]
    val <- theta[col]
    
    if (b$type == "uni") {
      if (b$r == 1L) {
        beta[b$j] <- val
      } else if (b$r == 2L) {
        H[b$j, b$j] <- H[b$j, b$j] + 2 * val
      }
    } else if (b$r == 1L && b$s == 1L) {
      H[b$j, b$k] <- H[b$j, b$k] + val
      H[b$k, b$j] <- H[b$k, b$j] + val
    }
  }
  
  H <- 0.5 * (H + t(H))
  ee <- eigen(H, symmetric = TRUE)
  ee$values <- pmin(ee$values, -eps)
  H_proj <- ee$vectors %*% diag(ee$values, d, d) %*% t(ee$vectors)
  
  theta_proj <- numeric(length(theta))
  for (col in seq_along(basis)) {
    b <- basis[[col]]
    if (b$type == "uni") {
      if (b$r == 1L) theta_proj[col] <- beta[b$j]
      if (b$r == 2L) theta_proj[col] <- 0.5 * H_proj[b$j, b$j]
    } else if (b$r == 1L && b$s == 1L) {
      theta_proj[col] <- H_proj[b$j, b$k]
    }
  }
  
  list(
    theta = theta_proj,
    beta = beta,
    H = H_proj,
    eigenvalues = ee$values
  )
}

# ------------------------------------------------------------
# (5) Fit
# ------------------------------------------------------------

fit_score_matching_mv_basic <- function(x,
                                        m = 2,
                                        include_interactions = TRUE,
                                        standardize = TRUE,
                                        ridge = 1e-8,
                                        solver = c("solve", "qr"),
                                        log_concave = c("off", "grid", "m2"),
                                        lc_grid_size = 5L,
                                        lc_box_expand = 0,
                                        lc_max_points = 500L,
                                        lc_penalty = 1e4,
                                        lc_tol = 1e-8,
                                        lc_optim_method = c("BFGS", "Nelder-Mead"),
                                        lc_control = list(maxit = 500),
                                        lc_seed = NULL,
                                        lc_m2_eps = 1e-8) {
  solver <- match.arg(solver)
  log_concave <- match.arg(log_concave)
  lc_optim_method <- match.arg(lc_optim_method)
  
  x <- as_obs_matrix(x)
  keep <- apply(x, 1, function(row) all(is.finite(row)))
  x <- x[keep, , drop = FALSE]
  if (nrow(x) < 2L) stop("x must have at least two finite rows.")
  
  if (isTRUE(standardize)) {
    sc <- scale_matrix_cols(x)
    z <- sc$x_scaled
    scaling <- list(
      center = sc$center,
      scale = sc$scale,
      standardize = TRUE
    )
  } else {
    z <- x
    scaling <- list(
      center = rep(0, ncol(x)),
      scale = rep(1, ncol(x)),
      standardize = FALSE
    )
  }
  
  m_fit <- if (identical(log_concave, "m2")) 2L else m
  
  basis_obj <- build_poly_basis_mv_basic(
    d = ncol(z),
    m = m_fit,
    include_interactions = include_interactions
  )
  
  design <- build_sm_design_mv_basic(z, basis_obj)
  
  K <- matrix(0, nrow = basis_obj$p, ncol = basis_obj$p)
  for (k in seq_len(basis_obj$d)) {
    K <- K + crossprod(design$D[[k]])
  }
  K <- K / nrow(z)
  
  ell <- colMeans(design$Lap)
  
  if (!is.numeric(ridge) || length(ridge) != 1L || ridge < 0) {
    stop("ridge must be a nonnegative scalar.")
  }
  
  K_reg <- K + ridge * diag(basis_obj$p)
  
  theta_init <- switch(
    solver,
    solve = tryCatch(
      solve(K_reg, ell),
      error = function(e) qr.solve(K_reg, ell)
    ),
    qr = qr.solve(K_reg, ell)
  )
  theta_init <- as.numeric(theta_init)
  
  theta <- theta_init
  lc_points <- NULL
  lc_diag <- NULL
  opt_details <- NULL
  m2_details <- NULL
  
  if (identical(log_concave, "grid")) {
    lc_points <- make_logconcavity_grid_mv_basic(
      z = z,
      grid_size = lc_grid_size,
      box_expand = lc_box_expand,
      max_points = lc_max_points,
      seed = lc_seed
    )
    lc_hess_design <- build_hessian_design_mv_basic(lc_points, basis_obj)
    
    opt <- optim(
      par = theta_init,
      fn = penalized_sm_objective_mv_basic,
      K_reg = K_reg,
      ell = ell,
      hess_design = lc_hess_design,
      d = basis_obj$d,
      penalty = lc_penalty,
      tol = lc_tol,
      method = lc_optim_method,
      control = lc_control
    )
    
    theta <- as.numeric(opt$par)
    lc_diag <- logconcavity_violation_mv_basic(theta, lc_hess_design, d = basis_obj$d, tol = lc_tol)
    opt_details <- list(
      convergence = opt$convergence,
      message = opt$message,
      value = opt$value,
      counts = opt$counts,
      method = lc_optim_method
    )
  }
  
  if (identical(log_concave, "m2")) {
    m2_proj <- project_theta_logconcave_m2_mv_basic(
      theta = theta_init,
      basis = basis_obj$basis,
      d = basis_obj$d,
      eps = lc_m2_eps
    )
    theta <- as.numeric(m2_proj$theta)
    m2_details <- list(
      H = m2_proj$H,
      beta = m2_proj$beta,
      max_eigenvalue = max(eigen(m2_proj$H, symmetric = TRUE, only.values = TRUE)$values),
      eps = lc_m2_eps
    )
  }
  
  structure(
    list(
      theta = theta,
      theta_unconstrained = theta_init,
      scaling = scaling,
      basis = basis_obj$basis,
      basis_names = basis_obj$basis_names,
      m = basis_obj$m,
      m_requested = m,
      d = basis_obj$d,
      include_interactions = include_interactions,
      ridge = ridge,
      log_concave = log_concave,
      lc_points = lc_points,
      lc_diagnostics = lc_diag,
      optimizer = opt_details,
      m2_details = m2_details
    ),
    class = "score_matching_mv_basic_fit"
  )
}

# ------------------------------------------------------------
# (6) Predict score
# ------------------------------------------------------------

predict_score_mv_basic <- function(newx, fit) {
  if (!inherits(fit, "score_matching_mv_basic_fit")) {
    stop("fit must be of class 'score_matching_mv_basic_fit'.")
  }
  
  newx <- as_obs_matrix(newx)
  keep <- apply(newx, 1, function(row) all(is.finite(row)))
  newx <- newx[keep, , drop = FALSE]
  
  if (ncol(newx) != fit$d) stop("Dimension mismatch.")
  
  if (isTRUE(fit$scaling$standardize)) {
    z <- apply_scaling_matrix(newx, fit$scaling)$z
  } else {
    z <- newx
  }
  
  basis_obj <- list(
    basis = fit$basis,
    basis_names = fit$basis_names,
    p = length(fit$basis),
    d = fit$d,
    m = fit$m
  )
  
  design <- build_sm_design_mv_basic(z, basis_obj)
  
  score_z <- do.call(
    cbind,
    lapply(seq_len(fit$d), function(k) design$D[[k]] %*% fit$theta)
  )
  
  # r_z(z) = diag(scale) r_x(x)  =>  r_x(x) = r_z(z) / scale
  score_x <- sweep(score_z, 2, fit$scaling$scale, FUN = "/")
  colnames(score_x) <- paste0("dim", seq_len(fit$d))
  score_x
}
