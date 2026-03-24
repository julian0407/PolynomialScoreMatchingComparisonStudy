# ============================================================
# Multivariate_Polynomial_Score_Matching_Basic.R
# Minimal multivariate polynomial score matching
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
# (4) Fit
# ------------------------------------------------------------

fit_score_matching_mv_basic <- function(x,
                                        m = 2,
                                        include_interactions = TRUE,
                                        standardize = TRUE,
                                        ridge = 1e-8,
                                        solver = c("solve", "qr")) {
  solver <- match.arg(solver)
  
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
  
  basis_obj <- build_poly_basis_mv_basic(
    d = ncol(z),
    m = m,
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
  
  theta <- switch(
    solver,
    solve = tryCatch(
      solve(K_reg, ell),
      error = function(e) qr.solve(K_reg, ell)
    ),
    qr = qr.solve(K_reg, ell)
  )
  
  structure(
    list(
      theta = as.numeric(theta),
      scaling = scaling,
      basis = basis_obj$basis,
      basis_names = basis_obj$basis_names,
      m = m,
      d = basis_obj$d,
      include_interactions = include_interactions,
      ridge = ridge
    ),
    class = "score_matching_mv_basic_fit"
  )
}

# ------------------------------------------------------------
# (5) Predict score
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