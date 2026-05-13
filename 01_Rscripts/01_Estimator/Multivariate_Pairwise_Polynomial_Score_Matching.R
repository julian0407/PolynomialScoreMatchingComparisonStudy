# ============================================================
# Multivariate_Pairwise_Polynomial_Score_Matching.R
# with optional log-concavity switch:
#   - log_concave = "off"  : Just pairiwse polynomial basis
#   - log_concave = "grid" : original basis + grid penalty on Hessian
#
# API:
#   - fit_score_matching_mv()
#   - predict_score_mv()
#
# Convention:
#   score r(x) = - grad log f(x)
# ============================================================

# ------------------------------------------------------------
# (1) Basis construction
# ------------------------------------------------------------

build_poly_basis_mv <- function(d, m, include_interactions = TRUE) {
  # Check if d and m are positive
  if (!is.numeric(d) || length(d) != 1L || d < 1L) stop("d must be >= 1.")
  if (!is.numeric(m) || length(m) != 1L || m < 1L) stop("m must be >= 1.")
  
  # initialize empty list that will save meta data for each basis term
  basis <- vector("list", 0L)
  
  # for each dimension j and each degree r save meta data for 
  # univariate bases in this list as a nested list object
  for (j in seq_len(d)) {
    for (r in seq_len(m)) {
      basis[[length(basis) + 1L]] <- list(
        type = "uni",
        j = j,
        r = r
      )
    }
  }
  
  # Build interaction terms only if parameter inlcude_interaction = True
  # and d and m >=2
  if (isTRUE(include_interactions) && d >= 2L && m >= 2L) {
    # Add to the list also the metadata for pairwise bases in this list
    # such that 1<=j<k<=d for pairwise interaction of dimensions j and k
    for (j in seq_len(d - 1L)) {
      for (k in (j + 1L):d) {
        # Do not contsruct univariate bases again -> stop at m-1
        for (r in seq_len(m - 1L)) {
          # Only add bases with 'combined' degree less or equal than m
          for (s in seq_len(m - r)) {
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
  
  # Create and save names for each basis term in variable basis_names
  basis_names <- vapply(basis, function(b) {
    if (b$type == "uni") {
      sprintf("z%d^%d", b$j, b$r)
    } else {
      sprintf("z%d^%d*z%d^%d", b$j, b$r, b$k, b$s)
    }
  }, character(1))
  
  # Return created basis and basis_names, the length p of the basis,
  # the dimension d and degree m which has been used to create the basis
  # and the boolean if the basis includes interaction terms
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

# Create cache of 'powered' data for efficient calculation of derivatives and Hessian later
build_power_cache_mv <- function(z, max_deg) {
  # Ensure fitting data is in matrix format
  z <- as_obs_matrix(z)
  # Get dimension d and sample size n of data
  n <- nrow(z)
  d <- ncol(z)
  
  # Initialize a output list which will save for each entry a matrix that
  # contains for each data point (rows) the power of this sample (columns)
  # up to max_deg
  out <- vector("list", d)
  for (j in seq_len(d)) {
    # First column equals 1 (power to zero)
    P <- matrix(1, nrow = n, ncol = max_deg + 1L)
    if (max_deg >= 1L) {
      # subsequent columns uses recursive calculation of power to r
      for (r in 1:max_deg) {
        P[, r + 1L] <- P[, r] * z[, j]
      }
    }
    out[[j]] <- P
  }
  out
}

# Helper function to receive powered data in dimension j
# with power degree deg
get_power_mv <- function(P, j, deg) {
  P[[j]][, deg + 1L]
}

# ------------------------------------------------------------
# (3) Build derivative design + Laplacian
# ------------------------------------------------------------

# Building matrices D and Laplacian based on the basis functions in basis_obj
build_sm_design_mv <- function(z, basis_obj) {
  # Ensure fitting data is in matrix format
  z <- as_obs_matrix(z)
  # Check if basis functions were created on the same dimension as data sample z
  if (ncol(z) != basis_obj$d) stop("Dimension mismatch.")
  
  # Get sample size n, dimension d and number of basis functions p
  n <- nrow(z)
  d <- basis_obj$d
  p <- basis_obj$p
  
  # Build cache for powered sample with helper function
  P <- build_power_cache_mv(z, max_deg = basis_obj$m)
  
  # Initialize D as list of empty matrices of site n*p
  D <- lapply(seq_len(d), function(k) matrix(0, nrow = n, ncol = p))
  # Initialize Laplacian term as empty matrix of size n*p
  Lap <- matrix(0, nrow = n, ncol = p)
  
  # Iterate over all basis terms (advantage of this representation
  # -> do not iterate over sample size n)
  for (col in seq_len(p)) {
    # Get basis function at position col
    b <- basis_obj$basis[[col]]
    
    # Check if basis function is univariate or a pairwise polynomial
    if (b$type == "uni") {
      if (b$r >= 1L) {
        # if degree of univariate polynomial is at least one save derivative r*z^{r-1}
        # for dimension j in list D
        D[[b$j]][, col] <- b$r * get_power_mv(P, b$j, b$r - 1L)
      }
      if (b$r >= 2L) {
        # if degree of univariate polynomial is at least two than the second derivative 
        # of this basis term must be saved in the Laplacian
        Lap[, col] <- b$r * (b$r - 1L) * get_power_mv(P, b$j, b$r - 2L)
      }
      
    } else {
      # if the basis function is pairwise polynomial we first derive the corresponding
      # cached data for each dimesnion j,k and degree r,s
      zj_r <- get_power_mv(P, b$j, b$r)
      zk_s <- get_power_mv(P, b$k, b$s)
      
      # Similar to the univariate case we save te first partial derivative into structure D
      # dependent on which dimension j,k is derivated
      D[[b$j]][, col] <- b$r * get_power_mv(P, b$j, b$r - 1L) * zk_s
      D[[b$k]][, col] <- b$s * zj_r * get_power_mv(P, b$k, b$s - 1L)
      
      # For the laplacian we must calculate the second partial derivative in regard
      # to the participating dimensions j and k of the pairwise polynomial basis
      # They only exists of the "degree of dimension j or k" in the basis function
      # is at least two. 
      term_j <- if (b$r >= 2L) {
        b$r * (b$r - 1L) * get_power_mv(P, b$j, b$r - 2L) * zk_s
      } else {
        0
      }
      
      term_k <- if (b$s >= 2L) {
        b$s * (b$s - 1L) * zj_r * get_power_mv(P, b$k, b$s - 2L)
      } else {
        0
      }
      # The sum of the two second partial derivatives is the contribution to the Laplacian
      Lap[, col] <- term_j + term_k
    }
  }
  
  # Set coulmn names in each matrix in D representing contribution of dimension k
  for (k in seq_len(d)) colnames(D[[k]]) <- basis_obj$basis_names
  colnames(Lap) <- basis_obj$basis_names
  
  # Return list structure D and matrix Lap
  list(D = D, Lap = Lap)
}


# ------------------------------------------------------------
# (4) grid penalty: log-concavity helpers
# ------------------------------------------------------------

# Build Hessian matrix for data z and a given basis_obj used later for
# enforce H being strictly positive definite at this points
build_hessian_design_mv <- function(z, basis_obj) {
  # Ensure fitting data is in matrix format
  z <- as_obs_matrix(z)
  # Check if basis functions were created on the same dimension as data sample z
  if (ncol(z) != basis_obj$d) stop("Dimension mismatch.")
  
  # Get sample size n, dimension d and number of basis functions p
  n <- nrow(z)
  d <- basis_obj$d
  p <- basis_obj$p
  
  # Build cache for powered sample with helper function
  P <- build_power_cache_mv(z, max_deg = basis_obj$m)
  
  # Initialize H as list of size d
  H <- vector("list", d)
  # Iterate through this list and initialize an empty list of size d at each entry
  for (a in seq_len(d)) {
    H[[a]] <- vector("list", d)
    for (b in seq_len(d)) {
      # The nested list structure can be interpreted as matrix H and at each entry
      # we initlaize an empty matrix of size n*p
      H[[a]][[b]] <- matrix(0, nrow = n, ncol = p)
    }
  }
  
  # Now iterate through all basis function
  for (col in seq_len(p)) {
    b <- basis_obj$basis[[col]]
    
    # Univariate basis functions only contributes to the diagonal of H
    # Calculate the contribution to the Hessian of this basis to
    # the second partial derivative in dimension j
    if (b$type == "uni") {
      if (b$r >= 2L) {
        H[[b$j]][[b$j]][, col] <- b$r * (b$r - 1L) * get_power_mv(P, b$j, b$r - 2L)
      }
      
    } else {
      # if basis in a pairwise polynomial function we must calculate the
      # contribution of this basis to the Hessian to the diagonal and
      # to the cross derivations not on the diagonal
      # Get first the evaluated basis vectors for the data vector z
      zj_r <- get_power_mv(P, b$j, b$r)
      zk_s <- get_power_mv(P, b$k, b$s)
      
      if (b$r >= 2L) {
        H[[b$j]][[b$j]][, col] <- b$r * (b$r - 1L) * get_power_mv(P, b$j, b$r - 2L) * zk_s
      }
      if (b$s >= 2L) {
        H[[b$k]][[b$k]][, col] <- b$s * (b$s - 1L) * zj_r * get_power_mv(P, b$k, b$s - 2L)
      }
      if (b$r >= 1L && b$s >= 1L) {
        cross_term <- b$r * b$s * get_power_mv(P, b$j, b$r - 1L) * get_power_mv(P, b$k, b$s - 1L)
        H[[b$j]][[b$k]][, col] <- cross_term
        H[[b$k]][[b$j]][, col] <- cross_term
      }
    }
  }
  
  # For every matrix saved in this nested list structure apply the 
  # naming of the basis_obj
  for (a in seq_len(d)) {
    for (b in seq_len(d)) {
      colnames(H[[a]][[b]]) <- basis_obj$basis_names
    }
  }
  H
}

# Generate a grid in matrix format that is used to enforce log-concave assumption at these points
# (Hessian > 0)
make_logconcavity_points_mv <- function(z,
                                              method = c("grid", "data"),
                                              grid_size = 5L,
                                              box_expand = 0) {
  # Check if argument method is valid
  method <- match.arg(method)
  # Ensure z is in matrix format
  z <- as_obs_matrix(z)
  # Get dimension d of data
  d <- ncol(z)
  
  # If method equals data,just retzrn z
  if (method == "data") {
    return(z)
  }
  
  # error handling for not valid input parameters
  if (!is.numeric(grid_size) || length(grid_size) != 1L || grid_size < 2L) {
    stop("grid_size must be an integer >= 2 when lc_method = 'grid'.")
  }
  if (!is.numeric(box_expand) || length(box_expand) != 1L || box_expand < 0) {
    stop("box_expand must be a nonnegative scalar.")
  }

  # Get for each column the max and min
  mins <- apply(z, 2, min)
  maxs <- apply(z, 2, max)
  # Caluclate difference and ensure the difference is greater than zero
  spans <- pmax(maxs - mins, 1e-8)
  # Optionally expand box bei input parameter and span
  mins <- mins - box_expand * spans
  maxs <- maxs + box_expand * spans
  
  # Derive grid in each dimension
  grid_list <- lapply(seq_len(d), function(j) seq(mins[j], maxs[j], length.out = grid_size))
  # Apply cartesian product and save values in df
  grid_df <- do.call(expand.grid, grid_list)
  grid <- as.matrix(grid_df)
  # Name the columns
  colnames(grid) <- paste0("z", seq_len(d))
  
  grid
}

# Calculate actual Hessian based on the model parameter theta
hessian_array_from_theta_mv <- function(theta, hess_design, d) {
  # Get sample size based on number of rows at an arbitrary matrix entry in the
  # Hessian design matrix (Here at (1,1))
  n <- nrow(hess_design[[1L]][[1L]])
  # Initialize an array containing zeros and for each sample a "d*d matrix"
  out <- array(0, dim = c(n, d, d))
  for (a in seq_len(d)) {
    for (b in seq_len(d)) {
      # Calculate the linear combination of the weighted basis function based on theta
      # for each sample
      out[, a, b] <- as.numeric(hess_design[[a]][[b]] %*% theta)
    }
  }
  out
}


logconcavity_violation_mv <- function(theta, hess_design, d, tol = 1e-8) {
  # Derive actual hessian for each sample
  H_arr <- hessian_array_from_theta_mv(theta, hess_design, d = d)
  # Get sample size n
  n <- dim(H_arr)[1L]
  # Initialize vector that saves the minimum eigenvalue of the Hessian of each sample
  min_eigs <- numeric(n)
  
  # Get for each sample the Hessian
  for (i in seq_len(n)) {
    Hi <- H_arr[i, , , drop = TRUE]
    # Ensure Hessian is indeed symmteric
    Hi <- 0.5 * (Hi + t(Hi))
    # Calculate minimal eigenvalue
    min_eigs[i] <- min(eigen(Hi, symmetric = TRUE, only.values = TRUE)$values)
  }
  
  # Negative less than the tolernac contributes as violation
  violations <- pmax(-min_eigs - tol, 0)
  # Return violations and meta data to analyze
  list(
    min_eigenvalues = min_eigs,
    violations = violations,
    max_violation = max(violations),
    mean_violation = mean(violations),
    n_violated = sum(violations > 0)
  )
}

# objective function for log concave penalized multivariate sm with grid
penalized_sm_objective_mv <- function(theta,
                                            K_reg,
                                            ell,
                                            hess_design = NULL,
                                            d = NULL,
                                            penalty = 0,
                                            tol = 1e-8) {
  # basic pairwise polynomial sm objective theta^T*K*theta - theta^T*l
  base_obj <- 0.5 * drop(crossprod(theta, K_reg %*% theta)) - sum(theta * ell)
  
  # if no hessian_design object is provided return standard objective
  if (is.null(hess_design) || penalty <= 0) {
    return(base_obj)
  }
  
  # Get penalize object from based on hess design and theta
  lc_penalize <- logconcavity_violation_mv(theta, hess_design, d = d, tol = tol)
  # Derive penalty and add it to objective function
  pen <- penalty * mean(lc_penalize$violations^2)
  base_obj + pen
}

# ------------------------------------------------------------
# (5) Fit
#   log_concave = FALSE/TRUE
#   lc_method   = "grid" / "data" / "m2"
# ------------------------------------------------------------

fit_score_matching_mv <- function(x,
                                        m = 2,
                                        include_interactions = TRUE,
                                        standardize = TRUE,
                                        ridge = 1e-8,
                                        solver = c("qr"),
                                        log_concave = FALSE,
                                        lc_method = c("grid", "data", "m2"),
                                        lc_grid_size = 5L,
                                        lc_box_expand = 0,
                                        lc_penalty = 1e4,
                                        lc_tol = 1e-8,
                                        lc_optim_method = c("BFGS"),
                                        lc_control = list(maxit = 500),
                                        lc_seed = NULL,
                                        lc_m2_eps = 1e-8) {
  # check if input parameters are valid
  solver <- match.arg(solver)
  lc_method <- match.arg(lc_method)
  lc_optim_method <- match.arg(lc_optim_method)
  
  # Ensure that fitting data is in matrix format and drop infinite data
  x <- as_obs_matrix(x)
  keep <- apply(x, 1, function(row) all(is.finite(row)))
  x <- x[keep, , drop = FALSE]
  # Check if at least two samples still exist
  if (nrow(x) < 2L) stop("x must have at least two finite rows.")
  
  # Standardize data if standardize = True
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
  
  # Ensure that m=2 if lc_method="2m"
  m_fit <- if (isTRUE(log_concave) && identical(lc_method, "m2")) 2L else m
  
  # Build basis
  basis_obj <- build_poly_basis_mv(
    d = ncol(z),
    m = m_fit,
    include_interactions = include_interactions
  )
  
  # Build structure D and laplacian
  design <- build_sm_design_mv(z, basis_obj)
  
  # Initialze matrix K as empty matrix
  K <- matrix(0, nrow = basis_obj$p, ncol = basis_obj$p)
  # Caluclate sum of crossprodukt of each matrix D_k in structure D
  for (k in seq_len(basis_obj$d)) {
    K <- K + crossprod(design$D[[k]])
  }
  K <- K / nrow(z)
  
  # ell equals Laplacian
  ell <- colMeans(design$Lap)
  
  # Check if ridge factor is valid
  if (!is.numeric(ridge) || length(ridge) != 1L || ridge < 0) {
    stop("ridge must be a nonnegative scalar.")
  }
  
  # Add ridge term to K
  K_reg <- K + ridge * diag(basis_obj$p)
  
  # Get initial theta
  # In case standard solve leads to error automatic fallback to qr
  theta_init <- switch(
    solver,
    solve = tryCatch(
      solve(K_reg, ell),
      error = function(e) qr.solve(K_reg, ell)
    ),
    qr = qr.solve(K_reg, ell)
  )
  theta_init <- as.numeric(theta_init)
  
  # Initialize log-concave constraints/details as NULL
  lc_points <- NULL
  lc_hess_design <- NULL
  lc_diag <- NULL
  opt_details <- NULL
  m2_details <- NULL
  theta <- theta_init
  
  # Apply loc-concave methods if selected
  if (isTRUE(log_concave)) {
    # if lc_method = grid or lc_method = data
    # Build grid dependent on lc_method and input parameters
    lc_points <- make_logconcavity_points_mv(
      z = z,
      method = lc_method,
      grid_size = lc_grid_size,
      box_expand = lc_box_expand
    )
    # Build Hessian_design on grid points
    lc_hess_design <- build_hessian_design_mv(lc_points, basis_obj)
    
    # Iteratively minimize score matching loss with additional
    # Hessian-based penalty enforcing log-concavity
    opt <- optim(
      par = theta_init,
      fn = penalized_sm_objective_mv,
      K_reg = K_reg,
      ell = ell,
      hess_design = lc_hess_design,
      d = basis_obj$d,
      penalty = lc_penalty,
      tol = lc_tol,
      method = lc_optim_method,
      control = lc_control
    )
    
    # Update theta
    theta <- as.numeric(opt$par)
    # Compute diagnostics from current solution
    lc_diag <- logconcavity_violation_mv(theta, lc_hess_design, d = basis_obj$d, tol = lc_tol)
    # Save optimization details
    opt_details <- list(
      convergence = opt$convergence,
      message = opt$message,
      value = opt$value,
      counts = opt$counts,
      method = lc_optim_method
    )
  }

  # Return fitting and optimization details
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
      log_concave = isTRUE(log_concave),
      lc_method = if (isTRUE(log_concave)) lc_method else NULL,
      lc_points = lc_points,
      lc_settings = if (isTRUE(log_concave)) {
        list(
          method = lc_method,
          grid_size = lc_grid_size,
          box_expand = lc_box_expand,
          penalty = lc_penalty,
          tol = lc_tol,
          optim_method = lc_optim_method,
          control = lc_control,
          m2_eps = lc_m2_eps
        )
      } else {
        NULL
      },
      lc_diagnostics = lc_diag,
      optimizer = opt_details,
      m2_details = m2_details
    ),
    class = "score_matching_mv_fit"
  )
}

# ------------------------------------------------------------
# (6) Predict score
# ------------------------------------------------------------

predict_score_mv <- function(newx, fit) {
  # Enforce newx to be in matrix format and drop infinite samples
  newx <- as_obs_matrix(newx)
  keep <- apply(newx, 1, function(row) all(is.finite(row)))
  newx <- newx[keep, , drop = FALSE]
  
  # Check if fitting and estimating dimensions match
  if (ncol(newx) != fit$d) stop("Dimension mismatch.")
  
  # Apply scaling if standardize = True
  if (isTRUE(fit$scaling$standardize)) {
    z <- apply_scaling_matrix(newx, fit$scaling)$z
  } else {
    z <- newx
  }
  
  # Get basis_obj from fit
  basis_obj <- list(
    basis = fit$basis,
    basis_names = fit$basis_names,
    p = length(fit$basis),
    d = fit$d,
    m = fit$m
  )
  
  # Get Score design from z and basis_obj
  design <- build_sm_design_mv(z, basis_obj)
  
  # Calculate score for each transformed sample z
  score_z <- do.call(
    cbind,
    lapply(seq_len(fit$d), function(k) design$D[[k]] %*% fit$theta)
  )
  
  # Back transformation to x
  score_x <- sweep(score_z, 2, fit$scaling$scale, FUN = "/")
  colnames(score_x) <- paste0("dim", seq_len(fit$d))
  score_x
}
