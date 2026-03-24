build_pairwise_score_matching <- function(
    x, m,
    build_Phi = TRUE,
    drop_constant = TRUE,
    standardize = TRUE,
    center = NULL,
    scale = NULL
) {
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("x must be numeric.")
  n <- nrow(x)
  d <- ncol(x)
  if (m < 1) stop("m must be >= 1.")
  
  # ----------------------------------------------------------
  # Standardisierung
  # ----------------------------------------------------------
  if (standardize) {
    if (is.null(center)) center <- colMeans(x)
    if (is.null(scale))  scale  <- apply(x, 2, sd)
    
    # konstante / fast konstante Spalten absichern
    scale[!is.finite(scale) | scale <= 0] <- 1
    
    x_work <- sweep(x, 2, center, FUN = "-")
    x_work <- sweep(x_work, 2, scale, FUN = "/")
  } else {
    center <- rep(0, d)
    scale  <- rep(1, d)
    x_work <- x
  }
  
  # ----------------------------------------------------------
  # Basis:
  # univariat: x_l^r, r = 0, ..., 2m-2
  # pairwise:  x_l^i x_u^j, 1 <= l < u <= d, i,j = 0, ..., m-1
  # ----------------------------------------------------------
  uni_deg <- 0:(2 * m - 2)
  pair_deg <- 0:(m - 1)
  
  basis <- vector("list", 0)
  
  for (l in seq_len(d)) {
    for (r in uni_deg) {
      if (drop_constant && r == 0) next
      basis[[length(basis) + 1L]] <- list(
        type = "uni",
        l = l,
        r = r
      )
    }
  }
  
  for (l in seq_len(d - 1L)) {
    for (u in (l + 1L):d) {
      for (i in pair_deg) {
        for (j in pair_deg) {
          if (drop_constant && i == 0 && j == 0) next
          basis[[length(basis) + 1L]] <- list(
            type = "pair",
            l = l,
            u = u,
            i = i,
            j = j
          )
        }
      }
    }
  }
  
  p <- length(basis)
  
  # ----------------------------------------------------------
  # Potenzen auf standardisierten Daten vorkalkulieren
  # ----------------------------------------------------------
  max_deg <- 2 * m - 2
  Xpow <- vector("list", d)
  
  for (l in seq_len(d)) {
    P <- matrix(1, nrow = n, ncol = max_deg + 1L)
    if (max_deg >= 1) {
      P[, 2] <- x_work[, l]
      if (max_deg >= 2) {
        for (r in 2:max_deg) {
          P[, r + 1L] <- P[, r] * x_work[, l]
        }
      }
    }
    Xpow[[l]] <- P
  }
  
  get_pow <- function(l, deg) {
    if (deg < 0) return(rep(0, n))
    Xpow[[l]][, deg + 1L]
  }
  
  Phi <- if (build_Phi) matrix(0, nrow = n, ncol = p) else NULL
  D <- lapply(seq_len(d), function(k) matrix(0, nrow = n, ncol = p))
  Lap <- matrix(0, nrow = n, ncol = p)
  
  for (col in seq_len(p)) {
    b <- basis[[col]]
    
    if (b$type == "uni") {
      l <- b$l
      r <- b$r
      
      if (build_Phi) {
        Phi[, col] <- get_pow(l, r)
      }
      
      if (r >= 1) {
        D[[l]][, col] <- r * get_pow(l, r - 1L)
      }
      
      if (r >= 2) {
        Lap[, col] <- r * (r - 1) * get_pow(l, r - 2L)
      }
      
    } else {
      l <- b$l
      u <- b$u
      i <- b$i
      j <- b$j
      
      xl_i <- get_pow(l, i)
      xu_j <- get_pow(u, j)
      
      if (build_Phi) {
        Phi[, col] <- xl_i * xu_j
      }
      
      if (i >= 1) {
        D[[l]][, col] <- i * get_pow(l, i - 1L) * xu_j
      }
      
      if (j >= 1) {
        D[[u]][, col] <- j * xl_i * get_pow(u, j - 1L)
      }
      
      term_l <- if (i >= 2) i * (i - 1) * get_pow(l, i - 2L) * xu_j else 0
      term_u <- if (j >= 2) j * (j - 1) * xl_i * get_pow(u, j - 2L) else 0
      Lap[, col] <- term_l + term_u
    }
  }
  
  K <- matrix(0, nrow = p, ncol = p)
  for (k in seq_len(d)) {
    K <- K + crossprod(D[[k]])
  }
  K <- K / n
  
  ell <- colMeans(Lap)
  
  basis_names <- vapply(basis, function(b) {
    if (b$type == "uni") {
      sprintf("z%d^%d", b$l, b$r)
    } else {
      sprintf("z%d^%d*z%d^%d", b$l, b$i, b$u, b$j)
    }
  }, character(1))
  
  colnames(K) <- rownames(K) <- basis_names
  names(ell) <- basis_names
  if (!is.null(Phi)) colnames(Phi) <- basis_names
  
  list(
    x_original = x,
    x_used = x_work,
    center = center,
    scale = scale,
    standardized = standardize,
    Phi = Phi,
    K = K,
    ell = ell,
    basis = basis,
    basis_names = basis_names,
    D = D,
    Lap = Lap,
    m = m
  )
}

transform_with_fit <- function(newx, fit_obj) {
  newx <- as.matrix(newx)
  sweep(sweep(newx, 2, fit_obj$center, "-"), 2, fit_obj$scale, "/")
}


# ============================================================
# Unconstrained multivariate pairwise polynomial score matching
# based on your build_pairwise_score_matching()
# ============================================================

fit_pairwise_score_matching_unconstrained <- function(
    x, m,
    standardize = TRUE,
    center = NULL,
    scale = NULL,
    drop_constant = TRUE,
    ridge = 0,
    use_ridge = TRUE,
    ridge_exclude_constant = TRUE,
    build_Phi = FALSE,
    solver = c("solve", "qr")
) {
  solver <- match.arg(solver)
  
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) stop("x must be numeric.")
  if (nrow(x) < 2) stop("Need at least 2 observations.")
  if (ncol(x) < 1) stop("Need at least 1 column.")
  if (m < 1) stop("m must be >= 1.")
  if (!is.numeric(ridge) || length(ridge) != 1 || ridge < 0) {
    stop("ridge must be a nonnegative scalar.")
  }
  
  prep <- build_pairwise_score_matching(
    x = x,
    m = m,
    build_Phi = build_Phi,
    drop_constant = drop_constant,
    standardize = standardize,
    center = center,
    scale = scale
  )
  
  K <- prep$K
  ell <- prep$ell
  p <- ncol(K)
  
  # ----------------------------------------------------------
  # Ridge matrix
  # optionally do not penalize constants
  # ----------------------------------------------------------
  pen_diag <- rep(1, p)
  
  if (ridge_exclude_constant) {
    is_const <- vapply(prep$basis, function(b) {
      if (b$type == "uni") {
        b$r == 0
      } else {
        b$i == 0 && b$j == 0
      }
    }, logical(1))
    pen_diag[is_const] <- 0
  }
  
  lambda_eff <- if (use_ridge) ridge else 0
  K_reg <- K + lambda_eff * diag(pen_diag, nrow = p, ncol = p)
  
  # ----------------------------------------------------------
  # Solve linear system:
  #   (K + lambda I) theta = ell
  # ----------------------------------------------------------
  theta <- switch(
    solver,
    solve = tryCatch(
      solve(K_reg, ell),
      error = function(e) {
        warning("solve() failed; falling back to qr.solve().")
        qr.solve(K_reg, ell)
      }
    ),
    qr = qr.solve(K_reg, ell)
  )
  
  fitted_score_train <- do.call(
    cbind,
    lapply(seq_along(prep$D), function(k) prep$D[[k]] %*% theta)
  )
  
  if (prep$standardized) {
    fitted_score_train <- sweep(fitted_score_train, 2, prep$scale, "/")
  }
  
  colnames(fitted_score_train) <- paste0("dim", seq_len(ncol(x)))
  
  objective_value <- as.numeric(
    0.5 * crossprod(theta, prep$K %*% theta) - crossprod(prep$ell, theta)
  )
  
  objective_value_reg <- as.numeric(
    objective_value + 0.5 * lambda_eff * sum(pen_diag * theta^2)
  )
  
  out <- list(
    theta = as.numeric(theta),
    K = prep$K,
    ell = prep$ell,
    K_reg = K_reg,
    ridge = ridge,
    use_ridge = use_ridge,
    lambda_eff = lambda_eff,
    ridge_penalty_diag = pen_diag,
    objective_value = objective_value,
    objective_value_reg = objective_value_reg,
    fitted_score_train = fitted_score_train,
    prep = prep,
    m = m,
    d = ncol(x),
    n = nrow(x),
    drop_constant = drop_constant,
    standardize = standardize,
    solver = solver,
    call = match.call()
  )
  
  class(out) <- "pairwise_score_matching_fit"
  out
}

# ------------------------------------------------------------
# Evaluate derivative design matrices on new data
# returns only what we need for score prediction
# ------------------------------------------------------------
evaluate_pairwise_derivatives <- function(newx, prep) {
  newx <- as.matrix(newx)
  if (!is.numeric(newx)) stop("newx must be numeric.")
  if (ncol(newx) != ncol(prep$x_original)) {
    stop("Dimension mismatch in newx.")
  }
  
  z <- transform_with_fit(newx, prep)
  n <- nrow(z)
  d <- ncol(z)
  p <- length(prep$basis)
  max_deg <- 2 * prep$m - 2
  
  Xpow <- vector("list", d)
  for (l in seq_len(d)) {
    P <- matrix(1, nrow = n, ncol = max_deg + 1L)
    if (max_deg >= 1) {
      P[, 2] <- z[, l]
      if (max_deg >= 2) {
        for (r in 2:max_deg) {
          P[, r + 1L] <- P[, r] * z[, l]
        }
      }
    }
    Xpow[[l]] <- P
  }
  
  get_pow <- function(l, deg) {
    if (deg < 0) return(rep(0, n))
    Xpow[[l]][, deg + 1L]
  }
  
  D_new <- lapply(seq_len(d), function(k) matrix(0, nrow = n, ncol = p))
  
  for (col in seq_len(p)) {
    b <- prep$basis[[col]]
    
    if (b$type == "uni") {
      l <- b$l
      r <- b$r
      
      if (r >= 1) {
        D_new[[l]][, col] <- r * get_pow(l, r - 1L)
      }
      
    } else {
      l <- b$l
      u <- b$u
      i <- b$i
      j <- b$j
      
      xl_i <- get_pow(l, i)
      xu_j <- get_pow(u, j)
      
      if (i >= 1) {
        D_new[[l]][, col] <- i * get_pow(l, i - 1L) * xu_j
      }
      if (j >= 1) {
        D_new[[u]][, col] <- j * xl_i * get_pow(u, j - 1L)
      }
    }
  }
  
  D_new
}

# ------------------------------------------------------------
# Predict negative gradient of log density:
# r_hat(x) = -∇ log p_hat(x) = ∇ S_hat(x)
# ------------------------------------------------------------
predict_pairwise_score_matching <- function(newx, fit) {
  if (!inherits(fit, "pairwise_score_matching_fit")) {
    stop("fit must be of class 'pairwise_score_matching_fit'.")
  }
  
  newx <- as.matrix(newx)
  D_new <- evaluate_pairwise_derivatives(newx, fit$prep)
  
  out_z <- do.call(
    cbind,
    lapply(seq_along(D_new), function(k) D_new[[k]] %*% fit$theta)
  )
  
  # chain rule back to x-scale
  out_x <- sweep(out_z, 2, fit$prep$scale, "/")
  
  colnames(out_x) <- paste0("dim", seq_len(ncol(out_x)))
  out_x
}

score_loss_direct_pairwise_mv <- function(
    x_test, fit, score_true,
    h = function(x) rep(1, nrow(x))
) {
  x_test <- as.matrix(x_test)
  r_hat  <- predict_pairwise_score_matching(x_test, fit)
  r_true <- score_true(x_test)
  
  if (!is.matrix(r_true)) r_true <- as.matrix(r_true)
  if (!all(dim(r_true) == dim(r_hat))) {
    stop("score_true must return a matrix with same dimension as x_test.")
  }
  
  w <- h(x_test)
  if (length(w) != nrow(x_test)) stop("h must return one weight per row.")
  w <- pmax(as.numeric(w), 0)
  
  sqerr <- rowSums((r_hat - r_true)^2)
  mean(0.5 * w * sqerr)
}

predict_pairwise_energy <- function(newx, fit) {
  newx <- as.matrix(newx)
  prep <- fit$prep
  
  z <- transform_with_fit(newx, prep)
  n <- nrow(z)
  p <- length(prep$basis)
  
  max_deg <- 2 * prep$m - 2
  d <- ncol(z)
  
  Xpow <- vector("list", d)
  for (l in seq_len(d)) {
    P <- matrix(1, nrow = n, ncol = max_deg + 1L)
    if (max_deg >= 1) {
      P[,2] <- z[,l]
      if (max_deg >= 2) {
        for (r in 2:max_deg) {
          P[,r+1] <- P[,r] * z[,l]
        }
      }
    }
    Xpow[[l]] <- P
  }
  
  get_pow <- function(l,deg){
    if(deg<0) return(rep(0,n))
    Xpow[[l]][,deg+1]
  }
  
  S <- numeric(n)
  
  for(col in seq_along(prep$basis)){
    
    b <- prep$basis[[col]]
    
    if(b$type=="uni"){
      S <- S + fit$theta[col] * get_pow(b$l,b$r)
      
    } else {
      S <- S + fit$theta[col] *
        get_pow(b$l,b$i) *
        get_pow(b$u,b$j)
    }
  }
  
  S
}

predict_pairwise_density <- function(newx, fit, mc_points) {
  
  S_new <- predict_pairwise_energy(newx, fit)
  S_mc  <- predict_pairwise_energy(mc_points, fit)
  
  Z_hat <- mean(exp(-S_mc))
  
  exp(-S_new) / Z_hat
}

# ============================================================
# Minimaler Vergleich: pairwise score matching vs KDE
# ============================================================

set.seed(1)

d <- 2
mu <- c(0, 0)
Sigma <- matrix(c(1, 0.4,
                  0.4, 1.5), 2, 2)

n_train <- 10000
n_test  <- 5000

x_train <- rmvnorm_base(n_train, mu, Sigma)
x_test  <- rmvnorm_base(n_test,  mu, Sigma)

# ------------------------------------------------------------
# Dein neuer pairwise score matching estimator
# ------------------------------------------------------------
fit_pair_ridge_small <- fit_pairwise_score_matching_unconstrained(
  x = x_train,
  m = 3,
  standardize = TRUE,
  ridge = 1e-3,
  use_ridge = TRUE
)

fit_pair_ridge_large <- fit_pairwise_score_matching_unconstrained(
  x = x_train,
  m = 3,
  standardize = TRUE,
  ridge = 1e-2,
  use_ridge = TRUE
)

# ------------------------------------------------------------
# KDE baselines
# ------------------------------------------------------------
fit_kde_pi   <- fit_kde_mv(x_train, H_method = "Hpi", diagonal = FALSE)
fit_kde_lscv <- fit_kde_mv(x_train, H_method = "Hlscv", diagonal = FALSE)

# ------------------------------------------------------------
# True score for comparison
# ------------------------------------------------------------
score_true_fun <- function(x) score_mvn_true(x, mu, Sigma)

# ------------------------------------------------------------
# Direct score-loss
# ------------------------------------------------------------
sl_pair_ridge_small <- score_loss_direct_pairwise_mv(
  x_test, fit_pair_ridge_small,
  score_true = score_true_fun
)

sl_pair_ridge_large <- score_loss_direct_pairwise_mv(
  x_test, fit_pair_ridge_large,
  score_true = score_true_fun
)


sl_kde_pi <- score_loss_direct_kde_mv(
  x_test, fit_kde_pi,
  score_true = score_true_fun
)

sl_kde_lscv <- score_loss_direct_kde_mv(
  x_test, fit_kde_lscv,
  score_true = score_true_fun
)

cat("Direct score-loss pairwise SM (ridge=1e-3): ", sl_pair_ridge_small, "\n")
cat("Direct score-loss pairwise SM (ridge=1e-2): ", sl_pair_ridge_large, "\n")
cat("Direct score-loss KDE (Hpi):                ", sl_kde_pi, "\n")
cat("Direct score-loss KDE (Hlscv):              ", sl_kde_lscv, "\n")

# ------------------------------------------------------------
# Optional: inspect first few score predictions
# ------------------------------------------------------------
head(
  cbind(
    x_test[1:6, ],
    predict_pairwise_score_matching(x_test[1:6, ], fit_pair_ridge_small),
    score_true_fun(x_test[1:6, ])
  )
)

mc_points <- rmvnorm_base(20000, mu, Sigma)

p_pair <- predict_pairwise_density(
  x_test,
  fit_pair_ridge_small,
  mc_points
)

p_kde <- predict_kde_density_mv(x_test, fit_kde_pi)

p_true <- dmvnorm_base(x_test, mu, Sigma)

L1_pair <- mean(abs(p_pair - p_true))
L1_kde  <- mean(abs(p_kde - p_true))

cat("L1 error pairwise SM: ", L1_pair, "\n")
cat("L1 error KDE:         ", L1_kde, "\n")
