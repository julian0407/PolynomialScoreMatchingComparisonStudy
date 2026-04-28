library(CVXR)
library(pracma)
library(extraDistr)
library(sn)

# ============================================================
# Robustness helpers: scaling + stable exp normalization
# ============================================================

scale_x <- function(x, clip = NULL) {
  mu <- mean(x)
  s  <- sd(x)
  if (!is.finite(s) || s <= 0) s <- 1
  z <- (x - mu) / s
  if (!is.null(clip)) z <- pmax(pmin(z, clip), -clip)
  list(z = z, mu = mu, s = s, clip = clip)
}

stable_exp_neg <- function(v) {
  vf <- v[is.finite(v)]
  if (length(vf) == 0) {
    # alles Inf/NaN -> gib Nullen zurück (caller muss damit umgehen)
    return(list(values = rep(0, length(v)), shift = 0))
  }
  shift <- min(vf)
  
  # exp(-(v - shift)), aber:
  #   v=+Inf -> exp(-Inf)=0
  #   v=NaN  -> 0
  out <- rep(0, length(v))
  idx <- is.finite(v)
  out[idx] <- exp(-(v[idx] - shift))
  
  list(values = out, shift = shift)
}

# ============================================================
# Generalized score matching helper: bounded, smooth, INCREASING h(|z|)
#   h(z)  = tanh(|z|/tau)^2   in [0,1], increases with |z|, saturates
#   h'(z) = 2*tanh(u)*sech(u)^2 * sign(z)/tau
# ============================================================

h_tanh_sq <- function(z, cap = 2) {
  tau <- cap
  u <- abs(z) / tau
  h  <- tanh(u)^2
  sech2 <- 1 / cosh(u)^2
  hp <- 2 * tanh(u) * sech2 * (sign(z) / tau)
  list(h = h, hp = hp)
}

# -----------------------------
# Build M(x) and N(x) matrices
# -----------------------------
make_MN <- function(x, m) {
  i0 <- 0:(m - 1)
  j0 <- 0:(m - 1)
  I  <- matrix(rep(i0, times = m), nrow = m, byrow = FALSE)
  J  <- matrix(rep(j0, each  = m), nrow = m, byrow = FALSE)
  
  N <- x^(I + J)
  M <- x^(I + J + 1) / (I + J + 1)
  
  list(M = M, N = N)
}

# ============================================================
# Robust (Generalized) score matching fit
#   Generalized branch:
#     - uses bounded increasing h(|z|)
#     - includes + (1/n) * h'(z_i) * s1_i term  (SIGN FIX)
#     - includes c0 (intercept) and c1 (linear term in z) explicitly
# ============================================================

fit_score_matching_matrixG_robust <- function(
    x,
    m,
    lambda_trace = 1e-2,
    lambda_frob  = 0,      # kept but off by default
    eps = 1e-8,
    standardize = TRUE,
    clip = NULL,           # IMPORTANT: default no clipping for tail behavior
    col_scale = c("maxabs", "sd"),
    solver = "SCS",
    verbose = TRUE,
    generalized = TRUE,
    h_cap = 4
) {
  col_scale <- match.arg(col_scale)
  
  # ---- (1) standardize x to avoid huge powers ----
  sc <- list(mu = 0, s = 1, clip = NULL)
  x_used <- x
  if (standardize) {
    sc <- scale_x(x, clip = clip)
    x_used <- sc$z
  }
  
  n <- length(x_used)
  p <- m * m
  
  # ---------- Precompute A and B ----------
  t_pre <- system.time({
    A <- matrix(0, nrow = n, ncol = p)
    B <- matrix(0, nrow = n, ncol = p)
    
    for (k in seq_len(n)) {
      MN <- make_MN(x_used[k], m)
      A[k, ] <- as.vector(MN$M)  # a_k = vec(Mk)
      B[k, ] <- as.vector(MN$N)  # b_k = vec(Nk)
    }
    
    if (generalized) {
      hh <- h_tanh_sq(x_used, cap = h_cap)
      h_vec  <- hh$h
      hp_vec <- hh$hp
      
      # b_bar = (1/n) sum_i h(x_i) b_i
      B_h   <- B * h_vec
      b_bar <- colMeans(B_h)
      
      h_sqrt <- sqrt(pmax(h_vec, 1e-12))
      hp_num <- hp_vec
    } else {
      b_bar  <- colMeans(B)
      h_sqrt <- rep(1, n)
      hp_num <- rep(0, n)
    }
  })
  
  # ---------- (2) Column scaling (preconditioning) ----------
  scale_vec <- switch(
    col_scale,
    maxabs = apply(abs(A), 2, function(v) max(v, 1e-12)),
    sd     = apply(A, 2, function(v) max(sd(v), 1e-12))
  )
  
  A_scaled    <- sweep(A, 2, scale_vec, "/")
  bbar_scaled <- b_bar / scale_vec
  
  # ---------- Build CVXR problem ----------
  t_build <- system.time({
    G  <- CVXR::Variable(m, m, PSD = TRUE)
    
    # intercept + linear term (in standardized coordinate z)
    c0 <- CVXR::Variable(1)
    c1 <- CVXR::Variable(1)
    
    # vec(G) robust across CVXR versions
    gvec <- tryCatch(CVXR::vec(G), error = function(e) CVXR::reshape(G, c(p, 1)))
    
    A_c   <- CVXR::Constant(A_scaled)                     # n x p
    ones  <- CVXR::Constant(matrix(1, nrow = n, ncol = 1)) # n x 1
    zcol  <- CVXR::Constant(matrix(x_used, ncol = 1))      # n x 1
    b_c   <- CVXR::Constant(matrix(bbar_scaled, ncol = 1)) # p x 1
    
    # s1_i = (A*g)_i + c0 + c1*z_i
    s1 <- A_c %*% gvec + c0 * ones + c1 * zcol
    
    w_c  <- CVXR::Constant(matrix(h_sqrt, ncol = 1))  # sqrt(h(z_i))
    hp_c <- CVXR::Constant(matrix(hp_num,  ncol = 1)) # h'(z_i)
    
    if (generalized) {
      obj <- (0.5 / n) * CVXR::sum_squares(CVXR::multiply(w_c, s1)) -
        t(b_c) %*% gvec +
        (1 / n) * t(hp_c) %*% s1 +
        lambda_trace * sum(CVXR::diag(G))     # <-- FIX HERE
    } else {
      obj <- (0.5 / n) * CVXR::sum_squares(s1) -
        t(b_c) %*% gvec +
        lambda_trace * sum(CVXR::diag(G))     # <-- FIX HERE
    }
    
    if (lambda_frob > 0) {
      obj <- obj + lambda_frob * CVXR::sum_squares(G)
    }
    
    prob <- CVXR::Problem(
      CVXR::Minimize(obj),
      constraints = list(
        G == t(G),
        CVXR::diag(G) >= eps
      )
    )
  })
  
  # ---------- Solve ----------
  t_solve <- system.time({
    sol <- CVXR::solve(prob, solver = solver)
  })
  
  # ---------- Back-transform G because we solved with scaled columns ----------
  gvec_sol  <- as.numeric(sol$getValue(gvec))
  gvec_orig <- gvec_sol / scale_vec
  G_orig <- matrix(gvec_orig, nrow = m, ncol = m)
  G_orig <- 0.5 * (G_orig + t(G_orig))
  
  if (verbose) {
    cat("\n--- Timing ---\n")
    cat(sprintf("Precompute (A,B):   %.3f sec\n", t_pre[["elapsed"]]))
    cat(sprintf("Build (CVXR):       %.3f sec\n", t_build[["elapsed"]]))
    cat(sprintf("Solve (solver):     %.3f sec\n", t_solve[["elapsed"]]))
    cat("Status:", sol$status, "\n")
    cat("--------------\n\n")
    
    if (standardize) {
      cat(sprintf("Standardization: mu=%.4f, sd=%.4f, clip=%s\n",
                  sc$mu, sc$s, ifelse(is.null(sc$clip), "NULL", as.character(sc$clip))))
    }
    cat(sprintf("Regularization: lambda_trace=%.2e, lambda_frob=%.2e, eps=%.2e\n",
                lambda_trace, lambda_frob, eps))
    cat(sprintf("Generalized: %s (h_cap=%.3f)\n", generalized, h_cap))
  }
  
  list(
    G  = G_orig,
    c0 = as.numeric(sol$getValue(c0)),
    c1 = as.numeric(sol$getValue(c1)),
    solution = sol,
    timing = list(pre = t_pre, build = t_build, solve = t_solve),
    scaling = list(
      standardize = standardize, mu = sc$mu, s = sc$s, clip = sc$clip,
      col_scale = col_scale, scale_vec = scale_vec
    ),
    generalized = list(enabled = generalized, h_cap = h_cap)
  )
}

# -----------------------------------
# Reconstruct s(z) from G and (c0,c1)  in standardized coordinate z
# -----------------------------------
s_function <- function(z, G, c0, c1) {
  m <- nrow(G)
  
  # Iterativ Potenzen aufbauen: p[k] = z^k, k=0..(2m)
  # Wir brauchen bis exponent = (2m) mindestens, da i+j+2 <= 2m
  max_pow <- 2 * m
  pows <- numeric(max_pow + 1)
  pows[1] <- 1  # z^0
  for (k in 2:(max_pow + 1)) {
    pows[k] <- pows[k - 1] * z
    # frühes Abbrechen bei Overflow
    if (!is.finite(pows[k])) return(Inf)
  }
  
  val <- 0
  for (i in 0:(m - 1)) {
    for (j in 0:(m - 1)) {
      expo <- i + j + 2  # exponent for z^(i+j+2)
      denom <- (i + j + 1) * (i + j + 2)
      term <- G[i + 1, j + 1] * pows[expo + 1] / denom
      val <- val + term
      if (!is.finite(val)) return(Inf)
    }
  }
  
  out <- val + c0 + c1 * z
  if (!is.finite(out)) Inf else out
}

# -----------------------------------
# Density based on score matching fit (stable + inverse transform)
# If fit was on z = (x - mu) / s, then p_X(x) = (1/s) * p_Z(z)
# -----------------------------------
density_from_fit <- function(xgrid, G, c0, c1, mu = 0, s = 1) {
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)
  
  zgrid <- (xgrid - mu) / s
  svals <- vapply(zgrid, s_function, numeric(1), G = G, c0 = c0, c1 = c1)
  
  tmp <- stable_exp_neg(svals)
  unnorm <- tmp$values
  
  # Sicherheitschecks: ersetze NA durch 0
  unnorm[!is.finite(unnorm)] <- 0
  
  Z <- pracma::trapz(zgrid, unnorm)
  
  if (!is.finite(Z) || Z <= 0) {
    # Diagnose helfen: gib sinnvolle Fehlermeldung statt NA
    stop(
      "Normalization failed: integral Z is non-finite or <= 0.\n",
      "Likely cause: s(z) overflow (Inf) on your grid. Try smaller grid range,\n",
      "or reduce polynomial degree m, or keep zgrid within e.g. [-6,6].\n",
      "Quick check: range(svals, finite=TRUE) = ",
      paste(range(svals[is.finite(svals)]), collapse = ", ")
    )
  }
  
  pz <- unnorm / Z
  (1 / s) * pz
}