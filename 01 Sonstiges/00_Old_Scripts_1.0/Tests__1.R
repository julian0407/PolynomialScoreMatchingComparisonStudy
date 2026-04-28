library(CVXR)
library(pracma)

make_MN <- function(x, m) {
  i0 <- 0:(m - 1)
  j0 <- 0:(m - 1)
  I  <- matrix(rep(i0, times = m), nrow = m, byrow = FALSE)
  J  <- matrix(rep(j0, each  = m), nrow = m, byrow = FALSE)
  N <- x^(I + J)                     # basis for s''(x)
  M <- x^(I + J + 1) / (I + J + 1)   # integrated basis for s'(x)
  list(M = M, N = N)
}

scale_x <- function(x) {
  mu <- mean(x)
  s  <- sd(x)
  if (!is.finite(s) || s <= 0) s <- 1
  z <- (x - mu) / s
  list(z = z, mu = mu, s = s)
}

fit_score_matching_matrixG_R <- function(
    x,
    m,
    h       = function(x) rep(1, length(x)),
    h_prime = function(x) rep(0, length(x)),
    standardize = TRUE,
    # coercivity + ridge amplifier:
    delta_lead = 1e-6,     # enforce G[m,m] >= delta_lead
    ridge      = 1e-8,     # diagonal amplifier (ridge on theta)
    # numerical:
    col_scale = c("maxabs", "sd", "none"),
    solver = "SCS"
) {
  col_scale <- match.arg(col_scale)
  
  sc <- list(mu = 0, s = 1)
  x_used <- x
  if (standardize) {
    sc <- scale_x(x)
    x_used <- sc$z
  }
  
  n <- length(x_used)
  p <- m * m
  
  # Precompute design matrices A (for s') and B (for s'')
  A <- matrix(0, nrow = n, ncol = p)
  B <- matrix(0, nrow = n, ncol = p)
  for (k in seq_len(n)) {
    MN <- make_MN(x_used[k], m)
    A[k, ] <- as.vector(MN$M)
    B[k, ] <- as.vector(MN$N)
  }
  
  hx  <- as.numeric(h(x_used))
  hpx <- as.numeric(h_prime(x_used))
  stopifnot(length(hx) == n, length(hpx) == n)
  
  # b_bar = E_n[h(X) * basis_for_s'']
  B_h   <- B * hx
  b_bar <- colMeans(B_h)  # length p
  
  # Column scaling for numerical conditioning
  if (col_scale == "none") {
    scale_vec <- rep(1, p)
  } else if (col_scale == "maxabs") {
    scale_vec <- apply(abs(A), 2, function(v) max(v, 1e-12))
  } else {
    scale_vec <- apply(A, 2, function(v) max(sd(v), 1e-12))
  }
  A_scaled    <- sweep(A, 2, scale_vec, "/")
  bbar_scaled <- b_bar / scale_vec
  
  # Build CVXR problem
  G  <- Variable(m, m, PSD = TRUE)
  c1 <- Variable(1)
  
  # vec(G) (robust across CVXR versions)
  gvec <- tryCatch(vec(G), error = function(e) reshape(G, c(p, 1)))
  
  A_c   <- Constant(A_scaled)                             # n x p
  ones  <- Constant(matrix(1, nrow = n, ncol = 1))        # n x 1
  w_c   <- Constant(matrix(sqrt(pmax(hx, 1e-12)), n, 1))  # n x 1
  hp_c  <- Constant(matrix(hpx, n, 1))                    # n x 1
  b_c   <- Constant(matrix(bbar_scaled, p, 1))            # p x 1
  
  # s'(x) = A g + c1
  s1 <- A_c %*% gvec + c1 * ones
  
  # generalized score matching objective (up to constants)
  obj_base <- (0.5 / n) * sum_squares(multiply(w_c, s1)) -
    (1 / n) * t(hp_c) %*% s1 -
    t(b_c) %*% gvec
  
  # ridge / diagonal amplifier -> strong convexity in theta
  obj_ridge <- (ridge / 2) * (sum_squares(gvec) + sum_squares(c1))
  
  # Constraints:
  # - symmetry (numerical)
  # - coercive tails: leading coefficient positive
  constr <- list(
    G == t(G),
    G[m, m] >= delta_lead
  )
  
  prob <- Problem(Minimize(obj_base + obj_ridge), constraints = constr)
  sol  <- solve(prob, solver = solver)
  
  # Back-transform from column scaling
  gvec_sol  <- as.numeric(sol$getValue(gvec))
  gvec_orig <- gvec_sol / scale_vec
  G_orig    <- matrix(gvec_orig, nrow = m, ncol = m)
  
  list(
    G  = G_orig,
    c1 = as.numeric(sol$getValue(c1)),
    solution = sol,
    scaling  = list(
      standardize = standardize,
      mu = sc$mu, s = sc$s,
      col_scale = col_scale,
      scale_vec = scale_vec
    ),
    hyper = list(delta_lead = delta_lead, ridge = ridge)
  )
}

# reconstruct s(x) from G and c1 on standardized scale
s_function <- function(x, G, c1) {
  m <- nrow(G)
  val <- 0
  for (i in 0:(m - 1)) {
    for (j in 0:(m - 1)) {
      val <- val +
        G[i + 1, j + 1] *
        x^(i + j + 2) / ((i + j + 1) * (i + j + 2))
    }
  }
  val + c1 * x
}

stable_exp_neg <- function(v) {
  shift <- min(v)
  list(values = exp(-(v - shift)), shift = shift)
}

density_from_fit <- function(xgrid, G, c1, mu = 0, s = 1) {
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)
  zgrid <- (xgrid - mu) / s
  svals <- sapply(zgrid, s_function, G = G, c1 = c1)
  tmp <- stable_exp_neg(svals)
  unnorm <- tmp$values
  Z <- pracma::trapz(zgrid, unnorm)
  pz <- unnorm / Z
  (1 / s) * pz
}

h_tanh_sq <- function(z) {
  tau <- log(length(z))
  u <- abs(z) / tau
  tanh(u)^2
}
h_tanh_sq_prime <- function(z) {
  tau <- log(length(z))
  u <- abs(z) / tau
  sech2 <- 1 / cosh(u)^2
  2 * tanh(u) * sech2 * (sign(z) / tau)
}

set.seed(1)
x <- rnorm(100000, mean = 50, sd = 10)
m <- 5

fit <- fit_score_matching_matrixG_R(
  x = x, m = m,
  h = h_tanh_sq,
  h_prime = h_tanh_sq_prime,
  standardize = TRUE,
  delta_lead = 1e-6,
  ridge = 1e-8
)

xgrid <- seq(min(x) - 1, max(x) + 1, length.out = 500)
p_hat <- density_from_fit(xgrid, fit$G, fit$c1, mu = fit$scaling$mu, s = fit$scaling$s)
p_ref <- dnorm(xgrid, 50, 10)

hist(x, breaks = 30, freq = FALSE, col = "grey90", border = "grey70")
lines(xgrid, p_hat, lwd = 2, col = "yellow")
lines(xgrid, p_ref, lwd = 2, col = "red", lty = 2)

