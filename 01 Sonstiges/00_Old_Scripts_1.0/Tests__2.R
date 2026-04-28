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

# -------- Robust Drton-style generalized score matching on R --------
fit_score_matching_logconcave_R <- function(
    x,
    m,
    h       = function(z) rep(1, length(z)),
    h_prime = function(z) rep(0, length(z)),
    standardize = TRUE,
    col_scale = c("sd", "maxabs", "none"),
    delta_lead = 1e-3,
    ridge = 1e-6,
    solver = c("SCS", "MOSEK"),
    scs_control = list(max_iters = 200000, eps = 1e-5, alpha = 1.8, verbose = FALSE),
    retry_if_inaccurate = TRUE
) {
  col_scale <- match.arg(col_scale)
  solver <- match.arg(solver)
  
  sc <- list(mu = 0, s = 1)
  z <- x
  if (standardize) {
    sc <- scale_x(x)
    z <- sc$z
  }
  
  n <- length(z)
  p <- m * m
  
  A <- matrix(0, nrow = n, ncol = p)
  B <- matrix(0, nrow = n, ncol = p)
  for (k in seq_len(n)) {
    MN <- make_MN(z[k], m)
    A[k, ] <- as.vector(MN$M)
    B[k, ] <- as.vector(MN$N)
  }
  
  hx  <- as.numeric(h(z))
  hpx <- as.numeric(h_prime(z))
  stopifnot(length(hx) == n, length(hpx) == n)
  hx <- pmax(hx, 0)
  
  b_bar <- colMeans(B * hx)
  
  if (col_scale == "none") {
    scale_vec <- rep(1, p)
  } else if (col_scale == "maxabs") {
    scale_vec <- apply(abs(A), 2, function(v) max(v, 1e-12))
  } else {
    scale_vec <- apply(A, 2, function(v) max(sd(v), 1e-12))
  }
  A_sc <- sweep(A, 2, scale_vec, "/")
  b_sc <- b_bar / scale_vec
  
  G  <- CVXR::Variable(m, m, PSD = TRUE)
  c1 <- CVXR::Variable(1)
  gvec <- tryCatch(CVXR::vec(G), error = function(e) CVXR::reshape(G, c(p, 1)))
  
  A_c  <- CVXR::Constant(A_sc)
  ones <- CVXR::Constant(matrix(1, n, 1))
  w_c  <- CVXR::Constant(matrix(sqrt(pmax(hx, 1e-12)), n, 1))
  hp_c <- CVXR::Constant(matrix(hpx, n, 1))
  b_c  <- CVXR::Constant(matrix(b_sc, p, 1))
  
  s1 <- A_c %*% gvec + c1 * ones
  
  obj <- (0.5 / n) * CVXR::sum_squares(CVXR::multiply(w_c, s1)) -
    (1 / n) * t(hp_c) %*% s1 -
    t(b_c) %*% gvec
  
  obj <- obj + (ridge / 2) * (CVXR::sum_squares(gvec) + CVXR::sum_squares(c1))
  
  constr <- list(
    G == t(G),
    G[m, m] >= delta_lead
  )
  
  prob <- CVXR::Problem(CVXR::Minimize(obj), constraints = constr)
  
  solve_one <- function(eps_override = NULL, it_override = NULL) {
    if (solver == "MOSEK") {
      CVXR::solve(prob, solver = "MOSEK")
    } else {
      ctl <- scs_control
      if (!is.null(eps_override)) ctl$eps <- eps_override
      if (!is.null(it_override)) ctl$max_iters <- it_override
      
      # explicit call avoids base::solve()
      CVXR::solve(
        prob, solver = "SCS",
        max_iters = ctl$max_iters,
        eps = ctl$eps,
        alpha = ctl$alpha,
        verbose = ctl$verbose
      )
    }
  }
  
  sol <- solve_one()
  
  if (retry_if_inaccurate && identical(sol$status, "optimal_inaccurate") && solver == "SCS") {
    sol <- solve_one(eps_override = 1e-6, it_override = 400000)
  }
  
  g_sc  <- as.numeric(sol$getValue(gvec))
  g_org <- g_sc / scale_vec
  G_org <- matrix(g_org, nrow = m, ncol = m)
  
  list(
    G = G_org,
    c1 = as.numeric(sol$getValue(c1)),
    status = sol$status,
    scaling = list(mu = sc$mu, s = sc$s, standardize = standardize,
                   scale_vec = scale_vec, col_scale = col_scale),
    hyper = list(delta_lead = delta_lead, ridge = ridge, solver = solver),
    solution = sol
  )
}

# reconstruction of s(z)
s_function <- function(z, G, c1) {
  m <- nrow(G)
  val <- 0
  for (i in 0:(m - 1)) {
    for (j in 0:(m - 1)) {
      val <- val +
        G[i + 1, j + 1] * z^(i + j + 2) / ((i + j + 1) * (i + j + 2))
    }
  }
  val + c1 * z
}

stable_exp_neg <- function(v) {
  shift <- min(v)
  list(values = exp(-(v - shift)), shift = shift)
}

density_from_fit <- function(xgrid, fit) {
  G <- fit$G
  c1 <- fit$c1
  mu <- fit$scaling$mu
  s  <- fit$scaling$s
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)
  
  zgrid <- (xgrid - mu) / s
  svals <- sapply(zgrid, s_function, G = G, c1 = c1)
  
  tmp <- stable_exp_neg(svals)
  unnorm <- tmp$values
  Z <- pracma::trapz(zgrid, unnorm)
  pz <- unnorm / Z
  (1 / s) * pz
}

density_from_fit_m1_closedform <- function(xgrid, fit) {
  G11 <- fit$G[1, 1]
  c1  <- fit$c1
  mu  <- fit$scaling$mu
  s   <- fit$scaling$s
  
  mu_z <- -c1 / G11
  sd_z <- 1 / sqrt(G11)
  
  mu_x <- mu + s * mu_z
  sd_x <- s * sd_z
  
  dnorm(xgrid, mean = mu_x, sd = sd_x)
}

L1_error_gauss <- function(mu1, sd1, mu2, sd2) {
  sd_max <- max(sd1, sd2)
  lo <- min(mu1, mu2) - 12 * sd_max
  hi <- max(mu1, mu2) + 12 * sd_max
  xg <- seq(lo, hi, length.out = 20001)
  f1 <- dnorm(xg, mu1, sd1)
  f2 <- dnorm(xg, mu2, sd2)
  pracma::trapz(xg, abs(f1 - f2))
}

fit_score_matching_m1_closedform <- function(
    x,
    standardize = TRUE,
    ridge = 0,          # set e.g. 1e-8 if you want
    delta_lead = 1e-8   # enforce G >= delta
) {
  sc <- list(mu = 0, s = 1)
  z <- x
  if (standardize) {
    sc$mu <- mean(x)
    sc$s  <- sd(x)
    if (!is.finite(sc$s) || sc$s <= 0) sc$s <- 1
    z <- (x - sc$mu) / sc$s
  }
  
  n <- length(z)
  m1 <- mean(z)
  m2 <- mean(z^2)
  
  # Normal equations for minimizing 0.5*E[(G z + c)^2] - G + (ridge/2)(G^2 + c^2)
  # Gradient:
  # (m2+ridge) G + m1 c = 1
  # m1 G + (1+ridge) c = 0
  A <- matrix(c(m2 + ridge, m1,
                m1, 1 + ridge), 2, 2, byrow = TRUE)
  b <- c(1, 0)
  sol <- solve(A, b)
  Ghat <- sol[1]
  c1hat <- sol[2]
  
  # enforce constraints (should be unnecessary for normal, but safe)
  if (Ghat < delta_lead) Ghat <- delta_lead
  
  # implied Gaussian on x-scale
  mu_z <- -c1hat / Ghat
  sd_z <- 1 / sqrt(Ghat)
  mu_x <- sc$mu + sc$s * mu_z
  sd_x <- sc$s * sd_z
  
  list(
    G = matrix(Ghat, 1, 1),
    c1 = c1hat,
    scaling = list(mu = sc$mu, s = sc$s, standardize = standardize),
    implied = list(mu_x = mu_x, sd_x = sd_x, mu_z = mu_z, sd_z = sd_z)
  )
}

density_from_fit_m1_closedform <- function(xgrid, fit) {
  dnorm(xgrid, mean = fit$implied$mu_x, sd = fit$implied$sd_x)
}


set.seed(1)
mu0 <- 50; sd0 <- 10
ns <- c(200, 500, 1000, 2000, 5000, 10000, 20000)
R <- 30

L1_error_gauss <- function(mu1, sd1, mu2, sd2) {
  sd_max <- max(sd1, sd2)
  lo <- min(mu1, mu2) - 12 * sd_max
  hi <- max(mu1, mu2) + 12 * sd_max
  xg <- seq(lo, hi, length.out = 20001)
  f1 <- dnorm(xg, mu1, sd1)
  f2 <- dnorm(xg, mu2, sd2)
  pracma::trapz(xg, abs(f1 - f2))
}

# Speicher für Fehler
errs <- matrix(NA_real_, nrow = length(ns), ncol = R)

# Speicher für eine repräsentative Stichprobe + Fit pro n (für Plot)
x_store   <- vector("list", length(ns))
fit_store <- vector("list", length(ns))

for (i in seq_along(ns)) {
  n <- ns[i]
  
  # Eine feste Stichprobe speichern (für den Density-Plot)
  x_rep <- rnorm(n, mu0, sd0)
  x_store[[i]] <- x_rep
  fit_store[[i]] <- fit_score_matching_m1_closedform(
    x_rep, standardize = TRUE, ridge = 0
  )
  
  # L1 über R Wiederholungen (statistische Stabilität)
  for (r in 1:R) {
    x <- rnorm(n, mu0, sd0)
    fit <- fit_score_matching_m1_closedform(x, standardize = TRUE, ridge = 0)
    
    errs[i, r] <- L1_error_gauss(
      fit$implied$mu_x, fit$implied$sd_x,
      mu0, sd0
    )
  }
  
  cat("n=", n, " median L1=", median(errs[i, ]), "\n")
}

# =========================
# Plot 1: L1-Konvergenz
# =========================
plot(ns, apply(errs, 1, median), log="x", type="b",
     pch = 19, col = "blue",
     xlab = "Sample size n",
     ylab = "Median L1 error",
     main = "L1-Konvergenz des Score-Matching Schätzers (m=1)")

grid()

# =========================
# Plot 2: Density für jedes n
# =========================
# Layout: mehrere Panels
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))

for (i in seq_along(ns)) {
  n   <- ns[i]
  x   <- x_store[[i]]
  fit <- fit_store[[i]]
  
  # adaptives Grid (wichtig für faire Darstellung der Tails)
  sd_hat <- fit$implied$sd_x
  mu_hat <- fit$implied$mu_x
  sd_max <- max(sd_hat, sd0)
  
  xgrid <- seq(
    min(mu_hat, mu0) - 4 * sd_max,
    max(mu_hat, mu0) + 4 * sd_max,
    length.out = 600
  )
  
  # Geschätzte Dichte (closed-form, statistisch korrekt)
  p_hat <- dnorm(xgrid, mean = mu_hat, sd = sd_hat)
  p_ref <- dnorm(xgrid, mean = mu0, sd = sd0)
  
  hist(x, breaks = 30, freq = FALSE,
       col = "grey90", border = "grey70",
       main = paste0("n = ", n),
       xlab = "", ylab = "")
  
  lines(xgrid, p_hat, lwd = 2, col = "darkgreen")
  lines(xgrid, p_ref, lwd = 2, col = "red", lty = 2)
  
  legend("topright",
         legend = c("Score Matching (m=1)", "True N(50,10)"),
         col = c("darkgreen", "red"),
         lwd = 2, lty = c(1, 2), bty = "n", cex = 0.8)
}

# Reset Layout
par(mfrow = c(1, 1))