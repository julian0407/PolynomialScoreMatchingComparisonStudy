
# source("Fourth_Script_(Robust).R")

# source("Fith_Script_(Generalized).R")

# source("Sixth_Script_(Skewed).R")

source("Polynomial_Score_Matching.R")






# ============================================================
# Example run
# ============================================================
set.seed(1)

x <- rnorm(1000, mean = 50)
# x <- rnorm(10000, mean = 0, sd = 1 / sqrt(2))
# x <- rlogis(1000)
# x <- rlaplace(5000)
# x <- rsn(5000, xi = 0, omega = 1, alpha = 5)

# x2 <- sample(x, size = 1000)

# positive examples
# x <- rexp(300, rate = 1)
# x <- rchisq(500, df = 5)
# x <- rgamma(5000, shape = 2, rate = 1)
# x <- rbeta(300, shape1 = 2, shape2 = 3)


m <- 10

# ============================================================
# Robust
# ============================================================
fit <- fit_score_matching_matrixG_robust(
  x, m,
  lambda_trace = 1e-2,
  lambda_frob  = 1e-3,
  eps = 1e-6,
  standardize = TRUE,
  clip = 6,
  col_scale = "maxabs",
  solver = "SCS",
  verbose = TRUE
)

# ============================================================
# Generalized with
# ============================================================
fit <- fit_score_matching_matrixG_robust(
  x, m,
  lambda_trace = 1e-5,
  lambda_frob  = 1e-9,
  eps = 1e-6,
  standardize = TRUE,
  clip = NULL,
  col_scale = "maxabs",
  solver = "SCS",
  verbose = TRUE,
  generalized = TRUE,
  h_cap = 4
)

# ============================================================
# For Skewed with
# ============================================================
fit <- fit_score_matching_yj(
  x, m,
  lambda_yj = 0.85,         # rsn: 0.9; gamma: 0.85
  lambda_trace = 1e-2,
  lambda_frob  = 1e-3,
  eps = 1e-6,
  standardize = TRUE,
  clip = NULL,
  col_scale = "maxabs",
  solver = "SCS",
  verbose = TRUE
)

# ============================================================
# Final
# ============================================================
fit <- fit_score_matching_matrixG_robust(
  x = x,
  m = m,
  h = function(x) h_tanh_sq(x),
  h_prime = function(x) h_tanh_sq_prime(x),
  eps = 1e-6,
  standardize = TRUE
)

xgrid <- seq(min(x) - 1, max(x) + 1, length.out = 500)
p_hat <- density_from_fit(xgrid, fit$G, fit$c1, mu = fit$scaling$mu, s = fit$scaling$s)

p_norm    <- dnorm(xgrid, 50, 1)
p_squared <- dnorm(xgrid, 0, 1 / sqrt(2))
# p_logis   <- dlogis(xgrid)
# p_laplace <- dlaplace(xgrid)
p_skewed_norm <- dsn(xgrid, xi = 0, omega = 1, alpha = 5)

# p_exp   <- dexp(xgrid, rate = 1)
# p_chisq <- dchisq(xgrid, df = 5)
p_gamma <- dgamma(xgrid, shape = 2, rate = 1)

hist(x, breaks = 30, freq = FALSE, col = "grey90", border = "grey70",
     main = "Score-matching log-concave density", xlab = "x")
lines(xgrid, p_hat,  lwd = 2, col = "blue")
lines(xgrid, p_norm, lwd = 2, col = "red", lty = 2)
legend("topright",
       legend = c("Score matching estimate", "Reference"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n")

