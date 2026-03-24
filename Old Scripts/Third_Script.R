library(CVXR)
library(pracma)

library(extraDistr)

# -----------------------------
# Build M(x) and N(x) matrices
# -----------------------------
make_MN <- function(x, m) {
  i0 <- 0:(m-1)
  j0 <- 0:(m-1)
  I  <- matrix(rep(i0, times = m), nrow = m, byrow = FALSE)
  J  <- matrix(rep(j0, each  = m), nrow = m, byrow = FALSE)
  
  N <- x^(I + J)
  M <- x^(I + J + 1) / (I + J + 1)
  
  list(M = M, N = N)
}

# -----------------------------
# Score matching objective (vectorized) + trace regularization
# -----------------------------
fit_score_matching_matrixG_fast <- function(x, m, lambda_trace = 0, eps = 1e-6,
                                            solver = "SCS", verbose = TRUE) {
  n <- length(x)
  p <- m * m
  
  # ---------- Precompute A und b_bar in Base R ----------
  t_pre <- system.time({
    A <- matrix(0, nrow = n, ncol = p)
    B <- matrix(0, nrow = n, ncol = p)
    
    for (k in seq_len(n)) {
      MN <- make_MN(x[k], m)
      A[k, ] <- as.vector(MN$M)  # a_k = vec(Mk)
      B[k, ] <- as.vector(MN$N)  # b_k = vec(Nk)
    }
    
    b_bar <- colMeans(B)         # (1/n) sum_k b_k
  })
  
  # ---------- Build small CVXR problem ----------
  t_build <- system.time({
    G  <- Variable(m, m, PSD = TRUE)
    c1 <- Variable(1)
    
    # vec(G) (versionsrobust)
    gvec <- tryCatch(vec(G), error = function(e) reshape(G, c(p, 1)))
    
    A_c   <- Constant(A)                                 # n x p
    ones  <- Constant(matrix(1, nrow = n, ncol = 1))      # n x 1
    b_c   <- Constant(matrix(b_bar, ncol = 1))            # p x 1
    
    # s1 = A %*% gvec + c1  (c1 als Skalar wird ggf. nicht gebroadcastet -> c1*ones)
    s1 <- A_c %*% gvec + c1 * ones
    
    # Objective: (0.5/n) * sum_squares(s1) - b_bar' gvec + lambda * tr(G)
    obj <- (0.5 / n) * sum_squares(s1) - t(b_c) %*% gvec + lambda_trace * sum(diag(G))
    
    prob <- Problem(Minimize(obj), constraints = list(G[m, m] >= eps))
  })
  
  # ---------- Solve ----------
  t_solve <- system.time({
    sol <- solve(prob, solver = solver)
  })
  
  if (verbose) {
    cat("\n--- Timing ---\n")
    cat(sprintf("Precompute (A,B):   %.3f sec\n", t_pre[["elapsed"]]))
    cat(sprintf("Build (CVXR):       %.3f sec\n", t_build[["elapsed"]]))
    cat(sprintf("Solve (solver):     %.3f sec\n", t_solve[["elapsed"]]))
    cat("Status:", sol$status, "\n")
    cat("--------------\n\n")
  }
  
  list(
    G = sol$getValue(G),
    c1 = as.numeric(sol$getValue(c1)),
    solution = sol,
    timing = list(pre = t_pre, build = t_build, solve = t_solve)
  )
}





# -----------------------------------
# Reconstruct s(x) from G and c1
# -----------------------------------
s_function <- function(x, G, c1) {
  m <- nrow(G)
  val <- 0
  for (i in 0:(m-1)) {
    for (j in 0:(m-1)) {
      val <- val +
        G[i+1, j+1] *
        x^(i + j + 2) / ((i + j + 1) * (i + j + 2))
    }
  }
  val + c1 * x
}


# -----------------------------------
# Density based on score matching fit
# -----------------------------------
density_from_fit <- function(xgrid, G, c1) {
  if (is.null(dim(G))) G <- matrix(G, ncol = 1)  # <- wichtig für m=1
  svals <- sapply(xgrid, s_function, G = G, c1 = c1)
  unnorm <- exp(-svals)
  Z <- pracma::trapz(xgrid, unnorm)
  unnorm / Z
}



# -----------------------------
# Example run
# -----------------------------
set.seed(1)


x <- rnorm(1000, mean = 0)
x <- rnorm(300, mean = 0, sd = 1 / sqrt(2))
x <- rlogis(1000)
x <- rlaplace(300)


# -------------------------
# Just positive denisties+
# ------------------------
# x <- rexp(300, rate = 1)
x <- rchisq(5000, df = 5)
x <- rgamma(3000, shape = 2, rate = 1)
# x <- rbeta(300, shape1 = 2, shape2 = 3)

m <- 3

fit <- fit_score_matching_matrixG_fast(x, m, lambda_trace = 0, verbose = TRUE)

xgrid  <- seq(min(x) - 1, max(x) + 1, length.out = 500)
p_hat  <- density_from_fit(xgrid, fit$G, fit$c1)

p_norm <- dnorm(xgrid, 0, 1)
p_squared <- dnorm(xgrid, 0, 1 / sqrt(2))
p_logis <- dlogis(xgrid)
p_laplace <- dlaplace(xgrid)

p_exp <- dexp(xgrid, rate = 1)
p_chisq<- dchisq(xgrid, df = 5)
p_gamma <- dgamma(xgrid, shape = 2, rate = 1)
p_beta <- dbeta(xgrid, shape1 = 2, shape2 = 3)


hist(x, breaks = 30, freq = FALSE, col = "grey90", border = "grey70",
     main = "Score-Matching log-concave density vs. Normal", xlab = "x")
lines(xgrid, p_hat,  lwd = 2, col = "blue")
lines(xgrid, p_norm, lwd = 2, col = "red", lty = 2)
legend("topright",
       legend = c("Score matching estimate", "Normal N(0,1)"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n")

