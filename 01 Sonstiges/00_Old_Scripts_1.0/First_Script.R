library(CVXR)
library(pracma)


# -----------------------------
# Build M(x) and N(x) matrices
# -----------------------------
make_MN <- function(x, m) {
  # M_ij(x) = x^(i+j+1) / (i+j+1)   for s'(x) inner product
  # N_ij(x) = x^(i+j)              for s''(x) inner product
  #
  # Here i,j in {0,...,m-1}. In R we use 1..m but compute powers with (i0+j0)
  i0 <- 0:(m-1)
  j0 <- 0:(m-1)
  I  <- matrix(rep(i0, times = m), nrow = m, byrow = FALSE)
  J  <- matrix(rep(j0, each  = m), nrow = m, byrow = FALSE)
  
  N <- x^(I + J)
  M <- x^(I + J + 1) / (I + J + 1)
  
  list(M = M, N = N)
}

# -----------------------------
# Score matching objective
# -----------------------------
fit_score_matching_matrixG <- function(x, m, verbose=TRUE) {
  n <- length(x)
  
  # Build time
  t_build <- system.time({
    # Decision variables
    G  <- Variable(m, m, PSD = TRUE)
    c1 <- Variable(1)
    
    obj_terms <- vector("list", n)
    
    for (k in seq_len(n)) {
      MN <- make_MN(x[k], m)
      Mk <- MN$M
      Nk <- MN$N
      
      # s'(x_k) = <G, M(x_k)> + c1
      s1k <- sum_entries(G * Mk) + c1
      
      # s''(x_k) = <G, N(x_k)>
      s2k <- sum_entries(G * Nk)
      
      obj_terms[[k]] <- 0.5 * square(s1k) - s2k
    }
    
    objective <- Minimize(Reduce(`+`, obj_terms) / n)
    
    constraints <- list(G[m, m] >= 1e-6)
    prob <- Problem(objective, constraints)
  })
  
  # solve time
  t_solve <- system.time(
    sol <- solve(prob, solver = "SCS"))  # or "MOSEK" if you have it
  
  if (verbose) {
    cat("\n--- Timing ---\n")
    cat(sprintf("Build: %.3f sec (user=%.3f, sys=%.3f)\n",
                t_build[["elapsed"]], t_build[["user.self"]], t_build[["sys.self"]]))
    cat(sprintf("Solve: %.3f sec (user=%.3f, sys=%.3f)\n",
                t_solve[["elapsed"]], t_solve[["user.self"]], t_solve[["sys.self"]]))
    cat("--------------\n\n")
  }
  
  list(
    G = sol$getValue(G),
    c1 = sol$getValue(c1),
    solution = sol,
    timing = list(build = t_build, solve = t_solve)
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
  svals <- sapply(xgrid, s_function, G = G, c1 = c1)
  unnorm <- exp(-svals)
  Z <- pracma::trapz(xgrid, unnorm)
  unnorm / Z
}


# -----------------------------
# Example run
# -----------------------------
set.seed(1)
x <- rnorm(200)     # your samples here
m <- 4            # basis size -> degree roughly 2m for s(x)

print(fit_score_matching_matrixG)
t <- system.time({
  fit <- fit_score_matching_matrixG(x, m)
})
print(t)
fit$c1
fit$G

# -----------------------------------
# Plot result
# -----------------------------------
xgrid <- seq(min(x) - 1, max(x) + 1, length.out = 500)

p_hat <- density_from_fit(xgrid, fit$G, fit$c1)
p_norm <- dnorm(xgrid, mean = 0, sd = 1)

hist(x,
     breaks = 30,
     freq = FALSE,
     col = "grey90",
     border = "grey70",
     main = "Score-Matching log-concave density vs. Normal",
     xlab = "x")

lines(xgrid, p_hat, lwd = 2, col = "blue")
lines(xgrid, p_norm, lwd = 2, col = "red", lty = 2)

legend("topright",
       legend = c("Score matching estimate", "Normal N(0,1)"),
       col = c("blue", "red"),
       lwd = 2,
       lty = c(1, 2),
       bty = "n")

