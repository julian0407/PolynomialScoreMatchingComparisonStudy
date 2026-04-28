library(CVXR)
library(pracma)

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
fit_score_matching_matrixG <- function(x, m, lambda_trace = 0, verbose = TRUE) {
  n <- length(x)
  
  # -------------------------
  # Build phase
  # -------------------------
  t_build <- system.time({
    G  <- Variable(m, m, PSD = TRUE)
    c1 <- Variable(1)
    
    obj_terms <- vector("list", n)
    
    for (k in seq_len(n)) {
      MN <- make_MN(x[k], m)
      Mk <- MN$M
      Nk <- MN$N
      
      s1k <- sum_entries(G * Mk) + c1
      s2k <- sum_entries(G * Nk)
      
      obj_terms[[k]] <- 0.5 * square(s1k) - s2k
    }
    
    base_obj <- Reduce(`+`, obj_terms)
    
    # Trace-Regularisierung versionssicher:
    trG <- sum(diag(G))
    
    objective <- Minimize(base_obj / n + lambda_trace * trG)
    constraints <- list(G[m, m] >= 1e-6)
    
    prob <- Problem(objective, constraints)
  })
  
  # -------------------------
  # Solve phase
  # -------------------------
  t_solve <- system.time({
    sol <- solve(prob, solver = "SCS")
  })
  
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
x <- rnorm(200)
m <- 4

t <- system.time({
  fit <- fit_score_matching_matrixG(x, m)
})
print(t)

# CVXR gibt c1 oft als 1x1 zurück -> numerisch machen
c1_hat <- as.numeric(fit$c1)
G_hat  <- fit$G

cat("c1 =", c1_hat, "\n")
print(G_hat)

# -----------------------------------
# Plot result
# -----------------------------------
xgrid <- seq(min(x) - 1, max(x) + 1, length.out = 500)

# Option A: nutze deine density_from_fit direkt
p_hat <- density_from_fit(xgrid, G_hat, c1_hat)

# Option B (robuster gegen Overflow): stabilisierte Dichte aus s_function
# (Wenn Option A bei dir NaNs/Infs produziert, nimm diese Version stattdessen)
# svals  <- sapply(xgrid, s_function, G = G_hat, c1 = c1_hat)
# svals  <- svals - min(svals)          # shift für Stabilität
# unnorm <- exp(-svals)
# Z <- pracma::trapz(xgrid, unnorm)
# p_hat <- unnorm / Z

# Vergleich: Normal
p_norm <- dnorm(xgrid, mean = 0, sd = 1)

hist(x,
     breaks = 30,
     freq   = FALSE,
     col    = "grey90",
     border = "grey70",
     main   = "Score-Matching log-concave density vs. Normal",
     xlab   = "x")

lines(xgrid, p_hat,  lwd = 2, col = "blue")
lines(xgrid, p_norm, lwd = 2, col = "red", lty = 2)

legend("topright",
       legend = c("Score matching estimate", "Normal N(0,1)"),
       col    = c("blue", "red"),
       lwd    = 2,
       lty    = c(1, 2),
       bty    = "n")



