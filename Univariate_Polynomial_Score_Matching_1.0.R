library(CVXR)
library(pracma)

# ------------------------------------------------------------
# (1) Precompute matrix K and vector l that are used for the
#     estimated score loss J_hat = (1/2)*y^T*K*y - l^T*y
# ------------------------------------------------------------

# Precompute matrices that are later used for entry-wise 
# matrix operations instead of loop-wise executions
precompute_exponents_univariate <- function(m) {
  # Precompute a mxm matrix containing i+j at the (i,j)-th entry
  # Needed for second devirative
  expoN <- as.vector(outer(0:(m - 1), 0:(m - 1), `+`))
  # Add one to each entry for matrix M
  # Needed for first devirative
  expoM <- expoN + 1
  # (Redundant) Save same values for the denominators in M
  denomM <- expoM
  # Return list of matrices for exponents used in M,N and
  # denominators in M
  list(
    expoN = expoN,
    expoM = expoM,
    denomM = denomM
  )
}

# Use precomputed exponents and denominators to build matrix
# K and vector l and a weighting function h and its derivative
build_K_and_l_univariate <- function(
    z,
    m,
    h = function(z) rep(1, length(z)),
    h_prime = function(z) rep(0, length(z))
    ) {
      n <- length(z)
      z <- as.numeric(z)
      # Build exponents
      exps <- precompute_exponents_univariate(m)
      # For every entry in z apply exponents in expoM (treated as vector)
      # Result is a matrix containing z_i^expoM_j as (i,j)-th element
      # --> Matrix with n rows and mxm columns
      A <- outer(z, exps$expoM, `^`)
      # Do the same with expoN
      B <- outer(z, exps$expoN, `^`)
      # For the first derivative we need to apply on each column j the
      # denominator the debominator denomM[j]
      A <- sweep(A, 2, exps$denomM, "/")
      
      # Get weights by applying function h and h_prime
      hx  <- as.numeric(h(z))
      hpx <- as.numeric(h_prime(z))
      # Only allow positive values for h
      # TODO: Check if this make sense
      hx <- pmax(hx, 0)
      
      # Use column scaling for numerical stability (std for each column)
      # as big m can be very unstable for exponentials
      # TODO: Mention in Report
      scale_vec <- apply(A, 2, function(v) max(stats::sd(v), 1e-12))
      # Apply coulmn scaling
      A_sc <- sweep(A, 2, scale_vec, "/")
      B_sc <- sweep(B, 2, scale_vec, "/")
      
      # Calculate aggregated moments
      Aw <- A_sc * sqrt(hx)
      S  <- crossprod(Aw) / n
      t  <- as.vector(crossprod(A_sc, hx) / n)
      u  <- sum(hx) / n
      r  <- as.vector(crossprod(A_sc, hpx) / n)
      q  <- sum(hpx) / n
      b_bar <- as.vector(crossprod(B_sc, hx) / n)
      
      # Use the moments to derive matrix K and vector l
      K <- rbind(
        cbind(S, matrix(t, ncol = 1)),
        cbind(matrix(t, nrow = 1), matrix(u, 1, 1))
      )
      K <- 0.5 * (K + t(K))
      l <- c(r + b_bar, q)
      
      # Return K, l und the column scaling vector
      list(K = K, l = l, scale_vec = scale_vec)
}


# ------------------------------------------------------------
# (2) Core fit: univariate Score Matching
# ------------------------------------------------------------
# - Use identity function as standard for weighting function h
fit_score_matching_univariate <- function(
    x,
    m,
    h = function(z) rep(1, length(z)),
    h_prime = function(z) rep(0, length(z)),
    standardize = TRUE,
    solver = "SCS",
    scs_control = list(max_iters = 100000, eps = 1e-5, alpha = 1.8, verbose = FALSE)
) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2) stop("Data vector x must contain at least 
                          two finite observations.")
  if (m < 1) stop("Degree m must be >= 1.")
  
  # optional standardization
  if (standardize) {
    sc <- scale_vector_1d(x)
    z <- sc$x_scaled
  } else {
    z <- x
    sc <- list(center = 0, scale = 1)
  }
  
  # Precompute matrix K and vector l which determine the estimated score loss
  score_loss_input <- build_K_and_l_univariate(z, m, h, h_prime)
  
  # Specify optimization variables G, c1 and set PSD constraint to True 
  G  <- Variable(m, m, PSD = TRUE)
  c1 <- Variable(1)
  # Transform these variables to vector representation
  p = m*m
  gvec <- tryCatch(vec(G), error = function(e) reshape(G, c(p, 1)))
  y <- vstack(gvec, c1)
  
  # Define Score loss as objective function by using the precomputed
  # matrix K and vector l
  obj <- 0.5 * quad_form(y, Constant(score_loss_input$K)) -
    t(Constant(matrix(score_loss_input$l, ncol = 1))) %*% y
  
  # Specify CVXR problem
  prob <- Problem(Minimize(obj))
  
  # Solve CVXR problem using provided solver settings
  sol <- solve(
    prob,
    solver = solver,
    max_iters = scs_control$max_iters,
    eps = scs_control$eps,
    alpha = scs_control$alpha,
    verbose = scs_control$verbose
  )
  
  # Back-transform column scaled solution for g_vec
  g_sc  <- as.numeric(sol$getValue(gvec))
  g_org <- g_sc / score_loss_input$scale_vec
  # Back-transform vector representation of g to Matrix G
  G_org <- matrix(g_org, nrow = m, ncol = m)
  
  # TODO: why using structure and class?
  structure(
    list(
      G = G_org,
      c1 = as.numeric(sol$getValue(c1)),
      status = sol$status,
      solution = sol,
      scaling = list(
        center = sc$center,
        scale = sc$scale,
        standardize = standardize
      ),
      column_scaling = list(
        scale_vec = score_loss_input$scale_vec
      ),
      z_train = z,
      m = m
    ),
    class = "score_matching_univariate_fit"
  )
}


# ------------------------------------------------------------
# (3) Evaluation metrics
# ------------------------------------------------------------

# ------------------------------------------------------------
# (3.1) Evaluate Score on original scale of x
# ------------------------------------------------------------
predict_score_z_univariate <- function(z, fit) {
  # TODO: Desccribe procedure
  z <- as.numeric(z)
  m <- fit$m
  exps <- precompute_exponents_univariate(m)
  A <- outer(z, exps$expoM, `^`)
  A <- sweep(A, 2, exps$denomM, "/")
  gvec <- as.vector(fit$G)
  as.numeric(A %*% gvec) + fit$c1
}

predict_score_univariate <- function(x, fit) {
  z <- apply_scaling_1d(x, fit$scaling)$z
  score_z <- predict_score_z_univariate(z, fit)
  # TODO: Check this
  # r_z(z) = s * r_x(x)  =>  r_x(x) = r_z(z) / s
  score_z / fit$scaling$scale
}

