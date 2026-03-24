library(CVXR)

# vech of a symmetric matrix (upper triangle, including diagonal)
vech_upper <- function(A) {
  A[upper.tri(A, diag = TRUE)]
}

# Weights for converting vech-dot to Frobenius inner product:
# <G, M> = sum_i G_ii M_ii + 2 * sum_{i<j} G_ij M_ij
vech_weights <- function(m) {
  W <- matrix(2, m, m)
  diag(W) <- 1
  W[lower.tri(W)] <- NA
  vech_upper(W)
}

# Build a symmetric CVXR matrix expression from g = vech(G)
# Order matches vech_upper(): column-wise upper triangle extraction in R
vech_to_sym_expr <- function(g, m) {
  Gexpr <- matrix(list(), m, m)
  idx <- 1
  for (j in 1:m) {
    for (i in 1:j) {
      Gexpr[[i, j]] <- g[idx]
      Gexpr[[j, i]] <- g[idx]
      idx <- idx + 1
    }
  }
  # Convert list-matrix to CVXR Expression matrix via bmat
  rows <- lapply(1:m, function(i) lapply(1:m, function(j) Gexpr[[i, j]]))
  bmat(rows)
}



# Zweiter Part

make_MN <- function(x, m) {
  i0 <- 0:(m-1)
  j0 <- 0:(m-1)
  I  <- matrix(rep(i0, times = m), nrow = m, byrow = FALSE)
  J  <- matrix(rep(j0, each  = m), nrow = m, byrow = FALSE)
  
  N <- x^(I + J)
  M <- x^(I + J + 1) / (I + J + 1)
  
  list(M = M, N = N)
}

fit_score_matching_vectorG <- function(x, m) {
  n <- length(x)
  d <- m * (m + 1) / 2
  
  g  <- Variable(d)   # vech(G)
  c1 <- Variable(1)
  
  # weights so that g' * (w * vech(M)) equals <G, M>
  w <- vech_weights(m)
  
  obj_terms <- vector("list", n)
  
  for (k in seq_len(n)) {
    MN <- make_MN(x[k], m)
    mk <- vech_upper(MN$M)
    nk <- vech_upper(MN$N)
    
    # <G, M> using vech + weights:
    # inner = sum_{i<=j} g_ij * (w_ij * m_ij)
    inner_M <- sum(g * (w * mk))
    inner_N <- sum(g * (w * nk))
    
    s1k <- inner_M + c1
    s2k <- inner_N
    
    obj_terms[[k]] <- 0.5 * square(s1k) - s2k
  }
  
  objective <- Minimize(sum(obj_terms) / n)
  
  # PSD constraint: build G(g) and constrain it
  Gexpr <- vech_to_sym_expr(g, m)
  constraints <- list(Gexpr >> 0)
  
  prob <- Problem(objective, constraints)
  sol <- solve(prob, solver = "SCS")  # or MOSEK
  
  list(g = sol$getValue(g),
       c1 = sol$getValue(c1),
       G = sol$getValue(Gexpr),
       solution = sol)
}

# Example run
set.seed(1)
x <- rnorm(200)
m <- 4

fit2 <- fit_score_matching_vectorG(x, m)
fit2$c1
fit2$G

