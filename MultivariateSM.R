build_pairwise_score_matching <- function(x, m, build_Phi = TRUE,
                                          drop_constant = TRUE) {
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("x must be numeric.")
  n <- nrow(x)
  d <- ncol(x)
  if (m < 1) stop("m must be >= 1.")
  
  # ----------------------------------------------------------
  # Basis:
  #   univariate: x_l^r, r = 0, ..., 2m-2
  #   pairwise:   x_l^i x_u^j, 1 <= l < u <= d, i,j = 0, ..., m-1
  #
  # Optional: Konstante weglassen, da sie fuer Score Matching
  # identifizierungsmaessig nutzlos ist (Gradient/Laplacian = 0).
  # ----------------------------------------------------------
  
  uni_deg <- 0:(2 * m - 2)
  pair_deg <- 0:(m - 1)
  
  basis <- vector("list", 0)
  
  # univariate basis
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
  
  # pairwise basis
  for (l in seq_len(d - 1L)) {
    for (u in (l + 1L):d) {
      for (i in pair_deg) {
        for (j in pair_deg) {
          # konstante 1 und reine univariate Doppelungen vermeiden
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
  # Potenzen vorkalkulieren:
  # Xpow[[l]][, r+1] = x[,l]^r, fuer r = 0, ..., 2m-2
  # ----------------------------------------------------------
  max_deg <- 2 * m - 2
  Xpow <- vector("list", d)
  for (l in seq_len(d)) {
    P <- matrix(1, nrow = n, ncol = max_deg + 1L)
    if (max_deg >= 1) {
      P[, 2] <- x[, l]
      if (max_deg >= 2) {
        for (r in 2:max_deg) {
          P[, r + 1L] <- P[, r] * x[, l]
        }
      }
    }
    Xpow[[l]] <- P
  }
  
  # Hilfsfunktion: x^k mit Konvention 0 falls k < 0
  get_pow <- function(l, deg) {
    if (deg < 0) return(rep(0, n))
    Xpow[[l]][, deg + 1L]
  }
  
  # ----------------------------------------------------------
  # Optional Phi sowie Ableitungsmatrizen D_k und Laplacian-Matrix
  # D[[k]][n,j] = d/dx_k phi_j(x_n)
  # Lap[n,j]    = Delta phi_j(x_n)
  # ----------------------------------------------------------
  Phi <- if (build_Phi) matrix(0, nrow = n, ncol = p) else NULL
  D <- lapply(seq_len(d), function(k) matrix(0, nrow = n, ncol = p))
  Lap <- matrix(0, nrow = n, ncol = p)
  
  for (col in seq_len(p)) {
    b <- basis[[col]]
    
    if (b$type == "uni") {
      l <- b$l
      r <- b$r
      
      # Phi
      if (build_Phi) {
        Phi[, col] <- get_pow(l, r)
      }
      
      # gradient
      if (r >= 1) {
        D[[l]][, col] <- r * get_pow(l, r - 1L)
      }
      
      # laplacian
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
      
      # Phi
      if (build_Phi) {
        Phi[, col] <- xl_i * xu_j
      }
      
      # gradient wrt x_l
      if (i >= 1) {
        D[[l]][, col] <- i * get_pow(l, i - 1L) * xu_j
      }
      
      # gradient wrt x_u
      if (j >= 1) {
        D[[u]][, col] <- j * xl_i * get_pow(u, j - 1L)
      }
      
      # laplacian = second wrt l + second wrt u
      term_l <- if (i >= 2) i * (i - 1) * get_pow(l, i - 2L) * xu_j else 0
      term_u <- if (j >= 2) j * (j - 1) * xl_i * get_pow(u, j - 2L) else 0
      Lap[, col] <- term_l + term_u
    }
  }
  
  # ----------------------------------------------------------
  # K = (1/n) sum_k crossprod(D_k)
  # ell = colMeans(Lap)
  # ----------------------------------------------------------
  K <- matrix(0, nrow = p, ncol = p)
  for (k in seq_len(d)) {
    K <- K + crossprod(D[[k]])
  }
  K <- K / n
  ell <- colMeans(Lap)
  
  # sprechende Namen
  basis_names <- vapply(basis, function(b) {
    if (b$type == "uni") {
      sprintf("x%d^%d", b$l, b$r)
    } else {
      sprintf("x%d^%d*x%d^%d", b$l, b$i, b$u, b$j)
    }
  }, character(1))
  
  colnames(K) <- rownames(K) <- basis_names
  names(ell) <- basis_names
  if (!is.null(Phi)) colnames(Phi) <- basis_names
  
  list(
    Phi = Phi,
    K = K,
    ell = ell,
    basis = basis,
    basis_names = basis_names,
    D = D
  )
}