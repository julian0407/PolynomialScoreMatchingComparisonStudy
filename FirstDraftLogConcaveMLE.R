# ============================================================
# LOG-CONCAVE MLE (1D) für die Vergleichsstudie
# basiert praktisch auf dem Active-Set-Ansatz via logcondens
# ============================================================

# Paket laden / installieren
load_logcondens <- function() {
  if (!requireNamespace("logcondens", quietly = TRUE)) {
    install.packages("logcondens")
  }
  invisible(TRUE)
}

# ------------------------------------------------------------
# (1) Fit: univariater log-concave MLE
#     optional auch smoothed version
# ------------------------------------------------------------
fit_logconcave_mle_1d <- function(x, smoothed = FALSE, print = FALSE) {
  stopifnot(is.numeric(x), length(x) >= 2)
  x <- x[is.finite(x)]
  
  if (length(unique(x)) < 2) {
    stop("Für den log-concave MLE braucht man mindestens zwei verschiedene Beobachtungen.")
  }
  
  load_logcondens()
  
  # logConDens ist laut CRAN die Hauptfunktion für MLE + smoothed MLE
  # und berechnet den univariaten log-concave Schätzer. :contentReference[oaicite:5]{index=5}
  fit <- logcondens::logConDens(
    x = sort(x),
    smoothed = smoothed,
    print = print
  )
  
  structure(
    list(
      fit = fit,
      smoothed = smoothed,
      support = range(fit$x),
      class_name = "logconcave_mle_1d"
    ),
    class = "logconcave_mle_1d"
  )
}

# ------------------------------------------------------------
# (2) Dichteauswertung
#     MLE: log f ist stückweise linear auf [x_(1), x_(n)],
#          außerhalb = 0
# ------------------------------------------------------------
density_from_logconcave_fit <- function(xgrid, obj, which = c("MLE", "smoothed")) {
  which <- match.arg(which)
  stopifnot(inherits(obj, "logconcave_mle_1d"))
  
  fit <- obj$fit
  xgrid <- as.numeric(xgrid)
  
  if (which == "MLE") {
    xmin <- min(fit$x)
    xmax <- max(fit$x)
    
    inside <- (xgrid >= xmin) & (xgrid <= xmax)
    out <- numeric(length(xgrid))
    
    if (any(inside)) {
      # phi = log-density an den support points
      phi_eval <- approx(
        x = fit$x,
        y = fit$phi,
        xout = xgrid[inside],
        method = "linear",
        rule = 2
      )$y
      
      out[inside] <- exp(phi_eval)
    }
    
    return(out)
  }
  
  if (!isTRUE(obj$smoothed)) {
    stop("Dieser Fit enthält keine smoothed-Version. Bitte mit smoothed = TRUE fitten.")
  }
  
  # smoothed estimator wird vom Paket auf einem Grid xs/f.smoothed geliefert
  # Wir interpolieren numerisch auf xgrid.
  xs <- fit$xs
  fs <- fit$f.smoothed
  
  if (is.null(xs) || is.null(fs)) {
    stop("Smoothed-Ausgabe im Fit nicht gefunden.")
  }
  
  approx(
    x = xs,
    y = fs,
    xout = xgrid,
    method = "linear",
    rule = 2
  )$y
}

# ------------------------------------------------------------
# (3) Out-of-sample Log-Likelihood
# ------------------------------------------------------------
oos_loglik_logconcave <- function(x_test, obj, which = c("MLE", "smoothed"), eps = 1e-300) {
  which <- match.arg(which)
  dens <- density_from_logconcave_fit(x_test, obj, which = which)
  mean(log(pmax(dens, eps)))
}

# ------------------------------------------------------------
# (4) Fehlermaße
# ------------------------------------------------------------
hellinger_error <- function(xg, p, q) {
  # H^2 = 1/2 \int (sqrt(p)-sqrt(q))^2 dx
  sqrt(0.5 * pracma::trapz(xg, (sqrt(pmax(p, 0)) - sqrt(pmax(q, 0)))^2))
}

ISE_error <- function(xg, p, q) {
  pracma::trapz(xg, (p - q)^2)
}

L1_error_num <- function(xg, p, q) {
  pracma::trapz(xg, abs(p - q))
}

# ------------------------------------------------------------
# (5) Ein einzelner Benchmark-Lauf für log-concave MLE
# ------------------------------------------------------------
benchmark_logconcave_mle <- function(
    x_train,
    x_test,
    xg_eval,
    p_true_eval,
    smoothed = FALSE
) {
  t_fit <- system.time({
    fit <- fit_logconcave_mle_1d(x_train, smoothed = smoothed, print = FALSE)
  })[["elapsed"]]
  
  t_dens <- system.time({
    p_hat <- density_from_logconcave_fit(
      xgrid = xg_eval,
      obj = fit,
      which = if (smoothed) "smoothed" else "MLE"
    )
  })[["elapsed"]]
  
  list(
    fit = fit,
    p_hat = p_hat,
    L1 = L1_error_num(xg_eval, p_hat, p_true_eval),
    ISE = ISE_error(xg_eval, p_hat, p_true_eval),
    Hellinger = hellinger_error(xg_eval, p_hat, p_true_eval),
    loglik_oos = oos_loglik_logconcave(
      x_test,
      fit,
      which = if (smoothed) "smoothed" else "MLE"
    ),
    t_fit = t_fit,
    t_dens = t_dens
  )
}


# Beispiel: Logistic
set.seed(1)

n <- 20000
x_train <- rlogis(n, location = 0, scale = 1)
x_test  <- rlogis(10000, location = 0, scale = 1)

xg <- seq(-8, 8, length.out = 2001)
p_true <- dlogis(xg, location = 0, scale = 1)

res_lc <- benchmark_logconcave_mle(
  x_train = x_train,
  x_test = x_test,
  xg_eval = xg,
  p_true_eval = p_true,
  smoothed = FALSE
)

res_lc_s <- benchmark_logconcave_mle(
  x_train = x_train,
  x_test = x_test,
  xg_eval = xg,
  p_true_eval = p_true,
  smoothed = TRUE
)

res_lc$L1
res_lc$ISE
res_lc$Hellinger
res_lc$loglik_oos

plot(xg, p_true, type = "l", lwd = 2, col = 2,
     main = "Logistic: truth vs log-concave estimators",
     ylab = "density", xlab = "x")
lines(xg, res_lc$p_hat, lwd = 2, col = 4)
lines(xg, res_lc_s$p_hat, lwd = 2, col = 3, lty = 2)
legend("topright",
       legend = c("True", "Log-concave MLE", "Smoothed log-concave MLE"),
       col = c(2, 4, 3), lty = c(1, 1, 2), lwd = 2, bty = "n")