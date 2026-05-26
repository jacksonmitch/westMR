# Observed-data log-likelihood for Gaussian mixture regression

obs_loglik_gmr <- function(A, B, y, beta_g, beta, sigma_g, pi_g) {
  A <- as.matrix(A)
  if (is.null(B)) {
    B <- matrix(
      numeric(0),
      nrow = nrow(A),
      ncol = 0
    )
  }
  else{
    B <- as.matrix(B)
  }
  y <- as.numeric(y)
  beta_g <- as.matrix(beta_g)
  beta <- as.numeric(beta)
  sigma_g <- as.numeric(sigma_g)
  pi_g <- as.numeric(pi_g)

  mu <- sweep(A %*% t(beta_g), 1, common_eta(B, beta),"+"
    )

  res2 <- (y - mu)^2

  log_w <- sweep(-0.5 * res2, 2, sigma_g^2, "/")

  log_w <- sweep(
    log_w,
    2,
    log(pi_g + 1e-16) - log(sigma_g + 1e-16),
    "+"
  ) - 0.5 * log(2 * pi)

  sum(row_logsumexp(log_w))
}
