# E-step for Gaussian mixture regression

e_step_fmr <- function(dat, em_state, family) {
  A <- dat$X_het
  B <- dat$X_com
  y <- dat$y
  n <- dat$n

  beta_g <- em_state$beta_g
  beta <- em_state$beta
  pi_g <- em_state$pi_g
  sigma_g <- em_state$sigma_g

  G <- nrow(beta_g)

  eta <- linear_predictor_matrix(
    A = A,
    B = B,
    beta_g = beta_g,
    beta = beta
  )

  eta <- as.matrix(eta)
  storage.mode(eta) <- "double"


  log_pi <- log(pmax(pi_g, 1e-16))
  log_w <- matrix(NA_real_, nrow = n, ncol = G)

  # Family specific E-step quantities

  if (family == "gaussian") {
    sigma_g <- as.numeric(sigma_g)

    if (length(sigma_g) != G) {
      stop("length(sigma_g) must equal nrow(beta_g).")
    }

    for (g in seq_len(G)) {
      log_w[, g] <- log_pi[g] + stats::dnorm(
        x = y,
        mean = eta[, g],
        sd = pmax(sigma_g[g], 1e-16),
        log = TRUE
      )
    }
  } else if (family == "poisson") {
    mu <- exp(eta)

    for (g in seq_len(G)) {
      log_w[, g] <- log_pi[g] + stats::dpois(
        x = y,
        lambda = mu[, g],
        log = TRUE
      )
    }
  } else if (family == "binomial") {
    binomial_size <- dat$binomial_size

    mu <- stats::plogis(eta)
    mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)

    for (g in seq_len(G)) {
      log_w[, g] <- log_pi[g] + stats::dbinom(
        x = y,
        size = binomial_size,
        prob = mu[, g],
        log = TRUE
      )
    }
  }

  log_denom <- row_logsumexp(log_w)

  loglik <- sum(log_denom)
  tau <- exp(log_w - log_denom)

  pi_new <- colMeans(tau)
  pi_new <- pmax(pi_new, 1e-16)
  pi_new <- pi_new / sum(pi_new)

  em_state$loglik <- loglik
  em_state$tau <- tau
  em_state$pi_g <- pi_new
  em_state$eta <- eta
}
