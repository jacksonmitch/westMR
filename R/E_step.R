# E-step for Gaussian mixture regression

e_step_fmr <- function(dat, em_state, family) {

  A <- dat$X_het
  B <- dat$X_com
  y <- dat$y
  n <- dat$n

  beta_g <- em_state[["beta_g"]]
  beta <- em_state[["beta"]]
  pi_g <- em_state[["pi_g"]]
  sigma_g <- em_state[["sigma_g"]]
  eta <- em_state[["mu"]]

  G <- nrow(beta_g)


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
  }
  else if (family == "poisson") {
    mu <- exp(eta)

    for (g in seq_len(G)) {
      log_w[, g] <- log_pi[g] + stats::dpois(
        x = y,
        lambda = mu[, g],
        log = TRUE
      )
    }
  }
  else if (family == "binomial") {
    trials <- dat$trials
    if (is.null(trials)) {
      trials <- rep(1, n)
    }

    trials <- as.numeric(trials)

    if (length(trials) != n) {
      stop("length(trials) must equal length(y).")
    }

    mu <- stats::plogis(eta)
    mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)

    for (g in seq_len(G)) {
      log_w[, g] <- log_pi[g] + stats::dbinom(
        x = y,
        size = trials,
        prob = mu[, g],
        log = TRUE
      )
    }
  }

  row_max <- apply(log_w, 1, max) # big time sink
  w <- exp(log_w - row_max)
  w_sum <- rowSums(w)

  loglik <- sum(row_max + log(w_sum))
  tau <- w / w_sum
  colnames(tau) <- paste0("g", seq_len(G)) #is this necesary?

  em_state[["loglik"]] <- loglik
  em_state[["tau"]] <- tau
  em_state
}
