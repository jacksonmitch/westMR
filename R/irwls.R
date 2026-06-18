
irwls_fmr <- function(dat, em_state, control,
                      family = c("poisson", "binomial")) {
  family <- match.arg(family)
  A <- dat$X_het
  B <- dat$X_com
  stopifnot(!is.null(dat$B))
  y <- dat$y

  tau <- em_state[["tau"]]
  beta_g <- em_state[["beta_g"]]
  beta <- em_state[["beta"]]

  n <- nrow(A)
  G <- ncol(tau)

  maxit <- control$maxit
  tol <- control$tol
  weight_floor <- control$weight_floor

  # I think trials should be saved somewhere but don't know where
  trials <- dat$trials
  # this will keep it being generated until we decide where it goes, but
  # calculations stay safe
  if (is.null(trials)) {
    trials <- rep(1, n)
  }

  trials <- as.numeric(trials)

  # Initial values for beta_g and beta
  if (is.null(beta_g) || is.null(beta)) {
    if (family == "poisson") {
      z0 <- log(pmax(y, 0) + 0.25)
    }

    if (family == "binomial") {
      p0 <- (y + 0.5) / (trials + 1)
      p0 <- pmin(pmax(p0, 1e-8), 1 - 1e-8)
      z0 <- stats::qlogis(p0)
    }
    out <- wls_sqr(
      dat = dat,
      z = z0,
      weights = tau,
      control = control
    )
    beta_g <- out[["beta_g"]]
    beta <- out[["beta"]]
  }

  converged <- FALSE

  for (r in seq_len(maxit)) {

    beta_g_old <- beta_g
    beta_old <- beta

    eta <- linear_predictor_matrix(A, B, beta_g, beta)

    if (family == "poisson") {

      mu <- exp(eta)

      z <- eta + (y - mu) / pmax(mu, 1e-10)
      w <- tau * mu
    }

    if (family == "binomial") {
      eta_safe <- pmin(pmax(eta, -30), 30)
      mu <- stats::plogis(eta_safe)
      mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)

      y_bar <- y / trials
      y_bar_mat <- matrix(y_bar, nrow = n, ncol = G)
      trials_mat <- matrix(trials, nrow = n, ncol = G)

      z <- eta_safe + (y_bar_mat - mu) / pmax(mu * (1 - mu), 1e-10)
      w <- tau * trials_mat * mu * (1 - mu)

    }

    w <- pmax(w, weight_floor)


    out <- wls_sqr(
      dat = dat,
      z = z,
      weights = w,
      control = control
    )
    beta_g <- out[["beta_g"]]
    beta <- out[["beta"]]

    diff <- max(abs(c(beta_g - beta_g_old, beta - beta_old)))

    if (is.finite(diff) && diff < tol) {
      converged <- TRUE
      break
    }
  }

  em_state[["beta_g"]] <- beta_g
  em_state[["beta"]] <- beta
  em_state[["irwls_iterations"]] <- r
  em_state[["irwls_converged"]] <- converged
  em_state
}
