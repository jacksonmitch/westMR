irwls_gmr <- function(A, B, y, tau,
                      beta_g = NULL, beta = NULL, trials = NULL,
                      maxit = 50, tol = 1e-8, weight_floor = 1e-10,
                      family = c("poisson", "binomial")
                      ){

  family <- match.arg(family)

  A <- as.matrix(A)

  if (is.null(B)) {
    B <- matrix(numeric(0), nrow = nrow(A), ncol = 0)
  } else {
    B <- as.matrix(B)
  }

  y <- as.numeric(y)
  tau <- as.matrix(tau)

  n <- nrow(A)
  G <- ncol(tau)

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

    init <- wls_sqr(
      A = A,
      B = B,
      z = z0,
      w = tau
    )

    beta_g <- init$beta_g
    beta <- init$beta

  } else {
    beta_g <- as.matrix(beta_g)
    beta <- as.numeric(beta)
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

    fit <- wls_sqr(
      A = A,
      B = B,
      z = z,
      w = w
    )

    beta_g <- fit$beta_g
    beta <- fit$beta

    diff <- max(abs(c(beta_g - beta_g_old, beta - beta_old)))

    if (is.finite(diff) && diff < tol) {
      converged <- TRUE
      break
    }
  }

  list(
    beta_g = beta_g,
    beta = beta,
    iterations = r,
    converged = converged
  )
}
