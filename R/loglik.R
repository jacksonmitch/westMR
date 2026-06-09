# Observed-data log-likelihood for Gaussian, Poisson, or binomial mixture regression

obs_loglik_gmr <- function(A, B, y, beta_g, beta, pi_g,
                           sigma_g = NULL,
                           family = c("gaussian", "poisson", "binomial"),
                           trials = NULL) {
  family <- match.arg(family)
  
  A <- as.matrix(A)
  
  if (is.null(B)) {
    B <- matrix(numeric(0), nrow = nrow(A), ncol = 0)
  } else {
    B <- as.matrix(B)
  }
  
  y <- as.numeric(y)
  beta_g <- as.matrix(beta_g)
  beta <- as.numeric(beta)
  pi_g <- as.numeric(pi_g)
  
  n <- nrow(A)
  G <- nrow(beta_g)
  
  eta <- linear_predictor_matrix(A, B, beta_g, beta)
  
  log_pi <- log(pmax(pi_g, 1e-16))
  log_w <- matrix(NA_real_, nrow = n, ncol = G)
  
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
  
  if (family == "poisson") {
    mu <- exp(eta)
    
    for (g in seq_len(G)) {
      log_w[, g] <- log_pi[g] + stats::dpois(
        x = y,
        lambda = mu[, g],
        log = TRUE
      )
    }
  }
  
  if (family == "binomial") {
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
  
  sum(row_logsumexp(log_w))
}
