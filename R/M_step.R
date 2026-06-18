# Wrapper function for M-Steps accross families and method
m_step <- function(dat, em_state, family, control) {
  if (dat$p_com > 0) {
    em_state <- switch(family,
                       gaussian = m_step_sqr_gaussian(dat, em_state, control),
                       poisson  = m_step_sqr_poisson(dat, em_state, control),
                       binomial = m_step_sqr_binomial(dat, em_state, control),
                       stop(sprintf("No sQR M-step for family '%s'", family))
    )

  } else {
    em_state <- switch(family,
                       gaussian = m_step_gaussian(dat, em_state, control),
                       poisson  = m_step_poisson(dat, em_state, control),
                       binomial = m_step_binomial(dat, em_state, control),
                       stop(sprintf("No M-step for family '%s'", family))
    )
  }
}

m_step_gaussian <- function(dat, em_state, control) {
  tau <- em_state[["tau"]]
  G <- ncol(tau)
  X <- dat$X_het
  y <- dat$y
  p <- dat$p_het

  beta_g <- matrix(NA_real_, nrow = G, ncol = p)

  for (g in seq_len(G)) {
    w <- tau[, g]          # n-vector
    Xw <- X * w            # n x p, broadcasts column-wise
    XtWX <- t(Xw) %*% X   # p x p
    XtWy <- t(Xw) %*% y   # p x 1

    beta_g[g, ] <- solve(XtWX + 1e-8 * diag(p), XtWy)
  }


  # TODO: look at this line and think of alternatives
  mu <- linear_predictor_matrix(A = dat$X_het,
                                B = dat$X_com,
                                beta_g = beta_g,
                                beta = NULL)
  res2 <- (dat$y - mu)^2

  # Would it be faster to do these with matrix form instead of col sums?
  # TODO: test with scripts (maybe do this when in C++ land I imagine)

  den <- pmax(colSums(tau), 1e-8)
  sigma_g <- sqrt(colSums(tau * res2) / den)

  stopifnot(!is.null(control$sigma_floor))
  sigma_g <- pmax(sigma_g, control$sigma_floor)

  em_state[["beta_g"]] <- beta_g
  em_state[["mu"]] <- mu
  em_state[["sigma_g"]] <- sigma_g
  em_state
}

m_step_poisson <- function(dat, em_state, control) {
  stop("poisson not done yet")
}

m_step_binomial <- function(dat, em_state, control) {
  stop("binomial not done yet")
}

# Gaussian M-step using structured QR

m_step_sqr_gaussian <- function(dat, em_state, control) {
  tau <- em_state[["tau"]]
  sigma_g <- em_state[["sigma_g"]]
  # Updated beta, beta_g and QR parts
  if (!is.null(sigma_g)) {
    sigma_g <- as.numeric(sigma_g)

    stopifnot(length(sigma_g) == ncol(tau))

    sigma_g <- pmax(sigma_g, 1e-16)

    w_beta <- sweep(
      tau,
      MARGIN = 2,
      STATS = sigma_g^2,
      FUN = "/"
    )
  } else {
    w_beta <- tau
  }

  out <- wls_sqr(dat = dat,
                 z = dat$y,
                 weights = w_beta,
                 control = control)
  beta <- out$beta
  beta_g <- out$beta_g

  # TODO: look at this line and think of alternatives
  mu <- linear_predictor_matrix(A = dat$X_het,
                                B = dat$X_com,
                                beta_g = beta_g,
                                beta = beta)

  res2 <- (dat$y - mu)^2
  # Would it be faster to do these with matrix form instead of col sums?
  # TODO: test with scripts (maybe do this when in C++ land I imagine)

  den <- pmax(colSums(tau), 1e-8)
  sigma_g <- sqrt(colSums(tau * res2) / den)

  stopifnot(!is.null(control$sigma_floor))
  sigma_g <- pmax(sigma_g, control$sigma_floor)

  em_state[["beta_g"]] <- beta_g
  em_state[["beta"]] <- beta
  em_state[["mu"]] <- mu
  em_state[["sigma_g"]] <- sigma_g
  em_state
}

# Poisson M-step using inner IRWLS and structured QR

m_step_sqr_poisson <- function(dat, em_state, control) {
  em_state <- irwls_fmr(dat, em_state, control, "poisson")
  em_state
}

# Binomial M-step using inner IRWLS and structured QR

m_step_sqr_binomial <- function(dat, em_state, control) {
  em_state <- irwls_fmr(dat, em_state, control, "binomial")
  em_state
}



