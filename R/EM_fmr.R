# EM algorithm for Gaussian mixture regression

em_fmr <- function(prepared_data,
                   G,
                   tau,
                   family,
                   control) {

  # Normalize mixing proportions just in case
  # pi_g <- pi_g / sum(pi_g)
  het <- prepared_data$X_het
  common <- prepared_data$X_com
  response <- prepared_data$y

  loglik_trace <- numeric(control$max_iter)
  converged <- FALSE
  sigma_floor <- control$sigma_floor

  # Initial log-likelihood before any EM updates
  loglik_old <- 0

  # TODO: pass prepared_data into the steps once they're updated
  for (iter in seq_len(control$max_iter)) {
    # M-step: update parameters using selected method
    m_fit <- m_step_sqr(
      A = het,
      B = common,
      y = response,
      tau = tau,
      sigma_floor = sigma_floor
    )

    beta_g <- m_fit$beta_g
    beta <- m_fit$beta
    sigma_g <- m_fit$sigma
    pi_g <- m_fit$pi

    # Observed-data log-likelihood after M-step
    loglik_new <- obs_loglik_gmr(
      A = het,
      B = common,
      y = response,
      beta_g = beta_g,
      beta = beta,
      sigma_g = sigma_g,
      pi_g = pi_g
    )

    loglik_trace[iter] <- loglik_new

    # E-step: compute posterior probabilities
    tau <- e_step_gmr(
      A = het,
      B = common,
      y = response,
      beta_g = beta_g,
      beta = beta,
      sigma_g = sigma_g,
      pi_g = pi_g
    )

    # Convergence check
    if (abs(loglik_new - loglik_old) < control$tol) {
      converged <- TRUE
      break
    }

    loglik_old <- loglik_new
  }

  loglik_trace <- loglik_trace[seq_len(iter)]

  list(
    beta_g = beta_g,
    beta = beta,
    sigma_g = sigma_g,
    pi_g = pi_g,
    tau = tau,
    loglik = loglik_trace[length(loglik_trace)],
    loglik_trace = loglik_trace,
    iterations = iter,
    converged = converged,
    G = G
  )
}
