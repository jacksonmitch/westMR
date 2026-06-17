# EM algorithm for finite mixture regression

em_fmr <- function(prepared_data,
                   G,
                   tau,
                   family = c("gaussian", "poisson", "binomial"),
                   control,
                   max_iter = control$max_iter,
                   tol = control$tol,
                   beta_g_start = NULL,
                   beta_start = NULL,
                   sigma_g_start = NULL,
                   pi_g_start = NULL) {

  family <- match.arg(family)

  het <- prepared_data$X_het
  common <- prepared_data$X_com
  response <- prepared_data$y

  loglik_trace <- numeric(max_iter)
  converged <- FALSE
  sigma_floor <- control$sigma_floor

  # Initial values carried through EM iterations
  beta_g <- beta_g_start
  beta <- beta_start
  sigma_g <- sigma_g_start

  if (is.null(pi_g_start)) {
    pi_g <- colMeans(tau)
  } else {
    pi_g <- as.numeric(pi_g_start)
  }

  irwls_iterations <- integer(max_iter)
  irwls_converged <- logical(max_iter)

  # Initial log-likelihood before any EM updates
  loglik_old <- -Inf

  # TODO: pass prepared_data into the steps once they're updated
  for (iter in seq_len(max_iter)) {

    # M-step: update parameters using the selected family
    m_fit <- m_step_gmr(
      A = het,
      B = common,
      y = response,
      tau = tau,
      family = family,
      beta_g = beta_g,
      beta = beta,
      sigma_g = sigma_g,
      sigma_floor = sigma_floor,
      irwls_maxit = control$irwls_max_iter,
      irwls_tol = control$irwls_tol,
      weight_floor = control$weight_floor
    )

    beta_g <- m_fit$beta_g
    beta <- m_fit$beta
    sigma_g <- m_fit$sigma_g
    pi_g <- m_fit$pi_g

    if (!is.null(m_fit$irwls_iterations)) {
      irwls_iterations[iter] <- m_fit$irwls_iterations
    }

    if (!is.null(m_fit$irwls_converged)) {
      irwls_converged[iter] <- isTRUE(m_fit$irwls_converged)
    }

    # Observed-data log-likelihood after M-step
    loglik_new <- obs_loglik_gmr(
      A = het,
      B = common,
      y = response,
      beta_g = beta_g,
      beta = beta,
      sigma_g = sigma_g,
      pi_g = pi_g,
      family = family
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
      pi_g = pi_g,
      family = family
    )

    if (is.finite(loglik_old) && tol > 0) {
      ll_diff <- abs(loglik_new - loglik_old)
      ll_scale <- 1 + abs(loglik_old)

      if (ll_diff < tol * ll_scale) {
        converged <- TRUE
        break
      }
    }

    loglik_old <- loglik_new
  }

  loglik_trace <- loglik_trace[seq_len(iter)]
  irwls_iterations <- irwls_iterations[seq_len(iter)]
  irwls_converged <- irwls_converged[seq_len(iter)]

  out <- list(
    beta_g = beta_g,
    beta = beta,
    sigma_g = sigma_g,
    pi_g = pi_g,
    tau = tau,

    loglik = loglik_trace[length(loglik_trace)],
    loglik_trace = loglik_trace,

    iterations = iter,
    converged = converged,

    irwls_iterations = if (family == "gaussian") NULL else irwls_iterations,
    irwls_converged = if (family == "gaussian") NULL else irwls_converged,

    family = family,
    G = G
  )

  out
}
