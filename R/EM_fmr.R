# EM algorithm for finite mixture regression

em_fmr <- function(prepared_data,
                   G,
                   tau,
                   family = c("gaussian", "poisson", "binomial"),
                   control) {
  
  family <- match.arg(family)
  
  het <- prepared_data$X_het
  common <- prepared_data$X_com
  response <- prepared_data$y
  
  loglik_trace <- numeric(control$max_iter)
  converged <- FALSE
  sigma_floor <- control$sigma_floor
  
  # Initial values carried through EM iterations
  beta_g <- NULL
  beta <- NULL
  sigma_g <- NULL
  pi_g <- colMeans(tau)
  
  irwls_iterations <- integer(control$max_iter)
  irwls_converged <- logical(control$max_iter)
  
  # Initial log-likelihood before any EM updates
  loglik_old <- -Inf

  # TODO: pass prepared_data into the steps once they're updated
  for (iter in seq_len(control$max_iter)) {
    
    # M-step: update parameters using the selected family
    m_fit <- m_step_gmr(
      A = het,
      B = common,
      y = response,
      tau = tau,
      family = family,
      beta_g = beta_g,
      beta = beta,
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
    
    # Convergence check
    if (is.finite(loglik_old) && abs(loglik_new - loglik_old) < control$tol) {
      converged <- TRUE
      break
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
