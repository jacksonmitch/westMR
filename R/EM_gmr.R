# EM algorithm for Gaussian mixture regression

em_gmr <- function(A,
                   B,
                   y,
                   beta_g,
                   beta,
                   sigma_g,
                   pi_g,
                   method = c("sqr", "qr"),
                   max_iter = 2000,
                   tol = 1e-6,
                   sigma_floor = NULL,
                   verbose = FALSE) {
  
  method <- match.arg(method)
  
  A <- as.matrix(A)
  B <- as.matrix(B)
  y <- as.numeric(y)
  beta_g <- as.matrix(beta_g)
  beta <- as.numeric(beta)
  sigma_g <- as.numeric(sigma_g)
  pi_g <- as.numeric(pi_g)
  
  n <- nrow(A)
  G <- nrow(beta_g)
  
  if (length(y) != n) {
    stop("length(y) must equal nrow(A).")
  }
  
  if (nrow(B) != n) {
    stop("nrow(B) must equal nrow(A).")
  }
  
  if (ncol(A) != ncol(beta_g)) {
    stop("ncol(A) must equal ncol(beta_g).")
  }
  
  if (length(beta) != ncol(B)) {
    stop("length(beta) must equal ncol(B).")
  }
  
  if (length(sigma_g) != G) {
    stop("length(sigma_g) must equal nrow(beta_g).")
  }
  
  if (length(pi_g) != G) {
    stop("length(pi_g) must equal nrow(beta_g).")
  }
  
  if (any(sigma_g <= 0)) {
    stop("All initial sigma_g values must be positive.")
  }
  
  if (any(pi_g < 0)) {
    stop("Initial pi_g values must be nonnegative.")
  }
  
  if (sum(pi_g) <= 0) {
    stop("Initial pi_g values must sum to a positive value.")
  }
  
  # Normalize mixing proportions just in case
  pi_g <- pi_g / sum(pi_g)
  
  loglik_trace <- numeric(max_iter)
  converged <- FALSE
  
  # Initial log-likelihood before any EM updates
  loglik_old <- obs_loglik_gmr(
    A = A,
    B = B,
    y = y,
    beta_g = beta_g,
    beta = beta,
    sigma_g = sigma_g,
    pi_g = pi_g
  )
  
  for (iter in seq_len(max_iter)) {
    
    # E-step: compute posterior probabilities
    tau <- e_step_gmr(
      A = A,
      B = B,
      y = y,
      beta_g = beta_g,
      beta = beta,
      sigma_g = sigma_g,
      pi_g = pi_g
    )
    
    # M-step: update parameters using selected method
    if (method == "qr") {
      m_fit <- m_step_qr(
        A = A,
        B = B,
        y = y,
        tau = tau,
        sigma_floor = sigma_floor
      )
    } else if (method == "sqr") {
      m_fit <- m_step_sqr(
        A = A,
        B = B,
        y = y,
        tau = tau,
        sigma_floor = sigma_floor
      )
    }
    
    beta_g <- m_fit$beta_g
    beta <- m_fit$beta
    sigma_g <- m_fit$sigma
    pi_g <- m_fit$pi
    
    # Observed-data log-likelihood after M-step
    loglik_new <- obs_loglik_gmr(
      A = A,
      B = B,
      y = y,
      beta_g = beta_g,
      beta = beta,
      sigma_g = sigma_g,
      pi_g = pi_g
    )
    
    loglik_trace[iter] <- loglik_new
    
    if (verbose) {
      message(
        "iter = ", iter,
        ", loglik = ", round(loglik_new, 6),
        ", diff = ", round(abs(loglik_new - loglik_old), 8)
      )
    }
    
    # Convergence check
    if (abs(loglik_new - loglik_old) < tol) {
      converged <- TRUE
      break
    }
    
    loglik_old <- loglik_new
  }
  
  loglik_trace <- loglik_trace[seq_len(iter)]
  
  # Recompute final tau using final parameter estimates
  tau <- e_step_gmr(
    A = A,
    B = B,
    y = y,
    beta_g = beta_g,
    beta = beta,
    sigma_g = sigma_g,
    pi_g = pi_g
  )
  
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
    method = method,
    G = G
  )
}