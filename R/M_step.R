# Wrapper for Gaussian, Poisson, and binomial M-steps

m_step_gmr <- function(A,
                       B,
                       y,
                       tau,
                       family = c("gaussian", "poisson", "binomial"),
                       beta_g = NULL,
                       beta = NULL,
                       sigma_floor = NULL,
                       trials = NULL,
                       irwls_maxit = 50,
                       irwls_tol = 1e-8,
                       weight_floor = 1e-10,
                       return_qr_parts = FALSE) {
  
  family <- match.arg(family)
  
  if (family == "gaussian") {
    return(m_step_gaussian(
      A = A,
      B = B,
      y = y,
      tau = tau,
      sigma_floor = sigma_floor,
      return_qr_parts = return_qr_parts
    ))
  }
  
  if (family == "poisson") {
    return(m_step_poisson(
      A = A,
      B = B,
      y = y,
      tau = tau,
      beta_g = beta_g,
      beta = beta,
      irwls_maxit = irwls_maxit,
      irwls_tol = irwls_tol,
      weight_floor = weight_floor
    ))
  }
  
  if (family == "binomial") {
    return(m_step_binomial(
      A = A,
      B = B,
      y = y,
      tau = tau,
      beta_g = beta_g,
      beta = beta,
      trials = trials,
      irwls_maxit = irwls_maxit,
      irwls_tol = irwls_tol,
      weight_floor = weight_floor
    ))
  }
}


# Gaussian M-step using structured QR

m_step_gaussian <- function(A,
                            B,
                            y,
                            tau,
                            sigma_floor = NULL,
                            return_qr_parts = FALSE) {
  
  A <- as.matrix(A)
  
  if (is.null(B)) {
    B <- matrix(numeric(0), nrow = nrow(A), ncol = 0)
  } else {
    B <- as.matrix(B)
  }
  
  y <- as.numeric(y)
  tau <- as.matrix(tau)
  
  pi_g <- colMeans(tau)
  
  fit <- wls_sqr(
    A = A,
    B = B,
    z = y,
    w = tau,
    return_qr_parts = return_qr_parts
  )
  
  beta_g <- fit$beta_g
  beta <- fit$beta
  
  mu <- linear_predictor_matrix(A, B, beta_g, beta)
  res2 <- (y - mu)^2
  
  den <- pmax(colSums(tau), 1e-8)
  
  if (is.null(sigma_floor)) {
    sigma_floor <- 0.05 * stats::sd(y)
  }
  
  sigma_g <- sqrt(colSums(tau * res2) / den)
  sigma_g <- pmax(sigma_g, sigma_floor)
  
  out <- list(
    beta_g = beta_g,
    beta = beta,
    pi_g = pi_g,
    sigma_g = sigma_g,
    irwls_iterations = NA_integer_,
    irwls_converged = NA
  )
  
  if (isTRUE(return_qr_parts)) {
    out$R1_list <- fit$R1_list
    out$E_list <- fit$E_list
    out$R2 <- fit$R2
  }
  
  out
}


# Poisson M-step using inner IRWLS and structured QR

m_step_poisson <- function(A,
                           B,
                           y,
                           tau,
                           beta_g = NULL,
                           beta = NULL,
                           irwls_maxit = 50,
                           irwls_tol = 1e-8,
                           weight_floor = 1e-10) {
  
  tau <- as.matrix(tau)
  pi_g <- colMeans(tau)
  
  fit <- irwls_gmr(
    A = A,
    B = B,
    y = y,
    tau = tau,
    beta_g = beta_g,
    beta = beta,
    trials = NULL,
    maxit = irwls_maxit,
    tol = irwls_tol,
    weight_floor = weight_floor,
    family = "poisson"
  )
  
  list(
    beta_g = fit$beta_g,
    beta = fit$beta,
    pi_g = pi_g,
    sigma_g = NULL,
    irwls_iterations = fit$iterations,
    irwls_converged = fit$converged
  )
}

# Binomial M-step using inner IRWLS and structured QR

m_step_binomial <- function(A,
                            B,
                            y,
                            tau,
                            beta_g = NULL,
                            beta = NULL,
                            trials = NULL,
                            irwls_maxit = 50,
                            irwls_tol = 1e-8,
                            weight_floor = 1e-10) {
  
  tau <- as.matrix(tau)
  pi_g <- colMeans(tau)
  
  if (is.null(trials)) {
    trials <- rep(1, length(y))
  }
  
  fit <- irwls_gmr(
    A = A,
    B = B,
    y = y,
    tau = tau,
    beta_g = beta_g,
    beta = beta,
    trials = trials,
    maxit = irwls_maxit,
    tol = irwls_tol,
    weight_floor = weight_floor,
    family = "binomial"
  )
  
  list(
    beta_g = fit$beta_g,
    beta = fit$beta,
    pi_g = pi_g,
    sigma_g = NULL,
    irwls_iterations = fit$iterations,
    irwls_converged = fit$converged
  )
}