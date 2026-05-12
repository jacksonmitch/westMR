# Helper used by the augmented QR M-step

build_bigX <- function(A, B, G) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  
  n <- nrow(A)
  p1 <- ncol(A)
  q <- ncol(B)
  
  X <- matrix(0, n * G, G * p1 + q)
  
  idx_g <- function(g) {
    ((g - 1) * n + 1):(g * n)
  }
  
  col_g <- function(g) {
    ((g - 1) * p1 + 1):(g * p1)
  }
  
  for (g in seq_len(G)) {
    ii <- idx_g(g)
    jj <- col_g(g)
    
    X[ii, jj] <- A
    
    if (q > 0) {
      common_cols <- (G * p1 + 1):(G * p1 + q)
      X[ii, common_cols] <- B
    }
  }
  
  X
}

m_step_qr <- function(A, B, y, tau, sigma_floor = NULL) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  y <- as.numeric(y)
  tau <- as.matrix(tau)
  
  n <- nrow(A)
  G <- ncol(tau)
  p1 <- ncol(A)
  q <- ncol(B)
  
  if (length(y) != n) {
    stop("length(y) must equal nrow(A).")
  }
  
  if (nrow(B) != n) {
    stop("nrow(B) must equal nrow(A).")
  }
  
  if (nrow(tau) != n) {
    stop("nrow(tau) must equal nrow(A).")
  }
  
  # Update mixing proportions
  pi_g <- colMeans(tau)
  
  # Stack y and weights in the group-by-group order
  y_stack <- rep(y, times = G)
  w_stack <- as.numeric(tau)
  
  # Build full augmented design matrix
  X_big <- build_bigX(A, B, G)
  
  # Weighted least squares fit
  fit <- lm.wfit(
    x = X_big,
    y = y_stack,
    w = w_stack
  )
  
  c_hat <- as.numeric(fit$coefficients)
  
  # Extract component-specific coefficients
  beta_g <- matrix(
    c_hat[seq_len(G * p1)],
    nrow = G,
    byrow = TRUE
  )
  
  # Extract common-effect coefficients
  if (q > 0) {
    beta <- c_hat[(G * p1 + 1):(G * p1 + q)]
  } else {
    beta <- numeric(0)
  }
  
  # Update sigma_g after coefficient update.
  mu <- sweep(
    A %*% t(beta_g),
    1,
    common_eta(B, beta),
    "+"
  )
  
  res2 <- (y - mu)^2
  den <- pmax(1e-8, colSums(tau))
  
  if (is.null(sigma_floor)) {
    sigma_floor <- 0.05 * stats::sd(y)
  }
  
  sigma_g <- pmax(
    sigma_floor,
    sqrt(colSums(tau * res2) / den)
  )
  
  list(
    beta = beta,
    beta_g = beta_g,
    pi = pi_g,
    sigma = sigma_g,
    fitted = fit$fitted.values
  )
}