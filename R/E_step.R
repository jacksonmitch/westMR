# E-step for Gaussian mixture regression

e_step_gmr <- function(A, B, y, beta_g, beta, sigma_g, pi_g){
  
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
  
  mu <- sweep(
    A %*% t(beta_g), 1, common_eta(B, beta), "+"
  )
  
  res2 <- (y - mu)^2
  
  # log pi_g + log normal density
  log_w <- sweep(-0.5 * res2, 2, sigma_g^2, "/")
  
  log_w <- sweep(
    log_w,
    2,
    log(pi_g + 1e-16) - log(sigma_g + 1e-16),
    "+"
  ) - 0.5 * log(2 * pi)
  
  tau <- row_softmax(log_w)
  
  colnames(tau) <- paste0("g", seq_len(G))
  
  return(tau)
}