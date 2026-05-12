# Model selection helper functions

count_params_gmr <- function(A, B, G) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  
  p1 <- ncol(A)
  q <- ncol(B)
  k <- G * p1 + q + G + (G - 1)
  
  k
}


compute_bic <- function(loglik, n, k) {
  -2 * loglik + k * log(n)
}
