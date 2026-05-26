# Model selection helper functions

count_params_gmr <- function(ncol_het, ncol_common, G) {
  k <- G * ncol_het + ncol_common + G + (G - 1)

  k
}


compute_bic <- function(loglik, n, k) {
  -2 * loglik + k * log(n)
}
