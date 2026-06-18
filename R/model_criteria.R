# Model selection helper functions

count_params_gmr <- function(ncol_het,
                             ncol_common,
                             G,
                             family = c("gaussian", "poisson", "binomial")) {
  family <- match.arg(family)

  k <- G * ncol_het + ncol_common + (G - 1)

  if (family == "gaussian") {
    k <- k + G
  }

  k
}


compute_bic <- function(loglik, n, k) {
  -2 * loglik + k * log(n)
}
