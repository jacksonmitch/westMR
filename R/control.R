#' Create a Control Object for westMR Models
#'
#' Creates and validates a control object of tuning parameters for
#' \code{westMR()}. Most users only need \code{alpha}, \code{direction}, and
#' \code{max_iter}/\code{tol} (EM convergence). The remaining parameters goveren
#' initialization strategies, inner IRWLS loops, and numerical safeguards. They
#' are available for tuning but the defaults are reasonable for most fits.
#'
#' @param alpha A numeric value between 0 and 1 specifying the significance
#'  level. Default is 0.05.
#' @param max_iter An integer specifying the maximum number of iterations
#'  allowed for the algorithm. Must be at least 1. Default is 300.
#' @param n_init An integer specifying the total number of initializations
#'  to use. Must be at least 1. Default is 10.
#' @param n_best_init An integer specifying the total number of initializations
#'  to carry into full convergence. Must be at least 1. Default is 1.
#' @param direction A character string specifying the testing direction.
#'   Defaults to 'forward'.
#' @param verbose A logical flag. If \code{TRUE}, detailed execution logs are
#'  printed to the console during execution Default is \code{FALSE}.
#' @param tol A numeric value specifying the iteration tolerance threshold.
#'  Must be greater than or equal to 0. Default is 1e-6.
#' @param n_kmeans_init An integer specifying how many of the total
#'  initializations (\code{n_init}) should be seeded using K-means clustering.
#'  Cannot exceed \code{n_init}. Default is 2.
#' @param kmeans_starts An integer specifying the number of random starts to use
#'  within the K-means algorithm itself. Must be at least 1. Default is 20.
#' @param sigma_floor An optional numeric value establishing a lower bound for
#'  variance estimates to prevent numerical instabilities. Must be greater than
#'  or equal to 0. If \code{NULL}, defaults to a percentage of the response's
#'  variance.
#' @param irwls_max_iter An integer specifying the maximum number of inner
#'  IRWLS iterations used in the Poisson/binomial M-step. Must be at least
#'  1. Default is 50.
#' @param irwls_tol A numeric value specifying the convergence tolerance for
#'  the inner IRWLS loop used in the Poisson/binomial M-step; iteration
#'  stops when the largest absolute change in coefficients falls below this
#'  value. Must be greater than or equal to 0. Default is 1e-8.
#' @param weight_floor A numeric value establishing a lower bound for the
#'  working weights used in the IRWLS loop and the plain Poisson M-step, to
#'  prevent numerical instabilities. Must be greater than or equal to 0.
#'  Default is 1e-10.
#' @param return_qr_parts A logical flag. If \code{TRUE}, the structured-QR
#'  M-step retains its intermediate QR decomposition components in its
#'  output instead of discarding them. Default is \code{FALSE}.
#' @param init_burnin An integer specifying the number of EM iterations run
#'  for each candidate initialization during the burn-in stage used to pick
#'  the best starting point. Must be at least 1. Default is 10.
#' @param init_eps A numeric value controlling how close an initial
#'  responsibility matrix (\code{tau}) built from a cluster or quantile
#'  partition is to a hard assignment: the assigned component gets
#'  probability \code{1 - init_eps} and the rest share \code{init_eps}.
#'  Must be between 0 and 0.5. Default is 1e-6.
#' @param init_min_size An integer specifying a minimum group size to
#'  enforce on initializations. Currently validated but not used by the
#'  estimation code. Must be at least 1. Default is \code{NULL}.
#' @param use_mclust A logical flag reserved for enabling model-based
#'  (mclust) initialization. Currently validated but not used by the
#'  estimation code, which uses k-means and quantile-based initializations
#'  only. Default is \code{TRUE}.
#' @param parallel A logical flag for enabling user-friendly parallel
#'  computation. If \code{TRUE}, a general \code{future::multisession} plan overrides
#'  the current one, which is restored on exit. A \code{future} can be
#'  specified manually while setting this flag as \code{FALSE} for
#'  hardware-specific customization. Default is \code{FALSE}
#'
#' @return A structured list of class \code{"WMRControl"} containing all
#'   validated control arguments.
#'
#' @export
#'
#' @examples
#' # Generate default control options
#' ctrl <- build_control()
#'
#' # Customize specific hyper-parameters
#' custom_ctrl <- build_control(max_iter = 500, verbose = TRUE)
build_control <- function(
  alpha = 0.05,
  max_iter = 300,
  n_init = 10,
  n_best_init = 2,
  direction = "forward",
  verbose = FALSE,
  tol = 1e-6,
  n_kmeans_init = 3,
  kmeans_starts = 20,
  sigma_floor = NULL,
  irwls_max_iter = 50,
  irwls_tol = 1e-8,
  weight_floor = 1e-10,
  init_burnin = 10,
  init_eps = 1e-6,
  init_min_size = NULL,
  use_mclust = TRUE,
  return_qr_parts = FALSE,
  parallel = FALSE
) {
  # User Input checks

  collection <- checkmate::makeAssertCollection()

  checkmate::assert_number(alpha, lower = 0, upper = 1, add = collection)
  checkmate::assert_int(max_iter, lower = 1, add = collection)
  checkmate::assert_int(n_init, lower = 1, add = collection)
  checkmate::assert_int(n_best_init, lower = 1, add = collection)
  if (n_best_init > n_init) {
    collection$push(sprintf(
      "n_best_init (number of initializations carried into full EM) cannot
      exceed n_init. (received %s and %s)", n_best_init, n_init
    ))
  }
  checkmate::assert_choice(direction,
    choices = c("forward", "backward"),
    add = collection
  )
  checkmate::assert_flag(verbose, add = collection)
  checkmate::assert_number(tol, lower = 0, add = collection) # TODO: allows 0
  checkmate::assert_int(n_kmeans_init, lower = 0, add = collection)
  if (n_kmeans_init > n_init) {
    collection$push(sprintf(
      "n_kmeans_init (subset of inits using kmeans) cannot be larger than
      n_init. (received %s and %s)", n_kmeans_init, n_init
    ))
  }
  checkmate::assert_int(kmeans_starts, lower = 1, add = collection)
  checkmate::assert_number(sigma_floor,
    lower = 0,
    add = collection, null.ok = TRUE
  )
  checkmate::assert_int(irwls_max_iter, lower = 1, add = collection)
  checkmate::assert_number(irwls_tol, lower = 0, add = collection)
  checkmate::assert_number(weight_floor, lower = 0, add = collection)
  checkmate::assert_int(init_burnin, lower = 1, add = collection)
  checkmate::assert_number(init_eps, lower = 0, upper = 0.5, add = collection)
  checkmate::assert_int(init_min_size,
    lower = 1, add = collection,
    null.ok = TRUE
  )
  checkmate::assert_flag(use_mclust, add = collection)
  checkmate::assert_flag(parallel, add = collection)

  max_iter <- as.integer(max_iter)
  n_init <- as.integer(n_init)
  n_best_init <- as.integer(n_best_init)
  n_kmeans_init <- as.integer(n_kmeans_init)
  kmeans_starts <- as.integer(kmeans_starts)
  irwls_max_iter <- as.integer(irwls_max_iter)
  init_burnin <- as.integer(init_burnin)

  control <- as.list(environment())
  class(control) <- "WMRControl"
  control
}
