#' Create a Control Object for westMR Models
#'
#' This helper function structures and validates tuning parameters,
#' initialization settings, and convergence criteria for the estimation algorithm.
#'
#' @param alpha A numeric value between 0 and 1 specifying the significance level.
#'   Default is 0.05.
#' @param max_iter An integer specifying the maximum number of iterations allowed
#'   for the algorithm. Must be at least 1. Default is 300.
#' @param n_init An integer specifying the total number of random initializations
#'   to attempt. Must be at least 1. Default is 10.
#' @param verbose A logical flag. If \code{TRUE}, detailed execution logs are printed
#'   to the console during optimization. Default is \code{FALSE}.
#' @param tol A numeric value specifying the convergence tolerance threshold.
#'   Must be greater than or equal to 0. Default is 1e-6.
#' @param n_kmeans_init An integer specifying how many of the total initializations
#'   (\code{n_init}) should be seeded using K-means clustering. Cannot exceed
#'   \code{n_init}. Default is 2.
#' @param kmeans_starts An integer specifying the number of random starts to use
#'   within the K-means algorithm itself. Must be at least 1. Default is 20.
#' @param sigma_floor An optional numeric value establishing a lower bound for variance
#'   estimates to prevent numerical instabilities. Must be greater than or equal to 0.
#'   Defaults to \code{NULL} (no floor).
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
  verbose = FALSE,
  tol = 1e-6,
  n_kmeans_init = 2,
  kmeans_starts = 20,
  sigma_floor = NULL
) {
  # User Input checks

  collection <- checkmate::makeAssertCollection()

  checkmate::assert_number(alpha, lower = 0, upper = 1, add = collection)
  checkmate::assert_int(max_iter, lower = 1, add = collection)
  checkmate::assert_int(n_init, lower = 1, add = collection)
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

  max_iter <- as.integer(max_iter)
  n_init <- as.integer(n_init)
  n_kmeans_init <- as.integer(n_kmeans_init)
  kmeans_starts <- as.integer(kmeans_starts)

  control <- as.list(environment())
  class(control) <- "WMRControl"
  control
}
