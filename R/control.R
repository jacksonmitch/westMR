
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
      "n_kmeans_init (subset of inits using kmeans) cannot be larger than n_init. (received %s and %s)", n_kmeans_init, n_init))
  }
  checkmate::assert_int(kmeans_starts, lower = 1, add = collection)
  checkmate::assert_number(sigma_floor, lower = 0, add = collection)

  max_iter <- as.integer(max_iter)
  n_init <- as.integer(n_init)
  n_kmeans_init <- as.integer(n_kmeans_init)
  kmeans_starts <- as.integer(kmeans_starts)

  control <- as.list(environment())
  class(control) <- "WMRControl"
  control
}
