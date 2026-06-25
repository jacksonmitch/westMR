EmState <- R6::R6Class("EmState",
  private = list(
    .tau              = NULL,
    .pi_g             = NULL,
    .beta_g           = NULL,
    .beta             = NULL,
    .sigma_g          = NULL,
    .eta              = NULL,
    .loglik           = NULL,
    .irwls_iterations = NULL,
    .irwls_converged  = NULL
  ),
  active = list(
    tau = function(value) {
      if (missing(value)) {
        return(private$.tau)
      }
      stopifnot(is.matrix(value))
      private$.tau <- value
    },
    G = function(value) {
      if (missing(value)) {
        return(ncol(private$.tau))
      }
      stop("G is read-only; it is derived from tau.")
    },
    pi_g = function(value) {
      if (missing(value)) {
        return(private$.pi_g)
      }
      stopifnot(is.numeric(value))
      private$.pi_g <- value
    },
    # These parameters can be set/start as NULL
    beta_g = function(value) {
      if (missing(value)) {
        return(private$.beta_g)
      }
      stopifnot(is.null(value) || is.matrix(value))
      private$.beta_g <- value
    },
    beta = function(value) {
      if (missing(value)) {
        return(private$.beta)
      }
      stopifnot(is.null(value) || is.numeric(value))
      private$.beta <- value
    },
    sigma_g = function(value) {
      if (missing(value)) {
        return(private$.sigma_g)
      }
      stopifnot(is.null(value) || is.numeric(value))
      private$.sigma_g <- value
    },
    eta = function(value) {
      if (missing(value)) {
        return(private$.eta)
      }
      stopifnot(is.null(value) || is.numeric(value))
      private$.eta <- value
    },
    loglik = function(value) {
      if (missing(value)) {
        return(private$.loglik)
      }
      stopifnot(is.numeric(value), length(value) == 1)
      private$.loglik <- value
    },
    irwls_iterations = function(value) {
      if (missing(value)) {
        return(private$.irwls_iterations)
      }
      stopifnot(is.numeric(value), length(value) == 1)
      private$.irwls_iterations <- value
    },
    irwls_converged = function(value) {
      if (missing(value)) {
        return(private$.irwls_converged)
      }
      stopifnot(is.logical(value), length(value) == 1)
      private$.irwls_converged <- value
    }
  ),
  public = list(
    initialize = function(tau,
                          pi_g = NULL,
                          beta_g = NULL,
                          beta = NULL,
                          sigma_g = NULL,
                          eta = NULL,
                          loglik = -Inf,
                          irwls_iterations = NA_integer_,
                          irwls_converged = FALSE) {
      self$tau <- tau
      self$pi_g <- if (is.null(pi_g)) colMeans(tau) else pi_g

      self$beta_g <- beta_g
      self$beta <- beta
      self$sigma_g <- sigma_g
      self$eta <- eta
      self$loglik <- loglik

      self$irwls_iterations <- irwls_iterations
      self$irwls_converged <- irwls_converged
    }
  )
)
