library(R6)

# Class for storing all constant variables (variables, data, settings)
GMRModel <- R6Class(
  "GMRModel",

  public = list(

    formula = NULL,
    response = NULL,
    predictors = NULL,
    G_values = NULL,
    data = NULL,
    method = NULL,
    alpha = NULL,
    n_init = NULL,
    n_kmeans_init = NULL,
    maxit = NULL,
    tol = NULL,
    verbose = NULL,
    seed = NULL,
    fit_function = NULL,

    initialize = function(formula = NULL,
                          G_values = NULL,
                          data = NULL,
                          method = NULL,
                          alpha = 0.05,
                          n_init = NULL,
                          n_kmeans_init = NULL,
                          maxit = NULL,
                          tol = NULL,
                          verbose = NULL,
                          seed = NULL,
                          fit_function = fit_gmr) {
      # Important checks that could be hard to debug
      checkmate::assert_function(fit_function)

      self$formula <- formula
      self$response <- get_response(formula)
      self$predictors <- get_predictors(formula)
      self$G_values <- G_values
      self$data <- data
      self$method <- method
      self$alpha <- alpha
      self$n_init <- n_init
      self$n_kmeans_init <- n_kmeans_init
      self$maxit <- maxit
      self$tol <- tol
      self$verbose <- verbose
      self$seed <- seed
      self$fit_function <- fit_function
    },

    # TODO: @included = predictors included in the model
    # @common = predictors that will be fir under an homogeneous effect
    # @G_values = all candidate orders fitted
    fit = function(included = NULL, common = NULL, ...) {
      tryCatch(
        {
          fit_result <- self$fit_function(
            formula = self$formula,
            common = common,
            data = self$data,
            G_values = self$G_values,
            method = self$method,
            ...)
        },
        error = function(e) {
          stop("Model fitting failed: ",e$message, call. = FALSE)
        }
      )
      fit_result
    }
  )
)


