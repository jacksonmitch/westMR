#' @importFrom R6 R6Class
WMRModel <- R6::R6Class(
  "WMRModel",
  public = list(
    data = NULL,
    G_values = NULL,
    family = NULL,
    control = NULL,
    initialize = function(formula,
                          data,
                          G_values,
                          family,
                          control) {
      self$formula <- formula # goes through active binding
      self$data <- data
      self$G_values <- G_values
      self$family <- family

      if (is.null(control$sigma_floor)) {
        mf <- stats::model.frame(formula, data, na.action = stats::na.fail)
        resp <- as.numeric(stats::model.response(mf))
        control$sigma_floor <- 0.05 * stats::sd(resp)
      }
      self$control <- control
    },
    print = function(...) {
      cat(
        "<WMRModel>  formula =", format(self$formula),
        " predictors =", paste(self$predictors, collapse = ", "),
        " G_values =", paste(self$G_values, collapse = ", "),
        " family =", self$family, "\n"
      )
      invisible(self)
    }
  ),
  private = list(
    formula_ = NULL,
    predictors_ = NULL,
    response_ = NULL
  ),
  active = list(
    formula = function(value) {
      if (missing(value)) {
        return(private$formula_)
      }
      private$formula_ <- value
      private$predictors_ <- get_predictors(value)
      private$response_ <- get_response(value)
      invisible(value)
    },
    predictors = function(value) {
      if (!missing(value)) stop("predictors is derived from formula;
                                set formula instead.")
      private$predictors_
    },
    response = function(value) {
      if (!missing(value)) stop("response is derived from formula;
                                set formula instead.")
      private$response_
    }
  )
)
