#' @importFrom R6 R6Class
WMRModel <- R6::R6Class(
  "WMRModel",
  public = list(
    formula = NULL,
    response = NULL,
    predictors = NULL,
    data = NULL,
    G_values = NULL,
    family = NULL,
    control = NULL,
    initialize = function(formula,
                          data,
                          G_values,
                          family,
                          control) {
      self$formula <- formula
      self$response <- get_response(formula)
      self$predictors <- get_predictors(formula)
      self$data <- data
      self$G_values <- G_values
      self$family <- family
      self$control <- control
    },
    print = function(...) {
      cat(
        "<WMRModel>  formula =", format(self$formula),
        " predictors =", self$predictos,
        " G_values =", self$G_values,
        " family =", self$family, "\n"
      )
      invisible(self)
    }
  )
)
