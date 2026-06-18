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

      if (is.null(control$sigma_floor)) {
        mf <- stats::model.frame(formula, data, na.action = stats::na.fail)
        response <- as.numeric(stats::model.response(mf))
        control$sigma_floor <- 0.05 * stats::sd(response)
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
  )
)




