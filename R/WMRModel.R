#' Container for One westMR Model Specification
#'
#' An R6 class holding the formula, data, candidate \code{G} values, family,
#' and control settings for a \code{westMR} run. \code{predictors} and
#' \code{response} are derived from \code{formula} via active bindings and
#' cannot be set directly.
#'
#' @param formula A formula object specifying the model.
#' @param data A data.frame containing the variables in \code{formula}.
#' @param G_values An integer vector of candidate numbers of mixture
#'   components.
#' @param family A character string specifying the error distribution:
#'   \code{"gaussian"}, \code{"poisson"}, or \code{"binomial"}.
#' @param control A \code{WMRControl} object (from \code{build_control()}).
#'   If \code{control$sigma_floor} is \code{NULL}, it is set here to
#'   \code{0.05 * sd(response)}.
#'
#' @return A new \code{WMRModel} object.
#' @importFrom R6 R6Class
#' @noRd
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
