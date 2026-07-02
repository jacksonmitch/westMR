# determine_effects.R

#' Stepwise Determination of Effect Types
#'
#' Runs the WEST stepwise procedure to decide, for each included predictor,
#' whether its effect is heterogeneous (component-specific) or homogeneous
#' (shared across mixture components). At each step, remaining predictors
#' are tested against a shared baseline fit via \code{test_predictors()}.
#'
#' @param model A \code{WMRModel} object.
#' @param direction A character string, either \code{"forward"} or
#'   \code{"backward"}, specifying the stepwise search direction.
#' @param predictors A character vector of predictor names to test for
#'   heterogeneity. Defaults to all predictors in \code{model}.
#'
#' @return A list of class \code{determine_effects} with elements:
#'   \code{direction}, \code{alpha}, \code{G_values}, \code{best_fit} (the
#'   fit at the BIC-optimal \code{G} for the final model),
#'   \code{heterogeneous} and \code{homogeneous} predictor names,
#'   \code{steps} (the step-by-step search log), \code{final_fits} (fits
#'   across \code{G_values} for the final model), \code{final_formula}, and
#'   \code{call}.
#' @noRd
determine_effects <- function(
  model,
  direction = "forward",
  predictors = model$predictors
) {
  stopifnot(direction %in% c("forward", "backward"))

  if (direction == "forward") {
    heterogeneous <- character(0)
  } else {
    heterogeneous <- predictors
  }

  result <- test_predictors(
    model           = model,
    predictors      = predictors,
    included        = predictors, # effects: included is always fixed
    heterogeneous   = heterogeneous,
    update_included = function(s, x) s, # effects: included never changes
    direction       = direction
  )

  homogeneous <- setdiff(predictors, result$heterogeneous)

  best_fit <- select_best_G(result$final_fits, criterion = "bic")

  out <- list(
    direction = direction,
    alpha = model$control$alpha,
    G_values = model$G_values,
    best_fit = best_fit,
    heterogeneous = result$heterogeneous,
    homogeneous = homogeneous,
    steps = result$steps,
    final_fits = result$final_fits,
    final_formula = model$formula,
    call = match.call()
  )

  class(out) <- "determine_effects"
  out
}
