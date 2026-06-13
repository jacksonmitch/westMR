# determine_effects.R

determine_effects <- function(
  model,
  direction = "forward",
  predictors = model$predictors
) {
  stopifnot(direction %in% c("forward", "backward"))

  if (direction == "forward") {
    common_predictors <- predictors
    is_eligible <- function(p) is.finite(p) & p < alpha
    find_best <- which.min
    update_common <- function(common, predictor) {
      setdiff(common, predictor)
    }
  } else {
    common_predictors <- character(0)
    is_eligible <- function(p) is.finite(p) & p >= alpha
    find_best <- which.max
    update_common <- function(common, predictor) {
      unique(c(common, predictor))
    }
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

  out <- list(
    direction = direction,
    alpha = model$control$alpha,
    G_values = model$G_values,
    heterogeneous = result$heterogeneous,
    homogeneous = homogeneous,
    steps = result$steps,
    final_fits = result$final_fits,
    final_formula = model$formula,
    final_common = make_formula(homogeneous),
    call = match.call()
  )

  class(out) <- "determine_effects"
  out
}
