# determine_effects.R

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
