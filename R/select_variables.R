# select_variables.R

select_variables <- function(
  model,
  direction = "forward",
  predictors = model$predictors
) {
  stopifnot(direction %in% c("forward", "backward"))

  if (direction == "forward") {
    included <- character(0)
    update_included <- c
  } else {
    included <- predictors
    update_included <- setdiff
  }

  result <- test_predictors(
    model           = model,
    predictors      = predictors,
    included        = included,
    heterogeneous   = included, # selection: het always mirrors included
    update_included = update_included,
    direction       = direction
  )

  out <- list(
    direction = direction,
    alpha = model$control$alpha,
    G_values = model$G_values,
    all_predictors = predictors,
    selected = result$included,
    steps = result$steps,
    final_fits = result$final_fits,
    final_formula = make_formula(result$included, response = model$response),
    call = match.call()
  )

  class(out) <- "select_variables"
  out
}
