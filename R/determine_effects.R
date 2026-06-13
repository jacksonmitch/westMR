# Effect-type determination

determine_effects <- function(
  model,
  direction = "forward",
  predictors = model$predictors
) {
  stopifnot(direction %in% c("forward", "backward"))

  control <- model$control
  alpha <- control$alpha

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

  # Fit the initial shared side once before the loop
  shared_data <- prepare_data(model,
    included = predictors,
    common   = common_predictors
  )
  shared_fits <- fit_across_G(model, shared_data)

  steps <- list()
  step_id <- 1L

  remaining_predictors <- predictors
  while (length(remaining_predictors) > 0) {
    tests <- lapply(remaining_predictors, function(predictor) {
      candidate_common <- update_common(common_predictors,
                                        predictor)

      candidate_data <- prepare_data(model,
        included = predictors,
        common = candidate_common
      )

      west_procedure(
        model = model,
        direction = direction,
        shared_fits = shared_fits,
        candidate_data = candidate_data
      )
    })
    names(tests) <- remaining_predictors

    p0 <- vapply(tests, function(t) t$p0, numeric(1))

    steps[[step_id]] <- list(
      step = step_id,
      direction = direction,
      p0 = p0,
      chosen = NA_character_,
      tests = tests
    )

    eligible_p <- p0[is_eligible(p0)]
    if (length(eligible_p) > 0) {
      chosen <- names(eligible_p)[find_best(eligible_p)]
      steps[[step_id]]$chosen <- chosen
    } else {
      break
    }
    # The chosen candidate's fits become the next step's shared fits
    shared_fits <- tests[[chosen]]$candidate_fits
    common_predictors <- update_common(common_predictors, chosen)
    remaining_predictors <- setdiff(remaining_predictors, chosen)

    step_id <- step_id + 1L
  }

  # The shared fits at loop exit are the final model fits across G
  final_fits <- shared_fits

  out <- list(
    direction = direction,
    alpha = alpha,
    G_values = model$G_values,
    heterogeneous = setdiff(predictors, common_predictors),
    homogeneous = common_predictors,
    steps = steps,
    final_fits = final_fits,
    final_formula = model$formula,
    final_common = make_formula(common_predictors),
    call = match.call()
  )

  class(out) <- "determine_effects"
  out
}
