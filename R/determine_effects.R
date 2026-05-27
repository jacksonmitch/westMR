# Forward effect-type determination

determine_effects<- function(
    model,
    direction = "forward"
){

  predictors <- model$predictors
  control <- model$control
  alpha <- control$alpha

  if (direction == "forward"){
    common_predictors <- predictors
    is_eligible <- function(p) is.finite(p) & p < alpha
    find_best <- which.min
    update_common <- setdiff # removes predictor from list
  }
  else{ #backward
    common_predictors <- NULL
    is_eligible <- function(p) is.finite(p) & p >= alpha
    find_best <- which.max
    update_common <- c # adds predictors to list
  }

  steps <- list()
  step_id <- 1

  while(length(predictors) > 0) {
    tests <- lapply(predictors, function(predictor_tested){
      if (direction == "forward"){
        null_common <- common_predictors
        alt_common <- setdiff(common_predictors,predictor_tested)
      }
      else{    # backward
        null_common <- c(common_predictors, predictor_tested)
        alt_common <- common_predictors
      }
      test <- west_procedure(model = model,
                             null_common = null_common,
                             alt_common = alt_common)
    })

    names(tests) <- predictors
    p0 <- vapply(tests, function(u) u$p0, numeric(1))

    eligible_p <- p0[is_eligible(p0)]

    if (length(eligible_p) == 0) {
      steps[[step_id]] <- list(
        step = step_id,
        direction = direction,
        p0 = p0,
        chosen = NA_character_,
        tests = tests
      )
      break
    }
    chosen <- names(eligible_p)[find_best(eligible_p)]
    predictors <- setdiff(predictors, chosen)
    common_predictors <- update_common(common_predictors, chosen)

    steps[[step_id]] <- list(
      step = step_id,
      direction = direction,
      p0 = p0,
      chosen = chosen,
      tests = tests
    )
    step_id <- step_id + 1L # L for Literal integer
  }

  # final_fit <- fit_fmr(
  #   model,
  #   G,
  #   included = predictors,
  #   common = common_predictors
  # )

  out <- list(
    direction = direction,
    method = model$method,
    alpha = model$alpha,
    G_values = model$G_values,
    heterogeneous = setdiff(model$predictors, common_predictors),
    homogeneous = common_predictors,
    steps = steps,
    final_formula = model$formula,
    final_common = make_formula(common_predictors),
    #final_fit = final_fit,
    call = match.call()
  )

  class(out) <- "effect_selection"
  out
}
