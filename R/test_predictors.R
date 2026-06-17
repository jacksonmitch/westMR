# test_predictors.R

test_predictors <- function(
  model,
  predictors,
  included,
  heterogeneous,
  update_included,
  direction = c("forward", "backward")
) {
  direction <- match.arg(direction)
  alpha <- model$control$alpha

  included <- as.character(included)
  heterogeneous <- as.character(heterogeneous)
  predictors <- as.character(predictors)

  heterogeneous <- intersect(heterogeneous, included)

  # These three always co-vary with direction — not task-specific
  if (direction == "forward") {
    is_eligible <- function(p) is.finite(p) & p < alpha
    find_best <- which.min
    update_heterogeneous <- union
    remaining_predictors <- setdiff(predictors, heterogeneous)
  } else {
    is_eligible <- function(p) is.finite(p) & p >= alpha
    find_best <- which.max
    update_heterogeneous <- setdiff
    remaining_predictors <- intersect(predictors, heterogeneous)
  }

  steps <- list()
  step_id <- 1L

  common <- setdiff(included, heterogeneous)
  shared_data <- prepare_data(model, included = included, common = common)
  shared_fits <- fit_across_G(model, shared_data)

  while (length(remaining_predictors) > 0) {
    tests <- lapply(remaining_predictors, function(predictor) {
      cand_included <- update_included(included, predictor)
      cand_heterogeneous <- update_heterogeneous(heterogeneous, predictor)
      cand_common <- setdiff(cand_included, cand_heterogeneous)

      candidate_data <- prepare_data(model,
        included = cand_included,
        common   = cand_common
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

    shared_fits <- tests[[chosen]]$candidate_fits
    included <- update_included(included, chosen)
    heterogeneous <- update_heterogeneous(heterogeneous, chosen)
    remaining_predictors <- setdiff(remaining_predictors, chosen)
    step_id <- step_id + 1L
  }

  list(
    included = included,
    heterogeneous = heterogeneous,
    common = setdiff(included, heterogeneous),
    steps = steps,
    final_fits = shared_fits
  )
}
