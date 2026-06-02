test_that("determine_effects", {
  sim_data <- simulate_effect_data(n = 400, seed = 123)
  true_heterogeneous <- "x1"
  true_homogeneous <- c("x2", "x3")
  selected_covariates <- c("x1", "x2", "x3")
  formula <- make_formula(predictors = selected_covariates, response = "y")

  model <- WMRModel$new(
    formula = formula,
    data = sim_data,
    G_values = 2,
    family = "gaussian",
    control = build_control()
  )

  expect_no_error(determine_effects(model = model, direction = "forward"))
  expect_no_error(determine_effects(model = model, direction = "backward"))
})
