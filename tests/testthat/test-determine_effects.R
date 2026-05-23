test_that("determine_effects", {
  data <- simulate_effect_data(n = 400, seed = 123)
  print(class(data))
  true_heterogeneous <- 'x1'
  true_homogeneous <- c('x2','x3')
  selected_covariates <- c('x1','x2','x3')
  formula <- make_formula(predictors = selected_covariates, response = 'y')

  model <- WMRModel$new(
    formula = formula,
    data = data,
    G_values = 1:2,
    m_step_var = m_step_qr,
    m_step_eff = m_step_sqr,
    log_lik = obs_loglik_gmr,
    control = build_control()
  )

  expect_no_error(determine_effects(model = model,direction = "forward"))
  expect_no_error(determine_effects(model = model,direction = "backward"))
})
