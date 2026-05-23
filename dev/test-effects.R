rm(list = ls())
devtools::load_all()


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

# forward_fit <- determine_effects(
#   model = model,
#   direction = "forward"
# )
# print(forward_fit)
# summary(forward_fit$final_fit)
#
# steps_df <- effect_selection_table(forward_fit)
# steps_df
#
bench <- microbenchmark::microbenchmark(
  forward_fit = determine_effects(
    model = model,
    direction = "forward"
  ),
  backward_fit = determine_effects(
    model = model,
    direction = "backward"
  ),
  times = 1
)
print(bench)
#
# backward_fit <- determine_effects(
#   model = model,
#   direction = "backward"
# )
# print(backward_fit)
# summary(backward_fit$final_fit)

# profvis::profvis({
#   forward_fit = determine_effects(
#     model = model,
#     direction = "forward"
#   )
# })

# test_bench <- microbenchmark::microbenchmark(
#   model = model$fit(),
#   raw_call = fit_fmr(formula, data = data, G_values = model$G_values),
#   times = 100
# )
# print(test_bench)



