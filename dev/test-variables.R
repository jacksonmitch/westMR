rm(list = ls())
devtools::load_all()


test_data <- do.call(simulate_fmr,
                       c(scenarios$three_group_variables,
                         list(n = 500, seed = 1)))

model <- WMRModel$new(
  formula = formula(test_data),
  data = test_data,
  G_values = 1:3,
  family = "gaussian",
  control = build_control()
)

forward_fit <- select_variables(
  model = model,
  direction = "forward"
)
print(forward_fit)

backward_fit <- select_variables(
  model = model,
  direction = "backward"
)
print(backward_fit)
