rm(list = ls())
devtools::load_all()


test_data <- do.call(simulate_fmr,
                       c(scenarios$three_group_four_variables,
                         list(n = 500, seed = 1)))

result <- westMR(formula(test_data), test_data)
print(result)
print(summary(result))

model <- WMRModel$new(
  formula = formula(test_data),
  data = test_data,
  G_values = 2:3,
  family = "gaussian",
  control = build_control(n_init = 5)
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

# bench <- microbenchmark::microbenchmark(
#   forward_fit = determine_effects(
#     model = model,
#     direction = "forward"
#   ),
#   backward_fit = determine_effects(
#     model = model,
#     direction = "backward"
#   ),
#   times = 10
# )
# print(bench)

# profvis::profvis({
#   forward_fit = determine_effects(
#     model = model,
#     direction = "forward"
#   )
# })
