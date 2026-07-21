rm(list = ls())
devtools::load_all()


test_data <- do.call(simulate_fmr,
                     c(scenarios$two_group_effects,
                       list(n = 400, seed = 123)))$data


model <- WMRModel$new(
  formula = formula(test_data),
  data = test_data,
  G_values = 2:3,
  family = "gaussian",
  control = build_control(n_init = 5)
)

forward_fit <- determine_effects(
  model = model,
  direction = "forward"
)
print(forward_fit)
#
# backward_fit <- determine_effects(
#   model = model,
#   direction = "backward"
# )
# print(backward_fit)

bench <- microbenchmark::microbenchmark(
  forward_fit = determine_effects(
    model = model,
    direction = "forward"
  ),
  # backward_fit = determine_effects(
  #   model = model,
  #   direction = "backward"
  # ),
  # westMR_test = westMR(
  #   formula = formula,
  #   data = test_data,
  #   G_max = 3,
  #   family = "gaussian",
  #   task = "effects",
  #   control = build_control(n_init = 5)
  # ),
  # westMR_sim = westMR(
  #   formula = formula,
  #   data = sim_data,
  #   G_max = 3,
  #   family = "gaussian",
  #   task = "effects",
  #   control = build_control(n_init = 5)
  # ),
  times = 10
)
print(bench)

# profvis::profvis({
#   westMR_test = westMR(
#     formula = formula,
#     data = sim_data,
#     G_max = 3,
#     family = "gaussian",
#     task = "effects",
#     control = build_control(n_init = 5)
#   )
# })
