rm(list = ls())
devtools::load_all()

# 2. Benchmark the execution of a specific test file
bm <- microbenchmark::microbenchmark(
  test_file_time = testthat::test_file("tests/testthat/test-fits.R", reporter = "silent"),
  times = 100
)

# 3. View the results
print(bm)
saveRDS(bm, "matrixStats_speed.rds")


twelve_vars_data <- do.call(simulate_fmr,
                      c(scenarios$three_group_twelve_variables_gaussian,
                        list(n = 500, seed = 1)))
four_vars_data <- do.call(simulate_fmr,
                       c(scenarios$three_group_four_variables,
                         list(n = 500, seed = 1)))
two_g_3_vars_data <- do.call(simulate_fmr,
                       c(scenarios$two_group_effects,
                         list(n = 500, seed = 1)))
simulation_study_data_500 <- do.call(simulate_fmr,
                       c(scenarios$four_group_twelve_variables,
                         list(n = 500, seed = 123)))
simulation_study_data_1k <- do.call(simulate_fmr,
                       c(scenarios$four_group_twelve_variables,
                         list(n = 1000, seed = 123)))

test_data <- simulation_study_data_500
result_500 <- westMR(test_data$formula, test_data$data, G_max = 7,
    control = build_control(parallel = TRUE, direction = "forward"))

result_g3 <- westMR(four_vars_data$formula, four_vars_data$data, G_max = 4,
    control = build_control(parallel = TRUE, direction = "forward"))
result_g2 <- westMR(two_g_3_vars_data$formula, two_g_3_vars_data$data, G_max = 4,
    control = build_control(parallel = TRUE, direction = "forward"))

bench <- microbenchmark::microbenchmark(
  # sequential = westMR(test_data$formula, test_data$data,
  #   control = build_control()),
  parallel = westMR(test_data$formula, test_data$data,
    control = build_control(parallel = TRUE)),
  times = 1
)
print(bench)
