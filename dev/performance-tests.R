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
test_data <- twelve_vars_data
westMR(formula(test_data), test_data,
    control = build_control(parallel = TRUE))

bench <- microbenchmark::microbenchmark(
  sequential = westMR(formula(test_data), test_data,
    control = build_control()),
  parallel = westMR(formula(test_data), test_data,
    control = build_control(parallel = TRUE)),
  times = 1
)
print(bench)

profvis::profvis({
  westMR(formula(twelve_vars_data), twelve_vars_data,
    control = build_control())
})

westMR(formula(twelve_vars_data), twelve_vars_data,
    control = build_control(parallel = TRUE))

t0 <- Sys.time()
future::plan(future::multisession, workers = 10)
cat("Pool setup:", Sys.time() - t0, "\n")
# ... run test_predictors ...
t1 <- Sys.time()
future::plan(future::sequential)
cat("Teardown:", Sys.time() - t1, "\n")