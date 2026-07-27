# tests/testthat/test-parallel.R
test_that("westMR runs without error with control$parallel = TRUE", {
  test_data <- do.call(
    simulate_fmr,
    c(
      scenarios$two_group_effects,
      list(n = 400, seed = 123)
    )
  )$data
  expect_no_error(
    westMR(formula(test_data), test_data,
      control = build_control(parallel = TRUE)
    )
  )
})
