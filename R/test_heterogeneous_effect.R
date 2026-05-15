# Single heterogeneous-effect test

test_heterogeneous_effect <- function(
    model,
    predictor_tested,
    common,
    direction
){

  if (direction == "forward"){
    null_common <- common
    alt_common <- setdiff(common,predictor_tested)
  }
  else{    # backward
    null_common <- c(common, predictor_tested)
    alt_common <- common
  }

  test <- west_procedure(model = model,
                         null_common = null_common,
                         alt_common = alt_common)

  out <- list(p0 = test$p0,
              reject = test$reject,
              west = test,
              table = test$table)

  class(out) <- "heterogeneous_test"
  out
}
