# Single heterogeneous-effect test

test_heterogeneous_effect <- function(
    response,
    covariate,
    heterogeneous,
    homogeneous,
    data,
    G_values,
    alpha = 0.05,
    method = c("sqr", "qr"),
    ...
){
  
  method <- match.arg(method)
  
  null_formula <- make_gmr_formula(response=response, heterogeneous=heterogeneous)
  
  null_common <- make_common_formula(homogeneous=homogeneous)
  
  alt_formula <- make_gmr_formula(response=response, heterogeneous=c(heterogeneous, covariate))
  
  alt_common <- make_common_formula(homogeneous=setdiff(homogeneous, covariate))
  
  test <- west_test(null_formula = null_formula,
                    null_common = null_common,
                    alt_formula = alt_formula,
                    alt_common = alt_common,
                    data = data,
                    G_values = G_values,
                    alpha = alpha,
                    method = method,
                    ...)
  
  out <- list(covariate = covariate,
              test_type = 'heterogeneous',
              alpha = alpha,
              p0 = test$p0,
              reject = test$reject,
              heterogeneous_before = heterogeneous,
              homogeneous_before = homogeneous,
              null_formula = null_formula,
              null_common = null_common,
              alt_formula = alt_formula,
              alt_common = alt_common,
              west = test,
              table = test$table)
  
  class(out) <- "heterogeneous_test"
  out
}