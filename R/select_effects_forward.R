# Forward effect-type determination

select_effects_forward <- function(
    response,
    covariates,
    data,
    G_values,
    alpha = 0.05,
    method = c("sqr", "qr"),
    ...
){
  
  method <- match.arg(method)
  
  heterogeneous <- character(0)
  homogeneous <- unique(covariates)
  steps <- list()
  step_id <- 1
  continue <- TRUE
  
  while(length(homogeneous) > 0 && continue){
    
    tests <- lapply(homogeneous, function(v){
      test_heterogeneous_effect(
        response = response,
        covariate = v,
        heterogeneous = heterogeneous,
        homogeneous = homogeneous,
        data = data,
        G_values = G_values,
        alpha = alpha,
        method = method,
        ...)
    })
    
    names(tests) <- homogeneous
    p0 <- vapply(tests, function(u) u$p0, numeric(1))
    
    if(all(!is.finite(p0))){
      chosen <- NA_character_
      reject <- FALSE
    }
    else{
      best_idx <- which.min(p0)
      chosen <- names(p0)[best_idx]
      reject <- is.finite(p0[best_idx]) && p0[best_idx] < alpha
    }
    
    steps[[step_id]] <- list(
      step = step_id,
      direction = "forward",
      heterogeneous_before = heterogeneous,
      homogeneous_before = homogeneous,
      tested = homogeneous,
      p0 = p0,
      chosen = ifelse(reject, chosen, NA_character_),
      tests = tests
    )
    
    if(reject){
      heterogeneous <- c(heterogeneous, chosen)
      homogeneous <- setdiff(homogeneous, chosen)
      step_id <- step_id + 1
    }
    else{
      continue <- FALSE
    }
  }
  
  final_formula <- make_gmr_formula(response, heterogeneous)
  final_common <- make_common_formula(homogeneous)
  
  final_fit <- fit_gmr(
    formula = final_formula,
    common = final_common,
    data = data,
    G_values = G_values,
    method = method,
    ...
  )
  
  out <- list(
    direction = "forward",
    response = response,
    covariates = covariates,
    heterogeneous = heterogeneous,
    homogeneous = homogeneous,
    alpha = alpha,
    G_values = G_values,
    method = method,
    steps = steps,
    final_formula = final_formula,
    final_common = final_common,
    final_fit = final_fit,
    call = match.call()
  )
  
  class(out) <- "effect_selection"
  out
}