# Backward effect-type determination

select_effects_backward <- function(
    response,
    covariates,
    data,
    G_values,
    alpha = 0.05,
    method = c("sqr", "qr"),
    ...
){
  
  method <- match.arg(method)
  
  heterogeneous <- unique(covariates)
  homogeneous <- character(0)
  steps <- list()
  step_id <- 1
  continue <- TRUE
  
  while (length(heterogeneous) > 0 && continue) {
    
    tests <- lapply(heterogeneous, function(v) {
      null_heterogeneous <- setdiff(heterogeneous, v)
      null_homogeneous <- c(homogeneous, v)
      
      null_formula <- make_gmr_formula(response, null_heterogeneous)
      null_common <- make_common_formula(null_homogeneous)
      
      alt_formula <- make_gmr_formula(response, heterogeneous)
      alt_common <- make_common_formula(homogeneous)
      
      test <- west_test(
        null_formula = null_formula,
        null_common = null_common,
        alt_formula = alt_formula,
        alt_common = alt_common,
        data = data,
        G_values = G_values,
        alpha = alpha,
        method = method,
        ...
      )
      
      out <- list(
        covariate = v,
        test_type = "heterogeneous_effect",
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
        table = test$table
      )
      
      class(out) <- "heterogeneous_test"
      out
    })
    
    names(tests) <- heterogeneous
    
    p0 <- vapply(tests, function(u) u$p0, numeric(1))
    
    nonsig <- is.finite(p0) & p0 >= alpha
    
    if (any(nonsig)) {
      chosen <- names(which.max(ifelse(nonsig, p0, -Inf)))
      move_to_homogeneous <- TRUE
    } else {
      chosen <- NA_character_
      move_to_homogeneous <- FALSE
    }
    
    steps[[step_id]] <- list(
      step = step_id,
      direction = "backward",
      heterogeneous_before = heterogeneous,
      homogeneous_before = homogeneous,
      tested = heterogeneous,
      p0 = p0,
      chosen = ifelse(move_to_homogeneous, chosen, NA_character_),
      tests = tests
    )
    
    if (move_to_homogeneous) {
      homogeneous <- c(homogeneous, chosen)
      heterogeneous <- setdiff(heterogeneous, chosen)
      step_id <- step_id + 1L
    } else {
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
    direction = "backward",
    response = response,
    covariates = covariates,
    heterogeneous = heterogeneous,
    homogeneous = homogeneous,
    no_effect = character(0),
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