# S3 functions 

# Print GMR fit

print.gmr_fit <- function(x, ...) {
  cat("Gaussian mixture regression fit\n")
  cat("--------------------------------\n")
  
  cat("Formula:        ")
  print(x$formula)
  
  if (!is.null(x$common)) {
    cat("Common effects: ")
    print(x$common)
  }
  
  cat("\n")
  cat("Method:         ", x$method, "\n", sep = "")
  cat("Selection:      ", x$criterion, "\n", sep = "")
  cat("Best G:         ", x$best_G, "\n", sep = "")
  cat("Log-likelihood: ", round(x$loglik, 4), "\n", sep = "")
  cat("BIC:            ", round(x$bic, 4), "\n", sep = "")
  cat("Parameters:     ", x$k, "\n", sep = "")
  cat("Iterations:     ", x$iterations, "\n", sep = "")
  cat("Converged:      ", x$converged, "\n", sep = "")
  cat("Best init:      ", x$best_init, " of ", x$n_init, "\n", sep = "")
  cat("\n")
  
  cat("Candidate models:\n")
  print(x$results, row.names = FALSE)
  
  invisible(x)
}

# Summarize GMR fit

summary.gmr_fit <- function(object, ...) {
  
  component_coefficients <- object$beta_g
  common_coefficients <- object$beta
  
  if (!is.null(component_coefficients)) {
    colnames(component_coefficients) <- object$A_colnames
    rownames(component_coefficients) <- paste0("g", seq_len(nrow(component_coefficients)))
  }
  
  if (!is.null(common_coefficients)) {
    names(common_coefficients) <- object$B_colnames
  }
  
  out <- list(
    call = object$call,
    formula = object$formula,
    common = object$common,
    
    method = object$method,
    criterion = object$criterion,
    best_G = object$best_G,
    
    loglik = object$loglik,
    bic = object$bic,
    k = object$k,
    iterations = object$iterations,
    converged = object$converged,
    
    pi_g = object$pi_g,
    sigma_g = object$sigma_g,
    beta_g = component_coefficients,
    beta = common_coefficients,
    
    results = object$results,
    n_init = object$n_init,
    best_init = object$best_init,
    n_valid_init = object$n_valid_init
  )
  
  class(out) <- "summary.gmr_fit"
  
  out
}

# Print GMR fit summary

print.summary.gmr_fit <- function(x, ...) {
  cat("Summary of Gaussian mixture regression fit\n")
  cat("------------------------------------------\n")
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  cat("Model:\n")
  cat("  Formula:        ")
  print(x$formula)
  
  cat("  Common effects: ")
  print(x$common)
  
  cat("\n")
  cat("Fitting information:\n")
  cat("  Method:         ", x$method, "\n", sep = "")
  cat("  Criterion:      ", x$criterion, "\n", sep = "")
  cat("  Selected G:     ", x$best_G, "\n", sep = "")
  cat("  Log-likelihood: ", round(x$loglik, 4), "\n", sep = "")
  cat("  BIC:            ", round(x$bic, 4), "\n", sep = "")
  cat("  Parameters:     ", x$k, "\n", sep = "")
  cat("  Iterations:     ", x$iterations, "\n", sep = "")
  cat("  Converged:      ", x$converged, "\n", sep = "")
  cat("  Best init:      ", x$best_init, " of ", x$n_init, "\n", sep = "")
  cat("  Valid init:     ", x$n_valid_init, " of ", x$n_init, "\n", sep = "")
  
  cat("\n")
  cat("Mixing proportions:\n")
  print(round(x$pi_g, 4))
  
  cat("\n")
  cat("Component standard deviations:\n")
  print(round(x$sigma_g, 4))
  
  cat("\n")
  cat("Component-specific coefficients:\n")
  print(round(x$beta_g, 4))
  
  if (!is.null(x$beta) && length(x$beta) > 0) {
    cat("\n")
    cat("Common-effect coefficients:\n")
    print(round(x$beta, 4))
  }
  
  cat("\n")
  cat("Model selection table:\n")
  print(x$results, row.names = FALSE)
  
  invisible(x)
}

# Print effect selection

print.effect_selection <- function(x, ...) {
  cat("GMR effect-type selection\n")
  cat("-------------------------\n")
  cat("Direction: ", x$direction, "\n", sep = "")
  cat("Method:    ", x$method, "\n", sep = "")
  cat("Alpha:     ", x$alpha, "\n", sep = "")
  cat("G values:  ", paste(x$G_values, collapse = ", "), "\n", sep = "")
  cat("\n")
  
  cat("Heterogeneous effects:\n")
  if (length(x$heterogeneous) == 0) {
    cat("  None\n")
  } else {
    cat("  ", paste(x$heterogeneous, collapse = ", "), "\n", sep = "")
  }
  
  cat("\nHomogeneous effects:\n")
  if (length(x$homogeneous) == 0) {
    cat("  None\n")
  } else {
    cat("  ", paste(x$homogeneous, collapse = ", "), "\n", sep = "")
  }
  
  cat("\nFinal model:\n")
  cat("  Formula: ")
  print(x$final_formula)
  cat("  Common:  ")
  print(x$final_common)
  
  if (!is.null(x$final_fit)) {
    cat("\nSelected G: ")
    if (!is.null(x$final_fit$best_G)) {
      cat(x$final_fit$best_G, "\n", sep = "")
    } else {
      cat("NA\n")
    }
    
    cat("BIC:        ")
    if (is.numeric(x$final_fit$bic) &&
        length(x$final_fit$bic) == 1 &&
        is.finite(x$final_fit$bic)) {
      cat(round(x$final_fit$bic, 4), "\n", sep = "")
    } else {
      cat("NA\n")
    }
  }
  
  invisible(x)
}

# Print heterogeneous effect test

print.heterogeneous_test <- function(x, ...) {
  cat("GMR heterogeneous-effect test\n")
  cat("-----------------------------\n")
  cat("Covariate: ", x$covariate, "\n", sep = "")
  cat("Alpha:     ", x$alpha, "\n", sep = "")
  cat("p0:        ")
  
  if (is.numeric(x$p0) && length(x$p0) == 1 && is.finite(x$p0)) {
    cat(signif(x$p0, 4), "\n", sep = "")
  } else {
    cat("NA\n")
  }
  
  cat("Reject H0: ", x$reject, "\n", sep = "")
  cat("\nNull model:\n")
  cat("  Formula: ")
  print(x$null_formula)
  cat("  Common:  ")
  print(x$null_common)
  
  cat("\nAlternative model:\n")
  cat("  Formula: ")
  print(x$alt_formula)
  cat("  Common:  ")
  print(x$alt_common)
  
  invisible(x)
}