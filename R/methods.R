# S3 print methods

#' @export
print.fit_fmr <- function(x, ...) {
  cat("Finite mixture regression fit\n")
  cat("------------------------------\n")
  cat("G:              ", x$best_fit$G, "\n", sep = "")
  cat("Log-likelihood: ", round(x$loglik, 4), "\n", sep = "")
  cat("Converged:      ", x$converged, "\n", sep = "")
  cat("Iterations:     ", x$iterations, "\n", sep = "")
  cat("Best init:      ", x$best_init, "\n", sep = "")
  cat("Valid inits:    ", x$n_valid_init, "\n", sep = "")
  cat("\n")

  cat("Mixing proportions:\n")
  pi_named <- x$pi_g
  names(pi_named) <- paste0("g", seq_along(pi_named))
  print(round(pi_named, 4))

  cat("\nComponent standard deviations:\n")
  sigma_named <- x$sigma_g
  names(sigma_named) <- paste0("g", seq_along(sigma_named))
  print(round(sigma_named, 4))

  cat("\nComponent-specific coefficients:\n")
  beta_g <- x$beta_g
  rownames(beta_g) <- paste0("g", seq_len(nrow(beta_g)))
  print(round(beta_g, 4))

  if (length(x$beta) > 0) {
    cat("\nCommon-effect coefficients:\n")
    print(round(x$beta, 4))
  }

  invisible(x)
}

#' @export
print.determine_effects <- function(x, ...) {
  cat("GMR effect-type determination\n")
  cat("------------------------------\n")
  cat("Direction: ", x$direction, "\n", sep = "")
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
  cat("  Formula: ", format(x$final_formula), "\n", sep = "")
  cat("  Common:  ", format(x$final_common), "\n", sep = "")

  cat("\nStep summary:\n")
  step_df <- do.call(rbind, lapply(x$steps, function(s) {
    data.frame(
      step    = s$step,
      chosen  = if (is.na(s$chosen)) "-" else s$chosen,
      min_p0  = round(min(s$p0, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  }))
  print(step_df, row.names = FALSE)

  invisible(x)
}

#' @export
print.variable_selection <- function(x, ...) {
  cat("GMR variable selection\n")
  cat("----------------------\n")
  cat("Direction: ", x$direction, "\n", sep = "")
  cat("Alpha:     ", x$alpha, "\n", sep = "")
  cat("G values:  ", paste(x$G_values, collapse = ", "), "\n", sep = "")
  cat("\n")

  cat("Selected variables:\n")
  if (length(x$selected) == 0) {
    cat("  None\n")
  } else {
    cat("  ", paste(x$selected, collapse = ", "), "\n", sep = "")
  }

  cat("\nFinal model formula: ", format(x$final_formula), "\n", sep = "")

  cat("\nStep summary:\n")
  step_df <- do.call(rbind, lapply(x$steps, function(s) {
    data.frame(
      step    = s$step,
      chosen  = if (is.na(s$chosen)) "-" else s$chosen,
      min_p0  = round(min(s$p0, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  }))
  print(step_df, row.names = FALSE)

  invisible(x)
}
