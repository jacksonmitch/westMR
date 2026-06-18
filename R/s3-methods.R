#' @export
print.fit_fmr <- function(x, ...) {
  cat("westMR mixture regression fit\n")
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

print_steps <- function(steps, all_predictors, direction, label_state) {
  in_set <- if (direction == "forward") character(0) else all_predictors

  for (s in steps) {
    out_set <- setdiff(all_predictors, in_set)
    state <- label_state(in_set, out_set)

    cat(sprintf("Step %d  (%s)\n", s$step, state))
    cat(sprintf("  %-20s  %-12s\n", "Candidate", "weighted p-value"))
    cat(sprintf("  %-20s  %-12s\n", "---------", "-------"))

    p_fmt <- formatC(s$p0, format = "e", digits = 2)
    for (nm in names(s$p0)) {
      marker <- if (!is.na(s$chosen) && nm == s$chosen) " *" else "  "
      cat(sprintf("  %-20s  %-10s%s\n", nm, p_fmt[[nm]], marker))
    }

    if (!is.na(s$chosen)) {
      cat(sprintf("  Chosen: %s\n", s$chosen))
      in_set <- if (direction == "forward") {
        c(in_set, s$chosen)
      } else {
        setdiff(in_set, s$chosen)
      }
    } else {
      cat("  No eligible candidate -- stopping.\n")
    }
    cat("\n")
  }
}

#' Print the step-by-step selection table
#'
#' @param x An object of class \code{effect_determination} or
#'  \code{variable_selection}.
#' @param ... Currently unused.
#' @export
#' @export
steps_table <- function(x, ...) UseMethod("steps_table")

#' @export
steps_table.determine_effects <- function(x, ...) {
  label_state <- function(het, common) {
    if (length(het) == 0) {
      return("all common")
    }
    if (length(common) == 0) {
      return("all heterogeneous")
    }
    paste0(
      "het: ", paste(het, collapse = ", "),
      "  |  common: ", paste(common, collapse = ", ")
    )
  }
  print_steps(
    x$steps, union(x$heterogeneous, x$homogeneous),
    x$direction, label_state
  )
  invisible(x)
}

#' @export
steps_table.select_variables <- function(x, ...) {
  label_state <- function(included, excluded) {
    if (length(included) == 0) {
      return("none included")
    }
    if (length(excluded) == 0) {
      return("all included")
    }
    paste0(
      "included: ", paste(included, collapse = ", "),
      "  |  excluded: ", paste(excluded, collapse = ", ")
    )
  }
  print_steps(x$steps, x$all_predictors, x$direction, label_state)
  invisible(x)
}

#' @export
print.determine_effects <- function(x, ...) {
  summary(x)
  steps_table(x)
  invisible(x)
}

#' @export
print.select_variables <- function(x, ...) {
  summary(x)
  steps_table(x)
  invisible(x)
}

#' @export
summary.determine_effects <- function(object, ...) {
  cat("westMR effect-type determination results\n")
  cat("---------------------------------\n")
  cat("Direction: ", object$direction, "\n", sep = "")
  cat("Alpha:     ", object$alpha, "\n", sep = "")
  cat("G values:  ", paste(object$G_values, collapse = ", "), "\n", sep = "")
  cat("\n")

  cat("Heterogeneous: ")
  cat(if (length(object$heterogeneous) == 0) {
    "none"
  } else {
    paste(object$heterogeneous, collapse = ", ")
  })
  cat("\nHomogeneous:   ")
  cat(if (length(object$homogeneous) == 0) {
    "none"
  } else {
    paste(object$homogeneous, collapse = ", ")
  })
  cat("\n\n")

  cat("Final model:\n")
  cat("  Formula: ", format(object$final_formula), "\n", sep = "")
  cat("  Common:  ", format(object$final_common), "\n", sep = "")
  cat("\n")
}

#' @export
summary.select_variables <- function(object, ...) {
  cat("westMR variable selection results\n")
  cat("--------------------------\n")
  cat("Direction: ", object$direction, "\n", sep = "")
  cat("Alpha:     ", object$alpha, "\n", sep = "")
  cat("G values:  ", paste(object$G_values, collapse = ", "), "\n", sep = "")
  cat("\n")

  cat("Selected: ")
  cat(if (length(object$selected) == 0) {
    "none"
  } else {
    paste(object$selected, collapse = ", ")
  })
  cat("\n\n")

  cat("Final formula: ", format(object$final_formula), "\n\n", sep = "")
}
