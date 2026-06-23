format_formula_clean <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }

  if (inherits(x, "formula")) {
    return(paste(deparse(x), collapse = " "))
  }

  paste(deparse(x), collapse = " ")
}

format_none <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return("none")
  }

  paste(x, collapse = ", ")
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

####
# westMR methods
####

#' @export
print.westMR <- function(x, ...) {
  cat("westMR result\n")
  cat("-------------\n")
  cat("Family:     ", x$family, "\n", sep = "")
  cat("Task:       ", x$task, "\n", sep = "")
  cat("G values:   ", paste(x$G_values, collapse = ", "), "\n", sep = "")

  if (!is.null(x$selected_G)) {
    cat("Selected G: ", x$selected_G, "\n", sep = "")
  }

  cat("\n")

  if (!is.null(x$variable_selection)) {
    cat("Variable selection:\n")
    cat("  Selected predictors: ",
        format_none(x$variable_selection$selected),
        "\n",
        sep = "")
  }

  if (!is.null(x$effect_determination)) {
    cat("\nEffect-type determination:\n")
    cat("  Heterogeneous effects: ",
        format_none(x$effect_determination$heterogeneous),
        "\n",
        sep = "")
    cat("  Homogeneous effects:   ",
        format_none(x$effect_determination$homogeneous),
        "\n",
        sep = "")
  }

  if (is.null(x$effect_determination) &&
      !is.null(x$variable_selection) &&
      length(x$variable_selection$selected) == 0L) {
    cat("\nEffect-type determination was skipped because no predictors were selected.\n")
  }

  invisible(x)
}

#' @export
summary.westMR <- function(object, ...) {
  variable_summary <- NULL
  effect_summary <- NULL

  if (!is.null(object$variable_selection)) {
    variable_summary <- list(
      direction = object$variable_selection$direction,
      alpha = object$variable_selection$alpha,
      G_values = object$variable_selection$G_values,
      selected_G = object$variable_selection$selected_G,
      selected = object$variable_selection$selected,
      final_formula = object$variable_selection$final_formula,
      n_steps = length(object$variable_selection$steps)
    )
  }

  if (!is.null(object$effect_determination)) {
    effect_summary <- list(
      direction = object$effect_determination$direction,
      alpha = object$effect_determination$alpha,
      G_values = object$effect_determination$G_values,
      selected_G = object$effect_determination$selected_G,
      heterogeneous = object$effect_determination$heterogeneous,
      homogeneous = object$effect_determination$homogeneous,
      final_formula = object$effect_determination$final_formula,
      final_common = object$effect_determination$final_common,
      n_steps = length(object$effect_determination$steps)
    )
  }

  out <- list(
    call = object$call,
    formula = object$formula,
    family = object$family,
    task = object$task,
    G_values = object$G_values,
    selected_G = object$selected_G,
    best_fit = object$best_fit,
    variable_selection = variable_summary,
    effect_determination = effect_summary
  )

  class(out) <- "summary.westMR"
  out
}

#' @export
print.summary.westMR <- function(x, ...) {
  cat("Summary of westMR result\n")
  cat("------------------------\n")
  cat("Family:     ", x$family, "\n", sep = "")
  cat("Task:       ", x$task, "\n", sep = "")
  cat("G values:   ", paste(x$G_values, collapse = ", "), "\n", sep = "")

  if (!is.null(x$selected_G)) {
    cat("Selected G: ", x$selected_G, "\n", sep = "")
  }

  cat("\n")

  if (!is.null(x$variable_selection)) {
    cat("Variable selection\n")
    cat("------------------\n")
    cat("Direction:     ", x$variable_selection$direction, "\n", sep = "")
    cat("Alpha:         ", x$variable_selection$alpha, "\n", sep = "")
    cat("Steps:         ", x$variable_selection$n_steps, "\n", sep = "")
    cat("Selected G:    ", x$variable_selection$selected_G, "\n", sep = "")
    cat("Selected:      ", format_none(x$variable_selection$selected), "\n", sep = "")
    cat("Final formula: ",
        format_formula_clean(x$variable_selection$final_formula),
        "\n\n",
        sep = "")
  }

  if (!is.null(x$effect_determination)) {
    cat("Effect-type determination\n")
    cat("-------------------------\n")
    cat("Direction:     ", x$effect_determination$direction, "\n", sep = "")
    cat("Alpha:         ", x$effect_determination$alpha, "\n", sep = "")
    cat("Steps:         ", x$effect_determination$n_steps, "\n", sep = "")
    cat("Selected G:    ", x$effect_determination$selected_G, "\n", sep = "")
    cat("Heterogeneous: ",
        format_none(x$effect_determination$heterogeneous),
        "\n",
        sep = "")
    cat("Homogeneous:   ",
        format_none(x$effect_determination$homogeneous),
        "\n",
        sep = "")
    cat("Final formula: ",
        format_formula_clean(x$effect_determination$final_formula),
        "\n",
        sep = "")
  }

  invisible(x)
}

####
# fit_fmr methods
####

#' @export
print.fit_fmr <- function(x, ...) {
  cat("westMR mixture regression fit\n")
  cat("------------------------------\n")

  if (!is.null(x$family)) {
    cat("Family:         ", x$family, "\n", sep = "")
  }

  if (!is.null(x$G)) {
    cat("G:              ", x$G, "\n", sep = "")
  }

  if (!is.null(x$loglik)) {
    cat("Log-likelihood: ", round(x$loglik, 4), "\n", sep = "")
  }

  if (!is.null(x$bic)) {
    cat("BIC:            ", round(x$bic, 4), "\n", sep = "")
  }

  if (!is.null(x$converged)) {
    cat("Converged:      ", x$converged, "\n", sep = "")
  }

  if (!is.null(x$iterations)) {
    cat("Iterations:     ", x$iterations, "\n", sep = "")
  }

  if (!is.null(x$best_init)) {
    cat("Best init:      ", x$best_init, "\n", sep = "")
  }

  if (!is.null(x$n_valid_init)) {
    cat("Valid inits:    ", x$n_valid_init, "\n", sep = "")
  }

  if (!is.null(x$pi_g)) {
    cat("\nMixing proportions:\n")
    pi_named <- x$pi_g
    names(pi_named) <- paste0("g", seq_along(pi_named))
    print(round(pi_named, 4))
  }

  if (!is.null(x$sigma_g)) {
    cat("\nComponent standard deviations:\n")
    sigma_named <- x$sigma_g
    names(sigma_named) <- paste0("g", seq_along(sigma_named))
    print(round(sigma_named, 4))
  }

  if (!is.null(x$beta_g)) {
    cat("\nHeterogeneous coefficients:\n")
    print(round(x$beta_g, 4))
  }

  if (!is.null(x$beta) && length(x$beta) > 0L) {
    cat("\nHomogeneous coefficients:\n")
    print(round(x$beta, 4))
  }

  invisible(x)
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
      return("all homogeneous")
    }
    if (length(common) == 0) {
      return("all heterogeneous")
    }
    paste0(
      "het: ", paste(het, collapse = ", "),
      "  |  homo: ", paste(common, collapse = ", ")
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
