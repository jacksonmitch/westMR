# ---- select_variables -------------------------------------------------------

#' Print a Variable Selection Result
#'
#' @param x An object of class \code{select_variables}, as returned by
#'   \code{select_variables()}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @export
print.select_variables <- function(x, ...) {
  cat("Variable selection\n")
  cat("  Direction: ", x$direction,
    "  |  Alpha: ", x$alpha, "\n",
    sep = ""
  )
  cat("  Selected:      ", format_none(x$selected), "\n", sep = "")
  cat("  Final formula: ", format_formula(x$final_formula), "\n", sep = "")
  invisible(x)
}

# ---- determine_effects ------------------------------------------------------

#' Print an Effect-Type Determination Result
#'
#' @param x An object of class \code{determine_effects}, as returned by
#'   \code{determine_effects()}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @export
print.determine_effects <- function(x, ...) {
  cat("Effect-type determination\n")
  cat("  Direction: ", x$direction,
    "  |  Alpha: ", x$alpha, "\n",
    sep = ""
  )
  cat("  Heterogeneous: ", format_none(x$heterogeneous), "\n", sep = "")
  cat("  Homogeneous:   ", format_none(x$homogeneous), "\n", sep = "")
  cat("  Final formula: ", format_formula(x$final_formula), "\n", sep = "")
  invisible(x)
}

# ---- steps_table ------------------------------------------------------------

#' Print the step-by-step selection table
#'
#' @param x An object of class \code{determine_effects} or
#'   \code{select_variables}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @export
steps_table <- function(x, ...) UseMethod("steps_table")

#' Print the Step-by-Step Effect-Type Determination Table
#'
#' @param x An object of class \code{determine_effects}, as returned by
#'   \code{determine_effects()}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
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
  print_steps(x$steps, union(x$heterogeneous, x$homogeneous), x$direction, label_state)
  invisible(x)
}

#' Print the Step-by-Step Variable Selection Table
#'
#' @param x An object of class \code{select_variables}, as returned by
#'   \code{select_variables()}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
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
      in_set <- if (direction == "forward") c(in_set, s$chosen) else setdiff(in_set, s$chosen)
    } else {
      cat("  No eligible candidate -- stopping.\n")
    }
    cat("\n")
  }
}

# ---- fit_fmr ----------------------------------------------------------------

#' Print a Fitted Finite Mixture Regression Model
#'
#' @param x An object of class \code{fit_fmr}, as returned by
#'   \code{fit_fmr()}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @export
print.fit_fmr <- function(x, ...) {
  cat("Best fit (G: ", x$G, ")\n", sep = "")

  if (!is.null(x$parameter_values$pi_g)) {
    cat("\nMixing proportions:\n")
    print(round(x$parameter_values$pi_g, 4))
  }

  if (!is.null(x$parameter_values$sigma_g)) {
    cat("\nComponent standard deviations:\n")
    print(round(x$parameter_values$sigma_g, 4))
  }

  if (!is.null(x$parameter_values$beta_g)) {
    cat("\nHeterogeneous coefficients:\n")
    print(round(x$parameter_values$beta_g, 4))
  }

  if (!is.null(x$parameter_values$beta) && length(x$em_state$beta) > 0L) {
    cat("\nHomogeneous coefficients:\n")
    print(round(x$parameter_values$beta, 4))
  }

  invisible(x)
}

# ---- westMR -----------------------------------------------------------------

#' Print a westMR Result
#'
#' @param x An object of class \code{westMR}, as returned by
#'   \code{\link{westMR}}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @export
print.westMR <- function(x, ...) {
  cat("westMR  |  family: ", x$family,
    "  |  G: ", paste(x$G_values, collapse = ":"),
    "  |  task: ", x$task, "\n\n",
    sep = ""
  )

  if (!is.null(x$variable_selection)) {
    print(x$variable_selection)
    cat("\n")
  }

  if (!is.null(x$effect_determination)) {
    print(x$effect_determination)
    cat("\n")
  }

  if (is.null(x$effect_determination) &&
    !is.null(x$variable_selection) &&
    length(x$variable_selection$selected) == 0L) {
    cat("Effect-type determination skipped: no predictors selected.\n\n")
  }

  if (!is.null(x$best_fit)) {
    print(x$best_fit)
  }

  invisible(x)
}

#' Summarize a westMR Result
#'
#' @param object An object of class \code{westMR}, as returned by
#'   \code{\link{westMR}}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{summary.westMR}: a list containing the
#'   call, formula, family, task, candidate \code{G} values, the best fit,
#'   the variable-selection and effect-determination results (if run), and
#'   their step-by-step tables.
#' @export
summary.westMR <- function(object, ...) {
  out <- list(
    call = object$call,
    formula = object$formula,
    family = object$family,
    task = object$task,
    G_values = object$G_values,
    best_fit = object$best_fit,
    selected_G = object$best_fit$G,
    variable_selection = object$variable_selection,
    variable_steps = object$variable_selection$steps,
    effect_determination = object$effect_determination,
    effect_steps = object$effect_determination$steps
  )
  class(out) <- "summary.westMR"
  out
}

#' Print a westMR Summary
#'
#' @param x An object of class \code{summary.westMR}, as returned by
#'   \code{summary.westMR()}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @export
print.summary.westMR <- function(x, ...) {
  cat("westMR summary  |  family: ", x$family,
    "  |  G: ", paste(x$G_values, collapse = ","),
    "  |  task: ", x$task, "\n\n",
    sep = ""
  )

  if (!is.null(x$variable_selection)) {
    print(x$variable_selection)
    cat("\n")
    steps_table(x$variable_selection)
  }

  if (!is.null(x$effect_determination)) {
    print(x$effect_determination)
    cat("\n")
    steps_table(x$effect_determination)
  }

  if (is.null(x$effect_determination) &&
    !is.null(x$variable_selection) &&
    length(x$variable_selection$selected) == 0L) {
    cat("Effect-type determination skipped: no predictors selected.\n\n")
  }

  if (!is.null(x$best_fit)) {
    print(x$best_fit)
  }

  invisible(x)
}
