make_formula <- function(predictors, response = NULL) {
  if (length(predictors) == 0 || is.null(predictors)) {
    predictors <- "1"
  }
  stats::reformulate(
    termlabels = predictors,
    response = response
  )
}

get_response <- function(formula) {
  deparse(formula[[2]])
}

get_predictors <- function(formula) {
  attr(stats::terms(formula), "term.labels")
}

format_formula <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }

  paste(deparse(x), collapse = " ")
}

make_effects_formula <- function(heterogeneous, homogeneous, response = NULL) {
  common_part <- if (length(homogeneous) == 0) NULL else paste(homogeneous, collapse = " + ")
  het_part <- if (length(heterogeneous) == 0) NULL else paste0("(", paste(heterogeneous, collapse = " + "), " | group)")

  rhs_terms <- c(common_part, het_part)

  if (length(rhs_terms) == 0) rhs_terms <- "1"

  rhs <- paste(rhs_terms, collapse = " + ")

  if (!is.null(response)) {
    stats::as.formula(paste(response, "~", rhs))
  } else {
    stats::as.formula(paste("~", rhs))
  }
}