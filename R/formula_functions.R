
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
  as.character(formula[[2]])
}

get_predictors <- function(formula) {
  attr(stats::terms(formula), "term.labels")
}
