# Class for storing all constant variables (variables, data, control)

WMRModel <- R6::R6Class(
  "WMRModel",

  public = list(
    formula = NULL,
    response = NULL,
    predictors = NULL,
    data = NULL,
    G_values = NULL,
    m_step_var = NULL,
    m_step_eff = NULL,
    log_lik = NULL,
    control = NULL,

    initialize = function(formula,
                          data,
                          G_values,
                          m_step_var,
                          m_step_eff,
                          log_lik,
                          control) {
      self$formula <- formula
      self$response <- get_response(formula)
      self$predictors <- get_predictors(formula)
      self$data <- data
      self$G_values <- G_values
      self$m_step_var <- m_step_var
      self$m_step_eff <- m_step_eff
      self$log_lik <- log_lik
      self$control <- control
    }
  )
)


