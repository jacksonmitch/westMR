
rm(list = ls())
source("R/utils.R")
source("R/loglik.R")
source("R/E_step.R")
source("R/M_step_qr.R")
source("R/M_step_sqr.R")
source("R/initialization.R")
source("R/model_criteria.R")
source("R/EM_gmr.R")
source("R/fit_gmr_fixed.R")
source("R/fit_gmr.R")

source("R/formula_functions.R")
source("R/west_procedure.R")
source("R/test_heterogeneous_effect.R")
source("R/select_effects.R")
source("R/effect_selection_table.R")
source("R/methods.R")
library(microbenchmark)

simulate_effect_data <- function(n=400, seed=123){

  set.seed(seed)
  n <- n

  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)

  pi_true <- c(0.4,0.6)
  z <- sample.int(2, size = n, replace = TRUE, prob=pi_true)

  beta0 <- c(-1.5,1.5)
  beta1 <- c(-3,3)
  beta2 <- 2
  beta3 <- 0.5

  sigma <- c(0.5,0.5)

  mu <- beta0[z] + beta1[z]*x1 + beta2*x2 + beta3*x3
  y <- rnorm(n, mu, sigma[z])

  data.frame(
    y = y,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    true_group = z
  )
}

data <- simulate_effect_data()

true_heterogeneous <- 'x1'
true_homogeneous <- c('x2','x3')
selected_covariates <- c('x1','x2','x3')

model <- GMRModel$new(
  formula = make_formula(predictors = selected_covariates, response = 'y'),
  G_values = 1:2,
  data = data,
  method = 'sqr',
  alpha = 0.05,
  n_init = 5,
  n_kmeans_init = 5,
  maxit = 150,
  #seed = seed,
  fit_function = fit_gmr
)

forward_fit <- select_effects(
  model = model,
  direction = "forward"
)

bench <- microbenchmark(
  forward_fit = select_effects(
    model = model,
    direction = "forward"
  ),
  backward_fit = select_effects(
    model = model,
    direction = "backward"
  ),
  times = 1
)

print(bench)

print(forward_fit)
summary(forward_fit$final_fit)


print(backward_fit)
summary(backward_fit$final_fit)
steps_df <- effect_selection_table(forward_fit)
steps_df


