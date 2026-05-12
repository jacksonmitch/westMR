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

source("R/make_formula.R")
source("R/west_test.R")
source("R/test_heterogeneous_effect.R")
source("R/select_effects_forward.R")
source("R/select_effects_backward.R")
source("R/effect_selection_table.R")
source("R/methods.R")


rm(list = ls())
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
  y <- mu + rnorm(n, mu, sigma[z])
  
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

forward_fit <- select_effects_forward(
  response = data$y,
  covariates = selected_covariates,
  data = data,
  G_values = 1:2,
  alpha = 0.05,
  method = 'sqr',
  n_init = 5,
  n_kmeans_init = 5,
  maxit = 150
)


