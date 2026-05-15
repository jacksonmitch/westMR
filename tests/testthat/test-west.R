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

model <- GMRModel$new(
  formula = make_formula(predictors = selected_covariates, response = 'y'),
  G_values = 1:2,
  data = data,
  method = 'sqr',
  alpha = 0.05,
  n_init = 1,
  n_kmeans_init = 5,
  maxit = 150,
  #seed = seed,
  fit_function = fit_gmr
)

forward_fit <- select_effects(
  model = model,
  direction = "forward"
)
print(forward_fit)
summary(forward_fit$final_fit)

backward_fit <- select_effects(
  model = model,
  direction = "backward"
)
print(backward_fit)
summary(backward_fit$final_fit)
steps_df <- effect_selection_table(backward_fit)
steps_df


