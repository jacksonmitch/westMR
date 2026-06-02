simulate_effect_data <- function(n, seed) {
  set.seed(seed)
  n <- n

  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)

  pi_true <- c(0.4, 0.6)
  z <- sample.int(2, size = n, replace = TRUE, prob = pi_true)

  beta0 <- c(-1.5, 1.5)
  beta1 <- c(-3, 3)
  beta2 <- 2
  beta3 <- 0.5

  sigma <- c(0.5, 0.5)

  mu <- beta0[z] + beta1[z] * x1 + beta2 * x2 + beta3 * x3
  y <- rnorm(n, mu, sigma[z])

  data.frame(
    y = y,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    true_group = z
  )
}

#
# useful_thing <- readRDS(test_path("fixtures", "useful_thing1.rds"))
