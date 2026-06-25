simulate_fmr <- function(n, betas, pi = NULL, sigma = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  betas <- as.matrix(betas)
  G <- nrow(betas)
  if (is.null(colnames(betas))) stop("betas must have column names")

  if (is.null(pi)) pi <- rep(1 / G, G)
  sigma <- rep_len(sigma, G)
  z <- sample.int(G, size = n, replace = TRUE, prob = pi / sum(pi))

  pred_names <- setdiff(colnames(betas), "x0")
  X <- matrix(rnorm(n * length(pred_names)),
    nrow = n,
    dimnames = list(NULL, pred_names)
  )

  Xfull <- cbind("x0" = 1, X)

  mu <- rowSums(Xfull * betas[z, colnames(Xfull), drop = FALSE])
  y <- rnorm(n, mu, sigma[z])

  cbind(data.frame(y = y), as.data.frame(X)) # data.frame(true_group = z)
}

scenarios <- list(
  two_group_effects = list(
    betas = rbind(
      g1 = c("x0" = -1.5, x1 = -3, x2 = 2, x3 = 0.5),
      g2 = c("x0" = 1.5, x1 = 3, x2 = 2, x3 = 0.5)
    ),
    pi = c(0.4, 0.6),
    sigma = c(0.5, 0.5)
  ),
  two_group_subtle = list(
    betas = rbind(
      g1 = c("x0" = 0, x1 = -1, x2 = 1.5),
      g2 = c("x0" = 0, x1 = 1, x2 = 1.5)
    ),
    pi = c(0.5, 0.5),
    sigma = c(1.0, 1.0)
  ),
  three_group_subtle = list(
    betas = rbind(
      g1 = c("x0" = 34, x1 = 0.21, x2 = 0),
      g2 = c("x0" = 37, x1 = 0.01, x2 = 0),
      g3 = c("x0" = 28, x1 = 0.09, x2 = 0)
    ),
    pi = c(0.7, 0.1, 0.2),
    sigma = c(6, 4, 8)
  ),
  three_group_four_variables = list(
    betas = rbind(
      g1 = c("x0" = -3, x1 = 0.3, x2 = -0.4, x3 = 0.25, x4 = 0),
      g2 = c("x0" = 0, x1 = 0.6, x2 = -0.4, x3 = 0.25, x4 = 0),
      g3 = c("x0" = 3, x1 = 0.8, x2 = -0.4, x3 = 0.25, x4 = 0)
    ),
    pi = c(0.2, 0.3, 0.5),
    sigma = c(0.5, 0.5, 0.5)
  ),
  unbalanced = list(
    betas = rbind(
      g1 = c("x0" = -3, x1 = -2, x2 = 0.5),
      g2 = c("x0" = 1, x1 = 1, x2 = 0.5)
    ),
    pi = c(0.1, 0.9),
    sigma = c(0.5, 0.5)
  ),
  het_variance = list(
    betas = rbind(
      g1 = c("x0" = -2, x1 = -2, x2 = 1.0),
      g2 = c("x0" = 2, x1 = 2, x2 = 1.0)
    ),
    pi = c(0.5, 0.5),
    sigma = c(0.3, 1.5)
  ),
  four_group = list(
    betas = rbind(
      g1 = c("x0" = -4, x1 = -2, x2 = 1.0, x3 = -0.5),
      g2 = c("x0" = -1, x1 = -0.5, x2 = 1.0, x3 = -0.5),
      g3 = c("x0" = 1, x1 = 0.5, x2 = 1.0, x3 = -0.5),
      g4 = c("x0" = 4, x1 = 2, x2 = 1.0, x3 = -0.5)
    ),
    pi = c(0.2, 0.3, 0.3, 0.2),
    sigma = c(0.5, 0.5, 0.5, 0.5)
  )
)
