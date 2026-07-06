simulate_fmr <- function(n,
                         betas,
                         pi = NULL,
                         sigma = 1,
                         family = c("gaussian", "poisson", "binomial"),
                         size = 1,
                         seed = NULL) {
  family <- match.arg(family)

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

  eta <- rowSums(Xfull * betas[z, colnames(Xfull), drop = FALSE])

  y <- switch(family,
    gaussian = rnorm(n, mean = eta, sd = sigma[z]),
    poisson = rpois(n, lambda = exp(eta)),
    binomial = rbinom(n, size = size, prob = stats::plogis(eta))
  )

  out <- cbind(data.frame(y = y), as.data.frame(X))

  if (family == "binomial") {
    out$size <- rep_len(size, n)
  }

  out
}

match_group_order <- function(fitted_betas, true_betas) {
  fitted <- as.matrix(fitted_betas)
  truth <- as.matrix(true_betas)

  G <- nrow(truth)
  stopifnot(nrow(fitted) == G, ncol(fitted) == ncol(truth))

  perms <- expand.grid(rep(list(seq_len(G)), G))
  valid <- apply(perms, 1, function(r) length(unique(r)) == G)
  perms <- perms[valid, , drop = FALSE]

  resid <- apply(perms, 1, function(perm) {
    sum((fitted[perm, , drop = FALSE] - truth)^2)
  })

  as.integer(perms[which.min(resid), ])
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
  ),
  two_group_effects_poisson = list(
    betas = rbind(
      g1 = c("x0" = -0.5, x1 = -1.0, x2 = 0.3, x3 = 0.1),
      g2 = c("x0" = 0.5, x1 = 1.0, x2 = 0.3, x3 = 0.1)
    ),
    pi = c(0.4, 0.6),
    family = "poisson"
  ),
  two_group_effects_binomial = list(
    betas = rbind(
      g1 = c("x0" = -1.0, x1 = -1.5, x2 = 0.5, x3 = 0.2),
      g2 = c("x0" = 1.0, x1 = 1.5, x2 = 0.5, x3 = 0.2)
    ),
    pi = c(0.5, 0.5),
    family = "binomial",
    size = 25
  )
)
