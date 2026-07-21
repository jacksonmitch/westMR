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

  pred_names <- setdiff(colnames(betas), "Intercept")
  X <- matrix(rnorm(n * length(pred_names)),
    nrow = n,
    dimnames = list(NULL, pred_names)
  )

  Xfull <- cbind("Intercept" = 1, X)

  eta <- rowSums(Xfull * betas[z, colnames(Xfull), drop = FALSE])

  y <- switch(family,
    gaussian = rnorm(n, mean = eta, sd = sigma[z]),
    poisson = rpois(n, lambda = exp(eta)),
    binomial = rbinom(n, size = size, prob = stats::plogis(eta))
  )

  out <- list(
    data = cbind(data.frame(y = y), as.data.frame(X)),
    formula = make_formula(pred_names, "y"),
    true_group = z,
    size = if (family == "binomial") rep_len(size, n) else NULL
  )
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
      g1 = c("Intercept" = -1.5, x1 = -3, x2 = 2, x3 = 0.5),
      g2 = c("Intercept" = 1.5, x1 = 3, x2 = 2, x3 = 0.5)
    ),
    pi = c(0.4, 0.6),
    sigma = c(0.5, 0.5)
  ),
  three_group_four_variables = list(
    betas = rbind(
      g1 = c("Intercept" = -3, x1 = 0.3, x2 = -0.4, x3 = 0.25, x4 = 0),
      g2 = c("Intercept" = 0, x1 = 0.6, x2 = -0.4, x3 = 0.25, x4 = 0),
      g3 = c("Intercept" = 3, x1 = 0.8, x2 = -0.4, x3 = 0.25, x4 = 0)
    ),
    pi = c(0.2, 0.3, 0.5),
    sigma = c(0.5, 0.5, 0.5)
  ),
  two_group_effects_poisson = list(
    betas = rbind(
      g1 = c("Intercept" = -0.5, x1 = -1.0, x2 = 0.3, x3 = 0.1),
      g2 = c("Intercept" = 0.5, x1 = 1.0, x2 = 0.3, x3 = 0.1)
    ),
    pi = c(0.4, 0.6),
    family = "poisson"
  ),
  two_group_effects_binomial = list(
    betas = rbind(
      g1 = c("Intercept" = -1.0, x1 = -1.5, x2 = 0.5, x3 = 0.2),
      g2 = c("Intercept" = 1.0, x1 = 1.5, x2 = 0.5, x3 = 0.2)
    ),
    pi = c(0.5, 0.5),
    family = "binomial",
    size = 25
  ),
  three_group_twelve_variables_gaussian = list(
    betas = rbind(
      g1 = c(
        "Intercept" = -2, x1 = 0.4, x2 = 0.2, x3 = -0.5,
        x4 = -1.1, x5 = 1.0, x6 = 0.2, x7 = -0.3,
        x8 = 0.3, x9 = -0.4, x10 = 0.2, x11 = 0,
        x12 = 0
      ),
      g2 = c(
        "Intercept" = 0, x1 = 0.8, x2 = 0.5, x3 = 0,
        x4 = -0.5, x5 = 1.4, x6 = 0.3, x7 = -0.2,
        x8 = 0.3, x9 = -0.4, x10 = 0.2, x11 = 0,
        x12 = 0
      ),
      g3 = c(
        "Intercept" = 2, x1 = 1.1, x2 = 0.6, x3 = 0.4,
        x4 = 0.1, x5 = 1.8, x6 = 0.4, x7 = -0.1,
        x8 = 0.3, x9 = -0.4, x10 = 0.2, x11 = 0,
        x12 = 0
      )
    ),
    pi = c(0.2, 0.5, 0.3),
    sigma = c(0.5, 0.5, 0.5)
  ),
  four_group_twelve_variables = list(
    betas = rbind(
      g1 = c(
        "Intercept" = -2, het1 = 0.4, het2 = 0.2,
        het3 = -0.5, het4 = -1.1, het5 = 1.0,
        hom1 = 0.6, hom2 = -0.3, hom3 = 0.2,
        null1 = 0, null2 = 0, null3 = 0, null4 = 0
      ),
      g2 = c(
        "Intercept" = 0, het1 = 0.8, het2 = 0.5, 
        het3 = 0, het4 = -0.5, het5 = 1.4,
        hom1 = 0.6, hom2 = -0.3, hom3 = 0.2,
        null1 = 0, null2 = 0, null3 = 0, null4 = 0
      ),
      g3 = c(
        "Intercept" = 2, het1 = 1.1, het2 = 0.6, 
        het3 = 0.4, het4 = 0.1, het5 = 1.8,
        hom1 = 0.6, hom2 = -0.3, hom3 = 0.2,
        null1 = 0, null2 = 0, null3 = 0, null4 = 0
      ),
      g4 = c(
        "Intercept" = 1, het1 = 0.9, het2 = 0.3, 
        het3 = 0.7, het4 = -0.3, het5 = -1.6,
        hom1 = 0.6, hom2 = -0.3, hom3 = 0.2,
        null1 = 0, null2 = 0, null3 = 0, null4 = 0
      )
    ),
    pi = c(0.2, 0.3, 0.25, 0.25),
    sigma = c(0.6, 0.5, 0.4, 0.6)
  )
)
