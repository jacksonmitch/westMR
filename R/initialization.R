# Initialization helpers

make_init_clustering <- function(y,
                                 G,
                                 method = c("kmeans", "random", "random_balanced"),
                                 kmeans_starts = 20) {
  method <- match.arg(method)
  y <- as.numeric(y)
  n <- length(y)

  if (G == 1) {
    return(rep(1L, n))
  }

  if (method == "kmeans") {

    # kmeans cannot use more centers than distinct response values.
    # This is especially important for Bernoulli/binomial responses.
    if (length(unique(y)) < G) {
      method <- "random_balanced"
    } else {
      return(stats::kmeans(
        x = y,
        centers = G,
        nstart = kmeans_starts
      )$cluster)
    }
  }

  if (method == "random") {
    repeat {
      cl <- sample.int(G, n, replace = TRUE)

      if (length(unique(cl)) == G) {
        return(cl)
      }
    }
  }

  if (method == "random_balanced") {
    cl <- rep(seq_len(G), length.out = n)
    cl <- sample(cl, size = n, replace = FALSE)

    return(cl)
  }
}

# Converts a cluster assignment vector to a 0/1 tau matrix (n x G)

clustering_to_tau <- function(cl, G) {
  n <- length(cl)
  Z <- matrix(0, nrow = n, ncol = G)
  Z[cbind(seq_len(n), cl)] <- 1
  Z
}


# Returns a list of n_init tau matrices ready to pass to fit_fmr.
# The first n_kmeans_init are kmeans-seeded; the rest are random_balanced.


make_tau_list <- function(y, G, control, family = c("gaussian", "poisson", "binomial")) {
  family <- match.arg(family)

  n_init <- control$n_init
  n_kmeans_init <- control$n_kmeans_init
  kmeans_starts <- control$kmeans_starts
  n <- length(y)

  lapply(seq_len(n_init), function(i) {

    if (G == 1) {
      return(matrix(1, nrow = n, ncol = 1))
    }

    if (family == "binomial") {
      return(make_random_hard_tau(n = n, G = G))
    }

    method <- if (i <= n_kmeans_init) {
      "kmeans"
    } else {
      "random_balanced"
    }

    cl <- make_init_clustering(
      y = y,
      G = G,
      method = method,
      kmeans_starts = kmeans_starts
    )

    clustering_to_tau(cl, G)
  })
}

make_random_hard_tau <- function(n, G) {
  assignments <- sample.int(G, size = n, replace = TRUE)

  # Ensure every component is represented
  while (length(unique(assignments)) < G) {
    assignments <- sample.int(G, size = n, replace = TRUE)
  }

  tau <- matrix(0, nrow = n, ncol = G)
  tau[cbind(seq_len(n), assignments)] <- 1
  tau
}

# Create initial EM parameter values from an initial clustering

initialize_parameters <- function(A,
                                  B,
                                  y,
                                  G_length,
                                  method,
                                  kmeans_starts,
                                  cl_init,
                                  sigma_floor) {
  n <- nrow(A)

  # Convert hard clustering to indicator matrix Z
  Z <- matrix(0, nrow = n, ncol = G_length)
  Z[cbind(seq_len(n), cl_init)] <- 1

  pi_g <- colMeans(Z)

  # Use one M-step from the initial hard clustering to get starting coefficients
  if (method == "sqr") {
    fit0 <- m_step_sqr(
      A = A,
      B = B,
      y = y,
      tau = Z,
      sigma_floor = sigma_floor
    )
  } else if (method == "qr") {
    fit0 <- m_step_qr(
      A = A,
      B = B,
      y = y,
      tau = Z,
      sigma_floor = sigma_floor
    )
  }

  beta_g <- fit0$beta_g
  beta <- fit0$beta

  # Match the old code: compute sigma from the initial hard clustering
  mu <- sweep(
    A %*% t(beta_g),
    1,
    common_eta(B, beta),
    "+"
  )

  res2 <- (y - mu)^2

  sigma_g <- sqrt(
    pmax(
      1e-8,
      colSums(Z * res2) / pmax(1e-8, colSums(Z))
    )
  )

  if (is.null(sigma_floor)) {
    sigma_floor <- 0.05 * stats::sd(y)
  }

  sigma_g <- pmax(sigma_floor, sigma_g)

  list(
    beta_g = beta_g,
    beta = beta,
    pi = pi_g,
    sigma = sigma_g,
    tau = Z,
    cl_init = cl_init,
    method = method
  )
}
