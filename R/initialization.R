# Initialization helpers

make_init_clustering <- function(y, G_length,
                                 cluster_method = c("kmeans", "random", "random_balanced"),
                                 kmeans_start) {
  cluster_method <- match.arg(cluster_method)
  n <- length(y)

  if (G_length == 1) {
    return(rep(1L, n))
  }

  if (cluster_method == "kmeans") {
    return(stats::kmeans(y, centers = G_length, nstart = kmeans_start)$cluster)
  }

  if (cluster_method == "random") {
    repeat {
      cl <- sample.int(G_length, n, replace = TRUE)

      if (length(unique(cl)) == G_length) {
        return(cl)
      }
    }
  }

  if (cluster_method == "random_balanced") {
    cl <- rep(seq_len(G_length), length.out = n)
    cl <- sample(cl, size = n, replace = FALSE)

    return(cl)
  }
}

make_start_list <- function(y,
                            G_length,
                            control) {
  n_init <- control$n_init
  n_kmeans_init <- control$n_kmeans_init
  kmeans_start <- control$kmeans_start

  starts <- vector("list", n_init)

  for (s in seq_len(n_init)) {
    starts[[s]] <- if (s <= n_kmeans_init) {
      make_init_clustering(
        y = y,
        G_length = G_length,
        cluster_method = "kmeans",
        kmeans_start = kmeans_start
      )
    } else {
      make_init_clustering(
        y = y,
        G_length = G_length,
        cluster_method = "random_balanced",
        kmeans_start = kmeans_start
      )
    }
  }
  starts
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
