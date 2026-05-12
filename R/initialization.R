# Initialization helpers

make_init_clustering <- function(y, G,
                                 method = c("kmeans", "random", "random_balanced"),
                                 kmeans_starts = 20) {
  method <- match.arg(method)
  n <- length(y)
  
  if (G == 1) {
    return(rep(1L, n))
  }
  
  if (method == "kmeans") {
    return(stats::kmeans(y, centers = G, nstart = kmeans_starts)$cluster)
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

make_start_list <- function(y,
                            G,
                            n_init = 5,
                            n_kmeans_init = 2,
                            kmeans_starts = 20) {
  starts <- vector("list", n_init)
  
  for (s in seq_len(n_init)) {
    starts[[s]] <- if (s <= n_kmeans_init) {
      make_init_clustering(
        y = y,
        G = G,
        method = "kmeans",
        kmeans_starts = kmeans_starts
      )
    } else {
      make_init_clustering(
        y = y,
        G = G,
        method = "random_balanced",
        kmeans_starts = kmeans_starts
      )
    }
  }
  
  starts
}


# Create initial EM parameter values from an initial clustering

initialize_parameters <- function(A,
                                  B,
                                  y,
                                  G,
                                  method = c("sqr", "qr"),
                                  init = c("kmeans", "random", "random_balanced"),
                                  kmeans_starts = 20,
                                  cl_init = NULL,
                                  sigma_floor = NULL) {
  method <- match.arg(method)
  init <- match.arg(init)
  
  A <- as.matrix(A)
  B <- as.matrix(B)
  y <- as.numeric(y)
  
  n <- nrow(A)
  
  if (length(y) != n) {
    stop("length(y) must equal nrow(A).")
  }
  
  if (nrow(B) != n) {
    stop("nrow(B) must equal nrow(A).")
  }
  
  if (G < 1) {
    stop("G must be at least 1.")
  }
  
  if (is.null(cl_init)) {
    cl_init <- make_init_clustering(
      y = y,
      G = G,
      method = init,
      kmeans_starts = kmeans_starts
    )
  }
  
  if (length(cl_init) != n) {
    stop("length(cl_init) must equal length(y).")
  }
  
  if (any(cl_init < 1) || any(cl_init > G)) {
    stop("cl_init must only contain values from 1 to G.")
  }
  
  # Convert hard clustering to indicator matrix Z
  Z <- matrix(0, nrow = n, ncol = G)
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
    init = init,
    method = method
  )
}