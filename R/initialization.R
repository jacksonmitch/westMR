make_tau_from_partition <- function(partition, G, eps = 1e-6) {
  partition <- as.integer(partition)
  n <- length(partition)

  if (G == 1L) {
    return(matrix(1, nrow = n, ncol = 1L))
  }

  if (length(unique(partition)) != G) {
    stop("partition must contain exactly G groups.")
  }

  partition <- match(partition, sort(unique(partition)))

  tau <- matrix(eps / (G - 1L), nrow = n, ncol = G)
  tau[cbind(seq_len(n), partition)] <- 1 - eps

  colnames(tau) <- paste0("g", seq_len(G))

  tau
}

random_partition_balanced <- function(n, G) {
  partition <- rep(seq_len(G), length.out = n)
  sample(partition, size = n, replace = FALSE)
}

rank_partition <- function(score, G) {
  score <- as.numeric(score)
  n <- length(score)

  ord <- order(score)
  partition <- integer(n)

  group_sizes <- rep(floor(n / G), G)
  remainder <- n %% G

  if (remainder > 0L) {
    group_sizes[seq_len(remainder)] <- group_sizes[seq_len(remainder)] + 1L
  }

  start <- 1L

  for (g in seq_len(G)) {
    end <- start + group_sizes[g] - 1L
    partition[ord[start:end]] <- g
    start <- end + 1L
  }

  partition
}

make_initialization_features <- function(prepared_data,
                                         family = c("gaussian", "poisson", "binomial")) {
  family <- match.arg(family)

  A <- as.matrix(prepared_data$X_het)
  B <- as.matrix(prepared_data$X_com)
  y <- as.numeric(prepared_data$y)

  X <- cbind(A, B)
  n <- length(y)

  if (family == "gaussian") {
    fit0 <- stats::lm.fit(x = X, y = y)
    fitted0 <- as.numeric(fit0$fitted.values)
    resid0 <- as.numeric(fit0$residuals)

    return(cbind(
      response = y,
      fitted = fitted0,
      residual = resid0
    ))
  }

  if (family == "poisson") {
    fit0 <- stats::glm.fit(x = X, y = y, family = stats::poisson())

    eta0 <- as.numeric(fit0$linear.predictors)
    mu0 <- as.numeric(fit0$fitted.values)

    resid0 <- (y - mu0) / sqrt(pmax(mu0, 1e-8))

    return(cbind(
      response = log(y + 0.5),
      eta = eta0,
      residual = resid0
    ))
  }

  if (family == "binomial") {
    binomial_size <- prepared_data$binomial_size

    if (is.null(binomial_size)) {
      binomial_size <- rep(1, n)
    }

    binomial_size <- as.numeric(binomial_size)
    y_bin <- cbind(y, binomial_size - y)

    fit0 <- stats::glm.fit(
      x = X,
      y = y_bin,
      family = stats::binomial()
    )

    eta0 <- as.numeric(X %*% fit0$coefficients)
    mu0 <- as.numeric(fit0$fitted.values)

    p0 <- (y + 0.5) / (binomial_size + 1)
    p0 <- pmin(pmax(p0, 1e-8), 1 - 1e-8)

    logit <- stats::qlogis(p0)

    resid0 <- (y - binomial_size * mu0) /
      sqrt(pmax(binomial_size * mu0 * (1 - mu0), 1e-8))

    return(cbind(
      response = logit,
      eta = eta0,
      residual = resid0
    ))
  }
}

make_tau_list <- function(prepared_data,
                          G,
                          control,
                          family = c("gaussian", "poisson", "binomial")) {
  family <- match.arg(family)

  n <- prepared_data$n

  if (G == 1L) {
    tau <- matrix(1, nrow = n, ncol = 1L)
    colnames(tau) <- "g1"
    return(list(single = tau))
  }

  features <- make_initialization_features(
    prepared_data = prepared_data,
    family = family
  )

  features_scaled <- scale(features)
  features_scaled[!is.finite(features_scaled)] <- 0

  min_size <- max(2L, prepared_data$p_het + prepared_data$p_com + 1L)
  min_size <- min(min_size, floor(n / G))

  tau_list <- list()

  part_response <- rank_partition(features[, "response"], G)

  if (min(table(part_response)) >= min_size) {
    tau_list$quantile_response <- make_tau_from_partition(
      partition = part_response,
      G = G,
      eps = control$init_eps
    )
  }

  part_resid <- rank_partition(features[, "residual"], G)

  if (min(table(part_resid)) >= min_size) {
    tau_list$quantile_residual <- make_tau_from_partition(
      partition = part_resid,
      G = G,
      eps = control$init_eps
    )
  }

  for (s in seq_len(control$n_kmeans_init)) {
    km <- try(
      stats::kmeans(
        x = features_scaled,
        centers = G,
        nstart = control$kmeans_starts
      ),
      silent = TRUE
    )

    if (!inherits(km, "try-error") &&
      length(unique(km$cluster)) == G &&
      min(table(km$cluster)) >= min_size) {
      tau_list[[paste0("kmeans_", s)]] <- make_tau_from_partition(
        partition = km$cluster,
        G = G,
        eps = control$init_eps
      )
    }
  }

  n_random <- max(0L, control$n_init - control$n_kmeans_init)

  for (s in seq_len(n_random)) {
    part_random <- random_partition_balanced(n, G)

    tau_list[[paste0("random_", s)]] <- make_tau_from_partition(
      partition = part_random,
      G = G,
      eps = control$init_eps
    )
  }
  if (family == "binomial" && prepared_data$p_het > 0L) {
    A <- as.matrix(prepared_data$X_het)

    for (j in seq_len(ncol(A))) {
      part_x <- rank_partition(A[, j], G)

      if (min(table(part_x)) >= min_size) {
        tau_list[[paste0("quantile_A", j)]] <- make_tau_from_partition(
          partition = part_x,
          G = G,
          eps = control$init_eps
        )
      }

      part_neg_x <- rank_partition(-A[, j], G)

      if (min(table(part_neg_x)) >= min_size) {
        tau_list[[paste0("quantile_neg_A", j)]] <- make_tau_from_partition(
          partition = part_neg_x,
          G = G,
          eps = control$init_eps
        )
      }
    }
  }

  if (length(tau_list) == 0L) {
    stop("No valid initialization starts were created.")
  }

  tau_list
}

select_best_initialization <- function(tau_list,
                                       prepared_data,
                                       G,
                                       family,
                                       control) {
  best_fit <- NULL
  best_loglik <- -Inf
  best_name <- NA_character_

  logliks <- rep(NA_real_, length(tau_list))
  names(logliks) <- names(tau_list)

  failures <- character(0)

  control$tol <- 0
  control$max_iter <- control$init_burnin

  for (s in seq_along(tau_list)) {
    em_state <- EmState$new(tau = tau_list[[s]])
    fit_s <- try(
      em_fmr(
        prepared_data = prepared_data,
        G = G,
        em_state = em_state,
        family = family,
        control = control
      ),
      silent = TRUE
    )

    if (inherits(fit_s, "try-error")) {
      failures <- c(failures, names(tau_list)[s])
      next
    }

    logliks[s] <- fit_s$loglik

    if (is.finite(fit_s$loglik) && fit_s$loglik > best_loglik) {
      best_fit <- fit_s
      best_loglik <- fit_s$loglik
      best_name <- names(tau_list)[s]
    }
  }

  if (is.null(best_fit)) {
    stop("All initialization attempts failed.")
  }

  list(
    best_fit = best_fit,
    best_loglik = best_loglik,
    best_name = best_name,
    logliks = logliks,
    failures = failures
  )
}
