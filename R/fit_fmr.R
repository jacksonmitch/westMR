# Fit a single model specification across all G values.
# Returns a list of fit_fmr objects, one per G.

fit_across_G <- function(model, prepared_data) {
  G_values <- model$G_values

  lapply(G_values, function(G) {
    init_list <- make_tau_list(
      y = prepared_data$y,
      G = G,
      control = model$control
    )

    fit_fmr(
      model = model,
      G = G,
      init_list = init_list,
      prepared_data = prepared_data
    )
  })
}


# F
fit_fmr <- function(model,
                    G,
                    init_list,
                    prepared_data) {

  control <- model$control
  family <- model$family
  # Build response and component-specific design matrix A

  n_init <- length(init_list)
  fits <- vector("list", n_init)
  logliks <- rep(-Inf, n_init)

  best_fit <- NULL
  best_loglik <- -Inf
  best_init <- NA_integer_
  n_valid_init <- 0L

  for (i in seq_len(n_init)) {
    fit <- em_fmr(
      prepared_data = prepared_data,
      G = G,
      tau = init_list[[i]],
      family = model$family,
      control = control
    )

    fits[[i]] <- fit
    logliks[i] <- fit$loglik

    if (!is.finite(fit$loglik)) next

    n_valid_init <- n_valid_init + 1L

    if (fit$loglik > best_loglik) {
      best_fit <- fit
      best_loglik <- fit$loglik
      best_init <- i
    }
  }

  if (is.null(best_fit)) {
    stop("All fits failed for G = ", G)
  }

  out <- list(
    best_fit = best_fit,
    fits = fits,
    beta_g = best_fit$beta_g,
    beta = best_fit$beta,
    sigma_g = best_fit$sigma_g,
    pi_g = best_fit$pi_g,
    tau = best_fit$tau,
    loglik = best_fit$loglik,
    loglik_trace = best_fit$loglik_trace,
    iterations = best_fit$iterations,
    converged = best_fit$converged,
    best_init = best_init,
    n_valid_init = n_valid_init,

    family = family,
    G = G,
    n_init = n_init,
    n_valid_init = n_valid_init,
    logliks = logliks,

    irwls_iterations = best_fit$irwls_iterations,
    irwls_converged = best_fit$irwls_converged,

    call = match.call()
  )

  class(out) <- "fit_fmr"

  out
}
