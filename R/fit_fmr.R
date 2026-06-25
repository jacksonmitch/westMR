# Fit a single model specification across all G values.
# Returns a list of fit_fmr objects, one per G.

fit_across_G <- function(model, prepared_data, extra_tau_starts = NULL) {
  G_values <- model$G_values

  fits <- lapply(seq_along(G_values), function(i) {
    G <- G_values[[i]]

    init_list <- make_tau_list(
      prepared_data = prepared_data,
      G = G,
      control = model$control,
      family = model$family
    )

    if (!is.null(extra_tau_starts)) {
      extra_tau <- extra_tau_starts[[i]]

      init_list <- add_extra_tau_start(
        init_list = init_list,
        tau = extra_tau,
        n = prepared_data$n,
        G = G,
        name = "from_shared_fit"
      )
    }

    fit_fmr(
      model = model,
      G = G,
      init_list = init_list,
      prepared_data = prepared_data
    )
  })

  names(fits) <- paste0("G", G_values)

  fits
}

# allows shared fit to be passed as initialization

add_extra_tau_start <- function(init_list,
                                tau,
                                n,
                                G,
                                name = "extra_tau") {
  if (is.null(tau)) {
    return(init_list)
  }

  tau <- as.matrix(tau)

  if (!identical(dim(tau), c(n, G))) {
    return(init_list)
  }

  if (any(!is.finite(tau)) || any(tau < 0)) {
    return(init_list)
  }

  rs <- rowSums(tau)

  if (any(!is.finite(rs)) || any(rs <= 0)) {
    return(init_list)
  }

  tau <- tau / rs
  colnames(tau) <- paste0("g", seq_len(G))

  init_list[[name]] <- tau

  init_list
}


# Find the best fit given a specific G and a list of initializations
fit_fmr <- function(model,
                    G,
                    init_list,
                    prepared_data) {
  control <- model$control
  family <- model$family

  # Short burn-in stage
  # run a small number of EM iterations from each tau start,
  # then choose the start with the largest burn-in log-likelihood
  init_fit <- select_best_initialization(
    tau_list = init_list,
    prepared_data = prepared_data,
    G = G,
    family = family,
    control = control
  )

  em_state <- init_fit$best_fit$em_state

  best_fit <- em_fmr(
    prepared_data = prepared_data,
    G = G,
    em_state = em_state,
    family = family,
    control = control
  )

  k <- count_params_gmr(
    ncol_het = prepared_data$p_het,
    ncol_common = prepared_data$p_com,
    G = G,
    family = family
  )

  bic <- compute_bic(loglik = best_fit$loglik, n = prepared_data$n, k = k)

  add_parameter_labels(em_state, prepared_data)

  out <- list(
    em_state = em_state,
    loglik = best_fit$loglik,
    loglik_trace = best_fit$loglik_trace,
    iterations = best_fit$iterations,
    converged = best_fit$converged,
    irwls_iterations = best_fit$irwls_iterations,
    irwls_converged = best_fit$irwls_converged,
    best_init_name = init_fit$best_name,
    best_init_loglik = init_fit$best_loglik,
    n_valid_init = sum(is.finite(init_fit$logliks)),
    bic = bic,
    k = k,
    family = family,
    G = G,
    n_init = length(init_list),
    logliks = init_fit$logliks,
    init = list(
      strategy = "multistart_tau_burnin",
      burnin = control$init_burnin,
      n_starts = length(init_list),
      best_start = init_fit$best_name,
      burnin_logliks = init_fit$logliks,
      failed_starts = init_fit$failures
    ),
    model = model,
    call = match.call()
  )

  class(out) <- "fit_fmr"
  out
}

add_parameter_labels <- function(em_state, prepared_data) {
  if (!is.null(em_state$beta_g)) {
    bg <- em_state$beta_g
    rownames(bg) <- paste0("g", seq_len(nrow(bg)))
    if (!is.null(colnames(prepared_data$X_het))) {
      colnames(bg) <- colnames(prepared_data$X_het)
    }
    em_state$beta_g <- bg
  }

  if (!is.null(em_state$beta) && length(em_state$beta) > 0L) {
    b <- as.numeric(em_state$beta)
    if (!is.null(colnames(prepared_data$X_com))) {
      names(b) <- colnames(prepared_data$X_com)
    }
    em_state$beta <- b
  }

  if (!is.null(em_state$sigma_g)) {
    sg <- em_state$sigma_g
    names(sg) <- paste0("g", seq_along(sg))
    em_state$sigma_g <- sg
  }

  if (!is.null(em_state$pi_g)) {
    pg <- em_state$pi_g
    names(pg) <- paste0("g", seq_along(pg))
    em_state$pi_g <- pg
  }

  if (!is.null(em_state$tau)) {
    tau <- em_state$tau
    colnames(tau) <- paste0("g", seq_len(ncol(tau)))
    em_state$tau <- tau
  }
}
