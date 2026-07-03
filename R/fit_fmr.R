# Fit a single model specification across all G values.
# Returns a list of fit_fmr objects, one per G.

fit_across_G <- function(model, prepared_data, extra_tau_starts = NULL) {
  G_values <- model$G_values

  features <- make_initialization_features(
    prepared_data = prepared_data,
    family = model$family
  )

  fits <- lapply(seq_along(G_values), function(i) {
    G <- G_values[[i]]

    init_list <- make_tau_list(
      prepared_data = prepared_data,
      G = G,
      control = model$control,
      family = model$family,
      features = features
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
#' Fit One Model Specification at a Given G
#'
#' Fits a finite mixture regression model for a fixed number of components
#' \code{G} and a fixed design (heterogeneous/common predictor split). Burns
#' in each initialization in \code{init_list} for a few EM iterations,
#' selects the best-starting initialization by log-likelihood via
#' \code{select_best_initialization()}, then runs \code{em_fmr()} to
#' convergence from that state.
#'
#' @param model A \code{WMRModel} object.
#' @param G An integer number of mixture components.
#' @param init_list A named list of candidate \code{tau} initialization
#'   matrices (n x G), as produced by \code{make_tau_list()}.
#' @param prepared_data A \code{WMRData} object (from \code{prepare_data()})
#'   for this model specification.
#'
#' @return A list of class \code{fit_fmr} with elements including
#'   \code{parameter_values} (a labeled snapshot of the fitted EM state:
#'   \code{tau}, \code{pi_g}, \code{beta_g}, \code{beta}, \code{sigma_g},
#'   \code{eta}, \code{loglik}), \code{loglik}, \code{loglik_trace},
#'   \code{iterations}, \code{converged}, \code{irwls_iterations}/
#'   \code{irwls_converged} (non-Gaussian families only; dropped entirely
#'   for Gaussian), initialization diagnostics (\code{best_init_name},
#'   \code{best_init_loglik}, \code{n_valid_init}, \code{init}),
#'   \code{bic}, \code{num_parameters}, \code{family}, \code{G},
#'   \code{n_init}, and \code{control}. Any \code{NULL}-valued element is
#'   dropped from the result via \code{compact()}.
#' @noRd
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

  k <- count_params_fmr(
    ncol_het = prepared_data$p_het,
    ncol_common = prepared_data$p_com,
    G = G,
    family = family
  )

  bic <- compute_bic(loglik = best_fit$loglik, n = prepared_data$n, k = k)

  out <- compact(list(
    parameter_values = em_state$to_list(prepared_data),
    loglik = best_fit$loglik,
    loglik_trace = best_fit$loglik_trace,
    iterations = best_fit$iterations,
    converged = best_fit$converged,
    irwls_iterations = best_fit$irwls_iterations,
    irwls_converged = best_fit$irwls_converged,
    best_init_name = init_fit$best_name,
    best_init_loglik = init_fit$best_loglik,
    bic = bic,
    num_parameters = k,
    family = family,
    G = G,
    n_init = length(init_list),
    n_valid_init = sum(is.finite(init_fit$logliks)),
    init = list(
      strategy = "multistart_tau_burnin",
      burnin = control$init_burnin,
      n_starts = length(init_list),
      best_start = init_fit$best_name,
      burnin_logliks = init_fit$logliks,
      failed_starts = init_fit$failures
    ),
    control = control
  ))

  class(out) <- "fit_fmr"
  out
}
