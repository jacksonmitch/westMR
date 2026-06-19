# Run the WEST procedure for one candidate predictor.
#
# shared_fits : output of fit_across_G for the shared side (null for forward,
#               alt for backward).
# direction   : "forward" or "backward",determines which role shared_fits plays
# candidate_prepared_data : prepared_data for the candidate side (already built
#                           by the caller for this specific predictor).

west_procedure <- function(
  model,
  direction,
  shared_fits,
  candidate_data
) {
  G_values <- model$G_values
  alpha <- model$control$alpha
  n <- candidate_data$n

  stopifnot(length(shared_fits) == length(G_values))

  extra_tau_starts <- lapply(shared_fits, function(fit) fit$tau)

  candidate_fits <- fit_across_G(
    model = model,
    prepared_data = candidate_data,
    extra_tau_starts = extra_tau_starts
  )

  rows <- vector("list", length(G_values))

  for (i in seq_along(G_values)) {
    G <- G_values[[i]]

    if (direction == "forward") {
      fit_null <- shared_fits[[i]]
      fit_alt <- candidate_fits[[i]]
    } else {
      fit_null <- candidate_fits[[i]]
      fit_alt <- shared_fits[[i]]
    }

    k_null <- count_params_gmr(
      ncol_het = ncol(fit_null$beta_g),
      ncol_common = length(fit_null$beta),
      G = G,
      family = model$family
    )
    k_alt <- count_params_gmr(
      ncol_het = ncol(fit_alt$beta_g),
      ncol_common = length(fit_alt$beta),
      G = G,
      family = model$family
    )

    bic_null <- compute_bic(fit_null$loglik, n, k_null)
    bic_alt <- compute_bic(fit_alt$loglik, n, k_alt)

    lrt <- -2 * (fit_null$loglik - fit_alt$loglik)
    df <- k_alt - k_null

    p_value <- if (is.finite(lrt) && is.finite(df) && df > 0) {
      stats::pchisq(lrt, df = df, lower.tail = FALSE)
    } else {
      NA_real_
    }

    rows[[i]] <- data.frame(
      G = G,
      loglik_null = fit_null$loglik,
      loglik_alt = fit_alt$loglik,
      bic_null = bic_null,
      bic_alt = bic_alt,
      k_null = k_null,
      k_alt = k_alt,
      lrt = lrt,
      df = df,
      p_value = p_value,
      null_converged = fit_null$converged,
      alt_converged = fit_alt$converged,
      stringsAsFactors = FALSE
    )
  }

  tab <- do.call(rbind, rows)

  if (direction == "forward") {
    model_bic <- tab$bic_alt
  }
  else {
    model_bic <- tab$bic_null
  }
  valid <- is.finite(model_bic) & is.finite(tab$p_value)

  if (!any(valid)) {
    p0 <- NA_real_
    weights <- rep(NA_real_, nrow(tab))
    reject <- FALSE
  } else {
    bic_valid <- model_bic[valid]
    shifted <- -0.5 * (bic_valid - min(bic_valid))
    w_valid <- exp(shifted) / sum(exp(shifted))

    weights <- rep(NA_real_, nrow(tab))
    weights[valid] <- w_valid

    p0 <- sum(tab$p_value[valid] * w_valid)
    reject <- is.finite(p0) && p0 < alpha
  }

  tab$weight <- weights

  out <- list(
    p0 = p0,
    reject = reject,
    alpha = alpha,
    table = tab,
    candidate_fits = candidate_fits,
    G_values = G_values
  )

  class(out) <- "west_procedure"
  out
}
