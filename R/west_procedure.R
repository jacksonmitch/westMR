west_procedure <- function(
    model,
    null_predictors = model$predictors,
    null_common = NULL,
    alt_predictors = model$predictors,
    alt_common = NULL
){

  G_values <- model$G_values
  alpha <- model$alpha

  rows <- vector("list", length(G_values))
  null_fits <- vector("list", length(G_values))
  alt_fits <- vector("list", length(G_values))

  for(i in seq_along(G_values)){
    G <- G_values[[i]]

    if (model$control$verbose) {
      message("Fitting candidate model with G = ", G)
    }

    fit_null <- fit_fmr(model = model,
                        G = G,
                        included = null_predictors,
                        common = null_common)
    fit_alt <- fit_fmr(model = model,
                       G = G,
                       included = alt_predictors,
                       common = alt_common)

    # Add error checking that fit_null and fit_alt are correct

    lrt <- -2*(fit_null$loglik - fit_alt$loglik)

    df <- fit_alt$k - fit_null$k

    if (model$control$verbose) cat("lrt :", lrt, " df: ", df, "\n")
    if (is.finite(lrt) && is.finite(df) && df > 0) {
      p_value <- stats::pchisq(lrt, df = df, lower.tail = FALSE)
    }
    else {
      p_value <- NA_real_
    }

    null_fits[[i]] <- fit_null
    alt_fits[[i]] <- fit_alt

    rows[[i]] <- data.frame(
      G = G,
      loglik_null = fit_null$loglik,
      loglik_alt = fit_alt$loglik,
      bic_null = fit_null$bic,
      bic_alt = fit_alt$bic,
      k_null = fit_null$k,
      k_alt = fit_alt$k,
      lrt = lrt,
      df = df,
      p_value = p_value,
      null_converged = fit_null$converged,
      alt_converged = fit_alt$converged,
      stringsAsFactors = FALSE)
  }

  tab <- do.call(rbind, rows)

  valid <- is.finite(tab$bic_null) & is.finite(tab$p_value)

  if (!any(valid)) {
    p0 <- NA_real_
    weights <- rep(NA_real_, nrow(tab))
    reject <- FALSE
  } else {
    # Stable version of exp(-0.5 * BIC) weights.
    bic_valid <- tab$bic_null[valid]
    shifted <- -0.5 * (bic_valid - min(bic_valid))
    w_valid <- exp(shifted)
    w_valid <- w_valid / sum(w_valid)

    weights <- rep(NA_real_, nrow(tab))
    weights[valid] <- w_valid

    p0 <- sum(tab$p_value[valid] * weights[valid])
    reject <- is.finite(p0) && p0 < alpha
  }

  tab$weight <- weights

  out <- list(
    p0 = p0,
    reject = reject,
    alpha = alpha,
    table = tab,
    null_fits = null_fits,
    alt_fits = alt_fits,
    G_values = G_values
  )

  class(out) <- "west_procedure"
  out

}
