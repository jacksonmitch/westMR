# Fit Gaussian mixture regression and select number of components by BIC
# Controller of the numeric functions, sets methods to be used and G's used

# TODO: @included = predictors included in the model
# @common = predictors that will be fir under an homogeneous effect
# @G_values = all candidate orders fitted
fit_fmr <- function(model,
                    G,
                    included = model$predictors,
                    common = NULL) {

  control <- model$control
  # Build response and component-specific design matrix A
  mf <- stats::model.frame(
    formula = model$formula,
    data = model$data,
    na.action = stats::na.fail
  )
  response <- as.numeric(stats::model.response(mf))

  heterogeneous_formula <- make_formula(
    predictors = setdiff(model$predictors, common),
    response = model$response)

  het_predictors_matrix <- stats::model.matrix(
    object = heterogeneous_formula,
    data = mf
  )

  # Build common-effect design matrix B

  # TODO: better logic here?
  common_predictors_matrix <- NULL

  if (!is.null(common) && length(common) != 0) {
    common_formula <- stats::as.formula(
      paste("~", paste(common, collapse = " + "), "- 1")
    )
    common_predictors_matrix <- stats::model.matrix(
      object = common_formula,
      data = mf
    )
    if (nrow(common_predictors_matrix) != nrow(het_predictors_matrix)) {
      stop("nrow(B) must equal nrow(A).")
    }
  }

  # if (length(response) != nrow(het_predictors_matrix)) {
  #   stop("length(y) must equal nrow(A).")
  # }

  # if (anyNA(het_predictors_matrix) || anyNA(common_predictors_matrix) || anyNA(response)) {
  #   stop("A, B, and y cannot contain missing values.")
  # }

  k <- count_params_gmr(
    ncol_het = length(included) - length(common),
    ncol_common = length(common),
    G = G
  )

  start_list <- make_start_list(
    y = response,
    G_length = G,
    control = control
  )

  n <- nrow(het_predictors_matrix)
  n_init <- model$control$n_init
  fits <- vector("list", n_init)
  logliks <- rep(-Inf, n_init)

  best_fit <- NULL
  best_loglik <- -Inf
  best_init <- NA_integer_
  n_valid_init <- 0L
  for (i in seq_len(n_init)) {

    Z <- matrix(0, nrow = n, ncol = G)
    Z[cbind(seq_len(n), start_list[[i]])] <- 1

    fit <- em_fmr(
      het = het_predictors_matrix,
      common = common_predictors_matrix,
      response = response,
      G = G,
      tau = Z,
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

  bic <- compute_bic(
    loglik = best_fit$loglik,
    n = n,
    k = k
  )

  out <- list(
    best_fit = best_fit,
    fits = fits,

    # Direct access to selected model estimates
    beta_g = best_fit$beta_g,
    beta = best_fit$beta,
    sigma_g = best_fit$sigma_g,
    pi_g = best_fit$pi_g,
    tau = best_fit$tau,

    loglik = best_fit$loglik,
    loglik_trace = best_fit$loglik_trace,
    bic = bic,
    k = k,

    iterations = best_fit$iterations,
    converged = best_fit$converged,

    best_init = best_fit$best_init,
    n_valid_init = best_fit$n_valid_init,

    call = match.call()
  )

  class(out) <- "fit_fmr"

  out
}
