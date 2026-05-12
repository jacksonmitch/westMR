# Fit Gaussian mixture regression and select number of components by BIC

fit_gmr <- function(formula,
                    common = ~ 0,
                    data,
                    G_max = 5,
                    G_values = NULL,
                    method = c("sqr", "qr"),
                    maxit = 300,
                    tol = 1e-6,
                    n_init = 10,
                    n_kmeans_init = 2,
                    kmeans_starts = 20,
                    verbose = FALSE,
                    sigma_floor = NULL) {
  
  method <- match.arg(method)
  
  # Input checks
  
  if (missing(formula)) {
    stop("formula must be supplied.")
  }
  
  if (!inherits(formula, "formula")) {
    stop("formula must be a formula, such as y ~ x1 + x2.")
  }
  
  if (missing(data)) {
    stop("data must be supplied.")
  }
  
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  
  if (!inherits(common, "formula")) {
    stop("common must be a one-sided formula, such as ~ z1 + z2 or ~ 0.")
  }
  
  if (length(common) != 2) {
    stop("common must be a one-sided formula, such as ~ z1 + z2 or ~ 0.")
  }
  
  # Build response and component-specific design matrix A
  
  mf <- stats::model.frame(
    formula = formula,
    data = data,
    na.action = stats::na.fail
  )
  
  y <- stats::model.response(mf)
  
  if (!is.numeric(y)) {
    stop("The response variable must be numeric.")
  }
  
  A <- stats::model.matrix(
    object = formula,
    data = mf
  )
  
  y <- as.numeric(y)
  
  # Build common-effect design matrix B
  
  common_terms <- attr(stats::terms(common), "term.labels")
  
  if (length(common_terms) == 0) {
    
    B <- matrix(
      numeric(0),
      nrow = nrow(A),
      ncol = 0
    )
    
  } else {
    
    common_formula <- stats::as.formula(
      paste("~", paste(common_terms, collapse = " + "), "- 1")
    )
    
    B <- stats::model.matrix(
      object = common_formula,
      data = data
    )
  }
  
  n <- nrow(A)
  
  if (length(y) != n) {
    stop("length(y) must equal nrow(A).")
  }
  
  if (nrow(B) != n) {
    stop("nrow(B) must equal nrow(A).")
  }
  
  if (anyNA(A) || anyNA(B) || anyNA(y)) {
    stop("A, B, and y cannot contain missing values.")
  }
  
  # Determine candidate G values
  
  if (is.null(G_values)) {
    
    if (length(G_max) != 1 || G_max < 1) {
      stop("G_max must be a single integer greater than or equal to 1.")
    }
    
    G_values <- seq_len(as.integer(G_max))
    
  } else {
    
    G_values <- sort(unique(as.integer(G_values)))
    
    if (length(G_values) < 1) {
      stop("G_values must contain at least one value.")
    }
    
    if (any(G_values < 1)) {
      stop("All values in G_values must be at least 1.")
    }
  }
  
  # Check fitting controls
  
  if (length(n_init) != 1 || n_init < 1) {
    stop("n_init must be a single integer greater than or equal to 1.")
  }
  
  if (length(n_kmeans_init) != 1 || n_kmeans_init < 0) {
    stop("n_kmeans_init must be a single nonnegative integer.")
  }
  
  if (n_kmeans_init > n_init) {
    stop("n_kmeans_init cannot be larger than n_init.")
  }
  
  if (length(kmeans_starts) != 1 || kmeans_starts < 1) {
    stop("kmeans_starts must be a single integer greater than or equal to 1.")
  }
  
  if (length(maxit) != 1 || maxit < 1) {
    stop("maxit must be a single integer greater than or equal to 1.")
  }
  
  if (length(tol) != 1 || tol <= 0) {
    stop("tol must be a single positive number.")
  }
  
  n_init <- as.integer(n_init)
  n_kmeans_init <- as.integer(n_kmeans_init)
  kmeans_starts <- as.integer(kmeans_starts)
  maxit <- as.integer(maxit)
  
  # Fit each candidate G
  
  fits <- vector("list", length(G_values))
  names(fits) <- paste0("G", G_values)
  
  results <- data.frame(
    G = G_values,
    loglik = NA_real_,
    bic = NA_real_,
    k = NA_integer_,
    iterations = NA_integer_,
    converged = NA,
    n_valid_init = NA_integer_
  )
  
  for (i in seq_along(G_values)) {
    
    G <- G_values[i]
    
    if (verbose) {
      message("Fitting candidate model with G = ", G)
    }
    
    fit_i <- fit_gmr_fixed(
      A = A,
      B = B,
      y = y,
      G = G,
      method = method,
      maxit = maxit,
      tol = tol,
      n_init = n_init,
      n_kmeans_init = n_kmeans_init,
      kmeans_starts = kmeans_starts,
      start_list = NULL,
      verbose = verbose,
      sigma_floor = sigma_floor
    )
    
    fits[[i]] <- fit_i
    
    results$loglik[i] <- fit_i$loglik
    results$bic[i] <- fit_i$bic
    results$k[i] <- fit_i$k
    results$iterations[i] <- fit_i$iterations
    results$converged[i] <- fit_i$converged
    results$n_valid_init[i] <- fit_i$n_valid_init
  }
  
  # Select best model by BIC
  
  if (all(!is.finite(results$bic))) {
    
    best_idx <- NA_integer_
    best_fit <- NULL
    best_G <- NA_integer_
    
  } else {
    
    best_idx <- which.min(results$bic)
    best_fit <- fits[[best_idx]]
    best_G <- results$G[best_idx]
  }
  
  if (is.null(best_fit)) {
    
    out <- list(
      best_fit = NULL,
      best_G = NA_integer_,
      criterion = "BIC",
      results = results,
      fits = fits,
      G_values = G_values,
      method = method,
      
      beta_g = NULL,
      beta = NULL,
      sigma_g = NULL,
      pi_g = NULL,
      tau = NULL,
      
      loglik = -Inf,
      bic = Inf,
      k = NA_integer_,
      iterations = NA_integer_,
      converged = FALSE,
      
      call = match.call()
    )
    
    class(out) <- "gmr_fit"
    
    return(out)
  }

  out <- list(
    best_fit = best_fit,
    best_G = best_G,
    criterion = "BIC",
    results = results,
    fits = fits,
    G_values = G_values,
    method = method,
    
    # Direct access to selected model estimates
    beta_g = best_fit$beta_g,
    beta = best_fit$beta,
    sigma_g = best_fit$sigma_g,
    pi_g = best_fit$pi_g,
    tau = best_fit$tau,
    
    loglik = best_fit$loglik,
    loglik_trace = best_fit$loglik_trace,
    bic = best_fit$bic,
    k = best_fit$k,
    
    iterations = best_fit$iterations,
    converged = best_fit$converged,
    
    best_init = best_fit$best_init,
    n_init = n_init,
    n_kmeans_init = n_kmeans_init,
    kmeans_starts = kmeans_starts,
    n_valid_init = best_fit$n_valid_init,
    
    call = match.call()
  )
  
  class(out) <- "gmr_fit"
  
  out
}