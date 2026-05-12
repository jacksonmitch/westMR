# Fit Gaussian mixture regression for fixed G

fit_gmr_fixed <- function(A,
                    B,
                    y,
                    G,
                    method = c("sqr", "qr"),
                    maxit = 300,
                    tol = 1e-6,
                    n_init = 10,
                    n_kmeans_init = 2,
                    kmeans_starts = 20,
                    start_list = NULL,
                    verbose = FALSE,
                    sigma_floor = NULL){
  
  method <- match.arg(method)
  
  A <- as.matrix(A)
  B <- as.matrix(B)
  y <- as.numeric(y)
  
  n <- nrow(A)
  
  # Input checks
  
  if (length(y) != n) {
    stop("length(y) must equal nrow(A).")
  }
  
  if (nrow(B) != n) {
    stop("nrow(B) must equal nrow(A).")
  }
  
  if (G < 1) {
    stop("G must be at least 1.")
  }
  
  if (n_init < 1) {
    stop("n_init must be at least 1.")
  }
  
  if (n_kmeans_init < 0) {
    stop("n_kmeans_init cannot be negative.")
  }
  
  if (n_kmeans_init > n_init) {
    stop("n_kmeans_init cannot be larger than n_init.")
  }
  
  if (maxit < 1) {
    stop("maxit must be at least 1.")
  }
  
  if (tol <= 0) {
    stop("tol must be positive.")
  }
  
  if (anyNA(A) || anyNA(B) || anyNA(y)) {
    stop("A, B, and y cannot contain missing values.")
  }
  
  # Create initial clusters
  
  if (is.null(start_list)) {
    
    start_list <- make_start_list(
      y = y,
      G = G,
      n_init = n_init,
      n_kmeans_init = n_kmeans_init,
      kmeans_starts = kmeans_starts
    )
    
  } else {
    
    if (!is.list(start_list)) {
      stop("start_list must be a list of initial clustering vectors.")
    }
    
    if (length(start_list) != n_init) {
      stop("length(start_list) must equal n_init.")
    }
  }
  
  # Run EM from each initialization
  
  fits <- vector("list", n_init)
  logliks <- rep(-Inf, n_init)
  errors <- rep(NA_character_, n_init)
  
  best_fit <- NULL
  best_loglik <- -Inf
  best_init <- NA_integer_
  n_valid_init <- 0L
  
  for (s in seq_len(n_init)) {
    
    if (verbose) {
      message("Fitting G = ", G, ", initialization ", s, " of ", n_init)
    }
    
    fit_s <- tryCatch({
      
      # Convert one initial clustering into starting parameter values
      init_s <- initialize_parameters(
        A = A,
        B = B,
        y = y,
        G = G,
        method = method,
        cl_init = start_list[[s]],
        sigma_floor = sigma_floor
      )
      
      # Run EM from these starting values
      em_s <- em_gmr(
        A = A,
        B = B,
        y = y,
        beta_g = init_s$beta_g,
        beta = init_s$beta,
        sigma_g = init_s$sigma,
        pi_g = init_s$pi,
        method = method,
        max_iter = maxit,
        tol = tol,
        sigma_floor = sigma_floor,
        verbose = FALSE
      )
      
      em_s$error_msg <- NA_character_
      em_s
      
    }, error = function(e) {
      
      list(
        beta_g = NULL,
        beta = NULL,
        sigma_g = rep(NA_real_, G),
        pi_g = rep(NA_real_, G),
        tau = matrix(NA_real_, nrow = n, ncol = G),
        loglik = -Inf,
        loglik_trace = NA_real_,
        iterations = NA_integer_,
        converged = FALSE,
        error_msg = conditionMessage(e)
      )
    })
    
    fits[[s]] <- fit_s
    logliks[s] <- fit_s$loglik
    errors[s] <- fit_s$error_msg
    
    if (is.finite(fit_s$loglik)) {
      n_valid_init <- n_valid_init + 1L
      
      if (fit_s$loglik > best_loglik) {
        best_fit <- fit_s
        best_loglik <- fit_s$loglik
        best_init <- s
      }
    }
    
    if (verbose) {
      message(
        "Finished initialization ", s,
        ": loglik = ", round(fit_s$loglik, 6),
        ", converged = ", fit_s$converged
      )
    }
  }
  
  # If all initializations failed
  
  k <- count_params_gmr(
    A = A,
    B = B,
    G = G
  )
  
  if (is.null(best_fit)) {
    
    out <- list(
      G = G,
      method = method,
      loglik = -Inf,
      bic = Inf,
      k = k,
      beta_g = NULL,
      beta = NULL,
      sigma_g = rep(NA_real_, G),
      pi_g = rep(NA_real_, G),
      tau = matrix(NA_real_, nrow = n, ncol = G),
      iterations = NA_integer_,
      converged = FALSE,
      best_init = NA_integer_,
      n_init = n_init,
      n_kmeans_init = n_kmeans_init,
      kmeans_starts = kmeans_starts,
      n_valid_init = n_valid_init,
      all_logliks = logliks,
      error_msgs = errors,
      all_fits = fits,
      call = match.call()
    )
    
    class(out) <- "gmr_fit_fixed"
    
    return(out)
  }
    
  # Compute BIC for best fit
  
  bic <- compute_bic(
    loglik = best_fit$loglik,
    n = n,
    k = k
  )
  
  # Return clean fitted object
  
  out <- list(
    G = G,
    method = method,
    
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
    
    best_init = best_init,
    n_init = n_init,
    n_kmeans_init = n_kmeans_init,
    kmeans_starts = kmeans_starts,
    n_valid_init = n_valid_init,
    all_logliks = logliks,
    error_msgs = errors,
    all_fits = fits,
    
    call = match.call()
  )
  
  class(out) <- "gmr_fit_fixed"
  
  out
}